// #######################################################################################
// Main setup and loops for performing well-tempered metadynamics on a
// diblock copolymer field-theoretic simulation
// #######################################################################################

#pragma once
#include <stdlib.h>
#include <cuda.h>
#include "GPUerror.h"
#include <iostream>
#include <fstream>
#include "diblock.h"
#include "anderson.h"
#include "metadynamics.h"
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <curand.h>

// Langevin update of w-(r) on the GPU using symmetrised noise
__global__ void langevin_sym(double *w_gpu, double *noise_gpu_new, double *noise_gpu_prev, const double chi_b, const double dt, const int M)
{
	int const tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= M) return;
	w_gpu[tid] += -(w_gpu[tid+2*M]+2*w_gpu[tid]/chi_b)*dt+0.5*(noise_gpu_prev[tid]+noise_gpu_new[tid]);
}

class wtmd_simulation {
    int    TpB_;                // GPU threads per block (default: 512)
    int    NA_;                 // Length of polymer A-block
    int    N_;                  // Total polymer length
    int    equil_its_;          // Number of Langevin steps for equilibration
    int    sim_its_;            // Number of Langevin steps for statistics
    int    sample_freq_;        // Number of Langevin steps between samples
    int    save_freq_;          // Number of steps between saving statistics to file
    int    *m_;                 // Number of mesh points [mx,my,mz]
    int    M_;                  // Total number of field mesh points
    int    Mk_;                 // Total number of field mesh points in k-space
    double chi_b_;              // Bare chi of the simulation
    double *L_;                 // Dimensions of simulation box [Lx,Ly,Lz] in units of aN^0.5
    double V_;                  // Volume of the simulation box
    double dt_;                 // Langevin time step multiplied by N
    double C_;                  // Dimensionless concentration, Nbar^0.5
    double sigma_;              // Standard deviation of random noise multiplied by N
    double *w_gpu_;             // GPU array containing: N*w-(r), N*w+(r), phi-(r), phi+(r)
    double *noise_gpu_;         // Array holding random noise for current step and previous step
    double *noise_gpu_new_;     // Pointer to portion of memory for new noise in noise_gpu_[]
    double *noise_gpu_prev_;    // Pointer to portion of memory for previous noise in noise_gpu_[]
    
    diblock  *dbc_;             // Diblock object for calculating phi-(r) and phi+(r)
    anderson *AM_;              // Anderson mixing object to solve for w+(r)
    metadynamics *WTMD_;
    curandGenerator_t RNG_;     // Random number generator for the GPU

    public:
        wtmd_simulation(std::string inputFile="", std::string biasFile="", int TpB=512) {
            
            // Check that input files exist before proceeding
            if (!(isValidFile(inputFile) && isValidFile(biasFile))) {
                std::cout << "ERROR => Cannot open one of the input files." << std::endl;
                exit(1);
            }

            TpB_ = TpB;
            L_ = new double[3];
            m_ = new int[3];

            // read input file member parameters and fields for copying to gpu
            double *w = readInputFile(inputFile);

            // Set up a metadynamics object
            WTMD_ = new metadynamics(m_, L_, C_*V_, Mk_, M_, biasFile);

            // Set up random number generator
            curandCreateGenerator(&RNG_, CURAND_RNG_PSEUDO_DEFAULT);
            int iseed = time(NULL);
            iseed = 123456789;
            std::cout << "RNG seed: " << iseed << std::endl;
            curandSetPseudoRandomGeneratorSeed(RNG_, iseed);

            // Allocate memory for Gaussian random noise on the GPU
            GPU_ERR(cudaMalloc((void**)&noise_gpu_,2*M_*sizeof(double)));
            
            // Generate initial "previous" Gaussian random noise on the gpu
            curandGenerateNormalDouble(RNG_, noise_gpu_, M_, 0.0, sigma_);
            noise_gpu_prev_ = noise_gpu_;
            noise_gpu_new_ = noise_gpu_ + M_;

            // Allocate memory for field array on the GPU copy over w-(r) and w+(r) from host
            GPU_ERR(cudaMalloc((void**)&w_gpu_,4*M_*sizeof(double)));
            GPU_ERR(cudaMemcpy(w_gpu_,w,2*M_*sizeof(double),cudaMemcpyHostToDevice));
            delete[] w;

            // Create a new diblock object
            dbc_ = new diblock(NA_, N_-NA_, M_, Mk_, m_, L_);

            // Create a new anderson mixing object
            AM_ = new anderson(M_);

            // Perform an initial mix to get phi-(r) and phi+(r) from the input fields
            AM_ -> mix(dbc_,200,1e-4,w_gpu_);
        }

        // Destructor
        ~wtmd_simulation() {
            delete[] L_;
            delete[] m_;
            GPU_ERR(cudaFree(w_gpu_));
            GPU_ERR(cudaFree(noise_gpu_));
            delete WTMD_;
            delete dbc_;
            delete AM_;
        }

        // Equilibration loop, during which statistics are NOT sampled
        void equilibrate() {
            int it;
            for (it=1; it<=equil_its_; it++) {

                // Perform a Langevin step with symmetrised noise to update w-(r). Last argument dictates whether inputted bias field is used.
                langevin(w_gpu_, chi_b_, sigma_, dt_, &noise_gpu_new_, &noise_gpu_prev_, WTMD_->get_read_bias());
                std::cout << "lnQ = " << dbc_->calc_concs(w_gpu_) << std::endl;
                
                // Save to file every save_freq_ steps
                if (it%save_freq_==0) { 
                    saveStdOutputFile("w_eq_" + std::to_string(it));
                    saveGPUArray("phi_eq_" + std::to_string(it) , w_gpu_+2*M_, 2*M_);
                }
            }

            // Final save to file at end of equilibration period
            saveStdOutputFile("w_eq_" + std::to_string(it-1));
            saveGPUArray("phi_eq_" + std::to_string(it-1) , w_gpu_+2*M_, 2*M_);
        }

        // Statistics loop, during which statistics are sampled
        void statistics() {
            int it;
            for (it=1; it<=sim_its_; it++) {

                // Perform a Langevin step with symmetrised noise to update w-(r)
                langevin(w_gpu_, chi_b_, sigma_, dt_, &noise_gpu_new_, &noise_gpu_prev_, true);
                std::cout << "lnQ = " << dbc_->calc_concs(w_gpu_) << std::endl;
                
                // Sample statistics every sample_freq_ steps
                if (it%sample_freq_==0) {
                }

                // Update the bias potential every WTMD_->get_update_freq() steps
                if (it%WTMD_->get_update_freq() == 0) {
                    WTMD_->update_bias_field(w_gpu_);
                }

                // Save fields to file every save_freq_ steps
                if (it%save_freq_==0) { 
                    saveStdOutputFile("w_st_" + std::to_string(it));
                    saveGPUArray("phi_st_" + std::to_string(it), w_gpu_+2*M_, 2*M_);
                    WTMD_->save_bias_fields("bias_st_" + std::to_string(it));
                }
            }

            // Final save to file at end of equilibration period
            saveStdOutputFile("w_st_" + std::to_string(it-1));
            saveGPUArray("phi_st_" + std::to_string(it-1), w_gpu_+2*M_, 2*M_);
            WTMD_->save_bias_fields("bias_st_" + std::to_string(it-1));
        }

        // Print parameters that were read from the input file to standard output
        void outputParameters() {
            std::cout << "N = "             << N_           << std::endl;
            std::cout << "NA = "            << NA_          << std::endl;
            std::cout << "chi_b = "         << chi_b_       << std::endl;
            std::cout << "C = "             << C_           << std::endl;
            std::cout << "dt = "            << dt_          << std::endl;
            std::cout << "m[0] = "          << m_[0]        << std::endl;
            std::cout << "m[1] = "          << m_[1]        << std::endl;
            std::cout << "m[2] = "          << m_[2]        << std::endl;
            std::cout << "L[0] = "          << L_[0]        << std::endl;
            std::cout << "L[1] = "          << L_[1]        << std::endl;
            std::cout << "L[2] = "          << L_[2]        << std::endl;
            std::cout << "equilit_s = "     << equil_its_   << std::endl;
            std::cout << "simit_s = "       << sim_its_     << std::endl;
            std::cout << "sample_freq = "   << sample_freq_ << std::endl;
        }

        // Calculate the diblock copolymer Hamiltonian
        double getH() {
            // Calculate the natural log of the partition function
            double lnQ = dbc_->calc_concs(w_gpu_);

            // Create a Thrust device pointer to the GPU memory for the fields
            thrust::device_ptr<double> dp = thrust::device_pointer_cast(w_gpu_);

            // Calculate the sum of w+(r) on the GPU
            double w_sum = thrust::reduce(dp+M_, dp+2*M_, 0.0, thrust::plus<double>());

            // Calculate the sum of w-(r)^2 on the GPU
            double w2_sum = thrust::transform_reduce(dp, dp+M_, thrust::square<double>(), 0.0, thrust::plus<double>());

            // Return the Hamiltonian
            return -lnQ + (w2_sum/chi_b_ - w_sum)/M_;
        }

    private:
        // Perform a Langevin update of the fields using symmetrised noise
        void langevin(double* w_gpu, double chi_b, double sigma, double dt, double **noise_gpu_new, double **noise_gpu_prev, bool useBias = true)
        {
            double *ptr_tmp;

            // Create new random noise on the GPU for the call to langevin()
            curandGenerateNormalDouble(RNG_, noise_gpu_new_, M_, 0.0, sigma_);

            // Perform the Langevin step on the GPU. The call depends on whether bias is being used.
            if (useBias) {
                WTMD_->langevin(w_gpu, chi_b, sigma, dt, noise_gpu_new, noise_gpu_prev);
            } else {
                langevin_sym<<<(M_+TpB_-1)/TpB_,TpB_>>>(w_gpu, *noise_gpu_new, *noise_gpu_prev, chi_b, dt, M_);
            }

            // Swap noise pointer positions to avoid shifting (copying) data every step
            ptr_tmp = *noise_gpu_prev;
            *noise_gpu_prev = *noise_gpu_new;
            *noise_gpu_new = ptr_tmp;

            // Perform anderson mixing to solve for w+(r) and obtain phi-(r) and phi+(r) (all stored in w_gpu_[])
            AM_->mix(dbc_,200,1e-4,w_gpu);
        }

        // Check whether a file exists
        bool isValidFile(std::string fileName) {
            if (fileName == "") return false;
            std::ifstream instream(fileName);
            return !(instream.fail());
        }

        // Read input parameters/fields and return a pointer to the w field on the host
        double *readInputFile(std::string fileName) {
            double *w;
            std::ifstream instream;
            instream.open(fileName);

            // Read the simulation parameters
            instream >> N_ >> NA_ >> chi_b_ >> C_ >> dt_;
            instream >> m_[0] >> m_[1] >> m_[2] >> L_[0] >> L_[1] >> L_[2];
            instream >> equil_its_ >> sim_its_ >> sample_freq_ >> save_freq_;

            // Calculate derived parameters
            sim_its_ = (sim_its_/sample_freq_)*sample_freq_;  
            M_ = m_[0]*m_[1]*m_[2];
            Mk_ = m_[0]*m_[1]*(m_[2]/2+1);
            V_ = L_[0]*L_[1]*L_[2];
            sigma_ = sqrt(2.0*M_*dt_/(C_*V_));

            // Read the w-(r) and w+(r) fields
            w = new double[2*M_];
            for (int r=0; r<2*M_; r++) {
                instream >> w[r];
            }
            instream.close();
            return w;
        }

        // Save the simulation parameters to file
        void saveOutputParams(std::string fileName, bool append=false) {
            std::ofstream outstream;
            if (append) outstream.open(fileName,std::ios_base::app);
            else outstream.open(fileName);
            outstream << N_ << " " << NA_ << " " << chi_b_ << " " << C_ << " " << dt_ << std::endl;
            outstream << m_[0] << " " << m_[1] << " " << m_[2] << " " << L_[0] << " " << L_[1] << " " << L_[2] << std::endl;
            outstream << equil_its_ << " " << sim_its_ << " " << sample_freq_ << " " << save_freq_ << std::endl;
            outstream.close();
        }

        // Save a GPU array to file
        void saveGPUArray(std::string fileName, double *arr_gpu, int n, bool append=false) {
            double *arr_cpu = new double[n];
            GPU_ERR(cudaMemcpy(arr_cpu,arr_gpu,n*sizeof(double),cudaMemcpyDeviceToHost));
            saveArray(fileName, arr_cpu, n, append);
            delete[] arr_cpu;
        }

        // Save a host array to file
        void saveArray(std::string fileName, double *arr, int n, bool append=false) {
            std::ofstream outstream;
            if (append) outstream.open(fileName,std::ios_base::app);
            else outstream.open(fileName);
            for (int r=0; r<n; r++) outstream << arr[r] << std::endl;
            outstream.close();
        }

        // Save data in a standard format to be used as in input file
        void saveStdOutputFile(std::string fileName) {
            saveOutputParams(fileName);
            double *w = new double[2*M_];
            GPU_ERR(cudaMemcpy(w,w_gpu_,2*M_*sizeof(double),cudaMemcpyDeviceToHost));
            saveArray(fileName, w, 2*M_, true);
            delete[] w;
        }
        
};
