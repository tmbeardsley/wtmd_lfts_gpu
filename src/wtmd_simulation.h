// ######################################################################################
// Exposes public methods to perform an L-FTS simulation: equilibrate() and statistics()
// ######################################################################################

#pragma once
#include <stdlib.h>
#include <string>
#include <cuda.h>
#include "GPUerror.h"
#include <iostream>
#include <fstream>
#include "diblock.h"
#include "anderson.h"
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <curand.h>
#include "field_generator.h"
#include "langevin.h"
#include "lfts_params.h"
#include "file_IO.h"
#include "metadynamics.h"


class wtmd_simulation {
    double *w_gpu_;             // GPU array containing: N*w-(r), N*w+(r), phi-(r), phi+(r)
    diblock  *dbc_;             // Diblock object for calculating phi-(r) and phi+(r)
    anderson *AM_;              // Anderson mixing object to solve for w+(r)
    langevin *Langevin_;        // Langevin object to update w-(r) at each step
    metadynamics *WTMD_;        // Metadynamics object to collect the biasing potentials for the langevin step
    curandGenerator_t RNG_;     // Random number generator for the GPU
    lfts_params *P_;            // Object to hold the simulation parameters - automatically updates derived parameters

    int M_;                     // Total number of field mesh points (constant - contained in lfts_params object but copied for tidier code)

    public:
        wtmd_simulation(std::string inputFile, std::string biasFile, int TpB=512) {

            // Check that input files exist before proceeding
            std::string s = "";
            if (!file_IO::isValidFile(inputFile))   s += "ERROR => Cannot open the L-FTS input file.\n";
            if (!file_IO::isValidFile(biasFile))    s += "ERROR => Cannot open the WTMD bias file.\n";
            if (s != "") {
                std::cout << s << std::endl;
                exit(1);
            }

            // Read simulation parameters from the input file and allocate temporary host memory for fields
            P_ = new lfts_params(inputFile);
            P_->outputParameters();
            M_=P_->M();
            double *w = new double[2*M_];

            // Set up random number generator
            curandCreateGenerator(&RNG_, CURAND_RNG_PSEUDO_DEFAULT);
            int iseed = time(NULL);
            iseed = 123456789;
            curandSetPseudoRandomGeneratorSeed(RNG_, iseed);
            std::cout << "RNG seed: " << iseed << std::endl;

            // Allocate memory for field array on the GPU
            GPU_ERR(cudaMalloc((void**)&w_gpu_, 4*M_*sizeof(double)));

            // Create a new diblock object
            std::cout << "creating diblock object..." << std::endl;
            dbc_ = new diblock(P_->NA(), P_->NB(), P_->m(), P_->L(), M_, P_->Mk(), TpB);

            // Create a new anderson mixing class
            std::cout << "creating anderson object..." << std::endl;
            AM_ = new anderson(M_, 10, TpB);

            // Set up a langevin object to upstate w-(r) at each step
            std::cout << "creating langevin object..." << std::endl;
            Langevin_ = new langevin(RNG_, P_->sigma(), M_, TpB);

            // Set up a metadynamics object
            std::cout << "creating metadynamics object..." << std::endl;
            WTMD_ = new metadynamics(P_->m(), P_->L(), P_->n(), P_->Mk(), M_, biasFile, TpB);

            // Read w-[r] and w+[r] from the input file
            if (P_->loadType() == 1) { 
                std::cout << "loading input field..." << std::endl;
                file_IO::readArray(w, inputFile, 2*M_, 3);
            }
            else generate_field(w, P_->loadType());

            // Copy w-(r) and w+(r) from host to GPU
            GPU_ERR(cudaMemcpy(w_gpu_,w,2*M_*sizeof(double),cudaMemcpyHostToDevice));
            delete[] w;

            // Perform an initial mix to get phi-(r) and phi+(r) from the input fields
            std::cout << "Initial Anderson mix..." << std::endl;
            AM_ -> mix(dbc_,200,1e-4,w_gpu_);

            // Output initial fields
            saveStdOutputFile("w_0");
            file_IO::saveGPUArray(w_gpu_+2*M_, "phi_0", 2*M_);
        }

        // Destructor
        ~wtmd_simulation() {
            GPU_ERR(cudaFree(w_gpu_));
            delete dbc_;
            delete AM_;
            delete Langevin_;
            delete WTMD_;
            delete P_;
        }

        // Equilibration loop, during which statistics are NOT sampled
        void equilibrate() {
            int it;
            for (it=1; it<=P_->equil_its(); it++) {

                // Perform a Langevin step with symmetrised noise to update w-(r)
                Langevin_->step_wm(w_gpu_, RNG_, P_->XbN(), P_->sigma(), P_->dt(), WTMD_->get_fBias(w_gpu_));

                // Calculate saddle point value of w+(r), phi-(r) and phi+(r)
                AM_->mix(dbc_,200,1e-4,w_gpu_);
                std::cout << "Psi = " << WTMD_->get_Psi(w_gpu_) << std::endl;

                // Save to file every save_freq_ steps
                if (it%P_->save_freq()==0) { 
                    saveStdOutputFile("w_eq_" + std::to_string(it));
                    file_IO::saveGPUArray(w_gpu_+2*M_, "phi_eq_"+std::to_string(it), 2*M_);
                }
            }
            // Final save to file at end of equilibration period
            saveStdOutputFile("w_eq_" + std::to_string(it-1));
            file_IO::saveGPUArray(w_gpu_+2*M_, "phi_eq_"+std::to_string(it-1), 2*M_);
        }

        // Statistics loop, during which statistics are sampled
        void statistics() {
            int it;
            for (it=1; it<=P_->sim_its(); it++) {

                // Perform a Langevin step with symmetrised noise to update w-(r)
                Langevin_->step_wm(w_gpu_, RNG_, P_->XbN(), P_->sigma(), P_->dt(), WTMD_->get_fBias(w_gpu_));

                // Calculate saddle point value of w+(r), phi-(r) and phi+(r)
                AM_->mix(dbc_,200,1e-4,w_gpu_);
                std::cout << "Psi = " << WTMD_->get_Psi(w_gpu_) << std::endl;

                // Sample statistics every sample_freq_ steps
                if (it%P_->sample_freq()==0) {
                }

                // Update the bias potential every WTMD_->get_update_freq() steps
                if (it%WTMD_->get_update_freq() == 0) {
                    WTMD_->update_bias_field(w_gpu_);
                }

                // Save fields to file every save_freq_ steps
                if (it%P_->save_freq()==0) { 
                    saveStdOutputFile("w_st_" + std::to_string(it));
                    file_IO::saveGPUArray(w_gpu_+2*M_, "phi_st_" + std::to_string(it), 2*M_);
                    WTMD_->save_bias_std_output("bias_st_" + std::to_string(it));
                }
            }
            // Final save to file at end of equilibration period
            saveStdOutputFile("w_st_" + std::to_string(it-1));
            file_IO::saveGPUArray(w_gpu_+2*M_, "phi_st_" + std::to_string(it-1), 2*M_);
            WTMD_->save_bias_std_output("bias_st_" + std::to_string(it-1));
        }

        

        // Calculate the diblock copolymer Hamiltonian
        double getH() {
            // Calculate the natural log of the partition function
            double lnQ = dbc_->calc_concs(w_gpu_);

            // Create a Thrust device pointer to the GPU memory for the fields
            thrust::device_ptr<double> dp = thrust::device_pointer_cast(w_gpu_);

            // Calculate the sum of w+(r) and the sum of w-(r)^2 on the GPU
            double w_sum = thrust::reduce(dp+M_, dp+2*M_, 0.0, thrust::plus<double>());
            double w2_sum = thrust::transform_reduce(dp, dp+M_, thrust::square<double>(), 0.0, thrust::plus<double>());

            // Return the Hamiltonian
            return -lnQ + (w2_sum/P_->XbN() - w_sum)/M_;
        }



    private:
        // Save data in a standard format to be used as in input file
        void saveStdOutputFile(std::string fileName) {
            P_->saveOutputParams(fileName);
            file_IO::saveGPUArray(w_gpu_, fileName, 2*M_, true);
        }

        // Generate an initial w-(r) field (sets w+(r) = 0)
        void generate_field(double *w, int loadType) {
            switch (loadType) {
                case 2:
                    field_generator::create_lamellar(w, P_->XbN(), P_->m());
                    std::cout << "Generated lamellar initial configuration..." << std::endl;
                    break;
                default:
                    // Create a random field with noise of amplitude XN/2 
                    // Using curand to keep a single, consistent RNG for reproducibilty
                    double *rand_gpu, *rand_cpu;
                    rand_cpu = new double[M_];
                    GPU_ERR(cudaMalloc((void**)&rand_gpu, M_*sizeof(double)));
                    curandGenerateUniformDouble(RNG_, rand_gpu, M_);
                    GPU_ERR(cudaMemcpy(rand_cpu,rand_gpu, M_*sizeof(double), cudaMemcpyDeviceToHost));
                    for (int r=0; r<M_; r++) {
                        w[r] = P_->XbN()*(rand_cpu[r]-0.5);
                        w[r+M_] = 0.0;
                    }
                    delete[] rand_cpu;
                    GPU_ERR(cudaFree(rand_gpu));
                    std::cout << "Generated disordered initial configuration..." << std::endl;
                    break;
            }
        }
        
};