// #######################################################################################
// Provides public methods to update the bias potential and perform a modified (biased)
// Langevin step
// #######################################################################################

#pragma once
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cufft.h>
#include <fstream>
#include <limits>
#include <thrust/transform_reduce.h>
#include <thrust/complex.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/zip_function.h>
#include <iomanip>
#include "wtmd_params.h"

// GPU kernel to calculate (DPsi/dwk) on the GPU. Involves manipulation of complex numbers
__global__ void get_dPsi_dwk(cufftDoubleComplex *dPsi_dwk_gpu, cufftDoubleComplex *wk_gpu, double *fk_gpu, double const Psi, double const ell, int const Mk)
{
	int const tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= Mk) return;

    dPsi_dwk_gpu[tid].x = wk_gpu[tid].x*wk_gpu[tid].x + wk_gpu[tid].y*wk_gpu[tid].y;
	dPsi_dwk_gpu[tid].x = pow(dPsi_dwk_gpu[tid].x,ell/2.0-1.0);
	dPsi_dwk_gpu[tid].x = dPsi_dwk_gpu[tid].x * pow(Psi,1.0-ell);
	dPsi_dwk_gpu[tid].x = dPsi_dwk_gpu[tid].x * fk_gpu[tid];
	dPsi_dwk_gpu[tid].y = dPsi_dwk_gpu[tid].x;
	dPsi_dwk_gpu[tid].x *= wk_gpu[tid].x;
	dPsi_dwk_gpu[tid].y *= wk_gpu[tid].y;
}

// Multiply all elements of B[] by constant, C, and return in A[]
__global__ void aEqBc(double *A, double *B, const double c, const int M)
{
    int const tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= M) return;
    A[tid] = B[tid] * c;
}

// GPU kernel to update the bias fields, U(Psi), DU(Psi)/DPsi, I0(Psi) and I1(Psi)
__global__ void update_bias_fields_gpu(double *u_gpu, double *up_gpu, double *I0_gpu, double *I1_gpu, const double Psi_min, const double dPsi, const double Psi_hat, const double sigma_Psi, const double n, const double DT, const double w2_hat, const int mPsi)
{
    int const tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= mPsi) return;
    const double Psi = Psi_min + tid*dPsi;
    const double x = (Psi_hat-Psi)/sigma_Psi;
    const double A = exp(-n*u_gpu[tid]/DT)/n;
    const double G = exp(-0.5*x*x);
    u_gpu[tid]  += A*G;
    up_gpu[tid] += (x/sigma_Psi-n*up_gpu[tid]/DT)*A*G;
    I0_gpu[tid] += G;
    I1_gpu[tid] += w2_hat*G;
}

// Structure used in conjunction with thrust library functions to calculate the order parameter, Psi
struct Psi_calc
{
    const double ell_;
    const int M_;
    Psi_calc(double ell, int M):ell_(ell),M_(M) {}
	__device__ __host__
	double operator()(thrust::tuple<cufftDoubleComplex, double, int> t) {
        // wk = get<0>(t), fk = get<1>(t), wt = get<2>(t)
        cufftDoubleComplex wk = thrust::get<0>(t);
        thrust::complex<double> wkT;
        wkT.real(wk.x);
        wkT.imag(wk.y);
        double fk = thrust::get<1>(t);
        int wt = thrust::get<2>(t);
        return pow(thrust::norm(wkT),0.5*ell_)*fk*wt/M_;
	}
};




class metadynamics {
    int TpB_;

    // non-changing arrays (set in calc_wt_fk(...))
    int *wt_gpu_;           // Weighting of contribution from wavevector k
    double *fk_gpu_;        // Weighting function to reduce the contribution of large wavevectors in the order parameter

    // update_bias()
    typedef thrust::device_ptr<double> dpD_;

    // get_Psi()
    cufftHandle wr_to_wk_;
    typedef thrust::device_ptr<int> dpI_;
    typedef thrust::device_ptr<cufftDoubleComplex> dpC_;
    typedef thrust::tuple<dpC_,dpD_,dpI_> t_dpC_dpD_dpI_;
    typedef thrust::zip_iterator<t_dpC_dpD_dpI_> z_dpC_dpD_dpI_;
    dpC_ wk_tptr;           // Thrust device pointer to a complex array
    dpD_ fk_tptr;           // Thrust device pointer to a double array
    dpI_ wt_tptr;           // Thrust device pointer to an integer array
    t_dpC_dpD_dpI_ t_;      // Tuple used to set up zip iterator for calculation of order parameter
    z_dpC_dpD_dpI_ z_;      // Zip iterator used for thrust reductions in calculation of order parameter on the GPU

    // get_fBias()
    cufftHandle dPsi_wk_to_dPsi_wr_;    // Plan handle for transforming (DPsi/Dwk) -> (DPsi/Dwr) with cufft on the GPU

    // update_bias_field()
    double *u_gpu_;         // Bias potential, U(Psi)
    double *up_gpu_;        // Derivative of the bias potential, (DU(Psi)/DPsi)
    double *I0_gpu_;        // Array to which Gaussians are added for computing <w-(r)^2>
    double *I1_gpu_;        // Array to which Gaussians scaled by w-(r)^2 are added for computing <w-(r)^2>

    // langevin()
    double *fBias_gpu_;     // The biasing force used in the modified Langevin step

    // get_fBias()
    cufftDoubleComplex *dPsi_dwk_gpu_;      // (DPsi/Dwk) on the GPU
    double *dPsi_dwr_gpu_;                  // (DPsi/Dwr) on the GPU
    double *up_i_;                          // Adjacent values of (DU(Psi)/DPsi) for extrapolation

    // get_fBias(), get_Psi()
    cufftDoubleComplex *wk_gpu_;            // Composition field in reciprocal space, w-(k), on the GPU
    
    bool nonZeroBias_;                      // Indicates whether there is a non-zero bias field (for efficiency in Langevin step)
    wtmd_params *B_;                        // Object to deal with bias-related parameters read from file

    // Simulation constants derived from the input file (see lfts_params.h for details)
    int M_;
    int Mk_;
    double n_;


    public:
        metadynamics(int *m, double *L, double CV, int Mk, int M, std::string biasFile, int TpB=512) {
            TpB_ = TpB;
            Mk_ = Mk;
            M_ = M;
            n_ = CV;

            nonZeroBias_ = false;

            // Read first line of the input 
            B_ = new wtmd_params(biasFile);
            B_->printInputParams();

            // Allocate memory for lookup arrays, f(k) and wt(k) on the GPU
            GPU_ERR(cudaMalloc((void**)&fk_gpu_, Mk_*sizeof(double)));
            GPU_ERR(cudaMalloc((void**)&wt_gpu_, Mk_*sizeof(int)));
            calc_wt_fk(wt_gpu_, fk_gpu_, m, L);

            // langevin() - Allocate memory for the biasing force on the GPU and set to zero
            GPU_ERR(cudaMalloc((void**)&fBias_gpu_, M_*sizeof(double)));
            Array_init<<<(M+TpB_-1)/TpB_,TpB_>>>(fBias_gpu_, 0.0, M_);

            // get_fBias() - Set up the cufft plan and allocate GPU memory for derivatives
            GPU_ERR(cufftPlan3d(&dPsi_wk_to_dPsi_wr_, m[0], m[1], m[2], CUFFT_Z2D));
            GPU_ERR(cudaMalloc((void**)&dPsi_dwk_gpu_, Mk_*sizeof(cufftDoubleComplex)));
            GPU_ERR(cudaMalloc((void**)&dPsi_dwr_gpu_, M_*sizeof(double)));
            up_i_ = new double[2];

            // get_fBias(), get_Psi() - Allocate memory for w(k) on the GPU
            GPU_ERR(cudaMalloc((void**)&wk_gpu_, Mk_*sizeof(cufftDoubleComplex)));

            // get_Psi() - Set up cufft plan to transform w-(r) to reciprocal space. Set up thrust pointers and interators.
            GPU_ERR(cufftPlan3d(&wr_to_wk_, m[0], m[1], m[2], CUFFT_D2Z));
            wk_tptr = thrust::device_pointer_cast(wk_gpu_);
            fk_tptr = thrust::device_pointer_cast(fk_gpu_);
            wt_tptr = thrust::device_pointer_cast(wt_gpu_);
            t_ = thrust::make_tuple(wk_tptr, fk_tptr, wt_tptr);
            z_ = thrust::make_zip_iterator(t_);

            // update_bias() - Allocate memory and pointers for bias-related arrays and set to zero
            int mPsi = B_->mPsi();
            GPU_ERR(cudaMalloc((void**)&u_gpu_, 4*mPsi*sizeof(double)));
            Array_init<<<(4*mPsi+TpB_-1)/TpB_,TpB_>>>(u_gpu_, 0.0, 4*mPsi);
            up_gpu_ = u_gpu_ + mPsi;
            I0_gpu_ = u_gpu_ + 2*mPsi;
            I1_gpu_ = u_gpu_ + 3*mPsi;

            // Read the bias field if the flag was set in the input file
            if (B_->read_bias() != 0) read_Bias_Fields(biasFile, u_gpu_);
        }

        // Destructor
        ~metadynamics() {
            delete[] up_i_;
            GPU_ERR(cudaFree(fk_gpu_));
            GPU_ERR(cudaFree(wt_gpu_));
            GPU_ERR(cudaFree(fBias_gpu_));
            GPU_ERR(cudaFree(dPsi_dwk_gpu_));
            GPU_ERR(cudaFree(dPsi_dwr_gpu_));
            GPU_ERR(cudaFree(wk_gpu_));
            GPU_ERR(cudaFree(u_gpu_));
            GPU_ERR(cufftDestroy(dPsi_wk_to_dPsi_wr_));
            GPU_ERR(cufftDestroy(wr_to_wk_));
            delete B_;
        }

        // Calculate the order parameter 
        double get_Psi(double *w_gpu)
        {
            // Fourier transform w-(r) to get w-(k)
            GPU_ERR(cufftExecD2Z(wr_to_wk_, w_gpu, wk_gpu_));

            // Perform a thrust transform reduction on the gpu and calculate Psi
            double Psi = thrust::transform_reduce(z_, z_ + Mk_, Psi_calc(B_->ell(), M_), 0.0, thrust::plus<double>());
            Psi = pow(Psi/M_, 1.0/(B_->ell()));
            return Psi;
        }


        // Update u(Psi), up(Psi), I0(Psi) and I1(Psi) - function overloaded so psi doesn't have to be recalculated if already known
        void update_bias_field(double *w_gpu) { update_bias_field(get_Psi(w_gpu), w_gpu); }
        void update_bias_field(double Psi_hat, double *w_gpu)
        {
            // Bias field exists as the field has been updated
            nonZeroBias_ = true;

            // Perform a thrust transform reduction to calculate current value of w-^2 on the GPU
            dpD_ tPtr = thrust::device_pointer_cast(w_gpu);
            double w2_hat = thrust::transform_reduce(tPtr, tPtr + M_, thrust::square<double>(), 0.0, thrust::plus<double>());
            w2_hat /= M_;

            // Update the bias fields by calling a GPU kernel
            int mPsi = B_->mPsi();
            update_bias_fields_gpu<<<(mPsi+TpB_-1)/TpB_,TpB_>>>(u_gpu_, up_gpu_, I0_gpu_, I1_gpu_, B_->Psi_min(), B_->dPsi(), Psi_hat, B_->sigma_Psi(), n_, B_->DT(), w2_hat, mPsi);
        }

        // Make update_freq publicly accessible for wtmd_simulation.h
        int get_update_freq() {
            return B_->update_freq();
        }

        // Calculate the biasing force for the modified Langevin step.
        double* get_fBias(double *w_gpu)
        {
            // fBias_gpu_[]=0 if there is no bias, so don't waste resources computing it
            if (!nonZeroBias_) return fBias_gpu_;

            int    i;
            double Psi, x, up_hat;

            // Calculate current value of U'(Psi)
            Psi = get_Psi(w_gpu);
            x = (Psi - B_->Psi_min()) / B_->dPsi();
            i = floor(x);
            if (i < 0) {printf("Error: Psi = %lf < Psi_min\n",Psi); exit(1);}
            if (i >= B_->mPsi()) {printf("Error: Psi = %lf > Psi_max\n",Psi); exit(1); }
            x = x-i;

            // Linear interpolation of the (DU(Psi)/DPsi) due to discrete mesh
            GPU_ERR(cudaMemcpy(up_i_,up_gpu_+i,2*sizeof(double),cudaMemcpyDeviceToHost));
            up_hat = (1.0-x)*up_i_[0] + x*up_i_[1];

            // Calculate derivative of order parameter with respect to wk.
            // Note: wk_gpu was evaluated in the above call to get_Psi(wt, ell, fk)
            get_dPsi_dwk<<<(M_+TpB_-1)/TpB_, TpB_>>>(dPsi_dwk_gpu_, wk_gpu_, fk_gpu_, Psi, B_->ell(), Mk_);

            // Calculate derivative of order parameter with respect to w
            GPU_ERR(cufftExecZ2D(dPsi_wk_to_dPsi_wr_, dPsi_dwk_gpu_, dPsi_dwr_gpu_));

            // Multiply array elements by a constant on the GPU
            aEqBc<<<(M_+TpB_-1)/TpB_, TpB_>>>(fBias_gpu_, dPsi_dwr_gpu_, up_hat/M_, M_);

            return fBias_gpu_;
        }

        // Save standard bias output file
        void save_bias_std_output(std::string fileName) {
            B_->saveBiasParams(fileName);
            save_bias_fields(fileName, true);
        }

        // Save the bias parameters and fields to file in the same format as the bias input file
        void save_bias_fields(std::string fileName, bool append=false) 
        {
            std::ofstream outstream;
            if (append) outstream.open(fileName,std::ios_base::app);
            else outstream.open(fileName);
            outstream.precision(6);
            outstream << std::fixed;

            // Copy the bias fields from the GPU to the host
            int mPsi = B_->mPsi();
            double *biasFields = new double[4*mPsi];
            GPU_ERR(cudaMemcpy(biasFields,u_gpu_,4*mPsi*sizeof(double),cudaMemcpyDeviceToHost));
            double *u_tmp = biasFields, *up_tmp = biasFields + mPsi;
            double *I0_tmp = biasFields + 2*mPsi, *I1_tmp = biasFields + 3*mPsi;

            // Write the output fields to file
            for (int i=0; i<mPsi; i++) {
                outstream   << B_->Psi_min()+i*B_->dPsi()   << " " << std::scientific
                            << u_tmp[i]                     << " " 
                            << up_tmp[i]                    << " " 
                            << I0_tmp[i]                    << " " 
                            << I1_tmp[i]                    << " " << std::fixed
                            << std::endl;
            }
            outstream.close();
            delete[] biasFields;
        }

    private:
        // Read the bias fields from the input bias file
        void read_Bias_Fields(std::string fileName, double *u_gpu) {
            double Psi;
            int mPsi = B_->mPsi();
            double *u = new double[4*mPsi];
            double *up=u+mPsi, *I0=u+2*mPsi, *I1=u+3*mPsi;

            // Bias field exists as the field has been updated
            nonZeroBias_ = true;

            std::ifstream inFile;
            inFile.open(fileName);

            // Ignore parameters on first line (already read with read_Bias_Params(...))
            inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (int i=0; i<mPsi; i++) {
                inFile >> Psi >> u[i] >> up[i] >> I0[i] >> I1[i];
            }
            inFile.close();

            // Copy the input data to the GPU
            GPU_ERR(cudaMemcpy(u_gpu,u,4*mPsi*sizeof(double),cudaMemcpyHostToDevice));
            delete[] u;
        }

        // Calculate wt(k) and f(k) and copy them to the GPU
        void calc_wt_fk(int *wt_gpu, double *fk_gpu, int *m, double *L) {
            int K0, K1, k;
            double kx_sq, ky_sq, kz_sq, K;
            int *wt = new int[Mk_];
            double *fk = new double[Mk_];
            double kc = B_->kc();

            for (k=0; k<Mk_; k++) wt[k]=2;

            for (int k0=-(m[0]-1)/2; k0<=m[0]/2; k0++) {
                K0 = (k0<0)?(k0+m[0]):k0;
                kx_sq = k0*k0/(L[0]*L[0]);

                for (int k1=-(m[1]-1)/2; k1<=m[1]/2; k1++) {
                    K1 = (k1<0)?(k1+m[1]):k1;
                    ky_sq = k1*k1/(L[1]*L[1]);

                    for (int k2=0; k2<=m[2]/2; k2++) {
                        kz_sq = k2*k2/(L[2]*L[2]);
                        k = k2 + (m[2]/2+1)*(K1+m[1]*K0);
                        K = 2*M_PI*pow(kx_sq+ky_sq+kz_sq,0.5);
                        fk[k] = 1.0/(1.0 + exp(12.0*(K-kc)/kc));
                        if ((k2==0)||(k2==m[2]/2)) wt[k]=1;
                    }
                }
            }

            // Copy arrays from the host to the GPU
            GPU_ERR(cudaMemcpy(fk_gpu,fk,Mk_*sizeof(double),cudaMemcpyHostToDevice));
            GPU_ERR(cudaMemcpy(wt_gpu,wt,Mk_*sizeof(int),cudaMemcpyHostToDevice));
            delete[] wt;
            delete[] fk;
        }

};
