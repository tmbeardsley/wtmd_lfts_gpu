// #######################################################################################
// Provides the public method: void fwd(...), which takes the propagators of the previous
// monomer as input and returns the propagators of the next monomer as output
// #######################################################################################

#pragma once
#include <cuda.h>
#include <cufft.h>
#include "GPUerror.h"
#include <math.h>

// Element by element multiplication of complex array, a[], by double array, b[]
__global__ void Mult_self(cufftDoubleComplex *a, const double *b, int const M)
{
	int const tid = threadIdx.x + blockIdx.x * blockDim.x;

	if (tid >= M) return;
	a[tid].x *= b[tid];
	a[tid].y *= b[tid];
}

// Element by element multiplication of double array, a[], by double array, b[]
__global__ void Mult_self(double *a, const double *b, int const M)
{
	int const tid = threadIdx.x + blockIdx.x * blockDim.x;

	if (tid >= M) return;
	a[tid] *= b[tid];
}




class step {
    // Step-specific variables
    int TpB_;                           // GPU threads per block (default: 512)
    int M_;                             // Total number of field mesh points
    int Mk_;                            // Total number of field mesh points in k-space
    int NA_;                            // Length of polymer A-block
    int NB_;                            // Length of polymer B-block
    int *m_;                            // Number of mesh points [mx,my,mz]

    double *g_gpu_;                     // Bond potential Boltzmann weight, Fourier transformed and /M_ on the GPU
    cufftDoubleComplex *qk_gpu_;        // Fourier transforms of q1 and q2 on the GPU (for cufftPlanMany())
    cufftHandle qr_to_qk_;              // cufft plan to transform q1[r] and q2[r] to k-space
    cufftHandle qk_to_qr_;              // cufft plan to transform q1[k] and q2[k] to real-space

    public:
        // Constructor
        step(int *m, double *_L, int NA, int NB, int Mk, int M, int TpB=512) {
            TpB_ = TpB;
            M_ = M;
            Mk_ = Mk;
            NA_ = NA;
            NB_ = NB;
            m_ = new int[3];
            for (int i=0; i<3; i++) m_[i] = m[i];

            // Allocate memory for g on cpu and calculate it
            double *g = new double[Mk_];
            calc_g(g, _L);

            // Allocate memory for g_gpu and copy g from the host. 
            // g_gpu_ contains two copies of g[] so that q1[k] and q2[k] can be multiplied on the GPU at the same time
            GPU_ERR(cudaMalloc((void**)&g_gpu_,2*Mk_*sizeof(double)));
            GPU_ERR(cudaMemcpy(g_gpu_,g,Mk_*sizeof(double),cudaMemcpyHostToDevice));
            GPU_ERR(cudaMemcpy(g_gpu_+Mk_,g,Mk_*sizeof(double),cudaMemcpyHostToDevice));
            delete[] g;

            // Allocate memory for q1[k] and q2[k], stored in contigious memory
            GPU_ERR(cudaMalloc((void**)&qk_gpu_,2*Mk_*sizeof(cufftDoubleComplex)));

            // Configure cufft plans. cufftPlanMany used for batched processing
            GPU_ERR(cufftPlanMany(&qr_to_qk_,3,m_,NULL,1,0,NULL,1,0,CUFFT_D2Z,2));
            GPU_ERR(cufftPlanMany(&qk_to_qr_,3,m_,NULL,1,0,NULL,1,0,CUFFT_Z2D,2));
        }

        // Calculate propagators for the next monomer given propagators of previous monomer
        // q_in  = q{i}(r), q^{N+1-i}(r)
        // q_out = q{i+1}(r), q^{N-i}(r)
        // h_gpu = hA_gpu, hB_gpu
        void fwd(double* q_in, double* q_out, double *h_gpu, int i)
        {
            // Fourier transform q{i}(r),q^{N+1-i}(r) to k-space: qk_gpu_(k)
            GPU_ERR(cufftExecD2Z(qr_to_qk_, q_in, qk_gpu_));

            // Multiply qk_gpu_(k) by g_gpu_(k)
            Mult_self <<<(2*Mk_+TpB_-1)/TpB_, TpB_>>> (qk_gpu_, g_gpu_, 2*Mk_);

            // Fourier transform qk_gpu_(k) to real-space: q_out[r]
            GPU_ERR(cufftExecZ2D(qk_to_qr_, qk_gpu_, q_out));

            // Multiply q_out[r] by hA[r] or hB[r] (depending on monomer index) to get q{i+1}(r)
            if (i < NA_) Mult_self<<<(M_+TpB_-1)/TpB_,TpB_>>>(q_out,h_gpu,M_);
            else Mult_self<<<(M_+TpB_-1)/TpB_,TpB_>>>(q_out,h_gpu+M_,M_);

            // Multiply q_out[r+M_] by hA[r] or hB[r] (depending on monomer index) to get q^{N-i}(r)
            if (i < NB_) Mult_self<<<(M_+TpB_-1)/TpB_,TpB_>>>(q_out+M_,h_gpu+M_,M_);
            else Mult_self<<<(M_+TpB_-1)/TpB_,TpB_>>>(q_out+M_,h_gpu,M_);
        }

        // Destructor
        ~step() {
            GPU_ERR(cufftDestroy(qr_to_qk_));
            GPU_ERR(cufftDestroy(qk_to_qr_));
            GPU_ERR(cudaFree(g_gpu_));
            GPU_ERR(cudaFree(qk_gpu_));
            delete[] m_;
        }


    private:
        // Calculate the Boltzmann weight of the bond potential in k-space, _g[k]
        void calc_g(double *g, double *_L) {
            int K0, K1, k;
            double K, kx_sq, ky_sq, kz_sq;

            for (int k0=-(m_[0]-1)/2; k0<=m_[0]/2; k0++) {
                K0 = (k0<0)?(k0+m_[0]):k0;
                kx_sq = k0*k0/(_L[0]*_L[0]);

                for (int k1=-(m_[1]-1)/2; k1<=m_[1]/2; k1++) {
                    K1 = (k1<0)?(k1+m_[1]):k1;
                    ky_sq = k1*k1/(_L[1]*_L[1]);

                    for (int k2=0; k2<=m_[2]/2; k2++) {
                        kz_sq = k2*k2/(_L[2]*_L[2]);
                        k = k2 + (m_[2]/2+1)*(K1+m_[1]*K0);
                        K = 2*M_PI*pow(kx_sq+ky_sq+kz_sq,0.5); 
                        g[k] = exp(-K*K/(6*(NA_+NB_)))/M_; 
                    }
                }
            }
        }

};




