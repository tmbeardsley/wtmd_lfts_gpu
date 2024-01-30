// ######################################################
// Provides public method: calc_concs(double *w_gpu), 
// to calculate concentrations (used in Anderson mixing)
// ######################################################

#pragma once
#include <cuda.h>
#include "GPUkernels.h"
#include "GPUerror.h"
#include "step.h"
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>


// Calculate hA[r] and hB[r]: h -> hA[0], h+M -> hB[0] 
__global__ void prepare_h(double *h, double *w_gpu, const int N, const int M)
{
	double *wm=w_gpu, *wp=w_gpu+M;
	int const tid = threadIdx.x + blockIdx.x * blockDim.x;

	if (tid >= 2*M) return;
	if (tid < M) h[tid] = exp(-(wm[tid]+wp[tid])/N);
	else h[tid] = exp(-(-wm[tid-M]+wp[tid-M])/N);
}

// Normalise concentrations phi-(r) and phi+(r): w_gpu+2*M -> phim[0],  w_gpu+3*M -> phip[0].
__global__ void normalize_phi(double *w_gpu, double *h, const double Q, const double N, const int M)
{
    double *phiA=h, *phiB=h+M, *hA=h, *hB=h+M;
    double *phim=w_gpu+2*M, *phip=w_gpu+3*M;
	int const tid = threadIdx.x + blockIdx.x * blockDim.x;

	if (tid >= M) return;
	phiA[tid] = phim[tid]/(N*Q*hA[tid]);
	phiB[tid] = phip[tid]/(N*Q*hB[tid]);
	phim[tid] = phiA[tid] - phiB[tid];
	phip[tid] = phiA[tid] + phiB[tid];
}

// Multiply and sum propagators for calculating either phiA[r] or phiB[r]
__global__ void sum_phi(double *phi_gpu, double *q1_gpu, double *q2_gpu, const int M)
{
	int const tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= M) return;
	phi_gpu[tid] += q1_gpu[tid]*q2_gpu[tid];
}







class diblock {

    // Diblock-specific variables
    int TpB_;                       // GPU threads per block (default: 512)
    int NA_;                        // Length of polymer A-block
    int N_;                         // Total polymer length
    int M_;                         // Total number of field mesh points
    double *qr_gpu_;                // Pointer to GPU memory for propagators: q_{i}(r) and q^_{N+1-i}(r) are contigious in memory
    double *h_gpu_;                 // Pointer to GPU memory for hA(r) and hB(r)
    double **q1_;                   // Array of pointers to q_{j=i}(r), where j is the monomer index and i is array index
    double **q2_;                   // Array of pointers to q^_{j=N+1-i}(r), where j is the monomer index and i is array index
    step *Step_;                    // Step object to get propagators for the next monomer

    public:
        // Constructor
        diblock(int NA, int _NB, int M, int _Mk, int *_m, double *_L, int TpB=512) {
            TpB_ = TpB;
            NA_ = NA;
            N_ = NA_ + _NB;
            M_ = M;

            // Allocate gpu memory for h_gpu_ and qr_gpu_
            GPU_ERR(cudaMalloc((void**)&h_gpu_,2*M_*sizeof(double)));
            GPU_ERR(cudaMalloc((void**)&qr_gpu_,2*(N_+1)*M_*sizeof(double)));

            // Allocate arrays of pointers for q_{j=1...N}(r) and q^_{j=1...N}(r)
            q1_ = new double* [N_+1];
            q2_ = new double* [N_+1];

            // Assign pointers such that q_{1}(r) and q_{N}(r) are in contigious memory,
            // as are q_{2}(r) and q_{N-1}(r), q_{3}(r) and q_{N-2}(r)... etc. (required for cufftPlanMany())
            for (int i=1; i<=N_; i++) {
                q1_[i] = qr_gpu_ + 2*i*M_;
                q2_[N_+1-i] = qr_gpu_ + (2*i+1)*M_;
            }

            // New step object containing methods to get next monomer's propagators
            Step_ = new step(_m, _L, NA_, _NB, _Mk, M_);
        }


        // Calculates phi-(r) and phi+(r): w+2*M -> phi-(0), w+3*M -> phi+(0).
        // Returns ln(Q)
        double calc_concs(double *w_gpu) {
            int i;
            double Q;                           // Single-chain partition function
            double *phiA_gpu=w_gpu+2*M_;
            double *phiB_gpu=w_gpu+3*M_;

            // Calculate hA[r] and hB[r] on the GPU
            prepare_h<<<(2*M_+TpB_-1)/TpB_, TpB_>>>(h_gpu_,w_gpu,N_,M_);

            // Set initial conditions: q[1][r]=hA[r] and q^[N][r]=hB[r] for all r
            Array_copy<<<(2*M_+TpB_-1)/TpB_, TpB_>>>(q1_[1],h_gpu_,2*M_);

            // Step the propagators q1 and q2 for each subsequent monomer (note q[i],q^[N+1-i]... contigious in memory)
            for (i=1; i<N_; i++) Step_->fwd(q1_[i], q1_[i+1], h_gpu_, i);

            // Calculate single-chain partition function using a Thrust reduction sum
            thrust::device_ptr<double> dp = thrust::device_pointer_cast(q1_[N_]);
            Q = thrust::reduce(dp, dp+M_, 0.0, thrust::plus<double>());
            Q /= M_;

            // Calculate concentrations using custom CUDA kernels
            Array_init<<<(2*M_+TpB_-1)/TpB_, TpB_>>>(phiA_gpu, 0.0, 2*M_);
            for (i=1; i<=NA_; i++) sum_phi<<<(M_+TpB_-1)/TpB_, TpB_>>>(phiA_gpu, q1_[i], q2_[i], M_);
            for (i=NA_+1; i<=N_; i++) sum_phi<<<(M_+TpB_-1)/TpB_, TpB_>>>(phiB_gpu, q1_[i], q2_[i], M_);
            normalize_phi<<<(M_+TpB_-1)/TpB_, TpB_>>>(w_gpu, h_gpu_, Q, N_, M_);

            // Return ln(Q)
            return log(Q);
        }

        // Destructor
        ~diblock() {
            GPU_ERR(cudaFree(h_gpu_));
            GPU_ERR(cudaFree(qr_gpu_));
            delete[] q1_;
            delete[] q2_;
            delete Step_;
        }
};