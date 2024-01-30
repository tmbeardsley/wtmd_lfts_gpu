// ##########################
// Commonly used GPU kernels
// ##########################

#pragma once
#include <cuda.h>
#include "GPUerror.h"

// Set all elements of array, a, to the constant b
__global__ void Array_init(double *a, const double b, int const M)
{
	int const tid = threadIdx.x + blockIdx.x * blockDim.x;

	if (tid >= M) return;
	a[tid] = b;
}

// Set array, a, equal to array b. Add an optional constant, C, to each array element
__global__ void Array_copy(double *a, const double *b, int const M, double C = 0.0)
{
	int const tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= M) return;
	a[tid] = b[tid] + C;
}
