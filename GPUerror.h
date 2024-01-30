// #####################################################################
// Functions to manage problems arising from CUDA and the cuFFT library
// #####################################################################

#pragma once
#include <cuda.h>
#include <cufft.h>
#include <stdio.h>


void HandleError(cudaError_t err, char const * const file, int const line)
{
    if (err != cudaSuccess) {
        fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err), file, line);
        exit(EXIT_FAILURE);
    }
}

void HandleError(cufftResult err, char const * const file, int const line)
{
    if (err != CUFFT_SUCCESS) {
        switch(err) {
            case CUFFT_INVALID_PLAN:
                printf ("cufft %s in %s at line %d\n", "CUFFT_INVALID_PLAN", file, line);
                break;
            case CUFFT_INVALID_VALUE:
                printf ("cufft %s in %s at line %d\n", "CUFFT_INVALID_VALUE", file, line);
                break;
            case CUFFT_INTERNAL_ERROR:
                printf ("cufft %s in %s at line %d\n", "CUFFT_INTERNAL_ERROR", file, line);
                break;
            case CUFFT_EXEC_FAILED:
                printf ("cufft %s in %s at line %d\n", "CUFFT_EXEC_FAILED", file, line);
                break;
            case CUFFT_SETUP_FAILED:
                printf ("cufft %s in %s at line %d\n", "CUFFT_EXEC_FAILEDCUFFT_SETUP_FAILED", file, line);
                break;
            default:
                printf ("cufft %s in %s at line %d\n", "CUFFT_EXEC_FAILEDCUFFT_SETUP_FAILED", file, line);
        }
        exit(EXIT_FAILURE);
    }
}

#define GPU_ERR( err ) (HandleError( err, __FILE__, __LINE__ ))

