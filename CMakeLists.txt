cmake_minimum_required(VERSION 3.18)
project(WTMD-LFTS-GPU LANGUAGES CXX CUDA)

# Set CUDA standard
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Specify the GPU architectures to compile for
set(CMAKE_CUDA_ARCHITECTURES all)

# Find required packages
find_package(CUDAToolkit REQUIRED)
find_package(GSL REQUIRED)

# Add the source file for the CUDA code
add_executable(wtmd-lfts-gpu ./src/wtmd_lfts_gpu.cu)

# Link required libraries (cuFFT, GSL, GSLcblas, cuRAND)
target_link_libraries(wtmd-lfts-gpu
    CUDA::cufft       # CUDA FFT library
    GSL::gsl          # GSL (GNU Scientific Library)
    CUDA::curand      # CUDA random number generation library
)