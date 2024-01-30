#!/bin/bash
#$ -cwd

nvcc -O2 -std=c++14 -L$EBROOTCUDA/lib64 -lcufft -lgsl -lgslcblas -lcurand -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80 fts_gpu.cu -o fts_gpu
