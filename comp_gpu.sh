#!/bin/bash
#$ -cwd

mkdir ./build

nvcc ./src/fts_gpu.cu -o ./build/wtmd-fts-gpu -O2 -std=c++14 -lcufft -lgsl -lgslcblas -lm -lcurand -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80
