//
// Created by Elmar Westphal on 25.03.20.
//

// proof of concept for replacing the masked dot product in DrProbe with something more modular
// not prettified yet, but it works :)

// Added more kernels 20.04.2020, Juri Barthel

#ifndef DRPROBEREDUCTION_REDUCTION_CUH
#define DRPROBEREDUCTION_REDUCTION_CUH

#include <cuda_runtime.h>
#include <algorithm>
#include <cassert>

// different loaders can be used to preprocess data
// in this one, complex numbers are read and processed to replace the call of the cpow kernel
template<typename T, typename T2>
struct maskedComplexDotProduct {
    __device__ __forceinline__ static T load(unsigned int idx, const int* const mask, const T2* const in_1, const T* const in_2) {
        int i = mask[idx];
        T2 t2 = in_1[i];
        T t = t2.x * t2.x + t2.y * t2.y; // pow on in_1
        return t * in_2[i]; // dot
    }
};

template<typename T>
struct maskedTDotProduct {
    __device__ __forceinline__ static T load(unsigned int idx, const int* const mask, const T* const in_1, const T* const in_2) {
        int i = mask[idx];
        return in_1[i] * in_2[i]; // dot
    }
};

template<typename T, typename T2>
struct ComplexDotProduct {
    __device__ __forceinline__ static T load(unsigned int idx, const T2* const in_1, const T* const in_2) {
        T2 t2 = in_1[idx];
        T t = t2.x * t2.x + t2.y * t2.y; // pow on in_1
        return t * in_2[idx]; // dot
    }
};

template<typename T>
struct TDotProduct {
    __device__ __forceinline__ static T load(unsigned int idx, const T* const in_1, const T* const in_2) {
        return in_1[idx] * in_2[idx]; // dot
    }
};

template<typename T, typename T2>
struct ComplexPower {
    __device__ __forceinline__ static T load(unsigned int idx, const T2* const in_1) {
        T2 t2 = in_1[idx];
        T t = t2.x * t2.x + t2.y * t2.y; // pow on in_1
        return t;
    }
};

template<typename T, typename T2>
struct ComplexPowerFMul {
    __device__ __forceinline__ static T load(unsigned int idx, const T2* const in_1, const T* const in_2) {
        T2 t2 = in_1[idx];
        T t = t2.x * t2.x + t2.y * t2.y; // pow on in_1
        return t * in_2[idx];
    }
};

// simple reduction kernel capable of using data loader classes as template arguments
// current compilers do all the loop unrolling etc. by default

template<typename Loader, int BlockSize, typename Out, typename ... In>
__global__ void addingReductionKernel(unsigned int n, Out *out, const In* const ... in) {
    Out t=0;
    for(unsigned int i=threadIdx.x+blockDim.x*blockIdx.x;i<n;i+=gridDim.x*blockDim.x) {
        t+=Loader::load(i,in...);
    }
    __shared__ Out c[BlockSize];
    int ti=threadIdx.x;
    c[ti]=t;
    for (int stage=BlockSize/2;stage;stage/=2) {
        __syncthreads();
        if (ti<stage)
            c[ti]+=c[ti+stage];
    }
    out[blockIdx.x]=c[0];
}

// wrapper class for reduction kernel, allocates output buffer in pinned host memory to avoid cudaMemcpy
// applies Kahan summation at the end

template<typename T>
class addingReduction {
    T* buffer;
    unsigned int max_blocks;
    static const int blockSize=256;

    public:

    addingReduction(int device=0) : buffer(nullptr), max_blocks(0) {
		cudaDeviceProp prop{}; // initialize to zero
        cudaError_t err = cudaGetDeviceProperties(&prop,device);
        max_blocks=prop.multiProcessorCount*(prop.maxThreadsPerMultiProcessor / blockSize);
        cudaMallocHost(&buffer,max_blocks*sizeof(T));
    }
    ~addingReduction() {
        if (buffer!= nullptr)
            cudaFreeHost(buffer);
    }
    template<typename Op, typename ... In>
    T perform(unsigned int n, In* ... in) {
        int gridSize=std::min(max_blocks,1+ (n-1) / blockSize);
        addingReductionKernel<Op,blockSize><<<gridSize, blockSize>>>(n,buffer,in...);
        cudaDeviceSynchronize();
        T dy,dc=0,dt,dsum=0;
        unsigned int nr = (unsigned int)gridSize;
        for (unsigned int ir = 0; ir < nr; ir++) {
            dy = buffer[ir] - dc; // next value including previous correction
            dt = dsum + dy; // intermediate new sum value
            dc = (dt - dsum) - dy; // new correction
            dsum = dt; // update result
        }
        return dsum;
    }


};

#endif //DRPROBEREDUCTION_REDUCTION_CUH
