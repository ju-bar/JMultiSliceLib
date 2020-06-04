//
// CUDA C header file: ArrayOps.cuh
// declaration of array operations to be performed on the current GPU device
// (implementation, see ArrayOps.cu)
//
// Author: Juri Barthel (juribarthel@gmail.com)
// Date  : 2018-03-18
//
/*
This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <https://www.gnu.org/licenses/>
*/
//
#pragma once
//
#ifdef __INTELLISENSE__
//
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
//
#define __CUDACC__
//
#include <device_functions.h>
//
#endif
#include <cufft.h>

struct ArrayOpStats1 { // 1D kernel call stats
	unsigned int uSize; // data size
	int nBlockSize; // use block size
	int nGridSize; // use grid size
	ArrayOpStats1() : uSize(0), nBlockSize(0), nGridSize(0) {} // init constructor
};

// -----------------------------------------------------------------------------
// CUDA calculation kernels
// -----------------------------------------------------------------------------
//
// For clarity, use function and kernel names that are alike. Prefix "ArrayOp"
// is used for API functions. Kernels have a Suffix "Kernel".
// Maker interfaces as in the following example:
// ( *out_1..k , *in_1..l, size)
// k = # output pointers on device
// l = # input pointers on device
//
// Important: 
// - Any array pointer provided should be on device memory.
// - Never pass pointers to host memory, this will cause an illegal access.
// - Scalars may be passed from host memory as their value will be copied to the
//   device cache on call.
//
// -----------------------------------------------------------------------------

//__global__ void MetKernel(int *out_1, int *out_2, int * out_3, int * out_4, unsigned int size);
//cudaError_t ArrayOpMet(int *out_1, int *out_2, int * out_3, int * out_4, ArrayOpStats1 stats);

// sets in_1 as the real part of out_1: out_1[i].x = in_1[i]
__global__ void SetReKernel(cuComplex* out_1, float* in_1, unsigned int size);

// sets in_1 as the imaginary part of out_1: out_1[i].y = in_1[i]
__global__ void SetImKernel(cuComplex* out_1, float* in_1, unsigned int size);

// gets out_1 as the real part of in_1: out_1[i] = in_1[i].x
__global__ void GetReKernel(float* out_1, cuComplex* in_1, unsigned int size);

// gets out_1 as the imaginary part of in_1: out_1[i] = in_1[i].y
__global__ void GetImKernel(float *out_1, cuComplex* in_1, unsigned int size);


// calculates a linear combination of two complex arrays with offset: out_1[i] = in_1[i] * a + in_2[i] * b + c
__global__ void AddKernel(cuComplex *out_1, cuComplex *in_1, cuComplex* in_2, cuComplex a, cuComplex b, cuComplex c, unsigned int size);

// calculates the sum of two complex arrays with offset: out_1[i] = in_1[i] + in_2[i]
__global__ void AddKernel0(cuComplex *out_1, cuComplex *in_1, cuComplex* in_2, unsigned int size);

// calculates a linear combination of two float arrays with offset: out_1[i] = in_1[i] * a + in_2[i] * b + c
__global__ void FAddKernel(float *out_1, float *in_1, float* in_2, float a, float b, float c, unsigned int size);

// calculates an atomic add of two float arrays: out_1[i] = a * in_1[i] + in_2[i]
__global__ void FAdd1Kernel(float *out_1, float *in_1, float* in_2, float a, unsigned int size);

// calculates an atomic add of two float arrays: out_1[i] = in_1[i] + a * in_2[i]
__global__ void FAdd2Kernel(float *out_1, float *in_1, float* in_2, float a, unsigned int size);

// calculates an atomic add of two float arrays: out_1[i] = in_1[i] + in_2[i]
__global__ void FAdd0Kernel(float *out_1, float *in_1, float* in_2, unsigned int size);

// calculates complex rescale: out_1[i] = in_1[i] * sca
// - use this to rescale wave
__global__ void ScaKernel(cuComplex *out_1, cuComplex *in_1, float sca, unsigned int size);

// calculates real rescale: out_1[i] = in_1[i] * sca
// - use this to rescale intensities
__global__ void FScaKernel(float *out_1, float *in_1, float sca, unsigned int size);

// calculates complex*complex elemet-wise: out_1[i] = in_1[i] * in_2[i]
// - use this to apply complex transmission functions
__global__ void MulKernel(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, unsigned int size);

// calculates complex*complex elemet-wise: out_1[j1][i1] = in_1[j1][i1] * in_2[j2][i2]
// - use this to apply complex transmission functions smaller than the wavefunction
// - uses periodic boundary conditions on in_2
__global__ void MulSub2dKernel(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, unsigned int n0, unsigned int n1, unsigned int m0, unsigned int m1, unsigned int size);

// calculates complex*float elemet-wise: out_1[i] = in_1[i] * in_2[i]
// - use this to apply real transmission functions
__global__ void FMulKernel(cuComplex *out_1, cuComplex *in_1, float *in_2, unsigned int size);

// calculates float*float elemet-wise: out_1[i] = in_1[i] * in_2[i]
// - use this to multiply up real transmission functions
__global__ void FFMulKernel(float *out_1, float *in_1, float *in_2, unsigned int psize);

// calculates complex*complex array multiplication: out_1[i] = in_1[i] * in_2[i] * sca
// - use this to apply complex transmission functions
__global__ void MulScaKernel(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, float sca, unsigned int size);

// calculates complex*float array multiplication: out_1[i] = in_1[i] * in_2[i] * sca
// - use this to apply real transmission functions to wave
__global__ void FMulScaKernel(cuComplex *out_1, float *in_1, cuComplex *in_2, float sca, unsigned int size);

// calculates complex array arguments: out_1[i] = atan2( in_1[i].y, in_1[i].x )
__global__ void CArgKernel(float *out_1, cuComplex *in_1, unsigned int size);

// calculates complex array absolute square with scale: out_1[i] = sqrt( in_1[i]*conjg(in_1[i]) )
__global__ void CAbsKernel(float *out_1, cuComplex *in_1, unsigned int size);

// calculates complex array absolute square with scale: out_1[i] = sqrt( in_1[i]*conjg(in_1[i]) ) * sca
__global__ void CAbsScaKernel(float *out_1, cuComplex *in_1, float sca, unsigned int size);

// calculates complex array absolute square with scale: out_1[i] = in_1[i]*conjg(in_1[i])
// - use this to calculate probability distributions from wave
__global__ void CPowKernel(float *out_1, cuComplex *in_1, unsigned int size);

// calculates complex array absolute square with scale: out_1[i] = in_1[i]*conjg(in_1[i]) * sca
// - use this to calculate probability distributions from wave
__global__ void CPowScaKernel(float *out_1, cuComplex *in_1, float sca, unsigned int size);

// adds the complex array absolute square with scale to the output: out_1[i] += in_1[i]*conjg(in_1[i]) * sca
// - use this to calculate probability distributions from wave
__global__ void AddCPowScaKernel(float* out_1, cuComplex* in_1, float sca, unsigned int size);

// copies from in_1 to out_1 using a cyclic 2d shift of sh0 and sh1 positive along dimensions n0 and n1
// - out_1 -> shifted output data
// - in_1 -> input data
// - sh0, sh1 -> positive (right) shift values of dimension 0 and 1
// - n0, n1 -> grid size on dimensions 0 and 1
__global__ void CShift2dKernel(cuComplex *out_1, cuComplex *in_1, unsigned int sh0, unsigned int sh1, unsigned int n0, unsigned int n1);

// applies a shift and defocus offset to a wave-function: out_1[i] = in_[1]*Exp{ -I *2*Pi * [ dx*in_2[i] + dy*in_3[i] ] }
// - in_1 -> input wave function (Fourier space)
// - in_2 -> kx field [1/nm]
// - in_3 -> ky field [1/nm]
// - dx, dy = shifts x and y [nm]
__global__ void MulPhasePlate00Kernel(cuComplex *out_1, cuComplex *in_1, float *in_2, float *in_3, float dx, float dy, unsigned int size);

// applies a shift and defocus offset to a wave-function: out_1[i] = in_[1]*Exp{ -I *2*Pi * [ dx*in_2[i] + dy*in_3[i] + dz*(in_2[i]*in_2[i]+in_3[i]*in_3[i]) ] }
// - in_1 -> input wave function (Fourier space)
// - in_2 -> kx field [1/nm]
// - in_3 -> ky field [1/nm]
// - dx, dy = shifts x and y [nm]
// - dz = defocus * wavelength / 2 [nm^2]
__global__ void MulPhasePlate01Kernel(cuComplex *out_1, cuComplex *in_1, float *in_2, float *in_3, float dx, float dy, float dz, unsigned int size);

//// calculates the total sum of a float array using a reduction scheme: out_1[blockIdx.x] = SUM( in_1[i] , i=1 ... N-1)
//template <unsigned int blockSize> __global__ void FAddReduceKernel(float *out_1, float *in_1, unsigned int size);
//
//// calculates the total sum of a float array using a reduction scheme: out_1[blockIdx.x] = SUM( in_1[imask] , imask=1 ... NMASK-1)
//// but using an access mask defining which i is to be used of in_1.
//// It is assumed that the indices in mask are not causing overflow in in_1.
//template <unsigned int blockSize> __global__ void MaskFAddReduceKernel(float *out_1, int* mask, float *in_1, unsigned int size);
//
//// calculates the total dot product of two float arrays using a reduction scheme: out_1[blockIdx.x] = SUM( in_1[imask]*in_2[imask] , imask=1 ... NMASK-1)
//// but using an access mask defining which i is to be used of in_1 and in_2.
//// It is assumed that the indices in mask are not causing overflow in in_1 and in_2.
//template <unsigned int blockSize> __global__ void MaskFDotReduceKernel(float *out_1, int* mask, float *in_1, float *in_2, unsigned int size);



// -----------------------------------------------------------------------------
// API helper functions
// -----------------------------------------------------------------------------

// determines optimized threading size for parallel 1 array ops on current device with given array size
// - checks with MulScaKernel calls
cudaError_t GetOptimizedMultStats(unsigned int *size, int *blockSize, int *gridSize);

// -----------------------------------------------------------------------------
// Kernel wrapper functions
// -----------------------------------------------------------------------------

// sets the real part of out_1 from in_1: out_1[i].x = in_1[i] on device 
cudaError_t ArrayOpSetRe(cuComplex *out_1, float *in_1, ArrayOpStats1 stats);

// sets the imaginary part of out_1 from in_1: out_1[i].y = in_1[i] on device 
cudaError_t ArrayOpSetIm(cuComplex *out_1, float *in_1, ArrayOpStats1 stats);

// gets out_1 as the real part of in_1: out_1[i] = in_1[i].x on device 
cudaError_t ArrayOpGetRe(float *out_1, cuComplex *in_1, ArrayOpStats1 stats);

// gets out_1 as the imaginary part of in_1: out_1[i] = in_1[i].y on device 
cudaError_t ArrayOpGetIm(float *out_1, cuComplex *in_1, ArrayOpStats1 stats);

// calculates complex linear combination out_1[i] = in_1[i] * a + in_2[i] * b + c on device 
cudaError_t ArrayOpAdd(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, cuComplex a, cuComplex b, cuComplex c, ArrayOpStats1 stats);

// calculates complex sum out_1[i] = in_1[i] + in_2[i] on device 
cudaError_t ArrayOpAdd0(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, ArrayOpStats1 stats);

// calculates float linear combination out_1[i] = in_1[i] * a + in_2[i] * b + c on device 
cudaError_t ArrayOpFAdd(float *out_1, float *in_1, float *in_2, float a, float b, float c, ArrayOpStats1 stats);

// calculates float linear combination out_1[i] = in_1[i] + in_2[i] * a  on device 
cudaError_t ArrayOpFAdd2(float *out_1, float *in_1, float *in_2, float a, ArrayOpStats1 stats);

// calculates float linear combination out_1[i] = in_1[i] * a + in_2[i]  on device 
cudaError_t ArrayOpFAdd1(float *out_1, float *in_1, float *in_2, float a, ArrayOpStats1 stats);

// calculates float addition out_1[i] = in_1[i] + in_2[i] on device 
cudaError_t ArrayOpFAdd0(float *out_1, float *in_1, float *in_2, ArrayOpStats1 stats);

// calculates out_1[i] = in_1[i] * sca  on device 
cudaError_t ArrayOpSca(cuComplex *out_1, cuComplex *in_1, float sca, ArrayOpStats1 stats);

// calculates out_1[i] = in_1[i] * sca  on device 
cudaError_t ArrayOpFSca(float *out_1, float *in_1, float sca, ArrayOpStats1 stats);

// calculates out_1[i] = in_1[i] * in_2[i]  on device 
cudaError_t ArrayOpMul(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, ArrayOpStats1 stats);

// calculates out_1[j][i] = in_1[j][i] * in_2[j][i]  on device  where out_1 and in_1 are of the same dimension (n0,n1) and in_2 can be of different dimension (nsub0,nsub1), periodic wrap around is used
cudaError_t ArrayOpMulSub2d(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, unsigned int n0, unsigned int n1, unsigned int nsub0, unsigned int nsub1, ArrayOpStats1 stats);

// calculates out_1[i] = in_1[i] * in_2[i]  on device 
cudaError_t ArrayOpFMul(cuComplex *out_1, cuComplex *in_1, float *in_2, ArrayOpStats1 stats);

// calculates out_1[i] = in_1[i] * in_2[i]  on device 
cudaError_t ArrayOpFFMul(float *out_1, float *in_1, float *in_2, ArrayOpStats1 stats);


// calculates out_1[i] = in_1[i] * in_2[i] * sca  on device 
cudaError_t ArrayOpMulSca(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, float sca, ArrayOpStats1 stats);

// calculates out_1[i] = in_1[i] * in_2[i] * sca  on device 
cudaError_t ArrayOpFMulSca(cuComplex *out_1, float *in_1, cuComplex *in_2, float sca, ArrayOpStats1 stats);

// calculates out_1[i] = atan2( in_1[i].y, in_1[i].x )  on device 
cudaError_t ArrayOpCArg(float *out_1, cuComplex *in_1, ArrayOpStats1 stats);

// calculates out_1[i] = sqrt( in_1[i] * conjg(in_1[i]) )  on device 
cudaError_t ArrayOpCAbs(float *out_1, cuComplex *in_1, ArrayOpStats1 stats);

// calculates out_1[i] = sqrt( in_1[i] * conjg(in_1[i]) ) * sca  on device 
cudaError_t ArrayOpCAbsSca(float *out_1, cuComplex *in_1, float sca, ArrayOpStats1 stats);

// calculates out_1[i] = in_1[i] * conjg(in_1[i])  on device 
cudaError_t ArrayOpCPow(float *out_1, cuComplex *in_1, ArrayOpStats1 stats);

// calculates out_1[i] = in_1[i] * conjg(in_1[i]) * sca  on device 
cudaError_t ArrayOpCPowSca(float *out_1, cuComplex *in_1, float sca, ArrayOpStats1 stats);

// calculates out_1[i] += in_1[i] * conjg(in_1[i]) * sca  on device 
cudaError_t ArrayOpAddCPowSca(float* out_1, cuComplex* in_1, float sca, ArrayOpStats1 stats);

// calculates out_1[(j+sh1)%n1][(i+sh0)%n0] = in_1[j][i] on device
cudaError_t ArrayOpCShift2d(cuComplex *out_1, cuComplex *in_1, unsigned int sh0, unsigned int sh1, unsigned int n0, unsigned int n1, ArrayOpStats1 stats);


// multiplies a shift phase plate to a complex array
// - out_1 is the modified wave function
// - in_1 is the input wave function
// - in_2 is the kx-array [1/nm] (full field)
// - in_3 is the ky-array [1/nm] (full field)
// - dx and dy are the shift coordinates [nm]
cudaError_t ArrayOpMulPP00(cuComplex *out_1, cuComplex *in_1, float *in_2, float *in_3, float dx, float dy, ArrayOpStats1 stats);

// multiplies a shift and defocus phase plate to a complex array
// - out_1 is the modified wave function
// - in_1 is the input wave function
// - in_2 is the kx-array [1/nm] (full field)
// - in_3 is the ky-array [1/nm] (full field)
// - dx and dy are the shift coordinates [nm]
// - dz is the defocus multiplied by wave length / 2 [nm^2]
cudaError_t ArrayOpMulPP01(cuComplex *out_1, cuComplex *in_1, float *in_2, float *in_3, float dx, float dy, float dz, ArrayOpStats1 stats);



//// calculates the sum of a float array on device: out_1 = SUM( in_1[i] )
//// - CPU_threshold determines from which number of items on CPU summation will be carried out, set 0 to do all on GPU
//cudaError_t ArrayOpFSum(float &out_1, float * in_1, ArrayOpStats1 stats, int CPU_threshold = 0);
//
//// calculates the sum of a float array on device using a mask: out_1 = SUM( in_1[imask] )_imask[i]
//// - set stats.uSize to the size of the mask array
//// - set stats.nBlockSize to as much as many threads you want, it will be capped to 1024.
//// - CPU_threshold determines from which number of items on CPU summation will be carried out, set 0 to do all on GPU
//cudaError_t ArrayOpMaskFSum(float &out_1, int * mask, float * in_1, ArrayOpStats1 stats, int CPU_threshold = 0);

// calculates the dot product of two float arrays on device using a mask: out_1 = SUM( in_1[imask]*in_2[imask] )_imask[i]
// - set stats.uSize to the size of the mask array
// - set stats.nBlockSize to as much as many threads you want, it will be capped to 1024.
// - use this to calculate integrated detector results with mask and detector sensitivity (in_2) from a intensity distrib. (in_1)
// - CPU_threshold determines from which number of items on CPU summation will be carried out, set 0 to do all on GPU
cudaError_t ArrayOpMaskFDot(float &out_1, int * mask, float * in_1, float * in_2, ArrayOpStats1 stats, int CPU_threshold = 0);
// change additional masked dot product to eliminate cpowkernel
cudaError_t ArrayOpMaskFDot(float& out_1, int* mask, float2* in_1, float* in_2, ArrayOpStats1 stats, int CPU_threshold = 0);
cudaError_t ArrayOpFDot(float& out_1, float* in_1, float* in_2, ArrayOpStats1 stats, int CPU_threshold = 0);
cudaError_t ArrayOpFDot(float& out_1, float2* in_1, float* in_2, ArrayOpStats1 stats, int CPU_threshold = 0);
// end change

// calculates out_1 = SUM( in_1[i]*conjg( in_1[i] ) ) * sca on device
// - CPU_threshold determines from which number of items on CPU summation will be carried out, set 0 to do all on GPU
cudaError_t ArrayOpCPowSum(float &out_1, float2 *in_1, float sca, ArrayOpStats1 stats, int CPU_threshold = 0);

// calculates out_1 = SUM( in_1[i]*conjg( in_1[i] ) * in_2[i] ) * sca on device
// - CPU_threshold determines from which number of items on CPU summation will be carried out, set 0 to do all on GPU
cudaError_t ArrayOpCPowSumFMul(float &out_1, float2 *in_1, float *in_2, float sca, ArrayOpStats1 stats, int CPU_threshold = 0);

