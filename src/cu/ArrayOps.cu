// file ArrayOps.cu
// author: Juri Barthel, ju.barthel@fz-juelich.de
// 2018
//
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
// Implementation for CUDA array operations
//

#include <stdio.h>
#include <cooperative_groups.h>
#include "ArrayOps.cuh"

namespace cg = cooperative_groups;

//__global__ void MetKernel(int *out_1, int *out_2, int * out_3, int * out_4, unsigned int size)
//{
//	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
//	if (idx < size)
//	{
//		out_1[idx] = (int)idx;
//		out_2[idx] = (int)threadIdx.x;
//		out_3[idx] = (int)blockIdx.x;
//		out_4[idx] = (int)blockDim.x;
//	}
//}
//cudaError_t ArrayOpMet(int *out_1, int *out_2, int * out_3, int * out_4, ArrayOpStats1 stats)
//{
//	cudaError_t cudaStatus;
//	unsigned int size = stats.uSize;	// input size
//	int blockSize = stats.nBlockSize;	// block size 
//	int gridSize = stats.nGridSize;		// grid size needed, based on input size
//
//	// Launch the parallel kernel operation
//	MetKernel<<<gridSize, blockSize>>>(out_1, out_2, out_3, out_4, size);
//
//	// Check for any errors launching the kernel
//	cudaStatus = cudaGetLastError();
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "- MetKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
//		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
//		goto Error;
//	}
//
//	// synchronize threads and wait for all to be finished
//	cudaStatus = cudaDeviceSynchronize();
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaDeviceSynchronize failed after launching MetKernel: %s\n", cudaGetErrorString(cudaStatus));
//		goto Error;
//	}
//
//Error:
//	return cudaStatus;
//}


// sets in_1 as the real part of out_1: out_1[i].x = in_1[i]
__global__ void SetReKernel(cuComplex* out_1, float* in_1, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx].x = in_1[idx];
	}
	return;
}

// sets in_1 as the imaginary part of out_1: out_1[i].y = in_1[i]
__global__ void SetImKernel(cuComplex* out_1, float* in_1, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx].y = in_1[idx];
	}
	return;
}

// gets out_1 as the real part of in_1: out_1[i] = in_1[i].x
__global__ void GetReKernel(float* out_1, cuComplex* in_1, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = in_1[idx].x;
	}
	return;
}

// gets out_1 as the imaginary part of in_1: out_1[i] = in_1[i].y
__global__ void GetImKernel(float *out_1, cuComplex* in_1, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = in_1[idx].y;
	}
	return;
}


// calculates a linear combination of two complex arrays with offset: out_1[i] = in_1[i] * a + in_2[i] * b + c
__global__ void AddKernel(cuComplex *out_1, cuComplex *in_1, cuComplex* in_2, cuComplex a, cuComplex b, cuComplex c, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		cuComplex c1a = cuCmulf(in_1[idx], a);
		cuComplex c2b = cuCmulf(in_2[idx], b);
		out_1[idx].x = c1a.x + c2b.x + c.x;
		out_1[idx].y = c1a.y + c2b.y + c.y;
	}
	return;
}

// calculates the sum of two complex arrays: out_1[i] = in_1[i] + in_2[i]
__global__ void AddKernel0(cuComplex *out_1, cuComplex *in_1, cuComplex* in_2, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = cuCaddf(in_1[idx],in_2[idx]);
	}
	return;
}

// calculates a linear combination of two float arrays with offset: out_1[i] = in_1[i] * a + in_2[i] * b + c
__global__ void FAddKernel(float *out_1, float *in_1, float* in_2, float a, float b, float c, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = in_1[idx] * a + in_2[idx] * b + c;
	}
	return;
}

// calculates an atomic add of two float arrays: out_1[i] = in_1[i] + in_2[i]
__global__ void FAdd0Kernel(float *out_1, float *in_1, float* in_2, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = in_1[idx] + in_2[idx];
	}
	return;
}

// calculates an atomic add of two float arrays: out_1[i] = in_1[i] * a + in_2[i]
__global__ void FAdd1Kernel(float *out_1, float *in_1, float* in_2, float a, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = in_1[idx] * a + in_2[idx];
	}
	return;
}

// calculates an atomic add of two float arrays: out_1[i] = in_1[i] + in_2[i] * a
__global__ void FAdd2Kernel(float *out_1, float *in_1, float* in_2, float a, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = in_1[idx] + in_2[idx] * a;
	}
	return;
}

// calculates complex rescale: out_1[i] = in_1[i] * sca
// - use this to rescale wave
__global__ void ScaKernel(cuComplex *out_1, cuComplex *in_1, float sca, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx].x = in_1[idx].x*sca;
		out_1[idx].y = in_1[idx].y*sca;
	}
	return;
}

// calculates real rescale: out_1[i] = in_1[i] * sca
// - use this to rescale intensities
__global__ void FScaKernel(float *out_1, float *in_1, float sca, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = in_1[idx] * sca;
	}
}

// calculates complex*complex element-wise: out_1[i] = in_1[i] * in_2[i]
// - use this to apply complex transmission functions
__global__ void MulKernel(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		float x = in_1[idx].x * in_2[idx].x - in_1[idx].y * in_2[idx].y;
		float y = in_1[idx].x * in_2[idx].y + in_1[idx].y * in_2[idx].x;
		//cuComplex ctemp = cuCmulf(in_1[idx], in_2[idx]);
		//out_1[idx] = ctemp;
		out_1[idx].x = x;
		out_1[idx].y = y;
	}
}

// calculates complex*complex element-wise: out_1[j1][i1] = in_1[j1][i1] * in_2[j2][i2]
// - use this to apply complex transmission functions smaller than the wavefunction
// - uses periodic boundary conditions on in_2 on periods (m0,m1)
__global__ void MulSub2dKernel(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, unsigned int n0, unsigned int n1, unsigned int m0, unsigned int m1, unsigned int size)
{
	unsigned int idx1 = threadIdx.x + blockIdx.x * blockDim.x; // pixel index in data stream 1
	if (idx1 < size) {
		unsigned int i1 = idx1 % n0; // i1 = column index in frame 1 (n0,n1) - from thread id + block set
		unsigned int i2 = i1 % m0; // i2 = column index in frame 2 (m0,m1)
		unsigned int j1 = (idx1 - i1) / n0; // j1 = row index in frame 1 (n0,n1) - from thread id + block set
		unsigned int j2 = j1 % m1; // j2 = row index in frame 2 (m0,m1)
		unsigned int idx2 = i2 + j2*m0; // pixel index in data stream 2
		float x = in_1[idx1].x * in_2[idx2].x - in_1[idx1].y * in_2[idx2].y; // complex product, real part
		float y = in_1[idx1].x * in_2[idx2].y + in_1[idx1].y * in_2[idx2].x; // cpmplex product, imaginary part
		//cuComplex ctemp = cuCmulf(in_1[idx], in_2[idx]);
		//out_1[idx] = ctemp;
		out_1[idx1].x = x;
		out_1[idx1].y = y;
	}
}

// calculates complex*float element-wise: out_1[i] = in_1[i] * in_2[i]
// - use this to apply real transmission functions
__global__ void FMulKernel(cuComplex *out_1, cuComplex *in_1, float *in_2, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx].x = in_1[idx].x * in_2[idx];
		out_1[idx].y = in_1[idx].y * in_2[idx];
	}
}

// calculates float*float element-wise: out_1[i] = in_1[i] * in_2[i]
// - use this to multiply up real transmission functions
__global__ void FFMulKernel(float *out_1, float *in_1, float *in_2, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = in_1[idx] * in_2[idx];
	}
}

// calculates complex*complex array multiplication: out_1[i] = in_1[i] * in_2[i] * sca
// - use this to apply complex transmission functions
__global__ void MulScaKernel(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, float sca, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size){
		cuComplex ctemp = cuCmulf(in_1[idx], in_2[idx]);
		out_1[idx].x = sca * ctemp.x;
		out_1[idx].y = sca * ctemp.y;
	}
}

// calculates complex*float array multiplication: out_1[i] = in_1[i] * in_2[i] * sca
// - use this to apply real transmission functions to wave
__global__ void FMulScaKernel(cuComplex *out_1, float *in_1, cuComplex *in_2, float sca, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx].x = sca * in_1[idx] * in_2[idx].x;
		out_1[idx].y = sca * in_1[idx] * in_2[idx].y;
	}
}

// calculates complex array arguments: out_1[i] = atan2( in_1[i].y, in_1[i].x )
__global__ void CArgKernel(float *out_1, cuComplex *in_1, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = atan2f( in_1[idx].y , in_1[idx].x );
	}
}


// calculates complex array absolute square with scale: out_1[i] = sqrt( in_1[i]*conjg(in_1[i]) )
__global__ void CAbsKernel(float *out_1, cuComplex *in_1, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = sqrtf( in_1[idx].x * in_1[idx].x + in_1[idx].y * in_1[idx].y );
	}
}

// calculates complex array absolute square with scale: out_1[i] = sqrt( in_1[i]*conjg(in_1[i]) ) * sca
__global__ void CAbsScaKernel(float *out_1, cuComplex *in_1, float sca, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = sca * sqrtf(in_1[idx].x * in_1[idx].x + in_1[idx].y * in_1[idx].y);
	}
}


// calculates complex array absolute square with scale: out_1[i] = in_1[i]*conjg(in_1[i])
// - use this to calculate probability distributions from wave
__global__ void CPowKernel(float *out_1, cuComplex *in_1, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = in_1[idx].x * in_1[idx].x + in_1[idx].y * in_1[idx].y;
	}
}

// calculates complex array absolute square with scale: out_1[i] = in_1[i]*conjg(in_1[i]) * sca
// - use this to calculate probability distributions from wave
__global__ void CPowScaKernel(float *out_1, cuComplex *in_1, float sca, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < size) {
		out_1[idx] = sca * (in_1[idx].x * in_1[idx].x + in_1[idx].y * in_1[idx].y);
	}
}

// copies from in_1 to out_1 using a cyclic 2d shift of sh0 and sh1 positive along dimensions n0 and n1
__global__ void CShift2dKernel(cuComplex *out_1, cuComplex *in_1, unsigned int sh0, unsigned int sh1, unsigned int n0, unsigned int n1) {
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x; // source index from thread id
	unsigned int size = n0 * n1;
	unsigned int i0 = idx % n0;
	unsigned int j0 = (idx - i0) / n0;
	unsigned int i1 = (i0 + sh0) % n0;
	unsigned int j1 = (j0 + sh1) % n1;
	unsigned int idx1 = i1 + n0 * j1;
	if (idx < size && idx1 < size) {
		out_1[idx1] = in_1[idx];
	}
}

// applies a shift offset to a wave-function: out_1[i] = in_[1]*Exp{ -I *2*Pi * [ dx*in_2[i] + dy*in_3[i]) ] }
// - in_1 -> input wave function (Fourier space)
// - in_2 -> kx field [1/nm]
// - in_3 -> ky field [1/nm]
// - dx, dy = shifts x and y [nm]
__global__ void MulPhasePlate00Kernel(cuComplex *out_1, cuComplex *in_1, float *in_2, float *in_3, float dx, float dy, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuComplex tf;
	if (idx < size) {
		float chi = -6.28318531f*(dx * in_2[idx] + dy * in_3[idx]);
		tf.x = cosf(chi);
		tf.y = sinf(chi);
		out_1[idx] = cuCmulf(in_1[idx], tf);
	}
}

// applies a shift and defocus offset to a wave-function: out_1[i] = in_[1]*Exp{ -I *2*Pi * [ dx*in_2[i] + dy*in_3[i] + dz*(in_2[i]*in_2[i]+in_3[i]*in_3[i]) ] }
// - in_1 -> input wave function (Fourier space)
// - in_2 -> kx field [1/nm]
// - in_3 -> ky field [1/nm]
// - dx, dy = shifts x and y [nm]
// - dz = defocus * wavelength [nm^2]
__global__ void MulPhasePlate01Kernel(cuComplex *out_1, cuComplex *in_1, float *in_2, float *in_3, float dx, float dy, float dz, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuComplex tf;
	if (idx < size) {
		float ksqr = in_2[idx] * in_2[idx] + in_3[idx] * in_3[idx];
		float chi = -6.28318531f*(dx*in_2[idx] + dy*in_3[idx] + dz*ksqr);
		tf.x = cosf(chi);
		tf.y = sinf(chi);
		out_1[idx] = cuCmulf(in_1[idx] , tf);
	}
}

//
//template <unsigned int blockSize>
//__global__ void	reduce6(int * g_idata, int * g_odata, unsigned int n)
//{
//	extern __shared__ int sdata[];
//	unsigned int tid = threadIdx.x;
//	unsigned int i = blockIdx.x*(blockSize * 2) + tid;
//	unsigned int gridSize = blockSize * 2 *gridDim.x;
//	sdata[tid] = 0;
//	while(i < n) {
//		sdata[tid] += g_idata[i] + g_idata[i + blockSize];
//		i += gridSize;
//	}
//	__syncthreads();
//	if (blockSize >= 512) {
//		if (tid < 256) {
//			sdata[tid] += sdata[tid + 256];
//		} 
//		__syncthreads();
//	}
//	if (blockSize >= 256) {
//		if (tid < 128) {
//			sdata[tid] += sdata[tid + 128];
//		} 
//		__syncthreads();
//	}
//	if (blockSize >= 128) {
//		if (tid <   64) {
//			sdata[tid] += sdata[tid + 64];
//		} __syncthreads();
//	}
//	if (tid < 32)
//	{
//		if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
//		if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
//		if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
//		if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
//		if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
//		if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
//	}
//	if (tid == 0) g_odata[blockIdx.x] = sdata[0];
//}
//
// reduction data access by threads 0-N, blocks 0-M, and grids 0-K:
// | | | (  0,  0,0) ... (  0,  0,N-1) & (  0,  0,N) ... (  0,  0,2N-1) | 
//     | (  0,  1,0) ... (  0,  1,N-1) & (  0,  1,N) ... (  0,  1,2N-1) |
//     ...
//     | (  0,M-1,0) ... (  0,M-1,N-1) & (  0,M-1,N) ... (  0,M-1,2N-1) | |
//   | | (  1,  0,0) ... (  1,  0,N-1) & (  1,  0,N) ... (  1,  0,2N-1) |
//     | (  1,  1,0) ... (  1,  1,N-1) & (  1,  1,N) ... (  1,  1,2N-1) |
//     ...
//     | (  1,M-1,0) ... (  1,M-1,N-1) & (  1,M-1,N) ... (  1,M-1,2N-1) | |
//   ...
//     | (K-1,M-1,0) ... (K-1,M-1,N-1) & (K-1,M-1,N) ... (K-1,M-1,2N-1) | | |

// calculates the total sum of a float array using a reduction scheme: out_1[blockIdx.x] = SUM( in_1[i] , i=1 ... N-1)
// launch this with half the number of blocks (gridSize) since we load strided by 2 * blockSize
// setup:
//      threads = (n < maxThreads * 2) ? nextPow2((n + 1) / 2) : maxThreads;
//      blocks  = (n + (threads * 2 - 1)) / (threads * 2);
//      smemSize = (threads <= 32) ? 2 * threads * sizeof(float) : threads * sizeof(float);
// call FAddReduceKernel<threads><<<blocks,threads,smemSize>>>(out,in,n);
template <unsigned int blockSize> __global__ void FAddReduceKernel(float *out_1, float *in_1, unsigned int size)
{
	// declare and initialize
	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	extern __shared__ float sdata[]; // ? unclear how the size of this is set, it should be blockSize
	unsigned int tid = threadIdx.x; // thread index
	unsigned int i = blockIdx.x*(blockSize * 2) + tid; // 2 x stride index 1 (initial)
	unsigned int j = i + blockSize; // 2x stride index 2 (initial)
	unsigned int gridSize = blockSize * 2 *gridDim.x; // offset to the next grid
	float mySum = 0.f;
	sdata[tid] = 0.f; // init cache for each thread
	// load initial data and first level of reduction
	// with this loop each of the blockSize threads in the current block will access
	// - two elements i and j per grid add them up
	// - then the same elements i and j of the next grid will be added also
	// - the added results will end up in shared memory sdata[tid] for each thread %tid%
	while(i < size) { // make sure not to read out of array bounds
		mySum += in_1[i];
		if (j < size) { // ensure to avoid out of bound for strided access
			mySum += in_1[j];
		}
		i += gridSize; // stride to next grid
		j = i + blockSize; // address to second half of the block
	}
	sdata[tid] = mySum;
	cg::sync(cta); // finish loading for all threads
	// reduce on shared memory
	// Note: more higher levels could be inserted if future compute classes allow more threads per block
	// unroll cascade for #thread >=512
	if (blockSize >= 1024) { // adds second 256 items to first 256 using threads 0 - 255
		if (tid < 512) {
			sdata[tid] = mySum = mySum + sdata[tid + 512];
		}
		cg::sync(cta);
	}
	// unroll cascade for #thread >=512
	if (blockSize >= 512) { // adds second 256 items to first 256 using threads 0 - 255
		if (tid < 256) {
			sdata[tid] = mySum = mySum + sdata[tid + 256];
		} 
		cg::sync(cta);
	}
	// unroll cascade for #threads >= 256
	if (blockSize >= 256) { // adds second 128 items to first 128 using threads 0 - 127
		if (tid < 128) {
			sdata[tid] = mySum = mySum + sdata[tid + 128];
		} 
		cg::sync(cta);
	}
	// unroll cascade for #threads >= 128
	if (blockSize >= 128) { // adds second 64 items to first 64 using threads 0 - 63
		if (tid <   64) {
			sdata[tid] = mySum = mySum + sdata[tid + 64];
		}
		cg::sync(cta);
	}
	//
	//// simpler code that also works, but without warp
	//// unroll SIMD synchroneous reductions within a warp
	//if ( tid < 32 )
	//{
	//	if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
	//	__syncthreads();
	//	if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
	//	__syncthreads();
	//	if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
	//	__syncthreads();
	//	if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
	//	__syncthreads();
	//	if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
	//	__syncthreads();
	//	if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
	//	__syncthreads();
	//}
	//
#if (__CUDA_ARCH__ >= 300 )
	if (tid < 32)
	{
		cg::coalesced_group active = cg::coalesced_threads();

		// Fetch final intermediate sum from 2nd warp
		if (blockSize >= 64) mySum += sdata[tid + 32];
		// Reduce final warp using shuffle
		for (int offset = warpSize / 2; offset > 0; offset /= 2)
		{
			mySum += active.shfl_down(mySum, offset);
		}
	}
#else
	// fully unroll reduction within a single warp
	if ((blockSize >= 64) && (tid < 32))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 32];
	}

	cg::sync(cta);

	if ((blockSize >= 32) && (tid < 16))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 16];
	}

	cg::sync(cta);

	if ((blockSize >= 16) && (tid <  8))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 8];
	}

	cg::sync(cta);

	if ((blockSize >= 8) && (tid <  4))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 4];
	}

	cg::sync(cta);

	if ((blockSize >= 4) && (tid <  2))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 2];
	}

	cg::sync(cta);

	if ((blockSize >= 2) && (tid <  1))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 1];
	}

	cg::sync(cta);
#endif
	// this block is reduced to one value, write it to output
	if (tid == 0) out_1[blockIdx.x] = mySum;
}


// calculates the total sum of a float array using a reduction scheme: out_1[blockIdx.x] = SUM( in_1[imask] , imask=1 ... NMASK-1)
// but using an access mask defining which i is to be used of in_1.
// It is assumed that the indices in mask are not causing overflow in in_1.
// launch this with half the number of blocks (gridSize) since we load strided by 2 * blockSize
// setup:
//      threads = (nmask < maxThreads * 2) ? nextPow2((nmask + 1) / 2) : maxThreads;
//      blocks  = (n + (threads * 2 - 1)) / (threads * 2);
//      smemSize = (threads <= 32) ? 2 * threads * sizeof(float) : threads * sizeof(float);
// call MaskFAddReduceKernel<threads><<<blocks,threads,smemSize>>>(out,mask,in,nmask);
template <unsigned int blockSize> __global__ void MaskFAddReduceKernel(float *out_1, int* mask, float *in_1, unsigned int size)
{
	// declare and initialize
	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	extern __shared__ float sdata[]; // ? unclear how the size of this is set, it should be blockSize
	unsigned int tid = threadIdx.x; // thread index
	unsigned int i = blockIdx.x*(blockSize * 2) + tid; // 2 x stride index 1 (initial)
	unsigned int j = i + blockSize; // 2x stride index 2 (initial)
	unsigned int imask = 0;
	unsigned int gridSize = blockSize * 2 * gridDim.x; // offset to the next grid
	float mySum = 0.f;
	sdata[tid] = 0.f; // init cache for each thread
					  // load initial data and first level of reduction
					  // with this loop each of the blockSize threads in the current block will access
					  // - two elements i and j per grid add them up
					  // - then the same elements i and j of the next grid will be added also
					  // - the added results will end up in shared memory sdata[tid] for each thread %tid%
	while (i < size) { // make sure not to read out of array bounds
		imask = (unsigned)mask[i]; // get the first array index from mask
		mySum += in_1[imask]; // add
		if (j < size) { // ensure to avoid out of bound for strided access
			imask = (unsigned)mask[j]; // get the second array index from mask
			mySum += in_1[imask]; // add
		}
		i += gridSize; // stride to next grid
		j = i + blockSize; // address to second half of the block
	}
	sdata[tid] = mySum; // store in shared memory
	cg::sync(cta); // finish loading for all threads
				   // reduce on shared memory
				   // Note: more higher levels could be inserted if future compute classes allow more threads per block
				   // unroll cascade for #thread >=512
	if (blockSize >= 1024) { // adds second 256 items to first 256 using threads 0 - 255
		if (tid < 512) {
			sdata[tid] = mySum = mySum + sdata[tid + 512];
		}
		cg::sync(cta);
	}
	// unroll cascade for #thread >=512
	if (blockSize >= 512) { // adds second 256 items to first 256 using threads 0 - 255
		if (tid < 256) {
			sdata[tid] = mySum = mySum + sdata[tid + 256];
		}
		cg::sync(cta);
	}
	// unroll cascade for #threads >= 256
	if (blockSize >= 256) { // adds second 128 items to first 128 using threads 0 - 127
		if (tid < 128) {
			sdata[tid] = mySum = mySum + sdata[tid + 128];
		}
		cg::sync(cta);
	}
	// unroll cascade for #threads >= 128
	if (blockSize >= 128) { // adds second 64 items to first 64 using threads 0 - 63
		if (tid <   64) {
			sdata[tid] = mySum = mySum + sdata[tid + 64];
		}
		cg::sync(cta);
	}
	//
	//// simpler code that also works, but without warp
	//// unroll SIMD synchroneous reductions within a warp
	//if ( tid < 32 )
	//{
	//	if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
	//	__syncthreads();
	//	if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
	//	__syncthreads();
	//	if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
	//	__syncthreads();
	//	if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
	//	__syncthreads();
	//	if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
	//	__syncthreads();
	//	if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
	//	__syncthreads();
	//}
	//
#if (__CUDA_ARCH__ >= 300 )
	if (tid < 32)
	{
		cg::coalesced_group active = cg::coalesced_threads();

		// Fetch final intermediate sum from 2nd warp
		if (blockSize >= 64) mySum += sdata[tid + 32];
		// Reduce final warp using shuffle
		for (int offset = warpSize / 2; offset > 0; offset /= 2)
		{
			mySum += active.shfl_down(mySum, offset);
		}
	}
#else
	// fully unroll reduction within a single warp
	if ((blockSize >= 64) && (tid < 32))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 32];
	}

	cg::sync(cta);

	if ((blockSize >= 32) && (tid < 16))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 16];
	}

	cg::sync(cta);

	if ((blockSize >= 16) && (tid <  8))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 8];
	}

	cg::sync(cta);

	if ((blockSize >= 8) && (tid <  4))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 4];
	}

	cg::sync(cta);

	if ((blockSize >= 4) && (tid <  2))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 2];
	}

	cg::sync(cta);

	if ((blockSize >= 2) && (tid <  1))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 1];
	}

	cg::sync(cta);
#endif
	// this block is reduced to one value, write it to output
	if (tid == 0) out_1[blockIdx.x] = mySum;
}


// calculates the total dot product of two float arrays using a reduction scheme: out_1[blockIdx.x] = SUM( in_1[imask]*in_2[imask] , imask=1 ... NMASK-1)
// but using an access mask defining which i is to be used of in_1 and in_2.
// It is assumed that the indices in mask are not causing overflow in in_1 and in_2.
// launch this with half the number of blocks (gridSize) since we load strided by 2 * blockSize
// setup:
//      threads = (nmask < maxThreads * 2) ? nextPow2((nmask + 1) / 2) : maxThreads;
//      blocks  = (n + (threads * 2 - 1)) / (threads * 2);
//      smemSize = (threads <= 32) ? 2 * threads * sizeof(float) : threads * sizeof(float);
// call MaskFDotReduceKernel<threads><<<blocks,threads,smemSize>>>(out,mask,in_1,in_2,nmask);
template <unsigned int blockSize> __global__ void MaskFDotReduceKernel(float *out_1, int* mask, float *in_1, float *in_2, unsigned int size)
{
	// declare and initialize
	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	extern __shared__ float sdata[]; // ? unclear how the size of this is set, it should be blockSize
	unsigned int tid = threadIdx.x; // thread index
	unsigned int i = blockIdx.x*(blockSize * 2) + tid; // 2 x stride index 1 (initial)
	unsigned int j = i + blockSize; // 2x stride index 2 (initial)
	unsigned int imask = 0;
	unsigned int gridSize = blockSize * 2 * gridDim.x; // offset to the next grid
	float mySum = 0.f;
	sdata[tid] = 0.f; // init cache for each thread
					  // load initial data and first level of reduction
					  // with this loop each of the blockSize threads in the current block will access
					  // - two elements i and j per grid add them up
					  // - then the same elements i and j of the next grid will be added also
					  // - the added results will end up in shared memory sdata[tid] for each thread %tid%
	while (i < size) { // make sure not to read out of array bounds
		imask = (unsigned)mask[i]; // get the first array index from mask
		mySum += in_1[imask]*in_2[imask]; // add projected component
		if (j < size) { // ensure to avoid out of bound for strided access
			imask = (unsigned)mask[j]; // get the second array index from mask
			mySum += in_1[imask]* in_2[imask]; // add projected component
		}
		i += gridSize; // stride to next grid
		j = i + blockSize; // address to second half of the block
	}
	sdata[tid] = mySum; // store in shared memory
	cg::sync(cta); // finish loading for all threads
				   // reduce on shared memory
				   // Note: more higher levels could be inserted if future compute classes allow more threads per block
				   // unroll cascade for #thread >=512
	if (blockSize >= 1024) { // adds second 256 items to first 256 using threads 0 - 255
		if (tid < 512) {
			sdata[tid] = mySum = mySum + sdata[tid + 512];
		}
		cg::sync(cta);
	}
	// unroll cascade for #thread >=512
	if (blockSize >= 512) { // adds second 256 items to first 256 using threads 0 - 255
		if (tid < 256) {
			sdata[tid] = mySum = mySum + sdata[tid + 256];
		}
		cg::sync(cta);
	}
	// unroll cascade for #threads >= 256
	if (blockSize >= 256) { // adds second 128 items to first 128 using threads 0 - 127
		if (tid < 128) {
			sdata[tid] = mySum = mySum + sdata[tid + 128];
		}
		cg::sync(cta);
	}
	// unroll cascade for #threads >= 128
	if (blockSize >= 128) { // adds second 64 items to first 64 using threads 0 - 63
		if (tid <   64) {
			sdata[tid] = mySum = mySum + sdata[tid + 64];
		}
		cg::sync(cta);
	}
	//
	//// simpler code that also works, but without warp
	//// unroll SIMD synchroneous reductions within a warp
	//if ( tid < 32 )
	//{
	//	if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
	//	__syncthreads();
	//	if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
	//	__syncthreads();
	//	if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
	//	__syncthreads();
	//	if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
	//	__syncthreads();
	//	if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
	//	__syncthreads();
	//	if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
	//	__syncthreads();
	//}
	//
#if (__CUDA_ARCH__ >= 300 )
	if (tid < 32)
	{
		cg::coalesced_group active = cg::coalesced_threads();

		// Fetch final intermediate sum from 2nd warp
		if (blockSize >= 64) mySum += sdata[tid + 32];
		// Reduce final warp using shuffle
		for (int offset = warpSize / 2; offset > 0; offset /= 2)
		{
			mySum += active.shfl_down(mySum, offset);
		}
	}
#else
	// fully unroll reduction within a single warp
	if ((blockSize >= 64) && (tid < 32))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 32];
	}

	cg::sync(cta);

	if ((blockSize >= 32) && (tid < 16))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 16];
	}

	cg::sync(cta);

	if ((blockSize >= 16) && (tid <  8))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 8];
	}

	cg::sync(cta);

	if ((blockSize >= 8) && (tid <  4))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 4];
	}

	cg::sync(cta);

	if ((blockSize >= 4) && (tid <  2))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 2];
	}

	cg::sync(cta);

	if ((blockSize >= 2) && (tid <  1))
	{
		sdata[tid] = mySum = mySum + sdata[tid + 1];
	}

	cg::sync(cta);
#endif
	// this block is reduced to one value, write it to output
	if (tid == 0) out_1[blockIdx.x] = mySum;
}


// Helper API to determine optimum multiplication call scenario
cudaError_t GetOptimizedMultStats(unsigned int *size, int *blockSize, int *gridSize)
{
	cudaError_t cudaStatus = cudaSuccess;
	int minGridSize = 0;    // The minimum grid size needed to achieve the maximum occupancy for a full device launch 

	cudaStatus = cudaOccupancyMaxPotentialBlockSize(&minGridSize, blockSize, MulScaKernel);
	if (cudaSuccess != cudaStatus) {
		fprintf(stderr, "Failed to determine optimum block size, code: %i\n", cudaStatus);
		goto Error;
	}
	cudaStatus = cudaDeviceSynchronize();
	if (cudaSuccess != cudaStatus) {
		fprintf(stderr, "Failed to snychronize device, code: %i\n", cudaStatus);
		goto Error;
	}
	// Round up according to array size 
	*gridSize = (*size + *blockSize - 1) / *blockSize;
Error:
	return cudaStatus;
}


// sets the real part of out_1 from in_1: out_1[i].x = in_1[i] on device 
cudaError_t ArrayOpSetRe(cuComplex *out_1, float *in_1, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	SetReKernel<<<gridSize, blockSize>>>(out_1, in_1, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- SetReKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching SetReKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// sets the imaginary part of out_1 from in_1: out_1[i].y = in_1[i] on device 
cudaError_t ArrayOpSetIm(cuComplex *out_1, float *in_1, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	SetImKernel<<<gridSize, blockSize>>>(out_1, in_1, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- SetImKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching SetImKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// gets out_1 as the real part of in_1: out_1[i] = in_1[i].x on device 
cudaError_t ArrayOpGetRe(float *out_1, cuComplex *in_1, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	GetReKernel<<<gridSize, blockSize>>>(out_1, in_1, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- GetReKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching GetReKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// gets out_1 as the imaginary part of in_1: out_1[i] = in_1[i].y on device 
cudaError_t ArrayOpGetIm(float *out_1, cuComplex *in_1, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	GetImKernel<<<gridSize, blockSize>>>(out_1, in_1, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- GetImKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching GetImKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}


// calculates complex linear combination out_1[i] = in_1[i] * a + in_2[i] * b + c on device 
cudaError_t ArrayOpAdd(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, cuComplex a, cuComplex b, cuComplex c, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	AddKernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, a, b, c, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- AddKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching AddKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates complex sum out_1[i] = in_1[i] + in_2[i] on device 
cudaError_t ArrayOpAdd0(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	AddKernel0<<<gridSize, blockSize>>>(out_1, in_1, in_2, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- AddKernel0 launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching AddKernel0: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates float linear combination out_1[i] = in_1[i] * a + in_2[i] * b + c on device 
cudaError_t ArrayOpFAdd(float *out_1, float *in_1, float *in_2, float a, float b, float c, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	FAddKernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, a, b, c, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- FAddKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching FAddKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates float linear combination out_1[i] = in_1[i] * a + in_2[i] on device 
cudaError_t ArrayOpFAdd1(float *out_1, float *in_1, float *in_2, float a, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	FAdd1Kernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, a, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- FAdd1Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching FAdd1Kernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates float linear combination out_1[i] = in_1[i] + in_2[i] * a on device 
cudaError_t ArrayOpFAdd2(float *out_1, float *in_1, float *in_2, float a, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	FAdd2Kernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, a, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- FAdd2Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching FAdd2Kernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates float addition out_1[i] = in_1[i] + in_2[i] on device 
cudaError_t ArrayOpFAdd0(float *out_1, float *in_1, float *in_2, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

										// Launch the parallel kernel operation
	FAdd0Kernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- FAdd0Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching FAddKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}


// calculates out_1[i] = in_1[i] * scale  on device 
cudaError_t ArrayOpSca(cuComplex *out_1, cuComplex *in_1, float sca, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	ScaKernel<<<gridSize, blockSize>>>(out_1, in_1, sca, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "ScaKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching ScaKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates out_1[i] = in_1[i] * scale  on device 
cudaError_t ArrayOpFSca(float *out_1, float *in_1, float sca, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	FScaKernel<<<gridSize, blockSize>>>(out_1, in_1, sca, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "FScaKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching FScaKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates out_1[i] = in_1[i] * in_2[i]  on device 
cudaError_t ArrayOpMul(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size
	
	// Launch the parallel kernel operation
	MulKernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- MulKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching MulKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}


// calculates out_1[j][i] = in_1[j][i] * in_2[j][i]  on device  where out_1 and in_1 are of the same dimension (n0,n1) and in_2 can be of different dimension (nsub0,nsub1), periodic wrap around is used
cudaError_t ArrayOpMulSub2d(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, unsigned int n0, unsigned int n1, unsigned int nsub0, unsigned int nsub1, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	MulSub2dKernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, n0, n1, nsub0, nsub1, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "- MulSub2dKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching MulSub2dKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates out_1[i] = in_1[i] * in_2[i]  on device 
cudaError_t ArrayOpFMul(cuComplex *out_1, cuComplex *in_1, float *in_2, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	FMulKernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "FMulKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching FMulKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates out_1[i] = in_1[i] * in_2[i]  on device 
cudaError_t ArrayOpFFMul(float *out_1, float *in_1, float *in_2, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	FFMulKernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "FFMulKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching FFMulKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates out_1[i] = in_1[i] * in_2[i] * sca  on device 
cudaError_t ArrayOpMulSca(cuComplex *out_1, cuComplex *in_1, cuComplex *in_2, float sca, ArrayOpStats1 stats)
{
    cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	MulScaKernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, sca, size);
    
	// Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "MulScaKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
        goto Error;
    }

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching NulScaKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
    
Error:
    return cudaStatus;
}

// calculates out_1[i] = in_1[i] * in_2[i] * sca  on device 
cudaError_t ArrayOpFMulSca(cuComplex *out_1, float *in_1, cuComplex *in_2, float sca, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	FMulScaKernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, sca, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "FMulScaKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching FMulScaKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates out_1[i] = atan2( in_1[i].y, in_1[i].x )  on device 
cudaError_t ArrayOpCArg(float *out_1, cuComplex *in_1, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	CArgKernel<<<gridSize, blockSize>>>(out_1, in_1, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "CArgKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching CArgKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates out_1[i] = sqrt( in_1[i] * conjg(in_1[i]) ) on device 
cudaError_t ArrayOpCAbs(float *out_1, cuComplex *in_1, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	CAbsKernel<<<gridSize, blockSize>>>(out_1, in_1, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "CAbsKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching CAbsKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates out_1[i] = sqrt( in_1[i] * conjg(in_1[i]) ) * sca  on device 
cudaError_t ArrayOpCAbsSca(float *out_1, cuComplex *in_1, float sca, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	CAbsScaKernel<<<gridSize, blockSize>>>(out_1, in_1, sca, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "CAbsScaKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching CAbsScaKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}


// calculates out_1[i] = in_1[i] * conjg(in_1[i]) on device 
cudaError_t ArrayOpCPow(float *out_1, cuComplex *in_1, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	CPowKernel<<<gridSize, blockSize>>>(out_1, in_1, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "CPowKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching CPowKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculates out_1[i] = in_1[i] * conjg(in_1[i]) * sca  on device 
cudaError_t ArrayOpCPowSca(float *out_1, cuComplex *in_1, float sca, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	CPowScaKernel<<<gridSize, blockSize>>>(out_1, in_1, sca, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "CPowScaKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching CPowScaKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// calculate out_1[(j+sh1)%n1][(i+sh0)%n0] = in_1[j][i] on device
cudaError_t ArrayOpCShift2d(cuComplex *out_1, cuComplex *in_1, unsigned int sh0, unsigned int sh1, unsigned int n0, unsigned int n1, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
//	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	CShift2dKernel<<<gridSize, blockSize>>>(out_1, in_1, sh0, sh1, n0, n1);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "CShift2dKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching CShift2dKernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
Error:
	return cudaStatus;
}


cudaError_t ArrayOpMulPP00(cuComplex *out_1, cuComplex *in_1, float *in_2, float *in_3, float dx, float dy, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

	// Launch the parallel kernel operation
	MulPhasePlate00Kernel <<<gridSize, blockSize>>>(out_1, in_1, in_2, in_3, dx, dy, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "MulPhasePlate00Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching MulPhasePlate00Kernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}


cudaError_t ArrayOpMulPP01(cuComplex *out_1, cuComplex *in_1, float *in_2, float *in_3, float dx, float dy, float dz, ArrayOpStats1 stats)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	int blockSize = stats.nBlockSize;	// block size 
	int gridSize = stats.nGridSize;		// grid size needed, based on input size

										// Launch the parallel kernel operation
	MulPhasePlate01Kernel<<<gridSize, blockSize>>>(out_1, in_1, in_2, in_3, dx, dy, dz, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "MulPhasePlate01Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - gridSize: %i, blockSize: %i\n", gridSize, blockSize);
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize failed after launching MulPhasePlate01Kernel: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

Error:
	return cudaStatus;
}

// Calculates the sum of a float array on device.
cudaError_t ArrayOpFSum(float &out_1, float * in_1, ArrayOpStats1 stats, int CPU_threshold)
{
	cudaError_t cudaStatus = cudaSuccess;
	unsigned int size = stats.uSize;
	int blockSize = (1024 < stats.nBlockSize) ? 1024 : stats.nBlockSize;
	int gridSize = (size + (blockSize * 2 - 1)) / (blockSize * 2); // determine gridSize locally
	int smemSize = (blockSize <= 32) ? 2 * blockSize * sizeof(float) : blockSize * sizeof(float);
	unsigned int ns1 = size;
	unsigned int nstage = 0;
	float dsum = 0.0, dc = 0.0, dy = 0.0, dt = 0.0;
	// init output
	out_1 = 0.f;
	// allocate temp output array on device, one slot for each block and preset with zeroes
	float * d_odata = NULL;
	float * d_idata = NULL;
	cudaStatus = cudaMalloc((void**)&d_odata, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage=1;  goto Error;	}
	cudaStatus = cudaMalloc((void**)&d_idata, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage=2;  goto Error; }
	cudaStatus = cudaMemset(d_odata, 0, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage=3;  goto Error; }
	cudaStatus = cudaMemset(d_idata, 0, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage=4;  goto Error; }
	// initial reduction (max blockSize 1024)
	// this is separated as it runs on the large input array, which do not want to change
	switch (blockSize)
	{
	case 1024:
		FAddReduceKernel<1024><<<gridSize, blockSize, smemSize>>>(d_odata, in_1, size);
		break;

	case 512:
		FAddReduceKernel<512><<<gridSize, blockSize, smemSize>>>(d_odata, in_1, size);
		break;

	case 256:
		FAddReduceKernel<256><<<gridSize, blockSize, smemSize>>>(d_odata, in_1, size);
		break;

	case 128:
		FAddReduceKernel<128><<<gridSize, blockSize, smemSize>>>(d_odata, in_1, size);
		break;

	case 64:
		FAddReduceKernel<64><<<gridSize, blockSize, smemSize>>>(d_odata, in_1, size);
		break;

	case 32:
		FAddReduceKernel<32><<<gridSize, blockSize, smemSize>>>(d_odata, in_1, size);
		break;

	case 16:
		FAddReduceKernel<16><<<gridSize, blockSize, smemSize>>>(d_odata, in_1, size);
		break;

	case  8:
		FAddReduceKernel<8><<<gridSize, blockSize, smemSize>>>(d_odata, in_1, size);
		break;

	case  4:
		FAddReduceKernel<4><<<gridSize, blockSize, smemSize>>>(d_odata, in_1, size);
		break;

	case  2:
		FAddReduceKernel<2><<<gridSize, blockSize, smemSize>>>(d_odata, in_1, size);
		break;

	case  1:
		FAddReduceKernel<1><<<gridSize, blockSize, smemSize>>>(d_odata, in_1, size);
		break;
	}
	cudaStatus = cudaGetLastError();
	if (cudaSuccess != cudaStatus) { nstage=5;  goto Error; }
	cudaStatus = cudaDeviceSynchronize();
	if (cudaSuccess != cudaStatus) { nstage=6;  goto Error; }
	// here, data is summed up to d_odata
	//
	// sum further while number of items > 1
	// this can be repeated until full reduction to a single float since it runs on
	// local arrays (d_odata, d_idata) of reduced size
	ns1 = gridSize;
	while (ns1 > (unsigned int)(1 + abs(CPU_threshold))) {
		gridSize = (ns1 + (blockSize * 2 - 1)) / (blockSize * 2); // update gridSize
		// prepare new input data on device
		cudaStatus = cudaMemcpy(d_idata, d_odata, sizeof(float)*ns1, cudaMemcpyDeviceToDevice);
		if (cudaSuccess != cudaStatus) { nstage=7;  goto Error; }
		// initialize output on device
		cudaStatus = cudaMemset(d_odata, 0, sizeof(float)*ns1);
		if (cudaSuccess != cudaStatus) { nstage=8;  goto Error; }
		switch (blockSize)
		{
		case 1024:
			FAddReduceKernel<1024><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 512:
			FAddReduceKernel<512><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 256:
			FAddReduceKernel<256><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 128:
			FAddReduceKernel<128><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 64:
			FAddReduceKernel<64><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 32:
			FAddReduceKernel<32><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 16:
			FAddReduceKernel<16><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case  8:
			FAddReduceKernel<8><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case  4:
			FAddReduceKernel<4><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case  2:
			FAddReduceKernel<2><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case  1:
			FAddReduceKernel<1><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;
		}
		cudaStatus = cudaGetLastError();
		if (cudaSuccess != cudaStatus) { nstage=9;  goto Error; }
		cudaStatus = cudaDeviceSynchronize();
		if (cudaSuccess != cudaStatus) { nstage=10;  goto Error; }
		ns1 = gridSize; // update size of the problem
		// d_odata contains the reduced chain of length ns1
		// when ns1 == 1, we are done with reducing
	}
	if (ns1 == 1) {
		// final retrieval of the sum
		cudaStatus = cudaMemcpy(&out_1, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaSuccess != cudaStatus) { nstage = 11;  goto Error; }
	}
	else {
		// do a Kohan Summation on the remaining items
		float *frest = (float*)malloc(sizeof(float)*ns1);
		cudaStatus = cudaMemcpy(frest, d_odata, sizeof(float)*ns1, cudaMemcpyDeviceToHost);
		if (cudaSuccess != cudaStatus) { nstage = 12;  goto Error; }
		for (unsigned int ir = 0; ir < ns1; ir++) {
			dy = frest[ir] - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
		out_1 = dsum;
		free(frest);
	}
	//
Error:
	if (NULL != d_odata) cudaFree(d_odata);
	if (NULL != d_idata) cudaFree(d_idata);
	if (nstage > 0) {
		fprintf(stderr, "CUDA error in ArrayOpFSum: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - procedure stage: %d\n", nstage);
		goto Error;
	}
	return cudaStatus;
}


// Calculates the sum of a float array on device using a mask.
cudaError_t ArrayOpMaskFSum(float &out_1, int * mask, float * in_1, ArrayOpStats1 stats, int CPU_threshold)
{
	cudaError_t cudaStatus = cudaSuccess;
	unsigned int size = stats.uSize;
	int blockSize = (1024 < stats.nBlockSize) ? 1024 : stats.nBlockSize;
	int gridSize = (size + (blockSize * 2 - 1)) / (blockSize * 2); // determine gridSize locally
	int smemSize = (blockSize <= 32) ? 2 * blockSize * sizeof(float) : blockSize * sizeof(float);
	unsigned int ns1 = size;
	unsigned int nstage = 0;
	float dsum = 0.0, dc = 0.0, dy = 0.0, dt = 0.0;
	// init output
	out_1 = 0.f;
	// allocate temp output array on device, one slot for each block and preset with zeroes
	float * d_odata = NULL;
	float * d_idata = NULL;
	cudaStatus = cudaMalloc((void**)&d_odata, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage = 1;  goto Error; }
	cudaStatus = cudaMalloc((void**)&d_idata, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage = 2;  goto Error; }
	cudaStatus = cudaMemset(d_odata, 0, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage = 3;  goto Error; }
	cudaStatus = cudaMemset(d_idata, 0, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage = 4;  goto Error; }
	// initial reduction with mask (max blockSize 1024)
	// this is separated as it runs on the large input array, which do not want to change
	switch (blockSize)
	{
	case 1024:
		MaskFAddReduceKernel<1024><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, size);
		break;

	case 512:
		MaskFAddReduceKernel<512><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, size);
		break;

	case 256:
		MaskFAddReduceKernel<256><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, size);
		break;

	case 128:
		MaskFAddReduceKernel<128><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, size);
		break;

	case 64:
		MaskFAddReduceKernel<64><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, size);
		break;

	case 32:
		MaskFAddReduceKernel<32><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, size);
		break;

	case 16:
		MaskFAddReduceKernel<16><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, size);
		break;

	case  8:
		MaskFAddReduceKernel<8><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, size);
		break;

	case  4:
		MaskFAddReduceKernel<4><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, size);
		break;

	case  2:
		MaskFAddReduceKernel<2><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, size);
		break;

	case  1:
		MaskFAddReduceKernel<1><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, size);
		break;
	}
	cudaStatus = cudaGetLastError();
	if (cudaSuccess != cudaStatus) { nstage = 5;  goto Error; }
	cudaStatus = cudaDeviceSynchronize();
	if (cudaSuccess != cudaStatus) { nstage = 6;  goto Error; }
	// here, data is summed up to d_odata
	//
	// Sum further while number of items > 1.
	// This can be repeated until full reduction to a single float since it runs on
	// local arrays (d_odata, d_idata) of reduced size.
	// And this is also done now without mask and the standard float reduction kernel.
	ns1 = gridSize;
	while (ns1 > (unsigned int)(1 + abs(CPU_threshold))) {
		gridSize = (ns1 + (blockSize * 2 - 1)) / (blockSize * 2); // update gridSize
																  // prepare new input data on device
		cudaStatus = cudaMemcpy(d_idata, d_odata, sizeof(float)*ns1, cudaMemcpyDeviceToDevice);
		if (cudaSuccess != cudaStatus) { nstage = 7;  goto Error; }
		// initialize output on device
		cudaStatus = cudaMemset(d_odata, 0, sizeof(float)*ns1);
		if (cudaSuccess != cudaStatus) { nstage = 8;  goto Error; }
		switch (blockSize)
		{
		case 1024:
			FAddReduceKernel<1024><<<gridSize, blockSize, smemSize >>>(d_odata, d_idata, ns1);
			break;

		case 512:
			FAddReduceKernel<512><<<gridSize, blockSize, smemSize >>>(d_odata, d_idata, ns1);
			break;

		case 256:
			FAddReduceKernel<256><<<gridSize, blockSize, smemSize >>>(d_odata, d_idata, ns1);
			break;

		case 128:
			FAddReduceKernel<128><<<gridSize, blockSize, smemSize >>>(d_odata, d_idata, ns1);
			break;

		case 64:
			FAddReduceKernel<64><<<gridSize, blockSize, smemSize >>>(d_odata, d_idata, ns1);
			break;

		case 32:
			FAddReduceKernel<32><<<gridSize, blockSize, smemSize >>>(d_odata, d_idata, ns1);
			break;

		case 16:
			FAddReduceKernel<16><<<gridSize, blockSize, smemSize >>>(d_odata, d_idata, ns1);
			break;

		case  8:
			FAddReduceKernel<8><<<gridSize, blockSize, smemSize >>>(d_odata, d_idata, ns1);
			break;

		case  4:
			FAddReduceKernel<4><<<gridSize, blockSize, smemSize >>>(d_odata, d_idata, ns1);
			break;

		case  2:
			FAddReduceKernel<2><<<gridSize, blockSize, smemSize >>>(d_odata, d_idata, ns1);
			break;

		case  1:
			FAddReduceKernel<1><<<gridSize, blockSize, smemSize >>>(d_odata, d_idata, ns1);
			break;
		}
		cudaStatus = cudaGetLastError();
		if (cudaSuccess != cudaStatus) { nstage = 9;  goto Error; }
		cudaStatus = cudaDeviceSynchronize();
		if (cudaSuccess != cudaStatus) { nstage = 10;  goto Error; }
		ns1 = gridSize; // update size of the problem
						// d_odata contains the reduced chain of length ns1
						// when ns1 == 1, we are done with reducing
	}
	if (ns1 == 1) {
		// final retrieval of the sum
		cudaStatus = cudaMemcpy(&out_1, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaSuccess != cudaStatus) { nstage = 11;  goto Error; }
	}
	else {
		// do a Kohan Summation on the remaining items
		float *frest = (float*)malloc(sizeof(float)*ns1);
		cudaStatus = cudaMemcpy(frest, d_odata, sizeof(float)*ns1, cudaMemcpyDeviceToHost);
		if (cudaSuccess != cudaStatus) { nstage = 12;  goto Error; }
		for (unsigned int ir = 0; ir < ns1; ir++) {
			dy = frest[ir] - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
		out_1 = dsum;
		free(frest);
	}
	//
Error:
	if (NULL != d_odata) cudaFree(d_odata);
	if (NULL != d_idata) cudaFree(d_idata);
	if (nstage > 0) {
		fprintf(stderr, "CUDA error in ArrayOpMaskFSum: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - procedure stage: %d\n", nstage);
		goto Error;
	}
	return cudaStatus;
}


// calculates the dot product of two float arrays on device using a mask: out_1 = SUM( in_1[imask]*in_2[imask] )_imask[i]
// - set stats.uSize to the size of the mask array
// - set stats.nBlockSize to as much as many threads you want, it will be capped to 1024.
// - use this to calculate integrated detector results with mask and detector sensitivity (in_2) from a intensity distrib. (in_1)
cudaError_t ArrayOpMaskFDot(float &out_1, int * mask, float * in_1, float * in_2, ArrayOpStats1 stats, int CPU_threshold)
{
	cudaError_t cudaStatus = cudaSuccess;
	unsigned int size = stats.uSize;
	int blockSize = (1024 < stats.nBlockSize) ? 1024 : stats.nBlockSize;
	int gridSize = (size + (blockSize * 2 - 1)) / (blockSize * 2); // determine gridSize locally
	int smemSize = (blockSize <= 32) ? 2 * blockSize * sizeof(float) : blockSize * sizeof(float);
	unsigned int ns1 = size;
	unsigned int nstage = 0;
	float dsum = 0.0, dc = 0.0, dy = 0.0, dt = 0.0;
	// init output
	out_1 = 0.f;
	// allocate temp output array on device, one slot for each block and preset with zeroes
	float * d_odata = NULL;
	float * d_idata = NULL;
	cudaStatus = cudaMalloc((void**)&d_odata, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage = 1;  goto Error; }
	cudaStatus = cudaMalloc((void**)&d_idata, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage = 2;  goto Error; }
	cudaStatus = cudaMemset(d_odata, 0, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage = 3;  goto Error; }
	cudaStatus = cudaMemset(d_idata, 0, sizeof(float)*gridSize);
	if (cudaSuccess != cudaStatus) { nstage = 4;  goto Error; }
	// initial reduction with mask (max blockSize 1024)
	// this is separated as it runs on the large input array, which do not want to change
	switch (blockSize)
	{
	case 1024:
		MaskFDotReduceKernel<1024><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, in_2, size);
		break;

	case 512:
		MaskFDotReduceKernel<512><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, in_2, size);
		break;

	case 256:
		MaskFDotReduceKernel<256><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, in_2, size);
		break;

	case 128:
		MaskFDotReduceKernel<128><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, in_2, size);
		break;

	case 64:
		MaskFDotReduceKernel<64><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, in_2, size);
		break;

	case 32:
		MaskFDotReduceKernel<32><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, in_2, size);
		break;

	case 16:
		MaskFDotReduceKernel<16><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, in_2, size);
		break;

	case  8:
		MaskFDotReduceKernel<8><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, in_2, size);
		break;

	case  4:
		MaskFDotReduceKernel<4><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, in_2, size);
		break;

	case  2:
		MaskFDotReduceKernel<2><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, in_2, size);
		break;

	case  1:
		MaskFDotReduceKernel<1><<<gridSize, blockSize, smemSize >>>(d_odata, mask, in_1, in_2, size);
		break;
	}
	cudaStatus = cudaGetLastError();
	if (cudaSuccess != cudaStatus) { nstage = 5;  goto Error; }
	cudaStatus = cudaDeviceSynchronize();
	if (cudaSuccess != cudaStatus) { nstage = 6;  goto Error; }
	// here, data is summed up to d_odata
	//
	// Sum further while number of items > 1.
	// This can be repeated until full reduction to a single float since it runs on
	// local arrays (d_odata, d_idata) of reduced size.
	// And this is also done now without mask and the standard float reduction kernel.
	ns1 = gridSize;
	while (ns1 > (unsigned int)(1 + abs(CPU_threshold))) {
		gridSize = (ns1 + (blockSize * 2 - 1)) / (blockSize * 2); // update gridSize
																  // prepare new input data on device
		cudaStatus = cudaMemcpy(d_idata, d_odata, sizeof(float)*ns1, cudaMemcpyDeviceToDevice);
		if (cudaSuccess != cudaStatus) { nstage = 7;  goto Error; }
		// initialize output on device
		cudaStatus = cudaMemset(d_odata, 0, sizeof(float)*ns1);
		if (cudaSuccess != cudaStatus) { nstage = 8;  goto Error; }
		switch (blockSize)
		{
		case 1024:
			FAddReduceKernel<1024><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 512:
			FAddReduceKernel<512><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 256:
			FAddReduceKernel<256><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 128:
			FAddReduceKernel<128><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 64:
			FAddReduceKernel<64><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 32:
			FAddReduceKernel<32><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case 16:
			FAddReduceKernel<16><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case  8:
			FAddReduceKernel<8><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case  4:
			FAddReduceKernel<4><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case  2:
			FAddReduceKernel<2><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;

		case  1:
			FAddReduceKernel<1><<<gridSize, blockSize, smemSize>>>(d_odata, d_idata, ns1);
			break;
		}
		cudaStatus = cudaGetLastError();
		if (cudaSuccess != cudaStatus) { nstage = 9;  goto Error; }
		cudaStatus = cudaDeviceSynchronize();
		if (cudaSuccess != cudaStatus) { nstage = 10;  goto Error; }
		ns1 = gridSize; // update size of the problem
						// d_odata contains the reduced chain of length ns1
						// when ns1 == 1, we are done with reducing
	}
	if (ns1 == 1) {
		// final retrieval of the sum
		cudaStatus = cudaMemcpy(&out_1, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaSuccess != cudaStatus) { nstage = 11;  goto Error; }
	}
	else {
		// do a Kohan Summation on the remaining items
		float *frest = (float*)malloc(sizeof(float)*ns1);
		cudaStatus = cudaMemcpy(frest, d_odata, sizeof(float)*ns1, cudaMemcpyDeviceToHost);
		if (cudaSuccess != cudaStatus) { nstage = 12;  goto Error; }
		for (unsigned int ir = 0; ir < ns1; ir++) {
			dy = frest[ir] - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
		out_1 = dsum;
		free(frest);
	}
	//
Error:
	if (NULL != d_odata) cudaFree(d_odata);
	if (NULL != d_idata) cudaFree(d_idata);
	if (nstage > 0) {
		fprintf(stderr, "CUDA error in ArrayOpMaskFSum: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - procedure stage: %d\n", nstage);
		goto Error;
	}
	return cudaStatus;
}


// calculates out_1[0] = SUM( in_1[i]*conjg( in_1[i] )) * sca on device
cudaError_t ArrayOpCPowSum(float &out_1, cuComplex *in_1, float sca, ArrayOpStats1 stats, int CPU_threshold)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	unsigned int nstage = 0;			// internal error stage code
	// allocate temporary output size and allocate on device
	size_t sz_tmp_out = sizeof(float)*size;
	float * tmp_out_1 = NULL;
	cudaStatus = cudaMalloc((void**)&tmp_out_1, sz_tmp_out);
	if (cudaStatus != cudaSuccess) { nstage = 1; goto Error; }
	cudaStatus = cudaMemset(tmp_out_1, 0, sz_tmp_out);
	if (cudaStatus != cudaSuccess) { nstage = 2; goto Error; }
	// Launch the parallel CPow operation: tmp_out_1[i] -> in_1[i]*conjg( in_1[i] )
	cudaStatus = ArrayOpCPow(tmp_out_1, in_1, stats);
	if (cudaStatus != cudaSuccess) { nstage = 3; goto Error; }
	// Launch the parallel float array summation: out_1 -> SUM( tmp_out_1[i] )
	cudaStatus = ArrayOpFSum(out_1, tmp_out_1, stats, CPU_threshold);
	if (cudaStatus != cudaSuccess) { nstage = 4; goto Error; }
	// Final scaling
	out_1 *= sca;
Error:
	if (tmp_out_1 != NULL) cudaFree(tmp_out_1);
	if (nstage > 0) {
		fprintf(stderr, "CUDA error in ArrayOpCPowSum: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - procedure stage: %d\n", nstage);
		goto Error;
	}
	return cudaStatus;
}


// calculates out_1 = SUM( in_1[i]*conjg( in_1[i] ) * in_2[i] ) * sca on device
cudaError_t ArrayOpCPowSumFMul(float &out_1, cuComplex *in_1, float *in_2, float sca, ArrayOpStats1 stats, int CPU_threshold)
{
	cudaError_t cudaStatus;
	unsigned int size = stats.uSize;	// input size
	unsigned int nstage = 0;			// internal error stage code
	size_t sz_tmp_out = sizeof(float)*size; // size of float temporaray working array on device
	float * tmp_out_1 = NULL;			// temporary working array on device

	// allocate temporary output size and allocate on device
	cudaStatus = cudaMalloc((void**)&tmp_out_1, sz_tmp_out);
	if (cudaStatus != cudaSuccess) { nstage = 1; goto Error; }
	cudaStatus = cudaMemset(tmp_out_1, 0, sz_tmp_out);
	if (cudaStatus != cudaSuccess) { nstage = 2; goto Error; }

	// Launch the parallel CPow operation: tmp_out_1[i] -> in_1[i]*conjg( in_1[i] )
	cudaStatus = ArrayOpCPow(tmp_out_1, in_1, stats);
	if (cudaStatus != cudaSuccess) { nstage = 3; goto Error; }
	// Launch the parallel float array multiplication: tmp_out_1[i] -> tmp_out_1[i] * in_2[i]
	cudaStatus = ArrayOpFFMul(tmp_out_1, tmp_out_1, in_2, stats);
	if (cudaStatus != cudaSuccess) { nstage = 4; goto Error; }
	// Launch the parallel float array summation: out_1 -> SUM( tmp_out_1[i] )
	cudaStatus = ArrayOpFSum(out_1, tmp_out_1, stats, CPU_threshold);
	if (cudaStatus != cudaSuccess) { nstage = 5; goto Error; }
	// Final scaling
	out_1 *= sca;
Error:
	if (tmp_out_1 != NULL) cudaFree(tmp_out_1);
	if (nstage > 0) {
		fprintf(stderr, "CUDA error in ArrayOpCPowSumFMul: %s\n", cudaGetErrorString(cudaStatus));
		fprintf(stderr, "  - procedure stage: %d\n", nstage);
		goto Error;
	}
	return cudaStatus;
}
