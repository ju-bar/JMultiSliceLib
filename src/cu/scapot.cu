// file scapot.cu
// author: Juri Barthel, ju.barthel@fz-juelich.de
// 2018 - 2020
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
// Implementation for CUDA calculations of scattering potentials
//


#include "scapot.cuh"
#include <iostream>


__global__ void AddFormFactor2dKernel(cuComplex* pot, cuComplex* ff, float* qx, float* qy, float* ap, float x, float y, float occ, unsigned int nx, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuComplex tf;
	float pha;
	unsigned int ix, iy;
	if (idx < size) {
		ix = idx % nx;
		iy = (idx - ix) / nx;
		pha = -6.2831853f * (x * qx[ix] + y * qy[iy]);
		tf.x = cosf(pha) * ap[idx] * occ; tf.y = sinf(pha) * ap[idx] * occ;
		pot[idx] = cuCaddf(pot[idx], cuCmulf(ff[idx], tf));
	}
}

__global__ void AddIonFormFactor2dKernel(cuComplex* pot, cuComplex* ff, cuComplex* pio, float* qx, float* qy, float* ap, float x, float y, float crg, float occ, unsigned int nx, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuComplex tf;
	cuComplex ctmp;
	float pha;
	unsigned int ix, iy;
	if (idx < size) {
		ix = idx % nx;
		iy = (idx - ix) / nx;
		pha = -6.2831853f * (x * qx[ix] + y * qy[iy]);
		tf.x = cosf(pha) * ap[idx] * occ; tf.y = sinf(pha) * ap[idx] * occ;
		ctmp.x = ff[idx].x + pio[idx].x * crg; ctmp.y = ff[idx].y + pio[idx].y * crg;
		pot[idx] = cuCaddf(pot[idx], cuCmulf(ctmp, tf));
	}
}

__global__ void AddDampedIonFormFactor2dKernel(cuComplex* pot, cuComplex* ff, cuComplex* pio, float* qx, float* qy, float* ap, float x, float y, float crg, float biso, float occ, unsigned int nx, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuComplex tf;
	cuComplex ctmp;
	float pha, dwf;
	unsigned int ix, iy;
	if (idx < size) {
		ix = idx % nx;
		iy = (idx - ix) / nx;
		pha = -6.2831853f * (x * qx[ix] + y * qy[iy]);
		dwf = expf(-0.25f * biso * (qx[ix] * qx[ix] + qy[iy] * qy[iy]));
		tf.x = cosf(pha) * ap[idx] * occ; tf.y = sinf(pha) * ap[idx] * occ;
		ctmp.x = ff[idx].x + pio[idx].x * crg * dwf; ctmp.y = ff[idx].y + pio[idx].y * crg * dwf;
		pot[idx] = cuCaddf(pot[idx], cuCmulf(ctmp, tf));
	}
}

__global__ void PotToPgrKernel(cuComplex* pgr, cuComplex* pot, cuComplex cisig, unsigned int size)
{
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	cuComplex ctmp;
	if (idx < size) {
		ctmp = cuCmulf(pot[idx], cisig);
		pgr[idx] = cu_cexpf(ctmp);
	}
}





cudaError_t AddFormFactor2d(cuComplex* pot, cuComplex* ff, float* qx, float* qy, float* ap, float x, float y, float occ, unsigned int nx, unsigned int ny)
{
	cudaError_t cudaStatus;
	unsigned int size = nx * ny; // input size
	int blockSize = 256; // default block size 
	int gridSize = (size + blockSize - 1) / blockSize; // grid size needed, based on input size

	// Launch the parallel kernel operation
	AddFormFactor2dKernel<<<gridSize,blockSize>>>(pot, ff, qx, qy, ap, x, y, occ, nx, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "AddFormFactor2dKernel launch failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "cudaDeviceSynchronize failed after launching AddFormFactor2dKernel: " << cudaGetErrorString(cudaStatus) << std::endl;
		goto Error;
	}

Error:
	return cudaStatus;
}


cudaError_t AddIonFormFactor2d(cuComplex* pot, cuComplex* ff, cuComplex* pio, float* qx, float* qy, float* ap, float x, float y, float crg, float occ, unsigned int nx, unsigned int ny)
{
	cudaError_t cudaStatus;
	unsigned int size = nx * ny; // input size
	int blockSize = 256; // default block size 
	int gridSize = (size + blockSize - 1) / blockSize; // grid size needed, based on input size

	// Launch the parallel kernel operation
	AddIonFormFactor2dKernel<<<gridSize,blockSize>>> (pot, ff, pio, qx, qy, ap, x, y, crg, occ, nx, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "AddIonFormFactor2dKernel launch failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "cudaDeviceSynchronize failed after launching AddIonFormFactor2dKernel: " << cudaGetErrorString(cudaStatus) << std::endl;
		goto Error;
	}

Error:
	return cudaStatus;
}


cudaError_t AddDampedIonFormFactor2d(cuComplex* pot, cuComplex* ff, cuComplex* pio, float* qx, float* qy, float* ap, float x, float y, float crg, float biso, float occ, unsigned int nx, unsigned int ny)
{
	cudaError_t cudaStatus;
	unsigned int size = nx * ny; // input size
	int blockSize = 256; // default block size 
	int gridSize = (size + blockSize - 1) / blockSize; // grid size needed, based on input size

	// Launch the parallel kernel operation
	AddDampedIonFormFactor2dKernel<<<gridSize,blockSize>>>(pot, ff, pio, qx, qy, ap, x, y, crg, biso, occ, nx, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "AddDampedIonFormFactor2dKernel launch failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "cudaDeviceSynchronize failed after launching AddDampedIonFormFactor2dKernel: " << cudaGetErrorString(cudaStatus) << std::endl;
		goto Error;
	}

Error:
	return cudaStatus;
}

cudaError_t PotToPgr(cuComplex* pgr, cuComplex* pot, cuComplex* cisig, unsigned int nx, unsigned int ny)
{
	cudaError_t cudaStatus;
	unsigned int size = nx * ny; // input size
	int blockSize = 256; // default block size 
	int gridSize = (size + blockSize - 1) / blockSize; // grid size needed, based on input size
	cuComplex cs = *cisig;

	// Launch the parallel kernel operation
	PotToPgrKernel<<<gridSize,blockSize>>>(pgr, pot, cs, size);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "PotToPgrKernel launch failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		goto Error;
	}

	// synchronize threads and wait for all to be finished
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "cudaDeviceSynchronize failed after launching PotToPgrKernel: " << cudaGetErrorString(cudaStatus) << std::endl;
		goto Error;
	}

Error:
	return cudaStatus;
}