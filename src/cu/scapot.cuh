//
// CUDA C header file: scapot.cuh
// declaration of device function for the calculation of scattering potentials and object transmission functions operations to be performed on the current GPU device
// (implementation, see scapot.cu)
//
// Author: Juri Barthel (juribarthel@gmail.com)
// Date  : 2020-05-28
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
#ifndef DRPROBE_SCAPOT_CUH
#define DRPROBE_SCAPOT_CUH
//
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
//
#include <cufft.h>

__device__ __forceinline__ cuComplex cu_cexpf(cuComplex z)
{
    cuComplex res;
    float t = expf(z.x);
    res.x = cosf(z.y) * t;
    res.y = sinf(z.y) * t;
    return res;
}

__global__ void AddFormFactor2dKernel(cuComplex* pot, cuComplex* ff, float* qx, float* qy, float* ap, float x, float y, float occ, unsigned int nx, unsigned int size);
__global__ void AddIonFormFactor2dKernel(cuComplex* pot, cuComplex* ff, cuComplex* pio, float* qx, float* qy, float* ap, float x, float y, float crg, float occ, unsigned int nx, unsigned int size);
__global__ void AddDampedIonFormFactor2dKernel(cuComplex* pot, cuComplex* ff, cuComplex* pio, float* qx, float* qy, float* ap, float x, float y, float crg, float biso, float occ, unsigned int nx, unsigned int size);
__global__ void PotToPgrKernel(cuComplex* pgr, cuComplex* pot, cuComplex cisig, unsigned int size);

cudaError_t AddFormFactor2d(cuComplex* pot, cuComplex* ff, float* qx, float* qy, float* ap, float x, float y, float occ, unsigned int nx, unsigned int ny);
cudaError_t AddIonFormFactor2d(cuComplex* pot, cuComplex* ff, cuComplex* pio, float* qx, float* qy, float* ap, float x, float y, float crg, float occ, unsigned int nx, unsigned int ny);
cudaError_t AddDampedIonFormFactor2d(cuComplex* pot, cuComplex* ff, cuComplex* pio, float* qx, float* qy, float* ap, float x, float y, float crg, float biso, float occ, unsigned int nx, unsigned int ny);
cudaError_t PotToPgr(cuComplex* pgr, cuComplex* pot, cuComplex* cisig, unsigned int nx, unsigned int ny);

#endif // DRPROBE_SCAPOT_CUH

