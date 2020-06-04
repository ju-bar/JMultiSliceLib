//
// CUDA C header file: FFTCallbacks.cuh
// declaration of FFT callback functions used by the cuFFT callback API
// (implementation, see FFTCallbacks.cu)
//
// Author: Juri Barthel (juribarthel@gmail.com)
// Date  : 2020-04-21
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
//#ifdef __INTELLISENSE__
//
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufftXt.h>
//
//#endif
#include <cufft.h>

#ifndef __FFT_CALLBACK_DEF__
#define __FFT_CALLBACK_DEF__
// define parameter structure for cufft-callbacks
typedef struct {
	cuComplex* load; // pointer to a data array processed on load before the fft
	cuComplex* store; // pointer to a data array processed on store after the fft
} fftCallbackParams;
#endif
// applies complex multiplication during load
__device__ static cufftComplex MultCLoadCallback(void* in_ptr, size_t index, void* params_ptr, void*)
{
    cufftComplex* in_data = (cufftComplex*)in_ptr;
    cufftComplex in_tmp = in_data[index];
    cufftComplex in;
    fftCallbackParams* params = (fftCallbackParams*)params_ptr;
    cufftComplex* load_data = params->load;
    if (load_data != NULL) {
        cufftComplex dat = params->load[index];
        in.x = in_tmp.x * dat.x - in_tmp.y * dat.y;
        in.y = in_tmp.x * dat.y + in_tmp.y * dat.x;
    }
    else {
        in = in_tmp;
    }
    return in;
}
// applies complex multiplication during store
__device__ static void MultCStoreCallback(void* out_ptr, size_t index, cufftComplex in, void* params_ptr, void*)
{
    fftCallbackParams* params = (fftCallbackParams*)params_ptr;
    cufftComplex out;
    cufftComplex* store_data = params->store;
    // do propagator multiplication
    if (store_data != NULL) {
        cufftComplex store = store_data[index];
        out.x = in.x * store.x - in.y * store.y;
        out.y = in.x * store.y + in.y * store.x;
    }
    else {
        out = in;
    }
    cufftComplex* out_data = (cufftComplex*)out_ptr;
    out_data[index] = out;
}
