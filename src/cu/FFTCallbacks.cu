// file FFTCallbacks.cu
// author: Juri Barthel, ju.barthel@fz-juelich.de
// April 21, 2020
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
// Implementation for cufft callbacks
//

#include <stdio.h>
#include "FFTCallbacks.cuh"



__device__ cufftComplex MultCLoadCallback(void* in_ptr, size_t index, void* params_ptr, void*)
{
    cufftComplex* in_data = (cufftComplex*)in_ptr;
    cufftComplex in_tmp = in_data[index];
    cufftComplex in;
    fftCallbackParams* params = (fftCallbackParams*) params_ptr;
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

__device__ void MultCStoreCallback(void* out_ptr, size_t index, cufftComplex in, void* params_ptr, void*)
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

