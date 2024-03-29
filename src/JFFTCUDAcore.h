//
// C++ header file: JFFTCUDAcore.h
// declaration of class CJFFTCUDAcore (implementation see JFFTCUDAcore.cpp)
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
// class CJFFTCUDAcore implements routines of Fourier transformations using
// libraries cufft.lib; cudart_static.lib
// You need to provide a version of the cuda runtime libraries ...
// - CUDA version 9.2: cudart64_92.dll, cufft64_92.dll
//
// The code is supposed to handle one type of Fourier transform only
// for single precision (float) complex data. The transformation dimension
// is limited to 1d, 2d, and 3d.
//
// - setup interface:
//   Use functions Init(...) and Deinit() to setup and de-initialize the
//   transformation code.
// - operational interface:
//   You can call forward and backward transformations with the function
//   FT() and IFT(), respectively. The Transformation will be "in place".
// - data transfer interface:
//   Set the data using SetData...() routines.
//   Get the data using GetData...() routines.
//
// - Transformation parameters:
//   m_ndim = number of dimensions
//   m_pdims = list of ndim values defining the array sizes
//   m_pcw = complex data (input and output)
//
#pragma once
#include "fcomplex.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include "cu\ArrayOps.cuh"
//#include "cu\FFTCallbacks.cuh"
//#include <windows.h>
//
#ifndef JFFTCUDACORE_H
#define JFFTCUDACORE_H
#define JFFTCUDACORE_DIM_MIN		(1)
#define JFFTCUDACORE_DIM_MAX		(3)
#endif // JFFTCUDACORE_H
//
class CJFFTCUDAcore
{
public:
	CJFFTCUDAcore();
	~CJFFTCUDAcore();
	
protected:
	cudaError m_cuerrLast; // last cuda error code
	cufftResult m_cufftresLast; // last cuda fft result
	char* m_scufftdllname; // linked cufft dll
	int m_nstatus; // core status: 0: not initialized, 1: initialized
	int m_ndim; // core plan number of dimensions: <=0: invalid, 1: 1d, 2: 2d, 3: 3d, >=4: invalid
	int* m_pdims; // core plan size of the array
	// change added separate plans for forward and backward fft for different callback handling, add pointer to callback parameters
	cufftHandle m_fw_plan; // core forward plan object
	//cufftHandle m_bw_plan; // core inverse plan object
	//fftCallbackParams* m_d_forwardCallbackParams;
	// end change
	cufftComplex* m_pcw; // core complex working array
	ArrayOpStats1 m_aos1dMult; // stats for 1d cuda array complex multiplication operation

	void PostCUDAError(const char* smsg, cudaError code); // post cuda error to cerr

	
public:
	
	// Initialization interface:
	
	// de-initializes the object
	void Deinit(void);
	// initialize the plan and internal array on current cuda device
	int Init(int ndim, int * pdims); 

	// Transformation operation interface:

	//// forward FFT with the current forward plan and data in the core
	//// using load and store callbacks with data provided by cb_ld_data and cb_st_data
	//int FT(cuComplex* cb_ld_data, cuComplex* cb_st_data);

	// forward FFT with the current forward plan and data in the core
	int FT(void);
	// inverse FFT with the current inverse plan and data in the core
	int IFT(void); 
	// forward FFT with the current plan and device data src
	int FT_d(cufftComplex * src);
	// inverse FFT with the current plan and device data src
	int IFT_d(cufftComplex * src);
	// sets the core data values to zero
	int Zero(void);
	// scale the complex data by a given factor
	int Scale(float scale);
	// multiply complex values on device
	int MultiplyC_d(cuComplex * src);
	// multiply complex values on device from a 2d sub-frame
	int MultiplySub2dC_d(cuComplex * src, int nsub0, int nsub1);
	// multiply complex values given as re,im aligned list from host on device
	int MultiplyC(fcmplx * src); 
	// multiply float values on device
	int MultiplyF(float * src);
	// cyclic shift pixels
	int CShift2d(int nsh0, int nsh1);
	
	// Data transfer interface:
	
	// returns the core status: 0 = not ready, 1 = initialized
	int GetStatus(void);
	// returns the current device id
	int GetDevice(void);
	// return the number of core data dimensions
	int GetDimension(void);
	// returns the size along each dimensions
	int GetDimSize(int * pdims);
	// returns the block size used
	int GetBlockSize(void);
	// returns total number of data items involved in the transformation
	size_t GetDataSize(void);
	// returns device memory usage for fft
	size_t GetFFTMemUsage(void);
	// returns device memory usage for multiplication
	size_t GetMulMemUsage(void);
	// sets the cuda device id to be used for FFTs and ArrayOps
	int SetCUDADevice(int ndev = 0);
	// sets data by copying from device memory
	int SetDataC_d(cuComplex * src);
	// sets data, assuming float re,im aligned source items from host memory
	int SetDataC(fcmplx * src);
	// sets data (re,0) from host
	int SetDataRe(float * src);
	// sets data (0,im) from host
	int SetDataIm(float * src);
	// returns the device address of the core buffer (only for accessing data, do not change the allocation state)
	cuComplex* GetData(void);
	// gets data to host, assuming float re,im aligned destination items
	int GetDataC(fcmplx * dst);
	// gets data to host, assuming float (real) destination items
	int GetDataRe(float * dst);
	// gets data to host, assuming float (imag) destination items
	int GetDataIm(float * dst);
	// gets data to host, assuming float (abs) destination items
	int GetDataAbs(float * dst);
	// gets data to host, assuming float (arg) destination items
	int GetDataArg(float * dst);
	// gets a pointer to the core data on device memory
	// ! Warning ! Use carefully! Do not change allocation status!
	cuComplex* GetData_d();
	// gets a copy of the data on device memory
	int GetDataC_d(cuComplex * dst);
	// sets data (0,im) from device source
	int SetDataIm_d(float * src);
	// sets data (re,0) from device source
	int SetDataRe_d(float * src);
	// gets data (0,im) to device buffer
	int GetDataIm_d(float * dst);
	// gets data (re,0) to device buffer
	int GetDataRe_d(float * dst);
	// gets a copy of the data power to device buffer
	int GetDataPow_d(float * dst);
	// gets a copy of the data power to host buffer
	int GetDataPow(float * dst);
	// gets total data power
	int GetDataTotalPow(float &pow);
};

