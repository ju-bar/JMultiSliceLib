//
// C++ header file: JFFTWcore.h
// declaration of class CJFFTWcore (implementation see JFFTWcore.cpp)
//
// Author: Juri Barthel (juribarthel@gmail.com)
// Date  : 2018-03-17
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
// class CJFFTWcore implements routines of Fourier transformations using
// libraries libfftw3f-3.lib
//
// The code is supposed to handle one type of Fourier transform only
// for single precision (float) complex data.
//
// - setup interface:
//   Use functions Init(...) and Deinit() to setup and de-initialize the
//   transformation code.
// - operational interface:
//   You can call forward and backward transformations with the function
//   FT() and IFT(), respectively. The Transformation will be "in place".
//   Each call increases the power by the square of the number of items.
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
//
#include "fcomplex.h"
#include "fftw3.h"
//
#ifndef JFFTWCORE_H
#define JFFTWCORE_H
#define JFFTWCORE_DIM_MIN		(1)
#endif // JFFTWCORE_H
//
class CJFFTWcore
{
public:
	CJFFTWcore();
	~CJFFTWcore();
protected:
	int m_nstatus; // core status: 0: not initialized, 1: initialized
	int m_ndim; // core plan number of dimensions: <=0: invalid, 1: 1d, 2: 2d, 3: 3d, >=4: invalid
	int m_nplanflag; // core plan flags (see FFTW documentation), FFTW_ESTIMATE is used by default
	int* m_pdims; // core plan size of the array
	fftwf_plan m_planf; // core plan object for forward transformation
	fftwf_plan m_planb; // core plan object for backward transformation
	fftwf_complex* m_pcw; // core complex working array
public:
	// Initialization interface:
	void Deinit(void); // de-initializes the object
	void CleanupFFTW(void); // cleans the FFWT module (use with care)
	int Init(int ndim, int * pdims, int nplanflag = FFTW_ESTIMATE); // initialize the plan and internal array
	// Transformation operation interface:
	int FT(void); // does a forward FFT with the current plan and data in the core
	int IFT(void); // does a backward FFT with the current plan and data in the core
	int Zero(void); // sets the core data values to zero
	int Conjugate(void); // applies complex conjugate to the data
	int Scale(float sca); // scales the data
	int MultiplyReal(float * src); // multiply with a list of real values of the same size
	int MultiplyC(float * src); // multiply complex values given as re,im aligned list
	int MultiplyC(fcmplx * src); // multiply complex values
	int MultiplyC(fftwf_complex * src); // multiply complex values
	int MultiplySub2dC(float * src, int nsub0, int nsub1); // multiply complex values given as re,im aligned list from an array of different size. Periodic wrap around applied.
	int MultiplySub2dC(fcmplx * src, int nsub0, int nsub1); // multiply complex values given as re,im aligned list from an array of different size. Periodic wrap around applied.
	int MultiplySub2dC(fftwf_complex * src, int nsub0, int nsub1); // multiply complex values given as re,im aligned list from an array of different size. Periodic wrap around applied.
	// Data transfer interface:
	int GetStatus(void); // returns the core status
	int GetDimension(void); // return the number of core data dimensions
	int GetDimSize(int * pdims); // returns the size along each dimensions
	int GetFlag(void); // returns the transformation plan creation flag
	size_t GetDataSize(void); // returns total number of data items involved in the transformation
	int SetDataC(float * src); // sets data, assuming float re,im aligned source items
	int SetDataC(fcmplx * src); // sets data, assuming complex<float> source items
	int SetDataC(fftwf_complex * src); // sets data, assuming fftwf_complex source items
	int SetDataRe(float * src); // sets data, assuming float (real) source items
	int SetDataIm(float * src); // sets data, assuming float (imag) source items
	int SetDataReInt(int * src); // sets data, assuming int (real) source items
	int SetDataImInt(int * src); // sets data, assuming int (imag) source items
	int GetDataC(float * dst); // gets data, assuming float re,im aligned destination items
	int GetDataC(fcmplx * dst); // gets data, assuming complex<float> destination items
	int GetDataC(fftwf_complex * dst); // gets data, assuming fftwf_complex destination items
	int GetDataRe(float * dst); // gets data, assuming float (real) destination items
	int GetDataIm(float * dst); // gets data, assuming float (imag) destination items
	int GetDataAbs(float * dst); // gets data, assuming float (abs) destination items
	int GetDataArg(float * dst); // gets data, assuming float (arg) destination items
	int GetDataPow(float * dst); // gets data, assuming float (abs**2) destination items
	float GetDataTotalPow(void); // returns the total power of the data
	float GetDataTotalRe(void); // returns the total of the data real part
	float GetDataTotalIm(void); // returns the total of the data imaginary part
};

