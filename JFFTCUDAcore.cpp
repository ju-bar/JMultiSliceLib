//
// C++ source file: JFFTCUDAcore.cpp
// implementation of class JFFTCUDAcore (declaration see JFFTCUDAcore.h)
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
// class JFFTCUDAcore implements routines of Fourier transformations using
// libraries cufft.lib; cudart_static.lib
//
// The code is supposed to handle one type of Fourier transform only
// for single precision (float) complex data.
//
#include "stdafx.h"
#include "JFFTCUDAcore.h"
//
using namespace std;
//


CJFFTCUDAcore::CJFFTCUDAcore()
{
	m_nstatus = 0;
	m_ndim = 0;
	m_pdims = NULL;
	m_pcw = NULL;
	m_plan = NULL;
}


CJFFTCUDAcore::~CJFFTCUDAcore()
{
	Deinit();
}

void CJFFTCUDAcore::PostCUDAError(char* smsg, cudaError code)
{
	if (code > 0) {
		cerr << "Error: " << smsg << ", code: " << code << endl;
		cerr << "     - " << cudaGetErrorString(code) << endl;
	}
}


/* ----------------------------------------------------------------------------- */
/*                                                                               */
/* >>> CUFFT WRAPPER <<<                                                         */
/*                                                                               */
/* ----------------------------------------------------------------------------- */



/* ----------------------------------------------------------------------------- */
/*                                                                               */
/* >>> SETUP INTERFACE <<<                                                       */
/*                                                                               */
/* ----------------------------------------------------------------------------- */

void CJFFTCUDAcore::Deinit(void)
{
	if (m_plan != NULL) m_cufftresLast = cufftDestroy(m_plan);
	if (m_pdims !=NULL) free(m_pdims);
	if (m_pcw != NULL) cudaFree(m_pcw);
	m_pdims = NULL;
	m_pcw = NULL;
	m_plan = NULL;
	m_nstatus = 0;
	m_ndim = 0;
}

int CJFFTCUDAcore::Init(int ndim, int * pdims)
{
	int nconsistent = 0;
	int i = 0;
	int nconsistent_target = 1 + ndim; // this is what should exist already
	//
	if (ndim<JFFTCUDACORE_DIM_MIN || ndim>JFFTCUDACORE_DIM_MAX) {
		cerr << "Error(JFFTCUDAcore::Init): Cannot initialize, invalid dimension." << endl;
		return 1;
	}
	//
	if (m_nstatus > 0) { // object is initialized // check consistency with init
		if (ndim == m_ndim) nconsistent++; // same dimension ?
		if (nconsistent == 1) { // same size and device?
			for (i = 0; i < m_ndim; i++) { // check for same size on each dimension
				if (pdims[i] == m_pdims[i]) nconsistent++;
			}
		}
	}
	// In case the following check doesn't fail (nconsistent == nconsistent_target),
	// the current core is already initialized for the requested transformation parameters
	// and no further action is required. This will never happen if the core is not
	// initialized before (m_status==0) or if any of the parameters has changed.
	if (nconsistent != nconsistent_target) { // consistency target failed
		// re-initialize with new parameters
		// - de-init the current core parameters
		Deinit(); // this makes also sure that all allocations on a previous device are cleared
		// - new dimension
		m_ndim = ndim;
		// - prepare new arrays
		m_pdims = (int*)malloc((size_t)(m_ndim * sizeof(int)));
		// - new number of complex data items
		int nd = 1;
		for (i = 0; i < m_ndim; i++) { // loop over all dimensions
			m_pdims[i] = pdims[i]; // store dimension length
			nd = nd*m_pdims[i]; // multiply to total number of items
		}
		//// - set active device
		//if (0 < SetCUDADevice(ndev)) {
		//	cout << "Error(JFFTCUDAcore): Failed to initialize, code: " << m_cuerrLast << endl;
		//	return 2;
		//}
		// - allocate new transformation array
		m_cuerrLast = cudaMalloc((void**)&m_pcw, sizeof(cufftComplex)*nd);
		if (cudaSuccess != m_cuerrLast) {
			PostCUDAError("(JFFTCUDAcore::Init): Failed to allocate transformation array", m_cuerrLast);
			return 3;
		}
		//// - prepare new transformation plan
		//m_cufftresLast = cufftSetCompatibilityMode(m_plan, CUFFT_COMPATIBILITY_FFTW_PADDING);
		//if (CUFFT_SUCCESS != m_cufftresLast) {
		//	cout << "Error(JFFTCUDAcore): Failed to set fftw compatibility mode for plan, code: " << m_cufftresLast << endl;
		//	return 4;
		//}
		switch (ndim) {
		case 1:
			m_cufftresLast = cufftPlan1d(&m_plan, m_pdims[0], CUFFT_C2C, 1);
			break;
		case 2:
			m_cufftresLast = cufftPlan2d(&m_plan, m_pdims[1], m_pdims[0], CUFFT_C2C);
			break;
		case 3:
			m_cufftresLast = cufftPlan3d(&m_plan, m_pdims[2], m_pdims[1], m_pdims[0], CUFFT_C2C);
			break;
		}
		if (CUFFT_SUCCESS != m_cufftresLast) {
			cerr << "Error(JFFTCUDAcore::Init): Failed to initialize plan, code: " << m_cufftresLast << endl;
			return 5;
		}
		// - determine optimum cuda block size for linear multiplications
		m_aos1dMult.uSize = nd;
		if (0 < GetOptimizedMultStats(&m_aos1dMult.uSize, &m_aos1dMult.nBlockSize, &m_aos1dMult.nGridSize)) {
			PostCUDAError("(JFFTCUDAcore::Init): Failed to determine optimum block size", m_cuerrLast);
			return 6;
		}
		// - update status
		m_nstatus = 1;
		// - re-initialize the transformation array as the plan creation may have used it
		if (0 < Zero()) {
			cerr << "Error(JFFTCUDAcore::Init): Failed to initialize transformation array." << endl;
			return 7;
		}

	}
	return 0;
}


int CJFFTCUDAcore::SetCUDADevice(int ndev)
{
	int imax = 0;
	int idev = ndev;
	m_cuerrLast = cudaGetDeviceCount(&imax);
	if (cudaSuccess != m_cuerrLast) {
		PostCUDAError("(JFFTCUDAcore::SetCUDADevice): Failed to get current cuda device", m_cuerrLast);
		return 1;
	}
	if (imax > 0 && idev >= 0 && idev < imax) {
		m_cuerrLast = cudaSetDevice(idev);
		if (cudaSuccess != m_cuerrLast) {
			PostCUDAError("(JFFTCUDAcore::SetCUDADevice): Failed to set cuda device", m_cuerrLast);
			return 1;
		}
	}
	else {
		cerr << "Error(JFFTCUDAcore::SetCUDADevice): Requested cuda device is invalid." << endl;
		return 2;
	}
	return 0;
}



/* ----------------------------------------------------------------------------- */
/*                                                                               */
/* >>> OPERATION INTERFACE <<<                                                   */
/*                                                                               */
/* ----------------------------------------------------------------------------- */


int CJFFTCUDAcore::FT(void)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::FT): Cannot transform, not initialized." << endl;
		return 1;
	}
	m_cufftresLast = cufftExecC2C(m_plan, m_pcw, m_pcw, CUFFT_FORWARD);
	if (CUFFT_SUCCESS != m_cufftresLast) {
		cerr << "Error(JFFTCUDAcore::FT): Forward transformation failed, code: " << m_cufftresLast << endl;
		return 2;
	}
	m_cuerrLast = cudaDeviceSynchronize();
	if (cudaSuccess != m_cuerrLast) {
		PostCUDAError("(JFFTCUDAcore::FT): Failed to snychronize device", m_cuerrLast);
		return 3;
	}
	return 0;
}


int CJFFTCUDAcore::IFT(void)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::IFT): Cannot transform, not initialized." << endl;
		return 1;
	}
	m_cufftresLast = cufftExecC2C(m_plan, m_pcw, m_pcw, CUFFT_INVERSE);
	if (CUFFT_SUCCESS != m_cufftresLast) {
		cerr << "Error(JFFTCUDAcore::IFT): Inverse transformation failed, code: " << m_cufftresLast << endl;
		return 2;
	}
	m_cuerrLast = cudaDeviceSynchronize();
	if (cudaSuccess != m_cuerrLast) {
		PostCUDAError("(JFFTCUDAcore::IFT): Failed to snychronize device", m_cuerrLast);
		return 3;
	}
	return 0;
}


int CJFFTCUDAcore::FT_d(cufftComplex * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::FT_d): Cannot transform, not initialized." << endl;
		return 1;
	}
	m_cufftresLast = cufftExecC2C(m_plan, src, src, CUFFT_FORWARD);
	if (CUFFT_SUCCESS != m_cufftresLast) {
		cerr << "Error(JFFTCUDAcore::FT_d): Forward transformation failed, code: " << m_cufftresLast << endl;
		return 2;
	}
	m_cuerrLast = cudaDeviceSynchronize();
	if (cudaSuccess != m_cuerrLast) {
		PostCUDAError("(JFFTCUDAcore::FT_d): Failed to snychronize device", m_cuerrLast);
		return 3;
	}
	return 0;
}


int CJFFTCUDAcore::IFT_d(cufftComplex * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::IFT_d): Cannot transform, not initialized." << endl;
		return 1;
	}
	m_cufftresLast = cufftExecC2C(m_plan, src, src, CUFFT_INVERSE);
	if (CUFFT_SUCCESS != m_cufftresLast) {
		cerr << "Error(JFFTCUDAcore::IFT_d): Inverse transformation failed, code: " << m_cufftresLast << endl;
		return 2;
	}
	m_cuerrLast = cudaDeviceSynchronize();
	if (cudaSuccess != m_cuerrLast) {
		PostCUDAError("(JFFTCUDAcore::IFT_d): Failed to snychronize device", m_cuerrLast);
		return 3;
	}
	return 0;
}


int CJFFTCUDAcore::Zero(void)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::Zero): Cannot zero data, not initialized." << endl;
		return 1;
	}
	size_t nbytes = sizeof(cufftComplex)*GetDataSize();
	m_cuerrLast = cudaMemset(m_pcw, 0, nbytes);
	if (cudaSuccess != m_cuerrLast) {
		PostCUDAError("(JFFTCUDAcore::Zero): Failed to set device memory to zero", m_cuerrLast);
		return 2;
	}
	return 0;
}


int CJFFTCUDAcore::Scale(float scale)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::Scale): Cannot scale data, not initialized." << endl;
		return 1;
	}
	// multiply on device pointwise (m_pcw * pcsrvdev)
	m_cuerrLast = ArrayOpSca(m_pcw, m_pcw, scale, m_aos1dMult);
	if (cudaSuccess != m_cuerrLast) {
		PostCUDAError("(JFFTCUDAcore::Scale): Failed to scale on device memory", m_cuerrLast);
		return 2;
	}
	return 0;
}


int CJFFTCUDAcore::MultiplyF(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::MultiplyF): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize(); // get data size (number of items)
	size_t nbytes = nd * sizeof(float);
	// prepare to copy source from host to the device
	float * pfsrcdev = NULL;
	m_cuerrLast = cudaMalloc((void**)&pfsrcdev, nbytes);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::MultiplyF): Failed to allocate on device", m_cuerrLast);
		return 3;
	}
	// copy from host to device
	m_cuerrLast = cudaMemcpy(pfsrcdev, src, nbytes, ::cudaMemcpyHostToDevice);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::MultiplyF): Failed to copy from host to device", m_cuerrLast);
		return 4;
	}
	// multiply on device pointwise (m_pcw * pcsrvdev)
	m_cuerrLast = ArrayOpFMul(m_pcw, m_pcw, pfsrcdev, m_aos1dMult);
	// destroy source helper on device
	m_cuerrLast = cudaFree(pfsrcdev);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::MultiplyF): Failed to destroy helper array on device", m_cuerrLast);
		return 5;
	}
	return 0;
}


int CJFFTCUDAcore::MultiplyC_d(cuComplex * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	// multiply on device pointwise (m_pcw * src)
	m_cuerrLast = ArrayOpMul(m_pcw, src, m_pcw, m_aos1dMult);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore): Failed to multiply complex arrays on device", m_cuerrLast);
		return 6;
	}
	return 0;
}

int CJFFTCUDAcore::MultiplySub2dC_d(cuComplex * src, int nsub0, int nsub1)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	if (m_ndim != 2) // error: invalid core number of dimensions
		return 3;
	// multiply on device stripe-pointwise (m_pcw * src)
	m_cuerrLast = ArrayOpMulSub2d(m_pcw, m_pcw, src, (unsigned)m_pdims[0], (unsigned)m_pdims[0], (unsigned)nsub0, (unsigned)nsub1, m_aos1dMult);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore): Failed to multiply complex arrays on device", m_cuerrLast);
		return 6;
	}
	return 0;
}


int CJFFTCUDAcore::MultiplyC(fcmplx * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	size_t nbytes = nd * sizeof(cufftComplex);
	// prepare to copy source from host to the device
	cufftComplex * pcsrcdev = NULL;
	m_cuerrLast = cudaMalloc((void**)&pcsrcdev, nbytes);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore): Failed to allocate on device", m_cuerrLast);
		return 3;
	}
	// copy from host to device
	m_cuerrLast = cudaMemcpy(pcsrcdev, src, nbytes, ::cudaMemcpyHostToDevice);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore): Failed to copy from host to device", m_cuerrLast);
		return 4;
	}
	// multiply on device pointwise (m_pcw * pcsrvdev)
	m_cuerrLast = ArrayOpMul(m_pcw, pcsrcdev, m_pcw, m_aos1dMult);
	// destroy source helper on device
	m_cuerrLast = cudaFree(pcsrcdev);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore): Failed to destroy helper array on device", m_cuerrLast);
		return 6;
	}
	return 0;
}


/* ----------------------------------------------------------------------------- */
/*                                                                               */
/* >>> DATA INTERFACE <<<                                                        */
/*                                                                               */
/* ----------------------------------------------------------------------------- */

int CJFFTCUDAcore::GetStatus(void)
{
	return m_nstatus;
}

int CJFFTCUDAcore::GetDevice(void)
{
	int idev = -1;
	cudaGetDevice(&idev);
	return idev;
}

int CJFFTCUDAcore::GetDimension(void)
{
	return m_ndim;
}

size_t CJFFTCUDAcore::GetDataSize(void)
{
	size_t nd = 0;
	if (m_nstatus > 0) {
		nd = 1;
		for (int i = 0; i < m_ndim; i++) { // loop over all dimensions
			nd = nd*(size_t)m_pdims[i]; // multiply to total number of items
		}
	}
	return nd;
}

int CJFFTCUDAcore::GetDimSize(int * pdims)
{
	if (m_nstatus == 0) // error: not initialized
		return 1;
	if (pdims == NULL) // error: invalid destination pointer
		return 2;
	memcpy(pdims, m_pdims, (size_t)(sizeof(int)*m_ndim));
	return 0;
}

int CJFFTCUDAcore::GetBlockSize(void)
{
	return m_aos1dMult.nBlockSize;
}


size_t CJFFTCUDAcore::GetFFTMemUsage(void) {
	size_t nResult = 0;
	if (m_nstatus > 0) {
		m_cufftresLast = cufftGetSize(m_plan, &nResult);
		if (m_cufftresLast != CUFFT_SUCCESS) {
			cerr << "Error(JFFTCUDAcore): Failed to determine cufft memory requirement, code: " << m_cufftresLast << endl;
		}
	}
	return nResult;
}

size_t CJFFTCUDAcore::GetMulMemUsage(void) {
	size_t nResult = 0;
	if (m_nstatus > 0) {
		nResult = sizeof(cuComplex) * 2 * GetDataSize();
	}
	return nResult;
}


int CJFFTCUDAcore::SetDataC_d(cuComplex * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::SetDataC_d): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) { // error: invalid source pointer
		cerr << "Error(JFFTCUDAcore::SetDataC_d): Invalid source pointer." << endl;
		return 2;
	}
	// copy from device to device
	size_t nbytes = sizeof(cuComplex)*GetDataSize();
	m_cuerrLast = cudaMemcpy(m_pcw, src, nbytes, ::cudaMemcpyDeviceToDevice);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::SetDataC_d): Failed to copy from device to device", m_cuerrLast);
		return 4;
	}
	return 0;
}


int CJFFTCUDAcore::SetDataC(fcmplx * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::SetDataC): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) { // error: invalid source pointer
		cerr << "Error(JFFTCUDAcore::SetDataC): Invalid source pointer." << endl;
		return 2;
	}
	// copy from host to device
	size_t nbytes = sizeof(cuComplex)*GetDataSize();
	m_cuerrLast = cudaMemcpy(m_pcw, src, nbytes, ::cudaMemcpyHostToDevice);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::SetDataC): Failed to copy from host to device", m_cuerrLast);
		return 4;
	}
	return 0;
}

int CJFFTCUDAcore::SetDataRe(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::SetDataRe): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) { // error: invalid source pointer
		cerr << "Error(JFFTCUDAcore::SetDataRe): Invalid source pointer." << endl;
		return 2;
	}
	size_t nd = GetDataSize();
	size_t nbytes = sizeof(cuComplex)*nd;
	// prepare complex array and set re components
	cuComplex *srctmp = (cuComplex*)malloc(nbytes);
	for (size_t i = 0; i < nd; i++) {
		srctmp[i].x = src[i];
		srctmp[i].y = 0.f;
	}
	// copy from host to device
	m_cuerrLast = cudaMemcpy(m_pcw, srctmp, nbytes, ::cudaMemcpyHostToDevice);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::SetDataRe): Failed to copy from host to device", m_cuerrLast);
		return 4;
	}
	free(srctmp);
	return 0;
}

int CJFFTCUDAcore::SetDataIm(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::SetDataIm): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) { // error: invalid source pointer
		cerr << "Error(JFFTCUDAcore::SetDataIm): Invalid source pointer." << endl;
		return 2;
	}
	size_t nd = GetDataSize();
	size_t nbytes = sizeof(cuComplex)*nd;
	// prepare complex array and set re components
	cuComplex *srctmp = (cuComplex*)malloc(nbytes);
	for (size_t i = 0; i < nd; i++) {
		srctmp[i].x = 0.f;
		srctmp[i].y = src[i];
	}
	// copy from host to device
	m_cuerrLast = cudaMemcpy(m_pcw, srctmp, nbytes, ::cudaMemcpyHostToDevice);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::SetDataIm): Failed to copy from host to device", m_cuerrLast);
		return 4;
	}
	free(srctmp);
	return 0;
}

int CJFFTCUDAcore::SetDataRe_d(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::SetDataRe_d): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) { // error: invalid source pointer
		cerr << "Error(JFFTCUDAcore::SetDataRe_d): Invalid source pointer." << endl;
		return 2;
	}
	m_cuerrLast = ArrayOpSetRe(m_pcw, src, m_aos1dMult);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::SetDataRe_d): Set data real part on device", m_cuerrLast);
		return 3;
	}
	return 0;
}

int CJFFTCUDAcore::SetDataIm_d(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::SetDataIm_d): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) { // error: invalid source pointer
		cerr << "Error(JFFTCUDAcore::SetDataIm_d): Invalid source pointer." << endl;
		return 2;
	}
	m_cuerrLast = ArrayOpSetIm(m_pcw, src, m_aos1dMult);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::SetDataIm_d): Set data imaginary part on device", m_cuerrLast);
		return 3;
	}
	return 0;
}

int CJFFTCUDAcore::GetDataC(fcmplx * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::GetDataC): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) { // error: invalid destination pointer
		cerr << "Error(JFFTCUDAcore::GetDataC): Invalid destination pointer." << endl;
		return 2;
	}
	size_t nbytes = sizeof(fcmplx)*GetDataSize();
	m_cuerrLast = cudaMemcpy(dst, m_pcw, nbytes, ::cudaMemcpyDeviceToHost);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataC): Failed to copy from device to host", m_cuerrLast);
		return 4;
	}
	return 0;
}

int CJFFTCUDAcore::GetDataRe(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::GetDataRe): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) { // error: invalid destination pointer
		cerr << "Error(JFFTCUDAcore::GetDataRe): Invalid destination pointer." << endl;
		return 2;
	}
	size_t nd = GetDataSize();
	size_t nbytes = sizeof(cuComplex)*nd;
	cuComplex *dsttmp = (cuComplex*)malloc(nbytes);
	if (dsttmp == NULL) {
		cerr << "Error(JFFTCUDAcore::GetDataRe): Failed to allocate memory." << endl;
		return 3;
	}
	m_cuerrLast = cudaMemcpy(dsttmp, m_pcw, nbytes, ::cudaMemcpyDeviceToHost);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataRe): Failed to copy from device to host", m_cuerrLast);
		return 4;
	}
	for (size_t i = 0; i < nd; i++) {
		dst[i] = dsttmp[i].x;
	}
	free(dsttmp);
	return 0;
}

int CJFFTCUDAcore::GetDataIm(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::GetDataIm): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) { // error: invalid destination pointer
		cerr << "Error(JFFTCUDAcore::GetDataIm): Invalid destination pointer." << endl;
		return 2;
	}
	size_t nd = GetDataSize();
	size_t nbytes = sizeof(cuComplex)*nd;
	cuComplex *dsttmp = (cuComplex*)malloc(nbytes);
	if (dsttmp == NULL) {
		cerr << "Error(JFFTCUDAcore::GetDataIm): Failed to allocate memory." << endl;
		return 3;
	}
	m_cuerrLast = cudaMemcpy(dsttmp, m_pcw, nbytes, ::cudaMemcpyDeviceToHost);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataIm): Failed to copy from device to host", m_cuerrLast);
		return 4;
	}
	for (size_t i = 0; i < nd; i++) {
		dst[i] = dsttmp[i].y;
	}
	free(dsttmp);
	return 0;
}

int CJFFTCUDAcore::GetDataRe_d(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::GetDataRe_d): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) { // error: invalid source pointer
		cerr << "Error(JFFTCUDAcore::GetDataRe_d): Invalid source pointer." << endl;
		return 2;
	}
	m_cuerrLast = ArrayOpGetRe(src, m_pcw, m_aos1dMult);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataRe_d): Get data real part on device", m_cuerrLast);
		return 3;
	}
	return 0;
}

int CJFFTCUDAcore::GetDataIm_d(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::GetDataIm_d): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) { // error: invalid source pointer
		cerr << "Error(JFFTCUDAcore::GetDataIm_d): Invalid source pointer." << endl;
		return 2;
	}
	m_cuerrLast = ArrayOpGetIm(src, m_pcw, m_aos1dMult);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataIm_d): Set data imaginary part on device", m_cuerrLast);
		return 3;
	}
	return 0;
}

int CJFFTCUDAcore::GetDataAbs(float * dst)
{
	int nerr = 0;
	size_t nd =0;
	size_t nbytes = 0;
	float *d_tmp = NULL;
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::GetDataAbs): Cannot transfer data, not initialized." << endl;
		nerr = 1;
		goto Error;
	}
	if (dst == NULL) { // error: invalid destination pointer
		cerr << "Error(JFFTCUDAcore::GetDataAbs): Invalid destination pointer." << endl;
		nerr = 1;
		goto Error;
	}
	nd = GetDataSize();
	nbytes = sizeof(float)*nd;
	m_cuerrLast = cudaMalloc((void**)&d_tmp, nbytes);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataAbs): Failed to allocate on device", m_cuerrLast);
		nerr = 3;
		goto Error;
	}
	m_cuerrLast = ArrayOpCAbs(d_tmp, m_pcw, m_aos1dMult);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataAbs): Failed to calculate CAbs on device", m_cuerrLast);
		nerr = 4;
		goto Error;
	}
	m_cuerrLast = cudaMemcpy(dst, d_tmp, nbytes, ::cudaMemcpyDeviceToHost);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataAbs): Failed to copy from device to host", m_cuerrLast);
		nerr = 5;
		goto Error;
	}
Error:
	if (NULL != d_tmp) cudaFree(d_tmp);
	return nerr;
}

int CJFFTCUDAcore::GetDataArg(float * dst)
{
	int nerr = 0;
	size_t nd = 0;
	size_t nbytes = 0;
	float *d_tmp = NULL;
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::GetDataArg): Cannot transfer data, not initialized." << endl;
		nerr = 1;
		goto Error;
	}
	if (dst == NULL) { // error: invalid destination pointer
		cerr << "Error(JFFTCUDAcore::GetDataArg): Invalid destination pointer." << endl;
		nerr = 1;
		goto Error;
	}
	nd = GetDataSize();
	nbytes = sizeof(float)*nd;
	m_cuerrLast = cudaMalloc((void**)&d_tmp, nbytes);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataArg): Failed to allocate on device", m_cuerrLast);
		nerr = 3;
		goto Error;
	}
	m_cuerrLast = ArrayOpCArg(d_tmp, m_pcw, m_aos1dMult);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataArg): Failed to calculate CArg on device", m_cuerrLast);
		nerr = 4;
		goto Error;
	}
	m_cuerrLast = cudaMemcpy(dst, d_tmp, nbytes, ::cudaMemcpyDeviceToHost);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataArg): Failed to copy from device to host", m_cuerrLast);
		nerr = 5;
		goto Error;
	}
Error:
	if (NULL != d_tmp) cudaFree(d_tmp);
	return nerr;
}

cuComplex* CJFFTCUDAcore::GetData_d()
{
	return m_pcw;
}

int CJFFTCUDAcore::GetDataC_d(cuComplex * dst)
{
	int nerr = 0;
	size_t nd = 0;
	size_t nbytes = 0;
	float *d_tmp = NULL;
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::GetDataC_d): Cannot transfer data, not initialized." << endl;
		nerr = 1;
		goto Error;
	}
	nd = GetDataSize();
	nbytes = sizeof(float)*nd;
	m_cuerrLast = cudaMemcpy(dst, m_pcw, nbytes, ::cudaMemcpyDeviceToDevice);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataC_d): Failed to copy from device to device", m_cuerrLast);
		nerr = 2;
		goto Error;
	}
Error:
	return nerr;
}

int CJFFTCUDAcore::GetDataPow_d(float * dst)
{
	int nerr = 0;
	size_t nd = 0;
	size_t nbytes = 0;
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::GetDataPow_d): Cannot transfer data, not initialized." << endl;
		nerr = 1;
		goto Error;
	}
	nd = GetDataSize();
	nbytes = sizeof(float)*nd;
	m_cuerrLast = ArrayOpCPow(dst, m_pcw, m_aos1dMult);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataPow_d): Failed calculate complex power on device", m_cuerrLast);
		nerr = 2;
		goto Error;
	}
Error:
	return nerr;
}


int CJFFTCUDAcore::GetDataPow(float * dst)
{
	int nerr = 0;
	size_t nd = 0;
	size_t nbytes = 0;
	float *d_tmp = NULL;
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::GetDataPow_d): Cannot transfer data, not initialized." << endl;
		nerr = 1;
		goto Error;
	}
	nd = GetDataSize();
	nbytes = sizeof(float)*nd;
	m_cuerrLast = cudaMalloc((void**)&d_tmp, nbytes);
	m_cuerrLast = ArrayOpCPow(d_tmp, m_pcw, m_aos1dMult);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataPow_d): Failed calculate complex power on device", m_cuerrLast);
		nerr = 2;
		goto Error;
	}
	m_cuerrLast = cudaMemcpy(dst, d_tmp, nbytes, ::cudaMemcpyDeviceToHost);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataArg): Failed to copy from device to host", m_cuerrLast);
		nerr = 5;
		goto Error;
	}
Error:
	if (NULL != d_tmp) cudaFree(d_tmp);
	return nerr;
}


int CJFFTCUDAcore::GetDataTotalPow(float &pow)
{
	float fpow = 0.f;
	int nerr = 0;
	size_t nd = 0;
	size_t nbytes = 0;
	float *d_tmp = NULL;
	if (m_nstatus < 1) {
		cerr << "Error(JFFTCUDAcore::GetDataArg): Cannot calculate on data, not initialized." << endl;
		nerr = 1;
		goto Error;
	}
	m_cuerrLast = ArrayOpCPowSum(fpow, m_pcw, 1.0f, m_aos1dMult);
	if (m_cuerrLast != cudaSuccess) {
		PostCUDAError("(JFFTCUDAcore::GetDataArg): Failed to calculate CPowSum on device", m_cuerrLast);
		nerr = 2;
		goto Error;
	}
	pow = fpow;
Error:
	return nerr;
}