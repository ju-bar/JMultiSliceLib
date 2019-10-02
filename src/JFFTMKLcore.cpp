//
// C++ source file: JFFTMKLcore.cpp
// implementation of class CJFFTMKLcore (declaration see JFFTMKLcore.h)
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
// class CJFFTMKLcore implements routines of Fourier transformations using
// the Intel Math Kernel Libraries
//
// The code is supposed to handle one type of Fourier transform only
// for single precision (float) complex data.
//
#include "stdafx.h"
#include "JFFTMKLcore.h"
//
using namespace std;
//

CJFFTMKLcore::CJFFTMKLcore()
{
	m_mkl_status = 0;
	m_nstatus = 0;
	m_ndim = 0;
	m_pdims = NULL;
	m_pcw = NULL;
}


CJFFTMKLcore::~CJFFTMKLcore()
{
	Deinit();
}


/* ----------------------------------------------------------------------------- */
/*                                                                               */
/* >>> SETUP INTERFACE <<<                                                       */
/*                                                                               */
/* ----------------------------------------------------------------------------- */


void CJFFTMKLcore::Deinit(void)
{
	if (m_nstatus > 0) {
		m_mkl_status = DftiFreeDescriptor(&m_descr);
	}
	if (m_pdims != NULL) free(m_pdims);
	if (m_pcw != NULL) free(m_pcw);
	m_pdims = NULL;
	m_pcw = NULL;
	m_nstatus = 0;
	m_ndim = 0;
}

int CJFFTMKLcore::Init(int ndim, int * pdims)
{
	if (ndim<JFFTMKLCORE_DIM_MIN) {
		cerr << "Error(JFFTMKLcore): Cannot initialize, invalid dimension." << endl;
		return 1;
	}
	int nconsistent = 0;
	int i = 0;
	int nconsistent_target = 1 + ndim; // this is what should exist already
	if (m_nstatus > 0) { // object is initialized // check consistency with init
		if (ndim == m_ndim) nconsistent++; // same dimension ?
		if (nconsistent == 1) { // same dimension?
			for (i=0; i < m_ndim; i++) { // check for same dimensions
				if (pdims[i]==m_pdims[i]) nconsistent++;
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
		Deinit();
		// - new dimension
		m_ndim = ndim;
		// - prepare new arrays
		m_pdims = (int*)malloc((size_t)(m_ndim*sizeof(int)));
		int * pinidims = (int*)malloc((size_t)(m_ndim * sizeof(int)));
		// - new number of complex data items
		int nd = 1;
		for (i = 0; i < m_ndim; i++) { // loop over all dimensions
			m_pdims[i] = pdims[i]; // store dimension length
			pinidims[m_ndim - i - 1] = pdims[i]; // initialize reverse dimension list for MKL creator
			nd = nd*m_pdims[i]; // multiply to total number of items
		}
		// - allocate new transformation array
		m_pcw = (fcmplx*)malloc((size_t)nd*sizeof(fcmplx));
		if (NULL == m_pcw) {
			cerr << "Error(JFFTMKLcore): Failed to allocate transformation array." << endl;
			if (NULL != pinidims) { free(pinidims); }
			return 3;
		}
		// - prepare new transformation plan
		m_mkl_status = DftiCreateDescriptor(&m_descr, DFTI_SINGLE, DFTI_COMPLEX, ndim, pinidims);
		if (m_mkl_status > 0 && !DftiErrorClass(m_mkl_status, DFTI_NO_ERROR)) {
			cerr << "Error(JFFTMKLcore): Failed to create fft descriptor." << endl;
			cerr << "-> MKL: " << DftiErrorMessage(m_mkl_status) << endl;
			if (NULL != pinidims) { free(pinidims); }
			return 4;
		}
		//
		m_mkl_status = DftiCommitDescriptor(m_descr);
		if (m_mkl_status > 0 && !DftiErrorClass(m_mkl_status, DFTI_NO_ERROR)) {
			cerr << "Error(JFFTMKLcore): Failed to commit fft descriptor." << endl;
			cerr << "-> MKL: " << DftiErrorMessage(m_mkl_status) << endl;
			if (NULL != pinidims) { free(pinidims); }
			return 5;
		}
		// - update status
		m_nstatus = 1;
		if (NULL != pinidims) { free(pinidims); }
		// - re-initialize the transformation array as the plan creation may have used it
		if (0 < Zero()) {
			cerr << "Error(JFFTMKLcore): Failed to initialize transformation array." << endl;
			return 6;
		}
		
	}
	return 0;
}


/* ----------------------------------------------------------------------------- */
/*                                                                               */
/* >>> OPERATION INTERFACE <<<                                                   */
/*                                                                               */
/* ----------------------------------------------------------------------------- */


int CJFFTMKLcore::FT(void)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transform, not initialized." << endl;
		return 1;
	}
	m_mkl_status = DftiComputeForward(m_descr, m_pcw);
	if (m_mkl_status > 0 && !DftiErrorClass(m_mkl_status, DFTI_NO_ERROR))
	{
		cerr << "Error(JFFTMKLcore): Failed to compute forward fft." << endl;
		cerr << "-> MKL: " << DftiErrorMessage(m_mkl_status) << endl;
		return 2;
	}
	return 0;
}

int CJFFTMKLcore::FT_h(fcmplx* src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transform, not initialized." << endl;
		return 1;
	}
	if (NULL == src) {
		cerr << "Error(JFFTMKLcore): Cannot transform, invalid data pointer." << endl;
		return 2;
	}
	m_mkl_status = DftiComputeForward(m_descr, src);
	if (m_mkl_status > 0 && !DftiErrorClass(m_mkl_status, DFTI_NO_ERROR))
	{
		cerr << "Error(JFFTMKLcore): Failed to compute forward fft." << endl;
		cerr << "-> MKL: " << DftiErrorMessage(m_mkl_status) << endl;
		return 3;
	}
	return 0;
}


int CJFFTMKLcore::IFT(void)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transform, not initialized." << endl;
		return 1;
	}
	m_mkl_status = DftiComputeBackward(m_descr, m_pcw);
	if (m_mkl_status > 0 && !DftiErrorClass(m_mkl_status, DFTI_NO_ERROR))
	{
		cerr << "Error(JFFTMKLcore): Failed to compute backward fft." << endl;
		cerr << "-> MKL: " << DftiErrorMessage(m_mkl_status) << endl;
		return 2;
	}
	return 0;
}

int CJFFTMKLcore::IFT_h(fcmplx* src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transform, not initialized." << endl;
		return 1;
	}
	if (NULL == src) {
		cerr << "Error(JFFTMKLcore): Cannot transform, invalid data pointer." << endl;
		return 2;
	}
	m_mkl_status = DftiComputeBackward(m_descr, src);
	if (m_mkl_status > 0 && !DftiErrorClass(m_mkl_status, DFTI_NO_ERROR))
	{
		cerr << "Error(JFFTMKLcore): Failed to compute backward fft." << endl;
		cerr << "-> MKL: " << DftiErrorMessage(m_mkl_status) << endl;
		return 3;
	}
	return 0;
}

int CJFFTMKLcore::Zero(void)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot zero data, not initialized." << endl;
		return 1;
	}
	size_t nbytes = sizeof(fcmplx)*GetDataSize();
	memset(m_pcw,0,nbytes);
	return 0;
}

int CJFFTMKLcore::Conjugate(void)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot conjugate data, not initialized." << endl;
		return 1;
	}
	size_t nd = GetDataSize();
	float iv = 0.f;
	for (size_t i = 0; i < nd; i++) {
		iv = -m_pcw[i].imag();
		m_pcw[i].imag(iv); // invert the imaginary part
	}
	return 0;
}

int CJFFTMKLcore::Scale(float sca)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot scale data, not initialized." << endl;
		return 1;
	}
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		m_pcw[i] *= sca;
		//m_pcw[i][0] *= sca;
		//m_pcw[i][1] *= sca;
	}
	return 0;
}

int CJFFTMKLcore::MultiplyReal(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	//fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	for (size_t i = 0; i<nd; i++) {
		m_pcw[i] *= src[i];
		//re1 = m_pcw[i][0]; im1 = m_pcw[i][1];
		//m_pcw[i][0] = re1*src[2 * i] - im1*src[2 * i + 1];
		//m_pcw[i][1] = im1*src[2 * i] + re1*src[2 * i + 1];
	}
	return 0;
}

int CJFFTMKLcore::MultiplyC(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	//fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	fcmplx * pcsrc = reinterpret_cast<fcmplx*>(src);
	for (size_t i = 0; i<nd; i++) {
		m_pcw[i] *= pcsrc[i];
		//re1 = m_pcw[i][0]; im1 = m_pcw[i][1];
		//m_pcw[i][0] = re1*src[2 * i] - im1*src[2 * i + 1];
		//m_pcw[i][1] = im1*src[2 * i] + re1*src[2 * i + 1];
	}
	return 0;
}

int CJFFTMKLcore::MultiplyC(fcmplx * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	//fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	//fcmplx * pcsrc = src;
	for (size_t i=0; i<nd; i++){
		m_pcw[i] *= src[i];
		//re1 = m_pcw[i][0]; im1 = m_pcw[i][1];
		//m_pcw[i][0] = re1*src[i].real() - im1*src[i].imag();
		//m_pcw[i][1] = im1*src[i].real() + re1*src[i].imag();
	}
	return 0;
}

int CJFFTMKLcore::MultiplySub2dC(float * src, int nsub0, int nsub1)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	if (m_ndim != 2) // error: invalid core number of dimensions
		return 3;
	int idx1 = 0, idx2 = 0;
	int n0 = m_pdims[0];
	int n1 = m_pdims[1];
	int i1 = 0, j1 = 0;
	int i2 = 0, j2 = 0;
	//fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	fcmplx * pcsrc = reinterpret_cast<fcmplx*>(src);
	for (j1 = 0; j1 < n1; j1++) {
		j2 = j1%nsub1;
		for (i1 = 0; i1 < n0; i1++) {
			i2 = i1%nsub0;
			idx1 = i1 + j1*n0;
			idx2 = i2 + j2*nsub0;
			m_pcw[idx1] *= pcsrc[idx2];
		}
	}
	return 0;
}

int CJFFTMKLcore::MultiplySub2dC(fcmplx * src, int nsub0, int nsub1)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	if (m_ndim != 2) // error: invalid core number of dimensions
		return 3;
	int idx1 = 0, idx2 = 0;
	int n0 = m_pdims[0];
	int n1 = m_pdims[1];
	int i1 = 0, j1 = 0;
	int i2 = 0, j2 = 0;
	//fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	fcmplx * pcsrc = src;
	for (j1 = 0; j1 < n1; j1++) {
		j2 = j1%nsub1;
		for (i1 = 0; i1 < n0; i1++) {
			i2 = i1%nsub0;
			idx1 = i1 + j1*n0;
			idx2 = i2 + j2*nsub0;
			m_pcw[idx1] *= pcsrc[idx2];
		}
	}
	return 0;
}


/* ----------------------------------------------------------------------------- */
/*                                                                               */
/* >>> DATA INTERFACE <<<                                                        */
/*                                                                               */
/* ----------------------------------------------------------------------------- */

int CJFFTMKLcore::GetStatus(void)
{
	return m_nstatus;
}

int CJFFTMKLcore::GetDimension(void)
{
	return m_ndim;
}

size_t CJFFTMKLcore::GetDataSize(void)
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

int CJFFTMKLcore::GetDimSize(int * pdims)
{
	if (m_nstatus == 0) // error: not initialized
		return 1;
	if (pdims == NULL) // error: invalid destination pointer
		return 2;
	memcpy(pdims, m_pdims, (size_t)(sizeof(int)*m_ndim));
	return 0;
}

int CJFFTMKLcore::SetDataC(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nbytes = sizeof(fcmplx)*GetDataSize();
	memcpy(m_pcw, src, nbytes);
	return 0;
}

int CJFFTMKLcore::SetDataC(fcmplx * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nbytes = sizeof(fcmplx)*GetDataSize();
	memcpy(m_pcw, src, nbytes);
	return 0;
}

int CJFFTMKLcore::SetDataRe(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		m_pcw[i].real(src[i]);
		m_pcw[i].imag( 0.f );
	}
	return 0;
}

int CJFFTMKLcore::SetDataIm(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		m_pcw[i].real(0.f);
		m_pcw[i].imag(src[i]);
	}
	return 0;
}

int CJFFTMKLcore::SetDataReInt(int * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		m_pcw[i].real((float)src[i]);
		m_pcw[i].imag(0.f);
	}
	return 0;
}

int CJFFTMKLcore::SetDataImInt(int * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		m_pcw[i].real(0.f);
		m_pcw[i].imag((float)src[i]);
	}
	return 0;
}

int CJFFTMKLcore::GetDataC(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid destination pointer
		return 2;
	size_t nbytes = sizeof(fcmplx)*GetDataSize();
	memcpy(dst, m_pcw, nbytes);
	return 0;
}

int CJFFTMKLcore::GetDataC(fcmplx * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid destination pointer
		return 2;
	size_t nbytes = sizeof(fcmplx)*GetDataSize();
	memcpy(dst, m_pcw, nbytes);
	return 0;
}

int CJFFTMKLcore::GetDataRe(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		dst[i] = m_pcw[i].real();
	}
	return 0;
}

int CJFFTMKLcore::GetDataIm(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		dst[i] = m_pcw[i].imag();
	}
	return 0;
}

int CJFFTMKLcore::GetDataAbs(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	float re, im;
	for (size_t i = 0; i < nd; i++) {
		re = m_pcw[i].real();
		im = m_pcw[i].imag();
		dst[i] = sqrtf(re*re + im*im);
	}
	return 0;
}

int CJFFTMKLcore::GetDataArg(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	float re, im;
	for (size_t i = 0; i < nd; i++) {
		re = m_pcw[i].real();
		im = m_pcw[i].imag();
		dst[i] = atan2(im,re);
	}
	return 0;
}

int CJFFTMKLcore::GetDataPow(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTMKLcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	float re, im;
	for (size_t i = 0; i < nd; i++) {
		re = m_pcw[i].real();
		im = m_pcw[i].imag();
		dst[i] = re*re + im*im;
	}
	return 0;
}



float CJFFTMKLcore::GetDataTotalPow(void)
{
	/*       Performs a Kahan sum on the absolute square of the data              */
	/*       https://en.wikipedia.org/wiki/Kahan_summation_algorithm              */
	double dsum = 0.0;
	double dabs = 0.0;
	double dc = 0.0;
	double dy = 0.0;
	double dt = 0.0;
	double re = 0.0;
	double im = 0.0;
	size_t len = GetDataSize();
	if (len > 0 && m_pcw != NULL && m_nstatus > 0 ) {
		for (size_t i = 0; i < len; i++) {
			re = (double)m_pcw[i].real();
			im = (double)m_pcw[i].imag();
			dabs = re * re + im * im;
			dy = dabs - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
	}
	return (float)dsum;
}


float CJFFTMKLcore::GetDataTotalRe(void)
{
	/*       Performs a Kahan sum on the data real part                           */
	/*       https://en.wikipedia.org/wiki/Kahan_summation_algorithm              */
	double dsum = 0.0;
	double dc = 0.0;
	double dy = 0.0;
	double dt = 0.0;
	double re = 0.0;
	size_t len = GetDataSize();
	if (len > 0 && m_pcw != NULL && m_nstatus > 0) {
		for (size_t i = 0; i < len; i++) {
			re = (double)m_pcw[i].real();
			dy = re - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
	}
	return (float)dsum;
}

float CJFFTMKLcore::GetDataTotalIm(void)
{
	/*       Performs a Kahan sum on the data imaginary part                      */
	/*       https://en.wikipedia.org/wiki/Kahan_summation_algorithm              */
	double dsum = 0.0;
	double dc = 0.0;
	double dy = 0.0;
	double dt = 0.0;
	double im = 0.0;
	size_t len = GetDataSize();
	if (len > 0 && m_pcw != NULL && m_nstatus > 0) {
		for (size_t i = 0; i < len; i++) {
			im = (double)m_pcw[i].imag();
			dy = im - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
	}
	return (float)dsum;
}