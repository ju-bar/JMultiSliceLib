//
// C++ source file: JFFTWcore.cpp
// implementation of class CJFFTWcore (declaration see JFFTWcore.h)
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
// class CJFFTWcore implements routines of Fourier transformations using
// libraries libfftw3f-3.lib
//
// The code is supposed to handle one type of Fourier transform only
// for single precision (float) complex data.
//
#include "stdafx.h"
#include "JFFTWcore.h"
//
using namespace std;
//

CJFFTWcore::CJFFTWcore()
{
	m_nstatus = 0;
	m_ndim = 0;
	m_nplanflag = FFTW_ESTIMATE;
	m_pdims = NULL;
	m_pcw = NULL;
	m_planf = NULL;
	m_planb = NULL;
}


CJFFTWcore::~CJFFTWcore()
{
	Deinit();
}


/* ----------------------------------------------------------------------------- */
/*                                                                               */
/* >>> SETUP INTERFACE <<<                                                       */
/*                                                                               */
/* ----------------------------------------------------------------------------- */


void CJFFTWcore::Deinit(void)
{
	if (m_planf != NULL) fftwf_destroy_plan(m_planf);
	if (m_planb != NULL) fftwf_destroy_plan(m_planb);
	if (m_pdims != NULL) free(m_pdims);
	if (m_pcw != NULL) fftwf_free(m_pcw);
	m_planf = NULL;
	m_planb = NULL;
	m_pdims = NULL;
	m_pcw = NULL;
	m_nstatus = 0;
	m_ndim = 0;
	m_nplanflag = FFTW_ESTIMATE;
}

void CJFFTWcore::CleanupFFTW(void)
{
	fftwf_cleanup();
}

int CJFFTWcore::Init(int ndim, int * pdims, int nplanflag)
{
	if (ndim<JFFTWCORE_DIM_MIN) {
		cerr << "Error(JFFTWcore): Cannot initialize, invalid dimension." << endl;
		return 1;
	}
	if (nplanflag != FFTW_MEASURE && nplanflag != FFTW_ESTIMATE &&
		nplanflag != FFTW_WISDOM_ONLY && nplanflag != FFTW_PATIENT) {
		cerr << "Error(JFFTWcore): Cannot initialize, invalid plan flag." << endl;
		return 2;
	}
	int nconsistent = 0;
	int i = 0;
	int nconsistent_target = 2 + ndim; // this is what should exist already
	if (m_nstatus > 0) { // object is initialized // check consistency with init
		if (ndim == m_ndim) nconsistent++; // same dimension ?
		if (nplanflag = m_nplanflag) nconsistent++; // same plan flag
		if (nconsistent == 2) { // same size and flag?
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
		// - new number of complex data items
		int nd = 1;
		for (i = 0; i < m_ndim; i++) { // loop over all dimensions
			m_pdims[i] = pdims[i]; // store dimension length
			nd = nd*m_pdims[i]; // multiply to total number of items
		}
		// - allocate new transformation array
		m_pcw = fftwf_alloc_complex((size_t)nd);
		if (NULL == m_pcw) {
			cerr << "Error(JFFTWcore): Failed to allocate transformation array." << endl;
			return 3;
		}
		// - prepare new transformation plan
		switch (ndim) {
		case 1:
			m_planf = fftwf_plan_dft_1d(m_pdims[0], m_pcw, m_pcw, FFTW_FORWARD, nplanflag);
			m_planb = fftwf_plan_dft_1d(m_pdims[0], m_pcw, m_pcw, FFTW_BACKWARD, nplanflag);
			break;
		case 2:
			m_planf = fftwf_plan_dft_2d(m_pdims[1], m_pdims[0], m_pcw, m_pcw, FFTW_FORWARD, nplanflag);
			m_planb = fftwf_plan_dft_2d(m_pdims[1], m_pdims[0], m_pcw, m_pcw, FFTW_BACKWARD, nplanflag);
			break;
		case 3:
			m_planf = fftwf_plan_dft_3d(m_pdims[2], m_pdims[1], m_pdims[0], m_pcw, m_pcw, FFTW_FORWARD, nplanflag);
			m_planb = fftwf_plan_dft_3d(m_pdims[2], m_pdims[1], m_pdims[0], m_pcw, m_pcw, FFTW_BACKWARD, nplanflag);
			break;
		default: // anything > 3
			m_planf = fftwf_plan_dft(m_ndim, m_pdims, m_pcw, m_pcw, FFTW_FORWARD, nplanflag);
			m_planb = fftwf_plan_dft(m_ndim, m_pdims, m_pcw, m_pcw, FFTW_BACKWARD, nplanflag);
		}
		if (NULL == m_planf || NULL == m_planb) {
			cerr << "Error(JFFTWcore): Failed to initialize plans." << endl;
			return 4;
		}
		m_nplanflag = nplanflag;
		// - update status
		m_nstatus = 1;
		// - re-initialize the transformation array as the plan creation may have used it
		if (0 < Zero()) {
			cerr << "Error(JFFTWcore): Failed to initialize transformation array." << endl;
			return 5;
		}
		
	}
	return 0;
}


/* ----------------------------------------------------------------------------- */
/*                                                                               */
/* >>> OPERATION INTERFACE <<<                                                   */
/*                                                                               */
/* ----------------------------------------------------------------------------- */


int CJFFTWcore::FT(void)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transform, not initialized." << endl;
		return 1;
	}
	fftwf_execute(m_planf);
	return 0;
}


int CJFFTWcore::IFT(void)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transform, not initialized." << endl;
		return 1;
	}
	fftwf_execute(m_planb);
	return 0;
}

int CJFFTWcore::Zero(void)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot zero data, not initialized." << endl;
		return 1;
	}
	size_t nbytes = sizeof(fftwf_complex)*GetDataSize();
	memset(m_pcw,0,nbytes);
	return 0;
}

int CJFFTWcore::Conjugate(void)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot conjugate data, not initialized." << endl;
		return 1;
	}
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		m_pcw[i][1] = -m_pcw[i][1]; // invert the imaginary part
	}
	return 0;
}

int CJFFTWcore::Scale(float sca)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot scale data, not initialized." << endl;
		return 1;
	}
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		m_pcw[i][0] *= sca;
		m_pcw[i][1] *= sca;
	}
	return 0;
}

int CJFFTWcore::MultiplyReal(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	for (size_t i = 0; i<nd; i++) {
		pcdst[i] *= src[i];
		//re1 = m_pcw[i][0]; im1 = m_pcw[i][1];
		//m_pcw[i][0] = re1*src[2 * i] - im1*src[2 * i + 1];
		//m_pcw[i][1] = im1*src[2 * i] + re1*src[2 * i + 1];
	}
	return 0;
}

int CJFFTWcore::MultiplyC(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	fcmplx * pcsrc = reinterpret_cast<fcmplx*>(src);
	for (size_t i = 0; i<nd; i++) {
		pcdst[i] *= pcsrc[i];
		//re1 = m_pcw[i][0]; im1 = m_pcw[i][1];
		//m_pcw[i][0] = re1*src[2 * i] - im1*src[2 * i + 1];
		//m_pcw[i][1] = im1*src[2 * i] + re1*src[2 * i + 1];
	}
	return 0;
}

int CJFFTWcore::MultiplyC(fcmplx * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	fcmplx * pcsrc = src;
	for (size_t i=0; i<nd; i++){
		pcdst[i] *= pcsrc[i];
		//re1 = m_pcw[i][0]; im1 = m_pcw[i][1];
		//m_pcw[i][0] = re1*src[i].real() - im1*src[i].imag();
		//m_pcw[i][1] = im1*src[i].real() + re1*src[i].imag();
	}
	return 0;
}

int CJFFTWcore::MultiplyC(fftwf_complex * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot multiply data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	fcmplx * pcsrc = reinterpret_cast<fcmplx*>(src);
	for (size_t i=0; i<nd; i++){
		pcdst[i] *= pcsrc[i];
		//re1 = m_pcw[i][0]; im1 = m_pcw[i][1];
		//m_pcw[i][0] = re1*src[i][0] - im1*src[i][1];
		//m_pcw[i][1] = im1*src[i][0] + re1*src[i][1];
	}
	return 0;
}


int CJFFTWcore::MultiplySub2dC(float * src, int nsub0, int nsub1)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot multiply data, not initialized." << endl;
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
	fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	fcmplx * pcsrc = reinterpret_cast<fcmplx*>(src);
	for (j1 = 0; j1 < n1; j1++) {
		j2 = j1%nsub1;
		for (i1 = 0; i1 < n0; i1++) {
			i2 = i1%nsub0;
			idx1 = i1 + j1*n0;
			idx2 = i2 + j2*nsub0;
			pcdst[idx1] *= pcsrc[idx2];
		}
	}
	return 0;
}

int CJFFTWcore::MultiplySub2dC(fcmplx * src, int nsub0, int nsub1)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot multiply data, not initialized." << endl;
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
	fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	fcmplx * pcsrc = src;
	for (j1 = 0; j1 < n1; j1++) {
		j2 = j1%nsub1;
		for (i1 = 0; i1 < n0; i1++) {
			i2 = i1%nsub0;
			idx1 = i1 + j1*n0;
			idx2 = i2 + j2*nsub0;
			pcdst[idx1] *= pcsrc[idx2];
		}
	}
	return 0;
}

int CJFFTWcore::MultiplySub2dC(fftwf_complex * src, int nsub0, int nsub1)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot multiply data, not initialized." << endl;
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
	fcmplx * pcdst = reinterpret_cast<fcmplx*>(m_pcw);
	fcmplx * pcsrc = reinterpret_cast<fcmplx*>(src);
	for (j1 = 0; j1 < n1; j1++) {
		j2 = j1%nsub1;
		for (i1 = 0; i1 < n0; i1++) {
			i2 = i1%nsub0;
			idx1 = i1 + j1*n0;
			idx2 = i2 + j2*nsub0;
			pcdst[idx1] *= pcsrc[idx2];
		}
	}
	return 0;
}


/* ----------------------------------------------------------------------------- */
/*                                                                               */
/* >>> DATA INTERFACE <<<                                                        */
/*                                                                               */
/* ----------------------------------------------------------------------------- */

int CJFFTWcore::GetStatus(void)
{
	return m_nstatus;
}

int CJFFTWcore::GetDimension(void)
{
	return m_ndim;
}

size_t CJFFTWcore::GetDataSize(void)
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

int CJFFTWcore::GetDimSize(int * pdims)
{
	if (m_nstatus == 0) // error: not initialized
		return 1;
	if (pdims == NULL) // error: invalid destination pointer
		return 2;
	memcpy(pdims, m_pdims, (size_t)(sizeof(int)*m_ndim));
	return 0;
}

int CJFFTWcore::GetFlag(void)
{
	return m_nplanflag;
}

int CJFFTWcore::SetDataC(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nbytes = sizeof(fftwf_complex)*GetDataSize();
	memcpy(m_pcw, src, nbytes);
	return 0;
}

int CJFFTWcore::SetDataC(fcmplx * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nbytes = sizeof(fftwf_complex)*GetDataSize();
	memcpy(m_pcw, src, nbytes);
	return 0;
}

int CJFFTWcore::SetDataC(fftwf_complex * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nbytes = sizeof(fftwf_complex)*GetDataSize();
	memcpy(m_pcw, src, nbytes);
	return 0;
}

int CJFFTWcore::SetDataRe(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		m_pcw[i][0] = src[i];
		m_pcw[i][1] = 0.f;
	}
	return 0;
}

int CJFFTWcore::SetDataIm(float * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		m_pcw[i][0] = 0.f;
		m_pcw[i][1] = src[i];
	}
	return 0;
}

int CJFFTWcore::SetDataReInt(int * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		m_pcw[i][0] = (float)src[i];
		m_pcw[i][1] = 0.f;
	}
	return 0;
}

int CJFFTWcore::SetDataImInt(int * src)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (src == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		m_pcw[i][0] = 0.f;
		m_pcw[i][1] = (float)src[i];
	}
	return 0;
}

int CJFFTWcore::GetDataC(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid destination pointer
		return 2;
	size_t nbytes = sizeof(fftwf_complex)*GetDataSize();
	memcpy(dst, m_pcw, nbytes);
	return 0;
}

int CJFFTWcore::GetDataC(fcmplx * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid destination pointer
		return 2;
	size_t nbytes = sizeof(fftwf_complex)*GetDataSize();
	memcpy(dst, m_pcw, nbytes);
	return 0;
}

int CJFFTWcore::GetDataC(fftwf_complex * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid destination pointer
		return 2;
	size_t nbytes = sizeof(fftwf_complex)*GetDataSize();
	memcpy(dst, m_pcw, nbytes);
	return 0;
}

int CJFFTWcore::GetDataRe(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		dst[i] = m_pcw[i][0];
	}
	return 0;
}

int CJFFTWcore::GetDataIm(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	for (size_t i = 0; i < nd; i++) {
		dst[i] = m_pcw[i][1];
	}
	return 0;
}

int CJFFTWcore::GetDataAbs(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	float re, im;
	for (size_t i = 0; i < nd; i++) {
		re = m_pcw[i][0];
		im = m_pcw[i][1];
		dst[i] = sqrtf(re*re + im*im);
	}
	return 0;
}

int CJFFTWcore::GetDataArg(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	float re, im;
	for (size_t i = 0; i < nd; i++) {
		re = m_pcw[i][0];
		im = m_pcw[i][1];
		dst[i] = atan2(im,re);
	}
	return 0;
}

int CJFFTWcore::GetDataPow(float * dst)
{
	if (m_nstatus < 1) {
		cerr << "Error(JFFTWcore): Cannot transfer data, not initialized." << endl;
		return 1;
	}
	if (dst == NULL) // error: invalid source pointer
		return 2;
	size_t nd = GetDataSize();
	float re, im;
	for (size_t i = 0; i < nd; i++) {
		re = m_pcw[i][0];
		im = m_pcw[i][1];
		dst[i] = re*re + im*im;
	}
	return 0;
}



float CJFFTWcore::GetDataTotalPow(void)
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
			re = (double)m_pcw[i][0];
			im = (double)m_pcw[i][1];
			dabs = re * re + im * im;
			dy = dabs - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
	}
	return (float)dsum;
}


float CJFFTWcore::GetDataTotalRe(void)
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
			re = (double)m_pcw[i][0];
			dy = re - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
	}
	return (float)dsum;
}

float CJFFTWcore::GetDataTotalIm(void)
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
			im = (double)m_pcw[i][1];
			dy = im - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
	}
	return (float)dsum;
}