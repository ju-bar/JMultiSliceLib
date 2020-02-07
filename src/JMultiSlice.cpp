//
// C++ source file: JMultiSlice.cpp
// implementation for library JMultislice.lib (declarations see JMultislice.h)
//
//
// Copyright (C) 2018, 2019 - Juri Barthel (juribarthel@gmail.com)
// Copyright (C) 2018, 2019 - Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
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

//#include "stdafx.h"
#include "JMultiSlice.h"
#include "NatureConstants.h"
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <thread>
/*
#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

using namespace std;
*/

inline int mod(int a, int b) { int ret = a % b; return ret >= 0 ? ret : ret + b; }

void ffstrsum(float *a, size_t n, float *s)
{
	*s = 0.f;
	float ftmp = 0.f;
	if (n > 0) {
		for (size_t i = 0; i < n; i++) {
			ftmp += a[i];
		}
		*s = ftmp;
	}
	return;
}

void fdstrsum(float *a, size_t n, float *s)
{
	*s = 0.f;
	double dtmp = 0.0;
	if (n > 0) {
		for (size_t i = 0; i < n; i++) {
			dtmp += (double)a[i];
		}
		*s = (float)dtmp;
	}
	return;
}


void fkahan(float *a, size_t n, float *s)
{
	*s = 0.f;
	if (n > 0) {
		float c = 0.f, t = 0.f, y = 0.f;
		for (size_t i = 0; i < n; i++) {
			y = a[i] - c;
			t = *s + y; // temp result
			c = (t - *s) - y; // new error
			*s = t; // updated sum
		}
	}
}


void fdncs2m(float *a, size_t n, float *s)
{
	// roll out to 5
	if (n <= 0) { *s = 0.; return; }
	if (n == 1) { *s = a[0]; return; }
	if (n == 2) { *s = a[0] + a[1]; return; }
	if (n < (int)_JMS_SUMMATION_BTF_THR) { ffstrsum(a, n, s); return; }
	// recurse shift 2
	*s = 0.f;
	size_t n0 = 0, n1 = 1, n2 = 2, nc = 0, idx = 0, itmp = 0;
	// calculate number of strides through the buffer
	size_t ntmp = (size_t)ceil((double)n / (double)_JMS_SUMMATION_BUFFER);
	float r = 0.0f, t = 0.0f;
	float* dst = s; // preset destination with single result number
	if (ntmp > 1) { // there will be more than one stride -> more than one number
		dst = (float*)calloc(ntmp, sizeof(float)); // allocate new destination buffer
	}
	else { // butterfly on "a" directly, using s as output target
		r = 0.0f;
		n2 = 2; n1 = 1;
		while (n2 <= n) {
			for (idx = n2 - 1; idx < n; idx += n2) {
				a[idx] += a[idx - n1];
			}
			if (n1 <= n % n2) { // handle left-over component
				r += a[idx - n1];
			}
			n1 = n2;
			n2 = n2 << 1;
		}
		*s += (a[n1 - 1] + r);
		return;
	}
	float* tmp = a; // pre-link to beginning of input temp buffer
	while (n0 < n) { // loop over strides, using dst as output target
		nc = __min(_JMS_SUMMATION_BUFFER, n - n0); // number of copied items
		if (nc == (int)_JMS_SUMMATION_BUFFER) { // butterfly on full 2^M buffer length. This is repeated ntmp-1 times.
			tmp = &a[n0]; // link to offset in a
			n2 = 2; n1 = 1;
			while (n2 < (int)_JMS_SUMMATION_BUFFER) {
				for (idx = n2 - 1; idx < (int)_JMS_SUMMATION_BUFFER; idx += n2) {
					tmp[idx] += tmp[idx - n1];
				}
				n1 = n2;
				n2 = n2 << 1;
			}
			dst[itmp] += (tmp[n1 - 1] + tmp[n2 - 1]); // store intermediate result in stride slot of destination
		}
		else { // butterfly on remaining buffer length (not power of two), this happens only once!
			// roll-out to 5
			if (1 == nc) { dst[itmp] = a[n0]; }
			else if (2 == nc) { dst[itmp] = a[n0] + a[1 + n0]; }
			else if (_JMS_SUMMATION_BTF_THR > nc) { ffstrsum(&a[n0], nc, &dst[itmp]); }
			else if (_JMS_SUMMATION_BTF_THR <= nc) {
				tmp = &a[n0]; // link to offset in a
				n2 = 2; n1 = 1;
				while (n2 <= nc) {
					for (idx = n2 - 1; idx < nc; idx += n2) {
						tmp[idx] += tmp[idx - n1];
					}
					if (n1 <= nc % n2) { // handle left-over component
						r += tmp[idx - n1];
					}
					n1 = n2;
					n2 = n2 << 1;
				}
				dst[itmp] += (tmp[n1 - 1] + r);
			}
		}
		n0 += (int)_JMS_SUMMATION_BUFFER;
		itmp++;
	}
	if (ntmp > 1) { // recurse dst if more than one buffer stride happened
		fdncs2m(dst, ntmp, s); // sum on dst buffer to s
		free(dst); // release dst buffer memory
	}
	return; // s received the result
}


void fdncs2(float *a, size_t n, float *s)
{
	// roll out to 5
	if (n <= 0) { *s = 0.; return; }
	if (n == 1) { *s = a[0]; return; }
	if (n == 2) { *s = a[0] + a[1]; return; }
	if (n < (int)_JMS_SUMMATION_BTF_THR) { ffstrsum(a, n, s); return; }
	// recurse shift 2
	*s = 0.f;
	size_t n0 = 0, n1 = 1, n2 = 2, nc = 0, idx = 0, itmp = 0;
	// calculate number of strides through the buffer
	size_t ntmp = (size_t)ceil((double)n / (double)_JMS_SUMMATION_BUFFER);
	size_t nbb = sizeof(float)*_JMS_SUMMATION_BUFFER;
	size_t nbc = 0;
	float r = 0.0f, t = 0.0f;
	float* tmp = (float*)malloc(nbb); // alloc working buffer
	float* dst = s; // preset destination with single result number
	if (ntmp > 1) { // there will be more than one stride -> more than one number
		dst = (float*)calloc(ntmp, sizeof(float)); // allocate new destination buffer
	}
	else { // small n -> butterfly on tmp (copy of a), using s as output target
		memcpy(tmp, a, sizeof(float)*n); // prepare summation buffer
		r = 0.0f;
		n2 = 2; n1 = 1;
		while (n2 <= n) {
			for (idx = n2 - 1; idx < n; idx += n2) {
				tmp[idx] += tmp[idx - n1];
			}
			if (n1 <= n % n2) { // handle left-over component
				r += tmp[idx - n1];
			}
			n1 = n2;
			n2 = n2 << 1;
		}
		*s += (tmp[n1 - 1] + r);
		free(tmp);
		return;
	}
	while (n0 < n) { // loop over strides, using dst as output target
		nc = __min(_JMS_SUMMATION_BUFFER, n - n0); // number of copied items
		if (nc == (int)_JMS_SUMMATION_BUFFER) { // butterfly on full 2^M buffer length. This is repeated ntmp-1 times.
			memcpy(tmp, &a[n0], nbb); // prepare summation buffer
			n2 = 2; n1 = 1;
			while (n2 < (int)_JMS_SUMMATION_BUFFER) {
				for (idx = n2 - 1; idx < (int)_JMS_SUMMATION_BUFFER; idx += n2) {
					tmp[idx] += tmp[idx - n1];
				}
				n1 = n2;
				n2 = n2 << 1;
			}
			dst[itmp] += (tmp[n1 - 1] + tmp[n2 - 1]); // store intermediate result in stride slot of destination
		}
		else { // butterfly on remaining buffer length (not power of two), this happens only once!
			// roll-out to 5
			if (1 == nc) { dst[itmp] = a[n0]; }
			else if (2 == nc) { dst[itmp] = a[n0] + a[1 + n0]; }
			else if (_JMS_SUMMATION_BTF_THR > nc) { ffstrsum(&a[n0], nc, &dst[itmp]); }
			else if (_JMS_SUMMATION_BTF_THR <= nc) {
				memcpy(tmp, &a[n0], sizeof(float)*nc); // prepare summation buffer
				n2 = 2; n1 = 1;
				while (n2 <= nc) {
					for (idx = n2 - 1; idx < nc; idx += n2) {
						tmp[idx] += tmp[idx - n1];
					}
					if (n1 <= nc % n2) { // handle left-over component
						r += tmp[idx - n1];
					}
					n1 = n2;
					n2 = n2 << 1;
				}
				dst[itmp] += (tmp[n1 - 1] + r);
			}
		}
		n0 += (int)_JMS_SUMMATION_BUFFER;
		itmp++;
	}
	free(tmp); // release tmp buffer memory
	if (ntmp > 1) { // recurse dst if more than one buffer stride happened
		fdncs2m(dst, ntmp, s); // sum on dst buffer to s
		free(dst); // release dst buffer memory
	}
	return; // s received the result
}


CJMultiSlice::CJMultiSlice()
{
	// initialize the class members
	m_dbg = 0;
	m_ht = 0.f;
	m_wl = 0.f;
	m_a0[0] = 0.f;
	m_a0[1] = 0.f;
	m_a0[2] = 0.f;
	m_nscslc = 0;
	m_nslcvar = 0;
	m_nobjslc = 0;
	m_nscx = 0;
	m_nscy = 0;
	m_npgx = 0;
	m_npgy = 0;
	m_npro = 0;
	m_nprotype = 0;
	m_ndet = 0;
	m_ndetper = 0;
	m_imagedet = (int)_JMS_DETECT_NONE;
	m_ndetslc = 0;
	m_status_setup_CPU = 0;
	m_status_setup_GPU = 0;
	m_threads_CPU_out = 0;
	m_ncputhreads = 0;
	m_status_calc_CPU = NULL;
	m_status_calc_GPU = 0;
	m_objslc = NULL;
	m_det_objslc = NULL;
	m_nvarslc = NULL;
	m_slcoffset = NULL;
	m_slcthick = NULL;
	m_proidx = NULL;
	m_detmask_len = NULL;
	m_h_fnx = NULL;
	m_h_fny = NULL;
	m_dif_descan_flg = 0;
	m_plasmonmc_flg = 0;
	m_d_knx = NULL;
	m_d_kny = NULL;
	m_h_wav = NULL;
	m_h_wav0 = NULL;
	m_h_pgr = NULL;
	m_h_pro = NULL;
	m_h_det = NULL;
	m_h_detmask = NULL;
	m_h_det_int = NULL;
	m_h_det_img = NULL;
	m_h_det_dif = NULL;
	m_h_det_wfr = NULL;
	m_h_det_wff = NULL;
	m_h_dif_ndescanx = NULL;
	m_h_dif_ndescany = NULL;
	m_d_wav = NULL;
	m_d_wav0 = NULL;
	m_d_pgr = NULL;
	m_d_pro = NULL;
	m_d_det = NULL;
	m_d_detmask = NULL;
	m_d_det_int = NULL;
	m_d_det_img = NULL;
	m_d_det_dif = NULL;
	m_d_det_wfr = NULL;
	m_d_det_wff = NULL;
	m_d_det_tmp = NULL;
	m_d_det_tmpwav = NULL;
	m_d_dif_ndescanx = 0;
	m_d_dif_ndescany = 0;
	m_d_pgr_src = 0;
	m_jcpuco = NULL;
	m_prng = &m_lrng;
}

CJMultiSlice::~CJMultiSlice()
{
	// clean up to be safe
	Cleanup();
}

int CJMultiSlice::SetDebugLevel(int dbg)
{
	int olddbg = m_dbg;
	m_dbg = dbg;
	return olddbg;
}

void CJMultiSlice::SetRng(CRng *prng) 
{
	if (prng) {
		m_prng = prng;
	}
	else {
		m_prng = &m_lrng;
	}
}

CRng* CJMultiSlice::GetRng()
{
	return m_prng;
}

void CJMultiSlice::SeedRngEx(int nseed)
{
	if (m_prng) {
		m_prng->seed(nseed);
	}
	else {
		m_lrng.seed(nseed);
	}
}

float CJMultiSlice::SetHighTension(float ht)
{
	float oldht = m_ht;
	m_ht = ht;
	m_wl = (float)(_WLELKEV / sqrt(ht*(ht + 2.*_EEL0KEV)));
	return oldht;
}

float CJMultiSlice::GetHighTension(void)
{
	return m_ht;
}

float CJMultiSlice::GetWaveLength(void)
{
	return m_wl;
}

void CJMultiSlice::SetSupercellSize(float *a0)
{
	if (NULL != a0) {
		m_a0[0] = a0[0];
		m_a0[1] = a0[1];
		m_a0[2] = a0[2];
	}
}

void CJMultiSlice::SetSupercellSize(float a, float b, float c)
{
	m_a0[0] = a;
	m_a0[1] = b;
	m_a0[2] = c;
}


void CJMultiSlice::SetGridSize(int nx, int ny)
{
	m_nscx = nx;
	m_nscy = ny;
}

int CJMultiSlice::GetDetetionSliceNum(void)
{
	return m_ndetslc;
}

int CJMultiSlice::GetDetNum(void)
{
	return m_ndet;
}

void CJMultiSlice::GetGridSize(int &nx, int &ny)
{
	nx = m_nscx;
	ny = m_nscy;
}

void CJMultiSlice::GetSupercellSize(float* a0)
{
	if (NULL != a0) {
		a0[0] = m_a0[0];
		a0[1] = m_a0[1];
		a0[2] = m_a0[2];
	}
}

void CJMultiSlice::GetSupercellSize(float &a, float &b, float &c)
{
	a = m_a0[0];
	b = m_a0[1];
	c = m_a0[2];
}



int CJMultiSlice::PhaseGratingSetup(int whichcode, int nx, int ny, int nslc, int nvarmax, int* nslcvar)
{
	// Assumes that SetGridSize has been set before.
	// Allocates m_nvarslc, m_h_pgr, m_d_pgr, m_slcoffset
	int nerr = 0;
	int i = 0;
	int64_t memoffset = 0;
	int64_t nitems = (int64_t)nx*ny;
	m_npgx = 0;
	m_npgy = 0;
	m_nscslc = 0;
	m_nslcvar = 1;
	
	// cleanup old setups, they are not valid anymore
	DeallocMem_h((void**)&m_nvarslc);
	DeallocMem_h((void**)&m_slcthick);
	DeallocMem_h((void**)&m_h_pgr);
	DeallocMem_h((void**)&m_slcoffset);
	DeallocMem_d((void**)&m_d_pgr);
	if (nslc < 1 || NULL == nslcvar) goto _Exit; // no phase gratings

	// setup of phasegrating data common to both codes
	m_npgx = nx;
	m_npgy = ny;
	m_nscslc = nslc; // number of phase grating types
	m_nslcvar = nvarmax; // max. number of variants per phase grating
	if (0<AllocMem_h((void**)&m_nvarslc, sizeof(int)*m_nscslc, "PhaseGratingSetup", "#variants per slice")) { nerr = 1; goto _Exit; }
	memcpy(m_nvarslc, nslcvar, sizeof(int)*m_nscslc); // copy number of variants
	if (0<AllocMem_h((void**)&m_slcthick, sizeof(float)*m_nscslc, "PhaseGratingSetup", "slice thicknesses", true)) { nerr = 2; goto _Exit; }
	if (0<AllocMem_h((void**)&m_slcoffset, sizeof(int64_t)*m_nscslc, "PhaseGratingSetup", "slice offset on GPU", true)) { nerr = 3; goto _Exit; }
	for (i = 0; i < m_nscslc; i++) { // do helper array presets for each slice
		if (m_nvarslc[i] > m_nslcvar) m_nvarslc[i] = m_nslcvar; // limit the number of variants to the maximum set by m_nslcvar
		m_slcoffset[i] = memoffset; // set the access pointers to slice data in device memory (number of array item offset)
		memoffset += ((int64_t)m_nvarslc[i]*nitems); // set offset for next slice data chunk
	}
	// Note: memoffset holds now the total number of phase grating items to be allocated

	// CPU code - phase grating setup -> Memory prepared but no data set.
	if (whichcode & (int)_JMS_CODE_CPU) {
		if (0<AllocMem_h((void**)&m_h_pgr, sizeof(fcmplx*)*m_nscslc, "PhaseGratingSetup", "phase grating addresses", true)) { nerr = 4; goto _Exit; }
	}

	// GPU code - phase grating setup -> Memory prepared but no data set.
	if (whichcode & (int)_JMS_CODE_GPU) {
		if (memoffset > 0) { // make sure to allocate only if data is expected.
			if (0 == m_d_pgr_src) { // try allocation for pre-loading of all phase gratings to device
				if (0 < AllocMem_d((void**)&m_d_pgr, sizeof(cuComplex)*memoffset, "PhaseGratingSetup", "phase gratings", true)) {
					m_d_pgr_src = 1; // continue trying to use load-on-demand
				}
			}
			if (1 == m_d_pgr_src) { // setup load-on-demand to device
				if (0 < AllocMem_d((void**)&m_d_pgr, sizeof(cuComplex)*nitems, "PhaseGratingSetup", "single phase grating", true)) {
					nerr = 104; goto _Exit; // cannot use GPU
				}
				if (NULL == m_h_pgr) { // single phase grating mode needs phase grating pointers on host
					// allocate host pointer list // allocate this list here, make sure to check setup in SetPhaseGratingData
					if (0 < AllocMem_h((void**)&m_h_pgr, sizeof(fcmplx*)*m_nscslc, "PhaseGratingSetup", "phase grating addresses", true)) { nerr = 105; goto _Exit; }
				}
			}
		}
	}

_Exit:
	if (nerr == 0) {
		if (whichcode & (int)_JMS_CODE_CPU) m_status_setup_CPU |= (int)_JMS_STATUS_PGR;
		if (whichcode & (int)_JMS_CODE_GPU) m_status_setup_GPU |= (int)_JMS_STATUS_PGR;
	}
	return nerr;
}



int CJMultiSlice::ObjectSliceSetup(int nobjslc, int* objslc)
{
	// No cross-check is done whether the structure slice IDs in objslc are valid
	int nerr = 0;
	m_nobjslc = 0;
	DeallocMem_h((void**)&m_objslc);
	if (nobjslc < 1 || NULL == objslc) goto _Exit; // no object slices

	m_nobjslc = nobjslc; // max. number of object slices
	if (0 < AllocMem_h((void**)&m_objslc, sizeof(int)*m_nobjslc, "ObjectSliceSetup", "object slice IDs")) { nerr = 1; goto _Exit; }
	memcpy(m_objslc, objslc, sizeof(int)*m_nobjslc); // copy slice IDs
_Exit:
	if (nerr == 0) {
		m_status_setup_CPU |= (int)_JMS_STATUS_OBJ;
		m_status_setup_GPU |= (int)_JMS_STATUS_OBJ;
	}
	return nerr;
}


int CJMultiSlice::PropagatorSetup(int whichcode, int npro, int *proidx)
{
	// No cross-check is done whether the structure slice IDs in proidx are valid
	// the length of array proidx is assumed to be m_nscslc
	int nerr = 0;
	int64_t nproitems = 0;
	int64_t nitems = (int64_t)m_nscx*m_nscy;
	m_npro = 0;
	if (npro < 1 || NULL == proidx) goto _Exit;

	// Setup relevant for both codes
	m_npro = npro; // number of propagator functions
	if (0 < AllocMem_h((void**)&m_proidx, sizeof(int)*m_nscslc, "PropagatorSetup", "propagator access indices")) { nerr = 1; goto _Exit; }
	memcpy(m_proidx, proidx, sizeof(int)*m_nscslc); // copy propagator access indices

	// CPU code - propagator setup -> Memory prepared but no data set.
	if (whichcode & (int)_JMS_CODE_CPU) {
		if (0 < AllocMem_h((void**)&m_h_pro, sizeof(fcmplx*)*m_npro, "PropagatorSetup", "propagator address list", true)) { nerr = 2; goto _Exit; }
	}

	// GPU code - propagator setup -> Memory prepared but no data set.
	if (whichcode & (int)_JMS_CODE_GPU) {
		nproitems = (int64_t)m_npro*nitems;
		if (nproitems > 0) { // make sure to allocate only if data is expected.
			if (0 < AllocMem_d((void**)&m_d_pro, sizeof(cuComplex*)*nproitems, "PropagatorSetup", "propagators")) { nerr = 102; goto _Exit; }
		}
	}

_Exit:
	if (nerr == 0) {
		if (whichcode & (int)_JMS_CODE_CPU) m_status_setup_CPU |= (int)_JMS_STATUS_PRO;
		if (whichcode & (int)_JMS_CODE_GPU) m_status_setup_GPU |= (int)_JMS_STATUS_PRO;
	}
	return nerr;
}


int CJMultiSlice::SetDetectionPlanes(int ndetper, int nobjslc, int * det_objslc)
{
	int islc = 0;
	int ndetslc = 0;
	int noslc = __max(0, nobjslc);
	if (det_objslc && nobjslc > 0) {
		for (islc = 0; islc <= nobjslc; islc++) { // preset all slices for no detection
			det_objslc[islc] = -1;
		}
	}
	// detection plane flagging logics:
	// - ndetper == 0 -> only last slice, no entrance plane
	//            > 0 -> at least entrance plane and exit plane
	//                   plus planes corresponding to multiples of ndetper
	// - nobjslc == 0 -> only this plane, regardless of ndetper
	//            > 0 -> always last plane, other planes possible depending on ndetper
	if (ndetper > 0) { // handle intermediate detections
		if (noslc > 0) { // object slices are present
			for (islc = 0; islc <= noslc; islc++) { // loop from entrance plane to exit plane
				if (0 == islc % ndetper) { // readout at intermetiate periodic planes
					if (det_objslc) {
						det_objslc[islc] = ndetslc;
					}
					ndetslc++;
				}
			}
			if (0 != noslc%ndetper) { // add out-of-period exit-plane
				if (det_objslc) {
					det_objslc[noslc] = ndetslc;
				}
				ndetslc++;
			}
		}
		else { // readout at entrance plane only
			if (det_objslc) {
				det_objslc[0] = 0;
			}
			ndetslc = 1;
		}
	}
	else { // readout at exit plane only
		if (det_objslc) {
			det_objslc[noslc] = 0;
		}
		ndetslc = 1;
	}
	return ndetslc;
}


int CJMultiSlice::SetDetectionPlanes(std::vector<int> slc_det, int nobjslc, int * det_objslc)
{
	int islc = 0;
	int ndetslc = 0;
	int ndetlist = (int)slc_det.size();
	int noslc = __min(__max(0, nobjslc), ndetlist);
	if (det_objslc && nobjslc > 0) {
		for (islc = 0; islc <= nobjslc; islc++) { // preset all slices for no detection
			det_objslc[islc] = -1;
		}
	}
	if (ndetlist > 0) { // handle intermediate detections
		if (noslc > 0) { // object slices are present
			for (islc = 0; islc <= noslc; islc++) { // loop from entrance plane to exit plane
				if (0 <= slc_det[islc]) { // readout at intermetiate periodic planes
					if (det_objslc) {
						det_objslc[islc] = slc_det[islc];
					}
					ndetslc++;
				}
			}
		}
		else { // readout at entrance plane only
			if (det_objslc) {
				det_objslc[0] = 0;
			}
			ndetslc = 1;
		}
	}
	else { // readout at exit plane only
		if (det_objslc) {
			det_objslc[nobjslc] = 0;
		}
		ndetslc = 1;
	}
	return ndetslc;
}


int CJMultiSlice::DetectorSetup(int whichcode, int ndetper, int ndetint, int imagedet, int nthreads_CPU)
{
	// Assumes prior object slice setup done and successful via ObjectSliceSetup
	int nerr = 0;
	int islc = 0;
	int64_t ndetitems = 0;
	size_t nitems = (size_t)m_nscx*m_nscy;
	int nthreads = __max(1, nthreads_CPU);
	size_t nbytes = 0;
	// clean previous setup (should usually not be necessary)
	m_ndet = 0;
	m_ndetper = 0;
	m_imagedet = (int)_JMS_DETECT_INTEGRATED;
	m_ndetslc = 0;
	m_threads_CPU_out = 0;

	// check for prerequisites
	if (((whichcode & (int)_JMS_CODE_CPU) > 0) && ((m_status_setup_CPU&_JMS_STATUS_OBJ) == 0)) {
		std::cerr << "Error: (DetectorSetup): no previous object slice setup for CPU code." << std::endl;
		nerr = 1;
	}
	if (((whichcode & (int)_JMS_CODE_GPU) > 0) && ((m_status_setup_GPU&_JMS_STATUS_OBJ) == 0)) {
		std::cerr << "Error: (DetectorSetup): no previous object slice setup for GPU code." << std::endl;
		nerr = 101;
	}
	if (nerr > 0) goto _Exit;

	// Preparations common to both codes
	m_imagedet = imagedet;
	// - determine the number of slices with detection and set detection hash table
	nbytes = sizeof(int)*(m_nobjslc + 1);
	if (0 < AllocMem_h((void**)&m_det_objslc, nbytes, "DetectorSetup", "slice -> detector hash")) { nerr = 2; goto _Exit; }
	for (islc = 0; islc <= m_nobjslc; islc++) { // preset for no detection
		m_det_objslc[islc] = -1;
	}
	m_ndetper = __max(0, ndetper);
	m_ndetslc = SetDetectionPlanes(m_ndetper, m_nobjslc, m_det_objslc);
	if (ndetint > 0) { // use integrating diffraction plane detectors
		m_ndet = ndetint; // set number of detectors
		if (0 < AllocMem_h((void**)&m_detmask_len, sizeof(int)*m_ndet, "DetectorSetup", "detector mask lengths", true)) { nerr = 7; goto _Exit; }
	}

	// CPU code - detector setup -> Memory prepared but no data set.
	if ((whichcode & (int)_JMS_CODE_CPU) && (m_imagedet > (int)_JMS_DETECT_NONE)) {
		m_threads_CPU_out = nthreads;
		if ((m_ndet > 0) && (m_imagedet & (int)_JMS_DETECT_INTEGRATED)) { // Prepare integrating detector arrays
			nbytes = sizeof(float*)*(size_t)m_ndet;
			if (0 < AllocMem_h((void**)&m_h_det, nbytes, "DetectorSetup", "detector function addresses", true))  { nerr = 3; goto _Exit; }
			nbytes = sizeof(int*)*(size_t)m_ndet;
			if (0 < AllocMem_h((void**)&m_h_detmask, nbytes, "DetectorSetup", "detector mask addresses", true)) { nerr = 8; goto _Exit; }
			nbytes = sizeof(float)*(size_t)m_ndet*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_int, nbytes, "DetectorSetup", "integrated detector channels", true)) { nerr = 4; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_IMAGE) { // prepare image plane readout arrays
			nbytes = sizeof(float)*nitems*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_img, nbytes, "DetectorSetup", "image planes", true)) { nerr = 5; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_DIFFRACTION) { // prepare diffraction plane readout arrays
			nbytes = sizeof(float)*nitems*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_dif, nbytes, "DetectorSetup", "diffraction planes", true)) { nerr = 6; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_WAVEREAL) { // prepare psi(r) readout arrays
			nbytes = sizeof(fcmplx)*nitems*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_wfr, nbytes, "DetectorSetup", "real-space wavefunction", true)) { nerr = 7; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_WAVEFOURIER) { // prepare psi(k) readout arrays
			nbytes = sizeof(fcmplx)*nitems*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_wff, nbytes, "DetectorSetup", "Fourier-space wavefunction", true)) { nerr = 9; goto _Exit; }
		}
	}
	
	// GPU code - detector setup -> Memory prepared but no data set.
	if ((whichcode & (int)_JMS_CODE_GPU) && (m_imagedet > (int)_JMS_DETECT_NONE)) {
		ndetitems = (int64_t)m_ndet*nitems;
		if ((ndetitems > 0) && (m_imagedet & (int)_JMS_DETECT_INTEGRATED)) { // make sure to allocate integrated detectors only if data is expected.
			nbytes = sizeof(float)*(size_t)ndetitems;
			if (0 < AllocMem_d((void**)&m_d_det, nbytes, "DetectorSetup", "detector functions", true)) { nerr = 103; goto _Exit; }
			nbytes = sizeof(int)*(size_t)ndetitems;
			if (0 < AllocMem_d((void**)&m_d_detmask, nbytes, "DetectorSetup", "detector masks", true)) { nerr = 108; goto _Exit; }
			nbytes = sizeof(float)*(size_t)m_ndet*m_ndetslc;
			if (0 < AllocMem_h((void**)&m_d_det_int, nbytes, "DetectorSetup", "integrated detector channels, GPU", true)) { nerr = 104; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_IMAGE) { // prepare image plane readout arrays
			nbytes = sizeof(float)*nitems*m_ndetslc; // calculate number of bytes to allocate
			if (0 < AllocMem_d((void**)&m_d_det_img, nbytes, "DetectorSetup", "image planes", true)) { nerr = 105; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_DIFFRACTION) { // prepare diffraction plane readout arrays
			nbytes = sizeof(float)*nitems*m_ndetslc; // calculate number of bytes to allocate
			if (0 < AllocMem_d((void**)&m_d_det_dif, nbytes, "DetectorSetup", "diffraction planes", true)) { nerr = 106; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_WAVEREAL) { // prepare psi(r) readout arrays
			nbytes = sizeof(cuComplex)*nitems*m_ndetslc; // calculate number of bytes to allocate
			if (0 < AllocMem_d((void**)&m_d_det_wfr, nbytes, "DetectorSetup", "real-space wavefunction", true)) { nerr = 107; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_WAVEFOURIER) { // prepare psi(k) readout arrays
			nbytes = sizeof(cuComplex)*nitems*m_ndetslc; // calculate number of bytes to allocate
			if (0 < AllocMem_d((void**)&m_d_det_wff, nbytes, "DetectorSetup", "Fourier-space wavefunction", true)) { nerr = 109; goto _Exit; }
		}
	}
	
_Exit:
	if (nerr == 0) {
		if (whichcode & (int)_JMS_CODE_CPU) m_status_setup_CPU |= (int)_JMS_STATUS_DET;
		if (whichcode & (int)_JMS_CODE_GPU) m_status_setup_GPU |= (int)_JMS_STATUS_DET;
	}
	return nerr;
}


int CJMultiSlice::DetectorSetup(int whichcode, std::vector<int> slc_det, int ndetint, int imagedet, int nthreads_CPU)
{
	// Assumes prior object slice setup done and successful via ObjectSliceSetup
	int nerr = 0;
	int islc = 0;
	int64_t ndetitems = 0;
	size_t nitems = (size_t)m_nscx*m_nscy;
	int nthreads = __max(1, nthreads_CPU);
	size_t nbytes = 0;
	// clean previous setup (should usually not be necessary)
	m_ndet = 0;
	m_ndetper = 0;
	m_imagedet = (int)_JMS_DETECT_INTEGRATED;
	m_ndetslc = 0;
	m_threads_CPU_out = 0;

	// check for prerequisites
	if (((whichcode & (int)_JMS_CODE_CPU) > 0) && ((m_status_setup_CPU&_JMS_STATUS_OBJ) == 0)) {
		std::cerr << "Error: (DetectorSetup): no previous object slice setup for CPU code." << std::endl;
		nerr = 1;
	}
	if (((whichcode & (int)_JMS_CODE_GPU) > 0) && ((m_status_setup_GPU&_JMS_STATUS_OBJ) == 0)) {
		std::cerr << "Error: (DetectorSetup): no previous object slice setup for GPU code." << std::endl;
		nerr = 101;
	}
	if (nerr > 0) goto _Exit;

	// Preparations common to both codes
	m_imagedet = imagedet;
	// - determine the number of slices with detection and set detection hash table
	nbytes = sizeof(int)*(m_nobjslc + 1);
	if (0 < AllocMem_h((void**)&m_det_objslc, nbytes, "DetectorSetup", "slice -> detector hash")) { nerr = 2; goto _Exit; }
	for (islc = 0; islc <= m_nobjslc; islc++) { // preset for no detection
		m_det_objslc[islc] = -1;
	}
	m_ndetslc = SetDetectionPlanes(slc_det, m_nobjslc, m_det_objslc);
	if (ndetint > 0) { // use integrating diffraction plane detectors
		m_ndet = ndetint; // set number of detectors
		if (0 < AllocMem_h((void**)&m_detmask_len, sizeof(int)*m_ndet, "DetectorSetup", "detector mask lengths", true)) { nerr = 7; goto _Exit; }
	}

	// CPU code - detector setup -> Memory prepared but no data set.
	if ((whichcode & (int)_JMS_CODE_CPU) && (m_imagedet > (int)_JMS_DETECT_NONE)) {
		m_threads_CPU_out = nthreads;
		if ((m_ndet > 0) && (m_imagedet & (int)_JMS_DETECT_INTEGRATED)) { // Prepare integrating detector arrays
			nbytes = sizeof(float*)*(size_t)m_ndet;
			if (0 < AllocMem_h((void**)&m_h_det, nbytes, "DetectorSetup", "detector function addresses", true)) { nerr = 3; goto _Exit; }
			nbytes = sizeof(int*)*(size_t)m_ndet;
			if (0 < AllocMem_h((void**)&m_h_detmask, nbytes, "DetectorSetup", "detector mask addresses", true)) { nerr = 8; goto _Exit; }
			nbytes = sizeof(float)*(size_t)m_ndet*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_int, nbytes, "DetectorSetup", "integrated detector channels", true)) { nerr = 4; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_IMAGE) { // prepare image plane readout arrays
			nbytes = sizeof(float)*nitems*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_img, nbytes, "DetectorSetup", "image planes", true)) { nerr = 5; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_DIFFRACTION) { // prepare diffraction plane readout arrays
			nbytes = sizeof(float)*nitems*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_dif, nbytes, "DetectorSetup", "diffraction planes", true)) { nerr = 6; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_WAVEREAL) { // prepare psi(r) readout arrays
			nbytes = sizeof(fcmplx)*nitems*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_wfr, nbytes, "DetectorSetup", "real-space wavefunction", true)) { nerr = 7; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_WAVEFOURIER) { // prepare psi(k) readout arrays
			nbytes = sizeof(fcmplx)*nitems*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_wff, nbytes, "DetectorSetup", "Fourier-space wavefunction", true)) { nerr = 9; goto _Exit; }
		}
	}

	// GPU code - detector setup -> Memory prepared but no data set.
	if ((whichcode & (int)_JMS_CODE_GPU) && (m_imagedet > (int)_JMS_DETECT_NONE)) {
		ndetitems = (int64_t)m_ndet*nitems;
		if ((ndetitems > 0) && (m_imagedet & (int)_JMS_DETECT_INTEGRATED)) { // make sure to allocate integrated detectors only if data is expected.
			nbytes = sizeof(float)*(size_t)ndetitems;
			if (0 < AllocMem_d((void**)&m_d_det, nbytes, "DetectorSetup", "detector functions", true)) { nerr = 103; goto _Exit; }
			nbytes = sizeof(int)*(size_t)ndetitems;
			if (0 < AllocMem_d((void**)&m_d_detmask, nbytes, "DetectorSetup", "detector masks", true)) { nerr = 108; goto _Exit; }
			nbytes = sizeof(float)*(size_t)m_ndet*m_ndetslc;
			if (0 < AllocMem_h((void**)&m_d_det_int, nbytes, "DetectorSetup", "integrated detector channels, GPU", true)) { nerr = 104; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_IMAGE) { // prepare image plane readout arrays
			nbytes = sizeof(float)*nitems*m_ndetslc; // calculate number of bytes to allocate
			if (0 < AllocMem_d((void**)&m_d_det_img, nbytes, "DetectorSetup", "image planes", true)) { nerr = 105; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_DIFFRACTION) { // prepare diffraction plane readout arrays
			nbytes = sizeof(float)*nitems*m_ndetslc; // calculate number of bytes to allocate
			if (0 < AllocMem_d((void**)&m_d_det_dif, nbytes, "DetectorSetup", "diffraction planes", true)) { nerr = 106; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_WAVEREAL) { // prepare psi(r) readout arrays
			nbytes = sizeof(cuComplex)*nitems*m_ndetslc; // calculate number of bytes to allocate
			if (0 < AllocMem_d((void**)&m_d_det_wfr, nbytes, "DetectorSetup", "real-space wavefunction", true)) { nerr = 107; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_WAVEFOURIER) { // prepare psi(k) readout arrays
			nbytes = sizeof(cuComplex)*nitems*m_ndetslc; // calculate number of bytes to allocate
			if (0 < AllocMem_d((void**)&m_d_det_wff, nbytes, "DetectorSetup", "Fourier-space wavefunction", true)) { nerr = 109; goto _Exit; }
		}
	}

_Exit:
	if (nerr == 0) {
		if (whichcode & (int)_JMS_CODE_CPU) m_status_setup_CPU |= (int)_JMS_STATUS_DET;
		if (whichcode & (int)_JMS_CODE_GPU) m_status_setup_GPU |= (int)_JMS_STATUS_DET;
	}
	return nerr;
}



void CJMultiSlice::SetSliceThickness(int islc, float fthickness)
{
	if (islc >= 0 && islc < m_nscslc && NULL != m_slcthick) {
		m_slcthick[islc] = fthickness;
	}
}

float CJMultiSlice::GetSliceThickness(int islc)
{
	if (islc >= 0 && islc < m_nscslc && NULL != m_slcthick) {
		return m_slcthick[islc];
	}
	return 0.0f;
}

void CJMultiSlice::SetPlasmonMC(bool use, float q_e, float q_c, float mfp, unsigned int nexmax)
{
	m_plasmonmc_flg = 0;
	if (use) m_plasmonmc_flg = 1;
	m_jplmc.m_q_e = abs(q_e);
	m_jplmc.m_q_c = abs(q_c);
	m_jplmc.m_meanfreepath = abs(mfp);
	m_jplmc.m_sca_num_max = nexmax;
}

int CJMultiSlice::CalculateProbeWaveFourier(CJProbeParams* prm, fcmplx *wav)
{
	// Assumes that general parameters are set up.
	if (NULL == prm || NULL == wav) { return 1; }
	if (m_wl <= 0.f || m_a0[0] * m_a0[1] <= 0.f || m_nscx*m_nscy <= 0 || prm->m_alpha <= 0.f) {
		return 1; // error, invalid setup
	}
	// uses m_a0, m_nscx, m_nscy, m_wl from CJMultiSlice
	CJProbeParams lprm;
	lprm = *prm;
	lprm.m_wl = m_wl;
	int nerr = m_jpg.CalculateProbeWaveFourier(&lprm, m_nscx, m_nscy, m_a0[0], m_a0[1], wav);
	return nerr;
}

int CJMultiSlice::CalculatePropagator(float fthick, float otx, float oty, fcmplx *pro, int ntype)
{
	if (pro == NULL || ntype < 0 || ntype > 1) {
		return 1;
	}
	int nitems = m_nscx*m_nscy; // number of pixels
	int nyqx = (m_nscx - (m_nscx % 2)) >> 1; // x nyquist number
	int nyqy = (m_nscy - (m_nscy % 2)) >> 1; // y nyquist number
	double itogx = 1. / m_a0[0]; // Fourier space sampling rate along x
	double itogy = 1. / m_a0[1]; // Fourier space sampling rate along y
	double gmax = __min(itogx*nyqx, itogy*nyqy); // smaller of the x,y gmax
	double gthr = gmax * _JMS_RELAPERTURE; // band-width limitation of g
	double gthr2 = gthr*gthr; // band-width limitation of g^2
	double otxr = (double)otx * _D2R; // object tilt x in rad
	double otyr = (double)oty * _D2R; // object tilt y in rad
	double wl = (double)GetWaveLength(); // electron wavelength [nm]
	double ot = sqrt(otxr*otxr + otyr*otyr); // object tilt magnitude in rad
	double chi = 0.; // phase-shift to diffracted beams in rad
	float pronorm = 1.f/((float)m_nscx*m_nscy); // amplitude normalization for JMS
	int i = 0, j = 0, i1 = 0, j1 = 0, idy = 0; // iterators
	fcmplx c0(0.f, 0.f), c1(0.f, 0.f);

	// fork to different propagator versions

	if (ntype == 0) {
		// planar propagator version
		m_nprotype = 0;
		double pfac = 2. * _PI / wl; // chi pre-factor
		double cot = cos(ot);
		double sot = sin(ot);
		double od = 0.; // object tilt direction in rad
		if (ot > 0.0) {
			od = atan2(otyr, otxr);
		}
		// loop through frequencies (scrambled, not transposed here!)
		double gx = 0., gy = 0., gy2 = 0., g2 = 0.; // spatial frequencies [1/nm]
		double dt = 0., dd = 0.; // theta modulus and orientation in rad
		for (j = 0; j < m_nscy; j++) { // loop rows, y frequencies
			j1 = ((j + nyqy) % m_nscy) - nyqy; // scrambled y frequency index
			gy = itogy * j1;
			gy2 = gy*gy;
			idy = j*m_nscx;
			for (i = 0; i < m_nscx; i++) { // loop columns, x frequencies
				i1 = ((i + nyqx) % m_nscx) - nyqx; // scrambled x frequency index
				gx = itogx * i1;
				g2 = gx*gx + gy2;
				if (g2 > gthr2) { // out of band-width limit
					pro[i + idy] = c0; // set to zero
				}
				else { // in band width limit
					dt = asin(wl*sqrt(g2)); // theta in rad
					dd = 0.;
					if (dt > 0.) {
						dd = atan2(gy, gx);
					}
					chi = pfac*(double)fthick*(1. / (cot*cos(dt) + sot*sin(dt)*cos(od - dd)) - 1. / cot);
					c1.real((float)cos(chi)*pronorm);
					c1.imag((float)(-sin(chi)*pronorm));
					pro[i + idy] = c1;
				}
			}
		}
	}

	if (ntype == 1) {
		// fresnel propagator version
		m_nprotype = 1;
		double pfac = _PI / wl * cos(ot); // chi pre-factor with thickness reduction 
		double pt0 = pfac * (otxr*otxr + otyr*otyr); // chi0 - reference phase ! Pi/lambda * Cos[tilt] * t^2
		double itowx = itogx * wl; // theta x per pixel in rad
		double itowy = itogy * wl; // theta y per pixel in rad
		// loop through frequencies (scrambled, not transposed here!)
		double wx = 0., wy = 0., wy2 = 0., w2 = 0., ux = 0., uy = 0., uy2 = 0., u2 = 0.; // spatial frequencies [rad]
		double wthr2 = gthr2*wl*wl; // band-width limitation threshold on theta^2
		double dt = 0., dd = 0.; // theta modulus and orientation in rad
		for (j = 0; j < m_nscy; j++) { // loop rows, y frequencies
			j1 = ((j + nyqy) % m_nscy) - nyqy; // scrambled y frequency index
			wy = itowy * j1;
			wy2 = wy*wy;
			uy = wy - otyr; // tilt y
			uy2 = uy*uy;
			idy = j*m_nscx;
			for (i = 0; i < m_nscx; i++) { // loop columns, x frequencies
				i1 = ((i + nyqx) % m_nscx) - nyqx; // scrambled x frequency index
				wx = itowx * i1;
				w2 = wx*wx + wy2;
				ux = wx - otxr; // tilt x
				u2 = ux*ux + uy2;
				if (w2 > wthr2) {
					// out of band-width limt
					pro[i + idy] = c0; // set to zero
				}
				else {
					// in band width limit
					// chi(k) = dZ * Pi / lambda * Cos[tilt] * (w ^ 2 - 2 w.t)
					//                                          ^Prop     ^Tilt
					// The tilted propagator effectively shifts the slices
					// against each other along the tilt direction.
					chi = (double)fthick*(pfac*u2 - pt0);
					c1.real((float)cos(chi)*pronorm);
					c1.imag((float)(-sin(chi)*pronorm));
					pro[i + idy] = c1;
				}
			}
		}
	}

	return 0;
}

void CJMultiSlice::DiffractionDescan(bool bActivate)
{
	if (bActivate) {
		m_dif_descan_flg = 1;
	}
	else {
		m_dif_descan_flg = 0;
	}
}

void CJMultiSlice::SetDiffractionDescanN(int whichcode, int ndescanx, int ndescany, int iThread)
{
	if ((0 < (whichcode & (int)_JMS_CODE_CPU)) && (0 < (whichcode & (int)_JMS_CODE_GPU))) {
		goto _Exit; // invalid whichcode
	}
	if ((0 == (whichcode & (int)_JMS_CODE_CPU)) && (0 == (whichcode & (int)_JMS_CODE_GPU))) {
		goto _Exit; // invalid whichcode
	}
	if (whichcode & (int)_JMS_CODE_CPU && (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC) > 0 &&
			iThread >= 0 && iThread < m_ncputhreads) {
		m_h_dif_ndescanx[iThread] = ndescanx; // 1/nm -> pixel
		m_h_dif_ndescany[iThread] = ndescany; // 1/nm -> pixel
	}
	if (whichcode & (int)_JMS_CODE_GPU && (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC) > 0) {
		m_d_dif_ndescanx = ndescanx; // 1/nm -> pixel;
		m_d_dif_ndescany = ndescany; // 1/nm -> pixel;
	}
_Exit:
	return;
}

void CJMultiSlice::SetDiffractionDescan(int whichcode, float descanx, float descany, int iThread)
{
	if ((0 < (whichcode & (int)_JMS_CODE_CPU)) && (0 < (whichcode & (int)_JMS_CODE_GPU))) {
		goto _Exit; // invalid whichcode
	}
	if ((0 == (whichcode & (int)_JMS_CODE_CPU)) && (0 == (whichcode & (int)_JMS_CODE_GPU))) {
		goto _Exit; // invalid whichcode
	}
	if (whichcode & (int)_JMS_CODE_CPU && (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC) > 0 &&
		iThread >= 0 && iThread < m_ncputhreads) {
		m_h_dif_ndescanx[iThread] = (int)round(descanx * m_a0[0]); // 1/nm -> pixel
		m_h_dif_ndescany[iThread] = (int)round(descany * m_a0[1]); // 1/nm -> pixel
	}
	if (whichcode & (int)_JMS_CODE_GPU && (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC) > 0) {
		m_d_dif_ndescanx = (int)round(descanx * m_a0[0]); // 1/nm -> pixel;
		m_d_dif_ndescany = (int)round(descany * m_a0[1]); // 1/nm -> pixel;
	}
_Exit:
	return;
}

void CJMultiSlice::SetDiffractionDescanMRad(int whichcode, float descanx, float descany, int iThread)
{
	// mrad -> 1/nm
	return SetDiffractionDescan(whichcode, descanx*0.001f/m_wl, descany*0.001f / m_wl, iThread);
}

int CJMultiSlice::LoadSTEMDetectorProfile(std::string sfile, int &len, float &refpix, float** profile)
{
	// init
	int nResult = 0;
	float ftmp = 0.f;
	std::string sline;
	std::string ssub;
	std::string ssep = "' ";
	len = 0;
	refpix = 0.f;
	*profile = NULL;
	int idx = 0, nline = 0;
	int i0 = 0, i1 = 0;
	// try opening the file for reading list directed data
	std::ifstream instr( sfile );
	if (instr.is_open()) {
		while (getline(instr, sline)) {
			nline++;
			if (nline == 1) { // header line
				// parse two numbers separated by comma and white spaces
				i0 = (int)sline.find_first_not_of(ssep, 0); // get first character index which is no separator
				i1 = (int)sline.find_first_of(ssep, i0); // get the next character index which is a separator
				// the part of the string from i0 to i1 is the first (int) number
				if (i0 != std::string::npos && i1 != std::string::npos && i0 < i1) {
					ssub = sline.substr(i0, i1 - i0);
					len = stoi(ssub);
				}
				else {
					nResult = 2;
					goto _exit;
				}
				if (len > 0) {
					i0 = (int)sline.find_first_not_of(ssep, i1); // get first index of the next item
					i1 = (int)sline.find_first_of(ssep, i0); // get the next separator index
					if (i1 == std::string::npos) { // in case no more separators are found ...
						i1 = (int)sline.length(); // set the length of the string as terminating index
					}
					// the part of the string from i0 to i1 is the second (float) number
					if (i0 != std::string::npos && i0 < i1) {
						ssub = sline.substr(i0, i1 - i0);
						refpix = stof(ssub);
					}
					else {
						nResult = 3;
						goto _exit;
					}
				}
				// length>0 and reference pixel index have been read in
				// allocate the profile buffer (slightly larger)
				*profile = (float*)malloc(sizeof(float)*(len+10));
				if (NULL == *profile) {
					nResult = 4;
					goto _exit;
				}
				memset(*profile, 0, sizeof(float)*(len+10)); // preset buffer with zeroes
				//
				// ... ready for reading more data
			}
			else {
				// to float (expecting no other stuff in the line)
				(*profile)[idx] = stof(sline);
				idx++;
			}
		}
		instr.close();
	}
	else {
		nResult = 1;
		goto _exit;
	}
_exit:
	if (0 < nResult) { // handle allocations on error
		if (NULL != *profile) {
			free(*profile); *profile = NULL;
		}
		len = 0; refpix = 0.f;
		if (instr.is_open()) { // handle open input stream
			instr.close();
		}
	}
	return nResult;
}

float CJMultiSlice::GetRadialDetectorSensitivity(float theta, float* profile, int len, float refpix, float beta1)
{
	// This function returns the detector sensitivity for a given diffraction angle.
	// It assumes that the first pixel of the profile corresponds to theta = 0.
	// The angular calibration is done via refpix and beta1 which identify the pixel+1 and diffraction angle
	// providing scale for a linear dependency. Linear interpolation is used.
	float sens = 1.0f;
	float tpix = 0.f; // target pixel of the profile
	float tfrac = 0.f;
	float rpix = __max(0.f, refpix - 1.f); // reduce refpix by 1 (this value is the sample number, starting to count with 1)
	int ipix = 0;
	if (NULL != profile && len > 0 && rpix >= 0.f) {
		tpix = __max(0.f, theta * rpix / beta1); // get target pixel
		ipix = (int)tpix; // get integer part for linear interpolation
		tfrac = tpix - ipix; // get fractional part for linear interpolation
		if (ipix > (len - 2)) { // special case at the upper bound
			sens = profile[len-1]; // ... set the sensitivity to that of the last profile value
		}
		else { // default case inside the range
			sens = (1.f-tfrac)*profile[ipix] + tfrac*profile[ipix + 1];
		}
	}
	return sens;
}

int CJMultiSlice::CalculateRingDetector(float beta0, float beta1, float phi0, float phi1, float theta0x, float theta0y, std::string sdsprofile, float *det, int &msklen, int *msk)
{
	int nuseprofile = 0; // preset the use of a sensitivity profile to 'no' by default
	int nitems = m_nscx * m_nscy; // set number of detector function items
	bool bmakemask = false; // preset the creation of a mask to 'no' by default
	bool busesprof = false; // preset the sensitivity profile usage to 'no' by default
	msklen = 0; // reset mask length to 0
	if (NULL != msk) {
		bmakemask = true;
	}
	memset(det, 0, sizeof(float)*nitems); // reset detector function
	if (bmakemask) {
		memset(msk, 0, sizeof(int)*nitems); // reset detector mask
	}
	int nDetProfileLen = 0; // number of items of the detector profile
	float *pfDetProfile = NULL; // detector profile data buffer
	float fRefPixelBeta1 = 0.0f; // pixel of the detector profile corresponding to beta1 of the detector parameter set
	int nyqx = (m_nscx - (m_nscx % 2)) >> 1; // x nyquist number
	int nyqy = (m_nscy - (m_nscy % 2)) >> 1; // y nyquist number
	double itogx = 1. / m_a0[0]; // Fourier space sampling rate along x
	double itogy = 1. / m_a0[1]; // Fourier space sampling rate along y
	double wl = (double)GetWaveLength(); // electron wavelength [nm]
	double k0 = 1. / wl; // k0 [1/nm]
	float sens = 1.f; // sensitivity o fthe detector (1. by default)
	double dr1 = 0., dr2 = 0., da1 = 0., da2 = 0., dc1 = 0., dc2 = 0.; // internal detector parameters [rad]
	double dgc1 = 0., dgc2 = 0.; // detector center [1/nm]
	double gx = 0., gy = 0., gy2 = 0., g = 0., phi = 0.; // reciprocal space (projected)
	int i = 0, j = 0, i1 = 0, j1 = 0, idx = 0, idy = 0; // iterators
	int nring = 0; // full ring flag
	int nerr = 0; // error code
	//
	// Detector sensitivity modification (optional)
	if (0 < sdsprofile.length()) { // there is a string in the profile name parameter
		// try loading the detector profile data
		nerr = LoadSTEMDetectorProfile(sdsprofile, nDetProfileLen, fRefPixelBeta1, &pfDetProfile);
		if ( 0 < nerr) { // error while reading the detector profile
			std::cout << "(CJMultiSlice::CalculateRingDetector) Warning: Failed to read detector sensitivity profile." << std::endl;
			std::cout << "  file name: " << sdsprofile << std::endl;
			std::cout << "  proceeding with default sensitivity of 1." << std::endl;
			nDetProfileLen = 0;
		}
		else {
			busesprof = true; // activate to use the sensitivity profile
		}
	}
	//
	// Re-arrange the detector paramters to the 1/nm units of the calculation grid.
	// This assumes that the angle paramters are determined for an axial, planar detector,
	// which may be shifted in the projected plane.
	// - translate parameters to radian units
	dr1 = (double)beta0 * 0.001;
	dr2 = (double)beta1 * 0.001;
	da1 = (double)phi0 * _D2R;
	da2 = (double)phi1 * _D2R;
	dc1 = (double)theta0x * 0.001;
	dc2 = (double)theta0y * 0.001;
	// - translate the detector shift to 1/nm units
	if (m_nprotype == 1) { // fresnel propagator (small angle approximations)
		dr1 *= k0;
		dr2 *= k0;
		dgc1 *= k0;
		dgc2 *= k0;
	}
	else { // default large angle propagator
		dr1 = k0 * sin(dr1);
		dr2 = k0 * sin(dr2);
		phi = 0.;
		gy2 = k0 * sin( sqrt(dc1*dc1 + dc2*dc2) ); // theta0 to 1/nm
		if (gy2 > 0.) {
			phi = atan2(dc2, dc1); // azimuth in rad
		}
		dgc1 = gy2 * cos(phi); // gx
		dgc2 = gy2 * sin(phi); // gy
	}
	// - check the angular ranges -> determine whether to use a full ring or a segment
	if (abs(da1 - da2) < 0.0001) { // extremely small azimuth parameter difference
		nring++;
	}
	if (abs(abs(da1 - da2) - _TPI) < 0.0001) { // azimuth difference is two pi
		nring++;
	}
	if (nring > 0) {
		nring = 1;
	}
	if (nring == 0) { // ring segment setup
		// setup aziumth range from 0 up to 2 * pi in case of a segment detector with azimuthal range
		// add 2pi to the azimuths until both are definitely positive
		while (da1 < 0.0 || da2 < 0.0) {
			da1 += _TPI;
			da2 += _TPI;
		}
		// da2 should always be larger than da1, if not, than add 2pi to da2
		if (da1 >= da2) {
			da2 += _TPI;
		}
		// now bring both angles back to the first rotation cycle, at least for da1
		while (da1 >= _TPI) {
			da1 -= _TPI;
			da2 -= _TPI;
		}
	}
	//
	// Go through the Fourier space area used by the multislice and
	// set the detector sensistivity to each active pixel
	// Remember that in JMS qy is on the outer loop other than in the Fortran code!
	for (j = 0; j < m_nscy; j++) { // loop rows, y frequencies
		j1 = ((j + nyqy) % m_nscy) - nyqy; // scrambled y frequency index
		gy = itogy * j1 - dgc2; // gy of this row in 1/nm relatve to the detector y-center
		gy2 = gy * gy;
		idy = j * m_nscx;
		for (i = 0; i < m_nscx; i++) { // loop columns, x frequencies
			i1 = ((i + nyqx) % m_nscx) - nyqx; // scrambled x frequency index
			gx = itogx * i1 - dgc1; // gx of this pixel column in 1/nm relative to the detector x-center
			idx = i + idy; // pixel index
			g = sqrt( gy2 + gx * gx );
			phi = 0.;
			if (g > 0.) {
				phi = atan2(gy, gx);
			}
			// range checks on the following general interval : min <= x < max
			if (g < dr1 || g >= dr2) { // no detector pixel, try next pixel
				continue;
			}
			if (nring == 0) { // segment checks
				while (phi < da1) { // try to find azimuthal cycle beyond da1
					phi += _TPI;
				}
				if (phi < da1 || phi >= da2) { // no segment pixel, try next pixel
					continue;
				}
			}
			// all position checks were passed
			// THIS IS A DETECTOR PIXEL
			if (busesprof) { // get detector sensitivity from profile
				sens = 1.0f;
				if (dr1 > 0.) { // require the reference angle to be non-zero and positive
					sens = GetRadialDetectorSensitivity((float)g, pfDetProfile, nDetProfileLen, fRefPixelBeta1, (float)dr1);
				}
				else if (dr2 > 0.) { // handle possible case of a disk detector, where the outer diameter is used for reference
					sens = GetRadialDetectorSensitivity((float)g, pfDetProfile, nDetProfileLen, fRefPixelBeta1, (float)dr2);
				}
				// If none of the above two options is used, the default sensitivity of 1 is used.
			}
			det[idx] = sens; // set pixel sensitivity
			if (bmakemask) { // add pixel index to the end of the mask list
				msk[msklen] = idx;
				msklen++;
			}
		} // loop columns
	} // loop rows
	if (NULL != pfDetProfile) {
		free(pfDetProfile);
		pfDetProfile = NULL;
	}
	return 0; // exit with success code
}

int CJMultiSlice::SetPhaseGratingData(int whichcode, int islc, int nvar, fcmplx* pgr)
{
	// Assumes previous call of PhaseGratingSetup
	int nerr = 0;
	int nitems = m_npgx*m_npgy;
	int nlvar = 0; // local and used number of variants
	cudaError cuerr;
	size_t nbytes = 0;
	cuComplex* pgrdst = NULL;
	if (islc >= 0 && islc < m_nscslc && NULL != pgr && NULL != m_nvarslc) { // valid structure slice ID and pointers
		nlvar = __min(m_nvarslc[islc],nvar); // update number of variants
		if (nlvar < m_nvarslc[islc]) {
			m_nvarslc[islc] = nlvar;
		}
		nbytes = sizeof(fcmplx)*(size_t)nlvar*nitems;
		if (nbytes == 0) { nerr = 1; goto _Exit; } // something is wrong with the setup
		if (whichcode & (int)_JMS_CODE_CPU || m_d_pgr_src==1) { // CPU code or phase-grating host sourceing is used
			if (NULL == m_h_pgr) { nerr = 2; goto _Exit; } // cannot link, pgr array not ready
			m_h_pgr[islc] = pgr; // just link the pointer to host address, the rest is handled by lib parameters
		}
		if (whichcode & (int)_JMS_CODE_GPU && m_d_pgr_src==0) { // copy from host to device only on default pre-copy mode
			if (NULL == m_d_pgr || NULL == m_slcoffset) { nerr = 101; goto _Exit; } // cannot copy, pgr array not ready
			pgrdst = m_d_pgr + m_slcoffset[islc];
			cuerr = cudaMemcpy(pgrdst, pgr, nbytes, cudaMemcpyHostToDevice); // copy to device
			if (cuerr != cudaSuccess) { // handle device error
				PostCUDAError("(SetPhaseGratingData): Failed to copy phase grating data to GPU devices", cuerr);
				nerr = 102; goto _Exit;
			}
		}
	}
_Exit:
	return nerr;
}


int CJMultiSlice::SetPropagatorData(int whichcode, int ipro, fcmplx* pro)
{
	// Assumes previous call of PropagatorSetup
	int nerr = 0;
	int nitems = m_nscx*m_nscy;
	cudaError cuerr;
	size_t nbytes = 0;
	cuComplex* prodst = NULL;
	if (ipro >= 0 && ipro < m_npro && NULL != pro ) { // valid structure slice ID and pointers
		nbytes = sizeof(fcmplx)*(size_t)nitems; // bytes per propagator
		if (nbytes == 0) { nerr = 1; goto _Exit; } // something is wrong with the setup
		if (whichcode & (int)_JMS_CODE_CPU) {
			if (NULL == m_h_pro) { nerr = 2; goto _Exit; } // cannot link, internal pro array not ready
			m_h_pro[ipro] = pro; // just link the pointer to host address, the rest is handled by lib parameters
		}
		if (whichcode & (int)_JMS_CODE_GPU) {
			if (NULL == m_d_pro) { nerr = 101; goto _Exit; } // cannot copy, internal pro array not ready
			prodst = m_d_pro + (size_t)ipro*nitems;
			cuerr = cudaMemcpy(prodst, pro, nbytes, cudaMemcpyHostToDevice); // copy to device
			if (cuerr != cudaSuccess) { // handle device error
				PostCUDAError("(SetPropagatorData): Failed to copy propagator data to GPU devices", cuerr);
				nerr = 102;	goto _Exit;
			}
		}
	}
_Exit:
	return nerr;
}


int CJMultiSlice::SetDetectorData(int whichcode, int idet, float* det, int msklen, int* msk)
{
	// Assumes previous call of DetectorSetup
	int nerr = 0;
	int nitems = m_nscx*m_nscy;
	cudaError cuerr;
	size_t nbytes = 0; // detector function # bytes
	size_t nbytes2 = 0; // detector mask max. # bytes
	size_t nbytes3 = 0; // detector mask used # bytes
	float* detdst = NULL;
	int *mskdst = NULL;
	if (idet >= 0 && idet < m_ndet && m_ndet > 0 && NULL != det) { // valid structure slice ID and pointers
		nbytes = sizeof(float)*(size_t)nitems; // bytes per detector function
		nbytes2 = sizeof(int)*(size_t)nitems; // max. bytes per mask
		if (nbytes == 0) { nerr = 1; goto _Exit; } // something is wrong with the setup
		if (NULL != m_detmask_len) { // detector mask length can be set
			m_detmask_len[idet] = 0; // default preset
			if ((NULL != msk) && (msklen > 0)) { // use input detector mask length
				m_detmask_len[idet] = msklen; // set
				nbytes3 = sizeof(int)*(size_t)msklen; // bytes per mask used
			}
		}
		if (whichcode & (int)_JMS_CODE_CPU) {
			if (NULL == m_h_det) { nerr = 2; goto _Exit; } // cannot link, internal det array not ready
			m_h_det[idet] = det; // just link the pointer to host address, the rest is handled by lib parameters
			if (NULL != m_h_detmask) { // detector mask can be linked
				m_h_detmask[idet] = NULL; // preset to none
				if (NULL != msk && 0 < msklen) { // use input mask
					m_h_detmask[idet] = msk; // set
				}
			}
		}
		if (whichcode & (int)_JMS_CODE_GPU) {
			if (NULL == m_d_det) { nerr = 101; goto _Exit; } // cannot copy, internal pro array not ready
			detdst = m_d_det + (size_t)idet*nitems;
			cuerr = cudaMemcpy(detdst, det, nbytes, cudaMemcpyHostToDevice); // copy to device
			if (cuerr != cudaSuccess) { // handle device error
				PostCUDAError("(SetDetectorData): Failed to copy detector function data to GPU devices", cuerr);
				nerr = 102;
				goto _Exit;
			}
			if (NULL != m_d_detmask) { // detector mask can be linked
				mskdst = m_d_detmask + (size_t)idet*nitems;
				cuerr = cudaMemset(mskdst, 0, nbytes2); // preset to zero mask
				if (cuerr != cudaSuccess) { // handle device error
					PostCUDAError("(SetDetectorData): Failed to preset detector mask on GPU devices", cuerr);
					nerr = 103;
					goto _Exit;
				}
				if (NULL != msk && 0 < msklen) { // use input mask
					cuerr = cudaMemcpy(mskdst, msk, nbytes3, cudaMemcpyHostToDevice); // copy to device
					if (cuerr != cudaSuccess) { // handle device error
						PostCUDAError("(SetDetectorData): Failed to copy detector mask data to GPU devices", cuerr);
						nerr = 104;
						goto _Exit;
					}
				}
			}
		}
	}
_Exit:
	return nerr;
}



int CJMultiSlice::GetRandomVariantSequence(int *out)
{
	// requires a previous setup of some basic module data
	// m_nscslc = number of slices
	// m_nvarslc = number of variants for each slice
	// m_objslc = object slice sequence
	
	// local variables
	int nerr = 0;
	int *q = NULL;
	int *qcur = NULL;
	int *qlen = NULL;
	int *qlast = NULL;
	int qlenmax = 0;
	int i = 0, j = 0, k = 0, k0 = 0;
	int islc = 0;
	int ivar = 0;

	// initial check for presence and setup of required data
	if (NULL == out) { nerr = 1; goto _Exit; }
	if (NULL == m_nvarslc || NULL == m_objslc || m_nscslc < 1 || m_nobjslc < 1) { nerr = 2; goto _Exit; }

	// allocate local memory used for sequence generation
	qlen = (int*)malloc(sizeof(int)*m_nscslc); // number of available variants for each slice ID
	qlast = (int*)malloc(sizeof(int)*m_nscslc); // last accessed index for each slice ID

	// preset length and access arrays
	for (i = 0; i<m_nscslc; i++) {
		qlen[i] = 0;
		qlast[i] = -1;
		qlenmax = __max(qlenmax, m_nvarslc[i]);
	}

	// allocate the index access queues
	q = (int*)malloc(sizeof(int)*m_nscslc*qlenmax);
	
	// Find a random variant index for each object slice by drawing without returning.
	// Also take care to avoid that the next drawing round starts with the same index
	// drawn at the end of the previous round.

	// Loop over all object slices
	for (i = 0; i<m_nobjslc; i++) {
		// get slice index
		islc = m_objslc[i];
		// quick-handle case with 1 variant
		if (m_nvarslc[islc] < 2) {
			out[i] = 0;
			continue;
		}
		// get pointer to the index access queue for this index
		qcur = q + islc*qlenmax;
		// check for an empty access list
		if (qlen[islc] <= 0) { // list is empty
			// prepare new list of numbers as ramp sequence
			for (j = 0; j < m_nvarslc[islc]; j++) {
				qcur[j] = j;
			}
			// set new length of the access list
			qlen[islc] = m_nvarslc[islc];
		}
		
		// draw a random number from the current access list
		//k = (int)((float)m_nvarslc[islc] * GetRand() / (float)(RAND_MAX + 1));
		k = (int)m_prng->unirand(0., (double)m_nvarslc[islc]);
		k = k%m_nvarslc[islc]; // periodic wrap around max queue length
		k0 = k; // remember this draw

		// get a valid draw from the access list
		while (qcur[k] < 0 || qcur[k] == qlast[islc]) { // the current draw is invalid. It has been drawn already or is equal to the last draw.
			k++; // increase to next item
			k = k%m_nvarslc[islc]; // periodic wrap around max queue length
			if (k == k0) { // oops, this should never happen
				break; // stop the while loop as it will get infinite
			}
		}
		// get preliminary variant index
		ivar = qcur[k];

		// handle unexpected problems
		if (ivar < 0) { // problem case -> redraw with returning
			//ivar = (int)((float)m_nvarslc[islc] * GetRand() / (float)(RAND_MAX + 1));
			ivar = (int)m_prng->unirand(0.0, (double)m_nvarslc[islc]);
			// reset access queue
			qlen[islc] = 0;
		}
		else { // default case -> draw without returning
			qcur[k] = -1;
			qlen[islc]--;
		}
		
		// Remember last draw.
		// This assumes, that number of variants is > 1, which is valid as
		// this line will only be executed if a slice passes the quick handling
		// in the beginning of the routine.
		qlast[islc] = ivar;

		// store variant index
		out[i] = ivar;

	}
_Exit:
	if (q != NULL) { free(q); q = NULL; }
	if (qlen != NULL) { free(qlen); qlen = NULL; }
	if (qlast != NULL) { free(qlast); qlast = NULL; }
	return nerr;
}


int CJMultiSlice::InitCore(int whichcode, int nCPUthreads)
{
	// The code assumes a previous setup of all calculation box parameters
	// via the routines SetGridSize, SetHighTension, SetSupercellSize
	// also all components of the setup chain need to be successfully called
	// before this function works: (m_status_setup_CPU & _JMS_THRESHOLD_CORE) > 0
	int nerr = 0;
	int i = 0, j = 0, i1 = 0, nx2 = 0, ny2 = 0;
	int icore = 0;
	int pdims[2];
	int ncore_ready = 0;
	pdims[0] = m_nscx;
	pdims[1] = m_nscy;
	int nitems = m_nscx*m_nscy;
	size_t nbytes = 0;
	float kscax = 1.f / m_a0[0], kscay = 1.f / m_a0[1], kx = 0.f, ky = 0.f;
	cudaError cuerr;
	float* ftmpx = NULL;
	float* ftmpy = NULL;
	
	if (nitems > 0) {
		nx2 = (m_nscx - (m_nscx % 2)) >> 1;
		nbytes = sizeof(int)*(size_t)m_nscx;
		if (0<AllocMem_h((void**)&m_h_fnx, nbytes, "InitCore", "horizontal frequency index")) { nerr = 7; goto _Exit; } // failed to allocate horizontal frequency index helpers
		for (i = 0; i < m_nscx; i++) { // setup horizontal index helper
			m_h_fnx[i] = ((i + nx2) % m_nscx) - nx2;
		}
		ny2 = (m_nscy - (m_nscy % 2)) >> 1;
		nbytes = sizeof(int)*(size_t)m_nscy;
		if (0<AllocMem_h((void**)&m_h_fny, nbytes, "InitCore", "vertical frequency index")) { nerr = 8; goto _Exit; } // failed to allocate vertical frequency index helpers
		for (j = 0; j < m_nscy; j++) {
			m_h_fny[j] = ((j + ny2) % m_nscy) - ny2;
		}
	}
	
	if (whichcode & (int)_JMS_CODE_CPU && (m_status_setup_CPU & (int)_JMS_THRESHOLD_CORE) > 0 && nitems > 0 ) {
		if (NULL != m_jcpuco && m_ncputhreads>0 ) { // there seem to be old cores, get rid of them
			for (icore = 0; icore < m_ncputhreads; icore++) {
				m_jcpuco[icore].Deinit();
			}
			delete[] m_jcpuco;
			m_jcpuco = NULL;
			m_ncputhreads = 0;
		}
		if (nCPUthreads > m_threads_CPU_out) { nerr = 1; goto _Exit; }
		m_ncputhreads = __max(1, nCPUthreads); // new number of cores
		//m_jcpuco = new CJFFTWcore[m_ncputhreads]; // allocate new FFTW cores
		m_jcpuco = new CJFFTMKLcore[m_ncputhreads]; // allocate new FFTMKL cores
		if (NULL == m_jcpuco) {
			std::cerr << "Error: (InitCore): Failed to create " << m_ncputhreads << " new CPU cores" << std::endl;
			m_ncputhreads = 0;
			nerr = 2; goto _Exit; 
		} // failed to allocate FFTW objects
		if (0 < AllocMem_h((void**)&m_status_calc_CPU, sizeof(int)*m_ncputhreads, "InitCore", "CPU calculation status", true)) { nerr = 10; goto _Exit; }
		for (icore = 0; icore < m_ncputhreads; icore++) { // initialize all FFTW cores
			//if (0 == m_jcpuco[icore].Init(2, pdims, (int)_JMS_FFTW_PLANFLAG)) {
			if (0 == m_jcpuco[icore].Init(2, pdims)) {
				ncore_ready++;
			}
			else {
				std::cerr << "Error: (InitCore): Failed to initialize CPU core #" << icore+1 << std::endl;
			}
		}
		if (ncore_ready < m_ncputhreads) { nerr = 3; goto _Exit; } // not all cores are initialized
		nbytes = sizeof(fcmplx)*(size_t)nitems;
		if (0 < AllocMem_h((void**)&m_h_wav0, nbytes, "InitCore", "wave function backup", true)) { nerr = 4; goto _Exit; }
		nbytes = sizeof(fcmplx)*(size_t)nitems*m_ncputhreads;
		if (0 < AllocMem_h((void**)&m_h_wav, nbytes, "InitCore", "wave functions", true)) { nerr = 5; goto _Exit; }
		nbytes = sizeof(int)*(size_t)m_ncputhreads;
		if (0 < AllocMem_h((void**)&m_h_dif_ndescanx, nbytes, "InitCore", "de-scan x", true)) { nerr = 6; goto _Exit; }
		if (0 < AllocMem_h((void**)&m_h_dif_ndescany, nbytes, "InitCore", "de-scan y", true)) { nerr = 7; goto _Exit; }
		m_status_setup_CPU |= (int)_JMS_STATUS_CORE; // mark CPU core setup as completed
	}

	if (whichcode & (int)_JMS_CODE_GPU && (m_status_setup_GPU & (int)_JMS_THRESHOLD_CORE) > 0 && nitems > 0) {
		m_jgpuco.Deinit();
		if (0 < m_jgpuco.Init(2, pdims)) {
			std::cerr << "Error: (InitCore): Failed to initialize GPU core." << std::endl;
			nerr = 103; goto _Exit; 
		} // gpu core init failure?
		nbytes = sizeof(cuComplex)*(size_t)nitems;
		if (0 < AllocMem_d((void**)&m_d_wav0, nbytes, "InitCore", "wave function backup", true)) { nerr = 104; goto _Exit; }
		if (0 < AllocMem_d((void**)&m_d_wav, nbytes, "InitCore", "wave function", true)) { nerr = 105; goto _Exit; }
		if (0 < AllocMem_d((void**)&m_d_det_tmpwav, nbytes, "InitCore", "temp wave readout", true)) { nerr = 111; goto _Exit; }
		nbytes = sizeof(float)*(size_t)nitems;
		if (0 < AllocMem_d((void**)&m_d_det_tmp, nbytes, "InitCore", "temp readout", true)) { nerr = 106; goto _Exit; }
		if (0 < AllocMem_d((void**)&m_d_knx, nbytes, "InitCore", "kx field", true)) { nerr = 107; goto _Exit; }
		if (0 < AllocMem_d((void**)&m_d_kny, nbytes, "InitCore", "ky field", true)) { nerr = 108; goto _Exit; }
		ftmpx = (float*)malloc(nbytes);
		ftmpy = (float*)malloc(nbytes);
		for (j = 0; j < m_nscy; j++) {
			ky = kscay * m_h_fny[j];
			i1 = j*m_nscx;
			for (i = 0; i < m_nscx; i++) {
				kx = kscax * m_h_fnx[i];
				ftmpx[i + i1] = kx;
				ftmpy[i + i1] = ky;
			}
		}
		cuerr = cudaMemcpy(m_d_knx, ftmpx, nbytes, cudaMemcpyHostToDevice);
		if (cuerr != cudaSuccess) {
			PostCUDAError("(InitCore): Failed to copy kx data to device.", cuerr);
			nerr = 109; goto _Exit;
		}
		cuerr = cudaMemcpy(m_d_kny, ftmpy, nbytes, cudaMemcpyHostToDevice);
		if (cuerr != cudaSuccess) {
			PostCUDAError("(InitCore): Failed to copy ky data to device.", cuerr);
			nerr = 110; goto _Exit;
		}
		m_status_setup_GPU |= (int)_JMS_STATUS_CORE; // mark GPU core setup as completed
	}
_Exit:
	if (NULL != ftmpx) { free(ftmpx); ftmpx = NULL; }
	if (NULL != ftmpy) { free(ftmpy); ftmpx = NULL; }
	return nerr;
}


int CJMultiSlice::SetIncidentWave(int whichcode, fcmplx* wav, bool bTranspose)
{
	int nerr = 0;
	cudaError cuerr;
	int nitems = m_nscx*m_nscy;
	int i = 0, j = 0, i1 = 0, j1 = 0, idx = 0, idx1 = 0, idy = 0, idy1 = 0; // iterator
	size_t nbytes = sizeof(fcmplx)*(size_t)nitems;
	fcmplx* wavuse = wav;
	fcmplx* wavtmp = NULL; // temp buffer only used when bTranspose == true
	// float ftmp1 = 0.f, ftmp2 = 0.f;
	// ftmp1 = GetAbsTotal(wav, nitems);
	if (bTranspose) { // transpose the input to wavtmp and relink wavuse
		wavtmp = (fcmplx*)malloc(nbytes);
		if (NULL == wavtmp) { // allocation failed
			nerr = 1; goto _Exit;
		}
		// copy a transposed version of the input to the temporary buffer
		// - assumes that the input wav has m_nscx rows of length m_nscy
		// - writes data transposed to wavtmp in m_nscy rows of length m_nscx
		// This can be slow due to massive memory jumps
		for (j = 0; j < m_nscx; j++) { // loop over input m_nscx rows -> tmp columns (j)
			idy = j * m_nscy;
			for (i = 0; i < m_nscy; i++) { // loop over input m_nscy columns -> tmp rows (i)
				idx = i + idy;
				idx1 = j + i * m_nscx;
				wavtmp[idx1] = wav[idx];
			}
		}
		wavuse = wavtmp; // link the temporary buffer to used
	}
	if (whichcode & (int)_JMS_CODE_CPU && (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC) > 0 &&
		nitems > 0 && m_h_wav0 != NULL && wavuse != NULL) {
		if (NULL == memcpy(m_h_wav0, wavuse, nbytes)) {
			nerr = 2; goto _Exit;
		}
		// debug
		// ftmp2 = GetAbsTotal(m_h_wav0, nitems);
		// end debug
	}
	if (whichcode & (int)_JMS_CODE_GPU && (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC) > 0 &&
		nitems > 0 && m_d_wav0 != NULL && wavuse != NULL) {
		cuerr = cudaMemcpy(m_d_wav0, wavuse, nbytes, cudaMemcpyHostToDevice);
		if (cuerr != cudaSuccess) {
			PostCUDAError("(SetIncidentWave): Failed to copy wave function to GPU devices", cuerr);
			nerr = 102; goto _Exit;
		}
	}
_Exit:
	if (NULL != wavtmp) { free(wavtmp); wavtmp = NULL; } // free wavtmp if used
	return nerr;
}


int CJMultiSlice::GetUnscrambleHash(unsigned int* phash)
{
	int nerr = 0;
	int nx2 = 0, ny2 = 0;
	int i = 0, j = 0, idy = 0, i1 = 0, j1 = 0, idy1 = 0;
	unsigned int idx = 0, idx1 = 0;
	if (NULL == phash) { // I/O buffer needs to be allocated
		nerr = 1; goto _exit_point;
	}
	if ((0 < (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC)) || (0 < (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC))) { // core needs to be initialized
		nx2 = (m_nscx - (m_nscx % 2)) >> 1;
		ny2 = (m_nscy - (m_nscy % 2)) >> 1;
		for (j = 0; j < m_nscy; j++) {
			idy = j * m_nscx;
			j1 = (m_h_fny[j] + ny2) % m_nscy;
			idy1 = j1 * m_nscx;
			for (i = 0; i < m_nscx; i++) {
				idx = i + idy;
				i1 = (m_h_fnx[i] + nx2) % m_nscx;
				idx1 = i1 + idy1;
				phash[idx1] = idx; // a_unscrambled[idx] = a_scrambled[phash[idx]]
			}
		}
	}
_exit_point:
	return nerr;
}


int CJMultiSlice::OffsetIncomingWave(int whichcode, float dx, float dy, float dz, int iThread)
{
	// Requires completed setup: (int)_JMS_THRESHOLD_CALC
	// Assumes that m_h_wav0 and m_d_wav0 is set in Fourier space representation
	// Leaves an offset wave in m_h_wav and m_d_wav
	int nerr = 0;
	int i = 0, j = 0, jj = 0, idx = 0;
	int nitems = m_nscx*m_nscy;
	float ca = 1.f / m_a0[0];
	float cb = 1.f / m_a0[1];
	float chi = 1.f;
	float kx = 0.f, ky = 0.f, ppy = 0.f, fz = dz*m_wl*0.5f;

	size_t nbytes = sizeof(fcmplx)*(size_t)nitems;
	fcmplx * _h_wav = NULL;
	fcmplx * _h_pp = NULL;
	cuComplex * _d_pp = NULL;
	cudaError cuerr;
	ArrayOpStats1 stats;

	// float ftmp1 = 0.f, ftmp2 = 0.f;

	if (nitems <= 0) { // prepare phase plate Exp[ - 2*Pi*I ( kx*dx + ky*dy + (kx*kx + ky*ky)*dz*Lambda/2 ) ] 
		nerr = 1; goto _Exit; // invalid number of grid points
	}

	// ftmp1 = GetAbsTotal(m_h_wav0, nitems);

	if (whichcode & (int)_JMS_CODE_CPU && (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC) > 0 &&
		iThread >= 0 && iThread < m_ncputhreads) {
		_h_wav = m_h_wav + iThread*nitems; // get thread wave slot
		if (0<AllocMem_h((void**)&_h_pp, nbytes, "OffsetIncomingWave", "phase plate")) { nerr = 2; goto _Exit; } // allocation of helper phase-plate failed
		
		for (j = 0; j < m_nscy; j++) {
			jj = j*m_nscx;
			ky = (float)m_h_fny[j] * cb;
			ppy = ky*(dy + ky*fz);
			for (i = 0; i < m_nscx; i++) {
				idx = i + jj;
				kx = (float)m_h_fnx[i] * ca;
				chi = (float)(_TPI*(ppy + kx*(dx + kx*fz)));
				_h_pp[idx].real(cos(chi)); // cos(chi)
				_h_pp[idx].imag(-sin(chi)); // - I * Sin(chi)
			}
		}
		// multiply the translation phase plate to the original incoming wave and store in wave slot
		for (i = 0; i < nitems; i++) {
			_h_wav[i] = m_h_wav0[i] * _h_pp[i];
		}
		// ftmp2 = GetAbsTotal(_h_wav, nitems);
	}

	if (whichcode & (int)_JMS_CODE_GPU && (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC) > 0 ) {
		// multiply the translation phase plate to the original incoming wave and store in wave slot
		//if (0<AllocMem_d((void**)&_d_pp, nbytes, "OffsetIncomingWave", "phase plate")) { nerr = 103; goto _Exit; } // allocation of helper phase-plate failed
		//cuerr = cudaMemcpy(_d_pp, _h_pp, nbytes, cudaMemcpyHostToDevice);
		//if (cuerr != cudaSuccess) {
		//	PostCUDAError("(OffsetIncomingWave): failed to copy phase-plate to device", cuerr);
		//	nerr = 104; goto _Exit;
		//}
		stats.uSize = nitems;
		stats.nBlockSize = m_jgpuco.GetBlockSize();
		if (stats.nBlockSize == 0) stats.nBlockSize = 512; // some common sense block size in case of failed setup?
		stats.nGridSize = (nitems + stats.nBlockSize - 1) / stats.nBlockSize;
		cuerr = ArrayOpMulPP01(m_d_wav, m_d_wav0, m_d_knx, m_d_kny, dx, dy, fz, stats);
		//cuerr = ArrayOpMul(m_d_wav, m_d_wav0, _d_pp, stats);
		if (cuerr != cudaSuccess) {
			PostCUDAError("(OffsetIncomingWave): failed to multiply phase-plate to wave function on device", cuerr);
			nerr = 105; goto _Exit;
		}
	}

_Exit:
	DeallocMem_h((void**)&_h_pp);
	DeallocMem_d((void**)&_d_pp);
	return nerr;
}


// There is something fishy with this version in the GPU code. The power of the wavefunction
// decreases afer re-inserting the backup. Check it, before using this again!
//
//int CJMultiSlice::DescanDif(int whichcode, float dx, float dy, float* dif, int iThread)
//{
//	// Requires completed setup: (int)_JMS_THRESHOLD_CALC
//	// Assumes that m_h_wav[iThread] and m_d_wav is set in Fourier space representation
//	// Leaves a shifted wave in m_h_wav[iThread] and m_d_wav
//	int nerr = 0;
//	int i = 0, j = 0, jj = 0, idx = 0;
//	int nitems = m_nscx * m_nscy;
//	float ca = m_a0[0] / m_nscx; // x steps [nm/pixel]
//	float cb = m_a0[1] / m_nscy; // y steps [nm/pixel]
//	float chi = 1.f;
//	float fscale = 1.f;
//	float rx = 0.f, ry = 0.f, ppy = 0.f;
//	float fdescan = sqrtf(dx * dx + dy * dy);
//	size_t nbytes = sizeof(fcmplx)*(size_t)nitems;
//	fcmplx* _h_pp = NULL;
//	fcmplx* _h_wav_bk = NULL; // wave function backup on host
//	cuComplex* _d_wav_bk = NULL; // wave function backup on device
//
//	if (0 == m_dif_descan_flg) { goto _Exit; } // descan not activated, skip
//	if ((0 < (whichcode & (int)_JMS_CODE_CPU)) && (0 < (whichcode & (int)_JMS_CODE_GPU))) {
//		nerr = 666;	goto _Exit; // invalid whichcode
//	}
//	if ((0 == (whichcode & (int)_JMS_CODE_CPU)) && (0 == (whichcode & (int)_JMS_CODE_GPU))) {
//		nerr = 667;	goto _Exit; // invalid whichcode
//	}
//
//	if (nitems <= 0 || NULL == dif) {
//		nerr = 1; goto _Exit; // invalid number of grid points or input array
//	}
//	fscale = 1.f / (float)nitems; // rescaling for preserving wave function norm
//	if (fdescan == 0.f) { // zero descan speed up
//		if (whichcode & (int)_JMS_CODE_CPU && (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC) > 0 &&
//			iThread >= 0 && iThread < m_nfftwthreads) {
//			m_jcpuco[iThread].GetDataPow(dif);
//		}
//		if (whichcode & (int)_JMS_CODE_GPU && (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC) > 0) {
//			m_jgpuco.GetDataPow_d(dif);
//		}
//		goto _Exit;
//	}
//	// prepare a phase ramp for shifting in diffraction space
//	// Exp[ -I * 2*Pi * (x * dx + y * dy) ]
//	if (0<AllocMem_h((void**)&_h_pp, nbytes, "DescanDif", "phase plate")) { nerr = 2; goto _Exit; } // allocation of helper phase-plate failed
//	for (j = 0; j < m_nscy; j++) { // loop over rows
//		jj = j * m_nscx; // row offset
//		ry = (float)m_h_fny[j] * cb; // y - coodinate [nm]
//		ppy = ry * dy;
//		for (i = 0; i < m_nscx; i++) {
//			idx = i + jj;
//			rx = (float)m_h_fnx[i] * ca;
//			chi = (float)(_TPI*(ppy + rx * dx));
//			_h_pp[idx] = fcmplx(fscale * cos(chi), -fscale * sin(chi)); // fscale*Exp[-I*chi]
//		}
//	}
//	if (whichcode & (int)_JMS_CODE_CPU && (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC) > 0 &&
//		iThread >=0 && iThread < m_nfftwthreads) {
//		AllocMem_h((void**)&_h_wav_bk, nbytes, "DescanDif", "_h_wav_bk", false);
//		CJFFTWcore* jco = &m_jcpuco[iThread]; // link to the FFT core of the thread
//		jco->GetDataC(_h_wav_bk); // backup the current wave function
//		jco->IFT(); // transform the wave function to real space
//		jco->MultiplyC(_h_pp); // multiply the translation phase plate
//		jco->FT(); // transform back to Fourier-space
//		jco->GetDataPow(dif); // real part -> host memory diffraction pattern 
//		jco->SetDataC(_h_wav_bk); // reset the core data to the backup
//	}
//	if (whichcode & (int)_JMS_CODE_GPU && (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC) > 0) {
//		AllocMem_d((void**)&_d_wav_bk, nbytes, "DescanDif", "_d_wav_bk", false);
//		nbytes = sizeof(cuComplex)*(size_t)nitems;
//		CJFFTCUDAcore* jco = &m_jgpuco; // link the FFT CUDA core
//		float fpow[7] = { 0.f,0.f,0.f,0.f,0.f,0.f };
//		jco->GetDataTotalPow(fpow[0]);
//		jco->GetDataC_d(_d_wav_bk); // backup the current wave function on device
//		jco->IFT(); // transform the wave function to real space
//		jco->GetDataTotalPow(fpow[1]);
//		jco->MultiplyC(_h_pp); // multiply the translation phase plate
//		jco->GetDataTotalPow(fpow[2]);
//		jco->FT(); // transform back to Fourier-space
//		jco->GetDataTotalPow(fpow[3]);
//		jco->GetDataPow_d(dif); // real part -> device memory diffraction pattern 
//		jco->GetDataTotalPow(fpow[4]);
//		jco->SetDataC_d(_d_wav_bk); // reset the core data to the backup on device
//		jco->GetDataTotalPow(fpow[5]);
//		fpow[6] = fpow[5];
//	}
//
//_Exit:
//	DeallocMem_h((void**)&_h_pp);
//	DeallocMem_h((void**)&_h_wav_bk);
//	DeallocMem_d((void**)&_d_wav_bk);
//	return nerr;
//}


int CJMultiSlice::DescanDifN(int whichcode, int dnx, int dny, float* dif, int iThread)
{
	// Requires completed setup: (int)_JMS_THRESHOLD_CALC
	// Assumes that m_h_wav[iThread] and m_d_wav is set in Fourier space representation
	// Leaves a shifted wave in m_h_wav[iThread] and m_d_wav
	int nerr = 0;
	int i = 0, j = 0, jj = 0, idx = 0;
	int j1 = 0, jj1 = 0, idx1 = 0;
	size_t nitems = m_nscx * m_nscy;
	size_t nbytes = 0;
	cudaError_t cuerr = cudaSuccess;
	float* _h_desc = NULL;
	float* _h_cp = NULL;
	float dsum = 0.f;
	//
	// check for invalid parameters
	if ((0 < (whichcode & (int)_JMS_CODE_CPU)) && (0 < (whichcode & (int)_JMS_CODE_GPU))) {
		nerr = 666;	goto _Exit; // invalid whichcode
	}
	if ((0 == (whichcode & (int)_JMS_CODE_CPU)) && (0 == (whichcode & (int)_JMS_CODE_GPU))) {
		nerr = 667;	goto _Exit; // invalid whichcode
	}
	if (nitems <= 0 || NULL == dif) {
		nerr = 1; goto _Exit; // invalid number of grid points or input array
	}
	//
	if ((dnx == 0 && dny == 0) || (0 == m_dif_descan_flg)) { // zero descan or no-descan flag -> speed up
		if (whichcode & (int)_JMS_CODE_CPU && (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC) > 0 &&
			iThread >= 0 && iThread < m_ncputhreads) {
			m_jcpuco[iThread].GetDataPow(dif); // direct readout to host buffer
		}
		if (whichcode & (int)_JMS_CODE_GPU && (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC) > 0) {
			m_jgpuco.GetDataPow_d(dif); // direct readout to device buffer
		}
		goto _Exit;
	}
	//
	// prepare a host buffer receiving the shifted diffraction pattern
	nbytes = sizeof(float)*((size_t)nitems);
	if (0<AllocMem_h((void**)&_h_cp, nbytes, "DescanDif", "_h_cp")) { nerr = 2; goto _Exit; } // allocation of array failed
	if (0<AllocMem_h((void**)&_h_desc, nbytes, "DescanDif", "_h_desc")) { nerr = 3; goto _Exit; } // allocation of array failed
	// fill data from cores
	if (whichcode & (int)_JMS_CODE_CPU && (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC) > 0 &&
		iThread >= 0 && iThread < m_ncputhreads) {
		nerr = m_jcpuco[iThread].GetDataPow(_h_cp); // direct readout to host buffer
		dsum = GetTotalF(_h_cp, nitems);
	}
	if (whichcode & (int)_JMS_CODE_GPU && (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC) > 0) {
		nerr = m_jgpuco.GetDataPow(_h_cp); // direct readout to host buffer
		dsum = GetTotalF(_h_cp, nitems);
	}
	if (nerr > 0) {
		sprintf_s(m_msg, "(CJMultiSlice::DescanDifN): Failed to readout diffraction power (code type: %d, thread: %d, error code: %d)", whichcode, iThread, nerr);
		nerr = 4; goto _Exit;
	}
	// pixel-shift de-scan loop with periodic wrap-around
	for (j = 0; j < m_nscy; j++) { // loop over rows
		j1 = mod(j + dny, m_nscy); // shifted row index // use inline mod since numbers can be negative
		jj = j * m_nscx; // row offset
		jj1 = j1 * m_nscx; // shifted row offset
		for (i = 0; i < m_nscx; i++) {
			idx = i + jj; // index
			idx1 = mod(i + dnx, m_nscx) + jj1; // shifted index // use inline mod since numbers can be negative
			_h_desc[idx1] = _h_cp[idx]; // copy shifted
		}
	}
	if (whichcode & (int)_JMS_CODE_CPU && (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC) > 0 &&
		iThread >= 0 && iThread < m_ncputhreads) {
		dsum = GetTotalF(_h_desc, nitems);
		if (dif != memcpy(dif, _h_desc, nbytes)) { // copy to host output buffer
			nerr = 5; goto _Exit; // result copy error
		}
	}
	if (whichcode & (int)_JMS_CODE_GPU && (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC) > 0) {
		dsum = GetTotalF(_h_desc, nitems);
		cuerr = cudaMemcpy(dif, _h_desc, nbytes, ::cudaMemcpyHostToDevice);
		if (cuerr != cudaSuccess) {
			nerr = 5; goto _Exit; // result copy error
		}
	}

_Exit:
	DeallocMem_h((void**)&_h_cp);
	DeallocMem_h((void**)&_h_desc);
	return nerr;
}



int CJMultiSlice::DescanDifWavN_h(int dnx, int dny, fcmplx* wav, int iThread)
{
	// Requires completed setup: (int)_JMS_THRESHOLD_CALC
	// Assumes that m_h_wav[iThread] is set in Fourier space representation
	// Leaves the original wave in m_h_wav[iThread] and ouputs shifted wave to wav
	int nerr = 0;
	int i = 0, j = 0, jj = 0, idx = 0;
	int j1 = 0, jj1 = 0, idx1 = 0;
	size_t nitems = m_nscx * m_nscy;
	size_t nbytes = 0;
	fcmplx* _h_cp = NULL;
	//
	// check for invalid parameters
	if (nitems <= 0 || NULL == wav) {
		nerr = 1; goto _Exit; // invalid number of grid points or input array
	}
	if (0 == (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC) || 0 > iThread || iThread >= m_ncputhreads) {
		nerr = 1; goto _Exit; // cpu core is not ready or invalid thread ID
	}
	//
	if ((dnx == 0 && dny == 0) || (0 == m_dif_descan_flg)) { // zero descan or no-descan flag -> speed up
		m_jcpuco[iThread].GetDataC(wav); // direct readout to host buffer
		goto _Exit;
	}
	//
	// prepare a host buffer receiving the shifted diffraction pattern
	nbytes = sizeof(fcmplx)*((size_t)nitems);
	if (0<AllocMem_h((void**)&_h_cp, nbytes, "DescanDifWavN_h", "_h_cp")) { nerr = 2; goto _Exit; } // allocation of array failed
	// fill data from cores
	nerr = m_jcpuco[iThread].GetDataC(_h_cp); // direct readout to host buffer
	if (nerr > 0) {
		sprintf_s(m_msg, "(CJMultiSlice::DescanDifWavN_h): Failed to readout diffraction power (thread: %d, error code: %d)", iThread, nerr);
		nerr = 4; goto _Exit;
	}
	// pixel-shift de-scan loop with periodic wrap-around
	for (j = 0; j < m_nscy; j++) { // loop over rows
		j1 = mod(j + dny, m_nscy); // shifted row index // use inline mod since numbers can be negative
		jj = j * m_nscx; // row offset
		jj1 = j1 * m_nscx; // shifted row offset
		for (i = 0; i < m_nscx; i++) {
			idx = i + jj; // index
			idx1 = mod(i + dnx, m_nscx) + jj1; // shifted index // use inline mod since numbers can be negative
			wav[idx1] = _h_cp[idx]; // copy shifted
		}
	}
_Exit:
	DeallocMem_h((void**)&_h_cp);
	return nerr;
}



int CJMultiSlice::DescanDifWavN_d(int dnx, int dny, cuComplex* wav)
{
	// Requires completed setup: (int)_JMS_THRESHOLD_CALC
	// Assumes that m_d_wav is set in Fourier space representation
	// Leaves a original wave in m_d_wav and shifted wave in wav
	int nerr = 0;
	int i = 0, j = 0, jj = 0, idx = 0;
	int j1 = 0, jj1 = 0, idx1 = 0;
	size_t nitems = m_nscx * m_nscy;
	size_t nbytes = 0;
	cudaError_t cuerr = cudaSuccess;
	fcmplx* _h_desc = NULL;
	fcmplx* _h_cp = NULL;
	float dsum = 0.f;
	//
	// check for invalid parameters
	if (nitems <= 0 || NULL == wav) {
		nerr = 1; goto _Exit; // invalid number of grid points or input array
	}
	if (0 == (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC)) {
		nerr = 1; goto _Exit; // gpu core is not ready
	}
	//
	if ((dnx == 0 && dny == 0) || (0 == m_dif_descan_flg)) { // zero descan or no-descan flag -> speed up
		m_jgpuco.GetDataC_d(wav); // direct readout to device buffer
		goto _Exit;
	}
	//
	// prepare a host buffer receiving the shifted diffraction pattern
	nbytes = sizeof(fcmplx)*((size_t)nitems);
	if (0<AllocMem_h((void**)&_h_cp, nbytes, "DescanDifWavN_d", "_h_cp")) { nerr = 2; goto _Exit; } // allocation of array failed
	if (0<AllocMem_h((void**)&_h_desc, nbytes, "DescanDifWavN_d", "_h_desc")) { nerr = 3; goto _Exit; } // allocation of array failed
	// fill data from core
	nerr = m_jgpuco.GetDataC(_h_cp); // direct readout to host buffer
	if (nerr > 0) {
		sprintf_s(m_msg, "(CJMultiSlice::DescanDifWavN_d): Failed to readout wave function from device (error code: %d)", nerr);
		nerr = 4; goto _Exit;
	}
	// pixel-shift de-scan loop with periodic wrap-around
	for (j = 0; j < m_nscy; j++) { // loop over rows
		j1 = mod(j + dny, m_nscy); // shifted row index // use inline mod since numbers can be negative
		jj = j * m_nscx; // row offset
		jj1 = j1 * m_nscx; // shifted row offset
		for (i = 0; i < m_nscx; i++) {
			idx = i + jj; // index
			idx1 = mod(i + dnx, m_nscx) + jj1; // shifted index // use inline mod since numbers can be negative
			_h_desc[idx1] = _h_cp[idx]; // copy shifted
		}
	}
	cuerr = cudaMemcpy(wav, _h_desc, nbytes, ::cudaMemcpyHostToDevice);
	if (cuerr != cudaSuccess) {
		nerr = 5; goto _Exit; // result copy error
	}
_Exit:
	DeallocMem_h((void**)&_h_cp);
	DeallocMem_h((void**)&_h_desc);
	return nerr;
}




int CJMultiSlice::GetResult(int whichcode, int whichresult, float *dst, int iThread)
{
	int nerr = 0;
	size_t nbytes = 0;
	size_t nitems = 0;
	float *detcur = 0;
	cudaError cuerr;

	if (NULL == dst) { // invalid destination, do not even try
		goto _Exit;
	}

	if ((whichcode & (int)_JMS_CODE_CPU) > 0 && iThread >=0 && iThread < m_ncputhreads && m_ndetslc>0 ) {
		if (whichresult == (int)_JMS_DETECT_INTEGRATED) {
			nitems = (size_t)m_ndet*m_ndetslc;
			nbytes = sizeof(float)*(size_t)nitems;
			detcur = m_h_det_int + iThread * nitems;
			memcpy(dst, detcur, nbytes);
		}
		if (whichresult == (int)_JMS_DETECT_IMAGE) {
			nitems = (size_t)m_nscx*m_nscy*m_ndetslc;
			nbytes = sizeof(float)*(size_t)nitems;
			detcur = m_h_det_img + iThread * nitems;
			memcpy(dst, detcur, nbytes);
		}
		if (whichresult == (int)_JMS_DETECT_DIFFRACTION) {
			nitems = (size_t)m_nscx*m_nscy*m_ndetslc;
			nbytes = sizeof(float)*(size_t)nitems;
			detcur = m_h_det_dif + iThread * nitems;
			memcpy(dst, detcur, nbytes);
		}
	}

	if ((whichcode & (int)_JMS_CODE_GPU) > 0 && m_ndetslc>0) {
		if (whichresult == (int)_JMS_DETECT_INTEGRATED) {
			nitems = (size_t)m_ndet*m_ndetslc;
			nbytes = sizeof(float)*(size_t)nitems;
			detcur = m_d_det_int;
			memcpy(dst, detcur, nbytes);
		}
		if (whichresult == (int)_JMS_DETECT_IMAGE) {
			nitems = (size_t)m_nscx*m_nscy*m_ndetslc;
			nbytes = sizeof(float)*(size_t)nitems;
			detcur = m_d_det_img;
			cuerr = cudaMemcpy(dst, detcur, nbytes, cudaMemcpyDeviceToHost);
		}
		if (whichresult == (int)_JMS_DETECT_DIFFRACTION) {
			nitems = (size_t)m_nscx*m_nscy*m_ndetslc;
			nbytes = sizeof(float)*(size_t)nitems;
			detcur = m_d_det_dif;
			cuerr = cudaMemcpy(dst, detcur, nbytes, cudaMemcpyDeviceToHost);
		}
	}

_Exit:
	return nerr;
}


int CJMultiSlice::Cleanup(void)
{
	int nerr = 0;
	int icore = 0;

	// (int)_JMS_STATUS_CORE
	DeallocMem_h((void**)&m_h_wav);
	DeallocMem_h((void**)&m_h_wav0);
	DeallocMem_d((void**)&m_d_wav);
	DeallocMem_d((void**)&m_d_wav0);
	DeallocMem_d((void**)&m_d_det_tmpwav);
	DeallocMem_d((void**)&m_d_det_tmp);
	DeallocMem_h((void**)&m_status_calc_CPU);
	if (NULL != m_jcpuco && m_ncputhreads >0) { // there seem to be old cores, get rid of them
		for (icore = 0; icore < m_ncputhreads; icore++) {
			m_jcpuco[icore].Deinit();
		}
		delete[] m_jcpuco;
		m_jcpuco = NULL;
		m_ncputhreads = 0;
	}
	m_jgpuco.Deinit();
	DeallocMem_h((void**)&m_h_fnx);
	DeallocMem_h((void**)&m_h_fny);
	DeallocMem_h((void**)&m_h_dif_ndescanx);
	DeallocMem_h((void**)&m_h_dif_ndescany);
	DeallocMem_d((void**)&m_d_knx);
	DeallocMem_d((void**)&m_d_kny);
	if ((m_status_setup_CPU & (int)_JMS_STATUS_CORE) > 0) m_status_setup_CPU -= (int)_JMS_STATUS_CORE;
	if ((m_status_setup_GPU & (int)_JMS_STATUS_CORE) > 0) m_status_setup_GPU -= (int)_JMS_STATUS_CORE;
	
	// (int)_JMS_STATUS_DET
	DeallocMem_h((void**)&m_h_det_wfr);
	DeallocMem_h((void**)&m_h_det_wff);
	DeallocMem_h((void**)&m_h_det_dif);
	DeallocMem_h((void**)&m_h_det_img);
	DeallocMem_h((void**)&m_h_det_int);
	DeallocMem_h((void**)&m_h_det);
	DeallocMem_h((void**)&m_h_detmask);
	DeallocMem_d((void**)&m_d_det_wff);
	DeallocMem_d((void**)&m_d_det_wfr);
	DeallocMem_d((void**)&m_d_det_dif);
	DeallocMem_d((void**)&m_d_det_img);
	DeallocMem_h((void**)&m_d_det_int); // ! This array is really on host memory
	DeallocMem_d((void**)&m_d_det);
	DeallocMem_d((void**)&m_d_detmask);
	DeallocMem_h((void**)&m_det_objslc);
	DeallocMem_h((void**)&m_detmask_len);
	m_ndet = 0;
	m_ndetper = 0;
	m_imagedet = (int)_JMS_DETECT_NONE;
	m_ndetslc = 0;
	m_threads_CPU_out = 0;
	if ((m_status_setup_CPU & (int)_JMS_STATUS_DET) > 0) m_status_setup_CPU -= (int)_JMS_STATUS_DET;
	if ((m_status_setup_GPU & (int)_JMS_STATUS_DET) > 0) m_status_setup_GPU -= (int)_JMS_STATUS_DET;

	// (int)_JMS_STATUS_PRO
	DeallocMem_h((void**)&m_h_pro);
	DeallocMem_d((void**)&m_d_pro);
	DeallocMem_h((void**)&m_proidx);
	m_npro = 0;
	if ((m_status_setup_CPU & (int)_JMS_STATUS_PRO) > 0) m_status_setup_CPU -= (int)_JMS_STATUS_PRO;
	if ((m_status_setup_GPU & (int)_JMS_STATUS_PRO) > 0) m_status_setup_GPU -= (int)_JMS_STATUS_PRO;

	// (int)_JMS_STATUS_OBJ
	DeallocMem_h((void**)&m_objslc);
	m_nobjslc = 0;
	if ((m_status_setup_CPU & (int)_JMS_STATUS_OBJ) > 0) m_status_setup_CPU -= (int)_JMS_STATUS_OBJ;
	if ((m_status_setup_GPU & (int)_JMS_STATUS_OBJ) > 0) m_status_setup_GPU -= (int)_JMS_STATUS_OBJ;

	// (int)_JMS_STATUS_PGR
	DeallocMem_h((void**)&m_nvarslc);
	DeallocMem_h((void**)&m_slcthick);
	DeallocMem_h((void**)&m_h_pgr);
	DeallocMem_h((void**)&m_slcoffset);
	DeallocMem_d((void**)&m_d_pgr);
	m_npgx = 0;
	m_npgy = 0;
	m_nscslc = 0;
	m_nslcvar = 0;
	if ((m_status_setup_CPU & (int)_JMS_STATUS_PGR) > 0) m_status_setup_CPU -= (int)_JMS_STATUS_PGR;
	if ((m_status_setup_GPU & (int)_JMS_STATUS_PGR) > 0) m_status_setup_GPU -= (int)_JMS_STATUS_PGR;

	return nerr;
}


void CJMultiSlice::FreeLibMem(void)
{
	//CJFFTWcore ctmp;
	CJFFTMKLcore ctmp;
	ctmp.FreeLibMem();
}

////////////////////////////////////////////////////////////////////////////////
//
// CPU CODE
//
////////////////////////////////////////////////////////////////////////////////

/*
// GetCPUNum by Dirk-Jan Kroon on 2010-06-09
int CJMultiSlice::GetCPUNum(void) {
#ifdef WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return sysinfo.dwNumberOfProcessors;
#elif MACOS
	int nm[2];
	size_t len = 4;
	uint32_t count;
	nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
	sysctl(nm, 2, &count, &len, NULL, 0);
	if (count < 1) {
		nm[1] = HW_NCPU;
		sysctl(nm, 2, &count, &len, NULL, 0);
		if (count < 1) { count = 1; }
	}
	return count;
#else
	return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
*/
int CJMultiSlice::GetCPUNum(void) {
	return (int)std::thread::hardware_concurrency();
}

int CJMultiSlice::AllocMem_h(void ** _h_a, size_t size, char* callfn, char* arrnam, bool zero)
{
	if (NULL != *_h_a) {
		free(*_h_a);
		*_h_a = NULL;
	}
	if ( (NULL == *_h_a) && (size > 0) ) {
		*_h_a = malloc(size);
		if (NULL == *_h_a) {
			std::cerr << "Error (" << callfn << "): Failed to allocate host memory (" << arrnam << ")" << std::endl;
			std::cerr << "- requested size [MB]: " << size / 1048576 << std::endl;
			return 1;
		}
		if (zero) {
			if (NULL==memset(*_h_a, 0, size)) {
				std::cerr << "Error (" << callfn << "): Failed to zeroe host memory (" << arrnam << ")" << std::endl;
				return 2;
			}
		}
	}
	return 0;
}

void CJMultiSlice::DeallocMem_h(void ** _h_a)
{
	if (NULL != *_h_a) {
		free(*_h_a);
		*_h_a = NULL;
	}
}




int CJMultiSlice::ClearDetMem_h(int iThread)
{
	int nerr = 0;
	int nitems = 0;
	size_t nbytes = 0;
	// integrating detector readouts
	if ((m_imagedet & (int)_JMS_DETECT_INTEGRATED) && m_ndet > 0 && m_ndetslc > 0 && m_h_det_int != NULL) {
		nitems = m_ndet*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		memset(&m_h_det_int[iThread*nitems], 0, nbytes);
	}
	// probe image readouts
	if ((m_imagedet & (int)_JMS_DETECT_IMAGE) > 0 && m_ndetslc > 0 && m_h_det_img != NULL) {
		nitems = m_nscx*m_nscy*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		memset(&m_h_det_img[iThread*nitems], 0, nbytes);
	}
	// probe diffraction readouts
	if ((m_imagedet & (int)_JMS_DETECT_DIFFRACTION) > 0 && m_ndetslc > 0 && m_h_det_dif != NULL) {
		nitems = m_nscx*m_nscy*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		memset(&m_h_det_dif[iThread*nitems], 0, nbytes);
	}
	// real-space wave function readouts
	if ((m_imagedet & (int)_JMS_DETECT_WAVEREAL) > 0 && m_ndetslc > 0 && m_h_det_wfr != NULL) {
		nitems = m_nscx*m_nscy*m_ndetslc;
		nbytes = sizeof(fcmplx)*nitems;
		memset(&m_h_det_wfr[iThread*nitems], 0, nbytes);
	}
	// Fourier-space wave function readouts
	if ((m_imagedet & (int)_JMS_DETECT_WAVEFOURIER) > 0 && m_ndetslc > 0 && m_h_det_wff != NULL) {
		nitems = m_nscx*m_nscy*m_ndetslc;
		nbytes = sizeof(fcmplx)*nitems;
		memset(&m_h_det_wff[iThread*nitems], 0, nbytes);
	}
//_Exit:
	// error handling
	return nerr;
}


float CJMultiSlice::GetCoreAbsTotal_h(int iThread)
{
	float ftot = 0.f;
	if (iThread >= 0 && iThread < m_ncputhreads && (0 < (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC))) {
		ftot = m_jcpuco[iThread].GetDataTotalPow();
	}
	return ftot;
}

float CJMultiSlice::GetCoreAbsTotal_d(void)
{
	float ftot = 0.f;
	if (0 < (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC)) {
		m_jgpuco.GetDataTotalPow(ftot);
	}
	return ftot;
}


float CJMultiSlice::GetAbsTotal(fcmplx* wav, size_t len)
{
	/*       Performs a Kahan sum on in_1[i]*in_2[i]                              */
	/*       https://en.wikipedia.org/wiki/Kahan_summation_algorithm              */
	double dsum = 0.0;
	double dabs = 0.0;
	double dc = 0.0;
	double dy = 0.0;
	double dt = 0.0;
	double re = 0.0;
	double im = 0.0;
	if (len > 0 && wav != NULL) {
		for (size_t i = 0; i < len; i++) {
			re = (double)wav[i].real();
			im = (double)wav[i].imag();
			dabs = re * re + im * im;
			dy = dabs - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
	}
	return (float)dsum;
}


float CJMultiSlice::GetTotalF(float* dat, size_t len)
{
	/*       Performs a Kahan sum on in_1[i]*in_2[i]                              */
	/*       https://en.wikipedia.org/wiki/Kahan_summation_algorithm              */
	double dsum = 0.0;
	double dc = 0.0;
	double dy = 0.0;
	double dt = 0.0;
	if (len > 0 && dat != NULL) {
		for (size_t i = 0; i < len; i++) {
			dy = dat[i] - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
	}
	return (float)dsum;
}

float CJMultiSlice::DotProduct_h(float *in_1, float *in_2, size_t len)
{
	float fsum = 0.0f;
	float *fdat = NULL;
	if (len > 0) {
		size_t idx = 0;
		size_t nodd = len % 2;
		size_t nlen2 = (len + nodd) / 2; // by 2 reduced length
		fdat = (float*)calloc(nlen2, sizeof(float));
		if (NULL != fdat) { // allocation successful
			for (size_t i = 0; i < len; i += 2) { // pre-multiply loop into fdat with x2 interleave (32 bit -> 64 bit)
				fdat[idx] = in_1[i] * in_2[i] + in_1[i+1] * in_2[i+1]; // products of intensity and detector sensitivity
				idx++; // increment fdat index
			}
			if (0 < nodd) { // there is one element left
				fdat[idx] = in_1[len - 1] * in_2[len - 1];
			}
			fdncs2(fdat, nlen2, &fsum); // do the butterfly sum
			free(fdat); // free temp. memory with products
		}
	}
	return fsum;
}

float CJMultiSlice::MaskedDotProduct_h(int *mask, float *in_1, float *in_2, size_t lenmask)
{
	float fsum = 0.0f;
	float *fdat = NULL;
	unsigned int imask0, imask1;
	if (lenmask > 0) {
		size_t idx = 0, i = 0;
		size_t nodd = lenmask % 2;
		size_t neve = lenmask - nodd;
		size_t nlen2 = 1 + neve / 2; // half length
		fdat = (float*)calloc(nlen2, sizeof(float));
		if (NULL != fdat) { // allocation successful
			for (i = 0; i < neve; i += 2) { // pre-multiply loop into fdat with x2 interleave (32 bit -> 64 bit)
				imask0 = (unsigned)mask[i];
				imask1 = (unsigned)mask[i + 1];
				fdat[idx] = in_1[imask0] * in_2[imask0] + in_1[imask1] * in_2[imask1]; // products of intensity and detector sensitivity
				idx++; // increment fdat index
			}
			if (0 < nodd) { // there is one element left
				imask0 = (unsigned)mask[lenmask-1];
				fdat[idx] = in_1[imask0] * in_2[imask0];
			}
			fdncs2m(fdat, nlen2, &fsum); // do the butterfly sum
			free(fdat); // free temp. memory with products
		}
	}
	return fsum;
}



fcmplx* CJMultiSlice::GetPhaseGrating_h(int iSlice, int* pVarID)
{
	// Assumes a proper setup - no sanity checks
	fcmplx * pgr = NULL;
	int jslc = 0;
	int jvar = 0;
	int nitems = m_npgx*m_npgy;
	
	jslc = m_objslc[iSlice] % m_nscslc; // get structure slice index (index modulo number of structure slices)
	if (NULL == pVarID) { // no prepared list of variant IDs
		//jvar = (int)((float)m_nvarslc[jslc] * (float)GetRand() / (RAND_MAX + 1)); // random number 0 ... m_nvarslc[jslc]-1
		jvar = (int)m_prng->unirand(0.0, (double)m_nvarslc[jslc]);
	}
	else { // take variant ID from prepared list
		jvar = pVarID[iSlice];
	}
	pgr = m_h_pgr[jslc]+jvar*nitems; // pointer to first item of slice variant # jvar
	return pgr;
}

fcmplx* CJMultiSlice::GetPropagator_h(int iSlice)
{
	fcmplx * pro = NULL;
	int jslc = 0;
	int proidx = 0;
	int nitems = m_nscx*m_nscy;
	if ((NULL != m_h_pro) && (NULL != m_objslc) && (NULL != m_proidx) && (iSlice >= 0) && (iSlice < m_nobjslc) && (nitems>0)) {
		jslc = m_objslc[iSlice] % m_nscslc; // get structure slice index (index modulo number of structure slices)
		proidx = m_proidx[jslc]; // get propagator index
		pro = m_h_pro[proidx]; // pointer to first item of propagator for structure slice # jslc
	}
	return pro;
}


int CJMultiSlice::ReadoutDifDet_h(int iSlice, int iThread, float weight)
{
	// Assumes that the data in jcpuco[iThread] is in Fourier space
	//              m_h_det_int is allocated on host
	//              m_h_det_dif is allocated on host
	int nerr = 0;
	//CJFFTWcore *jco = NULL;
	CJFFTMKLcore *jco = NULL;
	float *dif = NULL;
	float *det = NULL;
	float *out = NULL;
	fcmplx *wav = NULL;
	fcmplx *wavout = NULL;
	int *msk = NULL;
	int msklen = 0;
	int idetslc = -1;
	int nitems = m_nscx*m_nscy;
	int idx = 0;
	int idet = 0;
	if (m_ndetslc == 0) { goto _Exit; } // handle no detection
	if ((0 == m_ndet || 0 == (m_imagedet&_JMS_DETECT_INTEGRATED)) &&
		0 == (m_imagedet&_JMS_DETECT_DIFFRACTION) &&
		0 == (m_imagedet&_JMS_DETECT_WAVEFOURIER)) { goto _Exit; } // handle no diffraction detection requested
	if (NULL != m_jcpuco && NULL != m_det_objslc && (iSlice >= 0) && (iSlice < m_nobjslc+1) && (nitems > 0)) {
		idetslc = m_det_objslc[iSlice]; // get index of the slice in detection output arrays
		if (0 <= idetslc) { // there is detection registered for this slice
			if ((0 < m_ndet && 0 < (m_imagedet&_JMS_DETECT_INTEGRATED)) || 0 < (m_imagedet&_JMS_DETECT_DIFFRACTION)) { // we need a diffraction pattern
				if (0 < AllocMem_h((void**)&dif, sizeof(float)*nitems, "ReadoutDifDet_h", "diffraction pattern")) { nerr = 1; goto _Exit; }
				jco = &m_jcpuco[iThread]; // link cpu core
				if (0 < m_dif_descan_flg) { // apply diffraction descan on this thread and readout to dif
					if (0 < DescanDifN(_JMS_CODE_CPU, m_h_dif_ndescanx[iThread], m_h_dif_ndescany[iThread], dif, iThread)) { nerr = 2; goto _Exit; }
				}
				else { // readout directly, no de-scan
					jco->GetDataPow(dif); // copy the power of the current data in cpu core
				}
			}
			if (0 < m_ndet && 0 < (m_imagedet&_JMS_DETECT_INTEGRATED)) { // retrieve data for all integrating detectors
				for (idet = 0; idet < m_ndet; idet++) { // loop over detectors
					det = m_h_det[idet]; // pointer to current detector function
					msk = NULL; msklen = 0;
					if ((NULL != m_h_detmask) && (NULL != m_detmask_len)) {
						msk = m_h_detmask[idet]; // get detector mask if defined
						msklen = m_detmask_len[idet]; // get detector mask length
					}
					out = m_h_det_int + (iThread*m_ndetslc + idetslc)*m_ndet + idet; // pointer to integrating output channel
					// sum and add to current detector intensity
					if ((NULL!=msk)&&(0<msklen)) {
						*out += (weight*MaskedDotProduct_h(msk, dif, det, (size_t)msklen));
						//*out += (weight*DotProduct_h(dif, det, (size_t)nitems));
					}
					else {
						*out += (weight*DotProduct_h(dif, det, (size_t)nitems));
					}
				}
			}
			if (0 < (m_imagedet&_JMS_DETECT_DIFFRACTION)) { // retrieve data for diffraction pattern
				out = m_h_det_dif + (iThread*m_ndetslc + idetslc)*nitems; // pointer to integrating output channel
				for (idx = 0; idx < nitems; idx++) {
					out[idx] += (weight*dif[idx]); // accumulate to the output channel
				}
			}
			if (0 < (m_imagedet&_JMS_DETECT_WAVEFOURIER)) { // accumulate the current wave function
				wavout = m_h_det_wff + (iThread*m_ndetslc + idetslc)*nitems; // pointer to integrating output channel
				if (0 < AllocMem_h((void**)&wav, sizeof(fcmplx)*nitems, "ReadoutDifDet_h", "wave function")) { nerr = 3; goto _Exit; }
				if (0 < m_dif_descan_flg) { // apply diffraction descan on this thread and readout to dif
					if (0 < DescanDifWavN_h(m_h_dif_ndescanx[iThread], m_h_dif_ndescany[iThread], wav, iThread)) { nerr = 4; goto _Exit; }
				}
				else { // readout directly, no de-scan
					if (0 < jco->GetDataC(wav)) { nerr = 5; goto _Exit; };
				}
				for (idx = 0; idx < nitems; idx++) {
					wavout[idx] += (wav[idx] * weight); // accumulate to the output channel
				}
			}
			DeallocMem_h((void**)&dif);
			DeallocMem_h((void**)&wav);
		}
	}	
_Exit:
	DeallocMem_h((void**)&dif);
	return nerr;
}


int CJMultiSlice::ReadoutImgDet_h(int iSlice, int iThread, float weight)
{
	// Assumes that the data in jcpuco[iThread] is in Real space
	//              m_h_det_img is allocated on host
	int nerr = 0;
	//CJFFTWcore *jco = NULL;
	CJFFTMKLcore *jco = NULL;
	float *img = NULL;
	float *out = NULL;
	fcmplx *wav = NULL;
	fcmplx *wavout = NULL;
	int idetslc = -1;
	int nitems = m_nscx*m_nscy;
	int idx = 0;
	if (0 == m_ndetslc) { goto _Exit; } // handle no detection
	if (0 == (m_imagedet&_JMS_DETECT_IMAGE) && 
		0 == (m_imagedet&_JMS_DETECT_WAVEREAL)) { goto _Exit; } // handle no image detection requested
	if (NULL != m_jcpuco && NULL != m_det_objslc && (iSlice >= 0) && (iSlice < m_nobjslc+1) && (nitems > 0)) {
		idetslc = m_det_objslc[iSlice]; // get index of the slice in detection output arrays
		if (0 <= idetslc) { // there is detection registered for this slice
			jco = &m_jcpuco[iThread]; // link cpu core
			if (0 < (m_imagedet&_JMS_DETECT_IMAGE)) { // we need an image 
				if (0 < AllocMem_h((void**)&img, sizeof(float)*nitems, "ReadoutImgDet_h", "image pattern")) { nerr = 1; goto _Exit; }
				jco->GetDataPow(img); // copy the power of the current data in cpu core
			}
			if (0 < (m_imagedet&_JMS_DETECT_IMAGE)) { // retrieve data for image pattern // though already checked, we leave this 'if' for future extensions 
				out = m_h_det_img + (iThread*m_ndetslc + idetslc)*nitems; // pointer to integrating output channel
				for (idx = 0; idx < nitems; idx++) {
					out[idx] += (weight*img[idx]); // accumulate to the output channel
				}
			}
			if (0 < (m_imagedet&_JMS_DETECT_WAVEREAL)) { // we need a wave function
				if (0 < AllocMem_h((void**)&wav, sizeof(fcmplx)*nitems, "ReadoutImgDet_h", "wave function")) { nerr = 2; goto _Exit; }
				jco->GetDataC(wav); // copy the current data from cpu core
			}
			if (0 < (m_imagedet&_JMS_DETECT_WAVEREAL)) { // retrieve wave function data // though already checked, we leave this 'if' for future extensions 
				wavout = m_h_det_wfr + (iThread*m_ndetslc + idetslc)*nitems; // pointer to integrating output channel
				for (idx = 0; idx < nitems; idx++) {
					wavout[idx] += (wav[idx] * weight); // accumulate to the output channel
				}
			}
			DeallocMem_h((void**)&img);
			DeallocMem_h((void**)&wav);
		}
	}
_Exit:
	DeallocMem_h((void**)&img);
	return nerr;
}


int CJMultiSlice::SetIncomingWaveCPU(fcmplx* wav, bool bTranspose, int iThread)
{
	int nerr = 0;
	int nitems = m_nscx * m_nscy;
	int i = 0, j = 0, idx = 0, idx1 = 0, idy = 0; // iterator
	size_t nbytes = sizeof(fcmplx)*(size_t)nitems;
	fcmplx* wavuse = wav;
	fcmplx* _h_wav = NULL;
	fcmplx* wavtmp = NULL; // temp buffer only used when bTranspose == true
	if ( (m_status_setup_CPU & (int)_JMS_THRESHOLD_CALC) == 0 || nitems <= 0 || 
		m_h_wav == NULL || iThread < 0 || iThread >= m_ncputhreads) {
		nerr = 1; goto _Exit;
	}
	// float ftmp1 = 0.f, ftmp2 = 0.f;
	// ftmp1 = GetAbsTotal(wav, nitems);
	if (bTranspose) { // transpose the input to wavtmp and relink wavuse
		wavtmp = (fcmplx*)malloc(nbytes);
		if (NULL == wavtmp) { // allocation failed
			nerr = 1; goto _Exit;
		}
		// copy a transposed version of the input to the temporary buffer
		// - assumes that the input wav has m_nscx rows of length m_nscy
		// - writes data transposed to wavtmp in m_nscy rows of length m_nscx
		// This can be slow due to massive memory jumps
		for (j = 0; j < m_nscx; j++) { // loop over input m_nscx rows -> tmp columns (j)
			idy = j * m_nscy;
			for (i = 0; i < m_nscy; i++) { // loop over input m_nscy columns -> tmp rows (i)
				idx = i + idy;
				idx1 = j + i * m_nscx;
				wavtmp[idx1] = wav[idx];
			}
		}
		wavuse = wavtmp; // link the temporary buffer to used
	}
	_h_wav = m_h_wav + iThread * nitems; // get thread wave slot
	if (NULL == memcpy(_h_wav, wavuse, nbytes)) { // mem copy wavuse to incoming wave slot
		nerr = 2; goto _Exit;
	}
	// debug
	// ftmp2 = GetAbsTotal(_h_wav, nitems);
	// end debug
_Exit:
	if (NULL != wavtmp) { free(wavtmp); wavtmp = NULL; } // free wavtmp if used
	return nerr;
}



int CJMultiSlice::CPUMultislice(int islc0, int accmode, float weight, int iThread)
{
	int nerr = 0;
	//CJFFTWcore *jco = NULL;
	CJFFTMKLcore *jco = NULL;
	CJPlasmonMC jpl(m_jplmc); // copy plasmon parameters from template to thread instance
	jpl.SetRng(m_prng); // set rng
	fcmplx *wav = NULL;
	fcmplx *pgr = NULL;
	fcmplx *pro = NULL;
	int *var = NULL;
	int islc = islc0;
	int jslc = 0;
	int ish0 = 0, ish1 = 0; // shift indices used with plasmon scattering
	float fthick = 0.f; // slice thickness used with plasmon scattering
	int nitems = m_nscx*m_nscy;
	int npgitems = m_npgx*m_npgy;
	bool bsubframe = false;
	bool bdetect = false;
	if (nitems != npgitems) {
		bsubframe = true;
	}
	if (iThread<0 || iThread>= m_ncputhreads || NULL== m_jcpuco || NULL==m_h_wav || nitems <= 0 ) {
		//std::cout << "iThread :" << iThread << ", core: " << m_jcpuco << ", wav: " << m_h_wav << ", << n: " << nitems << std::endl;
		nerr = 1;
		goto _Exit;
	}
	if (accmode == 0) { // clear detection arrays
		if (0 < ClearDetMem_h(iThread)) {
			nerr = 2;
			goto _Exit;
		}
	}
	//
	// WARNING: No further error checks in the following code.
	//          Make sure to do all setups properly before calling this!
	//
	jco = &m_jcpuco[iThread];
	if (0 < jco->SetDataC(m_h_wav + iThread*nitems)) {
		nerr = 3; goto _Exit;
	}
	// Assuming the incident wave function in slice islc0 is prepared.
	// Assuming the current wave function is in Fourier space.
	if (islc < 0) islc = 0; // limit slice index to 0 = entrance plane
	if (islc >= m_nobjslc) islc = m_nobjslc; // limit slice index to m_nobjslc = exit plane
	
	// Readout will be at slice indices flagged in m_objslc_det only,
	// but always at the exit plane.
	// Readout will will always happen after the propagation to the next slice.
	//
	if (islc < m_nobjslc) {
		// prepare new random variant sequence
		var = (int*)malloc(sizeof(int)*m_nobjslc);
		GetRandomVariantSequence(var);
		if (m_plasmonmc_flg == 1) { // plasmon scattering
			jpl.Init();
			jpl.ResetMC();
		}
		// Default case, the current wave function is not the exit-plane wave function
		// Do the multislice
		for (jslc = islc; jslc < m_nobjslc; jslc++) {
			bdetect = (m_det_objslc[jslc] >= 0);
			// 1) readout (Fourier space)
			if (bdetect) {
				nerr = ReadoutDifDet_h(jslc, iThread, weight);
				if (nerr > 0) { nerr += 100; goto _CancelMS; }
			}
			// 2) scattering (real space)
			nerr = jco->IFT(); // inverse FFT
			if (nerr > 0) { nerr += 200; goto _CancelMS; }
			if (bdetect && ((m_imagedet&_JMS_DETECT_IMAGE)|| (m_imagedet&_JMS_DETECT_WAVEREAL))) {  // non-default real-space readout
				nerr = ReadoutImgDet_h(jslc, iThread, weight);
				if (nerr > 0) { nerr += 300; goto _CancelMS; }
			}
			pgr = GetPhaseGrating_h(jslc, var);
			if (bsubframe) { // sub-frame phase grating
				nerr = jco->MultiplySub2dC(pgr, m_npgx, m_npgy);
			}
			else { // identical phase-grating
				nerr = jco->MultiplyC(pgr);
			}
			if (nerr > 0) { nerr += 400; goto _CancelMS; }
			// 3) propagation (fourier space)
			nerr = jco->FT(); // forward FFT
			if (nerr > 0) { nerr += 500; goto _CancelMS; }
			// 3.1) plasmon scattering
			if (m_plasmonmc_flg == 1) {
				fthick = GetSliceThickness(m_objslc[jslc]);
				if (0 < jpl.ScatGridMC(fthick, m_a0[0], m_a0[1], &ish0, &ish1)) { // scattering happened
					if (ish0 != 0 || ish1 != 0) { // wave function needs shift
						jco->CShift2d(ish0, ish1);
					}
				}
			}
			pro = GetPropagator_h(jslc);
			nerr = jco->MultiplyC(pro);
			if (nerr > 0) { nerr += 600; goto _CancelMS; }
		}
	_CancelMS:
		free(var);
		var = NULL;
		if (nerr > 0) goto _Exit;
	}
	// Final readout for the exit plane (do this always even if there is no object slice)
	nerr = ReadoutDifDet_h(m_nobjslc, iThread, weight);
	if (nerr > 0) { nerr += 700; goto _Exit; }
	if ((m_imagedet&_JMS_DETECT_IMAGE)|| (m_imagedet&_JMS_DETECT_WAVEREAL)) { // non-default real-space readout
		nerr = jco->IFT(); // inverse FFT
		if (nerr > 0) { nerr += 800; goto _Exit; }
		nerr = ReadoutImgDet_h(m_nobjslc, iThread, weight);
		if (nerr > 0) { nerr += 900; goto _Exit; }
	}
	//
_Exit:
	return nerr;
}

////////////////////////////////////////////////////////////////////////////////
//
// CUDA CODE
//
////////////////////////////////////////////////////////////////////////////////

void CJMultiSlice::PostCUDAError(char* smsg, cudaError code)
{
	if (code > 0) {
		std::cerr << "Error: " << smsg << ", code: " << code << std::endl;
		std::cerr << "     - " << cudaGetErrorString(code) << std::endl;
	}
}

void CJMultiSlice::PostCUDAMemory(size_t nrequestedbytes)
{
	int64_t memavail = 0, memtotal = 0;
	int idev = 0, cc1 = 0, cc2 = 0, nmaxthread = 0;
	idev = GetCurrentGPU();
	if (idev >= 0) {
		if (0 == GetGPUStats(idev, cc1, cc2, nmaxthread, memtotal, memavail)) {
			std::cerr << "  - requested device memory [MB]: " << nrequestedbytes / 1048576;
			std::cerr << "  - available device memory [MB]: " << memavail / 1048576;
		}
	}
}

int CJMultiSlice::GetGPUNum(void)
{
	int ndev = 0;
	cudaError cuerr;
	cuerr = cudaGetDeviceCount(&ndev);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(GetGPUNum): Failed to retrieve number of GPU devices", cuerr);
		return 0;
	}
	return ndev;
}


int CJMultiSlice::GetGPUName(int idev, char* name)
{
	cudaError cuerr;
	cudaDeviceProp prop;
	cuerr = cudaGetDeviceProperties(&prop, idev);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(GetGPUName): Failed to retrieve GPU devices properties", cuerr);
		return 1;
	}
	memcpy(name, prop.name, 256);
	return 0;
}


int CJMultiSlice::GetCurrentGPU(void)
{
	int idev = -1;
	cudaError cuerr;
	cuerr = cudaGetDevice(&idev);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(GetCurrentGPU): Failed to retrieve current GPU device ID", cuerr);
		return -1;
	}
	return idev;
}


int CJMultiSlice::SetCurrentGPU(int idev)
{
	cudaError cuerr;
	cuerr = cudaSetDevice(idev);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(SetCurrentGPU): Failed to set current GPU device", cuerr);
		return 1;
	}
	return 0;
}

int CJMultiSlice::GetGPUStats(int idev, int &iCMajor, int &iCMinor, int &iMaxThread, int64_t &CUDAmemtotal, int64_t &CUDAmemfree)
{
	cudaError cuerr;
	cudaDeviceProp prop;
	cuerr = cudaGetDeviceProperties(&prop, idev);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(GetGPUStats): Failed to retrieve GPU devices properties", cuerr);
		return 1;
	}
	iCMajor = prop.major;
	iCMinor = prop.minor;
	iMaxThread = prop.maxThreadsPerBlock;
	size_t memAvail = 0, memDev = 0;
	cuerr = cudaMemGetInfo(&memAvail, &memDev);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(GetGPUStats): Failed to retrieve memory state of device", cuerr);
		return 2;
	}
	CUDAmemtotal = (int64_t)memDev;
	CUDAmemfree = (int64_t)memAvail;
	return 0;
}

int CJMultiSlice::GetGPUCores(int idev, int &nMultiProcs, int &nCores, int& nMaxThreadPerProc)
{
	cudaError cuerr;
	cudaDeviceProp prop;
	cuerr = cudaGetDeviceProperties(&prop, idev);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(GetGPUCores): Failed to retrieve GPU devices properties", cuerr);
		return 1;
	}
	nMultiProcs = prop.multiProcessorCount;
	nCores = 0;
	nMaxThreadPerProc = prop.maxThreadsPerMultiProcessor;
	switch (prop.major) {
	case 2: // Fermi
		if (prop.minor == 1) nCores = nMultiProcs * 48;
		else nCores = nMultiProcs * 32;
		break;
	case 3: // Kepler
		nCores = nMultiProcs * 192;
		break;
	case 5: // Maxwell
		nCores = nMultiProcs * 128;
		break;
	case 6: // Pascal
		if (prop.minor == 1) nCores = nMultiProcs * 128;
		else if (prop.minor == 0) nCores = nMultiProcs * 64;
		break;
	default:
		break;
	}
	return 0;
}


int CJMultiSlice::GetGPUMemInfo(size_t &memtotal, size_t &memfree)
{
	cudaError cuerr;
	cuerr = cudaMemGetInfo(&memfree, &memtotal);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(GetGPUMemInfo): Failed to retrieve memory state of device", cuerr);
		return 1;
	}
	return 0;
}


int CJMultiSlice::SetGPUPgrLoading(int npgrload)
{
	if (m_status_setup_GPU & (int)_JMS_STATUS_PGR) { // phase-gratings have been allocated already
		return 1; // return with error
	}
	m_d_pgr_src = 0;
	if (npgrload == 1) {
		m_d_pgr_src = 1;
	}
	return 0;
}


int CJMultiSlice::AllocMem_d(void ** _d_a, size_t size, char* callfn, char* arrnam, bool zero)
{
	cudaError cuerr;
	if (NULL != *_d_a) {
		cuerr = cudaFree(*_d_a);
		*_d_a = NULL;
	}
	if (size > 0) {
		cuerr = cudaMalloc(_d_a, size);
		if (cuerr != cudaSuccess) {
			sprintf_s(m_msg, "(%s): Failed to allocate device memory (%s)", callfn, arrnam);
			PostCUDAError(m_msg, cuerr);
			PostCUDAMemory(size);
			return 1;
		}
		if (zero) {
			cuerr = cudaMemset(*_d_a, 0, size);
			if (cuerr != cudaSuccess) {
				sprintf_s(m_msg, "(%s): Failed to zeroe device memory (%s)", callfn, arrnam);
				PostCUDAError(m_msg, cuerr);
				return 2;
			}
		}
	}
	return 0;
}

void CJMultiSlice::DeallocMem_d(void ** _d_a)
{
	cudaError cuerr;
	if (NULL != *_d_a) {
		cuerr = cudaFree(*_d_a);
		*_d_a = NULL;
	}
}



float CJMultiSlice::DotProduct_d(float *in_1, float *in_2, size_t len, int nBlockSize)
{
	/*       Performs a dot product, sum on in_1[i]*in_2[i]                       */
	/*       Using ArrayOps on Device                                             */
	float fsum = 0.f;
	ArrayOpStats1 stats;
	float* ftmp = NULL;
	cudaError cuerr;
	if (len < 1) goto _Exit;
	stats.uSize = (unsigned int)len;
	stats.nBlockSize = nBlockSize;
	stats.nGridSize = (int)(len + nBlockSize - 1) / nBlockSize;
	if (0 < AllocMem_d((void**)&ftmp, sizeof(float)*len, "DotProduct_d", "temp. array")) goto _Exit;
	cuerr = ArrayOpFFMul(ftmp, in_1, in_2, stats);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(DotProduct_d): Failed to multiply arrays on device (ArrayOpFFMul)", cuerr);
		goto _Exit;
	}
	cuerr = ArrayOpFSum(fsum, ftmp, stats, nBlockSize);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(DotProduct_d): Failed to summ array on device (ArrayOpFSum)", cuerr);
		goto _Exit;
	}
_Exit:
	DeallocMem_d((void**)&ftmp);
	return fsum;
}

float CJMultiSlice::MaskedDotProduct_d(int *mask, float *in_1, float *in_2, size_t lenmask, int nBlockSize)
{
	/*       Performs a dot product, sum on in_1[mask[i]]*in_2[mask[i]]           */
	/*       Using ArrayOps on Device                                             */
	float fsum = 0.f;
	ArrayOpStats1 stats;
	float* ftmp = NULL;
	cudaError cuerr;
	if (lenmask < 1) goto _Exit;
	stats.uSize = (unsigned int)lenmask;
	stats.nBlockSize = nBlockSize;
	stats.nGridSize = (int)(lenmask + nBlockSize - 1) / nBlockSize;
	cuerr = ArrayOpMaskFDot(fsum, mask, in_1, in_2, stats, nBlockSize);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(MaskedDotProduct_d): Failed to summ array on device (ArrayOpMaskFDot)", cuerr);
		goto _Exit;
	}
_Exit:
	return fsum;
}


int CJMultiSlice::ClearDetMem_d(void) {
	int nerr = 0;
	size_t nitems = 0;
	size_t nbytes = 0;
	cudaError cuerr;
	// integrating detector readouts (this is host memory, REALLY ! )
	if (m_ndet > 0 && m_ndetslc > 0 && m_d_det_int != NULL) {
		nitems = (size_t)m_ndet*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		memset(m_d_det_int, 0, nbytes);
	}
	// probe image readouts (this is device memory)
	if ((m_imagedet & (int)_JMS_DETECT_IMAGE) > 0 && m_ndetslc > 0 && m_d_det_img != NULL) {
		nitems = (size_t)m_nscx*m_nscy*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		cuerr = cudaMemset((void*)m_d_det_img, 0, nbytes);
		if (cuerr != cudaSuccess) {
			PostCUDAError("(ClearDetMem_d): Failed to reset memory m_d_det_img", cuerr);
			return 3;
		}
	}
	// probe diffraction readouts (this is device memory)
	if ((m_imagedet & (int)_JMS_DETECT_DIFFRACTION) > 0 && m_ndetslc > 0 && m_d_det_dif != NULL) {
		nitems = (size_t)m_nscx*m_nscy*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		cuerr = cudaMemset((void*)m_d_det_dif, 0, nbytes);
		if (cuerr != cudaSuccess) {
			PostCUDAError("(ClearDetMem_d): Failed to reset memory m_d_det_dif", cuerr);
			return 3;
		}
	}
	// real-space wave function readouts (this is device memory)
	if ((m_imagedet & (int)_JMS_DETECT_WAVEREAL) > 0 && m_ndetslc > 0 && m_d_det_wfr != NULL) {
		nitems = (size_t)m_nscx * m_nscy*m_ndetslc;
		nbytes = sizeof(cuComplex)*nitems;
		cuerr = cudaMemset((void*)m_d_det_wfr, 0, nbytes);
		if (cuerr != cudaSuccess) {
			PostCUDAError("(ClearDetMem_d): Failed to reset memory m_d_det_wfr", cuerr);
			return 3;
		}
	}
	// Fourier-space wave function readouts (this is device memory)
	if ((m_imagedet & (int)_JMS_DETECT_WAVEFOURIER) > 0 && m_ndetslc > 0 && m_d_det_wff != NULL) {
		nitems = (size_t)m_nscx * m_nscy*m_ndetslc;
		nbytes = sizeof(cuComplex)*nitems;
		cuerr = cudaMemset((void*)m_d_det_wff, 0, nbytes);
		if (cuerr != cudaSuccess) {
			PostCUDAError("(ClearDetMem_d): Failed to reset memory m_d_det_wff", cuerr);
			return 3;
		}
	}
	//_Exit:
	// error handling
	return nerr;
}

cuComplex* CJMultiSlice::GetPhaseGrating_d(int iSlice, int* pVarID)
{
	// Assumes proper setup and non-faulty calls - no sanity checks
	cuComplex * pgr = NULL;
	int jslc = m_objslc[iSlice] % m_nscslc; // get structure slice index (index modulo number of structure slices)
	int jvar = 0;
	int nitems = m_npgx*m_npgy;
	if (NULL == pVarID) { // no prepared list of variant IDs
		//jvar = (int)((float)m_nvarslc[jslc] * (float)GetRand() / (float)(RAND_MAX + 1)); // random number 0 ... m_nvarslc[jslc]-1
		jvar = (int)m_prng->unirand(0.0, (double)m_nvarslc[jslc]);
	}
	else { // take variant ID from prepared list
		jvar = pVarID[iSlice];
	}
	if (m_d_pgr_src == 0) { // default: phase gratings have been copied to device before
		pgr = m_d_pgr + m_slcoffset[jslc] + (int64_t)jvar*nitems; // pointer to first item of slice variant # jvar
	}
	else if (m_d_pgr_src == 1) { // v.0.15: low-GPU memory mode: phase gratings are on host and need to be copied now
		cudaError cuerr;
		fcmplx* pgrsrc = GetPhaseGrating_h(iSlice, pVarID); // get phase-grating address on host
		size_t nbytes = sizeof(fcmplx)*(size_t)nitems;
		cuerr = cudaMemcpy(m_d_pgr, pgrsrc, nbytes, cudaMemcpyHostToDevice); // copy to device
		if (cuerr != cudaSuccess) { // handle device error
			PostCUDAError("(SetPhaseGratingData): Failed to copy phase grating data to GPU devices", cuerr);
		}
		pgr = m_d_pgr; // only one phase-grating buffer on device
	}
	/*if (m_dbg >= 5) {
		std::cout << "debug(GetPhaseGrating_d): iSlice " << iSlice << std::endl;
		std::cout << "debug(GetPhaseGrating_d): jslc " << jslc << std::endl;
		std::cout << "debug(GetPhaseGrating_d): m_nvarslc[jslc] " << m_nvarslc[jslc] << std::endl;
		if (pVarID) std::cout << "debug(GetPhaseGrating_d): pVarID[iSlice] " << pVarID[iSlice] << std::endl;
		std::cout << "debug(GetPhaseGrating_d): jvar " << jvar << std::endl;
		std::cout << "debug(GetPhaseGrating_d): m_slcoffset[jslc] " << m_slcoffset[jslc] << std::endl;
		std::cout << "debug(GetPhaseGrating_d): jvar*nitems " << jvar*nitems << std::endl;
		std::cout << "debug(GetPhaseGrating_d): m_d_pgr 0x" << pgr << std::endl;
		std::cout << "debug(GetPhaseGrating_d): pgr 0x" << pgr << std::endl;
		std::cout << "debug(GetPhaseGrating_d): dif 0x" << pgr - m_d_pgr << std::endl;
	}*/
	return pgr;
}


cuComplex* CJMultiSlice::GetPropagator_d(int iSlice)
{
	// Assumes proper setup and non-faulty calls - no sanity checks
	cuComplex * pro = NULL;
	int jslc = m_objslc[iSlice] % m_nscslc; // get structure slice index (index modulo number of structure slices)
	int proidx = m_proidx[jslc]; // get propagator index
	int nitems = m_nscx*m_nscy;
	pro = m_d_pro + (int64_t)proidx*nitems; // pointer to first item of propagator for structure slice # jslc
	return pro;
}

int CJMultiSlice::ReadoutDifDet_d(int iSlice, float weight)
{
	// Assumes that the data in jgpuco is in Fourier space.
	//              m_d_det_int is allocated on host
	//              m_d_det_dif is allocated on device
	int nerr = 0;
	cudaError cuerr;
	CJFFTCUDAcore *jco = &m_jgpuco;
	ArrayOpStats1 stats;
	float *dif_d = m_d_det_tmp;
	float *det_d = NULL;
	float *out_d = NULL;
	cuComplex* wav_d = m_d_det_tmpwav;
	cuComplex* wavout_d = NULL;
	int *msk_d = NULL;
	int msklen = 0;
	int idetslc = -1;
	int nitems = m_nscx*m_nscy;
	int idx = 0;
	int idet = 0;
	if (m_ndetslc == 0) { goto _Exit; } // handle no detection planes set up
	if ((0 == m_ndet || 0 == (m_imagedet&_JMS_DETECT_INTEGRATED)) && 
		0 == (m_imagedet&_JMS_DETECT_DIFFRACTION) &&
		0 == (m_imagedet&_JMS_DETECT_WAVEFOURIER)) { goto _Exit; } // handle no diffraction detection requested
	if (NULL != m_det_objslc && (0 <= iSlice) && (iSlice < m_nobjslc+1) && (0 < nitems)) { // check valid data
		idetslc = m_det_objslc[iSlice]; // get index of the slice in detection output arrays
		stats.uSize = (unsigned int)nitems;
		stats.nBlockSize = jco->GetBlockSize();
		stats.nGridSize = (nitems + stats.nBlockSize - 1) / stats.nBlockSize;
		if (idetslc >= 0) { // there is detection registered for this slice
			if ((0 < m_ndet && 0 < (m_imagedet&_JMS_DETECT_INTEGRATED)) || (0 < (m_imagedet&_JMS_DETECT_DIFFRACTION))) { // we need a diffraction pattern ...
				if (0 < m_dif_descan_flg) { // apply diffraction descan and readout to dif_d
					if (0 < DescanDifN(_JMS_CODE_GPU, m_d_dif_ndescanx, m_d_dif_ndescany, dif_d)) { nerr = 2; goto _Exit; }
				}
				else { // no diffraction de-scan, readout directly
					jco->GetDataPow_d(dif_d); // get the power of the current data in gpu core
				}
			}
			if (0 < m_ndet && 0 < (m_imagedet&_JMS_DETECT_INTEGRATED)) { // retrieve data for all integrating detectors from the diffraction pattern
				for (idet = 0; idet < m_ndet; idet++) { // loop over detectors
					det_d = m_d_det + idet*nitems; // pointer to current detector function
					msk_d = NULL; msklen = 0;
					if ((NULL != m_d_detmask) && (NULL != m_detmask_len)) {
						msk_d = m_d_detmask + idet*nitems;
						msklen = m_detmask_len[idet];
					}
					out_d = m_d_det_int + idet + idetslc*m_ndet; // pointer to integrating output channel
					// ---> Here is the STEM detector readout!
					// sum and add to current detector intensity
					if ((NULL != msk_d) && (0 < msklen)) {
						*out_d += (weight*MaskedDotProduct_d(msk_d, dif_d, det_d, msklen, stats.nBlockSize));
						//*out_d += (weight*DotProduct_d(dif_d, det_d, nitems, stats.nBlockSize));
					}
					else {
						*out_d += (weight*DotProduct_d(dif_d, det_d, nitems, stats.nBlockSize));
					}
				}
			}
			if (0 < (m_imagedet&_JMS_DETECT_DIFFRACTION)) { // accumulate the diffraction pattern
				out_d = m_d_det_dif + idetslc*nitems; // pointer to integrating output channel
				// ---> Here is the diffraction pattern readout!
				if (fabsf(weight - 1.f) < (float)1.e-6) { // ignore weight and readout faster
					cuerr = ArrayOpFAdd0(out_d, out_d, dif_d, stats); // parallel add on GPU
				}
				else {
					cuerr = ArrayOpFAdd2(out_d, out_d, dif_d, weight, stats); // weighted parallel add on GPU, dif_d may be weighted
				}
				if (cuerr != cudaSuccess) {
					PostCUDAError("(ReadoutDifDet_d): Failed to integrate diffraction pattern on device", cuerr);
					nerr = 3;
					goto _Exit;
				}
			}
			if (0 < (m_imagedet&_JMS_DETECT_WAVEFOURIER)) { // accumulate the current wavefunction
				if (0 < m_dif_descan_flg) { // apply diffraction descan and readout to wav_d (m_d_det_tmpwav)
					if (0 < DescanDifWavN_d(m_d_dif_ndescanx, m_d_dif_ndescany, wav_d)) { nerr = 102; goto _Exit; }
				}
				else { // no diffraction de-scan, readout directly
					wav_d = jco->GetData_d(); // simply get the device pointer to the wave function (forget m_d_det_tmpwav)
					// Be careful later! The pointer wav_d is no longer m_d_det_tmpwav now.
				}
				wavout_d = m_d_det_wff + idetslc * nitems; // pointer to the integrating putput channel
				// ---> Here is the wave function readout in Fourier space!
				if (fabsf(weight - 1.f) < (float)1.e-6) { // ignore weight and readout faster
					cuerr = ArrayOpAdd0(wavout_d, wavout_d, wav_d, stats); // add from wav_d to wavout_d
				}
				else {
					cuComplex c0, c1, cw;
					c0.x = 0.f; c0.y = 0.f;
					c1.x = 1.f; c1.y = 0.f;
					cw.x = weight; cw.y = 0.f;
					cuerr = ArrayOpAdd(wavout_d, wavout_d, wav_d, c1, cw, c0, stats); // add wav_d to wavout_d with weight
				}
				if (cuerr != cudaSuccess) {
					PostCUDAError("(ReadoutDifDet_d): Failed to integrate Fourier-space wavefunction on device", cuerr);
					nerr = 3;
					goto _Exit;
				}
			}
		}
	}
	else {
		nerr = 1; // strange setup
	}
_Exit:
	return nerr;
}

int CJMultiSlice::ReadoutImgDet_d(int iSlice, float weight)
{
	// Assumes that the data in jgpuco is in Real space
	// Assumes that m_d_det_img is allocated on device.
	int nerr = 0;
	cudaError cuerr;
	CJFFTCUDAcore *jco = &m_jgpuco;
	ArrayOpStats1 stats;
	float *img_d = m_d_det_tmp;
	float *out_d = NULL;
	cuComplex *wav_d = NULL;
	cuComplex *wavout_d = NULL;
	int idetslc = -1;
	int nitems = m_nscx*m_nscy;
	int idx = 0;
	if (m_ndetslc == 0) { goto _Exit; } // handle no detection
	if (0 == (m_imagedet&_JMS_DETECT_IMAGE) && 
		0 == (m_imagedet&_JMS_DETECT_WAVEREAL)) { goto _Exit; } // handle no image detection requested
	if (NULL != m_det_objslc && (iSlice >= 0) && (iSlice < m_nobjslc+1) && (nitems > 0)) {
		idetslc = m_det_objslc[iSlice]; // get index of the slice in detection output arrays
		stats.uSize = (unsigned int)nitems;
		stats.nBlockSize = jco->GetBlockSize();
		stats.nGridSize = (nitems + stats.nBlockSize - 1) / stats.nBlockSize;
		if (0 <= idetslc) { // there is detection registered for this slice
			if (0 < (m_imagedet&_JMS_DETECT_IMAGE)) { // we need an image
				jco->GetDataPow_d(img_d); // get the power of the current data in gpu core
			}
			if (0 < (m_imagedet&_JMS_DETECT_IMAGE)) { // retrieve data for image pattern
				out_d = m_d_det_img + idetslc*nitems; // pointer to integrating output channel
				// ---> Here is the probe image readout!
				if (fabsf(weight-1.f) < (float)1.e-6) { // ignore weight
					cuerr = ArrayOpFAdd0(out_d, out_d, img_d, stats); // parallel add on GPU
				}
				else {
					cuerr = ArrayOpFAdd2(out_d, out_d, img_d, weight, stats); // weighted parallel add on GPU
				}
				if (cuerr != cudaSuccess) {
					PostCUDAError("(ReadoutImgDet_d): Failed to integrate image pattern on device", cuerr);
					nerr = 3;
					goto _Exit;
				}
			}
			if (0 < (m_imagedet&_JMS_DETECT_WAVEREAL)) { // retrieve data for wave function accumulation
				wavout_d = m_d_det_wfr + idetslc * nitems; // pointer to the integrating output channel
				// ---> Here is the wave function readout!
				wav_d = jco->GetData_d(); // link to the current wave function of the CUDA core
				if (fabsf(weight - 1.f) < (float)1.e-6) { // ignore weight
					cuerr = ArrayOpAdd0(wavout_d, wavout_d, wav_d, stats); // parallel add on GPU
				}
				else {
					cuComplex c0, c1, cw;
					c0.x = 0.f; c0.y = 0.f;
					c1.x = 1.f; c1.y = 0.f;
					cw.x = weight; cw.y = 0.f;
					cuerr = ArrayOpAdd(wavout_d, wavout_d, wav_d, c1, cw, c0, stats); // weighted parallel add on GPU
				}
				if (cuerr != cudaSuccess) {
					PostCUDAError("(ReadoutImgDet_d): Failed to integrate wave function on device", cuerr);
					nerr = 3;
					goto _Exit;
				}
			}
		}
	}
_Exit:
	return nerr;
}


int CJMultiSlice::SetIncomingWaveGPU(fcmplx* wav, bool bTranspose)
{
	int nerr = 0;
	cudaError cuerr;
	int nitems = m_nscx * m_nscy;
	int i = 0, j = 0, idx = 0, idx1 = 0, idy = 0; // iterator
	size_t nbytes = sizeof(fcmplx)*(size_t)nitems;
	fcmplx* wavuse = wav;
	fcmplx* wavtmp = NULL; // temp buffer only used when bTranspose == true
	if ( (m_status_setup_GPU & (int)_JMS_THRESHOLD_CALC) == 0 || nitems <= 0 || m_d_wav == NULL) {
		nerr = 101; goto _Exit; // module not ready
	}
	// float ftmp1 = 0.f, ftmp2 = 0.f;
	// ftmp1 = GetAbsTotal(wav, nitems);
	if (bTranspose) { // transpose the input to wavtmp and relink wavuse
		wavtmp = (fcmplx*)malloc(nbytes);
		if (NULL == wavtmp) { // allocation failed
			nerr = 1; goto _Exit;
		}
		// copy a transposed version of the input to the temporary buffer
		// - assumes that the input wav has m_nscx rows of length m_nscy
		// - writes data transposed to wavtmp in m_nscy rows of length m_nscx
		// This can be slow due to massive memory jumps
		for (j = 0; j < m_nscx; j++) { // loop over input m_nscx rows -> tmp columns (j)
			idy = j * m_nscy;
			for (i = 0; i < m_nscy; i++) { // loop over input m_nscy columns -> tmp rows (i)
				idx = i + idy;
				idx1 = j + i * m_nscx;
				wavtmp[idx1] = wav[idx];
			}
		}
		wavuse = wavtmp; // link the temporary buffer to used
	}
	cuerr = cudaMemcpy(m_d_wav, wavuse, nbytes, cudaMemcpyHostToDevice);
	if (cuerr != cudaSuccess) {
		PostCUDAError("(SetIncomingWaveGPU): Failed to copy wave function to GPU devices", cuerr);
		nerr = 102; goto _Exit;
	}
_Exit:
	if (NULL != wavtmp) { free(wavtmp); wavtmp = NULL; } // free wavtmp if used
	return nerr;
}



int CJMultiSlice::GPUMultislice(int islc0, int accmode, float weight)
{
	int nerr = 0;
	CJFFTCUDAcore *jco = &m_jgpuco;
	CJPlasmonMC jpl(m_jplmc); // copy plasmon MC template
	jpl.SetRng(m_prng); // set rng
	cuComplex *wav = NULL;
	cuComplex *pgr = NULL;
	cuComplex *pro = NULL;
	int *var = NULL;
	int islc = islc0;
	int jslc = 0;
	int ish0 = 0, ish1 = 0;
	int nitems = m_nscx*m_nscy;
	int npgitems = m_npgx*m_npgy;
	float fthick = 0.f;
	bool bsubframe = false;
	bool bdetect = false;
	if (nitems != npgitems) { // calculation grid and phase gratings are on different sizes
		bsubframe = true;
	}

	if ( NULL == m_d_wav || nitems <= 0) {
		nerr = 1;
		goto _Exit;
	}
	if (0 == jco->GetStatus()) { // CUDA FFT core not initialized
		nerr = 2;
		goto _Exit;
	}
	if (accmode == 0) { // clear detection arrays
		if (0 < ClearDetMem_d()) {
			nerr = 3;
			goto _Exit;
		}
	}
	//
	// WARNING: No further error checks in the following code.
	//          Make sure to do all setups properly before calling this!
	//
	if (0 < jco->SetDataC_d(m_d_wav)) {
		nerr = 4; goto _Exit;
	}
	// Assuming the incident wave function in slice islc0 is prepared.
	// Assuming the current wave function is in Fourier space.
	if (islc < 0) islc = 0; // limit slice index to 0 = entrance plane
	if (islc >= m_nobjslc) islc = m_nobjslc; // limit slice index to m_nobjslc = exit plane

	// Readout will be at slice indices flagged in m_objslc_det only,
	// but always at the exit plane.
	// Readout will always happen after the propagation to the next slice.
	//
	// Call DetectorSetup to generate a coherent setup for the multislice.
	// The same setup will be used by all threads.
	//
	if (islc < m_nobjslc) {
		
		// prepare new random variant sequence
		var = (int*)malloc(sizeof(int)*m_nobjslc);
		GetRandomVariantSequence(var);
		if (m_plasmonmc_flg == 1) {
			jpl.Init();
			jpl.ResetMC();
		}
		/*if (m_dbg >= 5) {
			std::cout << "debug(GPUMultislice): RandomVariantSequence" << std::endl;
			for (jslc = islc; jslc < m_nobjslc; jslc++) {
				std::cout << "debug(GPUMultislice): jslc = " << jslc << ", var[jslc] = " << var[jslc] << std::endl;
			}
		}*/
		// Default case, the current wave function is not the exit-plane wave function
		// Do the multislice
		for (jslc = islc; jslc < m_nobjslc; jslc++) {
			bdetect = (m_det_objslc[jslc] >= 0);
			// 1) readout (Fourier space)
			//if (1 == m_objslc_det[jslc]) {
			if (bdetect) {
				nerr = ReadoutDifDet_d(jslc, weight);
				if (0 < nerr) {	nerr += 100; goto _CancelMS; }
			}
			// 2) scattering (real space)
			nerr = jco->IFT(); // inverse FFT
			if (0 < nerr) { nerr += 200; goto _CancelMS; }
			if (bdetect && ((m_imagedet&_JMS_DETECT_IMAGE)||(m_imagedet&_JMS_DETECT_WAVEREAL))) {  // non-default real-space readout
				nerr = ReadoutImgDet_d(jslc, weight);
				if (0 < nerr) { nerr += 300; goto _CancelMS; }
			}
			pgr = GetPhaseGrating_d(jslc, var);
			if (bsubframe) { // different size phase grating multiplication
				nerr = jco->MultiplySub2dC_d(pgr, m_npgx, m_npgy);
			}
			else { // identical size phase grating multiplication
				nerr = jco->MultiplyC_d(pgr);
			}
			
			if (0 < nerr) { 
				/*std::cerr << jslc << std::endl;
				std::cerr << pgr << std::endl;
				std::cerr << m_d_pgr << std::endl;*/
				nerr += 400; goto _CancelMS; 
			}
			// 3) propagation (fourier space)
			nerr = jco->FT(); // forward FFT
			if (0 < nerr) { nerr += 500; goto _CancelMS; }
			// 3.1) plasmon scattering
			if (m_plasmonmc_flg == 1) {
				fthick = GetSliceThickness(m_objslc[jslc]);
				if (0 < jpl.ScatGridMC(fthick, m_a0[0], m_a0[1], &ish0, &ish1)) { // scattering happened
					if (ish0 != 0 || ish1 != 0) { // wave function needs shift
						jco->CShift2d(ish0, ish1);
					}
				}
			}
			pro = GetPropagator_d(jslc);
			nerr = jco->MultiplyC_d(pro);
			if (0 < nerr) { nerr += 600; goto _CancelMS; }
		}
	_CancelMS:
		free(var);
		var = NULL;
		if (nerr > 0) goto _Exit;
	}
	// Final readout for the exit plane (do this always)
	nerr = ReadoutDifDet_d(m_nobjslc, weight);
	if (0 < nerr) { nerr += 700;  goto _Exit; }
	if ((m_imagedet&_JMS_DETECT_IMAGE) || (m_imagedet&_JMS_DETECT_WAVEREAL)) { // non-default real-space readout
		nerr = jco->IFT(); // inverse FFT
		if (0 < nerr) {	nerr += 800;  goto _Exit; }
		nerr = ReadoutImgDet_d(m_nobjslc, weight);
		if (0 < nerr) { nerr += 900;  goto _Exit; }
	}
	//
_Exit:
	return nerr;
}