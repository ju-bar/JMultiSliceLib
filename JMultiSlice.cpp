//
// C++ source file: JMultiSlice.cpp
// implementation for library JMultislice.lib (declarations see JMultislice.h)
//
//
// Copyright (C) 2018 - Juri Barthel (juribarthel@gmail.com)
// Copyright (C) 2018 - RWTH Aachen University, 52074 Aachen, Germany
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

#include "stdafx.h"
#include "JMultiSlice.h"
#include "NatureConstants.h"
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include <fstream>
#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

using namespace std;


CJMultiSlice::CJMultiSlice()
{
	// initialize the class members
	m_dbg = 0;
	m_rngseed_ex = 0;
	m_rnd_last = 0;
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
	m_imagedet = 0;
	m_ndetslc = 0;
	m_status_setup_CPU = 0;
	m_status_setup_GPU = 0;
	m_threads_CPU_out = 0;
	m_nfftwthreads = 0;
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
	m_d_wav = NULL;
	m_d_wav0 = NULL;
	m_d_pgr = NULL;
	m_d_pro = NULL;
	m_d_det = NULL;
	m_d_detmask = NULL;
	m_d_det_int = NULL;
	m_d_det_img = NULL;
	m_d_det_dif = NULL;
	m_d_det_tmp = NULL;
	m_jcpuco = NULL;
	
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

int CJMultiSlice::SetRNGSeedEx(int seed)
{
	int oldseed = m_rngseed_ex;
	m_rngseed_ex = seed;
	return oldseed;
}

int CJMultiSlice::GetLastRand(void)
{
	return m_rnd_last;
}


void CJMultiSlice::SeedRngEx(int nseed)
{
	m_rngseed_ex = nseed;
	srand((unsigned)m_rngseed_ex);
}

int CJMultiSlice::GetRand(void)
{
	m_rnd_last = rand();
	return m_rnd_last;
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
	int nitems = nx*ny;
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
	if (whichcode&_JMS_CODE_CPU) {
		if (0<AllocMem_h((void**)&m_h_pgr, sizeof(fcmplx*)*m_nscslc, "PhaseGratingSetup", "phase grating addresses", true)) { nerr = 4; goto _Exit; }
	}

	// GPU code - phase grating setup -> Memory prepared but no data set.
	if (whichcode&_JMS_CODE_GPU) {
		if (memoffset > 0) { // make sure to allocate only if data is expected.
			if (0<AllocMem_d((void**)&m_d_pgr, sizeof(cuComplex)*memoffset, "PhaseGratingSetup", "phase gratings", true)) { nerr = 104; goto _Exit; }
		}
	}

_Exit:
	if (nerr == 0) {
		if (whichcode&_JMS_CODE_CPU) m_status_setup_CPU |= _JMS_STATUS_PGR;
		if (whichcode&_JMS_CODE_GPU) m_status_setup_GPU |= _JMS_STATUS_PGR;
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
		m_status_setup_CPU |= _JMS_STATUS_OBJ;
		m_status_setup_GPU |= _JMS_STATUS_OBJ;
	}
	return nerr;
}


int CJMultiSlice::PropagatorSetup(int whichcode, int npro, int *proidx)
{
	// No cross-check is done whether the structure slice IDs in proidx are valid
	// the length of array proidx is assumed to be m_nscslc
	int nerr = 0;
	int64_t nproitems = 0;
	int nitems = m_nscx*m_nscy;
	m_npro = 0;
	if (npro < 1 || NULL == proidx) goto _Exit;

	// Setup relevant for both codes
	m_npro = npro; // number of propagator functions
	if (0 < AllocMem_h((void**)&m_proidx, sizeof(int)*m_nscslc, "PropagatorSetup", "propagator access indices")) { nerr = 1; goto _Exit; }
	memcpy(m_proidx, proidx, sizeof(int)*m_nscslc); // copy propagator access indices

	// CPU code - propagator setup -> Memory prepared but no data set.
	if (whichcode&_JMS_CODE_CPU) {
		if (0 < AllocMem_h((void**)&m_h_pro, sizeof(fcmplx*)*m_npro, "PropagatorSetup", "propagator address list", true)) { nerr = 2; goto _Exit; }
	}

	// GPU code - propagator setup -> Memory prepared but no data set.
	if (whichcode&_JMS_CODE_GPU) {
		nproitems = (int64_t)m_npro*nitems;
		if (nproitems > 0) { // make sure to allocate only if data is expected.
			if (0 < AllocMem_d((void**)&m_d_pro, sizeof(cuComplex*)*nproitems, "PropagatorSetup", "propagators")) { nerr = 102; goto _Exit; }
		}
	}

_Exit:
	if (nerr == 0) {
		if (whichcode&_JMS_CODE_CPU) m_status_setup_CPU |= _JMS_STATUS_PRO;
		if (whichcode&_JMS_CODE_GPU) m_status_setup_GPU |= _JMS_STATUS_PRO;
	}
	return nerr;
}


int CJMultiSlice::DetectorSetup(int whichcode, int ndetper, int ndetint, int imagedet, int nthreads_CPU)
{
	// Assumes prior object slice setup done and successful via ObjectSliceSetup
	int nerr = 0;
	int islc = 0;
	int64_t ndetitems = 0;
	int nitems = m_nscx*m_nscy;
	int nthreads = max(1, nthreads_CPU);
	size_t nbytes = 0;
	// clean previous setup (should usually not be necessary)
	m_ndet = 0;
	m_ndetper = 0;
	m_imagedet = 0;
	m_ndetslc = 0;
	m_threads_CPU_out = 0;

	// check for prerequisites
	if (((whichcode&_JMS_CODE_CPU) > 0) && ((m_status_setup_CPU&_JMS_STATUS_OBJ) == 0)) {
		cerr << "Error: (DetectorSetup): no previous object slice setup for CPU code." << endl;
		nerr = 1;
	}
	if (((whichcode&_JMS_CODE_GPU) > 0) && ((m_status_setup_GPU&_JMS_STATUS_OBJ) == 0)) {
		cerr << "Error: (DetectorSetup): no previous object slice setup for GPU code." << endl;
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
	if (ndetper > 0) { // handle intermediate detections
		m_ndetper = ndetper;
		for (islc = 0; islc < m_nobjslc; islc++) {
			if (0 == islc%ndetper) { // intermediate slice detection
				m_det_objslc[islc] = m_ndetslc;
				m_ndetslc++;
			}
		}
	}
	m_det_objslc[m_nobjslc] = m_ndetslc; // hash to exit-plane detector arrays
	m_ndetslc++; // count with exit plane

	if (ndetint > 0) { // use integrating diffraction plane detectors
		m_ndet = ndetint; // set number of detectors
		if (0 < AllocMem_h((void**)&m_detmask_len, sizeof(int)*m_ndet, "DetectorSetup", "detector mask lengths", true)) { nerr = 7; goto _Exit; }
	}



	// CPU code - detector setup -> Memory prepared but no data set.
	if (whichcode&_JMS_CODE_CPU) {
		m_threads_CPU_out = nthreads;
		if (m_ndet > 0) { // Prepare integrating detector arrays
			nbytes = sizeof(float*)*(size_t)m_ndet;
			if (0 < AllocMem_h((void**)&m_h_det, nbytes, "DetectorSetup", "detector function addresses", true))  { nerr = 3; goto _Exit; }
			nbytes = sizeof(int*)*(size_t)m_ndet;
			if (0 < AllocMem_h((void**)&m_h_detmask, nbytes, "DetectorSetup", "detector mask addresses", true)) { nerr = 8; goto _Exit; }
			nbytes = sizeof(float)*(size_t)m_ndet*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_int, nbytes, "DetectorSetup", "integrated detector channels", true)) { nerr = 4; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_IMAGE) { // prepare image plane readout arrays
			nbytes = sizeof(float)*(size_t)nitems*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_img, nbytes, "DetectorSetup", "image planes", true)) { nerr = 5; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_DIFFRACTION) { // prepare diffraction plane readout arrays
			nbytes = sizeof(float)*(size_t)nitems*m_ndetslc*m_threads_CPU_out;
			if (0 < AllocMem_h((void**)&m_h_det_dif, nbytes, "DetectorSetup", "diffraction planes", true)) { nerr = 6; goto _Exit; }
		}
	}
	
	// GPU code - detector setup -> Memory prepared but no data set.
	if (whichcode&_JMS_CODE_GPU) {
		ndetitems = (int64_t)m_ndet*nitems;
		if (ndetitems > 0) { // make sure to allocate only if data is expected.
			nbytes = sizeof(float)*(size_t)ndetitems;
			if (0 < AllocMem_d((void**)&m_d_det, nbytes, "DetectorSetup", "detector functions", true)) { nerr = 103; goto _Exit; }
			nbytes = sizeof(int)*(size_t)ndetitems;
			if (0 < AllocMem_d((void**)&m_d_detmask, nbytes, "DetectorSetup", "detector masks", true)) { nerr = 108; goto _Exit; }
			nbytes = sizeof(float)*(size_t)m_ndet*m_ndetslc;
			if (0 < AllocMem_h((void**)&m_d_det_int, nbytes, "DetectorSetup", "integrated detector channels, GPU", true)) { nerr = 104; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_IMAGE) { // prepare image plane readout arrays
			nbytes = sizeof(float)*(size_t)nitems*m_ndetslc; // calculate number of bytes to allocate
			if (0 < AllocMem_d((void**)&m_d_det_img, nbytes, "DetectorSetup", "image planes", true)) { nerr = 105; goto _Exit; }
		}
		if (m_imagedet&_JMS_DETECT_DIFFRACTION) { // prepare diffraction plane readout arrays
			nbytes = sizeof(float)*(size_t)nitems*m_ndetslc; // calculate number of bytes to allocate
			if (0 < AllocMem_d((void**)&m_d_det_dif, nbytes, "DetectorSetup", "diffraction planes", true)) { nerr = 106; goto _Exit; }
		}
	}
	
_Exit:
	if (nerr == 0) {
		if (whichcode&_JMS_CODE_CPU) m_status_setup_CPU |= _JMS_STATUS_DET;
		if (whichcode&_JMS_CODE_GPU) m_status_setup_GPU |= _JMS_STATUS_DET;
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

int CJMultiSlice::CalculateProbeWaveFourier(CJProbeParams prm, fcmplx *wav)
{
	// Assumes that general parameters are set up.
	if (m_wl <= 0.f || m_a0[0] * m_a0[1] <= 0.f || m_nscx*m_nscy <= 0 || NULL == wav || prm.m_alpha <= 0.f) {
		return 1; // error, invalid setup
	}
	// uses m_a0, m_nscx, m_nscy, m_wl from CJMultiSlice
	CJProbeParams lprm;
	lprm = prm;
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
	double gmax = min(itogx*nyqx, itogy*nyqy); // smaller of the x,y gmax
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

int CJMultiSlice::LoadSTEMDetectorProfile(std::string sfile, int &len, float &refpix, float** profile)
{
	// init
	int nResult = 0;
	float ftmp = 0.f;
	string sline;
	string ssub;
	string ssep = "' ";
	len = 0;
	refpix = 0.f;
	*profile = NULL;
	int idx = 0, nline = 0;
	int i0 = 0, i1 = 0;
	// try opening the file for reading list directed data
	ifstream instr( sfile );
	if (instr.is_open()) {
		while (getline(instr, sline)) {
			nline++;
			if (nline == 1) { // header line
				// parse two numbers separated by comma and white spaces
				i0 = (int)sline.find_first_not_of(ssep, 0); // get first character index which is no separator
				i1 = (int)sline.find_first_of(ssep, i0); // get the next character index which is a separator
				// the part of the string from i0 to i1 is the first (int) number
				if (i0 != string::npos && i1 != string::npos && i0 < i1) {
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
					if (i1 == string::npos) { // in case no more separators are found ...
						i1 = (int)sline.length(); // set the length of the string as terminating index
					}
					// the part of the string from i0 to i1 is the second (float) number
					if (i0 != string::npos && i0 < i1) {
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
	// The angular calibration is done via refpix and beta1 which identify the pixel and diffraction angle
	// providing scale for a linear dependency. Linear interpolation is used.
	float sens = 1.0f;
	float tpix = 0.f; // target pixel of the profile
	float tfrac = 0.f;
	int ipix = 0;
	if (NULL != profile && len > 0 && refpix >= 0.f) {
		tpix = min( max(0.f, theta * refpix / beta1), (float)(len-1) ); // get target pixel clipped(!) to range 0 ... len-1
		ipix = (int)tpix; // get integer part for linear interpolation
		tfrac = tpix - ipix; // get fractional part for linear interpolation
		if (ipix == (len - 1)) { // special case at the upper bound
			sens = profile[ipix];
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
	BOOL bmakemask = FALSE; // preset the creation of a mask to 'no' by default
	BOOL busesprof = FALSE; // preset the sensitivity profile usage to 'no' by default
	msklen = 0; // reset mask length to 0
	if (NULL != msk) {
		bmakemask = TRUE;
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
			cout << "(CJMultiSlice::CalculateRingDetector) Warning: Failed to read detector sensitivity profile." << endl;
			cout << "  file name: " << sdsprofile << endl;
			cout << "  proceeding with default sensitivity of 1." << endl;
			nDetProfileLen = 0;
		}
		else {
			busesprof = TRUE; // activate to use the sensitivity profile
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
		nlvar = min(m_nvarslc[islc],nvar); // update number of variants
		if (nlvar < m_nvarslc[islc]) {
			m_nvarslc[islc] = nlvar;
		}
		nbytes = sizeof(fcmplx)*(size_t)nlvar*nitems;
		if (nbytes == 0) { nerr = 1; goto _Exit; } // something is wrong with the setup
		if (whichcode&_JMS_CODE_CPU) {
			if (NULL == m_h_pgr) { nerr = 2; goto _Exit; } // cannot link, pgr array not ready
			m_h_pgr[islc] = pgr; // just link the pointer to host address, the rest is handled by lib parameters
		}
		if (whichcode&_JMS_CODE_GPU) {
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
		if (whichcode&_JMS_CODE_CPU) {
			if (NULL == m_h_pro) { nerr = 2; goto _Exit; } // cannot link, internal pro array not ready
			m_h_pro[ipro] = pro; // just link the pointer to host address, the rest is handled by lib parameters
		}
		if (whichcode&_JMS_CODE_GPU) {
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
		if (whichcode&_JMS_CODE_CPU) {
			if (NULL == m_h_det) { nerr = 2; goto _Exit; } // cannot link, internal det array not ready
			m_h_det[idet] = det; // just link the pointer to host address, the rest is handled by lib parameters
			if (NULL != m_h_detmask) { // detector mask can be linked
				m_h_detmask[idet] = NULL; // preset to none
				if (NULL != msk && 0 < msklen) { // use input mask
					m_h_detmask[idet] = msk; // set
				}
			}
		}
		if (whichcode&_JMS_CODE_GPU) {
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
		qlenmax = max(qlenmax, m_nvarslc[i]);
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
		k = (int)((float)m_nvarslc[islc] * GetRand() / (float)(RAND_MAX + 1));
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
			ivar = (int)((float)m_nvarslc[islc] * GetRand() / (float)(RAND_MAX + 1));
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
	
	if (m_rngseed_ex==0) {
		srand((unsigned)time(NULL));
	}
	else {
		srand((unsigned)m_rngseed_ex);
	}

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
	
	if (whichcode&_JMS_CODE_CPU && (m_status_setup_CPU & _JMS_THRESHOLD_CORE) > 0 && nitems > 0 ) {
		if (NULL != m_jcpuco && m_nfftwthreads>0 ) { // there seem to be old cores, get rid of them
			for (icore = 0; icore < m_nfftwthreads; icore++) {
				m_jcpuco[icore].Deinit();
			}
			delete[] m_jcpuco;
			m_jcpuco = NULL;
			m_nfftwthreads = 0;
		}
		if (nCPUthreads > m_threads_CPU_out) { nerr = 1; goto _Exit; }
		m_nfftwthreads = max(1, nCPUthreads); // new number of cores
		m_jcpuco = new CJFFTWcore[m_nfftwthreads]; // allocate new FFTW cores
		if (NULL == m_jcpuco) {
			cerr << "Error: (InitCore): Failed to create " << m_nfftwthreads << " new CPU cores" << endl;
			m_nfftwthreads = 0; 
			nerr = 2; goto _Exit; 
		} // failed to allocate FFTW objects
		if (0 < AllocMem_h((void**)&m_status_calc_CPU, sizeof(int)*m_nfftwthreads, "InitCore", "CPU calculation status", true)) { nerr = 10; goto _Exit; }
		for (icore = 0; icore < m_nfftwthreads; icore++) { // initialize all FFTW cores
			if (0 == m_jcpuco[icore].Init(2, pdims, FFTW_MEASURE)) {
				ncore_ready++;
			}
			else {
				cerr << "Error: (InitCore): Failed to initialize CPU core #" << icore+1 << endl;
			}
		}
		if (ncore_ready < m_nfftwthreads) { nerr = 3; goto _Exit; } // not all cores are initialized
		nbytes = sizeof(fcmplx)*(size_t)nitems;
		if (0 < AllocMem_h((void**)&m_h_wav0, nbytes, "InitCore", "wave function backup", true)) { nerr = 4; goto _Exit; }
		nbytes = sizeof(fcmplx)*(size_t)nitems*m_nfftwthreads;
		if (0 < AllocMem_h((void**)&m_h_wav, nbytes, "InitCore", "wave functions", true)) { nerr = 5; goto _Exit; }
		m_status_setup_CPU |= _JMS_STATUS_CORE; // mark CPU core setup as completed
	}

	if (whichcode&_JMS_CODE_GPU && (m_status_setup_GPU & _JMS_THRESHOLD_CORE) > 0 && nitems > 0) {
		m_jgpuco.Deinit();
		if (0 < m_jgpuco.Init(2, pdims)) {
			cerr << "Error: (InitCore): Failed to initialize GPU core." << endl;
			nerr = 103; goto _Exit; 
		} // gpu core init failure?
		nbytes = sizeof(cuComplex)*(size_t)nitems;
		if (0 < AllocMem_d((void**)&m_d_wav0, nbytes, "InitCore", "wave function backup", true)) { nerr = 104; goto _Exit; }
		if (0 < AllocMem_d((void**)&m_d_wav, nbytes, "InitCore", "wave function", true)) { nerr = 105; goto _Exit; }
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
		m_status_setup_GPU |= _JMS_STATUS_CORE; // mark GPU core setup as completed
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
	if (whichcode&_JMS_CODE_CPU && (m_status_setup_CPU & _JMS_THRESHOLD_CALC) > 0 &&
		nitems > 0 && m_h_wav0 != NULL && wavuse != NULL) {
		if (NULL == memcpy(m_h_wav0, wavuse, nbytes)) {
			nerr = 2; goto _Exit;
		}
		// debug
		// ftmp2 = GetAbsTotal(m_h_wav0, nitems);
		// end debug
	}
	if (whichcode&_JMS_CODE_GPU && (m_status_setup_GPU & _JMS_THRESHOLD_CALC) > 0 &&
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


int CJMultiSlice::GetUnscrambleHash(UINT* phash)
{
	int nerr = 0;
	int nx2 = 0, ny2 = 0;
	int i = 0, j = 0, idy = 0, i1 = 0, j1 = 0, idy1 = 0;
	UINT idx = 0, idx1 = 0;
	if (NULL == phash) { // I/O buffer needs to be allocated
		nerr = 1; goto _exit_point;
	}
	if ((0 < (m_status_setup_CPU & _JMS_THRESHOLD_CALC)) || (0 < (m_status_setup_GPU & _JMS_THRESHOLD_CALC))) { // core needs to be initialized
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
	// Requires completed setup: _JMS_THRESHOLD_CALC
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

	if (whichcode&_JMS_CODE_CPU && (m_status_setup_CPU & _JMS_THRESHOLD_CALC) > 0 &&
		iThread >= 0 && iThread < m_nfftwthreads) {
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

	if (whichcode&_JMS_CODE_GPU && (m_status_setup_GPU & _JMS_THRESHOLD_CALC) > 0 ) {
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



int CJMultiSlice::GetResult(int whichcode, int whichresult, float *dst, int iThread)
{
	int nerr = 0;
	int nbytes = 0;
	int nitems = 0;
	float *detcur = 0;
	cudaError cuerr;

	if (NULL == dst) { // invalid destination, do not even try
		goto _Exit;
	}
	// TODO: check calculation status !

	if ((whichcode&_JMS_CODE_CPU) > 0 && iThread >=0 && iThread < m_nfftwthreads && m_ndetslc>0 ) {
		if (whichresult == _JMS_DETECT_INTEGRATED) {
			nitems = m_ndet*m_ndetslc;
			nbytes = sizeof(float)*(size_t)nitems;
			detcur = m_h_det_int + iThread * nitems;
			memcpy(dst, detcur, nbytes);
		}
		if (whichresult == _JMS_DETECT_IMAGE) {
			nitems = m_nscx*m_nscy*m_ndetslc;
			nbytes = sizeof(float)*(size_t)nitems;
			detcur = m_h_det_img + iThread * nitems;
			memcpy(dst, detcur, nbytes);
		}
		if (whichresult == _JMS_DETECT_DIFFRACTION) {
			nitems = m_nscx*m_nscy*m_ndetslc;
			nbytes = sizeof(float)*(size_t)nitems;
			detcur = m_h_det_dif + iThread * nitems;
			memcpy(dst, detcur, nbytes);
		}
	}

	if ((whichcode&_JMS_CODE_GPU) > 0 && m_ndetslc>0) {
		if (whichresult == _JMS_DETECT_INTEGRATED) {
			nitems = m_ndet*m_ndetslc;
			nbytes = sizeof(float)*(size_t)nitems;
			detcur = m_d_det_int;
			memcpy(dst, detcur, nbytes);
		}
		if (whichresult == _JMS_DETECT_IMAGE) {
			nitems = m_nscx*m_nscy*m_ndetslc;
			nbytes = sizeof(float)*(size_t)nitems;
			detcur = m_d_det_img;
			cuerr = cudaMemcpy(dst, detcur, nbytes, cudaMemcpyDeviceToHost);
		}
		if (whichresult == _JMS_DETECT_DIFFRACTION) {
			nitems = m_nscx*m_nscy*m_ndetslc;
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

	// _JMS_STATUS_CORE
	DeallocMem_h((void**)&m_h_wav);
	DeallocMem_h((void**)&m_h_wav0);
	DeallocMem_d((void**)&m_d_wav);
	DeallocMem_d((void**)&m_d_wav0);
	DeallocMem_d((void**)&m_d_det_tmp);
	DeallocMem_h((void**)&m_status_calc_CPU);
	if (NULL != m_jcpuco && m_nfftwthreads>0) { // there seem to be old cores, get rid of them
		for (icore = 0; icore < m_nfftwthreads; icore++) {
			m_jcpuco[icore].Deinit();
		}
		delete[] m_jcpuco;
		m_jcpuco = NULL;
		m_nfftwthreads = 0;
	}
	m_jgpuco.Deinit();
	DeallocMem_h((void**)&m_h_fnx);
	DeallocMem_h((void**)&m_h_fny);
	DeallocMem_d((void**)&m_d_knx);
	DeallocMem_d((void**)&m_d_kny);
	if ((m_status_setup_CPU & _JMS_STATUS_CORE) > 0) m_status_setup_CPU -= _JMS_STATUS_CORE;
	if ((m_status_setup_GPU & _JMS_STATUS_CORE) > 0) m_status_setup_GPU -= _JMS_STATUS_CORE;
	
	// _JMS_STATUS_DET
	DeallocMem_h((void**)&m_h_det_dif);
	DeallocMem_h((void**)&m_h_det_img);
	DeallocMem_h((void**)&m_h_det_int);
	DeallocMem_h((void**)&m_h_det);
	DeallocMem_h((void**)&m_h_detmask);
	DeallocMem_d((void**)&m_d_det_dif);
	DeallocMem_d((void**)&m_d_det_img);
	DeallocMem_h((void**)&m_d_det_int); // ! This array is really on host memory
	DeallocMem_d((void**)&m_d_det);
	DeallocMem_d((void**)&m_d_detmask);
	DeallocMem_h((void**)&m_det_objslc);
	DeallocMem_h((void**)&m_detmask_len);
	m_ndet = 0;
	m_ndetper = 0;
	m_imagedet = 0;
	m_ndetslc = 0;
	m_threads_CPU_out = 0;
	if ((m_status_setup_CPU & _JMS_STATUS_DET) > 0) m_status_setup_CPU -= _JMS_STATUS_DET;
	if ((m_status_setup_GPU & _JMS_STATUS_DET) > 0) m_status_setup_GPU -= _JMS_STATUS_DET;

	// _JMS_STATUS_PRO
	DeallocMem_h((void**)&m_h_pro);
	DeallocMem_d((void**)&m_d_pro);
	DeallocMem_h((void**)&m_proidx);
	m_npro = 0;
	if ((m_status_setup_CPU & _JMS_STATUS_PRO) > 0) m_status_setup_CPU -= _JMS_STATUS_PRO;
	if ((m_status_setup_GPU & _JMS_STATUS_PRO) > 0) m_status_setup_GPU -= _JMS_STATUS_PRO;

	// _JMS_STATUS_OBJ
	DeallocMem_h((void**)&m_objslc);
	m_nobjslc = 0;
	if ((m_status_setup_CPU & _JMS_STATUS_OBJ) > 0) m_status_setup_CPU -= _JMS_STATUS_OBJ;
	if ((m_status_setup_GPU & _JMS_STATUS_OBJ) > 0) m_status_setup_GPU -= _JMS_STATUS_OBJ;

	// _JMS_STATUS_PGR
	DeallocMem_h((void**)&m_nvarslc);
	DeallocMem_h((void**)&m_slcthick);
	DeallocMem_h((void**)&m_h_pgr);
	DeallocMem_h((void**)&m_slcoffset);
	DeallocMem_d((void**)&m_d_pgr);
	m_npgx = 0;
	m_npgy = 0;
	m_nscslc = 0;
	m_nslcvar = 0;
	if ((m_status_setup_CPU & _JMS_STATUS_PGR) > 0) m_status_setup_CPU -= _JMS_STATUS_PGR;
	if ((m_status_setup_GPU & _JMS_STATUS_PGR) > 0) m_status_setup_GPU -= _JMS_STATUS_PGR;

	return nerr;
}


void CJMultiSlice::CleanFFTW(void)
{
	CJFFTWcore ctmp;
	ctmp.CleanupFFTW();
}

////////////////////////////////////////////////////////////////////////////////
//
// CPU CODE
//
////////////////////////////////////////////////////////////////////////////////

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


int CJMultiSlice::AllocMem_h(void ** _h_a, size_t size, char* callfn, char* arrnam, bool zero)
{
	if (NULL != *_h_a) {
		free(*_h_a);
		*_h_a = NULL;
	}
	if ( (NULL == *_h_a) && (size > 0) ) {
		*_h_a = malloc(size);
		if (NULL == *_h_a) {
			cerr << "Error (" << callfn << "): Failed to allocate host memory (" << arrnam << ")" << endl;
			cerr << "- requested size [MB]: " << size / 1048576 << endl;
			return 1;
		}
		if (zero) {
			if (NULL==memset(*_h_a, 0, size)) {
				cerr << "Error (" << callfn << "): Failed to zeroe host memory (" << arrnam << ")" << endl;
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
	if (m_ndet > 0 && m_ndetslc > 0 && m_h_det_int != NULL) {
		nitems = m_ndet*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		memset(&m_h_det_int[iThread*nitems], 0, nbytes);
	}
	// probe image readouts
	if ((m_imagedet & _JMS_DETECT_IMAGE) > 0 && m_ndetslc > 0 && m_h_det_img != NULL) {
		nitems = m_nscx*m_nscy*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		memset(&m_h_det_img[iThread*nitems], 0, nbytes);
	}
	// probe diffraction readouts
	if ((m_imagedet & _JMS_DETECT_DIFFRACTION) > 0 && m_ndetslc > 0 && m_h_det_dif != NULL) {
		nitems = m_nscx*m_nscy*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		memset(&m_h_det_dif[iThread*nitems], 0, nbytes);
	}
//_Exit:
	// error handling
	return nerr;
}


float CJMultiSlice::GetCoreAbsTotal_h(int iThread)
{
	float ftot = 0.f;
	if (iThread >= 0 && iThread < m_nfftwthreads && (0 < (m_status_setup_CPU & _JMS_THRESHOLD_CALC))) {
		ftot = m_jcpuco[iThread].GetDataTotalPow();
	}
	return ftot;
}

float CJMultiSlice::GetCoreAbsTotal_d(void)
{
	float ftot = 0.f;
	if (0 < (m_status_setup_GPU & _JMS_THRESHOLD_CALC)) {
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
	/*       Performs a Kahan sum on in_1[i]*in_2[i]                              */
	/*       https://en.wikipedia.org/wiki/Kahan_summation_algorithm              */
	double dsum = 0.0;
	double dc = 0.0;
	double dy = 0.0;
	double dt = 0.0;
	if (len > 0) {
		for (size_t i = 0; i < len; i++) {
			dy = ((double)in_1[i])*(double)in_2[i] - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
	}
	return (float)dsum;
}

float CJMultiSlice::MaskedDotProduct_h(int *mask, float *in_1, float *in_2, size_t lenmask)
{
	/*       Performs a Kahan sum on in_1[i]*in_2[i]                              */
	/*       https://en.wikipedia.org/wiki/Kahan_summation_algorithm              */
	double dsum = 0.0;
	double dc = 0.0;
	double dy = 0.0;
	double dt = 0.0;
	unsigned int imask;
	if (lenmask > 0) {
		for (size_t i = 0; i < lenmask; i++) {
			imask = (unsigned)mask[i];
			dy = ((double)in_1[imask])*(double)in_2[imask] - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
	}
	return (float)dsum;
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
		jvar = (int)((float)m_nvarslc[jslc] * (float)GetRand() / (RAND_MAX + 1)); // random number 0 ... m_nvarslc[jslc]-1
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
	CJFFTWcore *jco = NULL;
	float *dif = NULL;
	float *det = NULL;
	float *out = NULL;
	int *msk = NULL;
	int msklen = 0;
	int idetslc = -1;
	int nitems = m_nscx*m_nscy;
	int idx = 0;
	int idet = 0;
	if (m_ndet == 0 && (m_imagedet&_JMS_DETECT_DIFFRACTION) == 0) { goto _Exit; } // handle no diffraction detection requested
	if (m_ndetslc == 0) { goto _Exit; } // handle no detection
	if (NULL != m_jcpuco && NULL != m_det_objslc && (iSlice >= 0) && (iSlice < m_nobjslc+1) && (nitems > 0)) {
		idetslc = m_det_objslc[iSlice]; // get index of the slice in detection output arrays
		if (idetslc >= 0) { // there is detection registered for this slice
			if (0 < AllocMem_h((void**)&dif, sizeof(float)*nitems, "ReadoutDifDet_h", "diffraction pattern")) { nerr = 1; goto _Exit; }
			jco = &m_jcpuco[iThread]; // link cpu core
			jco->GetDataPow(dif); // copy the power of the current data in cpu core
			if (m_ndet > 0) { // retrieve data for all integrating detectors
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
					}
					else {
						*out += (weight*DotProduct_h(dif, det, (size_t)nitems));
					}
				}
			}
			if (m_imagedet&_JMS_DETECT_DIFFRACTION) { // retrieve data for diffraction pattern
				out = m_h_det_dif + (iThread*m_ndetslc + idetslc)*nitems; // pointer to integrating output channel
				for (idx = 0; idx < nitems; idx++) {
					out[idx] += (weight*dif[idx]);
				}
			}
			DeallocMem_h((void**)&dif);
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
	CJFFTWcore *jco = NULL;
	float *img = NULL;
	float *out = NULL;
	int idetslc = -1;
	int nitems = m_nscx*m_nscy;
	int idx = 0;
	if ((m_imagedet&_JMS_DETECT_IMAGE) == 0) { goto _Exit; } // handle no image detection requested
	if (m_ndetslc == 0) { goto _Exit; } // handle no detection
	if (NULL != m_jcpuco && NULL != m_det_objslc && (iSlice >= 0) && (iSlice < m_nobjslc+1) && (nitems > 0)) {
		idetslc = m_det_objslc[iSlice]; // get index of the slice in detection output arrays
		if (idetslc >= 0) { // there is detection registered for this slice
			if (0 < AllocMem_h((void**)&img, sizeof(float)*nitems, "ReadoutImgDet_h", "image pattern")) { nerr = 1; goto _Exit; }
			jco = &m_jcpuco[iThread]; // link cpu core
			jco->GetDataPow(img); // copy the power of the current data in cpu core
			if (m_imagedet&_JMS_DETECT_IMAGE) { // retrieve data for image pattern // though already checked, we leave this 'if' for future extensions 
				out = m_h_det_img + (iThread*m_ndetslc + idetslc)*nitems; // pointer to integrating output channel
				for (idx = 0; idx < nitems; idx++) {
					out[idx] += (weight*img[idx]);
				}
			}
			DeallocMem_h((void**)&img);
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
	int i = 0, j = 0, i1 = 0, j1 = 0, idx = 0, idx1 = 0, idy = 0, idy1 = 0; // iterator
	size_t nbytes = sizeof(fcmplx)*(size_t)nitems;
	fcmplx* wavuse = wav;
	fcmplx* _h_wav = NULL;
	fcmplx* wavtmp = NULL; // temp buffer only used when bTranspose == true
	if ( (m_status_setup_CPU & _JMS_THRESHOLD_CALC) == 0 || nitems <= 0 || 
		m_h_wav == NULL || iThread < 0 || iThread >= m_nfftwthreads) {
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
	CJFFTWcore *jco = NULL;
	fcmplx *wav = NULL;
	fcmplx *pgr = NULL;
	fcmplx *pro = NULL;
	int *var = NULL;
	int islc = islc0;
	int jslc = 0;
	int nitems = m_nscx*m_nscy;
	int npgitems = m_npgx*m_npgy;
	bool bsubframe = false;
	if (nitems != npgitems) {
		bsubframe = true;
	}
	if (iThread<0 || iThread>=m_nfftwthreads || NULL== m_jcpuco || NULL==m_h_wav || nitems <= 0 ) {
		//cout << "iThread :" << iThread << ", core: " << m_jcpuco << ", wav: " << m_h_wav << ", << n: " << nitems << endl;
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
		// Default case, the current wave function is not the exit-plane wave function
		// Do the multislice
		for (jslc = islc; jslc < m_nobjslc; jslc++) {
			// 1) readout (Fourier space)
			//if (1 == m_objslc_det[jslc]) {
			if (m_det_objslc[jslc] >= 0) {
				nerr = ReadoutDifDet_h(jslc, iThread, weight);
				if (nerr > 0) { nerr += 100; goto _CancelMS; }
			}
			// 2) scattering (real space)
			nerr = jco->IFT(); // inverse FFT
			if (nerr > 0) { nerr += 200; goto _CancelMS; }
			//if (1 == m_objslc_det[jslc] && (m_imagedet&_JMS_DETECT_IMAGE)) {  // non-default real-space readout
			if ((m_det_objslc[jslc] >= 0) && (m_imagedet&_JMS_DETECT_IMAGE)) {  // non-default real-space readout
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
	if (m_imagedet&_JMS_DETECT_IMAGE) { // non-default real-space readout
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
		cerr << "Error: " << smsg << ", code: " << code << endl;
		cerr << "     - " << cudaGetErrorString(code) << endl;
	}
}

void CJMultiSlice::PostCUDAMemory(size_t nrequestedbytes)
{
	int64_t memavail = 0, memtotal = 0;
	int idev = 0, cc1 = 0, cc2 = 0, nmaxthread = 0;
	idev = GetCurrentGPU();
	if (idev >= 0) {
		if (0 == GetGPUStats(idev, cc1, cc2, nmaxthread, memtotal, memavail)) {
			cerr << "  - requested device memory [MB]: " << nrequestedbytes / 1048576;
			cerr << "  - available device memory [MB]: " << memavail / 1048576;
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
	int nitems = 0;
	size_t nbytes = 0;
	cudaError cuerr;
	// integrating detector readouts (this is host memory, REALLY ! )
	if (m_ndet > 0 && m_ndetslc > 0 && m_d_det_int != NULL) {
		nitems = m_ndet*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		memset(m_d_det_int, 0, nbytes);
	}
	// probe image readouts (this is device memory)
	if ((m_imagedet & _JMS_DETECT_IMAGE) > 0 && m_ndetslc > 0 && m_d_det_img != NULL) {
		nitems = m_nscx*m_nscy*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		cuerr = cudaMemset((void*)m_d_det_img, 0, nbytes);
		if (cuerr != cudaSuccess) {
			PostCUDAError("(ClearDetMem_d): Failed to reset memory m_d_det_img", cuerr);
			return 3;
		}
	}
	// probe diffraction readouts (this is device memory)
	if ((m_imagedet & _JMS_DETECT_DIFFRACTION) > 0 && m_ndetslc > 0 && m_d_det_dif != NULL) {
		nitems = m_nscx*m_nscy*m_ndetslc;
		nbytes = sizeof(float)*nitems;
		cuerr = cudaMemset((void*)m_d_det_dif, 0, nbytes);
		if (cuerr != cudaSuccess) {
			PostCUDAError("(ClearDetMem_d): Failed to reset memory m_d_det_dif", cuerr);
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
		jvar = (int)((float)m_nvarslc[jslc] * (float)GetRand() / (float)(RAND_MAX + 1)); // random number 0 ... m_nvarslc[jslc]-1
	}
	else { // take variant ID from prepared list
		jvar = pVarID[iSlice];
	}
	pgr = m_d_pgr + m_slcoffset[jslc] + (int64_t)jvar*nitems; // pointer to first item of slice variant # jvar
	/*if (m_dbg >= 5) {
		cout << "debug(GetPhaseGrating_d): iSlice " << iSlice << endl;
		cout << "debug(GetPhaseGrating_d): jslc " << jslc << endl;
		cout << "debug(GetPhaseGrating_d): m_nvarslc[jslc] " << m_nvarslc[jslc] << endl;
		if (pVarID) cout << "debug(GetPhaseGrating_d): pVarID[iSlice] " << pVarID[iSlice] << endl;
		cout << "debug(GetPhaseGrating_d): jvar " << jvar << endl;
		cout << "debug(GetPhaseGrating_d): m_slcoffset[jslc] " << m_slcoffset[jslc] << endl;
		cout << "debug(GetPhaseGrating_d): jvar*nitems " << jvar*nitems << endl;
		cout << "debug(GetPhaseGrating_d): m_d_pgr 0x" << pgr << endl;
		cout << "debug(GetPhaseGrating_d): pgr 0x" << pgr << endl;
		cout << "debug(GetPhaseGrating_d): dif 0x" << pgr - m_d_pgr << endl;
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
	int *msk_d = NULL;
	int msklen = 0;
	int idetslc = -1;
	int nitems = m_nscx*m_nscy;
	int idx = 0;
	int idet = 0;
	if (m_ndet == 0 && (m_imagedet&_JMS_DETECT_DIFFRACTION) == 0) { goto _Exit; } // handle no diffraction detection requested
	if (m_ndetslc == 0) { goto _Exit; } // handle no detection
	if (NULL != m_det_objslc && (iSlice >= 0) && (iSlice < m_nobjslc+1) && (nitems > 0)) {
		idetslc = m_det_objslc[iSlice]; // get index of the slice in detection output arrays
		stats.uSize = (unsigned int)nitems;
		stats.nBlockSize = jco->GetBlockSize();
		stats.nGridSize = (nitems + stats.nBlockSize - 1) / stats.nBlockSize;
		if (idetslc >= 0) { // there is detection registered for this slice
			jco->GetDataPow_d(dif_d); // get the power of the current data in gpu core
			if (m_ndet > 0) { // retrieve data for all integrating detectors
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
					}
					else {
						*out_d += (weight*DotProduct_d(dif_d, det_d, nitems, stats.nBlockSize));
					}
				}
			}
			if (m_imagedet&_JMS_DETECT_DIFFRACTION) { // retrieve data for diffraction pattern
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
	int idetslc = -1;
	int nitems = m_nscx*m_nscy;
	int idx = 0;
	if ((m_imagedet&_JMS_DETECT_IMAGE) == 0) { goto _Exit; } // handle no image detection requested
	if (m_ndetslc == 0) { goto _Exit; } // handle no detection
	if (NULL != m_det_objslc && (iSlice >= 0) && (iSlice < m_nobjslc+1) && (nitems > 0)) {
		idetslc = m_det_objslc[iSlice]; // get index of the slice in detection output arrays
		stats.uSize = (unsigned int)nitems;
		stats.nBlockSize = jco->GetBlockSize();
		stats.nGridSize = (nitems + stats.nBlockSize - 1) / stats.nBlockSize;
		if (idetslc >= 0) { // there is detection registered for this slice
			jco->GetDataPow_d(img_d); // get the power of the current data in gpu core
			if (m_imagedet&_JMS_DETECT_IMAGE) { // retrieve data for image pattern
				out_d = m_d_det_img + idetslc*nitems; // pointer to integrating output channel
				// Here is the probe image readout!
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
	int i = 0, j = 0, i1 = 0, j1 = 0, idx = 0, idx1 = 0, idy = 0, idy1 = 0; // iterator
	size_t nbytes = sizeof(fcmplx)*(size_t)nitems;
	fcmplx* wavuse = wav;
	fcmplx* wavtmp = NULL; // temp buffer only used when bTranspose == true
	if ( (m_status_setup_GPU & _JMS_THRESHOLD_CALC) == 0 || nitems <= 0 || m_d_wav == NULL) {
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
	cuerr = cudaMemcpy(m_d_wav0, wavuse, nbytes, cudaMemcpyHostToDevice);
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
	cuComplex *wav = NULL;
	cuComplex *pgr = NULL;
	cuComplex *pro = NULL;
	int *var = NULL;
	int islc = islc0;
	int jslc = 0;
	int nitems = m_nscx*m_nscy;
	int npgitems = m_npgx*m_npgy;
	bool bsubframe = false;
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
											 // Readout will be at slice indices flagged in m_objslc_det only
											 // Call DetectorSetup to generate a coherent setup for the multislice.
											 // The same setup will be used by all threads.
											 // Readout will never include the entrance plane, but may include the exit plane
											 // and will always happen after the propagation to the next slice.
											 //
	if (islc < m_nobjslc) {
		
		// prepare new random variant sequence
		var = (int*)malloc(sizeof(int)*m_nobjslc);
		GetRandomVariantSequence(var);
		
		/*if (m_dbg >= 5) {
			cout << "debug(GPUMultislice): RandomVariantSequence" << endl;
			for (jslc = islc; jslc < m_nobjslc; jslc++) {
				cout << "debug(GPUMultislice): jslc = " << jslc << ", var[jslc] = " << var[jslc] << endl;
			}
		}*/
		// Default case, the current wave function is not the exit-plane wave function
		// Do the multislice
		for (jslc = islc; jslc < m_nobjslc; jslc++) {
			// 1) readout (Fourier space)
			//if (1 == m_objslc_det[jslc]) {
			if (m_det_objslc[jslc] >= 0) {
				nerr = ReadoutDifDet_d(jslc, weight);
				if (0 < nerr) {	nerr += 100; goto _CancelMS; }
			}
			// 2) scattering (real space)
			nerr = jco->IFT(); // inverse FFT
			if (0 < nerr) { nerr += 200; goto _CancelMS; }
			//if (1 == m_objslc_det[jslc] && (m_imagedet&_JMS_DETECT_IMAGE)) {  // non-default real-space readout
			if ((m_det_objslc[jslc] >= 0) && (m_imagedet&_JMS_DETECT_IMAGE)) {  // non-default real-space readout
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
				/*cerr << jslc << endl;
				cerr << pgr << endl;
				cerr << m_d_pgr << endl;*/
				nerr += 400; goto _CancelMS; 
			}
			// 3) propagation (fourier space)
			nerr = jco->FT(); // forward FFT
			if (0 < nerr) { nerr += 500; goto _CancelMS; }
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
	if (m_imagedet&_JMS_DETECT_IMAGE) { // non-default real-space readout
		nerr = jco->IFT(); // inverse FFT
		if (0 < nerr) {	nerr += 800;  goto _Exit; }
		nerr = ReadoutImgDet_d(m_nobjslc, weight);
		if (0 < nerr) { nerr += 900;  goto _Exit; }
	}
	//
_Exit:
	return nerr;
}