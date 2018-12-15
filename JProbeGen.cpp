//
// C++ source file: JProbeGen.cpp
// implementation for library JMultislice.lib (declarations see JProbeGen.h)
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
#include "JProbeGen.h"
#include "NatureConstants.h"
#include <time.h>
#include <math.h>
//#include <stdlib.h>
//#include <fstream>
#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

using namespace std;


// ****************************************************************************
//
// Definition of class CJProbeParams
//

CJProbeParams::CJProbeParams()
{
	// initialize the class members
	m_abrr_num = 0;
	m_abrr_coeff = NULL;
	// aberration table config
	// - max expected index for the order _JMS_ABERRATION_ORDER_MAX
	int nomax = (int)_JPG_ABERRATION_ORDER_MAX;
	m_abrr_num = (int)(4 * nomax - nomax % 2 + nomax * nomax) / 4;
	int nac = 2 * m_abrr_num; // number of supported aberration coefficients
	m_abrr_coeff = (float*)malloc(sizeof(float)*nac);
	//
	m_wl = 0.001969f;
	m_alpha = 0.025f;
	m_alpha_x0 = 0.0f;
	m_alpha_y0 = 0.0f;
	m_btx = 0.0f;
	m_bty = 0.0f;
	m_source_width = 0.0f;
	m_source_shape = 0;
	m_fspread_width = 0.0f;
	m_fspread_kernel_width = 2.0f;
	m_fspread_kernel_samples = 7;
	memset(m_abrr_coeff, 0, sizeof(float)*nac);
}

CJProbeParams::CJProbeParams(const CJProbeParams & src)
{
	// initialize the class members
	m_abrr_num = 0;
	m_abrr_coeff = NULL;
	// aberration table config
	// - max expected index for the order _JMS_ABERRATION_ORDER_MAX
	int nomax = (int)_JPG_ABERRATION_ORDER_MAX;
	m_abrr_num = (int)(4 * nomax - nomax % 2 + nomax * nomax) / 4;
	int nac = 2 * m_abrr_num; // number of supported aberration coefficients
	m_abrr_coeff = (float*)malloc(sizeof(float)*nac);
	//
	m_wl = src.m_wl;
	m_alpha = src.m_alpha;
	m_alpha_x0 = src.m_alpha_x0;
	m_alpha_y0 = src.m_alpha_y0;
	m_btx = src.m_btx;
	m_bty = src.m_bty;
	m_source_width = src.m_source_width;
	m_source_shape = src.m_source_shape;
	m_fspread_width = src.m_fspread_width;
	m_fspread_kernel_width = src.m_fspread_kernel_width;
	m_fspread_kernel_samples = src.m_fspread_kernel_samples;
	int nacsrc = 2*(const_cast <CJProbeParams*>(&src)->GetAberrationNum());
	int naccpy = min(nacsrc, nac);
	memcpy(m_abrr_coeff, src.m_abrr_coeff, sizeof(float)*naccpy);
}

CJProbeParams::~CJProbeParams()
{
	if (NULL != m_abrr_coeff) {
		free(m_abrr_coeff);
	}
}

void CJProbeParams::operator=(const CJProbeParams &other)
{
	int nac = 2 * m_abrr_num; // number of supported aberration coefficients
	int nacsrc = 2 * (const_cast <CJProbeParams*>(&other)->GetAberrationNum()); // number of supported coefficients of the other object
	int naccpy = min(nacsrc, nac); // number of coefficients to copy
	m_wl = other.m_wl;
	m_alpha = other.m_alpha;
	m_alpha_x0 = other.m_alpha_x0;
	m_alpha_y0 = other.m_alpha_y0;
	m_btx = other.m_btx;
	m_bty = other.m_bty;
	m_source_width = other.m_source_width;
	m_source_shape = other.m_source_shape;
	m_fspread_width = other.m_fspread_width;
	m_fspread_kernel_width = other.m_fspread_kernel_width;
	m_fspread_kernel_samples = other.m_fspread_kernel_samples;
	memcpy(m_abrr_coeff, 0, sizeof(float)*naccpy);

}

bool CJProbeParams::operator==(const CJProbeParams &other) const
{
	bool bResult = true;
	int nac = 2 * m_abrr_num; // number of supported aberration coefficients
	int nacsrc = 2 * (const_cast <CJProbeParams*>(&other)->GetAberrationNum()); // number of supported coefficients of the other object
	bResult &= (bool)(nac == nacsrc);
	if (bResult) {
		bResult &= (m_wl == other.m_wl);
		bResult &= (m_alpha == other.m_alpha);
		bResult &= (m_alpha_x0 == other.m_alpha_x0);
		bResult &= (m_alpha_y0 == other.m_alpha_y0);
		bResult &= (m_btx == other.m_btx);
		bResult &= (m_bty == other.m_bty);
		bResult &= (m_source_width == other.m_source_width);
		bResult &= (m_source_shape == other.m_source_shape);
		bResult &= (m_fspread_width == other.m_fspread_width);
		bResult &= (m_fspread_kernel_width == other.m_fspread_kernel_width);
		bResult &= (m_fspread_kernel_samples == other.m_fspread_kernel_samples);
		bResult &= (0 == memcmp(m_abrr_coeff, other.m_abrr_coeff, sizeof(float)*nac));
	}
	return bResult;
}

int CJProbeParams::GetAberrationNum()
{
	return m_abrr_num;
}



// ****************************************************************************
//
// Definition of class CJProbeGen
//

CJProbeGen::CJProbeGen()
{
	m_abrr_num = 0;
	m_abrr_binom = NULL;
	m_abrr_widx = NULL;
	m_amorph_dim = 0;
	m_amorph_samp = 0.f;
	m_amorph_pgr = NULL;
	m_dLastProbeCalculationTime = 0.;
	// aberration table config
	// - max expected index for the order _JMS_ABERRATION_ORDER_MAX
	int nomax = (int)_JPG_ABERRATION_ORDER_MAX;
	m_abrr_num = (int)(4*nomax-nomax%2+nomax*nomax)/4;
	
	int nac = 2 * m_abrr_num; // number of supported aberration coefficients

	// - allocate and set the aberration index table
	//   This allocation exists until the object is destroyed.

	m_abrr_widx = (int*)malloc(nac * sizeof(int));
	memset(m_abrr_widx, 0, nac * sizeof(int));
	if (m_abrr_widx) {
		// the following setup defines the sequence of aberrations in the interface list
		// index of aberration (x,y) -> (order, rotational symmetry)
		int l = 0;
		for (int m = 1; m <= _JPG_ABERRATION_ORDER_MAX; m++) { // loop over aberration order 1 ... _JMS_ABERRATION_ORDER_MAX
			for (int n = 0; n <= m; n++) { // loop over aberration rotational symmetry 0 .. m
				if (0 == (m + n) % 2) { // valid combination of order and rotational symmetry (1,1), (2,0), (2,2), (3,1) ... sum is even
					m_abrr_widx[2 * l] = m; // store order
					m_abrr_widx[1 + 2 * l] = n; // store rotational symmetry
					l++; // next aberration index
				}
			}
		}
	}
	//
	int nol = 1 + 2 * (int)_JPG_ABERRATION_ORDER_MAX; // number of order indices
	int nsz = nol * nol; // total order items
						 // - allocate and set the aberration binomial table
						 //   This allocation exists until the object is destroyed.
	m_abrr_binom = (int*)malloc(nsz * sizeof(int));
	memset(m_abrr_binom, 0, nsz * sizeof(int));
	if (m_abrr_binom) {
		// the following calculates binomial factors (n over k) for n, k = aberration order
		int res = 1;
		int lk = 0;
		int i = 0;
		for (int k = 0; k < nol; k++) {
			for (int n = 0; n < nol; n++) {
				if (n < k) {
					res = 0;
				}
				else {
					res = 1;
					lk = k;
					// Since C(n, k) = C(n, n-k) 
					if (k > n - k) {
						lk = n - k;
					}
					// Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1] 
					for (i = 0; i < lk; ++i) {
						res *= (n - i);
						res /= (i + 1);
					}
				}
				m_abrr_binom[n + k * nol] = res;
			}
		}
	}
	// set the aberration names
	m_abrr_name.push_back("Image shift"); // 1
	m_abrr_name.push_back("Defocus"); // 2
	m_abrr_name.push_back("2-fold astigmatism"); // 3
	m_abrr_name.push_back("Coma"); // 4
	m_abrr_name.push_back("3-fold astigmatism"); // 5
	m_abrr_name.push_back("Spherical aberration (Cs)"); // 6
	m_abrr_name.push_back("Star aberration"); // 7
	m_abrr_name.push_back("4-fold astigmatism"); // 8
	m_abrr_name.push_back("Coma (5th order)"); // 9
	m_abrr_name.push_back("3-lobe aberration"); // 10
	m_abrr_name.push_back("5-fold astigmatism"); // 11
	m_abrr_name.push_back("Spherical aberration (C5)"); // 12
	m_abrr_name.push_back("Star aberration (6th order)"); // 13
	m_abrr_name.push_back("Rosette aberration"); // 14
	m_abrr_name.push_back("6-fold astigmatism"); // 15
	m_abrr_name.push_back("Coma (7th order)"); // 16
	m_abrr_name.push_back("Aberration a73"); // 17
	m_abrr_name.push_back("Aberration a75"); // 18
	m_abrr_name.push_back("7-fold astigmatism"); // 19
	m_abrr_name.push_back("Spherical aberration (C7)"); // 20
	m_abrr_name.push_back("Aberration a82"); // 21
	m_abrr_name.push_back("Aberration a84"); // 22
	m_abrr_name.push_back("Aberration a86"); // 23
	m_abrr_name.push_back("8-fold astigmatism"); // 24

	// set the aberration symbols
	char stmp[8];
	for (int j = 0; j < m_abrr_num; j++) {
		sprintf_s(stmp, 8, "a%i%i", m_abrr_widx[2 * j], m_abrr_widx[2 * j + 1]);
		m_abrr_symbol.push_back(stmp);
	}

	// set the probe function names
	m_pfunc_name.push_back("Coherent Probe Intensity"); // 1
	m_pfunc_name.push_back("Coherent Probe Phase"); // 2
	m_pfunc_name.push_back("Partially Spatial Coherent Probe Intensity"); // 3
	m_pfunc_name.push_back("Partially Temporal Coherent Probe Intensity"); // 4
	m_pfunc_name.push_back("Partially Coherent Probe Intensity"); // 5
	m_pfunc_name.push_back("Real Part of Probe Wave Function"); // 6
	m_pfunc_name.push_back("Imaginary Part of Probe Wave Function"); // 7
	m_pfunc_name.push_back("Aberration Phase Plate"); // 8
	m_pfunc_name.push_back("Ronchigram of thin amorphous sample"); // 9
}

CJProbeGen::~CJProbeGen()
{
	if (NULL != m_abrr_widx) free(m_abrr_widx);
	if (NULL != m_abrr_binom) free(m_abrr_binom);
	if (NULL != m_amorph_pgr) free(m_amorph_pgr);
	m_abrr_name.clear();
	m_abrr_symbol.clear();
	m_pfunc_name.clear();
}


int CJProbeGen::GetMaxAberrationOrder(void)
{
	return (int)_JPG_ABERRATION_ORDER_MAX;
}

int CJProbeGen::GetAberrationNum(void)
{
	return m_abrr_num;
}

int CJProbeGen::GetMaxProbeFunction(void)
{
	return (int)_JPG_PROBEFUNC_NUM;
}

int CJProbeGen::GetAberrationName(int idx, string * abrr_name)
{
	if (NULL == abrr_name) { // invalid parameter reference <abrr_name>
		return 1;
	}
	if (idx > m_abrr_num || idx < 0) { // invalid aberration index <idx>
		return 2;
	}
	size_t nnames = m_abrr_name.size();
	if (idx < nnames) {
		*abrr_name = m_abrr_name[idx];
	}
	else {
		*abrr_name = "Aberration " + m_abrr_symbol[idx];
	}
	return 0;
}

int CJProbeGen::GetAberrationSymbol(int idx, string * abrr_symbol)
{
	if (NULL == abrr_symbol) { // invalid parameter reference <abrr_symbol>
		return 1;
	}
	if (idx > m_abrr_num || idx < 0) { // invalid aberration index <idx>
		return 2;
	}
	*abrr_symbol = m_abrr_symbol[idx];
	return 0;
}

int CJProbeGen::GetProbeFunctionName(int idx, string * func_name)
{
	if (NULL == func_name) { // invalid parameter reference <func_name>
		return 1;
	}
	if (idx > GetMaxProbeFunction() || idx < 0) { // invalid aberration index <idx>
		return 2;
	}
	size_t nfuncs = m_abrr_name.size();
	if (idx < nfuncs) {
		*func_name = m_pfunc_name[idx];
	}
	else {
		*func_name = "Probe Function #" + std::to_string(idx);
	}
	return 0;
}

float CJProbeGen::GetTotalIntensity(size_t n, float* pdata)
{
	// Performs a Kahan sum on the data
	// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
	double dsum = 0.0;
	double dc = 0.0;
	double dy = 0.0;
	double dt = 0.0;
	double val = 0.0;
	if (n > 0 && pdata != NULL) {
		for (size_t i = 0; i < n; i++) {
			val = (double)pdata[i];
			dy = val - dc; // next value including previous correction
			dt = dsum + dy; // intermediate new sum value
			dc = (dt - dsum) - dy; // new correction
			dsum = dt; // update result
		}
	}
	return (float)dsum;
}

void CJProbeGen::ScaleIntensity(float fsca, size_t n, float* pdata)
{
	if (n > 0 && pdata != NULL) {
		for (size_t i = 0; i < n; i++) {
			pdata[i] *= fsca;
		}
	}
}


float CJProbeGen::ApertureFunction(float qx, float qy, float qcx, float qcy, float qlim)
{
	float rtmp = 0.f;
	if (qlim>0.f) {
		float dqx = qx - qcx;
		float dqy = qy - qcy;
		if (dqx * dqx + dqy * dqy <= qlim) {
			rtmp = 1.f;
		}
	}
	return rtmp;
}


float CJProbeGen::ApertureFunctionS(float qx, float qy, float qcx, float qcy, float qlim, float ax, float ay, float smooth)
{
	// needs physical grid size: ax and ay
	// A <smooth>-pixel wide smoothing is applied to the aperture edge
	float rtmp = 0.f;
	if (qlim>0.f) {
		float dqx = qx - qcx; // qx distance from aperture center
		float dqy = qy - qcy; // qy distance from aperture center
		float dqm = sqrtf(dqx * dqx + dqy * dqy); // magnitude of distance
		if (dqm > 0.f) { // q is not at center
			float dpx = dqx * ax; // x pixel distance
			float dpy = dqy * ay; // y pixel distance
			float dpm = sqrtf(dpx * dpx + dpy * dpy); // total pixel distance
			float dlr = qlim / dqm; // ratio between q and qlim
			float dlpx = dqx * dlr * ax; // rescale to aperture edge pixel distance along x
			float dlpy = dqy * dlr * ay; // rescale to aperture edge pixel distance along y
			float dpl = sqrtf(dlpx * dlpx + dlpy * dlpy); // aperture total pixel distance
			rtmp = (1.f - tanhf((dpm - dpl)*((float)_PI) / smooth)) *0.5f; // sigmoid edge
		}
		else { // at center of aperture set amplitude to 1.
			rtmp = 1.f;
		}
	}
	return rtmp;
}


float CJProbeGen::AberrationFunction(float qx, float qy, float wl, int nabrr, float* aberrcoeff)
{
	// This code assumes that one of the cores fulfills _JMS_THRESHOLD_CALC
	if (0 == nabrr || NULL == aberrcoeff || wl <= 0.f) {
		return 0.f; // no coefficients / invalid setup
	}
	int i = 0, j = 0, k = 0, k2 = 0, l = 0, n = 0, m = 0;
	int nol = (int)(_JPG_ABERRATION_ORDER_MAX)+1;
	int nol2 = (int)(_JPG_ABERRATION_ORDER_MAX) * 2 + 1;
	int nab = min(nabrr, m_abrr_num); // limit applied list length to what is supported
	float* wfield = NULL;
	float* wabs = NULL;
	float sgnprm[4] = { 1.f,0.f,-1.f,0.f };
	float wx = qx * wl; // 1/nm -> rad
	float wy = qy * wl; // 1/nm -> rad
	float w = sqrtf(wx * wx + wy * wy);
	float prefac = (float)_TPI / wl;
	float chi = 0.f;
	float pwx = 1.f, pwy = 1.f, pw = 1.f;
	float wax = 0.f, way = 0.f;
	float rtmp = 0.f, ttmp = 0.f, ttmp1 = 0.f, tsgn = 0.f;
	wfield = (float*)malloc(sizeof(float) * 2 * nol);
	wabs = (float*)malloc(sizeof(float) * nol);
	wfield[0] = pwx; wfield[1] = pwy; wabs[0] = pw;
	for (i = 1; i < nol; i++) { // prepare list of powres of wx, wy, and w
		pwx *= wx;
		pwy *= wy;
		pw *= w;
		wfield[2 * i] = pwx;
		wfield[1 + 2 * i] = pwy;
		wabs[i] = pw;
	}
	for (k = 0; k < nab; k++) { // loop k over all aberrations in the ordered list
								// possible to implement skips here for non-active aberrations, need info
		k2 = 2 * k;
		wax = aberrcoeff[k2];
		way = aberrcoeff[1 + k2];
		if ((wax*wax + way * way) >(float)(_JPG_ABERRATION_STRTHRESH)) { // skip insignificant aberration 
			m = m_abrr_widx[k2]; // aberration order
			n = m_abrr_widx[1 + k2]; // aberration rotational symmetry
			ttmp = 0.f;
			for (l = 0; l <= n; l++) { // loop l over all terms -> possible exponents of w*
				ttmp1 = 0.f;
				// get first term sign
				tsgn = sgnprm[l % 4];
				ttmp1 = ttmp1 + tsgn * wax; // add scaling with wax
											// get second term sign
				tsgn = sgnprm[(l + 3) % 4];
				ttmp1 = ttmp1 + tsgn * way; // add scaling with way
											// get exponent j of w
				j = n - l;
				//   ++ (s1*wax + s2*way) * (n over l)           * wx^j          * wy^l
				ttmp = ttmp + ttmp1 * m_abrr_binom[n + l * nol2] * wfield[2 * j] * wfield[1 + 2 * l];
			}

			j = m - n; // get prefactor scale
			ttmp = ttmp * wabs[j] / ((float)m); // ! final scaling for current term
			rtmp = rtmp + ttmp;
		}
	}
	chi = rtmp * prefac; // scale result with prefactor to a phase [rad]
						 //
	if (NULL != wabs) free(wabs);
	if (NULL != wfield) free(wfield);
	//
	return chi;
}


int CJProbeGen::GetSourceDistribution(int src_type, float src_width, int nx, int ny, float ax, float ay, float *krn)
{
	int nerr = 0;
	int i=0, ij=0, j=0, i1=0, j1 = 0; // iterators
	int nx2 = (nx - (nx % 2)) >> 1; // x nyquist
	int ny2 = (ny - (ny % 2)) >> 1; // y nyquist
	float sx = ax / (float)nx, sy = ay / (float)ny;
	float ktot = 0.f; // total intensity of the kernel
	float kprm = 1.f, kprm2 = 1.f; // kernel parameters
	float kval = 0.f; // kernel values
	float* xn = NULL; // coordinate hash: i -> x
	float* yn = NULL; // coordinate hash: j -> y
	float y2 = 0.f, r2 = 0.f, rlim = 0.f, rlim2 = 0.f;
	float kap = 1.f; // exponent of the modiefied Bessel function
	double ksum = 0.0; // summation variable (double for cheap error prevention)
	int nitems = nx * ny; // number of pixels
	size_t nsbytes = sizeof(float)*nitems; // bytes per kernel
	memset(krn, 0, nsbytes); // preset with zeroes
	switch (src_type) {
	case 1: // Gaussian
		/*	Exp[ -(x^2 + y^2) /HWHM^2 *Log[2]]
		*/
		rlim = 3.f * src_width; // cut-off radius
		rlim2 = rlim * rlim;
		kprm = -log(2.0f) / (src_width*src_width);
		xn = (float*)malloc(sizeof(float)*nx);
		yn = (float*)malloc(sizeof(float)*ny);
		for (i = 0; i < nx; i++) { // setup horizontal index -> x
			xn[i] = sx * (((i + nx2) % nx) - nx2);
		}
		for (j = 0; j < ny; j++) { // setup vertical index -> y
			yn[j] = sy * (((j + ny2) % ny) - ny2);
		}
		for (j = 0; j < ny; j++) {
			y2 = yn[j] * yn[j];
			ij = j * nx;
			for (i = 0; i < nx; i++) {
				r2 = y2 + xn[i] * xn[i];
				if (r2 <= rlim2) { // apply radial cut-off
					kval = exp(kprm * r2);
					ksum += (double)kval;
					krn[i + ij] = kval; // write to kernel array
				}
			}
		}
		ktot = 1.f / (float)ksum;
		// normalize the kernel
		for (i = 0; i < nitems; i++) {
			krn[i] *= ktot;
		}
		free(xn);
		free(yn);
		break;
	case 2: // Modified Cauchy (Lorentzian)
		/*	A * HWHM / (2 * Pi *((A * HWHM) ^ 2 + (x - x0) ^ 2 + (y - y0) ^ 2) ^ [Kappa])
			A = Sqrt[1 / 3 (1 + 2 * 2 ^ (1 / 3) + 2 ^ (2 / 3))] = 1.3047660265041066
			Kappa = 1.5
		*/
		rlim = 10.f * src_width; // cut-off radius
		rlim2 = rlim * rlim;
		kprm = 1.3047660265f * src_width;
		kprm2 = kprm * kprm;
		xn = (float*)malloc(sizeof(float)*nx);
		yn = (float*)malloc(sizeof(float)*ny);
		for (i = 0; i < nx; i++) { // setup horizontal index -> x
			xn[i] = sx * (((i + nx2) % nx) - nx2);
		}
		for (j = 0; j < ny; j++) { // setup vertical index -> y
			yn[j] = sy * (((j + ny2) % ny) - ny2);
		}
		for (j = 0; j < ny; j++) {
			y2 = yn[j] * yn[j];
			ij = j * nx;
			for (i = 0; i < nx; i++) {
				r2 = y2 + xn[i] * xn[i];
				if (r2 <= rlim2) { // apply radial cut-off
					kval = kprm / (float)_TPI / powf(kprm2 + r2, kap);
					ksum += (double)kval;
					krn[i + ij] = kval; // write to kernel array
				}
			}
		}
		ktot = 1.f / (float)ksum;
		// normalize the kernel
		for (i = 0; i < nitems; i++) {
			krn[i] *= ktot;
		}
		free(xn);
		free(yn);
		break;
	case 3: // Disk Kernel
		rlim = 1.5f * src_width; // cut-off radius
		rlim2 = rlim * rlim;
		kprm = src_width;
		xn = (float*)malloc(sizeof(float)*nx);
		yn = (float*)malloc(sizeof(float)*ny);
		for (i = 0; i < nx; i++) { // setup horizontal index -> x
			xn[i] = sx * (((i + nx2) % nx) - nx2);
		}
		for (j = 0; j < ny; j++) { // setup vertical index -> y
			yn[j] = sy * (((j + ny2) % ny) - ny2);
		}
		for (j = 0; j < ny; j++) {
			y2 = yn[j] * yn[j];
			ij = j * nx;
			for (i = 0; i < nx; i++) {
				r2 = y2 + xn[i] * xn[i];
				if (r2 <= rlim2) { // apply radial cut-off
					kval = ApertureFunctionS(xn[i],yn[j],0.f,0.f,kprm,ax,ay);
					ksum += (double)kval;
					krn[i + ij] = kval; // write to kernel array
				}
			}
		}
		ktot = 1.f / (float)ksum;
		// normalize the kernel
		for (i = 0; i < nitems; i++) {
			krn[i] *= ktot;
		}
		free(xn);
		free(yn);
		break;
	default: // delta function
		krn[0] = 1.f;
		ktot = 1.f;
		break;
	}
	return nerr;
}


int CJProbeGen::CalculateProbeWaveFourier(CJProbeParams* prm, int nx, int ny, float ax, float ay, fcmplx *wav)
{
	// Assumes that parameters are all valid, do checks before calling this routine.
	if (NULL == wav) {
		return 1;
	}
	int nitems = nx * ny;
	int i = 0, j = 0, ij = 0;
	size_t nbytesx = sizeof(float) * nx;
	size_t nbytesy = sizeof(float) * ny;
	int nx2 = (nx - (nx % 2)) >> 1;
	int ny2 = (ny - (ny % 2)) >> 1;
	int nab = 0;
	float *qnx = NULL;
	float *qny = NULL;
	float btx = 0.f, bty = 0.f;
	float qlcx = 0.f, qlcy = 0.f, qlim = 0.f;
	float chi = 0.f;
	float wap = 0.f;
	float wpow = 0.f, fsca = 1.f;
	nab = prm->GetAberrationNum();
	qlim = prm->m_alpha * 0.001f / prm->m_wl; // mrad -> 1/nm
	qlcx = prm->m_alpha_x0 * 0.001f / prm->m_wl; // mrad -> 1/nm
	qlcy = prm->m_alpha_y0 * 0.001f / prm->m_wl; // mrad -> 1/nm
	btx = prm->m_btx * 0.001f / prm->m_wl; // mrad -> 1/nm
	bty = prm->m_bty * 0.001f / prm->m_wl; // mrad -> 1/nm
	//
	qnx = (float*)malloc(nbytesx);
	qny = (float*)malloc(nbytesy);
	memset(qnx, 0, nbytesx);
	memset(qny, 0, nbytesy);
	memset(wav, 0, sizeof(fcmplx)*nx*ny);
	for (i = 0; i < nx; i++) { // setup horizontal index -> qx + btx
		qnx[i] = (float)(((i + nx2) % nx) - nx2) / ax + btx;
	}
	for (j = 0; j < ny; j++) { // setup vertical index -> qy + bty
		qny[j] = (float)(((j + ny2) % ny) - ny2) / ay + bty;
	}
	//
	for (j = 0; j < ny; j++) { // loop over grid rows
		ij = j * nx;
		for (i = 0; i < nx; i++) { // loop over grid columns
			wap = ApertureFunctionS(qnx[i], qny[j], qlcx, qlcy, qlim, ax, ay); // aperture with smoothing of 1 pixel
			if (wap >(float)_JPG_PROBE_APERTURE_THRESH) {
				if (nab > 0) {
					chi = AberrationFunction(qnx[i], qny[j], prm->m_wl, nab, prm->m_abrr_coeff);
				}
				wav[i + ij] = fcmplx(wap*cosf(chi), -wap * sinf(chi));
				wpow += (wap*wap);
			}
		}
	}
	// normalize
	fsca = 1.f / sqrtf(wpow);
	for (j = 0; j < ny; j++) { // loop over grid rows
		ij = j * nx;
		for (i = 0; i < nx; i++) { // loop over grid columns
			wav[i + ij] = wav[i + ij] * fsca; // apply normalization factor
		}
	}
	//
	if (NULL != qnx) free(qnx);
	if (NULL != qny) free(qny);
	return 0;
}


int CJProbeGen::CalculateProbeIntensityCoh(CJProbeParams* prm, int ndim, float s, float* pdata)
{
	// assuming that all parameters are valid. Check before calling!
	int nerr = 0; // error code
	if (NULL == pdata) { // always check the validity of the pointer
		return 1;
	}
	int nft[2] = { ndim,ndim };
	CJFFTWcore jfft;
	jfft.Init(2, nft);
	size_t nbytes = sizeof(float) * ndim; // number of bytes per row and column
	size_t nsbytes = sizeof(float) * ndim * ndim; // number of bytes in output
	memset(pdata, 0, nsbytes);
	CJProbeParams lprm = *prm; // get a copy of the probe parameters
	int nab = lprm.GetAberrationNum(); // number of aberrations
	int ndim2 = (ndim - (ndim % 2)) >> 1; // nyquist pixel index
	int i = 0, j = 0, ij = 0; // iterators
	float a = s * ndim;  // physical box size(x and y)
	float sq = 1.f / a; // q sampling = 1 / physical box size (x and y)
	float *qn = (float*)malloc(nbytes); // allocate frequency helper
	memset(qn, 0, nbytes);
	fcmplx *wav = (fcmplx*)malloc(nsbytes * 2); // allocate buffer holding the wave function (2 floats per pixel)
	memset(wav, 0, nsbytes * 2);
	float btx = lprm.m_btx * 0.001f / lprm.m_wl; // beam tilt X: mrad -> 1/nm
	float bty = lprm.m_bty * 0.001f / lprm.m_wl; // beam tilt Y: mrad -> 1/nm
	float qlim = lprm.m_alpha * 0.001f / lprm.m_wl; // aperture size: mrad -> 1/nm
	float qlcx = lprm.m_alpha_x0 * 0.001f / lprm.m_wl; // aperture shift X: mrad -> 1/nm
	float qlcy = lprm.m_alpha_y0 * 0.001f / lprm.m_wl; // aperture shift Y: mrad -> 1/nm
	lprm.m_abrr_coeff[0] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	lprm.m_abrr_coeff[1] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	float chi = 0.f; // aberration function value
	float wap = 0.f; // aperture value
	float wpow = 0.f; // accumulated power
	float fsca = 1.f; // normalization
	float qx = 0.f, qy = 0.f; // coordinates
	for (i = 0; i < ndim; i++) { // setup frequency hash: i -> q (not scrambled!!!)
		qn[i] = sq * (((i + ndim2) % ndim) - ndim2);
	}
	// set Fourier coefficients in wav
	for (j = 0; j < ndim; j++) { // loop over grid rows
		qy = qn[j] + bty;
		ij = j * ndim;
		for (i = 0; i < ndim; i++) { // loop over grid columns
			qx = qn[i] + btx;
			wap = ApertureFunctionS(qx, qy, qlcx, qlcy, qlim, a, a); // aperture with smoothing of 1 pixel
			if (wap >(float)_JPG_PROBE_APERTURE_THRESH) {
				if (nab > 0) {
					chi = AberrationFunction(qx, qy, lprm.m_wl, nab, lprm.m_abrr_coeff);
				}
				wav[i + ij] = fcmplx(wap * cosf(chi), -wap * sinf(chi));
			}
		}
	}
	// inverse fourier transform
	jfft.SetDataC(wav); // copy data to fft module
	jfft.IFT(); // make an inverse fft
	jfft.GetDataPow(pdata); // copy the result to the output array
	wpow = GetTotalIntensity(ndim*ndim, pdata); // get total probe intensity
	ScaleIntensity(1.f/wpow, ndim * ndim, pdata); // normalize to total intensity 1
	jfft.Deinit(); // de-init the fft module
	free(wav); // free memory of wave helper
	free(qn); // free memory of q hash
	//
	return nerr;
}



int CJProbeGen::CalculateProbeIntensityPSC(CJProbeParams* prm, int ndim, float s, float* pdata)
{
	// assuming that all parameters are valid. Check before calling!
	int nerr = 0; // error code
	if (NULL == pdata) { // always check the validity of the pointer
		return 1;
	}
	int nft[2] = { ndim,ndim };
	CJFFTWcore jfft;
	jfft.Init(2, nft, FFTW_MEASURE);
	size_t nbytes = sizeof(float) * ndim; // number of bytes per row and column
	size_t nsbytes = sizeof(float) * ndim * ndim; // number of bytes in output
	memset(pdata, 0, nsbytes);
	CJProbeParams lprm = *prm; // get a copy of the probe parameters
	int nab = lprm.GetAberrationNum(); // number of aberrations
	int ndim2 = (ndim - (ndim % 2)) >> 1; // nyquist pixel index
	int i = 0, j = 0, ij = 0; // iterators
	float a = s * ndim;  // physical box size(x and y)
	float sq = 1.f / a; // q sampling = 1 / physical box size (x and y)
	float *qn = (float*)malloc(nbytes); // allocate frequency helper
	memset(qn, 0, nbytes);
	fcmplx *wav = (fcmplx*)malloc(nsbytes * 2); // allocate buffer holding the wave function (2 floats per pixel)
	memset(wav, 0, nsbytes * 2);
	float srcw = lprm.m_source_width; // source radius [nm] (HWHM)
	if (fabs(srcw) < 0.001f) { // extremely small or invalid source size
		lprm.m_source_shape = 0; // turn the source convolution off
	}
	int srct = max(0, min(_JPG_SOURCE_TYPE_MAX, lprm.m_source_shape)); // source type, limited to supported range
	float btx = lprm.m_btx * 0.001f / lprm.m_wl; // beam tilt X: mrad -> 1/nm
	float bty = lprm.m_bty * 0.001f / lprm.m_wl; // beam tilt Y: mrad -> 1/nm
	float qlim = lprm.m_alpha * 0.001f / lprm.m_wl; // aperture size: mrad -> 1/nm
	float qlcx = lprm.m_alpha_x0 * 0.001f / lprm.m_wl; // aperture shift X: mrad -> 1/nm
	float qlcy = lprm.m_alpha_y0 * 0.001f / lprm.m_wl; // aperture shift Y: mrad -> 1/nm
	lprm.m_abrr_coeff[0] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	lprm.m_abrr_coeff[1] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	float chi = 0.f; // aberration function value
	float wap = 0.f; // aperture value
	float ptot = 0.f; // total probe intensity for normalization
	float fsca = 1.f; // normalization
	float qx = 0.f, qy = 0.f; // coordinates
	for (i = 0; i < ndim; i++) { // setup frequency hash: i -> q
		qn[i] = sq * (((i + ndim2) % ndim) - ndim2);
	}
	// set Fourier coefficients in wav
	for (j = 0; j < ndim; j++) { // loop over grid rows
		qy = qn[j] + bty;
		ij = j * ndim;
		for (i = 0; i < ndim; i++) { // loop over grid columns
			qx = qn[i] + btx;
			wap = ApertureFunctionS(qx, qy, qlcx, qlcy, qlim, a, a); // aperture with smoothing of 1 pixel
			if (wap >(float)_JPG_PROBE_APERTURE_THRESH) {
				if (nab > 0) {
					chi = AberrationFunction(qx, qy, lprm.m_wl, nab, lprm.m_abrr_coeff);
				}
				wav[i + ij] = fcmplx(wap * cosf(chi), -wap * sinf(chi));
			}
		}
	}
	// inverse fourier transform
	jfft.SetDataC(wav); // copy data to fft module
	jfft.IFT(); // make an inverse fft
	ptot = jfft.GetDataTotalPow();
	jfft.Scale(1.f / sqrtf(ptot)); // -> scale power of the coherent wave function to 1
	jfft.GetDataPow(pdata); // copy the result to the output array (still not normalized)
	// apply partial spatial coherence by convolving the coherent probe with the source distribution
	if (srct > 0) { // do source size convolution of the probe intensity if source type is not a delta function
		jfft.Zero(); // zeroe the data
		jfft.SetDataRe(pdata); // set the probe intensity as new data
		jfft.FT(); // transform to Fourier-space
		jfft.GetDataC(wav); // store the probe-ft in wav
		ptot = jfft.GetDataTotalPow();
		float* krn = (float*)malloc(nsbytes); // allocate kernel memory with ndim*ndim floats
		nerr = GetSourceDistribution(srct, srcw, ndim, ndim, a, a, krn);
		jfft.Zero(); // zeroe again
		jfft.SetDataRe(krn); // set kernel
		ptot = jfft.GetDataTotalRe();
		jfft.FT(); // transform to Fourier space
		jfft.Conjugate(); // take complex conjugate
		jfft.MultiplyC(wav); // multiply with probe-ft
		ptot = jfft.GetDataTotalPow();
		jfft.IFT(); // inverse FT back to real-space
		ptot = jfft.GetDataTotalPow();
		jfft.GetDataRe(pdata); // get the real part as result of the convolution
		free(krn); // free memory of the kernel
	}
	// Normalize the data
	ptot = GetTotalIntensity(ndim*ndim, pdata);
	ScaleIntensity(1.f / ptot, ndim*ndim, pdata);
	jfft.Deinit();
	free(wav); // free memory of wave helper
	free(qn); // free memory of q hash
	return nerr;
}


int CJProbeGen::CalculateProbeIntensityPTC(CJProbeParams* prm, int ndim, float s, float* pdata)
{
	// assuming that all parameters are valid. Check before calling!
	int nerr = 0; // error code
	if (NULL == pdata) { // always check the validity of the pointer
		return 1;
	}
	int nft[2] = { ndim,ndim };
	CJFFTWcore jfft;
	jfft.Init(2, nft, FFTW_MEASURE);
	size_t nbytes = sizeof(float) * ndim; // number of bytes per row and column
	size_t nsbytes = sizeof(float) * ndim * ndim; // number of bytes in output
	memset(pdata, 0, nsbytes);
	CJProbeParams lprm = *prm; // get a copy of the probe parameters
	int nab = lprm.GetAberrationNum(); // number of aberrations
	int ndim2 = (ndim - (ndim % 2)) >> 1; // nyquist pixel index
	int i = 0, j = 0, ij = 0, ij3 = 0, k = 0; // iterators
	float a = s * ndim;  // physical box size(x and y)
	float sq = 1.f / a; // q sampling = 1 / physical box size (x and y)
	float *qn = (float*)malloc(nbytes); // allocate frequency helper
	memset(qn, 0, nbytes);
	fcmplx *wav = (fcmplx*)malloc(nsbytes * 2); // allocate buffer holding the wave function (2 floats per pixel)
	memset(wav, 0, nsbytes * 2);
	float *tmp3 = (float*)malloc(nsbytes * 3); // allocate buffer holding the aberration function, focal term and aperture (3 floats per pixel)
	memset(tmp3, 0, nsbytes * 3);
	float *probe = (float*)malloc(nsbytes); // allocate buffer for intermediate probe inetnsities
	memset(probe, 0, nsbytes);
	float fsw = lprm.m_fspread_width;
	if (fabs(fsw) < 0.01f) { // extremely small or invalid focus spread
		lprm.m_fspread_kernel_samples = 1; // turn the focus spread convolution off
	}
	int fsnum = max(1, lprm.m_fspread_kernel_samples);
	float fsrw = fabs(lprm.m_fspread_kernel_width); // relative width of the focus spread kernel
	float btx = lprm.m_btx * 0.001f / lprm.m_wl; // beam tilt X: mrad -> 1/nm
	float bty = lprm.m_bty * 0.001f / lprm.m_wl; // beam tilt Y: mrad -> 1/nm
	float qlim = lprm.m_alpha * 0.001f / lprm.m_wl; // aperture size: mrad -> 1/nm
	float qlcx = lprm.m_alpha_x0 * 0.001f / lprm.m_wl; // aperture shift X: mrad -> 1/nm
	float qlcy = lprm.m_alpha_y0 * 0.001f / lprm.m_wl; // aperture shift Y: mrad -> 1/nm
	lprm.m_abrr_coeff[0] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	lprm.m_abrr_coeff[1] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	float fswgt = 1.f; // focal kernel weight
	float dzval = 0.f; // addional defocus in focal averaging
	float chi = 0.f; // aberration function value
	float wap = 0.f; // aperture value
	float zchi = 0.f; // focal phase term
	float ptot = 0.f; // total probe intensity for normalization
	float fsca = 1.f; // normalization
	float qx = 0.f, qy = 0.f; // coordinates
	float* ptmp3 = NULL;
	for (i = 0; i < ndim; i++) { // setup frequency hash: i -> q (not scrambled!!!)
		qn[i] = sq * (((i + ndim2) % ndim) - ndim2);
	}
	// set basis aberration function and the extra focal term in chi2
	for (j = 0; j < ndim; j++) { // loop over grid rows
		qy = qn[j] + bty;
		ij3 = j * ndim * 3;
		for (i = 0; i < ndim; i++) { // loop over grid columns
			ptmp3 = &tmp3[3 * i + ij3];
			qx = qn[i] + btx;
			wap = ApertureFunctionS(qx, qy, qlcx, qlcy, qlim, a, a); // aperture with smoothing of 1 pixel
			chi = 0.f;
			zchi = 0.f;
			if (wap >(float)_JPG_PROBE_APERTURE_THRESH) {
				if (nab > 0) {
					chi = AberrationFunction(qx, qy, lprm.m_wl, nab, lprm.m_abrr_coeff); // aberration function
				}
				zchi = (float)_PI*lprm.m_wl*(qx * qx + qy * qy); // Pi * Lambda * g^2
			}
			ptmp3[0] = wap; // store aperture
			ptmp3[1] = chi; // store aberration function
			ptmp3[2] = zchi; // store focal term
		}
	}
	// explicit focal averaging
	for (k = 0; k < fsnum; k++) {
		if (fsnum > 1) {
			dzval = 2.f * fsw * fsrw * ((float)k / (fsnum - 1) - 0.5f);
			fswgt = exp(-dzval * dzval / (fsw*fsw));
		}
		else {
			dzval = 0.f;
			fswgt = 1.f;
		}
		// calculate the coherent probe intensity for this defocus (dzval)
		for (j = 0; j < ndim; j++) { // loop over grid rows
			ij = j * ndim;
			ij3 = ij * 3;
			for (i = 0; i < ndim; i++) { // loop over grid columns
				ptmp3 = &tmp3[3 * i + ij3];
				wap = ptmp3[0]; // get aperture
				if (wap >(float)_JPG_PROBE_APERTURE_THRESH) {
					chi = ptmp3[1] + dzval * ptmp3[2]; // get aberration function with extra focal term
					wav[i + ij] = fcmplx(wap * cosf(chi), -wap * sinf(chi)); // set probe wave function coefficient
				}
			}
		}
		// inverse fourier transform
		jfft.SetDataC(wav); // copy data to fft module
		jfft.IFT(); // make an inverse fft
		jfft.GetDataPow(probe); // copy the result to the intermediate buffer
		for (size_t ls = 0; ls < ndim*ndim; ls++) { // accumulate the probe to the ouput array (not normalized)
			pdata[ls] += (fswgt*probe[ls]); // add probe intensities with focal kernel weight to pdata
		}
	}
	// Normalize the probe intensity
	ptot = GetTotalIntensity(ndim*ndim, pdata);
	ScaleIntensity(1.f / ptot, ndim*ndim, pdata); // with this call pdata contains the final result
	jfft.Deinit();
	free(probe); // free memory of the intermediate probe buffer
	free(tmp3); // free memory of the temp buffer
	free(wav); // free memory of wave helper
	free(qn); // free memory of q hash
	return nerr;
}



int CJProbeGen::CalculateProbeIntensity(CJProbeParams *prm, int ndim, float s, float* pdata)
{
	// assuming that all parameters are valid. Check before calling!
	int nerr = 0; // error code
	if (NULL == pdata) { // always check the validity of the pointer
		return 1;
	}
	int nft[2] = { ndim,ndim };
	CJFFTWcore jfft;
	jfft.Init(2, nft, FFTW_MEASURE);
	size_t nbytes = sizeof(float) * ndim; // number of bytes per row and column
	size_t nsbytes = sizeof(float) * ndim * ndim; // number of bytes in output
	memset(pdata, 0, nsbytes);
	CJProbeParams lprm = *prm; // get a copy of the probe parameters
	int nab = lprm.GetAberrationNum(); // number of aberrations
	int ndim2 = (ndim - (ndim % 2)) >> 1; // nyquist pixel index
	int i = 0, j = 0, ij = 0, ij3 = 0, k = 0; // iterators
	float a = s * ndim;  // physical box size(x and y)
	float sq = 1.f / a; // q sampling = 1 / physical box size (x and y)
	float *qn = (float*)malloc(nbytes); // allocate frequency helper
	memset(qn, 0, nbytes);
	fcmplx *wav = (fcmplx*)malloc(nsbytes * 2); // allocate buffer holding the wave function (2 floats per pixel)
	memset(wav, 0, nsbytes * 2);
	float *tmp3 = (float*)malloc(nsbytes * 3); // allocate buffer holding the aberration function, focal term and aperture (3 floats per pixel)
	memset(tmp3, 0, nsbytes * 3);
	float *probe = (float*)malloc(nsbytes); // allocate buffer for intermediate probe inetnsities
	memset(probe, 0, nsbytes);
	float srcw = lprm.m_source_width; // source radius [nm] (HWHM)
	if (fabs(srcw) < 0.001f) { // extremely small or invalid source size
		lprm.m_source_shape = 0; // turn the source convolution off
	}
	int srct = max(0, min(_JPG_SOURCE_TYPE_MAX, lprm.m_source_shape)); // source type, limited to supported range
	float fsw = lprm.m_fspread_width;
	if (fabs(fsw) < 0.01f) { // extremely small or invalid focus spread
		lprm.m_fspread_kernel_samples = 1; // turn the focus spread convolution off
	}
	int fsnum = max(1, lprm.m_fspread_kernel_samples);
	float fsrw = fabs(lprm.m_fspread_kernel_width); // relative width of the focus spread kernel
	float btx = lprm.m_btx * 0.001f / lprm.m_wl; // beam tilt X: mrad -> 1/nm
	float bty = lprm.m_bty * 0.001f / lprm.m_wl; // beam tilt Y: mrad -> 1/nm
	float qlim = lprm.m_alpha * 0.001f / lprm.m_wl; // aperture size: mrad -> 1/nm
	float qlcx = lprm.m_alpha_x0 * 0.001f / lprm.m_wl; // aperture shift X: mrad -> 1/nm
	float qlcy = lprm.m_alpha_y0 * 0.001f / lprm.m_wl; // aperture shift Y: mrad -> 1/nm
	lprm.m_abrr_coeff[0] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	lprm.m_abrr_coeff[1] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	float fswgt = 1.f; // focal kernel weight
	float dzval = 0.f; // addional defocus in focal averaging
	float chi = 0.f; // aberration function value
	float wap = 0.f; // aperture value
	float zchi = 0.f; // focal phase term
	float ptot = 0.f; // total probe intensity for normalization
	float fsca = 1.f; // normalization
	float qx = 0.f, qy = 0.f; // coordinates
	float* ptmp3 = NULL;
	for (i = 0; i < ndim; i++) { // setup frequency hash: i -> q (not scrambled!!!)
		qn[i] = sq * (((i + ndim2) % ndim) - ndim2);
	}
	// set basis aberration function and the extra focal term in chi2
	for (j = 0; j < ndim; j++) { // loop over grid rows
		qy = qn[j] + bty;
		ij3 = j * ndim * 3;
		for (i = 0; i < ndim; i++) { // loop over grid columns
			ptmp3 = &tmp3[3 * i + ij3];
			qx = qn[i] + btx;
			wap = ApertureFunctionS(qx, qy, qlcx, qlcy, qlim, a, a); // aperture with smoothing of 1 pixel
			chi = 0.f;
			zchi = 0.f;
			if (wap > (float)_JPG_PROBE_APERTURE_THRESH) {
				if (nab > 0) {
					chi = AberrationFunction(qx, qy, lprm.m_wl, nab, lprm.m_abrr_coeff); // aberration function
				}
				zchi = (float)_PI*lprm.m_wl*(qx * qx + qy * qy); // Pi * Lambda * g^2
			}
			ptmp3[0] = wap; // store aperture
			ptmp3[1] = chi; // store aberration function
			ptmp3[2] = zchi; // store focal term
		}
	}
	// explicit focal averaging
	for (k = 0; k < fsnum; k++) {
		if (fsnum > 1) {
			dzval = 2.f * fsw * fsrw * ( (float)k / (fsnum-1) - 0.5f );
			fswgt = exp(-dzval * dzval / (fsw*fsw));
		}
		else {
			dzval = 0.f;
			fswgt = 1.f;
		}
		// calculate the coherent probe intensity for this defocus (dzval)
		for (j = 0; j < ndim; j++) { // loop over grid rows
			ij = j * ndim;
			ij3 = ij * 3;
			for (i = 0; i < ndim; i++) { // loop over grid columns
				ptmp3 = &tmp3[3 * i + ij3];
				wap = ptmp3[0]; // get aperture
				if (wap >(float)_JPG_PROBE_APERTURE_THRESH) {
					chi = ptmp3[1] + dzval * ptmp3[2]; // get aberration function with extra focal term
					wav[i + ij] = fcmplx(wap * cosf(chi), -wap * sinf(chi)); // set probe wave function coefficient
				}
			}
		}
		// inverse fourier transform
		jfft.SetDataC(wav); // copy data to fft module
		jfft.IFT(); // make an inverse fft
		jfft.GetDataPow(probe); // copy the result to the intermediate buffer
		for (size_t ls = 0; ls < ndim*ndim; ls++) { // accumulate the probe to the ouput array (not normalized)
			pdata[ls] += (fswgt*probe[ls]); // add probe intensities with focal kernel weight to pdata
		}
	}
	// apply partial spatial coherence by convolving the coherent probe with the source distribution
	if (srct > 0) { // do source size convolution of the probe intensity if source type is not a delta function
		jfft.Zero(); // zeroe the data
		jfft.SetDataRe(pdata); // set the probe intensity as new data
		jfft.FT(); // transform to Fourier-space
		jfft.GetDataC(wav); // store the probe-ft in wav
		float* krn = (float*)malloc(nsbytes); // allocate kernel memory with ndim*ndim floats
		nerr = GetSourceDistribution(srct, srcw, ndim, ndim, a, a, krn);
		jfft.Zero(); // zeroe again
		jfft.SetDataRe(krn); // set kernel
		jfft.FT(); // transform to Fourier space
		jfft.Conjugate(); // take complex conjugate
		jfft.MultiplyC(wav); // multiply with probe-ft
		jfft.IFT(); // inverse FT back to real-space
		jfft.GetDataRe(pdata); // get the real part as result of the convolution
		free(krn); // free memory of the kernel
	}
	// Normalize the probe intensity
	ptot = GetTotalIntensity(ndim*ndim, pdata);
	ScaleIntensity(1.f / ptot, ndim*ndim, pdata); // with this call pdata contains the final result
	jfft.Deinit();
	free(probe); // free memory of the intermediate probe buffer
	free(tmp3); // free memory of the temp buffer
	free(wav); // free memory of wave helper
	free(qn); // free memory of q hash
	return nerr;
}



int CJProbeGen::CalculateProbePhase(CJProbeParams* prm, int ndim, float s, float* pdata)
{
	// assuming that all parameters are valid. Check before calling!
	int nerr = 0; // error code
	if (NULL == pdata) { // always check the validity of the pointer
		return 1;
	}
	int nft[2] = { ndim,ndim };
	CJFFTWcore jfft;
	jfft.Init(2, nft);
	size_t nbytes = sizeof(float) * ndim; // number of bytes per row and column
	size_t nsbytes = sizeof(float) * ndim * ndim; // number of bytes in output
	memset(pdata, 0, nsbytes);
	CJProbeParams lprm = *prm; // get a copy of the probe parameters
	int nab = lprm.GetAberrationNum(); // number of aberrations
	int ndim2 = (ndim - (ndim % 2)) >> 1; // nyquist pixel index
	int i = 0, j = 0, ij = 0; // iterators
	float a = s * ndim;  // physical box size(x and y)
	float sq = 1.f / a; // q sampling = 1 / physical box size (x and y)
	float *qn = (float*)malloc(nbytes); // allocate frequency helper
	memset(qn, 0, nbytes);
	fcmplx *wav = (fcmplx*)malloc(nsbytes * 2); // allocate buffer holding the wave function (2 floats per pixel)
	memset(wav, 0, nsbytes * 2);
	float btx = lprm.m_btx * 0.001f / lprm.m_wl; // beam tilt X: mrad -> 1/nm
	float bty = lprm.m_bty * 0.001f / lprm.m_wl; // beam tilt Y: mrad -> 1/nm
	float qlim = lprm.m_alpha * 0.001f / lprm.m_wl; // aperture size: mrad -> 1/nm
	float qlcx = lprm.m_alpha_x0 * 0.001f / lprm.m_wl; // aperture shift X: mrad -> 1/nm
	float qlcy = lprm.m_alpha_y0 * 0.001f / lprm.m_wl; // aperture shift Y: mrad -> 1/nm
	lprm.m_abrr_coeff[0] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	lprm.m_abrr_coeff[1] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	float chi = 0.f; // aberration function value
	float wap = 0.f; // aperture value
	// float wpow = 0.f; // accumulated power
	float fsca = 1.f; // normalization
	float qx = 0.f, qy = 0.f; // coordinates
	for (i = 0; i < ndim; i++) { // setup frequency hash: i -> q (not scrambled!!!)
		qn[i] = sq * (((i + ndim2) % ndim) - ndim2);
	}
	// set Fourier coefficients in wav
	for (j = 0; j < ndim; j++) { // loop over grid rows
		qy = qn[j] + bty;
		ij = j * ndim;
		for (i = 0; i < ndim; i++) { // loop over grid columns
			qx = qn[i] + btx;
			wap = ApertureFunctionS(qx, qy, qlcx, qlcy, qlim, a, a); // aperture with smoothing of 1 pixel
			if (wap >(float)_JPG_PROBE_APERTURE_THRESH) {
				if (nab > 0) {
					chi = AberrationFunction(qx, qy, lprm.m_wl, nab, lprm.m_abrr_coeff);
				}
				wav[i + ij] = fcmplx(wap * cosf(chi), -wap * sinf(chi));
			}
		}
	}
	// inverse fourier transform
	jfft.SetDataC(wav); // copy data to fft module
	jfft.IFT(); // make an inverse fft
	// The probe doesn't need to be normalized to measure the phase
	// wpow = jfft.GetDataTotalPow(); // get the total power of the result
	// jfft.Scale(1.f / sqrt(wpow)); // rescale the wave function to normalize to power = 1
	jfft.GetDataArg(pdata); // copy the result to the output array
	jfft.Deinit(); // de-init the fft module
	free(wav); // free memory of wave helper
	free(qn); // free memory of q hash
			  //
	return nerr;
}


int CJProbeGen::CalculateProbeRe(CJProbeParams *prm, int ndim, float s, float* pdata)
{
	// assuming that all parameters are valid. Check before calling!
	int nerr = 0; // error code
	if (NULL == pdata) { // always check the validity of the pointer
		return 1;
	}
	int nft[2] = { ndim,ndim };
	CJFFTWcore jfft;
	jfft.Init(2, nft);
	size_t nbytes = sizeof(float) * ndim; // number of bytes per row and column
	size_t nsbytes = sizeof(float) * ndim * ndim; // number of bytes in output
	memset(pdata, 0, nsbytes);
	CJProbeParams lprm = *prm; // get a copy of the probe parameters
	int nab = lprm.GetAberrationNum(); // number of aberrations
	int ndim2 = (ndim - (ndim % 2)) >> 1; // nyquist pixel index
	int i = 0, j = 0, ij = 0; // iterators
	float a = s * ndim;  // physical box size(x and y)
	float sq = 1.f / a; // q sampling = 1 / physical box size (x and y)
	float *qn = (float*)malloc(nbytes); // allocate frequency helper
	memset(qn, 0, nbytes);
	fcmplx *wav = (fcmplx*)malloc(nsbytes * 2); // allocate buffer holding the wave function (2 floats per pixel)
	memset(wav, 0, nsbytes * 2);
	float btx = lprm.m_btx * 0.001f / lprm.m_wl; // beam tilt X: mrad -> 1/nm
	float bty = lprm.m_bty * 0.001f / lprm.m_wl; // beam tilt Y: mrad -> 1/nm
	float qlim = lprm.m_alpha * 0.001f / lprm.m_wl; // aperture size: mrad -> 1/nm
	float qlcx = lprm.m_alpha_x0 * 0.001f / lprm.m_wl; // aperture shift X: mrad -> 1/nm
	float qlcy = lprm.m_alpha_y0 * 0.001f / lprm.m_wl; // aperture shift Y: mrad -> 1/nm
	lprm.m_abrr_coeff[0] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	lprm.m_abrr_coeff[1] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	float chi = 0.f; // aberration function value
	float wap = 0.f; // aperture value
	float wpow = 0.f; // accumulated power
	float fsca = 1.f; // normalization
	float qx = 0.f, qy = 0.f; // coordinates
	for (i = 0; i < ndim; i++) { // setup frequency hash: i -> q (not scrambled!!!)
		qn[i] = sq * (((i + ndim2) % ndim) - ndim2);
	}
	// set Fourier coefficients in wav
	for (j = 0; j < ndim; j++) { // loop over grid rows
		qy = qn[j] + bty;
		ij = j * ndim;
		for (i = 0; i < ndim; i++) { // loop over grid columns
			qx = qn[i] + btx;
			wap = ApertureFunctionS(qx, qy, qlcx, qlcy, qlim, a, a); // aperture with smoothing of 1 pixel
			if (wap >(float)_JPG_PROBE_APERTURE_THRESH) {
				if (nab > 0) {
					chi = AberrationFunction(qx, qy, lprm.m_wl, nab, lprm.m_abrr_coeff);
				}
				wav[i + ij] = fcmplx(wap * cosf(chi), -wap * sinf(chi));
			}
		}
	}
	// inverse fourier transform
	jfft.SetDataC(wav); // copy data to fft module
	jfft.IFT(); // make an inverse fft
	wpow = jfft.GetDataTotalPow(); // get the total power of the result
	jfft.Scale(1.f / sqrt(wpow)); // rescale the wave function to normalize to power = 1
	jfft.GetDataRe(pdata); // copy the result to the output array
	jfft.Deinit(); // de-init the fft module
	free(wav); // free memory of wave helper
	free(qn); // free memory of q hash
			  //
	return nerr;
}


int CJProbeGen::CalculateProbeIm(CJProbeParams *prm, int ndim, float s, float* pdata)
{
	// assuming that all parameters are valid. Check before calling!
	int nerr = 0; // error code
	if (NULL == pdata) { // always check the validity of the pointer
		return 1;
	}
	int nft[2] = { ndim,ndim };
	CJFFTWcore jfft;
	jfft.Init(2, nft);
	size_t nbytes = sizeof(float) * ndim; // number of bytes per row and column
	size_t nsbytes = sizeof(float) * ndim * ndim; // number of bytes in output
	memset(pdata, 0, nsbytes);
	CJProbeParams lprm = *prm; // get a copy of the probe parameters
	int nab = lprm.GetAberrationNum(); // number of aberrations
	int ndim2 = (ndim - (ndim % 2)) >> 1; // nyquist pixel index
	int i = 0, j = 0, ij = 0; // iterators
	float a = s * ndim;  // physical box size(x and y)
	float sq = 1.f / a; // q sampling = 1 / physical box size (x and y)
	float *qn = (float*)malloc(nbytes); // allocate frequency helper
	memset(qn, 0, nbytes);
	fcmplx *wav = (fcmplx*)malloc(nsbytes * 2); // allocate buffer holding the wave function (2 floats per pixel)
	memset(wav, 0, nsbytes * 2);
	float btx = lprm.m_btx * 0.001f / lprm.m_wl; // beam tilt X: mrad -> 1/nm
	float bty = lprm.m_bty * 0.001f / lprm.m_wl; // beam tilt Y: mrad -> 1/nm
	float qlim = lprm.m_alpha * 0.001f / lprm.m_wl; // aperture size: mrad -> 1/nm
	float qlcx = lprm.m_alpha_x0 * 0.001f / lprm.m_wl; // aperture shift X: mrad -> 1/nm
	float qlcy = lprm.m_alpha_y0 * 0.001f / lprm.m_wl; // aperture shift Y: mrad -> 1/nm
	lprm.m_abrr_coeff[0] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	lprm.m_abrr_coeff[1] += 0.5f * a; // apply an image shift of 1/2 of the box to place the probe in the center
	float chi = 0.f; // aberration function value
	float wap = 0.f; // aperture value
	float wpow = 0.f; // accumulated power
	float fsca = 1.f; // normalization
	float qx = 0.f, qy = 0.f; // coordinates
	for (i = 0; i < ndim; i++) { // setup frequency hash: i -> q (not scrambled!!!)
		qn[i] = sq * (((i + ndim2) % ndim) - ndim2);
	}
	// set Fourier coefficients in wav
	for (j = 0; j < ndim; j++) { // loop over grid rows
		qy = qn[j] + bty;
		ij = j * ndim;
		for (i = 0; i < ndim; i++) { // loop over grid columns
			qx = qn[i] + btx;
			wap = ApertureFunctionS(qx, qy, qlcx, qlcy, qlim, a, a); // aperture with smoothing of 1 pixel
			if (wap >(float)_JPG_PROBE_APERTURE_THRESH) {
				if (nab > 0) {
					chi = AberrationFunction(qx, qy, lprm.m_wl, nab, lprm.m_abrr_coeff);
				}
				wav[i + ij] = fcmplx(wap * cosf(chi), -wap * sinf(chi));
			}
		}
	}
	// inverse fourier transform
	jfft.SetDataC(wav); // copy data to fft module
	jfft.IFT(); // make an inverse fft
	wpow = jfft.GetDataTotalPow(); // get the total power of the result
	jfft.Scale(1.f / sqrt(wpow)); // rescale the wave function to normalize to power = 1
	jfft.GetDataIm(pdata); // copy the result to the output array
	jfft.Deinit(); // de-init the fft module
	free(wav); // free memory of wave helper
	free(qn); // free memory of q hash
			  //
	return nerr;
}


int CJProbeGen::CalculateAberrationPhasePlate(CJProbeParams* prm, int ndim, float s, float* pdata)
{
	// assuming that all parameters are valid. Check before calling!
	int nerr = 0; // error code
	if (NULL == pdata) { // always check the validity of the pointer
		return 1;
	}
	memset(pdata, 0, sizeof(float) * ndim * ndim);
	int nab = prm->GetAberrationNum(); // number of aberrations
	if (nab > 0) {
		size_t nbytes = sizeof(float) * ndim; // number of bytes per row and column
		int ndim2 = (ndim - (ndim % 2)) >> 1; // nyquist pixel index
		int i = 0, j = 0, ij = 0; // iterators
		float sq = 1.f / (s * ndim); // q sampling = 1 / physical box size (x and y)
		float *qn = (float*)malloc(nbytes); // allocate frequency helper
		float btx = prm->m_btx * 0.001f / prm->m_wl; // beam tilt X: mrad -> 1/nm
		float bty = prm->m_bty * 0.001f / prm->m_wl; // beam tilt Y: mrad -> 1/nm
		//
		for (i = 0; i < ndim; i++) { // setup frequency hash: i -> q (not scrambled!!!)
			qn[i] = sq * (i - ndim2);
		}
		// set aberration phases
		for (j = 0; j < ndim; j++) { // loop over grid rows
			ij = j * ndim;
			for (i = 0; i < ndim; i++) { // loop over grid columns
				pdata[i+ij] = AberrationFunction(qn[i]+btx, qn[j]+bty, prm->m_wl, nab, prm->m_abrr_coeff);
			}
		}
		free(qn);
	}
	//
	return nerr;
}


int CJProbeGen::SetAmorph(int ndim, float s, bool bForce)
{
	int i = 0, j = 0, ij = 0, ndim2 = 0;
	unsigned int urnd = 0;
	size_t nitems = 0, ls = 0;
	int nerr = 0;
	int newpgr = 0;
	if (ndim <= 0 || s <= 0.f) {
		return 1;
	}
	if (NULL == m_amorph_pgr || ndim != m_amorph_dim || fabs(s - m_amorph_samp) > 0.0001f || bForce) {
		// we are here because we need a new phase grating
		if (NULL != m_amorph_pgr && ndim != m_amorph_dim) { // pgr is allocated but with wrong size -> free
			free(m_amorph_pgr);
			m_amorph_pgr = NULL;
		}
		if (NULL == m_amorph_pgr) { // pgr is not allocated -> allocate for square size ndim
			m_amorph_pgr = (fcmplx*)malloc(sizeof(fcmplx)*ndim*ndim);
			memset(m_amorph_pgr, 0, sizeof(fcmplx)*ndim*ndim);
			m_amorph_dim = ndim;
		}
		if (NULL == m_amorph_pgr) { // allocation failed
			return 2;
		}
		// calculate on new array
		nitems = ndim * ndim;
		m_amorph_samp = s;
		ndim2 = (ndim - (ndim % 2)) >> 1; // nyquist pixel index
		float fd = (float)(1. / pow(_JPG_AMORPH_SCATANG1,2));
		float dwf = (float)(0.25*_JPG_AMORPH_DWPRM);
		float sq = 1.f / (s * ndim); // Fourier scale for atomic form factors
		float qy2 = 0.f, q2 = 0.f;
		// 1) setup atomic form factor in Fourier space (float* atff)
		float* qn = (float*)malloc(sizeof(float)*ndim);
		if (NULL == qn) return 2;
		for (i = 0; i < ndim; i++) { // setup frequency hash: i -> q (not scrambled!!!)
			qn[i] = sq * (((i + ndim2) % ndim) - ndim2);
		}
		float* atff = (float*)malloc(sizeof(float)*nitems);
		if (NULL == atff) return 2;
		for (j = 0; j < ndim; j++) {
			ij = j * ndim;
			qy2 = qn[j] * qn[j];
			for (i = 0; i < ndim; i++) {
				q2 = qy2 + qn[i] * qn[i];
				atff[i + ij] = (float)_JPG_AMORPH_PHASEMAX / (q2*fd + 1.0f)*exp(-dwf * q2);
			}
		}
		// 2) setup random phases in real-space (float* rndpha)
		float* rndpha = (float*)malloc(sizeof(float)*nitems);
		if (NULL == rndpha) return 2;
		for (ls = 0; ls < nitems; ls++) {
			rand_s(&urnd); // get random number
			rndpha[ls] = (float)(_TPI * (double)urnd / (double)(UINT_MAX-1) ); // set random phase values from 0 to 2*Pi
		}
		// 3) transform random phase to Fourier-space
		CJFFTWcore jfft;
		int nd2[2] = { ndim,ndim };
		jfft.Init(2, nd2);
		jfft.Zero();
		jfft.SetDataRe(rndpha);
		jfft.FT(); // Fourier Transform
		jfft.Scale(1.f / (float)ndim); // rescale
		// 4) multiply by atff -> (fcmplx* obj)
		jfft.MultiplyReal(atff);
		// 5) inverse Fourier transform
		jfft.IFT(); // Inverse Fourier Transform
		jfft.Scale(1.f / (float)ndim); // rescale
		// 6) m_amorph_pgr = exp( i * obj.re)
		fcmplx obj;
		jfft.GetDataC(m_amorph_pgr);
		for (ls = 0; ls < nitems; ls++) {
			obj = m_amorph_pgr[ls];
			m_amorph_pgr[ls] = fcmplx( cos(obj.real()), sin(obj.real()) ); // exp(i*obj.re) = cos(obj.re) + i*sin(obj.re)
		}
		// done.
		free(qn);
		free(atff);
		free(rndpha);
		//free(obj);
		jfft.Deinit();
	}
	return nerr;
}


int CJProbeGen::CalculateRonchigram(CJProbeParams* prm, int ndim, float s, float* pdata)
{
	// assuming that all parameters are valid. Check before calling!
	int nerr = 0; // error code
	if (NULL == pdata) { // always check the validity of the pointer
		return 1;
	}
	int nft[2] = { ndim,ndim };
	CJFFTWcore jfft;
	jfft.Init(2, nft, FFTW_MEASURE);
	size_t nbytes = sizeof(float) * ndim; // number of bytes per row and column
	size_t nsbytes = sizeof(float) * ndim * ndim; // number of bytes in output
	memset(pdata, 0, nsbytes);
	CJProbeParams lprm = *prm; // get a copy of the probe parameters
	int nab = lprm.GetAberrationNum(); // number of aberrations
	int ndim2 = (ndim - (ndim % 2)) >> 1; // nyquist pixel index
	int i = 0, j = 0, ij = 0, ij3 = 0, k = 0, i1 = 0, j1 = 0; // iterators
	float a = s * ndim;  // physical box size(x and y)
	float sq = 1.f / a; // q sampling = 1 / physical box size (x and y)
	int *usc = (int*)malloc(sizeof(int)*ndim); // allocate unscrambling hash
	memset(usc, 0, sizeof(int)*ndim);
	float *qn = (float*)malloc(nbytes); // allocate frequency helper
	memset(qn, 0, nbytes);
	fcmplx *wav = (fcmplx*)malloc(nsbytes * 2); // allocate buffer holding the wave function (2 floats per pixel)
	memset(wav, 0, nsbytes * 2);
	float *tmp3 = (float*)malloc(nsbytes * 3); // allocate buffer holding the aberration function, focal term and aperture (3 floats per pixel)
	memset(tmp3, 0, nsbytes * 3);
	float *probe = (float*)malloc(nsbytes); // allocate buffer for intermediate probe inetnsities
	memset(probe, 0, nsbytes);
	float fsw = lprm.m_fspread_width;
	if (fabs(fsw) < 0.01f) { // extremely small or invalid focus spread
		lprm.m_fspread_kernel_samples = 1; // turn the focus spread convolution off
	}
	int fsnum = max(1, lprm.m_fspread_kernel_samples);
	float fsrw = fabs(lprm.m_fspread_kernel_width); // relative width of the focus spread kernel
	float btx = lprm.m_btx * 0.001f / lprm.m_wl; // beam tilt X: mrad -> 1/nm
	float bty = lprm.m_bty * 0.001f / lprm.m_wl; // beam tilt Y: mrad -> 1/nm
	float qlim = lprm.m_alpha * 0.001f / lprm.m_wl; // aperture size: mrad -> 1/nm
	float qlcx = lprm.m_alpha_x0 * 0.001f / lprm.m_wl; // aperture shift X: mrad -> 1/nm
	float qlcy = lprm.m_alpha_y0 * 0.001f / lprm.m_wl; // aperture shift Y: mrad -> 1/nm
	float fswgt = 1.f; // focal kernel weight
	float dzval = 0.f; // addional defocus in focal averaging
	float chi = 0.f; // aberration function value
	float wap = 0.f; // aperture value
	float zchi = 0.f; // focal phase term
	float ptot = 0.f; // total probe intensity for normalization
	float psum = 0.f; // sum of probe intensities
	float fsca = 1.f; // normalization
	float qx = 0.f, qy = 0.f; // coordinates
	float *ptmp3 = NULL;
	for (i = 0; i < ndim; i++) { // setup frequency hash: i -> q and unscramble hash
		usc[i] = (i + ndim2) % ndim;
		qn[i] = sq * (usc[i] - ndim2);
	}
	// set amorphous sample data (uses a module backup if present)
	nerr = SetAmorph(ndim, s);
	// set basis aberration function and the extra focal term in chi2
	for (j = 0; j < ndim; j++) { // loop over grid rows
		qy = qn[j] + bty;
		ij3 = j * ndim * 3;
		for (i = 0; i < ndim; i++) { // loop over grid columns
			ptmp3 = &tmp3[3 * i + ij3];
			qx = qn[i] + btx;
			wap = ApertureFunctionS(qx, qy, qlcx, qlcy, qlim, a, a); // aperture with smoothing of 1 pixel
			chi = 0.f;
			zchi = 0.f;
			if (wap >(float)_JPG_PROBE_APERTURE_THRESH) {
				if (nab > 0) {
					chi = AberrationFunction(qx, qy, lprm.m_wl, nab, lprm.m_abrr_coeff); // aberration function
				}
				zchi = (float)_PI*lprm.m_wl*(qx * qx + qy * qy); // Pi * Lambda * g^2
			}
			ptmp3[0] = wap; // store aperture
			ptmp3[1] = chi; // store aberration function
			ptmp3[2] = zchi; // store focal term
		}
	}
	// explicit focal averaging
	for (k = 0; k < fsnum; k++) {
		if (fsnum > 1) {
			dzval = 2.f * fsw * fsrw * ((float)k / (fsnum - 1) - 0.5f);
			fswgt = exp(-dzval * dzval / (fsw*fsw));
		}
		else {
			dzval = 0.f;
			fswgt = 1.f;
		}
		// calculate the coherent probe intensity for this defocus (dzval)
		for (j = 0; j < ndim; j++) { // loop over grid rows
			ij = j * ndim;
			ij3 = ij * 3;
			for (i = 0; i < ndim; i++) { // loop over grid columns
				ptmp3 = &tmp3[3 * i + ij3];
				wap = ptmp3[0]; // get aperture
				if (wap >(float)_JPG_PROBE_APERTURE_THRESH) {
					chi = ptmp3[1] + dzval * ptmp3[2]; // get aberration function with extra focal term
					wav[i + ij] = fcmplx(wap * cosf(chi), -wap * sinf(chi)); // set probe wave function coefficient
				}
			}
		}
		// inverse fourier transform
		jfft.SetDataC(wav); // copy data to fft module
		jfft.IFT(); // make an inverse fft
		jfft.MultiplyC(m_amorph_pgr); // multiply the transmission function
		jfft.FT(); // transform to Fourier space again
		ptot = jfft.GetDataTotalPow();
		fsca = fswgt / ptot;
		psum += fswgt; // accumulate focal weights
		jfft.GetDataPow(probe); // copy the result to the intermediate buffer
		for (j = 0; j < ndim; j++) { // loop over rows
			j1 = usc[j]; // index in the unscrambled array
			ij = j * ndim;
			ij3 = j1 * ndim;
			for (i = 0; i < ndim; i++) { // loop over columns
				i1 = usc[i]; // index in the unscrambled array
				pdata[i1 + ij3] += (fsca*probe[i + ij]); // add to output, unscrambled
			}
		}
	}
	// Normalize the probe intensity
	ScaleIntensity(1.f / psum, ndim*ndim, pdata); // with this call pdata contains the final result
	jfft.Deinit();
	free(probe); // free memory of the intermediate probe buffer
	free(tmp3); // free memory of the temp buffer
	free(wav); // free memory of wave helper
	free(qn); // free memory of q hash
	free(usc); // free memory of unscramble hash
	return nerr;
}


int CJProbeGen::CalculateProbeFunction(int idx, CJProbeParams* prm, int ndim, float s, float* pdata)
{
	int nerr = 0;
	// requires sufficient memory allocated to receive sizeof(float)*ndim*ndim bytes in pdata
	if (NULL == pdata || idx < 0 || idx >= _JPG_PROBEFUNC_NUM || ndim < 0 || ndim > 8192 || s <= 0.f ||
		prm->m_alpha <= 0.f || prm->m_wl <= 0.f) {
		return 1; // invalid setup
	}
	// the code forks into sub-functions
	m_clock.now();
	auto cl_start = m_clock.now();
	switch (idx)
	{
	case 0: // Coherent Probe Intensity
		nerr = CalculateProbeIntensityCoh(prm, ndim, s, pdata);
		break;
	case 1: // Coherent Probe Phase
		nerr = CalculateProbePhase(prm, ndim, s, pdata);
		break;
	case 2: // Partially Spatial Coherent Probe Intensity 
		nerr = CalculateProbeIntensityPSC(prm, ndim, s, pdata);
		break;
	case 3: // Partially Temporal Coherent Probe Intensity
		nerr = CalculateProbeIntensityPTC(prm, ndim, s, pdata);
		break;
	case 4: // Partially Coherent Probe Intensity
		nerr = CalculateProbeIntensity(prm, ndim, s, pdata);
		break;
	case 5: // Real Part of Probe Wave Function
		nerr = CalculateProbeRe(prm, ndim, s, pdata);
		break;
	case 6: // Imaginary Part of Probe Wave Function
		nerr = CalculateProbeIm(prm, ndim, s, pdata);
		break;
	case 7: // Aberration Phase Plate
		nerr = CalculateAberrationPhasePlate(prm, ndim, s, pdata);
		break;
	case 8: // Ronchigram
		nerr = CalculateRonchigram(prm, ndim, s, pdata);
		break;
	default: // Partially Coherent Probe Intensity
		nerr = CalculateProbeIntensity(prm, ndim, s, pdata);
		break;
	}
	auto cl_stop = m_clock.now();
	auto cl_delta = cl_stop - cl_start;
	m_dLastProbeCalculationTime = (double)cl_delta.count() * 1.e-9;
	return nerr;
}

bool CJProbeGen::IsProbeFunctionDiffraction(int idx)
{
	if (idx == 7 || idx == 8) { // any implemented diffraction function
		return true;
	}
	return false;
}