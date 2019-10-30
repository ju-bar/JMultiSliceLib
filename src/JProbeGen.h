//
// C++ header file: JProbeGen.h
// class declaration for library JMultislice.lib (implementation see JProbeGen.cpp)
//
//
// Copyright (C) 2018 - Juri Barthel (juribarthel@gmail.com)
// Copyright (C) 2018 - Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
//
// Verions of JProbeGen: 0.12 (2019 - Sep - 16)
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
//
// JProbeGen implements routines for generating probe electron wave functions
// for the JMultiSliceLib project using CPU code.
//
// The implementation is with two classes:
// 1) CJProbeParams : implements a data class defining physical and
//    numerical parameters for calculations with the class CJProbeGen
// 2) CJProbeGen : implements probe calculation routines.
//
// Memory management & threading:
//
// The code includes a data class CJProbeParams which is used as
// function interface on calling the calculation routines of CJProbeGen.
// By this way no data is stored in objects of CJProbeGen.
//
//
#pragma once
//
#include <vector>
#include <chrono>
//#include "JFFTWcore.h"
#include "JFFTMKLcore.h"
using namespace std;
//
#ifndef __JPROBEGEN__
#define __JPROBEGEN__
#define _JPG_MESSAGE_LEN				2048 // default message length
#define _JPG_PROBE_APERTURE_THRESH		(1.E-6) // aperture strength threshold
#define _JPG_ABERRATION_ORDER_MAX		8 // maximum order of axial coherent aberrations
#define _JPG_ABERRATION_STRTHRESH		(1.E-15) // aberration strength threshold [nm]
#define _JPG_PROBEFUNC_NUM				9 // number of supported probe functions
#define _JPG_SOURCE_TYPE_MAX			3 // max index of source types
#define _JPG_AMORPH_PHASEMAX			(0.4) // max phase of the amorphous sample
#define _JPG_AMORPH_SCATANG1			(3.0) // scattering angle parameter 1 [1/nm]
#define _JPG_AMORPH_SCATANG2			(5.0) // scattering angle parameter 2 [1/nm]
#define _JPG_AMORPH_SCATTAIL			(0.2) // relative power of large angle tails
#define _JPG_AMORPH_DWPRM				(0.08) // vibration amplitude parameter for Debye - Waller factors
#endif
//
//
// Declaration of class class CJProbeParams
//
class CJProbeParams {

	// ----------------------------------------------------------------------------
	// ----------------------------------------------------------------------------
	// constructors & desctructor

public:

	// default constructor
	CJProbeParams(void);
	// copy constructore
	CJProbeParams(const CJProbeParams & src); 
	// destructor
	~CJProbeParams(void);

	// ----------------------------------------------------------------------------
	// ----------------------------------------------------------------------------
	// data member variables (all public for easy access)

	// electron wavelength [nm]
	float m_wl;
	// semi-angle of convergence [mrad], illumination aperture radius
	float m_alpha;
	// illumination aperture relative smoothness
	float m_alpha_rs;
	// illumination aperture center [mrad]
	float m_alpha_x0;
	float m_alpha_y0;
	// illumination aperture relative asymmetry
	float m_alpha_asym;
	// illumination aperture asymmetry direction [deg]
	float m_alpha_adir;
	// beam tilt [mrad]
	float m_btx;
	float m_bty;
	// effective source size (HWHM) [nm]
	float m_source_width;
	// effective source shape (index 0: point, 1: Gaussian, 2: Cauchy, 3: disk)
	int m_source_shape;
	// effective focus spread (1/e-half-width) [nm]
	float m_fspread_width;
	// relative focus spread kernel width (w.r.t. m_source_width)
	float m_fspread_kernel_width;
	// number of focus spread kernel samples
	int m_fspread_kernel_samples;
	// sorted list of aberration corefficients
	float* m_abrr_coeff;

protected:
	// processing member variables (all protected)
	// These members are not shared in copy constructors or opertors
	
	// default message string for output (not initialized)
	char m_msg[_JPG_MESSAGE_LEN];
	// max. length of accepted aberration table (number of aberrations)
	int m_abrr_num;

public:
	// Operators

	// copy data
	void operator = (const CJProbeParams &other);
	// compare data
	bool operator == (const CJProbeParams &other) const;
	
	// Member functions

	// returns number of supported aberrations
	int GetAberrationNum(void);

};
//
//
// Declaration of class class CJProbeGen
//
class CJProbeGen {

	// ----------------------------------------------------------------------------
	// ----------------------------------------------------------------------------
	// constructor

public:
	CJProbeGen();
	~CJProbeGen();

	// ----------------------------------------------------------------------------
	// ----------------------------------------------------------------------------
	// member variables

public:
	// information data (public)

	// last probe calculation time in seconds
	double m_dLastProbeCalculationTime;

protected:
	// processing member variables (all protected)
	// These members are not shared in copy constructors or operators

	// default message string for output (not initialized)
	char m_msg[_JPG_MESSAGE_LEN];
	
	// length of the aberration table (number of aberrations)
	int m_abrr_num;
	
	// list of aberration indices l, n
	int* m_abrr_widx;
	
	// list of aberration function binomial factors
	int* m_abrr_binom;
	
	// list of aberration names
	vector<string> m_abrr_name;
	
	// list of aberration symbols
	vector<string> m_abrr_symbol;
	
	// list of function name
	vector<string> m_pfunc_name;
	
	// dimension of amorphous backup data for ronchigram simulation
	int m_amorph_dim;
	
	// sampling rate of the amorphous backup data
	float m_amorph_samp;
	
	// amorphous backup data for ronchigram simulation
	fcmplx* m_amorph_pgr;
	
	// clock object with nanoseconds resolution (does this work on Linux or Mac?)
	chrono::high_resolution_clock m_clock;

public:

	// ----------------------------------------------------------------------------
	// ----------------------------------------------------------------------------
	// member functions

	// cleans memory used to calculate FFTs (use with utmost care)
	void FreeLibMem(void);

	// returns the sum of the first n values in pdata
	float GetTotalIntensity(size_t n, float* pdata);

	// multiplies the first n values of pdata by fsca
	void ScaleIntensity(float fsca, size_t n, float *pdata);

	// calculates the value of a round aperture function in the diffraction plane
	// - qx, qy: diffraction plane coordinate [1/nm]
	// - qcx, qcy: aperture center [1/nm]
	// - qlim: radius of the aperture [1/nm]
	float ApertureFunction(float qx, float qy, float qcx, float qcy, float qlim);
	
	// calculates the value of a round aperture function in the diffraction plane
	// The edge of the aperture is smoothed by about one pixel
	// - qx, qy: diffraction plane coordinate [1/nm]
	// - qcx, qcy: aperture center [1/nm]
	// - qlim: radius of the aperture [1/nm]
	// - smooth: relative smoothness of the edge
	float ApertureFunctionS(float qx, float qy, float qcx, float qcy, float qlim, float smooth = 0.05f);

	// calculates the value of an elliptical aperture function in the diffraction plane
	// The edge of the aperture is smoothed by about one pixel
	// - qx, qy: diffraction plane coordinate [1/nm]
	// - qcx, qcy: aperture center [1/nm]
	// - qlim: radius of the aperture [1/nm]
	// - alim: relative asymmetry of the aperture
	// - adir: direction of asymettry [rad]
	// - smooth: relative smoothness of the edge
	float ApertureFunctionA(float qx, float qy, float qcx, float qcy, float qlim, float alim, float adir, float smooth = 0.05f);
	
	// calculates the value of an aberration function for a given diffraction plane coordinate {qx,qy}
	// - qx, qy: diffraction plane coordinate [1/nm]
	// - wl: electron wavelength [nm]
	// - nabrr: highest index of the aberration list (length-1)
	// - aberrcoeff: ordered list of nabrr aberrations with 2 coefficients each
	//   The order is {a11x,a11y,a20x,a20y,a22x,a22y,a31x,a31y,a33x,a33y,a40x,a40y, ...)
	//   amnx,y are the coefficients [nm]
	float AberrationFunction(float qx, float qy, float wl, int nabrr, float* aberrcoeff);

	// calculates a map of the normalized source distribution with periodic wrap-around
	// The distribution is centered at pixel {0,0} and spans in positive directions
	// over {0,0}->{nx/2-1,ny/2-1} and negative directions {nx/2,ny/2}->{nx,ny}
	// - src_type: source type index (0: delta, 1: Gaussian, 2: Cauchy, 3: disk)
	// - src_width: source width [nm] (HWHM)
	// - nx, ny: number of samples along x and y
	// - ax, ay: size of the grid along x and y [nm]
	// - krn: pointer to a pre-allocated array receiving the result
	int GetSourceDistribution(int src_type, float src_width, int nx, int ny, float ax, float ay, float *krn);
	
	// returns the maximum aberration order supported by the module
	int GetMaxAberrationOrder(void);
	
	// returns the number of aberrations supported by the module
	int GetAberrationNum(void);
	
	// returns the name of the aberration of index idx in abrr_name
	int GetAberrationName(int idx, std::string * abrr_name);
	
	// returns the symbol of the aberration of index idx in abrr_symbol
	int GetAberrationSymbol(int idx, std::string * abrr_symbol);
	
	// returns the number of probe functions (max. index + 1) supported by the module
	int GetMaxProbeFunction(void);

	// calculates a STEM probe wavefunction in Fourier space using current parameters
	// - prm: address of a CJProbeParams object defining physical probe parameters
	// - nx, ny: number of pixels fo wav
	// - ax, ay: physical size of wav [nm]
	// - wav: address receiving the probe wave function
	int CalculateProbeWaveFourier(CJProbeParams* prm, int nx, int ny, float ax, float ay, fcmplx *wav);

	// calculates an aberration phase plate
	// - prm: aberration parameters
	// - ndim: size of the square phase plate
	// - s: physical sampling rate in real space [nm/pix]
	// - pdata: pointer to float array receiving data
	int CalculateAberrationPhasePlate(CJProbeParams* prm, int ndim, float s, float* pdata);

	// calculates the intensity distribution of a coherent STEM probe
	// - prm: probe parameters
	// - ndim: size of the square phase plate
	// - s: physical sampling rate in real space [nm/pix]
	// - pdata: pointer to float array receiving data
	int CalculateProbeIntensityCoh(CJProbeParams* prm, int ndim, float s, float* pdata);

	// calculates the intensity distribution of a STEM probe with partial spatial coherence
	// - prm: probe parameters
	// - ndim: size of the square phase plate
	// - s: physical sampling rate in real space [nm/pix]
	// - pdata: pointer to float array receiving data
	int CalculateProbeIntensityPSC(CJProbeParams* prm, int ndim, float s, float* pdata);

	// calculates the intensity distribution of a STEM probe with partial temporal coherence
	// - prm: probe parameters
	// - ndim: size of the square phase plate
	// - s: physical sampling rate in real space [nm/pix]
	// - pdata: pointer to float array receiving data
	int CalculateProbeIntensityPTC(CJProbeParams* prm, int ndim, float s, float* pdata);

	// calculates the intensity distribution of a STEM probe with partial coherence
	// - prm: probe parameters
	// - ndim: size of the square phase plate
	// - s: physical sampling rate in real space [nm/pix]
	// - pdata: pointer to float array receiving data
	int CalculateProbeIntensity(CJProbeParams* prm, int ndim, float s, float* pdata);

	// calculates the phase map of a coherent STEM probe
	// - prm: probe parameters
	// - ndim: size of the square phase plate
	// - s: physical sampling rate in real space [nm/pix]
	// - pdata: pointer to float array receiving data
	int CalculateProbePhase(CJProbeParams* prm, int ndim, float s, float* pdata);

	// calculates the real part of a coherent STEM probe
	// - prm: probe parameters
	// - ndim: size of the square phase plate
	// - s: physical sampling rate in real space [nm/pix]
	// - pdata: pointer to float array receiving data
	int CalculateProbeRe(CJProbeParams* prm, int ndim, float s, float* pdata);

	// calculates the imaginary part of a coherent STEM probe
	// - prm: probe parameters
	// - ndim: size of the square phase plate
	// - s: physical sampling rate in real space [nm/pix]
	// - pdata: pointer to float array receiving data
	int CalculateProbeIm(CJProbeParams* prm, int ndim, float s, float* pdata);

	// calculates a phase grating for a thin amorphous samples on a square grid
	// of size ndim and sampling rate s
	int SetAmorph(int ndim, float s, bool bForce = false);

	// calculates the ronchigram of a thin amorphous object with a STEM probe
	// - prm: probe parameters
	// - ndim: size of the square phase plate
	// - s: physical sampling rate in real space [nm/pix]
	// - pdata: pointer to float array receiving data
	int CalculateRonchigram(CJProbeParams* prm, int ndim, float s, float* pdata);

	// returns the name of the probe function of index idx in func_name
	int GetProbeFunctionName(int idx, string * func_name);

	// calculates a probe function of type idx with square size
	// - idx: function index
	// - prm: probe parameters
	// - ndim: size of the square phase plate
	// - s: physical sampling rate in real space [nm/pix]
	// - pdata: pointer to float array receiving data
	int CalculateProbeFunction(int idx, CJProbeParams* prm, int ndim, float s, float* pdata);

	// returns whether the probe function displays in diffraction space
	// - idx: function index
	bool IsProbeFunctionDiffraction(int idx);

};
//
