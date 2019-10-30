//
// C++ header file: JPlasmonMC.h
// class declaration for library JMultislice.lib (implementation see JPlasmonMC.cpp)
//
//
// Copyright (C) 2019 - Juri Barthel (juribarthel@gmail.com)
// Copyright (C) 2019 - Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
//
// Verions of JPlasmonMC: 1.0.0 (2019 - Oct - 15)
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
// JPlasmonMC implements routines for plasmon scattering Monte-Carlo
// used by the JMultiSliceLib project using CPU and GPU code.
//
// 1) CJPlasmonMC : implements plasmon scattering data and functions.
//
//
// Assumes that the random number generator of the thread has been
// initialized by the calling code.
//
#pragma once
//
//#include "JFFTWcore.h"
#include "JFFTMKLcore.h"
#include "JFFTCUDAcore.h"
using namespace std;
//
#ifndef __JPLASMC__
#define __JPLASMC__
#define _JPL_MESSAGE_LEN				2048 // default message length
#endif
//
//
// Declaration of class class CJPlasmonMC
//
class CJPlasmonMC {

	// ----------------------------------------------------------------------------
	// ----------------------------------------------------------------------------
	// constructors & desctructor

public:

	// default constructor
	CJPlasmonMC(void);
	// copy constructore
	CJPlasmonMC(const CJPlasmonMC & src);
	// destructor
	~CJPlasmonMC(void);

	// ----------------------------------------------------------------------------
	// ----------------------------------------------------------------------------
	// data member variables

	// number of allowed plasmon excitations per probing electron in the sample
	// - set to zero to deactive the Monte-Carlo
	int m_sca_num_max;
	// characteristic scattering vector [1/nm]
	float m_q_e;
	float m_q_e2; // squared
	// critical scattering vector [1/nm]
	float m_q_c;
	float m_q_c2; // squared
	// mean-free path for single scattering [nm]
	float m_meanfreepath;

	// number of plasmon excitations in current run
	int m_sca_num;
	// total accumulated scattering vector 1/nm, x-component
	float m_sca_tot_qx;
	// total accumulated scattering vector 1/nm, y-component
	float m_sca_tot_qy;
	
	// processing member variables
	// These members are not shared in copy constructors or operators
protected:
	// default message string for output (not initialized)
	char m_msg[_JPL_MESSAGE_LEN];
	

public:
	// Operators

	// copy data
	void operator = (const CJPlasmonMC &src);
	// compare data
	bool operator == (const CJPlasmonMC &src) const;
	
	// Member functions
	// initialized the MC parameters
	// - requires members m_meanfreepath, m_q_e, and m_q_c to be set before calling
	void Init(void);
	// reset scattering run
	void ResetMC(void);
	// tests for scattering and returns current total scattering angle,
	// doesn't modify any of the class members
	// - dt: slice thickness in [nm], input
	// - num_sca_max: max. number of excitations allowed in the slice, input
	// - num_sca: number of excitation occurring in the slice, output
	// - sca_qx: scattering vector, x-component, output
	// - sca_qy: scattering vector, y-component, output
	void ScatteringMC(float dt, UINT num_sca_max, UINT & num_sca, float & sca_qx, float & sca_qy);

	// tests for scattering and returns pixel shifts to be applied,
	// returns number of excitations occurring in the slice,
	// modifies class members
	// - dt: slice thickness in [nm], input
	// - a_x: horizontal grid size [nm], input
	// - a_y: vertical grid size [nm], input
	// - di: horizontal pixel shift, output
	// - dj: vertical pixel shift, output
	int ScatGridMC(float dt, float a_x, float a_y, int *di, int *dj);

protected:
	// returns uniform random variates between 0 and 1
	float UniRand(void);
	// returns Poissonian random variates of mean value m
	UINT PoissonRand(float m);
};
//
