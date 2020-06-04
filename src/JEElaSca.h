//
// C++ header file: JEElaSca.h
// class declaration for library JMultislice.lib (implementation see JEElaSca.cpp)
//
//
// Copyright (C) 2020 - Juri Barthel (juribarthel@gmail.com)
// Copyright (C) 2020 - Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
//
// Verions of JEElaSca: 1.0 (2020 - June - 02)
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
// JEElaSca implements routines for generating object transmission functions
// for the JMultiSliceLib project using GPU and CPU code. The transmission
// functions are for elastic scattering, quantum excitation of phonons and
// absorptive potentials.
//
//
// Memory management and setup:
//
// The class handles helper arrays for the calculation by itself. This requires
// a strict sequence of setup calls to be done:
// 1) ResetAll() : resets all variables into a state to receive new data
// 2) SetStructure(CStructure st) : copy a list of atom types to prepare form factors from that
// 2) SetGrid(size_t nx, size_t ny, size_t nz, unsigned int pot3d) : sets samplings
// 3) SetPhysics()
// 4) AllocHelpers()
// 5) SetHelpers()
// 6) Once these setup steps are done, use one of the following functions to
//    calculate phase gratings or potentials
//    a) CalculatePhasegratingCPU(CStructure st, fcmplx* pgr, int ithread)
//       from a projected potential of the given structure with CPU thread ithread.
//    b) CalculatePhasegratingGPU(CStructure st, fcmplx* pgr)
//       from a projected potential of the given structure with the currently set GPU.
//    The structure parameters st must work with the syme atom type list as that used
//    with the SetStructure call in step (2) !
//
#pragma once

#include "NatureConstants.h"
#include "Structure.h"
#include "rng.h"
//#include "JFFTWcore.h"
#include "JFFTMKLcore.h"
#include "JFFTCUDAcore.h"

#ifndef __EELSCA_MACROS
#define __EELSCA_MACROS

#define EELSCA_STATUS_NOTREADY		0
#define EELSCA_STATUS_STRUCTURE		1
#define EELSCA_STATUS_PHYSICS		2
#define EELSCA_STATUS_GRID			4
#define EELSCA_STATUS_ATOMICFF		8
#define EELSCA_STATUS_HELPER_ALLOC	16
#define EELSCA_STATUS_HELPER_SET	32

#define EELSCA_THRESHOLD_ATOMICFF	7
#define EELSCA_THRESHOLD_HELPERS	21
#define EELSCA_THRESHOLD_POTENTIALS	63

#endif // !__EELSCA_MACROS

#ifndef __EELSCA_CONST
#define __EELSCA_CONST

// inverse Yukawa range to handle divergence in ionic charge potentials. (This is a fudge.)
constexpr float EELSCA_CONST_IYR = 0.2f; // [1 / nm]

// 1.0 / sqrt( 8 * Pi^2 ) rescales from sqrt( Biso ) to <ux>
constexpr float EELSCA_CONST_RSQR8PI2 = 0.1125395f;

// projected potential scale: 0.03809982080 ! hbar ^ 2 / (2 * m0 * e) * (10 ^ 9) ^ 2 )[V nm^2]
constexpr float EELSCA_CONST_POT_SCA = 0.03809982f;

// relative band-width limit of phase gratings (recommendation: set equal to _JMS_RELAPERTURE, see "JMultislice.h")
constexpr float EELSCA_CONST_PGR_BWL  = 0.6666667f;

// interaction constant: 2.0886573708 ! m0 * e / (2*Pi * hbar^2) * (10^-9)^2  [ V^-1 nm^-2 ]
constexpr float EELSCA_CONST_PGR_SIG = 2.08865737f;

#endif // !__EELSCA_CONST


class CJEElaSca {

	// ----------------------------------------------------------------------------
	// ----------------------------------------------------------------------------
	// constructors & desctructor

public:

	// default constructor
	CJEElaSca(void);
	// destructor
	~CJEElaSca(void);

	// ----------------------------------------------------------------------------
	// ----------------------------------------------------------------------------
	// member variables

protected:

	// infrastructure

	unsigned int m_status; // tracks the status of the parameter setup

	// CPU core objects for FFTs and array ops on host using FFTW
	//CJFFTWcore* m_jcpuco;

	// CPU core objects for FFTs and array ops on host using the Intel MKL
	CJFFTMKLcore* m_jcpuco;

	// GPU core object for FFTs and array ops on device using the cufft library
	CJFFTCUDAcore m_jgpuco;

	// random number generator
	CRng m_lrng;
	CRng* m_prng;
	
	// structure data

	CV3D m_vCellSize; // (a,b,c) size of the structure super-cell [nm]
	std::vector<CAtomType> m_vTypes; // list of atom types

	// sampling (these variables control the class memory consumption)

	unsigned int m_grid_x; // grid size along x
	unsigned int m_grid_y; // grid size along y
	unsigned int m_grid_z; // grid size along z
	unsigned int m_pot3d; // swiches the calculation of 3D potentials
	
	// physics
	bool m_bwl_pot; // flag apply band width limit to potentials (ensures round cut-off of form factors close to spatial frequency limit determined by the chosen x,y sampling)
	bool m_bwl_pgr; // flag apply band width limit to phase gratings (ensures cut-off of alias frequencies beyond 2/3 of Nyquist due to the application of repeated scattering an propagation with Fourier transforms)
	unsigned int mod_atffacs; // model of atomic form factors (0: Weickenmeier & Kohl)
	unsigned int mod_thermal; // model of thermal vibration effects (0: none, 1: Debye-Waller factors, 2: QEP)
	unsigned int mod_absorb; // model of absorption effects (0: none, 1: Hashimoto et al., 2: Hall & Hirsch)
	float m_ekv; // electron energy in keV
	float m_vabf; // fix absorption factor (mod_absorb == 1)

	/* --------------- object calculation helper variables --------------- */
	/*                                                                     */
	/* These helpers are used by the objects calculation functions. The    */
	/* object handles the states of these variable by itself.              */
	/*                                                                     */
	
	size_t m_att_num; // number of atom types (used for allocation of helpers)
	int m_cpu_num; // number of cpu threads (used for allocation of helpers)
	int m_gpu_id; // id of the gpu
	float* m_h_qx; // x shifting helper on host
	float* m_h_qy; // y shifting helper on host
	float* m_h_qz; // z shifting helper on host
	float* m_d_qx; // x shifting helper on device
	float* m_d_qy; // y shifting helper on device
	float* m_d_qz; // z shifting helper on device
	float* m_h_q2; // q^2 2d helper array on host
	float* m_h_ap; // two 2d helper arrays with apertures for 1: potentials and 2: phase gratings on host
	float* m_d_ap; // two 2d helper arrays with apertures for 1: potentials and 2: phase gratings on device
	fcmplx* m_h_potio; // ionic point-charge potential term on host
	fcmplx* m_d_potio; // ionic point-charge potential term on device
	fcmplx* m_h_att_ff2d; // 2d projected form factors of all atom types on host (prepared to be copied by each thread or gpu call) 
	fcmplx* m_d_ff2d; // 2d projected atomic form factor on device

	/*                                                                     */
	/* ------------------------------------------------------------------- */


	// ----------------------------------------------------------------------------
	// ----------------------------------------------------------------------------
	// member functions

public:

	// Set the current random number generator
	// Calling with NULL pointer will reset m_prng to the local CRng object m_lrng
	// Use this, if you want to control the random number generation
	void SetRng(CRng* prng = NULL);

	// Returns a pointer to the current random number generator object
	CRng* GetRng(void);

	// Seed the current random number generator m_prng
	// Seeding with 0 will cause a seed depending on current system time
	void SeedRngEx(int nseed = 0);

	// Resets all information
	void ResetAll(void);

	// Set the super cell size and the atom type list from the given structure data
	int SetStructure(CStructure st);

	// Set the grid dimensions
	int SetGrid(unsigned int nx, unsigned int ny, unsigned int nz, unsigned int pot3d);

	// Set the physics of beam and scattering simulation
	// ekv : primary electron energy in keV
	// attfacs : atomic form factor table
	//           0 : Weickenmeier and Kohl, Acta Cryst. A 47 (1991) 590-597.
	// absorp : absorptive form factor model
	//          0 : no absorption
	//          1 : Hashimoto et al (fraction of the elastic form factor)
	//          2 : Hall & Hirsch (loss in the elastic channel due to thermal diffuse scattering)
	// vabf : absorption factor used when absorption model is used with absorb = 1
	// thermal : thermal diffuse scattering (TDS)
	//          0 : not simulated
	//          1 : Debye-Waller factors (TDS removed)
	//          2 : Quantum Excitation of Phonons model (alternatively, frozen lattice, frozen phonons)
	int SetPhysics(float ekv, unsigned int attfacs, unsigned int absorb, float vabf, unsigned int thermal);

	// Return true if potentials are dampened by Debye-Waller factors.
	bool UseDWF(void);

	// Returns the Debye-Waller factor for a given biso value (nm^2) and a squared spatial fequency ksqr (1/nm^2)
	float GetDWF(float biso, float ksqr);

	// Returns the relativistic correction for a given kinetic electron energy (keV)
	float GetRelCorr(float ekv);

	// Returns true if atomic potentials are displaced to simulate thermal vibrations in the QEP model
	bool UseQEP(void);

	// Runs the calculation of radial form factor tables for all atom types in m_vTypes structure data.
	// The form factors are stored in the object member items of m_vTypes and thus work only for
	// structure data (lists of CAtom) which refer to the same sequence of atom types.
	int CalculateFormFactors(bool bverbose = false);

	// Prepares scattering functions for using them with later calls of CalculatePhaseGrating or CalculatePotential.
	// The prepared data will be stored in the object helper arrays and remain on host and device memory until the
	// object is destructed or until you call ResetAll().
	int PrepareScatteringFunctions(int igpu, int ncpu, bool bverbose = false);

	// Calculates the mean inner potential of a structure
	int CalculateMeanInnerPotential(CStructure* pst, fcmplx& mip);

	// Calculate a phase grating from the projected potential of the structure st output to the buffer pgr
	// using the CPU thread ID ithread.
	int CalculatePhasegratingCPU(CStructure* pst, fcmplx* pgr, int ithread);

	// Calculate a phase grating from the projected potential of the structure st output to the buffer pgr
	// using the current GPU.
	int CalculatePhasegratingGPU(CStructure* pst, fcmplx* pgr);

protected:

	// Allocates memory for calculation helper arrays
	int AllocHelper(int igpu, int ncpu);

	// Deallocates memory of helper arrays
	void DeallocHelper(void);

	// Initializes helper array values
	int InitHelper(void);
	
	// returns uniform random variates between 0 and 1
	float UniRand(void);

	// returns normal random variates with mean 0 and standard deviation 1
	float NormRand(void);
};