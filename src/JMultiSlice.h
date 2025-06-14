//
// C++ header file: JMultiSlice.h
// declaration for library JMultislice.lib (implementation see JMultislice.cpp)
//
//
// Copyright (C) 2018 - 2020 - Juri Barthel (juribarthel@gmail.com)
// Copyright (C) 2018 - 2020 - Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
//
// Verions of JMultiSlice: 0.45b (2020 - June - 10)
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
// JMultiSlice implements routines of running the multislice algorithm on CPU and GPU.
// Uses libraries cufft.lib; cudart_static.lib;
//                mkl_intel_lp64.lib; mkl_sequential.lib; mkl_core.lib
//                (libfftw3f-3.lib is currently no longer used)
//
//
// Concerning Multithreading:
//
// The FFT-part of the class supports CPU multi-threading organized from
// external code (MKL). By initializing for a given number of threads,
// separate plans and arrays will be generated to allow independent
// calculation by different threads or processes.
// The CUDA-part of the class assumes that all parallel computing is done
// on the CUDA device.
//
// The header declares variables global in the library scope.
// The variables are parameters of the algorithm as well as pointers to large
// amount of data such as phase gratings, propagators, wave functions, and
// detector functions.
// Multislice data may be allocated on host or device memory depending on
// where the algorithm runs. Parallel usage of both, GPU and CPU is possible.
//
//
// Memory management:
//
// CPU code works on host arrays and outputs to host arrays. The input arrays
// can be quite large, so no copying is done here. The code accesses arrays
// prepared by external routines via lists of pointers (float**) or (fcmplx**).
// This concerns phase gratings, propagators, and detector functions. Make
// sure to allocate these arrays externally and keep them alive during the
// calculation. 
//
// Important: In order to preserve the wave function normalization (number
//            of particles), both, phase gratings and propagators, need to be
//            scaled by 1./sqrt( nx * ny ) !!!
//            Alternatively, one of the two gratings should be scaled by
//            1./(nx * ny), preferrably the propagators.
//            Still, image detector outputs are not properly scaled. This
//            should be done later by the calling code.
//
// GPU code works on device arrays and outputs to host as well as to device
// arrays. This code will allocate all arrays on the device and on the host.
//
// Interface routines are with public flag, internal routines are with
// protected or private flag.
// 
// The GetResult function allows to copy the results from the result arrays
// on host and device to be stored in host memory. All output arrays are
// allocated by the code itself to keep results separated for each thread.
//
// SETUP-CHAIN (minimum, without much physics)
// 1) JMS_SetGridSize, JMS_SetHighTension, JMS_SetSupercellSize
// 2) 2.1) JMS_PhaseGratingSetup,
//    2.2) JMS_ObjectSliceSetup,
//    2.3) JMS_PropagatorSetup,
//    2.4) JMS_DetectorSetup (to be called after JMS_ObjectSliceSetup)
// 3) 3.1) JMS_SetPhaseGratingData (for each slice, calculation possible via CJEElaSca)
//    3.2) JMS_SetSliceThickness (for each slice) (required for plasmon scattering)
// 4) JMS_SetPropagatorData (for each propagator)
// 5) JMS_SetDetectorData (for each detectors)
// 6) JMS_InitCore (once)
// 7) JMS_SetIncidentWave (once)
//
// For each multislice and as often as required
// 1) JMS_OffsetIncidentWave
// 2) JMS_CPUMultislice or JMS_GPUMultislice
// 3) JMS_GetResult
//
// DE-INITIALIZATION-CHAIN (minimum, without physics)
// 1) JMS_Cleanup
//
//
// MULTISLICE:
//
// (for each slice)
// - readout detectors
// - apply phase grating
// - propagation
// (for exit plane)
// - readout detectors
//
//
#pragma once
//
//#include "JFFTWcore.h"
#include "JFFTMKLcore.h"
#include "JFFTCUDAcore.h"
#include "JProbeGen.h"
#include "JPlasmonMC.h"
#include "rng.h"
//
#ifndef __JMS__
#define __JMS__
// VERSION NUMBERS
#define __JMS_VERSION__			0
#define __JMS_VERSION_SUB__		4
#define __JMS_VERSION_SUB_SUB__	5
#define __JMS_VERSION_BUILD__	20200610
// CODE IDs
#define _JMS_CODE_CPU			1
#define _JMS_CODE_GPU			2
// STATUS IDs
#define _JMS_STATUS_NONE		0
#define _JMS_STATUS_PGR			1
#define _JMS_STATUS_OBJ			2
#define _JMS_STATUS_PRO			4
#define _JMS_STATUS_DET			8
#define _JMS_STATUS_CORE		16
#define _JMS_THRESHOLD_CORE		15 // status threshold for core initialization
#define _JMS_THRESHOLD_CALC		31 // status threshold for calculation
// DETECTION AND ACCUMULATION FLAGS
#define _JMS_ACCMODE_NONE		0
#define _JMS_ACCMODE_INTEGRATE	1
#define _JMS_DETECT_NONE		0
#define _JMS_DETECT_INTEGRATED	1    // flag integrating detection
#define _JMS_DETECT_IMAGE		2    // flag image detection
#define _JMS_DETECT_DIFFRACTION	4    // flag diffraction detection
#define _JMS_DETECT_WAVEREAL	8    // flag real-space wave-function detection
#define _JMS_DETECT_WAVEFOURIER	16   // flag reciprocal-space wave-function detection
#define _JMS_DETECT_IMAGE_AVG	32   // flag average image detection
#define _JMS_DETECT_DIFFR_AVG	64   // flag average diffraction detection
#define _JMS_DETECT_WAVER_AVG	128  // flag average real-space wave-function detection
#define _JMS_DETECT_WAVEF_AVG	256  // flag average reciprocal-space wave-function detection
// OTHER PARAMETERS
#define _JMS_MESSAGE_LEN		2048 // max. length of message strings
#define _JMS_FFTW_PLANFLAG		FFTW_MEASURE // used FFTW planner flag, on typical sizes FFTW_MEASURE perfoms best, whereas FFTW_PATIENT requires much more initialization time
#define _JMS_SUMMATION_BUFFER	0x1000	// number of temporary buffer items for summation
#define _JMS_SUMMATION_BTF_THR	32	// min. number of items for a butterfly to work, otherwise fallback to straight sum
									// Using numbers larger than 100 may lead to increased computation time.

constexpr double _JMS_RELAPERTURE = 0.666666666666667; // relative band-width limit of propagators

//
// global functions

// straight float summation with float accumulator
void ffstrsum(float* a, size_t n, float *s);

// straight float summation with double accumulator
void fdstrsum(float* a, size_t n, float *s);

// kahan summation with float error correction
// see <https://en.wikipedia.org/wiki/Kahan_summation_algorithm>
void fkahan(float* a, size_t n, float *s);

// strided 2-fold butterfly summation
// ! Warning: This routine may modify the input a.
// ! Make a backup of the data or provide a copy if you still need it!
void fdncs2m(float* a, size_t n, float *s);
// strided 2-fold butterfly summation
void fdncs2(float* a, size_t n, float *s);

//
//
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Class CJMultiSlice makes the global scope for the calculations.
// It contains all shared data and functions to work on it.
// Declare one global (static) instance of the class and use this to run
// functions. Pass it to all threads!
class CJMultiSlice {

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// constructor

public:
	CJMultiSlice();
	~CJMultiSlice();

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// member variables

protected:

// ----------------------------------------------------------------------------
// library handling

	// debug level
	int m_dbg;
	// random number generator object
	CRng m_lrng;
	CRng *m_prng;

// ----------------------------------------------------------------------------
// major multislice parameters

	// high tension [kV] or elecron kinetic energy [keV]
	float m_ht;
	// electron wave lenght [nm] (linked to _JMS_ht)
	float m_wl;
	// super-cell dimensions [nm] (x, y, z)
	float m_a0[3];
	// number of slices per supercell
	int m_nscslc;
	// number of configurations per slice (max. used)
	// The library will not accept more variants than this number.
	// It can be used to limit memory consumption.
	int m_nslcvar;
	// number of slices in object
	int m_nobjslc;
	// number of pixels along a0[0]
	int m_nscx;
	// number of pixels along a0[1]
	int m_nscy;
	// number of phase grating pixels along a0[0] but same sampling rate
	int m_npgx;
	// number of phase grating pixels along a0[1] but same sampling rate
	int m_npgy;
	// number of propagator functions
	int m_npro;
	// used propagator type (0: geometrical (default), 1: Fresnel)
	int m_nprotype;
	// number of integrating detectors
	int m_ndet;
	// detector readout period (used also for wave extraction)
	int m_ndetper;
	// intensity distribution output flag (0: off, 1: integrating, 2: image, 4: diffraction, 8: wave real, 16: wave fourier, 32: averaged image, 64: averaged diffraction, 128: averaged wave real, 256: averaged wave fourierand all OR combinations of these numbers)
	int m_imagedet;
	// number of slices with detection (sum of _JMS_objslc_det)
	int m_ndetslc;
	// status flag for CPU multislice setup
	int m_status_setup_CPU;
	// status flag for GPU multislice setup
	int m_status_setup_GPU;
	// global diffraction de-scan flag (0: off, 1: on)
	int m_dif_descan_flg;
	// plasmon scattering Monte-Carlo flag (0: off, 1: on)
	int m_plasmonmc_flg;
public:
	// default message string for output (not initialized)
	char m_msg[_JPG_MESSAGE_LEN];

// ----------------------------------------------------------------------------
// multislice data objects
protected:
	// number of CPU thread output channels initialized
	int m_threads_CPU_out;
	// number of CPU threads initialized (should be less than or equal to m_threads_CPU_out)
	int m_ncputhreads;
	// status flags for CPU multislice calculations
	int* m_status_calc_CPU;
	// status flag for GPU multislice calculation
	int m_status_calc_GPU;
	// list of object slice sequence
	int *m_objslc;
	// hash list for object slice detection ( -1: no detection, otherwise idetslc = index of the receiving array in m_h_det_int, m_h_det_img,m _h_det_dif ... and m_d_* )
	int *m_det_objslc;
	// list of number of variants for each slice
	int *m_nvarslc;
	// list of offsets of structure slice data in device memory (item counts)
	int64_t *m_slcoffset;
	// list of slice thickness values [nm]
	float *m_slcthick;
	// list of propagators indices for each structure slice (length m_nscslc)
	int *m_proidx;
	// host memory holding detector mask lengths (pointers for each detector function to external numbers)
	int *m_detmask_len;
	// host memory 1D int horizontal frequency numbers
	int *m_h_fnx;
	// host memory 1D int vertical frequency numbers
	int *m_h_fny;
	// device memory 2D horizontal frequency numbers in [1/nm]
	float *m_d_knx;
	// device memory 2D vertical frequency numbers in [1/nm]
	float *m_d_kny;
	// list of host memory incident wave functions (one per thread)
	fcmplx *m_h_wav;
	// host memory backup of incident wave function (common for all threads)
	fcmplx *m_h_wav0;
	// host memory of a constant phase plate due to aberrations of the probe-forming optics
	float *m_h_cppl;
	// list of host memory phase gratings (pointers for each slice to external arrays, for each slice several variants via m_nvarslc)
	fcmplx **m_h_pgr;
	// list of host memory holding propagators (pointers for each propagator to external arrays, ID via m_proidx)
	fcmplx **m_h_pro;
	// host memory holding integrating diffraction detector functions (pointers for each detector function to external arrays)
	float **m_h_det;
	// host memory holding integrating diffraction detector masks (pointers for each detector function to external arrays)
	int **m_h_detmask;
	// per thread list of host memory holding integrated detector results for all detection thicknesses ( iThread, -> idet + idetslc*m_ndet )
	float *m_h_det_int;
	// per thread list of host memory holding probe intensity distributions for all detection thicknesses ( iThread, -> idetslc*m_nscx*m_nscy )
	float *m_h_det_img;
	// per thread list of host memory holding probe diffraction patterns all detection thicknesses ( iThread, -> idetslc*m_nscx*m_nscy )
	float *m_h_det_dif;
	// per thread list of host memory holding real-space wave functions for all detection thicknesses ( iThread, -> idetslc*m_nscx*m_nscy )
	fcmplx *m_h_det_wfr;
	// per thread list of host memory holding Fourier-space wave functions for all detection thicknesses ( iThread, -> idetslc*m_nscx*m_nscy )
	fcmplx *m_h_det_wff;
	// per thread list of host memory holding averaged probe intensity distributions for all detection thicknesses ( iThread, -> idetslc*m_nscx*m_nscy )
	float *m_h_det_img_avg;
	// per thread list of host memory holding averaged probe diffraction patterns all detection thicknesses ( iThread, -> idetslc*m_nscx*m_nscy )
	float *m_h_det_dif_avg;
	// per thread list of host memory holding averaged real-space wave functions for all detection thicknesses ( iThread, -> idetslc*m_nscx*m_nscy )
	fcmplx *m_h_det_wfr_avg;
	// per thread list of host memory holding averaged Fourier-space wave functions for all detection thicknesses ( iThread, -> idetslc*m_nscx*m_nscy )
	fcmplx *m_h_det_wff_avg;
	// per thread list of weights on host image averging channels ( iThread )
	float *m_h_weight_imgavg;
	// per thread list of weights on host wave function averging channels ( iThread )
	float* m_h_weight_wavavg;
	// per thread list of host memory diffraction de-scan in x direction [pixels]
	int *m_h_dif_ndescanx;
	// per thread list of host memory diffraction de-scan in y direction [pixels]
	int *m_h_dif_ndescany;
	// device memory holding the incident wave function (one)
	cuComplex *m_d_wav;
	// device memory holding a backup of the incident wave function (one)
	cuComplex *m_d_wav0;
	// device memory of a constant phase plate due to aberrations of the probe-forming optics
	float *m_d_cppl;
	// device memory holding phase gratings ( -> m_slcoffset[islc] + ivar*m_nscx*m_nscy )
	cuComplex *m_d_pgr;
	// device memory holding propagators ( -> ipro*m_nscx*m_nscy )
	cuComplex *m_d_pro;
	// device memory holding detector functions ( -> idet*m_nscx*m_nscy )
	float *m_d_det;
	// device memory holding detector masks (use m_detmask_len to get detector mask lengths)
	int *m_d_detmask;
	// memory holding integrated detector results for all detection thicknesses ( -> idet + idetslc*ndet )
	// This is for GPU calculations.
	// ! Warning ! The memory of this buffer is on host, since the integration routine output is on host scalars.
	float *m_d_det_int;
	// device memory holding probe intensity distributions for all detection thicknesses ( -> idetslc*m_nscx*m_nscy )
	float *m_d_det_img;
	// device memory holding probe diffraction patterns for all detection thicknesses ( -> idetslc*m_nscx*m_nscy )
	float *m_d_det_dif;
	// device memory holding real-space wave functions for all detection thicknesses ( -> idetslc*m_nscx*m_nscy )
	cuComplex *m_d_det_wfr;
	// device memory holding Fourier-space wave functions for all detection thicknesses ( -> idetslc*m_nscx*m_nscy )
	cuComplex *m_d_det_wff;
	// device memory holding averaged probe intensity distributions for all detection thicknesses ( -> idetslc*m_nscx*m_nscy )
	float *m_d_det_img_avg;
	// device memory holding averaged probe diffraction patterns for all detection thicknesses ( -> idetslc*m_nscx*m_nscy )
	float *m_d_det_dif_avg;
	// device memory holding averaged real-space wave functions for all detection thicknesses ( -> idetslc*m_nscx*m_nscy )
	cuComplex *m_d_det_wfr_avg;
	// device memory holding averaged Fourier-space wave functions for all detection thicknesses ( -> idetslc*m_nscx*m_nscy )
	cuComplex *m_d_det_wff_avg;
	// weight on device image averging channels
	float m_d_weight_imgavg;
	// weight on device wave function averging channels
	float m_d_weight_wavavg;
	// device memory used temporary for readout steps (managed by InitCore)
	float *m_d_det_tmp;
	// device memory used temporary for readout steps (managed by InitCore)
	float *m_d_det_tmp2;
	// device memory used temporary for readout steps (managed by InitCore)
	cuComplex *m_d_det_tmpwav;
	// device diffraction de-scan in x-direction pixels
	int m_d_dif_ndescanx;
	// device diffraction de-scan in y-direction pixels
	int m_d_dif_ndescany;
	// GPU phasegrating source flag (0: device memory, 1: host memory -> copy-on-demand)
	int m_d_pgr_src;

// ----------------------------------------------------------------------------
// multislice objects doing calculations

	// Plasmon scattering Monte-Carlo handler template
	CJPlasmonMC m_jplmc;
	// Probe function calculations
	CJProbeGen m_jpg;
	// CPU core objects for FFTs and array ops on host
	//CJFFTWcore* m_jcpuco;
	CJFFTMKLcore* m_jcpuco;
	// GPU core object for FFTs and array ops on device
	CJFFTCUDAcore m_jgpuco;


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// member functions

public:

// ----------------------------------------------------------------------------
// library handling

	// Set debug level
	int SetDebugLevel(int dbg);

// ----------------------------------------------------------------------------
// cleanup functions

	// Cleans the complete setup (CPU & GPU) !!!
	// - all arrays will be de-allocated on device and host
	// - scalar values remain if not linked to the array handling
	// - resets the setup status to NONE
	int Cleanup(void);

	// Cleans the FFTW module
	// - call this only if all FFTW/MKL routines are halted and no longer used.
	// - FFTW/MKL needs to be reinitialized afterwards.
	void FreeLibMem(void);

// ----------------------------------------------------------------------------
// multislice setup functions

	// Sets high tension in [keV] and returns previous high tension value
	// - also changes wave length
	float SetHighTension(float ht);
	
	// Sets size of the super cell [nm]
	void SetSupercellSize(float *a0);
	
	// Sets size of the super cell [nm]
	void SetSupercellSize(float a, float b, float c);
	
	// Sets calculation grid dimensions in global parameters
	// - nx, ny: grid plane size (wave functions, transmission functions, ...)
	void SetGridSize(int nx, int ny);

	// Set the slice thickness in nm
	// - islc: a slice index valid in the current phase grating setup
	// - fthickness: the thickness of slice islc in nm
	void SetSliceThickness(int islc, float fthickness);

	// Sets the detection planes (overload 1)
	// - ndetper: detection slice period
	// - nobjslc: number of slices in the object
	// - det_objslc: detection hash list, calling with NULL will avoid access
	//               if provided, det_objslc should point to an array of length nobjslc+1
	int SetDetectionPlanes(int ndetper, int nobjslc, int * det_objslc=NULL);

	// Sets the detection planes (overload 2)
	// - slc_det: detection slices listed as std::vector<int>
	// - nobjslc: number of slices in the object
	// - det_objslc: detection hash list, calling with NULL will avoid access
	//               if provided, det_objslc should point to an array of length nobjslc+1
	int SetDetectionPlanes(std::vector<int> slc_det, int nobjslc, int * det_objslc = NULL);

	// Set plasmon scattering Monte-Carlo parameters and switch
	// - use: flag the function on and off
	// - q_e: characteristic scattering vector [1/nm]
	// - q_c: critical scattering vector [1/nm]
	// - mfp: mean free path for single scattering [nm]
	// - nexmax: max. number of excitations per probing electron
	void SetPlasmonMC(bool use, float q_e, float q_c, float mfp, unsigned int nexmax);


	// Turns the diffraction de-scan on or off
	void DiffractionDescan(bool bActivate);

	// Set the diffraction descan values for a specific calculation thread.
	// Use the same beam tilts as in the probe tilt, the routine inverts it.
	// Requires successful InitCore.
	// - whichcode: flag code type, either _JMS_CODE_CPU or _JMS_CODE_GPU
	// - ndescanx, ndescany: horizontal and vertical descan [pixel]
	// - iThread: thread ID (for calls with whichcode == _JMS_CODE_CPU)
	void SetDiffractionDescanN(int whichcode, int ndescanx, int ndescany, int iThread = 0);

	// Set the diffraction descan values for a specific calculation thread.
	// Use the same beam tilts as in the probe tilt, the routine inverts it.
	// Requires successful InitCore.
	// - whichcode: flag code type, either _JMS_CODE_CPU or _JMS_CODE_GPU
	// - descanx, descany: horizontal and vertical descan [1/nm]
	// - iThread: thread ID (for calls with whichcode == _JMS_CODE_CPU)
	void SetDiffractionDescan(int whichcode, float descanx, float descany, int iThread = 0);

	// Set the diffraction descan values for a specific calculation thread.
	// Use the same beam tilts as in the probe tilt, the routine inverts it.
	// Requires successful InitCore.
	// - whichcode: flag code type, either _JMS_CODE_CPU or _JMS_CODE_GPU
	// - descanx, descany: horizontal and vertical descan [mrad]
	// - iThread: thread ID (for calls with whichcode == _JMS_CODE_CPU)
	void SetDiffractionDescanMRad(int whichcode, float descanx, float descany, int iThread = 0);

protected:
	// Load a detector sensitivity profile from a text file
	int LoadSTEMDetectorProfile(std::string sfile, int &len, float &refpix, float** profile);

	// Calculates the radial detector sensitivity for a given diffraction angle theta
	// depending on present profile data.
	// - theta: diffraction angle [rad]
	// - profile: list of sensitivity values
	// - len: length of the profile list
	// - refpix: 
	// - beta1: reference pixel collection angle [rad]
	float GetRadialDetectorSensitivity(float theta, float* profile, int len, float refpix, float beta1);

public:
	// Calculates a STEM probe wavefunction in Fourier space using current parameters
	// - prm: address of a CJProbeParams object defining physical probe parameters
	// - wav: address receiving the probe wave function
	int CalculateProbeWaveFourier(CJProbeParams* prm, fcmplx *wav);

	// Calculates a propagator function for a given thickness, object tilt
	// and current size and wavelength parameters. The propagator amplitudes
	// will be scaled by 1./(m_nscx * m_nsxy) to preserve wave function norms
	// during each multislice step.
	// - fthick: propagation length in nm
	// - otx, oty: object tilt parameters in deg
	// - pro: address receiving the propagator grid
	// - ntype: propagator type switch: 0 = planar, 1 = Fresnel
	int CalculatePropagator(float fthick, float otx, float oty, fcmplx *pro, int ntype=0);

	// Calculates a ring detector function and writes it to a pre-allocated output
	// array 'det' of size (m_nscx * m_nsxy). Also creates a mask if a pointer
	// is provided to an array 'msk'. The function returns the length of the mask
	// to 'msklen', but you should pre-allocate 'msk' to be able to hold (m_nscx * m_nsxy)
	// items.
	// - beta0, beta1: ring collection range in mrad
	// - phi0, phi1: ring segment range in deg
	// - theta0x, theta0y: ring center in mrad
	// - sdsprofile: sensitivity profile file name
	// - det: pointer to a float array receiving the detector function data (pre-allocated)
	// - msklen: length of an access pixel mask list
	// - msk: pointer to an access pixel mask list, which speeds up the detector readout, optional (pre-allocated)
	// - kmom1: first moment pattern index (0 by default)
	int CalculateRingDetector(float beta0, float beta1, float phi0, float phi1, float theta0x, float theta0y, std::string sdsprofile, float *det, int &msklen, int *msk=NULL, int kmom1=0);

	// Prepare phase grating data setup
	// This function is common for both codes, CPU and GPU.
	// Which code is prepared is triggered by parameter whichcode.
	// All arrays are allocated relevant to phase grating data and access to it.
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	// - nx, ny: horizontal and vertical number of phasegrating samples (may differ from wave function grid)
	// - nslc: number of phase gratings
	// - nvarmax: max. number of variants allowed for each slice
	// - nslcvar: list of number of configurations for each phase grating
	int PhaseGratingSetup(int whichcode, int nx, int ny, int nslc, int nvarmax, int* nslcvar);

	// Sets the object slice sequence
	// - nobjslc: max. number of object slices (max. thickness)
	// - objslc: list of structure slice IDs, one for each object slice
	int ObjectSliceSetup(int nobjslc, int* objslc);

	// Prepare propagator setup
	// This function is common for both codes, CPU and GPU.
	// Which code is prepared is triggered by parameter whichcode.
	// All arrays are allocated relevant to propagator data and access to it.
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	// - npro: number of propagator functions
	// - proidx: list of propagator indices handling the access to a propagator function for each structure slice, length _JMS_nscslc
	int PropagatorSetup(int whichcode, int npro, int *proidx);

	// Prepare detector setup (overload 1)
	// This function is common for both codes, CPU and GPU.
	// Which code is prepared is triggered by parameter whichcode.
	// All arrays are allocated relevant to detection data and access to it.
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	// - ndetper: detector readout periodicity in slices (0 = exit plane only, >0 many planes + final plane)
	// - ndetint: number of integrating diffraction plane detectors
	// - imagedet: image detection switches (_JMS_DETECT_IMAGE | _JMS_DETECT_DIFFRACTION )
	// - nthreads_CPU: number of CPU output threads to initialize
	int DetectorSetup(int whichcode, int ndetper, int ndetint, int imagedet, int nthreads_CPU = 1);

	// Prepare detector setup (overload 2)
	// This function is common for both codes, CPU and GPU.
	// Which code is prepared is triggered by parameter whichcode.
	// All arrays are allocated relevant to detection data and access to it.
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	// - ndetper: detector readout periodicity in slices (0 = exit plane only, >0 many planes + final plane)
	// - ndetint: number of integrating diffraction plane detectors
	// - imagedet: image detection switches (_JMS_DETECT_IMAGE | _JMS_DETECT_DIFFRACTION )
	// - nthreads_CPU: number of CPU output threads to initialize
	int DetectorSetup(int whichcode, std::vector<int> slc_det, int ndetint, int imagedet, int nthreads_CPU = 1);

	// Sets phase grating data from external source.
	// Call for each slice separately!
	// For CPU code, this will not allocate more memory, the external source will be used, so keep it alive!
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	// - islc: structure slice index
	// - nvar: number of variants contained by the input pgr
	// - pgr: pointer to external source
	int SetPhaseGratingData(int whichcode, int islc, int nvar, fcmplx* pgr);

	// Sets propagator data from external source.
	// Call for each propagator ID separately!
	// For CPU code, this will not allocate more memory, the external source will be used, so keep it alive!
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	// - ipro: propagator index
	// - pro: pointer to external source
	int SetPropagatorData(int whichcode, int ipro, fcmplx* pro);

	// Sets detector function data from external source.
	// Call for each detector function ID separately!
	// For CPU code, this will not allocate more memory, the external source will be used, so keep it alive!
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	// - idet: detector function index
	// - det: detector function (sensitivity) pointer to external source
	// - msklen: detector mask length (optional)
	// - msk: detector mask pointer to external source (optional)
	int SetDetectorData(int whichcode, int idet, float* det, int msklen = 0, int* msk = NULL);

	// Initializes calculation core objects
	// Also prepares arrays holding incident wave functions.
	// Assumes previous setup via JMS_SetGridSize, JMS_SetHighTension, JMS_SetSupercellSize
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	// - nCPUthreads: number of CPU threads possible (ignored if only _JMS_CODE_GPU)
	int InitCore(int whichcode, int nCPUthreads = 1);

	// Stores a copy backup of the incident wave function
	// Makes copies of the wave function on host _h_JMS_wav0 and device _d_JMS_wav0.
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	// - wav: pointer to wave function data (assumed to be in Fourier space)
	// - bTranspose: flag causing a 2d transpose from input to module buffer (default: false)
	int SetIncidentWave(int whichcode, fcmplx* wav, bool bTranspose=false);

	// Stores the phase plate to be appplied to the incident wave function
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	// - ppl: pointer to phase plate data (scrambled Fourier space)
	// - bTranspose: flag causing a 2d transpose from input to module buffer (default: false)
	int SetIncidentWavePhaseplate(int whichcode, float* ppl, bool bTranspose = false);

	// Sets the phase plate to be appplied to the incident wave function to zero
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	int ZeroIncidentWavePhaseplate(int whichcode);

	// Returns a hash table to unscramble the Fourier space by means of index access.
	// After unscrambling, the zero beam is on pixel (m_nscx/2, m_nscy/2).
	// Use this hash like this: a_unscrambled[idx] = a_scrambled[phash[idx]];
	// - phash: unsigned int list of indices in the order of their occurrance in
	//          an unscrambled array (preallocated by the calling code)
	int GetUnscrambleHash(unsigned int* phash);


// ----------------------------------------------------------------------------
// multislice Get setup functions

	// Returns number of registered detection slices
	int GetDetetionSliceNum(void);
	// Returns number of registered detectors
	int GetDetNum(void);
	// Returns the current high tension value in kV
	float GetHighTension(void);
	// Returns the current wave length in nm
	float GetWaveLength(void);
	// Returns the current grid size
	void GetGridSize(int &nx, int &ny);
	// Returns the current super-cell size in nm
	void GetSupercellSize(float* a0);
	// Returns the current super-cell size in nm
	void GetSupercellSize(float &a, float &b, float &c);
	// Returns the thickness of slice islc in nm (0.0 if not set)
	float GetSliceThickness(int islc);
	


// ----------------------------------------------------------------------------
// CUDA device interface functions
// - in case of I/O via the function interface, return values of 0 indicate
//   success and return values of >0 indicate failure of the routine, error
//   codes may be returned

protected:
	// Posts a CUDA error message to the standard error console
	void PostCUDAError(const char* smsg, cudaError code);
	// Posts a CUDA memory allocation problem
	void PostCUDAMemory(size_t nrequestedbytes);
public:
	// Returns number of available CPU devices
	int GetCPUNum(void);
	// Returns number of initiated CPU core objects
	int GetJCPUcorenum(void);
	// Returns number of available GPU devices
	int GetGPUNum(void);
	// Returns name of a GPU device to name (allocate 256 bytes for that)
	int GetGPUName(int idev, char *name);
	// Gets the current CUDA device used
	int GetCurrentGPU(void);
	// Sets the current CUDA device to be used
	int SetCurrentGPU(int idev);
	// Returns a few statistical numbers on the GPU device
	int GetGPUStats(int idev, int &iCMajor, int &iCMinor, int &iMaxThread, int64_t &CUDAmemtotal, int64_t &CUDAmemfree);
	//// Returns core numbers on the GPU device
	//int GetGPUCores(int idev, int &nMultiProcs, int &nCores, int& nMaxThreadPerProc);
	// Returns memory info on current GPU device
	int GetGPUMemInfo(size_t &memtotal, size_t &memfree);
	// Sets GPU phase grating loading scheme
	// - npgrload: 0 = pre-load all to device (default),
	//             1 = each phase-grating is loaded to device on demand
	int SetGPUPgrLoading(int npgrload = 0);


// ----------------------------------------------------------------------------
// multislice utility functions (usually not exposed to extern)

	// Sets the active random number generator object
	// Calling with NULL pointer resets to the local CRng member m_lrng.
	void SetRng(CRng *prng=NULL);

	// Returns a pointer to the current random number generator object
	CRng* GetRng(void);

	// Re-Initializes the local random number generator
	void SeedRngEx(int nseed = 0);

protected:
	// Returns a random sequence of variant numbers for a given sequence of
	// object slices based on the current slice variant data
	// - out: (output) list of variant indices for each slice in input objslc
	int GetRandomVariantSequence(int *out);

	// Allocates host memory to pointer _h_a
	// - _h_a: returned host memory address
	// - size: memory size to allocate in bytes
	// - callfn: calling function name (for error report)
	// - arrnam: array usage (for error report)
	// - zeroe: flag causing the routine to preset the array with zeroes
	int AllocMem_h(void ** _h_a, size_t size, const char* callfn, const char* arrnam, bool zero=false);

	// Allocates device memory to pointer _d_a
	// - _d_a: returned device memory address
	// - size: memory size to allocate in bytes
	// - callfn: calling function name (for error report)
	// - arrnam: array usage (for error report)
	// - zeroe: flag causing the routine to preset the array with zeroes
	int AllocMem_d(void ** _d_a, size_t size, const char* callfn, const char* arrnam, bool zero = false);

	// Deallocate host memory registered for address _h_a
	// - _h_a: host memory address to free (will be nulled afterwards)
	void DeallocMem_h(void ** _h_a);

	// Deallocate device memory registered for address _d_a
	// - _d_a: device memory address to free (will be nulled afterwards)
	void DeallocMem_d(void ** _d_a);

public:
	// clears host detector memory for a thread
	// - thread ID is not checked
	int ClearDetMem_h(int iThread);

	//// clears averaging detector memory for a thread
	//// - thread ID is not checked
	//int ClearDetAvgMem_h(int iThread);
	
	// clears device detector memory
	int ClearDetMem_d(void);

	//// clears device detector memory
	//int ClearDetAvgMem_d(void);
	
protected:
	// returns host memory phase grating for a slice
	fcmplx* GetPhaseGrating_h(int iSlice, int* pVarID);
	
	// returns device memory phase grating for a slice
	cuComplex* GetPhaseGrating_d(int iSlice, int* pVarID);
	
	// returns host memory propagator for a slice
	fcmplx* GetPropagator_h(int iSlice);
	
	// returns device memory propagator for a slice
	cuComplex* GetPropagator_d(int iSlice);
	
	// returns the absolute square total of a wave function array (free)
	// - wav: pointer to wave function data
	// - len: number of values
	float GetAbsTotal(fcmplx* wav, size_t len);
	
	// returns the absolute total of a float array (free)
	// - dat pointer to the array
	// - len: number of values
	float GetTotalF(float* dat, size_t len);
	
	// returns the absolute square total of the current CPU core wave function
	// ! Take care to only use from the correct thread to avoid synchroneous access !
	float GetCoreAbsTotal_h(int iThread);
	
	// returns the absolute square total of the current GPU core wave function
	// ! Take care to only use from the correct thread to avoid synchroneous access !
	float GetCoreAbsTotal_d(void);
	
	// returns the dot product of two vectors on host memory
	float DotProduct_h(float *in_1, float *in_2, size_t len);
	
	// returns the dot product of two vectors on device memory
	float DotProduct_d(float *in_1, float *in_2, size_t len, int nBlockSize);
	// change alternative dot product working from complex data directly
	float DotProduct_d(float2* in_1, float* in_2, size_t len, int nBlockSize);
	// end change

	// returns the masked dot product of two vectors on host memory
	float MaskedDotProduct_h(int *mask, float *in_1, float *in_2, size_t lenmask);

	// returns the masked dot product of two vectors on device memory
	float MaskedDotProduct_d(int *mask, float *in_1, float *in_2, size_t lenmask, int nBlockSize);
	// change alternative dot product for detector readout, working on complex data
	float MaskedDotProduct_d(int* mask, float2* in_1, float* in_2, size_t lenmask, int nBlockSize);
	// end change
	
	// performs diffraction detector readouts for CPU thread in given slice
	int ReadoutDifDet_h(int iSlice, int iThread, float weight = 1.0f);
	
	// performs diffraction detector readouts for GPU in given slice
	int ReadoutDifDet_d(int iSlice, float weight = 1.0f);
	
	// performs image detector readouts for CPU thread in given slice
	int ReadoutImgDet_h(int iSlice, int iThread, float weight = 1.0f);
	
	// performs image detector readouts for GPU in given slice
	int ReadoutImgDet_d(int iSlice, float weight = 1.0f);

	// Applies a descan and records the diffraction pattern, storing in dif.
	// Allows only integer pixel shifts.
	// - whichcode: flag signaling which code to work on, either _JMS_CODE_CPU or _JMS_CODE_GPU.
	//              DOESN'T ALLOW COMBINING MODES. RETURNS AN ERROR IF CALLED WITH BOTH MODES.
	// - dnx, dny: de-scan in the diffraction plane in [pixels]
	// - dif: address of the diffraction pattern for output (on device or host depending on whichcode)
	// - iThread: thread ID for CPU code, ignored for GPU code
	int DescanDifN(int whichcode, int dnx, int dny, float *dif, int iThread = 0);

	//int CJMultiSlice::DescanDif(int whichcode, float dx, float dy, float* dif, int iThread)

	// Applies a descan and records the Fourier space wave function, storing in wav on device.
	// Allows only integer pixel shifts.
	// - dnx, dny: de-scan in the diffraction plane in [pixels]
	// - wav: address of the wave function for output (on device)
	int DescanDifWavN_d(int dnx, int dny, cuComplex *wav);
	
	// Applies a descan and records the Fourier space wave function, storing in wav on host.
	// Allows only integer pixel shifts.
	// - dnx, dny: de-scan in the diffraction plane in [pixels]
	// - wav: address of the wave function for output (on device)
	// - iThread: thread ID for CPU code
	int DescanDifWavN_h(int dnx, int dny, fcmplx *wav, int iThread);

// ----------------------------------------------------------------------------
// multislice functions

public:
	// Takes a copy of the backup wave function m_h_wav0 / m_d_wav0,
	// applies offsets in (x, y, z), and stores the result in the the active
	// wave function channels m_h_wav / m_d_wav used for the multislice calculation.
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_CPU | _JMS_CODE_GPU)
	// - dx, dy, dz: offset distances in 3 dimensions and nm units.
	// - iThread: thread ID for CPU code, ignored for GPU code
	// - ibtx, ibty: offset tilt in Fourier pixels. Defaults to zero. (added 2022-Apr-07 for SPED simulations)
	// - bcppl: flag to use the constant phase plates stored in m_h_cppl / m_d_cppl. Defaults tu true. (added 2022-Apr-14 for SPED simulations)
	int OffsetIncomingWave(int whichcode, float dx, float dy, float dz, int iThread = 0, int ibtx = 0, int ibty = 0, bool bcppl = true);

	// Ignores the backup wave function m_d_wav0 and stores wav directly in
	// the wave function channel m_d_wav for the GPU multislice calculation.
	// - wav: wave function in Fourier space
	// - bTranspose: flag signalizing that the input wave function is transposed
	int SetIncomingWaveGPU(fcmplx* wav, bool bTranspose=false);

	// Ignores the backup wave function m_h_wav0 and stores wav directly in
	// the wave function channel m_h_wav of thread iThread for the CPU multislice calculation.
	// - wav: wave function in Fourier space
	// - bTranspose: flag signalizing that the input wave function is transposed
	// - iThread: CPU thread ID
	int SetIncomingWaveCPU(fcmplx* wav, bool bTranspose=false, int iThread=0);

	// Runs a multislice calculation on a CPU thread
	// Assumes incident wave function present in _h_JMS_wav
	// - islc0: index of the starting slice
	// - accmode = 0: new accumulation, 1: accumulate to previous detections
	// - weight = weighting factor for the accumulation of this run to previous data
	// - iThread: thread ID
	int CPUMultislice(int islc0, int accmode, float weight=1.0f, int iThread=0);

	// Runs a multislice calculation on GPU.
	// Assumes incident wave function present in _d_JMS_wav.
	// - islc0: index of the starting slice
	// - accmode = 0: new accumulation, 1: accumulate to previous detections
	// - weight = weighting factor for the accumulation of this run to previous data
	int GPUMultislice(int islc0, int accmode, float weight = 1.0f);

	// Copies a set of detection results to provided host address
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_GPU or _JMS_CODE_CPU)
	// - whichresult: flag signaling which result to retrieve, one of _JMS_DETECT_INTEGRATED, _JMS_DETECT_IMAGE, _JMS_DETECT_DIFFRACTION, ...
	// - dst: destination address recieving results.
	// - iThread: thread ID of CPU code
	int GetResult(int whichcode, int whichresult, float *dst, int iThread = 0);

	// Copies a set of averaged detection results to provided host address
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_GPU or _JMS_CODE_CPU)
	// - whichresult: flag signaling which result to retrieve, one of _JMS_DETECT_IMAGE_AVG, _JMS_DETECT_DIFFR_AVG, ...
	// - dst: destination address recieving results.
	// - weight: floatvariable receiving the weight of averging
	// - iThread: thread ID of CPU code
	int GetAvgResult(int whichcode, int whichresult, float* dst, float& weight, int iThread = 0);

	// Returns the image averaging weight of a given code channel
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_GPU or _JMS_CODE_CPU)
	// - iThread: thread ID of CPU code
	float GetImageAveragingWeight(int whichcode, int iThread=0);

	// Returns the wave function averaging weight of a given code channel
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_GPU or _JMS_CODE_CPU)
	// - iThread: thread ID of CPU code
	float GetWaveAveragingWeight(int whichcode, int iThread = 0);


	// Resets active image averaging channels and weight counters for a given code channel
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_GPU or _JMS_CODE_CPU)
	// - iThread: thread ID of CPU code
	// Returns an error code or 0 in case of success.
	int ResetImageAveraging(int whichcode, int iThread = 0);

	// Resets active wave function averaging channels and weight counters for a given code channel
	// - whichcode: flag signaling which code to prepare (_JMS_CODE_GPU or _JMS_CODE_CPU)
	// - iThread: thread ID of CPU code
	// Returns an error code or 0 in case of success.
	int ResetWaveAveraging(int whichcode, int iThread = 0);

	// Reads detector intensities from the previously calculated average wave functions
	// - whichcode: flag signaling which code to prepare (one of _JMS_CODE_GPU or _JMS_CODE_CPU)
	// - whichresult: flag signaling which detector is to be read out (one of _JMS_DETECT_INTEGRATED, _JMS_DETECT_IMAGE, or _JMS_DETECT_DIFFRACTION)
	// - iThread: thread ID of CPU code
	// Returns an error code or 0 in case of success.
	// Writes the intensities to buffer dst, assuming that a sufficiently large amount of memory is allocated.
	// Writes the weight of the data to weight.
	int ReadoutDetAvg(int whichcode, int whichresult, float *dst, float &weight, int iThread = 0);

};

#endif