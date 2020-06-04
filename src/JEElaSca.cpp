//
// C++ source file: JEElaSca.cpp
// implementation for library JMultislice.lib (declarations see JEElaSca.h)
//
//
// Copyright (C) 2020 - Juri Barthel (juribarthel@gmail.com)
// Copyright (C) 2020 - Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
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

#include "JEElaSca.h"
#include "cu\scapot.cuh"
#include <algorithm>

CJEElaSca::CJEElaSca()
{
	m_status = (unsigned int)EELSCA_STATUS_NOTREADY;

	m_vCellSize = CV3D();
	
	m_grid_x = 0;
	m_grid_y = 0;
	m_grid_z = 0;
	
	m_bwl_pot = false;
	m_bwl_pgr = true;
	mod_atffacs = 0;
	mod_thermal = 0;
	mod_absorb = 0;
	m_ekv = 0.f;
	m_vabf = 0.f;

	m_att_num = 0;
	m_cpu_num = 0;
	m_gpu_id = -1;
	
	m_h_qx = NULL;
	m_h_qy = NULL;
	m_h_qz = NULL;
	m_d_qx = NULL;
	m_d_qy = NULL;
	m_d_qz = NULL;
	m_h_q2 = NULL;
	m_h_ap = NULL;
	m_d_ap = NULL;
	m_h_potio = NULL;
	m_d_potio = NULL;
	m_h_att_ff2d = NULL;
	m_d_ff2d = NULL;

	m_jcpuco = NULL;

	m_prng = &m_lrng;
}

CJEElaSca::~CJEElaSca()
{
	ResetAll();
}


void CJEElaSca::SetRng(CRng* prng)
{
	if (prng) {
		m_prng = prng;
	}
	else {
		m_prng = &m_lrng;
	}
}

CRng* CJEElaSca::GetRng()
{
	return m_prng;
}

void CJEElaSca::SeedRngEx(int nseed)
{
	if (m_prng) {
		m_prng->seed(nseed);
	}
	else {
		m_lrng.seed(nseed);
	}
}

float CJEElaSca::UniRand(void)
{
	return (float)m_prng->unirand();
}

float CJEElaSca::NormRand(void)
{
	return (float)m_prng->normrand();
}

void CJEElaSca::ResetAll(void)
{
	DeallocHelper();
	
	m_vCellSize = CV3D(0.f, 0.f, 0.f);
	m_vTypes.clear();

	m_grid_z = 0;
	m_grid_y = 0;
	m_grid_x = 0;

	m_ekv = 0.f;
	m_vabf = 0.f;
	mod_absorb = 0;
	mod_thermal = 0;
	mod_atffacs = 0;

	m_status = (unsigned int)EELSCA_STATUS_NOTREADY;
}


void CJEElaSca::DeallocHelper(void)
{
	if (m_d_qx != NULL) { cudaFree(m_d_qx); m_d_qx = NULL; }
	if (m_d_qy != NULL) { cudaFree(m_d_qy); m_d_qy = NULL; }
	if (m_d_qz != NULL) { cudaFree(m_d_qz); m_d_qz = NULL; }
	if (m_h_qx != NULL) { free(m_h_qx); m_h_qx = NULL; }
	if (m_h_qy != NULL) { free(m_h_qy); m_h_qy = NULL; }
	if (m_h_qz != NULL) { free(m_h_qz); m_h_qz = NULL; }
	if (m_h_q2 != NULL) { free(m_h_q2); m_h_q2 = NULL; }
	if (m_d_ap != NULL) { cudaFree(m_d_ap); m_d_ap = NULL; }
	if (m_h_ap != NULL) { free(m_h_ap); m_h_ap = NULL; }
	if (m_d_potio != NULL) { cudaFree(m_d_potio); m_d_potio = NULL; }
	if (m_h_potio != NULL) { free(m_h_potio); m_h_potio = NULL; }
	if (m_h_att_ff2d != NULL) { free(m_h_att_ff2d); m_h_att_ff2d = NULL; }
	if (m_d_ff2d != NULL) { cudaFree(m_d_ff2d); m_d_ff2d = NULL; }
	if (m_jcpuco != NULL && m_cpu_num > 0) {
		for (int i = 0; i < m_cpu_num; i++) {
			m_jcpuco[i].Deinit();
		}
		delete[] m_jcpuco;
		m_jcpuco = NULL;
		m_cpu_num = 0;
	}
	m_jgpuco.Deinit();
	m_gpu_id = -1;

	m_att_num = 0;

	if (m_status & (unsigned int)EELSCA_STATUS_HELPER_SET) m_status -= (unsigned int)EELSCA_STATUS_HELPER_SET; // clear helper set status
	if (m_status & (unsigned int)EELSCA_STATUS_HELPER_ALLOC) m_status -= (unsigned int)EELSCA_STATUS_HELPER_ALLOC; // clear helper alloc status
}


int CJEElaSca::AllocHelper(int igpu, int ncpu)
{
	int nerr = 0;
	cudaError cuerr;
	size_t nitems = 0;

	if (m_status < (unsigned int)EELSCA_THRESHOLD_ATOMICFF) {
		nerr = 1;
		std::cerr << "Error (CJEElaSca::AllocHelper): Incomplete setup (" << m_status << ")." << std::endl;
		goto _exit;
	}

	DeallocHelper();

	m_att_num = m_vTypes.size();

	// allocate memory for projected 2d form factors on host
	nitems = (size_t)m_grid_x * m_grid_y * m_att_num;
	if (nitems > 0 && 0 == m_pot3d) { // 2d form factors for all atom types
		m_h_att_ff2d = (fcmplx*)malloc(sizeof(fcmplx) * nitems);
		if (NULL == m_h_att_ff2d) {
			std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate host memory (projected atomic form factors)." << std::endl;
			nerr = 202;
			goto _exit;
		}
		memset(m_h_att_ff2d, 0, sizeof(fcmplx) * nitems);
	}
	// allocate other host buffes which we want to have independent of cpu usage (it is not much)
	if (m_grid_x > 0) { // x-shift helper on host
		m_h_qx = (float*)malloc(sizeof(float) * m_grid_x);
		if (NULL == m_h_qx) {
			std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate host memory (shift-x)." << std::endl;
			nerr = 3;
			goto _exit;
		}
		memset(m_h_qx, 0, sizeof(float) * m_grid_x);
	}
	if (m_grid_y > 0) { // y-shift helper on host
		m_h_qy = (float*)malloc(sizeof(float) * m_grid_y);
		if (NULL == m_h_qy) {
			std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate host memory (shift-y)." << std::endl;
			nerr = 4;
			goto _exit;
		}
		memset(m_h_qy, 0, sizeof(float) * m_grid_y);
	}
	if (m_grid_z > 0 && m_pot3d > 0) { // z-shift helper on host
		m_h_qz = (float*)malloc(sizeof(float) * m_grid_z);
		if (NULL == m_h_qz) {
			std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate host memory (shift-z)." << std::endl;
			nerr = 5;
			goto _exit;
		}
		memset(m_h_qz, 0, sizeof(float) * m_grid_z);
	}
	nitems = (size_t)m_grid_x * m_grid_y;
	if (nitems > 0 && 0 == m_pot3d) { // q^2 helper
		m_h_q2 = (float*)malloc(sizeof(float) * nitems);
		if (NULL == m_h_q2) {
			std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate host memory (q^2 helper)." << std::endl;
			nerr = 8;
			goto _exit;
		}
		memset(m_h_q2, 0, sizeof(float) * nitems);
	}
	nitems = (size_t)m_grid_x * m_grid_y * 2;
	if (nitems > 0 && 0 == m_pot3d) { // aperture helper
		m_h_ap = (float*)malloc(sizeof(float) * nitems);
		if (NULL == m_h_ap) {
			std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate host memory (aperture helper)." << std::endl;
			nerr = 7;
			goto _exit;
		}
		memset(m_h_ap, 0, sizeof(float) * nitems);
	}
	nitems = (size_t)m_grid_x * m_grid_y;
	if (nitems > 0 && 0 == m_pot3d) { // ionic term
		m_h_potio = (fcmplx*)malloc(sizeof(fcmplx) * nitems);
		if (NULL == m_h_potio) {
			std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate host memory (ionic term)." << std::endl;
			nerr = 8;
			goto _exit;
		}
		memset(m_h_potio, 0, sizeof(fcmplx) * nitems);
	}

	if (igpu >= 0) { // helper allocations needed on device
		cuerr = cudaSetDevice(igpu);
		if (cudaSuccess != cuerr) {
			std::cerr << "Error (CJEElaSca::AllocHelper): Select CUDA device #" << igpu << ". " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
			nerr = 102;
			goto _exit;
		}
		if (m_grid_x > 0) { // x-shift helper on device
			cuerr = cudaMalloc(&m_d_qx, sizeof(float) * m_grid_x);
			if (cuerr != cudaSuccess) {
				std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate device memory (shift-x). " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
				nerr = 103;
				goto _exit;
			}
			cuerr = cudaMemset(m_d_qx, 0, sizeof(float) * m_grid_x);
		}
		if (m_grid_y > 0) { // y-shift helper on device
			cuerr = cudaMalloc(&m_d_qy, sizeof(float) * m_grid_y);
			if (cuerr != cudaSuccess) {
				std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate device memory (shift-y). " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
				nerr = 104;
				goto _exit;
			}
			cuerr = cudaMemset(m_d_qy, 0, sizeof(float) * m_grid_y);
		}
		if (m_grid_z > 0 && m_pot3d > 0) { // z-shift helper on device
			cuerr = cudaMalloc(&m_d_qz, sizeof(float) * m_grid_z);
			if (cuerr != cudaSuccess) {
				std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate device memory (shift-z). " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
				nerr = 105;
				goto _exit;
			}
			cuerr = cudaMemset(m_d_qz, 0, sizeof(float) * m_grid_z);
		}
		nitems = (size_t)m_grid_x * m_grid_y * 2;
		if (nitems > 0 && 0 == m_pot3d) { // aperture helper on device
			cuerr = cudaMalloc(&m_d_ap, sizeof(float) * nitems);
			if (cuerr != cudaSuccess) {
				std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate device memory (apertures). " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
				nerr = 107;
				goto _exit;
			}
			cuerr = cudaMemset(m_d_ap, 0, sizeof(float)* nitems);
		}
		nitems = (size_t)m_grid_x * m_grid_y;
		if (nitems > 0 && 0 == m_pot3d) { // ionic term on device
			cuerr = cudaMalloc(&m_d_potio, sizeof(fcmplx) * nitems);
			if (cuerr != cudaSuccess) {
				std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate device memory (ionic term). " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
				nerr = 108;
				goto _exit;
			}
			cuerr = cudaMemset(m_d_potio, 0, sizeof(fcmplx)* nitems);
		}
		nitems = (size_t)m_grid_x * m_grid_y;
		if (nitems > 0 && 0 == m_pot3d) { // form factor helper on device
			cuerr = cudaMalloc(&m_d_ff2d, sizeof(fcmplx) * nitems);
			if (cuerr != cudaSuccess) {
				std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate device memory (2d form factors). " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
				nerr = 110;
				goto _exit;
			}
			cuerr = cudaMemset(m_d_ff2d, 0, sizeof(fcmplx)* nitems);
		}
		m_gpu_id = igpu;
	}

	// * TODO * //
	// Rethink how to use many CPU threads for 3D potential calculations. Perhaps each running on a kz-plane of the 3d Fourier space.

	if (ncpu > 0) { // allocate host buffers depending on cpu usage
		//m_jcpuco = (CJFFTWcore*)malloc(sizeof(CJFFTWcore) * ncpu);
		m_jcpuco = new CJFFTMKLcore[ncpu]; // initialize the FFT cores
		if (NULL == m_jcpuco) {
			std::cerr << "Error (CJEElaSca::AllocHelper): Failed to allocate host memory (FFT cores)." << std::endl;
			nerr = 11;
			goto _exit;
		}
		m_cpu_num = ncpu;
	}

	m_status |= (unsigned int)EELSCA_STATUS_HELPER_ALLOC;

_exit:
	return nerr;
}

int CJEElaSca::SetStructure(CStructure st)
{
	int nerr = 0;
	if (m_status & (unsigned int)EELSCA_STATUS_STRUCTURE) { // there is already structure information set
		// return with error
		nerr = 1;
		std::cerr << "Error (CJEElaSca::SetStructure): Structure data is already set. Call ResetAll before setting new data." << std::endl;
		goto _exit;
	}
	m_vCellSize = st.m_vCellSize;
	m_vTypes = st.m_vTypes;
	m_status |= (unsigned int)EELSCA_STATUS_STRUCTURE;
_exit:
	return nerr;
}

int CJEElaSca::SetGrid(unsigned int nx, unsigned int ny, unsigned int nz, unsigned int pot3d)
{
	int nerr = 0;
	if (m_status & (unsigned int)EELSCA_STATUS_GRID) { // there is already grid information set
		// return with error
		nerr = 1;
		std::cerr << "Error (CJEElaSca::SetGrid): Grid sampling is already set. Call ResetAll before setting new data." << std::endl;
		goto _exit;
	}
	m_grid_x = nx;
	m_grid_y = ny;
	m_grid_z = 1;
	m_pot3d = pot3d;
	if (m_pot3d > 0) {
		m_grid_z = nz;
	}
	m_status |= (unsigned int)EELSCA_STATUS_GRID;
_exit:
	return nerr;
}

int CJEElaSca::SetPhysics(float ekv, unsigned int attfacs, unsigned int absorb, float vabf, unsigned int thermal)
{
	int nerr = 0;
	if (m_status & (unsigned int)EELSCA_STATUS_PHYSICS) { // there is already grid information set
		// return with error
		nerr = 1;
		std::cerr << "Error (CJEElaSca::SetPhysics): Physics options are already set. Call ResetAll before setting new data." << std::endl;
		goto _exit;
	}
	mod_atffacs = attfacs;
	mod_absorb = absorb;
	mod_thermal = thermal;
	m_vabf = vabf;
	m_ekv = ekv;
	m_status |= (unsigned int)EELSCA_STATUS_PHYSICS;
_exit:
	return nerr;
}

bool CJEElaSca::UseDWF(void)
{
	return (1 == mod_thermal);
}

float CJEElaSca::GetDWF(float biso, float ksqr)
{
	if (UseDWF()) {
		return std::exp(-0.25f * biso * ksqr);
	}
	return 1.f;
}

bool CJEElaSca::UseQEP(void)
{
	return (2 == mod_thermal);
}

float CJEElaSca::GetRelCorr(float ekv)
{
	float me0 = (float)(_EEL0KEV);
	return (me0 + ekv) / me0;
}

int CJEElaSca::CalculateFormFactors(bool bverbose)
{
	int nerr = 0, ierr = 0;
	size_t iatt = 0;
	size_t natt = m_vTypes.size();
	size_t nk = 0;
	unsigned int nyqx = (m_grid_x - (m_grid_x % 2)) >> 1;
	unsigned int nyqy = (m_grid_y - (m_grid_y % 2)) >> 1;
	unsigned int nyqz = (m_grid_z - (m_grid_z % 2)) >> 1; // x, y, z nyquist numbers
	float kmax = 0.f, sk = 0.f, bwlk = 0.f;
	if (m_status < (unsigned int)EELSCA_THRESHOLD_ATOMICFF) {
		nerr = 1;
		std::cerr << "Error (CJEElaSca::CalculateFormFactors): Incomplete setup (" << m_status << ")." << std::endl;
		goto _exit;
	}
	if (natt == 0) {
		if (bverbose) {
			std::cout << std::endl;
			std::cout << "  No form factors to calculate because there are no atom types." << std::endl;
		}
		m_status |= (unsigned int)EELSCA_STATUS_ATOMICFF;
		goto _exit;
	}
	if (bverbose) {
		std::cout << std::endl;
		std::cout << "  Calculating radial electron form factors for " << natt << " atom types ..." << std::endl;
		if (mod_atffacs == 0) {
			std::cout << "  - Using the parameterization of Weickenmeier and Kohl, Acta Cryst. A 47 (1991) 590-597." << std::endl;
		}
		if (UseDWF()) {
			std::cout << "  - Form factors will be attenuated by Debye-Waller factors." << std::endl;
		}
		if (mod_absorb == 1) {
			std::cout << "  - Calculating absorptive form factors as fraction " << m_vabf << " of the elastic form factors." << std::endl;
		}
		if (mod_absorb == 2) {
			if (UseDWF()) {
				std::cout << "  - Calculating absorptive form factors due to thermal diffuse scattering." << std::endl;
			}
			else {
				std::cout << "  - Calculating absorptive form factors due to the band-width limitation of the calculation." << std::endl;
			}
		}
	}
	// determine sampling
	nk = (size_t)12 * std::max(std::max(m_grid_x, m_grid_y), m_grid_z); // use an oversampling of 12
	// determine the k-space sampling rate as 2 times of the highest scattering angle = 2 * n/2 * 1/a = n / a
	sk = std::max(std::max((float)m_grid_x / m_vCellSize.x, (float)m_grid_y / m_vCellSize.y), (float)m_grid_z / m_vCellSize.z) / (float)nk;
	// determine the maximum scattering vector in the x-y plane from the smaller of the two Nyquist frequencies
	kmax = std::min((float)nyqx / m_vCellSize.x, (float)nyqy / m_vCellSize.y);
	// calculate the 2d band-width limit = 2/3 * kmax
	bwlk = kmax * EELSCA_CONST_PGR_BWL;
	// loop over all atom types
	for (iatt = 0; iatt < natt; iatt++) {
		ierr = m_vTypes[iatt].CalculateRadialFormFactors(nk, sk, bwlk, m_ekv, UseDWF(), mod_absorb, m_vabf, mod_atffacs, bverbose);
		if (0 < ierr) {
			nerr = 2;
			std::cerr << "Error (CJEElaSca::CalculateFormFactors): failed to calculate radial form fractor for " << m_vTypes[iatt].SpeciesName() << "." << std::endl;
			goto _exit;
		}
	}
	m_status |= (unsigned int)EELSCA_STATUS_ATOMICFF;
_exit:
	return nerr;
}


int CJEElaSca::InitHelper()
{
	int nerr = 0, ierr = 0;
	cudaError cuerr;
	float* qx = NULL;
	float* qy = NULL;
	float* qz = NULL;
	fcmplx* cg2d = NULL;
	float* fg2d = NULL;
	size_t npix = 0, iatt = 0, idx = 0;
	int i = 0, j = 0, ny2 = 0, nx2 = 0, nz2 = 0;
	int fi = 0, fj = 0;
	int icpu = 0, ndims = 0;
	int* pdims = NULL;
	float itogx = 1.f / m_vCellSize.x;
	float itogy = 1.f / m_vCellSize.y;
	float itogz = 1.f / m_vCellSize.z;
	float gx = 0.f, gy = 0.f, gy2 = 0.f, g2 = 0.f, gm = 0.f, apxy1 = 1.f, apxy2 = 1.f, kmax = 0.f, pfacio = 0.f, kthr23 = 0.f;
	float ftrn = 1.f / (float)m_grid_x / (float)m_grid_y; // re-normalization factor so that after an FT and an IFT the norm is preserved
	float iyr2 = 1.f / EELSCA_CONST_IYR / EELSCA_CONST_IYR;
	if (0 == (m_status & (unsigned int)EELSCA_STATUS_HELPER_ALLOC)) {
		nerr = 1;
		std::cerr << "Error (CJEElaSca::CalculateHelper): Incomplete setup (" << m_status << ")." << std::endl;
		goto _exit;
	}
	npix = (size_t)m_grid_x * m_grid_y;
	if (npix == 0) {
		nerr = 2;
		std::cerr << "Error (CJEElaSca::CalculateHelper): Invalid (x,y) grid size (" << m_grid_x << ", " << m_grid_y << ")." << std::endl;
		goto _exit;
	}
	// allocation of temp. working arrays
	qx = (float*)malloc(sizeof(float) * m_grid_x);
	qy = (float*)malloc(sizeof(float) * m_grid_y);
	if (m_pot3d > 0 && m_grid_z > 0) {
		qz = (float*)malloc(sizeof(float) * m_grid_z);
	}
	if (0 == m_pot3d && npix > 0) { // 2d projection mode
		cg2d = (fcmplx*)malloc(sizeof(fcmplx) * npix);
		fg2d = (float*)malloc(sizeof(float) * npix * 3);
	}
	// calculations
	// ---------------------------------------------------------------------------
	// - shift helpers
	nx2 = (m_grid_x - m_grid_x % 2) >> 1;
	for (i = 0; i < (int)m_grid_x; i++) {
		fi = (((i + nx2) % (int)m_grid_x) - nx2);
		gx = (float)fi * itogx;
		qx[i] = gx;
	}
	ny2 = (m_grid_y - m_grid_y % 2) >> 1;
	for (j = 0; j < (int)m_grid_y; j++) {
		fj = (((j + ny2) % (int)m_grid_y) - ny2);
		gy = (float)fj * itogy;
		qy[j] = gy;
	}
	if (m_pot3d > 0 && m_grid_z > 0) {
		nz2 = (m_grid_z - m_grid_z % 2) >> 1;
		for (j = 0; j < (int)m_grid_z; j++) {
			fj = (((j + nz2) % (int)m_grid_z) - nz2);
			gy = (float)fj * itogz;
			qz[j] = gy;
		}
	}
	//
	if (m_gpu_id >= 0) { // copy to device buffers
		//cuerr = cudaSetDevice(m_gpu_id); // (set the device again?)
		cuerr = cudaMemcpy(m_d_qx, qx, sizeof(float) * m_grid_x, cudaMemcpyHostToDevice);
		if (cuerr != cudaSuccess) {
			nerr = 103;
			std::cerr << "Error (CJEElaSca::CalculateHelper): Failed to copy x-shift helper to device. " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
		}
		cuerr = cudaMemcpy(m_d_qy, qy, sizeof(float) * m_grid_y, cudaMemcpyHostToDevice);
		if (cuerr != cudaSuccess) {
			nerr = 104;
			std::cerr << "Error (CJEElaSca::CalculateHelper): Failed to copy y-shift helper to device. " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
		}
		if (m_pot3d > 0 && m_grid_z > 0) {
			cuerr = cudaMemcpy(m_d_qz, qz, sizeof(float) * m_grid_z, cudaMemcpyHostToDevice);
			if (cuerr != cudaSuccess) {
				nerr = 105;
				std::cerr << "Error (CJEElaSca::CalculateHelper): Failed to copy z-shift helper to device. " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
			}
		}
	}
	// copy to host buffers (independent of cpu usage, we want this to exist)
	memcpy(m_h_qx, qx, sizeof(float) * m_grid_x);
	memcpy(m_h_qy, qy, sizeof(float) * m_grid_y);
	if (m_pot3d > 0 && m_grid_z > 0) {
		memcpy(m_h_qz, qz, sizeof(float) * m_grid_z);
	}
	
	//
	if (m_pot3d > 0) { // 3d potential mode, no need for 2d helpers
		goto _pot3d;
	}
	//

	// ---------------------------------------------------------------------------
	// * TODO * Investigate an analytical form of an absorptive form factor for the ionic term
	idx = 0;
	apxy1 = 1.; apxy2 = 1.;
	kmax = std::min((float)nx2 * itogx, (float)ny2 * itogy);
	kthr23 = kmax * EELSCA_CONST_PGR_BWL;
	pfacio = GetRelCorr(m_ekv) * (float)(_CFFA) * (float)(_FPI) * 10.f; // prefactor for the ionic rest charge potential : rel-corr * C * 4 Pi * 10.[->nm ^ -1]
	for (j = 0; j < (int)m_grid_y; j++) {
		gy2 = qy[j] * qy[j];
		for (i = 0; i < (int)m_grid_x; i++) {
			g2 = qx[i] * qx[i] + gy2;
			fg2d[idx] = g2; // <--------------------------- g^2 helper
			//
			gm = std::sqrt(g2);
			if (m_bwl_pot) {
				apxy1 = 0.5f - 0.5f * std::tanh((gm / kmax - 0.9f) * 30.f);
			}
			fg2d[idx + npix] = apxy1; // <----------------- potential aperture
			//
			if (m_bwl_pgr) {
				if (gm > kthr23) {
					apxy2 = 0.;
				}
				else {
					apxy2 = ftrn;
				}
			}
			fg2d[idx + 2 * npix] = apxy2; // <------------- phase grating aperture
			//
			// * TODO * : calculate an absorptive form factor for the ionic point-charge potential
			cg2d[idx] = fcmplx(apxy1 * pfacio / (0.25f * g2 + iyr2), 0.f); // <----------- ionic point-charge potential
			idx++;
		}
	}
	if (m_gpu_id >= 0) { // copy to device buffer
		cuerr = cudaMemcpy(m_d_ap, &fg2d[npix], sizeof(float) * npix * 2, cudaMemcpyHostToDevice);
		if (cuerr != cudaSuccess) {
			nerr = 112;
			std::cerr << "Error (CJEElaSca::CalculateHelper): Failed to copy aperture helpers to device. " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
		}
		cuerr = cudaMemcpy(m_d_potio, cg2d, sizeof(fcmplx) * npix, cudaMemcpyHostToDevice);
		if (cuerr != cudaSuccess) {
			nerr = 113;
			std::cerr << "Error (CJEElaSca::CalculateHelper): Failed to copy ionic term helper to device. " << cuerr << ": " << cudaGetErrorString(cuerr) << std::endl;
		}
	}
	// copy to host buffer (independent of cpu usage, we want this to exist)
	memcpy(m_h_q2, fg2d, sizeof(float) * npix);
	memcpy(m_h_ap, &fg2d[npix], sizeof(float) * npix * 2);
	memcpy(m_h_potio, cg2d, sizeof(fcmplx) * npix);
	//
	ndims = 2;
	pdims = new int[ndims];
	pdims[0] = m_grid_x; pdims[1] = m_grid_y;
	goto _cores;

_pot3d: // jump point for 3d potential mode

	ndims = 3;
	pdims = new int[ndims];
	pdims[0] = m_grid_x; pdims[1] = m_grid_y; pdims[2] = m_grid_z;
	goto _cores;

_cores: // jump point for core setup

	if (m_gpu_id >= 0) {
		m_jgpuco.Deinit();
		if (0 < m_jgpuco.Init(ndims, pdims)) {
			nerr = 120;
			std::cerr << "Error (CJEElaSca::InitHelper): Failed to initialize GPU FFT module." << std::endl;
			goto _exit;
		}
	}

	if (m_cpu_num > 0) {
		int ncpu_init = 0;
		for (icpu = 0; icpu < m_cpu_num; icpu++) {
			m_jcpuco[icpu].Deinit();
			if (0 == m_jcpuco[icpu].Init(ndims, pdims)) {
				ncpu_init++;
			}
		}
		if (ncpu_init < m_cpu_num) {
			nerr = 20;
			std::cerr << "Error (CJEElaSca::InitHelper): Failed to initialize CPU FFT module." << std::endl;
			goto _exit;
		}
	}
	
	m_status |= (unsigned int)EELSCA_STATUS_HELPER_SET;

_exit:
	if (NULL != pdims) delete[] pdims;
	if (NULL != qx) free(qx);
	if (NULL != qy) free(qy);
	if (NULL != qz) free(qz);
	if (NULL != cg2d) free(cg2d);
	if (NULL != fg2d) free(fg2d);
	return nerr;
}

int CJEElaSca::PrepareScatteringFunctions(int igpu, int ncpu, bool bverbose)
{
	int nerr = 0, ierr = 0;
	size_t nats = 0, iats = 0;
	size_t npix = 0, idx = 0;
	float q = 0.f, ap = 1.f;
	fcmplx ff(0.f, 0.f);
	fcmplx* ff2d = m_h_att_ff2d;
	float* pot_ap = NULL; // aperture buffer address

	if (m_status < (unsigned int)EELSCA_THRESHOLD_ATOMICFF) {
		nerr = 1;
		std::cerr << "Error (CJEElaSca::PrepareScatteringFunctions): Incomplete setup (" << m_status << ")." << std::endl;
		goto _exit;
	}

	if (bverbose) {
		std::cout << std::endl;
		std::cout << "  Preparing scattering functions ..." << std::endl;
	}

	if (bverbose) {
		std::cout << "  - setting helper data" << std::endl;
	}

	// allocate helpers
	ierr = AllocHelper(igpu, ncpu);
	if (0 < ierr) {
		nerr = 2;
		std::cerr << "Error (CJEElaSca::PrepareScatteringFunctions): failed to allocate helper arrays (" << ierr << ")." << std::endl;
		goto _exit;
	}

	// calculate helper data
	ierr = InitHelper();
	if (0 < ierr) {
		nerr = 3;
		std::cerr << "Error (CJEElaSca::PrepareScatteringFunctions): failed to calculate helper data (" << ierr << ")." << std::endl;
		goto _exit;
	}

	nats = m_vTypes.size();
	npix = (size_t)m_grid_x * m_grid_y;
	pot_ap = m_h_ap; // potential bwl
	if (0 == m_pot3d && m_att_num > 0 && m_att_num == nats && npix > 0) { // projection of full form factors to 2d planes
		if (bverbose) {
			std::cout << "  - calculating projected atomic form-factors" << std::endl;
		}
		for (iats = 0; iats < nats; iats++) { // loop over atom types set in the module
			ff2d = &m_h_att_ff2d[iats * npix]; // address to the 2d form factor helper of this atom type 
			for (idx = 0; idx < npix; idx++) { // loop over grid points
				q = std::sqrt(m_h_q2[idx]); // get q 1/nm from q^2 helper
				ff = m_vTypes[iats].GetFormFactor(q, 1); // get form factor for this grid point with linear q interpolation
				if (m_bwl_pot) ap = pot_ap[idx]; // get potential band width limit
				ff2d[idx] = ff * ap; // store projected potential of the atom type
			}
		}
	}
	else { // calculate a 3-d potential and then project.
		nerr = 1001;
		std::cerr << "Error (CJEElaSca::PrepareScatteringFunctions): 3d potential calculation is not yet implemented." << std::endl;
		goto _exit;
	}

_exit:
	return nerr;
}


int CJEElaSca::CalculateMeanInnerPotential(CStructure* pst, fcmplx& mip)
{
	int nerr = 0;
	mip = fcmplx(0.f, 0.f); // output: mean inner potential (volts)
	
	size_t natt = m_vTypes.size(); // get number of atom types prepared (must correspond to that of the input structure)
	size_t nat = 0; // number of atoms of the input structure
	size_t iats = 0, iat = 0, jat = 0, iatt = 0, nat2 = 0;
	float rel_corr = GetRelCorr(m_ekv); // relativistic correction
	float occ = 1.f; // atom site occupancy
	float vol = 0.f; // super-cell volume (nm^3)
	fcmplx ff_cur(0.f, 0.f); // current atomic form factor (nm)
	fcmplx sf_acc(0.f, 0.f); // accumulates structure factor (nm)

	if (NULL == pst) goto _exit;

	nat = pst->m_vAtoms.size(); // get the number of atoms of the input structure
	vol = pst->GetVolume(); // super-cell volume (nm^3)

	if (m_status < EELSCA_THRESHOLD_POTENTIALS) {
		nerr = 1;
		std::cerr << "Error (CJEElaSca::CalculateMeanInnerPotential): Incomplete setup (" << m_status << ")." << std::endl;
		goto _exit;
	}

	if (natt != pst->m_vTypes.size()) {
		nerr = 2;
		std::cerr << "Error (CJEElaSca::CalculateMeanInnerPotential): Inconsistent list of atom types (" << natt << " != " << pst->m_vTypes.size() << ")." << std::endl;
		goto _exit;
	}

	if (natt == 0 || nat == 0) goto _exit;

	// The procedure is structured in the following order:
	// - loop over atom types
	//   * loop over all atoms and for all atoms of the current type:
	//     + accumulate the dc values of the atomic form factors (ignoring ionic charges and biso)
	// - rescale to volts


	// calculate the projected potential
	for (iatt = 0; iatt < natt; iatt++) { // loop over atom types
		ff_cur = m_vTypes[iatt].GetFormFactor(0.f, 1); // projected atomic form factor for this type (nm)
		for (iat = 0; iat < nat; iat++) { // loop over atoms
			if ((size_t)pst->m_vAtoms[iat].m_nAtomTypeIndex == iatt) { // this is an atom of this type
				occ = pst->m_vAtoms[iat].m_fOccupancy; // atom occupancy
				sf_acc += ff_cur * occ; // add to the structure factor
			}
		}
	}

	mip = sf_acc * EELSCA_CONST_POT_SCA / vol / rel_corr; // + extpot_mean ---> V

_exit:
	return nerr;
}


int CJEElaSca::CalculatePhasegratingCPU(CStructure* pst, fcmplx* pgr, int ithread)
{
	int nerr = 0;
	size_t nats = 0; // number of atom sites of the input structure
	size_t natt = m_vTypes.size(); // get number of atom types prepared (must correspond to that of the input structure)
	size_t nat = 0; // get the number of atoms of the input structure
	size_t nx = (size_t)m_grid_x;
	size_t ny = (size_t)m_grid_y;
	size_t npix = nx * ny; // number of pixels in the (x,y) plane
	size_t ix = 0, iy = 0, idx = 0, iats = 0, iat = 0, jat = 0, iatt = 0, nat2 = 0;
	CV3D pos_cur(0.f, 0.f, 0.f);
	float ux = 0.f, phy = 0.f, pha = 0.f, p_x = 0.f, p_y = 0.f;
	float occ = 0.f, crg = 0.f, biso = 0.f;
	float vol = 0.f; // volume of the structure cell (nm^3)
	float wl = ((float)(_WLELKEV) / sqrt(m_ekv * (m_ekv + 2.f * (float)(_EEL0KEV)))); // electron wavelength in vacuum (nm)
	//float rel_corr = GetRelCorr(m_ekv); // relativistic correction
	CV3D* at_pos = NULL;
	CM3D cell_basis;
	fcmplx* pot_cur = NULL;
	fcmplx* ff_cur = NULL;
	fcmplx cst(0.f, 0.f), csig(0.f, 1.f);
	bool bion = false;
	bool bion_supported = false;

	if (NULL == pst) goto _exit;
	if (NULL == pgr) goto _exit;

	nats = pst->m_vSites.size(); // get number of atom sites of the input structure
	nat = pst->m_vAtoms.size(); // get the number of atoms of the input structure
	vol = pst->GetVolume(); // volume of the structure cell (nm^3)
	cell_basis = pst->m_mCellBasis;
	cell_basis.x = cell_basis.x * pst->m_vCellSize.x;
	cell_basis.y = cell_basis.y * pst->m_vCellSize.y;
	cell_basis.z = cell_basis.z * pst->m_vCellSize.z;

	if (m_status < EELSCA_THRESHOLD_POTENTIALS) {
		nerr = 1;
		std::cerr << "Error (CJEElaSca::CalculatePhasegratingCPU): Incomplete setup (" << m_status << ")." << std::endl;
		goto _exit;
	}

	if (natt != pst->m_vTypes.size()) {
		nerr = 2;
		std::cerr << "Error (CJEElaSca::CalculatePhasegratingCPU): Inconsistent list of atom types (" << natt << " != " << pst->m_vTypes.size() << ")." << std::endl;
		goto _exit;
	}

	if (ithread < 0 || ithread >= m_cpu_num) {
		nerr = 3;
		std::cerr << "Error (CJEElaSca::CalculatePhasegratingCPU): Invalid CPU thread id (" << ithread << " / " << m_cpu_num  << ")." << std::endl;
		goto _exit;
	}

	pot_cur = m_jcpuco[ithread].GetData(); // link to the thread fft core buffer and use it as potential buffer
	m_jcpuco[ithread].Zero(); // reset the potential to zero

	if (nats == 0 || natt == 0 || nat == 0) goto _phase_grating;

	// The following procedure is structured in the following order:
	// - loop over atom sites
	//   * update positions for all atoms linked to sites
	// - loop over atom types
	//   * loop over all atoms and for all atoms of the current type:
	//     + calculate the translational phase factor from the linked atom site position
	//     + multiply the 2d form factor (prepared before) with the 2d translational phase plate
	//     + add the result to the potential

	// allocate memory for a list of positions for all atoms
	at_pos = (CV3D*)malloc(nat * sizeof(CV3D));
	if (NULL == at_pos) {
		nerr = 4;
		std::cerr << "Error (CJEElaSca::CalculatePhasegratingCPU): position buffer allocation failed." << std::endl;
		goto _exit;
	}
	// preset positions with atom positions
	for (iat = 0; iat < nat; iat++) {
		at_pos[iat] = cell_basis.VProduct(pst->m_vAtoms[iat].m_vPos);
	}

	// set atom positions to positions of atomic sites
	// * for calculations in the QEP model, move atoms linked to one site coherently
	for (iats = 0; iats < nats; iats++) { // loop over atom sites
		pos_cur = cell_basis.VProduct(pst->m_vSites[iats].m_vPos); // physical position of the site in the structure cell
		if (mod_thermal == 2) { // QEP calculation: modify the site position according to thermal vibration amplitude
			ux = sqrt(pst->m_vSites[iats].m_fBiso) * EELSCA_CONST_RSQR8PI2;
			pos_cur.x += ux * NormRand();
			pos_cur.y += ux * NormRand();
		}
		nat2 = pst->m_vSites[iats].m_vatocc.size(); // number of linked atoms
		if (nat2 > 0) {
			for (jat = 0; jat < nat2; jat++) { // loop over all atoms linked to the current site
				iat = pst->m_vSites[iats].m_vatocc[jat];
				at_pos[iat] = pos_cur;  // store updated site position per linked atom
			}
		}
	}

	// calculate the projected potential
	for (iatt = 0; iatt < natt; iatt++) { // loop over atom types
		crg = m_vTypes[iatt].m_fIonicCharge; // ionic charge
		biso = m_vTypes[iatt].m_fBiso; // biso parameters
		ff_cur = &m_h_att_ff2d[iatt * npix]; // address of the projected atomic form factor for this type
		if (std::abs(crg) > 0.001f && bion_supported) { bion = true; } // switch ionic term on
		for (iat = 0; iat < nat; iat++) { // loop over atoms
			if ((size_t)pst->m_vAtoms[iat].m_nAtomTypeIndex == iatt) { // this is an atom of this type
				occ = pst->m_vAtoms[iat].m_fOccupancy; // atom occupancy
				idx = 0;
				for (iy = 0; iy < ny; iy++) {
					phy = at_pos[iat].y * m_h_qy[iy]; // y part of the shifting phase
					for (ix = 0; ix < nx; ix++) {
						pha = (float)_TPI * (phy + at_pos[iat].x * m_h_qx[ix]); // total phase shift
						cst.real(std::cos(pha)); cst.imag(-std::sin(pha)); // complex phase factor shifting the atomic form factor exp[-I 2*Pi*(x*qx + y*qy)]
						if (bion) {
							if (UseDWF()) { // neutral potential + dampened ionic potential 
								pot_cur[idx] += (ff_cur[idx] + m_h_potio[idx] * crg * GetDWF(biso, m_h_q2[idx])) * cst * m_h_ap[idx] * occ;
							}
							else { // neutral & ionic potential without DWF
								pot_cur[idx] += (ff_cur[idx] + m_h_potio[idx] * crg) * cst * m_h_ap[idx] * occ;
							}
						}
						else { // only neutral atom form factor is added
							// potential += form factor * shift factor * potential bwl * occupancy
							pot_cur[idx] += ff_cur[idx] * cst * m_h_ap[idx] * occ;
						}
						idx++; // increment pixel index
					} // loop columns
				} // loop rows
			}
		} // loop atoms
	} // loop atom types

	// to real space
	m_jcpuco[ithread].IFT();
	// * TODO * //
	// Check whether there is a re-normalization required after the IFT. (Usually 1/(nx*ny) or 1/sqrt(nx*ny)).
	m_jcpuco[ithread].Scale(EELSCA_CONST_POT_SCA / vol);

	// * TODO * //
	// Implement that an external projected potential can be added ! Relativistic correction needs to be applied to the external potential!
	// m_jcpuco[ithread].AddC(pot_ext);
	// Implement that the projected potential may be saved (remove relativistic correction before doing that).

_phase_grating:

	// Remark:
	// Simplification possible, since EELSCA_CONST_POT_SCA * EELSCA_CONST_PGR_SIG = 1 / (4 * PI).
	// However, the two constants are used to obtain the potential in volts at this stage (rel. corr. factor still included).

	csig = fcmplx(0.f, EELSCA_CONST_PGR_SIG * pst->m_vCellSize.z * wl);
	for (idx = 0; idx < npix; idx++) {
		pot_cur[idx] = std::exp(pot_cur[idx] * csig); // phase = exp( i Sigma V )
	}

	if (m_bwl_pgr) {
		m_jcpuco[ithread].FT();
		m_jcpuco[ithread].MultiplyF(&m_h_ap[npix]); // The aperture is precalculated with a renormalizion to account for FT summations.
		m_jcpuco[ithread].IFT();
	}

	m_jcpuco[ithread].GetDataC(pgr);

_exit:
	if (NULL != at_pos) { free(at_pos); }
	return nerr;
}

int CJEElaSca::CalculatePhasegratingGPU(CStructure* pst, fcmplx* pgr)
{
	int nerr = 0;
	cudaError cuerr;
	size_t nats = 0; // number of atom sites of the input structure
	size_t natt = m_vTypes.size(); // get number of atom types prepared (must correspond to that of the input structure)
	size_t nat = 0; // get the number of atoms of the input structure
	size_t nx = (size_t)m_grid_x;
	size_t ny = (size_t)m_grid_y;
	size_t npix = nx * ny; // number of pixels in the (x,y) plane
	size_t ix = 0, iy = 0, idx = 0, iats = 0, iat = 0, jat = 0, iatt = 0, nat2 = 0;
	CV3D pos_cur(0.f, 0.f, 0.f);
	float ux = 0.f, phy = 0.f, pha = 0.f, p_x = 0.f, p_y = 0.f;
	float occ = 0.f, crg = 0.f, biso = 0.f;
	float vol = 0.f; // volume of the structure cell (nm^3)
	float wl = ((float)(_WLELKEV) / sqrt(m_ekv * (m_ekv + 2.f * (float)(_EEL0KEV)))); // electron wavelength in vacuum (nm)
	//float rel_corr = GetRelCorr(m_ekv); // relativistic correction
	CV3D* at_pos = NULL;
	CM3D cell_basis;
	cuComplex* pot_cur = NULL;
	fcmplx cst(0.f, 0.f), csig(0.f, 1.f);
	bool bion = false;
	bool bion_supported = false;

	if (NULL == pst) goto _exit;
	if (NULL == pgr) goto _exit;

	nats = pst->m_vSites.size(); // get number of atom sites of the input structure
	nat = pst->m_vAtoms.size(); // get the number of atoms of the input structure
	vol = pst->GetVolume(); // volume of the structure cell (nm^3)
	cell_basis = pst->m_mCellBasis;
	cell_basis.x = cell_basis.x * pst->m_vCellSize.x;
	cell_basis.y = cell_basis.y * pst->m_vCellSize.y;
	cell_basis.z = cell_basis.z * pst->m_vCellSize.z;

	if (m_status < EELSCA_THRESHOLD_POTENTIALS) {
		nerr = 1;
		std::cerr << "Error (CJEElaSca::CalculatePhasegratingGPU): Incomplete setup (" << m_status << ")." << std::endl;
		goto _exit;
	}

	if (natt != pst->m_vTypes.size()) {
		nerr = 2;
		std::cerr << "Error (CJEElaSca::CalculatePhasegratingGPU): Inconsistent list of atom types (" << natt << " != " << pst->m_vTypes.size() << ")." << std::endl;
		goto _exit;
	}

	pot_cur = m_jgpuco.GetData(); // link to the thread fft core buffer and use it as potential buffer
	m_jgpuco.Zero(); // reset the potential to zero

	if (nats == 0 || natt == 0 || nat == 0) goto _phase_grating;

	// The following procedure is structured in the following order:
	// - loop over atom sites
	//   * update positions for all atoms linked to sites
	// - loop over atom types
	//   * loop over all atoms and for all atoms of the current type:
	//     + calculate the translational phase factor from the linked atom site position
	//     + multiply the 2d form factor (prepared before) with the 2d translational phase plate
	//     + add the result to the potential

	// allocate memory for a list of positions for all atoms
	at_pos = (CV3D*)malloc(nat * sizeof(CV3D));
	if (NULL == at_pos) {
		nerr = 4;
		std::cerr << "Error (CJEElaSca::CalculatePhasegratingGPU): position buffer allocation failed." << std::endl;
		goto _exit;
	}
	// preset positions with atom positions
	for (iat = 0; iat < nat; iat++) {
		at_pos[iat] = cell_basis.VProduct(pst->m_vAtoms[iat].m_vPos);
	}

	// set atom positions to positions of atomic sites
	// * for calculations in the QEP model, move atoms linked to one site coherently
	for (iats = 0; iats < nats; iats++) { // loop over atom sites
		pos_cur = cell_basis.VProduct(pst->m_vSites[iats].m_vPos); // physical position of the site in the structure cell
		if (mod_thermal == 2) { // QEP calculation: modify the site position according to thermal vibration amplitude
			ux = sqrt(pst->m_vSites[iats].m_fBiso) * EELSCA_CONST_RSQR8PI2;
			pos_cur.x += ux * NormRand();
			pos_cur.y += ux * NormRand();
		}
		nat2 = pst->m_vSites[iats].m_vatocc.size(); // number of linked atoms
		if (nat2 > 0) {
			for (jat = 0; jat < nat2; jat++) { // loop over all atoms linked to the current site
				iat = pst->m_vSites[iats].m_vatocc[jat];
				at_pos[iat] = pos_cur;  // store updated site position per linked atom
			}
		}
	}

	// calculate the projected potential
	for (iatt = 0; iatt < natt; iatt++) { // loop over atom types
		crg = m_vTypes[iatt].m_fIonicCharge; // ionic charge
		biso = m_vTypes[iatt].m_fBiso; // biso parameters
		cuerr = cudaMemcpy(m_d_ff2d, &m_h_att_ff2d[iatt * npix], sizeof(cuComplex) * npix, cudaMemcpyHostToDevice);
		if (cuerr != cudaSuccess) {
			nerr = 5;
			std::cerr << "Error (CJEElaSca::CalculatePhasegratingGPU): Failed to copy 2d form factor of atom type " << iatt + 1 << " of " << natt << "." << std::endl;
			goto _exit;
		}
		if (std::abs(crg) > 0.001f && bion_supported) { bion = true; } // switch ionic term on
		for (iat = 0; iat < nat; iat++) { // loop over atoms
			if ((size_t)pst->m_vAtoms[iat].m_nAtomTypeIndex == iatt) { // this is an atom of this type
				occ = pst->m_vAtoms[iat].m_fOccupancy; // atom occupancy
				idx = 0;
				if (bion) {
					if (UseDWF()) { // neutral potential + dampened ionic potential
						cuerr = AddDampedIonFormFactor2d(pot_cur, (cuComplex*)m_d_ff2d, (cuComplex*)m_d_potio, m_d_qx, m_d_qy, m_d_ap, at_pos[iat].x, at_pos[iat].y, crg, biso, occ, m_grid_x, m_grid_y);
					}
					else { // neutral potential + ionic potential without DWF
						cuerr = AddIonFormFactor2d(pot_cur, (cuComplex*)m_d_ff2d, (cuComplex*)m_d_potio, m_d_qx, m_d_qy, m_d_ap, at_pos[iat].x, at_pos[iat].y, crg, occ, m_grid_x, m_grid_y);
					}
				}
				else { // only neutral atom form factor is added
					// potential += form factor * shift factor * potential bwl * occupancy
					cuerr = AddFormFactor2d(pot_cur, (cuComplex*)m_d_ff2d, m_d_qx, m_d_qy, m_d_ap, at_pos[iat].x, at_pos[iat].y, occ, m_grid_x, m_grid_y);
				}
				if (cuerr != cudaSuccess) {
					nerr = 6;
					std::cerr << "Error (CJEElaSca::CalculatePhasegratingGPU): Failed to add 2d form factor of atom " << iat + 1 << " of " << nat << "." << std::endl;
					goto _exit;
				}
			}
		} // loop atoms
	} // loop atom types

	// to real space
	m_jgpuco.IFT();
	// * TODO * //
	// Check whether there is a re-normalization required after the IFT. (Usually 1/(nx*ny) or 1/sqrt(nx*ny)).
	m_jgpuco.Scale(EELSCA_CONST_POT_SCA / vol);

	// * TODO * //
	// Implement that an external projected potential can be added ! Relativistic correction needs to be applied to the external potential!
	// m_jcpuco[ithread].AddC(pot_ext);
	// Implement that the projected potential may be saved (remove relativistic correction before doing that).

_phase_grating:

	// Remark:
	// Simplification possible, since EELSCA_CONST_POT_SCA * EELSCA_CONST_PGR_SIG = 1 / (4 * PI).
	// However, the two constants are used to obtain the potential in volts at this stage (rel. corr. factor still included).

	csig = fcmplx(0.f, EELSCA_CONST_PGR_SIG * pst->m_vCellSize.z * wl);
	cuerr = PotToPgr(pot_cur, pot_cur, (cuComplex*)(&csig), m_grid_x, m_grid_y); // phase = exp( i Sigma V )
	if (cuerr != cudaSuccess) {
		nerr = 7;
		std::cerr << "Error (CJEElaSca::CalculatePhasegratingGPU): Transform the projected potential to a transmission function." << std::endl;
		goto _exit;
	}
	
	if (m_bwl_pgr) {
		m_jgpuco.FT();
		m_jgpuco.MultiplyF(&m_h_ap[npix]); // The aperture is precalculated with a renormalizion to account for FT summations.
		m_jgpuco.IFT();
	}

	m_jgpuco.GetDataC(pgr);

_exit:
	return nerr;
}