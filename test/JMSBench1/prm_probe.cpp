// file: 'prm_probe.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementation of the class prm_probe.
//
// Class prm_probe handles parameters and setup routines around the
// TEM probe / beam and lens parameters.
//
/* -----------------------------------------------------------------------

	This file is part of JMSBench1.

	JMSBench1 is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	JMSBench1 is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with JMSBench1.  If not, see <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------- */

#include "prm_probe.h"
#include "NatureConstants.h" // from JMS


prm_probe::prm_probe()
{
	ekv = 200.f;
	alpha = 25.f;
	aberrations = NULL;
	focus_spread = 5.f;
	size = 0.05f;
	size_distribution = 1;
	num_aberrations = 0;
}

prm_probe::~prm_probe()
{
	if (NULL != aberrations) free(aberrations);
}


unsigned int prm_probe::prepare_aberrations(unsigned int num)
{
	unsigned int nready = 0;
	if (NULL != aberrations) {
		free(aberrations);
		aberrations = NULL;
		num_aberrations = 0;
	}
	if (num > 0) {
		size_t sz_list = sizeof(float) * (size_t)(2 * num);
		aberrations = (float*)malloc(sz_list);
		if (NULL != aberrations) {
			nready = num;
			memset(aberrations, 0, sz_list);
		}
	}
	return nready;
}


void prm_probe::copy_data_from(prm_probe *psrc)
{
	if (psrc) {
		ekv = psrc->ekv;
		alpha = psrc->alpha;
		focus_spread = psrc->focus_spread;
		size = psrc->size;
		size_distribution = psrc->size_distribution;

		prepare_aberrations(psrc->get_num_aberrations());
		if (num_aberrations > 0 && NULL != aberrations) {
			psrc->copy_aberrations_to(aberrations);
		}
	}
}


int prm_probe::copy_aberrations_to(float *buf)
{
	int nerr = 0;
	size_t sz_list = sizeof(float) * 2 * num_aberrations;
	if (NULL == buf) {
		return 1;
	}
	if (sz_list > 0) {
		memcpy(buf, aberrations, sz_list);
	}
	return nerr;
}


unsigned int prm_probe::get_num_aberrations()
{
	return num_aberrations;
}


int prm_probe::get_aberration_coeff(unsigned int idx, float &x, float &y)
{
	if (idx < num_aberrations) {
		x = aberrations[2 * idx];
		y = aberrations[1 + 2 * idx];
	}
	else {
		return 1; // error: invalid index
	}
	return 0;
}


int prm_probe::set_aberration_coeff(unsigned int idx, float x, float y)
{
	if (idx < num_aberrations) {
		aberrations[2 * idx] = x;
		aberrations[1 + 2 * idx] = y;
	}
	else {
		return 1; // error: invalid index
	}
	return 0;
}


float prm_probe::get_wl()
{
	return (float)(_WLELKEV / sqrt(ekv*(ekv + 2.*_EEL0KEV)));
}


// user input and setup of electron energy (high tension)
int prm_probe::setup_ht(void)
{
	int ierr = 0, i = 0;
	std::string  stmp, sprm;
	
	if (binteractive) {
_repeat_input:
		std::cout << std::endl;
		std::cout << "  Set the probe electron energy in keV: ";
		std::cin >> ekv;
		if (ekv < 10.f || ekv > 10000.f) {
			std::cout << std::endl;
			std::cerr << "Error: Invalid probe electron energy, accepted range: 10 keV - 10000 keV.";
			goto _repeat_input;
		}
		v_str_ctrl.push_back("set_probe_ekv");
		v_str_ctrl.push_back(format("%8.2f", ekv));
	}
	else {
		i = ctrl_find_param("set_probe_ekv", &stmp);
		if (i >= 0) {
			read_param(0, &stmp, &sprm);
			ekv = to_float(sprm);
		}
	}

	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  - electron wavelength in vacuum: " << format("%6.4f", get_wl()*1000.f) << " pm" << std::endl;
	}

	return ierr;
}

// user input and setup of the probe-forming aperture
int prm_probe::setup_aperture(void)
{
	int ierr = 0, i = 0;
	std::string  stmp, sprm;

	if (binteractive) {
	_repeat_input:
		std::cout << std::endl;
		std::cout << "  Set the convergence semi-angle of the electron probe in mrad: ";
		std::cin >> alpha;
		if (alpha < 0.001f || alpha > 200.f) {
			std::cout << std::endl;
			std::cerr << "Error: Invalid probe convergence, accepted range: 0.001 mrad - 200 mrad.";
			goto _repeat_input;
		}
		v_str_ctrl.push_back("set_probe_convergence");
		v_str_ctrl.push_back(format("%8.2f", alpha));
	}
	else {
		i = ctrl_find_param("set_probe_convergence", &stmp);
		if (i >= 0) {
			read_param(0, &stmp, &sprm);
			alpha = to_float(sprm);
		}
	}

	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  - probe max. diffraction vector: " << format("%6.4f", alpha/get_wl()*0.001f) << " nm-1" << std::endl;
	}

	return ierr;
}

// user input and setup of aberrations
int prm_probe::setup_aberrations(void)
{
	int ierr = 0;

	// TODO: implement aberration setup

	return ierr;
}

int prm_probe::get_jpg_params(CJProbeParams *pjpp)
{
	int ierr = 0, nab = 0;
	if (NULL == pjpp) {
		ierr = 1;
		std::cerr << "Error: (set_jpg_params) invalid argument." << std::endl;
		goto _exit;
	}
	pjpp->m_wl = get_wl();
	pjpp->m_alpha = alpha;
	pjpp->m_alpha_adir = 0.f;
	pjpp->m_alpha_asym = 0.f;
	pjpp->m_alpha_rs = 0.03f;
	pjpp->m_alpha_x0 = 0.f;
	pjpp->m_alpha_y0 = 0.f;
	pjpp->m_btx = 0.f;
	pjpp->m_bty = 0.f;
	pjpp->m_fspread_kernel_samples = 7;
	pjpp->m_fspread_kernel_width = 2.f;
	pjpp->m_fspread_width = focus_spread;
	pjpp->m_source_shape = (int)size_distribution;
	pjpp->m_source_width = size;
	nab = __min((int)num_aberrations, pjpp->GetAberrationNum());
	if (nab > 0) {
		for (int i = 0; i < nab; i++) {
			pjpp->m_abrr_coeff[2 * i] = aberrations[2 * i];
			pjpp->m_abrr_coeff[1 + 2 * i] = aberrations[1 + 2 * i];
		}
	}

_exit:
	return ierr;
}