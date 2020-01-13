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
#include "string_format.h"
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


unsigned int prm_probe::get_num_aberrations()
{
	return num_aberrations;
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

	// to be implemented later

	return ierr;
}
