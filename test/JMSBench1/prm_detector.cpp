// file: 'prm_detector.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementation of the class prm_detector.
//
// Class prm_detector handles parameters and setup routines around the
// detection of signal for TEM image simulations.
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

#include "prm_detector.h"
#include "string_format.h"

/******************** --- prm_annular --- *******************************/

prm_annular::prm_annular()
{
	name = "";
	beta_inner = 0.f;
	beta_outer = 10.f;
	type = (unsigned int)ANNULAR_RING;
	phi_begin = 0.f;
	phi_end = 0.f;
}

prm_annular::~prm_annular()
{
}

void prm_annular::reset()
{
	name = "";
	beta_inner = 0.f;
	beta_outer = 10.f;
	type = (unsigned int)ANNULAR_RING;
	phi_begin = 0.f;
	phi_end = 0.f;
}

std::string prm_annular::get_str_info(bool header)
{
	std::string stmp = "";

	if (header) {
		stmp = "scatt. angles [mrad] | azimuth range [deg] | name";
	}
	else {
		if (type == ANNULAR_RING) {
			stmp = format(" %8.2f - %8.2f |                full | %s", beta_inner, beta_outer, name.c_str());
		}

		if (type & ANNULAR_SEGMENTED) {
			stmp = format(" %8.2f - %8.2f | %8.2f - %8.2f | %s", beta_inner, beta_outer, phi_begin, phi_end, name.c_str());
		}
	}

	return stmp;
}

std::string prm_annular::get_str_params()
{
	if (type == ANNULAR_RING) {
		return format("%g %g 0.0 0.0 %s", beta_inner, beta_outer, name.c_str());
	}
	return format("%g %g %g %g %s", beta_inner, beta_outer, phi_begin, phi_end, name.c_str());
}

int prm_annular::setup_data(std::string str_ctrl)
{
	if (binteractive) {
		int iseg = 0;
		reset();
		std::cout << std::endl;
		std::cout << "  Enter the inner detection angle (mrad): ";
		std::cin >> beta_inner;
		std::cout << std::endl;
		std::cout << "  Enter the outer detection angle (mrad): ";
		std::cin >> beta_outer;
		std::cout << std::endl;
		std::cout << "  Is this an azimuthal detector degment:  <0> No  <1> Yes ? ";
		std::cin >> iseg;
		if (iseg != 0) {
			type ^= (unsigned int)ANNULAR_SEGMENTED;
			std::cout << std::endl;
			std::cout << "  Enter the segment azimuth begin (deg): ";
			std::cin >> phi_begin;
			std::cout << std::endl;
			std::cout << "  Enter the segment azimuth end   (deg): ";
			std::cin >> phi_end;
		}
		std::cout << std::endl;
		std::cout << "  Enter a detector name: ";
		std::cin >> name;
	}
	else {
		if (str_ctrl.size() > 0) {
			std::string sprm = "";
			reset();
			int ipos = 0;
			ipos = read_param(ipos, &str_ctrl, &sprm);
			if (ipos < 0) { return 50; } // parsing error, no parameters found
			beta_inner = to_float(sprm);
			ipos = read_param(ipos, &str_ctrl, &sprm);
			if (ipos < 0) { return 51; } // parsing error, no outer angle found
			beta_outer = to_float(sprm);
			ipos = read_param(ipos, &str_ctrl, &sprm);
			if (ipos < 0) { return 52; } // parsing error, no azimuth begin found
			phi_begin = to_float(sprm);
			ipos = read_param(ipos, &str_ctrl, &sprm);
			if (ipos < 0) { return 53; } // parsing error, no azimuth end found
			phi_end = to_float(sprm);
			ipos = read_param(ipos, &str_ctrl, &sprm);
			if (ipos < 0) { return 54; } // parsing error, no name found
			name = to_string(sprm);
			//
			// analyze azimuth range for segmented flag
			if (phi_begin + 0.001f < phi_end) {
				type ^= (unsigned int)ANNULAR_SEGMENTED;
			}
		}
		else { // empty setup string
			return 99;
		}
	}
	return 0;
}




/******************** --- prm_detector --- ******************************/

prm_detector::prm_detector()
{
	b_annular = false;
	b_difpat = false;
	b_difpat_avg = false;
	b_image = false;
	b_wave = false;
	b_waveft = false;
	b_wave_avg = false;
	b_waveft_avg = false;
}

prm_detector::~prm_detector()
{
	v_annular.clear();
}

unsigned int prm_detector::get_num_det_annular()
{
	num_det_annular = (unsigned int)v_annular.size();
	return num_det_annular;
}

void prm_detector::print_setup_annular(std::vector<prm_annular> *pv_ad)
{
	int ndet = 0, i = 0;
	prm_annular adtmp;
	std::vector<prm_annular> *pv_ad_use = &v_annular;
	if (NULL != pv_ad) pv_ad_use = pv_ad;
	ndet = (int)pv_ad_use->size();
	std::cout << std::endl;
	std::cout << "  Annular detector setup (" << ndet << " detectors):" << std::endl;
	if (ndet > 0) {
		std::cout << "  ID | " << adtmp.get_str_info(true) << std::endl;
		for (i = 0; i < ndet; i++) {
			std::cout << format("  %2d | ", i) << pv_ad_use->at(i).get_str_info() << std::endl;
		}
	}
	return;
}

int prm_detector::setup_annular()
{
	int ierr = 0, ichk = 0, i = 0, idet = 0, j = 0;
	std::string  stmp, sprm;
	prm_annular adtmp;
	std::vector<prm_annular> v_ad_tmp;
	bool bsegmented = false;
	int ndet = 0;

	// initialize local control interface
	adtmp.btalk = btalk;
	adtmp.binteractive = binteractive;
	adtmp.ndebug = ndebug;
	v_ad_tmp = v_annular;

	if (binteractive) {
	_repeat_input: // display the current setup of annular detectors
		ndet = (int)v_ad_tmp.size();
		print_setup_annular(&v_ad_tmp); // print the current temporary setup
		// allow to modify the setup
		std::cout << std::endl;
		std::cout << "  <0> Accept  <1> Add  <2> Delete  <3> Change ? ";
		std::cin >> ichk;
		switch (ichk) {
		case 0: // jump out and use this list
			break;
		case 1: // add another detector
			adtmp.setup_data();
			v_ad_tmp.push_back(adtmp);
			goto _repeat_input;
			break;
		case 2: // delete
			std::cout << std::endl;
			std::cout << "  Select the ID of the detector to delete from the current setup: ";
			std::cin >> idet;
			if (idet >= 0 && idet < ndet) {
				v_ad_tmp.erase(v_ad_tmp.begin() + idet);
			}
			goto _repeat_input;
			break;
		case 3: // change
			std::cout << std::endl;
			std::cout << "  Select the ID of the detector to change in the current setup: ";
			std::cin >> idet;
			if (idet >= 0 && idet < ndet) {
				v_ad_tmp[idet].setup_data();
			}
			goto _repeat_input;
			break;
		default:
			goto _repeat_input;
			break;
		}
		// store setup in the object member
		v_annular = v_ad_tmp;
		num_det_annular = (unsigned int)v_annular.size();
		ndet = (int)num_det_annular;
		if (ndet > 0) { // write the setup to the list of control lines
			for (i = 0; i < ndet; i++) {
				v_str_ctrl.push_back("set_det_annular");
				v_str_ctrl.push_back(v_annular[i].get_str_params());
			}
		}
	}
	else { // read the setup from the control lines
		v_annular.clear(); i = 0; j = 0;
		while (i >= 0) {
			i = ctrl_find_param("set_det_annular", &stmp, j);
			if (i >= 0) {
				ierr = adtmp.setup_data(stmp);
				if (0 == ierr) {
					v_annular.push_back(adtmp);
				}
				j = i + 2;
			}
		}
		print_setup_annular(); // print the current temporary setup
	}

	return ierr;
}
