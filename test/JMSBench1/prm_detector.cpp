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
#include "JMultiSlice.h"

/******************** --- prm_annular --- *******************************/

prm_annular::prm_annular()
{
	name = "";
	beta_inner = 0.f;
	beta_outer = 10.f;
	type = (unsigned int)ANNULAR_RING;
	phi_begin = 0.f;
	phi_end = 0.f;
	grid_nx = 0;
	grid_ny = 0;
	msklen = 0;
	msk = NULL;
	det = NULL;
}

prm_annular::prm_annular(const prm_annular &src)
{
	name = src.name;
	beta_inner = src.beta_inner;
	beta_outer = src.beta_outer;
	type = src.type;
	phi_begin = src.phi_begin;
	phi_end = src.phi_end;
	grid_nx = 0;
	grid_ny = 0;
	msklen = 0;
	msk = NULL;
	det = NULL;
	if (src.grid_nx > 0 && src.grid_ny > 0 && src.det != NULL) {
		size_t sz_det = sizeof(float) * src.grid_nx * src.grid_ny;
		det = (float*)malloc(sz_det);
		memcpy(det, src.det, sz_det);
		grid_nx = src.grid_nx;
		grid_ny = src.grid_ny;
		if (msk != NULL) {
			sz_det = sizeof(int) * src.grid_nx * src.grid_ny;
			msk = (int*)malloc(sz_det);
			memcpy(msk, src.msk, sz_det);
			msklen = src.msklen;
		}
	}
}


prm_annular::~prm_annular()
{
	if (NULL != msk) { free(msk); }
	if (NULL != det) { free(det); }
}


prm_annular& prm_annular::operator= (const prm_annular &src)
{
	if (this != &src) {
		name = src.name;
		beta_inner = src.beta_inner;
		beta_outer = src.beta_outer;
		type = src.type;
		phi_begin = src.phi_begin;
		phi_end = src.phi_end;
		grid_nx = 0;
		grid_ny = 0;
		msklen = 0;
		msk = NULL;
		det = NULL;
		if (src.grid_nx > 0 && src.grid_ny > 0 && src.det != NULL) {
			size_t sz_det = sizeof(float) * src.grid_nx * src.grid_ny;
			det = (float*)malloc(sz_det);
			memcpy(det, src.det, sz_det);
			grid_nx = src.grid_nx;
			grid_ny = src.grid_ny;
			if (msk != NULL) {
				sz_det = sizeof(int) * src.grid_nx * src.grid_ny;
				msk = (int*)malloc(sz_det);
				memcpy(msk, src.msk, sz_det);
				msklen = src.msklen;
			}
		}
	}
	return *this;
}

void prm_annular::reset()
{
	name = "";
	beta_inner = 0.f;
	beta_outer = 10.f;
	type = (unsigned int)ANNULAR_RING;
	phi_begin = 0.f;
	phi_end = 0.f;
	grid_nx = 0;
	grid_ny = 0;
	msklen = 0;
	if (NULL != msk) { free(msk); msk = NULL; }
	if (NULL != det) { free(det); det = NULL; }
}


void prm_annular::copy_setup_from(prm_annular *psrc)
{
	if (NULL != psrc) {
		reset();
		name = psrc->name;
		beta_inner = psrc->beta_inner;
		beta_outer = psrc->beta_outer;
		type = psrc->type;
		phi_begin = psrc->phi_begin;
		phi_end = psrc->phi_end;
	}
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

int prm_annular::calculate_func_jms(CJMultiSlice *pjms, bool b_create_mask)
{
	int nerr = 0;
	int nx = 0, ny = 0;
	size_t nitems = 0;

	if (NULL == pjms) {
		nerr = 1;
		std::cerr << "Error: (calculate_func_jms) invalid pointer to CJMultiSlice object." << std::endl;
		goto _exit;
	}

	pjms->GetGridSize(nx, ny);
	nitems = (size_t)nx * ny;

	if (nitems > 0) {

		if (NULL != msk) { free(msk); msk = NULL; }
		if (NULL != det) { free(det); det = NULL; }
		msklen = 0;

		det = (float*)malloc(sizeof(float) * nitems);

		if (det == NULL) {
			nerr = 2;
			std::cerr << "Error: (calculate_func_jms) failed to allocate memory for detector function." << std::endl;
			goto _exit;
		}

		grid_nx = nx;
		grid_ny = ny;

		if (b_create_mask) {
			msk = (int*)malloc(sizeof(int) * nitems);
			if (msk == NULL) {
				nerr = 3;
				std::cerr << "Error: (calculate_func_jms) failed to allocate memory for detector access mask." << std::endl;
				goto _exit;
			}
			msklen = 0;
		}
		nerr = pjms->CalculateRingDetector(beta_inner, beta_outer, phi_begin, phi_end, 0.f, 0.f, sens_profile_file, det, msklen, msk);
	}

_exit:
	return nerr;
}


int prm_annular::set_func_jms(CJMultiSlice *pjms, int whichcode, int idx)
{
	int nerr = 0;

	if (NULL == pjms) {
		nerr = 1;
		std::cerr << "Error: (set_func_jms) invalid pointer to CJMultiSlice object." << std::endl;
		goto _exit;
	}

	if (idx < 0 || idx >= pjms->GetDetNum()) {
		nerr = 2;
		std::cerr << "Error: (set_func_jms) called with invalid detector index." << std::endl;
		goto _exit;
	}
	
	if (NULL == det) {
		nerr = 3;
		std::cerr << "Error: (set_func_jms) invalid pointer to detector function array." << std::endl;
		goto _exit;
	}

	nerr = pjms->SetDetectorData(whichcode, idx, det, msklen, msk);

_exit:
	return nerr;
}




/******************** --- prm_detector --- ******************************/

prm_detector::prm_detector()
{
	b_annular = false;
	b_difpat = false;
	b_difpat_avg = false;
	b_image = false;
	b_image_avg = false;
	b_wave = false;
	b_waveft = false;
	b_wave_avg = false;
	b_waveft_avg = false;
}

prm_detector::~prm_detector()
{
	v_annular.clear();
}

void prm_detector::copy_data_from(prm_detector *psrc)
{
	if (psrc) {
		b_annular = psrc->b_annular;
		b_difpat = psrc->b_difpat;
		b_difpat_avg = psrc->b_difpat_avg;
		b_image = psrc->b_image;
		b_image_avg = psrc->b_image_avg;
		b_wave = psrc->b_wave;
		b_waveft = psrc->b_waveft;
		b_wave_avg = psrc->b_wave_avg;
		b_waveft_avg = psrc->b_waveft_avg;
		v_annular = psrc->v_annular;
	}
}

void prm_detector::copy_setup_from(prm_detector *psrc)
{
	if (psrc) {
		b_annular = psrc->b_annular;
		b_difpat = psrc->b_difpat;
		b_difpat_avg = psrc->b_difpat_avg;
		b_image = psrc->b_image;
		b_image_avg = psrc->b_image_avg;
		b_wave = psrc->b_wave;
		b_waveft = psrc->b_waveft;
		b_wave_avg = psrc->b_wave_avg;
		b_waveft_avg = psrc->b_waveft_avg;
		size_t num_ann = psrc->v_annular.size();
		v_annular.clear();
		if (num_ann > 0) {
			prm_annular ann_tmp;
			for (int i = 0; i < (int)num_ann; i++) {
				ann_tmp.copy_setup_from(&psrc->v_annular[i]);
				v_annular.push_back(ann_tmp);
			}
		}
	}
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
	adtmp.set_ctrl(*this);
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
		ndet = (int)v_annular.size();
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
		if (btalk) {
			print_setup_annular(); // print the current temporary setup
		}
	}

	if (v_annular.size() > 0) {
		b_annular = true;
	}

	return ierr;
}


void prm_detector::print_setup_pixelated()
{
	std::vector<std::string> v_str_switch;
	if (b_difpat) {
		v_str_switch.push_back("on");
	}
	else {
		v_str_switch.push_back("off");
	}
	if (b_difpat_avg) {
		v_str_switch.push_back("on");
	}
	else {
		v_str_switch.push_back("off");
	}
	if (b_image) {
		v_str_switch.push_back("on");
	}
	else {
		v_str_switch.push_back("off");
	}
	if (b_image_avg) {
		v_str_switch.push_back("on");
	}
	else {
		v_str_switch.push_back("off");
	}
	std::cout << std::endl;
	std::cout << "  Pixelated detector settings         : status" << std::endl;
	std::cout << "   <1> scanned diffraction (4-D STEM) : " << v_str_switch[0] << std::endl;
	std::cout << "   <2> averaged diffraction (PACBED)  : " << v_str_switch[1] << std::endl;
	std::cout << "   <3> scanned imaging                : " << v_str_switch[2] << std::endl;
	std::cout << "   <4> averaged imaging (ISTEM)       : " << v_str_switch[3] << std::endl;
	v_str_switch.clear();
}


int prm_detector::setup_pixelated()
{
	int ierr = 0, ichk = 0, i = 0, idet = 0, j = 0;
	std::string  stmp, sprm;
	
	
	if (binteractive) {
	_repeat_input: // display the current setup of annular detectors
		// allow to modify the setup
		print_setup_pixelated();
		std::cout << "  Switch a detector on/off or select <0> to proceed? ";
		std::cin >> ichk;
		switch (ichk) {
		case 0: // jump out and use this list
			break;
		case 1: // switch scanned diffraction
			b_difpat |= true;
			goto _repeat_input;
			break;
		case 2: // switch PACBED
			b_difpat_avg |= true;
			goto _repeat_input;
			break;
		case 3: // switch scanned images
			b_image |= true;
			goto _repeat_input;
			break;
		case 4:
			b_image_avg |= true;
			goto _repeat_input;
			break;
		default:
			goto _repeat_input;
			break;
		}
		
		if (b_difpat) {
			v_str_ctrl.push_back("set_det_dif");
		}
		if (b_difpat_avg) {
			v_str_ctrl.push_back("set_det_padif");
		}
		if (b_image) {
			v_str_ctrl.push_back("set_det_img");
		}
		if (b_image_avg) {
			v_str_ctrl.push_back("set_det_paimg");
		}

	}
	else { // read the setup from the control lines
		if (0 <= ctrl_find_param("set_det_dif", &stmp, j)) {
			b_difpat = true;
		}
		if (0 <= ctrl_find_param("set_det_padif", &stmp, j)) {
			b_difpat_avg = true;
		}
		if (0 <= ctrl_find_param("set_det_img", &stmp, j)) {
			b_image = true;
		}
		if (0 <= ctrl_find_param("set_det_paimg", &stmp, j)) {
			b_image_avg = true;
		}

		if (btalk) {
			print_setup_pixelated(); // print the current temporary setup
		}
	}

	if (v_annular.size() > 0) {
		b_annular = true;
	}

	return ierr;
}


int prm_detector::get_jms_flags()
{
	int ndetflg = (int)_JMS_DETECT_NONE;
	if (b_annular) ndetflg += (int)_JMS_DETECT_INTEGRATED;
	if (b_difpat || b_difpat_avg) ndetflg += (int)_JMS_DETECT_DIFFRACTION;
	if (b_image || b_image_avg) ndetflg += (int)_JMS_DETECT_IMAGE;
	if (b_wave || b_waveft_avg) ndetflg += (int)_JMS_DETECT_WAVEREAL;
	if (b_waveft || b_waveft_avg) ndetflg += (int)_JMS_DETECT_WAVEFOURIER;
	return ndetflg;
}