// file: 'prm_sample.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementations of the class prm_sample.
//
// Class prm_sample handles parameters and setup routines around the
// TEM sample and numerical electron diffraction calculation.
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

#include "prm_sample.h"
#include "prm_slice.h"
#include "string_format.h"


prm_sample::prm_sample()
{
	input_form = 0;
	str_slc_file_pre = "slice";
	str_slc_file_suf = ".sli";
	slc_num = 0;
	slc_num_obj = 0;
	slc_num_det = 0;
	slc_det_per = 0;
	grid_nx = 0;
	grid_ny = 0;
	grid_a0 = 0.0f;
	grid_a1 = 0.0f;
	grid_ekv = 0.0f;
	slc_obj = NULL;
	slc_det = NULL;
	slc_dz = NULL;
	slc_var_num = NULL;
	slc_pgr_offset = NULL;
}


prm_sample::~prm_sample()
{
	if (NULL != slc_obj) free(slc_obj);
	if (NULL != slc_det) free(slc_det);
	if (NULL != slc_dz) free(slc_dz);
	if (NULL != slc_var_num) free(slc_var_num);
	if (NULL != slc_pgr_offset) free(slc_pgr_offset);
}


unsigned int prm_sample::setup_input_form()
{
	unsigned int iinf = input_form; // preset with current value
	int i;
	std::string stmp, sprm;

	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Supported sample data input forms: " << std::endl;
		std::cout << "  <0> No sample. *not supported yet*" << std::endl;
		std::cout << "  <1> Load structure file (CIF, CEL). *not supported yet*" << std::endl;
		std::cout << "  <2> Load transmission functions (SLI)." << std::endl;
		while (iinf < SAMPLE_INF_MIN || iinf > SAMPLE_INF_MAX) {
			std::cout << std::endl;
			std::cout << "  Select the input form of sample data: ";
			std::cin >> iinf;
		}
		v_str_ctrl.push_back("set_sample_input_form");
		v_str_ctrl.push_back(format("%d", iinf));
	}
	else {
		i = ctrl_find_param("set_sample_input_form", &stmp);
		if (i >= 0) {
			read_param(0, &stmp, &sprm);
			input_form = to_int(sprm);
		}
		iinf = input_form;
	}
	return iinf;
}

unsigned int prm_sample::setup_thickness()
{
	unsigned int ierr = 0, j = 0, k = 0;
	int i = -1, nobj = 0, nstep = 0;
	std::string stmp, sprm;
	float zmax = 0.f, zstep = 0.f;
	
	nobj = slc_num_obj;
	nstep = slc_det_per;

_repeat_input:
	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Set maximum sample thickness in number of slices: ";
		std::cin >> nobj;
		std::cout << std::endl;
		std::cout << "  * The program allows to calculate thickness series" << std::endl;
		std::cout << "    of images a periodic thickness steps. Setting this" << std::endl;
		std::cout << "    period to 0 slices means to calculate images at" << std::endl;
		std::cout << "    maximum sample thickness only." << std::endl;
		std::cout << std::endl;
		std::cout << "  Set detector readout period in number of slices: ";
		std::cin >> nstep;
	}
	else {
		i = ctrl_find_param("set_sample_thickness_max", &stmp);
		if (i >= 0) {
			read_param(0, &stmp, &sprm);
			slc_num_obj = to_int(sprm);
		}
		nobj = slc_num_obj;
		i = ctrl_find_param("set_sample_thickness_step", &stmp);
		if (i >= 0) {
			read_param(0, &stmp, &sprm);
			slc_det_per = to_int(sprm);
		}
		nstep = slc_det_per;
	}

	if (btalk || binteractive) {
		zmax = 0.f, zstep = 0.f;
		if (nobj > 0 && slc_num > 0) {
			for (i = 0; i < nobj; i++) {
				if (i == nstep) zstep = zmax;
				zmax += slc_dz[i%slc_num];
			}
		}
		std::cout << std::endl;
		std::cout << "  - maximum sample thickness: " << zmax << " nm" << std::endl;
		if (nstep > 0) {
			std::cout << "  - thickness series step size: " << zstep << " nm" << std::endl;
		}
		else {
			std::cout << "  - detection at maximum thickness." << std::endl;
		}
	}

	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Continue with these settings or change:  <0> Continue  <1> Change? ";
		std::cin >> i;
		if (i == 1) goto _repeat_input;
		slc_num_obj = nobj;
		slc_det_per = nstep;
		v_str_ctrl.push_back("set_sample_thickness_max");
		v_str_ctrl.push_back(format("%d",slc_num_obj));
		v_str_ctrl.push_back("set_sample_thickness_step");
		v_str_ctrl.push_back(format("%d", slc_det_per));
	}

	// this is where essential data is prepared
	// - allocate array of object slices
	if (NULL != slc_obj) {
		free(slc_obj);
		slc_obj = NULL;
	}
	slc_obj = (unsigned int*)malloc(sizeof(unsigned int)*slc_num_obj);
	// - allocate array of detection slices
	if (NULL != slc_det) {
		free(slc_det);
		slc_det = NULL;
	}
	slc_det = (int*)malloc(sizeof(int)*slc_num_obj);
	// - setup slice stack and detection planes
	slc_num_det = 0;
	for (i = 0; i < (int)slc_num_obj; i++) {
		j = (unsigned int)(i % slc_num);
		slc_obj[i] = j; // periodic structure stacking
		slc_det[i] = -1; // preset no detection for this slice
		if (slc_det_per > 0) { // periodic detection
			j = (unsigned int)(i % slc_det_per);
			if (0 == j) { // add plane to detection list
				slc_det[i] = (int)slc_num_det;
				slc_num_det++;
			}
			else {
				slc_det[i] = -1;
			}
		}
	}
	if (slc_det[slc_num_obj - 1] < 0) { // extra add exit plane
		slc_det[slc_num_obj - 1] = (int)slc_num_det;
		slc_num_det++;
	}

	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  Registered " << slc_num_det << " planes of signal detection." << std::endl;
	}

	return ierr;
}


unsigned int prm_sample::get_num_slc()
{
	return slc_num;
}

unsigned int prm_sample::get_num_slc_obj()
{
	return slc_num_obj;
}


unsigned int prm_sample::find_sli_files()
{
	unsigned int nslc = slc_num; // preset with current value
	unsigned int i = 0;
	int icmd;
	std::string sfile = "", stmp = "", sprm= "";
	std::string pre = str_slc_file_pre;
	str_slc_file_suf = ".sli";
	bool bfilefound = true;
repeat_input:
	if (binteractive) {
		nslc = 0; // reset number of slices to zero
		std::cout << std::endl;
		std::cout << "  Input a slice file name prefix: ";
		std::cin >> pre;
	}
	else { // get data from control file input
		icmd = ctrl_find_param("set_slc_num", &stmp);
		if (icmd >= 0) {
			read_param(0, &stmp, &sprm);
			slc_num = to_int(sprm); // we read this here though this input is not really used
		}
		nslc = slc_num;
		icmd = ctrl_find_param("set_slc_file_name", &stmp);
		if (icmd >= 0) {
			read_param(0, &stmp, &sprm);
			str_slc_file_pre = to_string(sprm);
		}
		pre = str_slc_file_pre;
	}
	i = 0;
	while (bfilefound) {
		i++;
		sfile = generate_ser_file_name(pre + "_", i, 3, str_slc_file_suf);
		bfilefound = file_exists(sfile);
	}
	if (i > 1) {
		nslc = i - 1;
	}
	if (nslc == 0) { // no file found
		sfile = generate_ser_file_name(pre + "_", 1, 3, str_slc_file_suf);
		if (binteractive) { // allow repeat of input
			std::cout << std::endl;
			std::cout << "  - Couldn't find file " << sfile << std::endl;
			goto repeat_input;
		}
		// exit
		std::cerr << "Input slice file not found: " << sfile << std::endl;
		goto exit;
	}
	else { // nslc files found
		if (btalk || binteractive) {
			std::cout << std::endl;
			std::cout << "  - found " << nslc << " slice files." << std::endl;
			std::cout << "  - file names: " << pre + "_###" + str_slc_file_suf << std::endl;
		}
		if (binteractive) {
			slc_num = nslc;
			v_str_ctrl.push_back("set_slc_num");
			v_str_ctrl.push_back(format("%d",slc_num));
			str_slc_file_pre = pre;
			v_str_ctrl.push_back("set_slc_file_name");
			v_str_ctrl.push_back("'" + str_slc_file_pre + "'");
		}
	}
	
exit:
	return nslc;
}

int prm_sample::load_sli_file_headers(void)
{
	int ierr = 0;
	if (slc_num > 0) { // get ready to receive slc_num data sets
		// free existing allocation
		if (NULL != slc_dz) free(slc_dz); slc_dz = NULL;
		if (NULL != slc_var_num) free(slc_var_num); slc_var_num = NULL;
		if (NULL != slc_pgr_offset) free(slc_pgr_offset); slc_pgr_offset = NULL;
		// allocate new memory
		slc_dz = (float*)malloc((size_t)slc_num * sizeof(float));
		if (NULL == slc_dz) {
			ierr = 1;
			std::cerr << "Failed to allocate memory (slc_dz)." << std::endl;
			goto exit;
		}
		slc_var_num = (unsigned int*)malloc((size_t)slc_num * sizeof(unsigned int));
		if (NULL == slc_var_num) {
			ierr = 2;
			std::cerr << "Failed to allocate memory (slc_var_num)." << std::endl;
			goto exit;
		}
		slc_pgr_offset = (unsigned long long*)malloc((size_t)slc_num * sizeof(unsigned long long));
		if (NULL == slc_pgr_offset) {
			ierr = 3;
			std::cerr << "Failed to allocate memory (slc_pgr_offset)." << std::endl;
			goto exit;
		}
		// fill content
		std::string sfile;
		prm_slice slc;
		for (unsigned int i = 0; i < slc_num; i++) {
			slc_dz[i] = 0.0f;
			slc_var_num[i] = 0;
			slc_pgr_offset[i] = 0;
			sfile = generate_ser_file_name(str_slc_file_pre + "_", i+1, 3, str_slc_file_suf);
			ierr = slc.load_ems_header(sfile);
			if (ierr > 0) {
				ierr = 100 + ierr;
				std::cerr << "Failed to read ems header from file " << sfile << std::endl;
				goto exit;
			}
			slc_dz[i] = slc.get_thickness();
			slc_var_num[i] = slc.get_variant_num();
			slc_pgr_offset[i] = slc.get_file_data_offset();
			if (i == 0) {
				slc.get_grid_dim(&grid_nx, &grid_ny);
				slc.get_grid_size(&grid_a0, &grid_a1);
				grid_ekv = slc.get_energy_kev();
				if (btalk) { // 
					std::cout << "    slice electron energy (keV): " << grid_ekv << std::endl;
					std::cout << "    slice (x, y) size (nm): (" << grid_a0 << ", " << grid_a1 << ")" << std::endl;
					std::cout << "    slice (x, y) grid size: (" << grid_nx << ", " << grid_ny << ")" << std::endl;
				}
			}
			if (btalk) {
				std::cout << "    slice #" << i << ": thickness (nm): " << slc_dz[i] << ", variants: " << slc_var_num[i] << std::endl;
			}
		}
	}
exit:
	return ierr;
}


unsigned int prm_sample::get_num_slc_det(void)
{
	return slc_num_det;
}

unsigned int* prm_sample::get_slc_obj(void)
{
	return slc_obj;
}

int* prm_sample::get_slc_det(void)
{
	return slc_det;
}

