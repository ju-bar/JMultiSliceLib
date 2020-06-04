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
#include "prime_numbers.h"
#include <thread>



// ----------------------------------------------------------------------------
//
// class prm_sample
//
// ----------------------------------------------------------------------------


prm_sample::prm_sample()
{
	input_form = 0;
	structure_file_format = 0;
	str_slc_file_pre = "slice";
	str_slc_file_suf = ".sli";
	str_structure_file = "";
	grid_nx = 0;
	grid_ny = 0;
	grid_nz = 0;
	grid_a0 = 0.0f;
	grid_a1 = 0.0f;
	grid_a2 = 0.0f;
	grid_ekv = 0.0f;
	tilt_x = 0.0f;
	tilt_y = 0.0f;
	pot3d = 0;
	mod_atffacs = 0;
	mod_thermal = 0;
	mod_absorb = 0;
	vabf = 0.f;
	num_var = 1;
	save_sli = 0;
	lis_exc_max = 0; // 0 = no low-loss inelastic scattering
	lis_qe = 0.f;
	lis_qc = 0.f;
	lis_mfp = 0.f;
}


prm_sample::~prm_sample()
{
}


void prm_sample::copy_data_from(prm_sample *psrc)
{
	if (psrc) {
		input_form = psrc->input_form;
		structure_file_format = psrc->structure_file_format;
		str_structure_file = psrc->str_structure_file;
		str_slc_file_pre = psrc->str_slc_file_pre;
		str_slc_file_suf = psrc->str_slc_file_suf;
		grid_a0 = psrc->grid_a0;
		grid_a1 = psrc->grid_a1;
		grid_a2 = psrc->grid_a2;
		grid_ekv = psrc->grid_ekv;
		grid_nx = psrc->grid_nx;
		grid_ny = psrc->grid_ny;
		grid_nz = psrc->grid_nz;
		tilt_x = psrc->tilt_x;
		tilt_y = psrc->tilt_y;
		pot3d = psrc->pot3d;
		mod_atffacs = psrc->mod_atffacs;
		mod_thermal = psrc->mod_thermal;
		mod_absorb = psrc->mod_absorb;
		vabf = psrc->vabf;
		num_var = psrc->num_var;
		save_sli = psrc->save_sli;
		lis_exc_max = psrc->lis_exc_max;
		lis_qe = psrc->lis_qe;
		lis_qc = psrc->lis_qc;
		lis_mfp = psrc->lis_mfp;

		v_slc = psrc->v_slc;
		
		v_slc_obj = psrc->v_slc_obj;
		
		slc_det_per = 0;
		v_slc_det = psrc->v_slc_det;

		st_Input = psrc->st_Input;
		st_Used = psrc->st_Used;
	}
}


void prm_sample::copy_setup_from(prm_sample *psrc)
{
	if (psrc) {
		input_form = psrc->input_form;
		structure_file_format = psrc->structure_file_format;
		str_structure_file = psrc->str_structure_file;
		str_slc_file_pre = psrc->str_slc_file_pre;
		str_slc_file_suf = psrc->str_slc_file_suf;
		grid_a0 = psrc->grid_a0;
		grid_a1 = psrc->grid_a1;
		grid_a2 = psrc->grid_a2;
		grid_ekv = psrc->grid_ekv;
		grid_nx = psrc->grid_nx;
		grid_ny = psrc->grid_ny;
		grid_nz = psrc->grid_nz;
		tilt_x = psrc->tilt_x;
		tilt_y = psrc->tilt_y;
		pot3d = psrc->pot3d;
		mod_atffacs = psrc->mod_atffacs;
		mod_thermal = psrc->mod_thermal;
		mod_absorb = psrc->mod_absorb;
		vabf = psrc->vabf;
		num_var = psrc->num_var;
		save_sli = psrc->save_sli;
		lis_exc_max = psrc->lis_exc_max;
		lis_qe = psrc->lis_qe;
		lis_qc = psrc->lis_qc;
		lis_mfp = psrc->lis_mfp;

		v_slc.clear();
		int nslc = (int)psrc->v_slc.size();
		if (nslc > 0) {
			prm_slice slc_tmp;
			for (int i = 0; i < nslc; i++) {
				slc_tmp.header = psrc->v_slc[i].header;
				v_slc.push_back(slc_tmp);
			}
		}

		v_slc_obj = psrc->v_slc_obj;

		slc_det_per = 0;
		v_slc_det = psrc->v_slc_det;

		st_Input = psrc->st_Input;
		st_Used = psrc->st_Used;
	}
}


float prm_sample::get_wl()
{
	return (float)(_WLELKEV / sqrt(grid_ekv * (grid_ekv + 2. * _EEL0KEV)));
}


unsigned int prm_sample::setup_input_form()
{
	unsigned int iinf = input_form; // preset with current value
	int i;
	std::string stmp, sprm;

	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Supported sample data input forms: " << std::endl;
		std::cout << "  <1> Load structure file (CEL, CIF, XTL)." << std::endl;
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

int prm_sample::setup_thickness()
{
	unsigned int ierr = 0, j = 0, k = 0;
	int i = -1, l = -1, nobj = 0, nstep = 0, nslc = 0;
	int slc_num_det = 0;
	std::string stmp, sprm;
	float zmax = 0.f, zstep = 0.f;
	
	nslc = (int)v_slc.size();
	if (nslc == 0) {
		ierr = 1;
		std::cerr << "Error: Cannot setup sample thickness without slice information." << std::endl;
		goto _exit;
	}

	nobj = (int)v_slc_obj.size();
	nstep = slc_det_per;

	// TODO: Insert other options to setup readout at multiple thickness.

_repeat_input:
	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Set maximum sample thickness in number of slices: ";
		std::cin >> nobj;
		std::cout << std::endl;
		std::cout << "  * The program allows to calculate thickness series" << std::endl;
		std::cout << "    of images in periodic thickness steps. Setting this" << std::endl;
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
			nobj = to_int(sprm);
		}
		i = ctrl_find_param("set_sample_thickness_step", &stmp);
		if (i >= 0) {
			read_param(0, &stmp, &sprm);
			slc_det_per = to_int(sprm);
		}
		nstep = slc_det_per;
	}

	if (btalk || binteractive) {
		zmax = 0.f, zstep = 0.f;
		if (nobj > 0 && nslc > 0) {
			for (i = 0; i < nobj; i++) {
				if (i == nstep) zstep = zmax;
				zmax += v_slc[i%nslc].header.sz;
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
		slc_det_per = nstep;
		v_str_ctrl.push_back("set_sample_thickness_max");
		v_str_ctrl.push_back(format("%d", nobj));
		v_str_ctrl.push_back("set_sample_thickness_step");
		v_str_ctrl.push_back(format("%d", slc_det_per));
	}

	// this is where essential data is prepared
	v_slc_obj.clear();
	v_slc_obj.resize((size_t)nobj, 0);
	// - allocate array of detection slices
	v_slc_det.clear();
	v_slc_det.resize((size_t)nobj+1, -1);
	// - setup slice stack and detection planes
	if (slc_det_per > 0) { // setup with periodic readout
		v_slc_det[0] = 0; // -> add entrance plane
		slc_num_det = 1;
	}
	else { // setup wihout periodic readout
		v_slc_det[0] = -1; // -> don't detect on entrance plane
		slc_num_det = 0;
	}
	for (i = 0; i < nobj; i++) {
		j = (unsigned int)(i % nslc);
		v_slc_obj[i] = j; // periodic structure stacking
		if (slc_det_per > 0) { // periodic detection
			j = (unsigned int)((i+1) % slc_det_per);
			if (0 == j) { // add plane to detection list
				v_slc_det[(size_t)i+1] = (int)slc_num_det;
				slc_num_det++;
			}
			else {
				v_slc_det[(size_t)i+1] = -1;
			}
		}
	}
	if (v_slc_det[nobj] < 0) { // extra add exit plane
		v_slc_det[nobj] = slc_num_det;
		slc_num_det++;
	}

	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  Registered " << slc_num_det << " planes of signal detection." << std::endl;
		float fthick = 0.f;
		for (i = 0; i < (int)v_slc_det.size(); i++) {
			if (i > 0) {
				fthick += get_slc_obj_thickness((unsigned int)i);
			}
			l = v_slc_det[i];
			if (l >= 0) {
				std::cout << "  - detection #" << format("%03d", l) << " at plane #" << format("%03d", i) << " ( " << format("%6.2f", fthick) << " nm)" << std::endl;
			}
		}
	}

_exit:
	return ierr;
}


int prm_sample::setup_tilt()
{
	int ierr = 0, i = 0, j = 0;
	std::string  stmp, sprm;

	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Set the sample tilt to the optics axis along x in mrad: ";
		std::cin >> tilt_x;
		std::cout << std::endl;
		std::cout << "  Set the sample tilt to the optics axis along y in mrad: ";
		std::cin >> tilt_y;
		v_str_ctrl.push_back("set_sample_tilt");
		v_str_ctrl.push_back(format("%8.2f %8.2f", tilt_x, tilt_y));
	}
	else {
		i = ctrl_find_param("set_sample_tilt", &stmp);
		if (i >= 0) {
			j = read_param(0, &stmp, &sprm);
			tilt_x = to_float(sprm);
			j = read_param(j, &stmp, &sprm);
			tilt_y = to_float(sprm);
		}
		if (btalk) {
			std::cout << std::endl;
			std::cout << "  - sample tilt vector: " << format("(%7.2f, %7.2f)", tilt_x, tilt_y) << " mrad" << std::endl;
		}
	}

	return ierr;
}


int prm_sample::setup_lis()
{
	int ierr = 0, i = 0, j = 0, lis_mod = 0;
	std::string  stmp, sprm;
	float ep = 0.f; // plasmon energy (eV)
	float ek = 1000.f * grid_ekv; // kinetic energy of probe electrons (eV)
	float e0 = (float)(_EEL0EV); // rest energy of electrons (eV)
	float wthr = 0.01f, tol = 0.f , wt = 1.f, facn = 1.f;
	lis_exc_max = 0;
	lis_mfp = 0.f;
	lis_qe = 0.f;
	lis_qc = 0.f;
	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Supported input for low-loss inelastic scattering: " << std::endl;
		std::cout << "  <0> No low-loss inelastic scattering simulation" << std::endl;
		std::cout << "  <1> Plasmonic model (free electron gas, metals)" << std::endl;
		std::cout << "  <2> Generalized model" << std::endl;
		std::cout << std::endl;
		std::cout << "  Select form of low-loss inelastic scattering: ";
		std::cin >> lis_mod;
		while (lis_mod < SAMPLE_LIS_INF_MIN || lis_mod > SAMPLE_LIS_INF_MAX) {
			std::cout << std::endl;
			std::cout << "  Select form of low-loss inelastic scattering: ";
			std::cin >> lis_mod;
		}
		switch (lis_mod) {
		case 0:
			lis_exc_max = 0;
			break;
		case 1:
			std::cout << std::endl;
			std::cout << "  Set the plasmon energy in eV: ";
			std::cin >> ep;
			std::cout << std::endl;
			std::cout << "  Set the mean free path for inelastic scattering in nm: ";
			std::cin >> lis_mfp;
			// calculate the characteristic angle
			lis_qe = ep / ek * (ek + e0) / (ek + 2.f * e0);
			// estimate the critical angle
			lis_qc = sqrt(ep / e0);
			// estimate max. number of excitations
			tol = get_thickness_max() / lis_mfp;
			wt = 1.f - exp(-tol);
			facn = 1.;
			while (wt > wthr) { // include more excitation levels until remaining total probability of all higher levels is below the threshold wthr
				lis_exc_max++;
				facn = facn * (float)lis_exc_max; // n!
				wt -= (float)pow(tol, lis_exc_max) * exp(-tol) / facn;
			}
			break;
		case 2:
			std::cout << std::endl;
			std::cout << "  Set the maximum number of excitations per electron in the sample: ";
			std::cin >> lis_exc_max;
			std::cout << std::endl;
			std::cout << "  Set the mean free path for inelastic scattering in nm: ";
			std::cin >> lis_mfp;
			std::cout << std::endl;
			std::cout << "  Set the characteristic angle for inelastic scattering in mrad: ";
			std::cin >> lis_qe;
			std::cout << std::endl;
			std::cout << "  Set the critical angle for inelastic scattering in mrad: ";
			std::cin >> lis_qc;
			break;
		}

		v_str_ctrl.push_back("set_lis");
		v_str_ctrl.push_back(format("%d %10.4f %10.4f %10.4f", lis_exc_max, lis_mfp, lis_qe, lis_qc));
	}
	else {
		i = ctrl_find_param("set_lis", &stmp);
		if (i >= 0) {
			// TODO: implement input error handling by checking sprm.length() > 0
			j = read_param(0, &stmp, &sprm);
			lis_exc_max = (unsigned int)to_int(sprm);
			j = read_param(j, &stmp, &sprm);
			lis_mfp = to_float(sprm);
			j = read_param(j, &stmp, &sprm);
			lis_qe = to_float(sprm);
			j = read_param(j, &stmp, &sprm);
			lis_qc = to_float(sprm);
		}
		if (btalk && lis_exc_max > 0) {
			std::cout << std::endl;
			std::cout << "  Low-loss inelastic scattering parameters:" << std::endl;
			std::cout << "  - max. number of excitations per electron: " << lis_exc_max << std::endl;
			std::cout << "  - inelastic mean free path in nm: " << lis_mfp << std::endl;
			std::cout << "  - characteristic angle in mrad: " << lis_qe << std::endl;
			std::cout << "  - critical angle in mrad: " << lis_qc << std::endl;
		}
	}

	return ierr;
}

unsigned int prm_sample::get_sli_file_num(std::string sfile_name_prefix)
{
	// * TODO * //
	// This routine finds slice files with 3 digit indices only.
	// Expand it to find slice files independent of the number of digits used for the index.
	//
	bool bfilefound = true;
	unsigned int nfiles = 0;
	unsigned int i = 0;
	std::string sfile = "";
	i = 0;
	while (bfilefound) {
		i++;
		sfile = generate_ser_file_name(sfile_name_prefix + "_", i, 3, str_slc_file_suf);
		bfilefound = file_exists(sfile);
	}
	if (i > 1) {
		nfiles = i - 1;
	}
	return nfiles;
}

unsigned int prm_sample::find_sli_files()
{
	unsigned int nslc = (unsigned int)v_slc.size(); // preset with current value
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
			nslc = to_int(sprm); // we read this here though this input is not really used
		}
		icmd = ctrl_find_param("set_slc_file_name", &stmp);
		if (icmd >= 0) {
			read_param(0, &stmp, &sprm);
			str_slc_file_pre = to_string(sprm);
		}
		pre = str_slc_file_pre;
	}

	// get number of files with given prefix and current suffix str_slc_file_suf
	nslc = get_sli_file_num(pre);
	
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
		str_slc_file_pre = pre;
		if (btalk || binteractive) {
			std::cout << std::endl;
			std::cout << "  - found " << nslc << " slice files." << std::endl;
			std::cout << "  - file names: " << str_slc_file_pre + "_###" + str_slc_file_suf << std::endl;
		}
		if (binteractive) {
			v_str_ctrl.push_back("set_slc_num");
			v_str_ctrl.push_back(format("%d",nslc));
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
	unsigned int nslc;
	std::vector<prm_slice> v_slc_tmp;
	
	nslc = get_sli_file_num(str_slc_file_pre);
	
	if (nslc > 0) { // get ready to receive slc_num data sets
		// fill content
		std::string sfile;
		prm_slice slc;
		slc.set_ctrl(*this);
		grid_a2 = 0.f;
		for (unsigned int i = 0; i < nslc; i++) {
			sfile = generate_ser_file_name(str_slc_file_pre + "_", i+1, 3, str_slc_file_suf);
			ierr = slc.load_ems_header(sfile);
			if (ierr > 0) {
				ierr = 100 + ierr;
				std::cerr << "Failed to read ems header from file " << sfile << std::endl;
				goto exit;
			}
			v_slc_tmp.push_back(slc);
			grid_a2 += slc.header.sz;
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
				std::cout << "    slice #" << i << ": thickness (nm): " << slc.header.sz << ", variants: " << slc.header.num_var << std::endl;
			}
		}
		v_slc = v_slc_tmp;
		if (btalk) {
			std::cout << "    slice data total thickness (nm): " << grid_a2 << std::endl;
		}
	}
exit:
	return ierr;
}


int prm_sample::load_sli_file_data(void)
{
	int nerr = 0, ierr = 0;
	unsigned int nslc = (unsigned int)v_slc.size();

	if (nslc > 0) {
		std::string sfile;
		for (unsigned int i = 0; i < nslc; i++) {
			sfile = generate_ser_file_name(str_slc_file_pre + "_", i + 1, 3, str_slc_file_suf);
			ierr = v_slc[i].load_ems_data(sfile);
			if (0 < ierr) {
				nerr = 1;
				std::cerr << "Error: (load_sli_file_data) failed to load data from file " << sfile << " (" << ierr << ")." << std::endl;
				goto _exit;
			}
		}
	}

_exit:
	return nerr;
}

int prm_sample::load_structure(void)
{
	int nerr = 0, ierr = 0;
	int isff = 0;
	int icmd = 0;
	std::string stmp = "", sprm = "";

	if (binteractive) {
	repeat_input_format: // input a valid structure file format
		std::cout << std::endl;
		std::cout << "  Supported structure file formats: " << std::endl;
		std::cout << "  <1> CEL" << std::endl;
		std::cout << "  <2> CIF *not implemented yet*" << std::endl;
		std::cout << "  <3> XTL *not implemented yet*" << std::endl;
		std::cout << std::endl;
		std::cout << "  Select the structure file format: ";
		std::cin >> isff;
		if (isff < 1 || isff > STRUCTURE_FILE_FORMAT_MAX) goto repeat_input_format;
		structure_file_format = (unsigned int)isff;
	repeat_input_file: // input a structure file name
		std::cout << std::endl;
		std::cout << "  Input a structure file name: ";
		std::cin >> str_structure_file;
		if (!file_exists(str_structure_file)) {
			std::cerr << "Error (): file " + str_structure_file + " could not be found." << std::endl;
			goto repeat_input_file;
		}
		// store commands
		v_str_ctrl.push_back("set_structure_file_format");
		v_str_ctrl.push_back(format("%d", structure_file_format));
		v_str_ctrl.push_back("set_structure_file_name");
		v_str_ctrl.push_back(format("%d", str_structure_file));
	}
	else { // get data from control file input
		icmd = ctrl_find_param("set_structure_file_format", &stmp);
		if (icmd >= 0) {
			read_param(0, &stmp, &sprm);
			isff = to_int(sprm); // structure file format pre-load
			if (isff < 1 || isff > STRUCTURE_FILE_FORMAT_MAX) {
				std::cerr << "Error: (load_structure): invalid structure file format (" << isff << ")" << std::endl;
				nerr = 1;
				goto _exit;
			}
			structure_file_format = (unsigned int)isff;
		}
		icmd = ctrl_find_param("set_structure_file_name", &stmp);
		if (icmd >= 0) {
			read_param(0, &stmp, &sprm);
			str_structure_file = to_string(sprm);
		}
	}
	// Now load the structure
	switch (structure_file_format) {
	case 1: // CEL
		ierr = st_Input.LoadCEL(str_structure_file);
		if (ierr != 0) {
			std::cerr << "Error (load_structure): Failed to load structure data from CEL file: " << str_structure_file << std::endl;
			nerr = 600 + ierr;
			goto _exit;
		}
		break;
	case 2: // CIF
		//
		std::cerr << "Error (load_structure): Input of CIF file format is not supported yet." << std::endl;
		nerr = 700;
		goto _exit;
		//
		break;
	case 3: // XTL
		std::cerr << "Error (load_structure): Input of XTL file format is not supported yet." << std::endl;
		nerr = 800;
		goto _exit;
		//
		break;
	default:
		std::cerr << "Error (load_structure): Unsupported atomic structure file format (" << structure_file_format << ")." << std::endl;
		nerr = 900;
		goto _exit;
	}
	//
	if (btalk || binteractive) { // post some general information on the loaded structure data
		std::cout << std::endl;
		std::cout << "  Loaded atomic structure from file: " << str_structure_file << std::endl;
		std::cout << "  - cell dimensions (a, b, c) (nm): " << format("%9.5f %9.5f %9.5f", st_Input.m_vCellSize.x, st_Input.m_vCellSize.y, st_Input.m_vCellSize.z) << std::endl;
		std::cout << "  - cell angles (alpha, beta, gamma) (deg): " << format("%7.2f %7.2f %7.2f", st_Input.m_vCellAngles.x, st_Input.m_vCellAngles.y, st_Input.m_vCellAngles.z) << std::endl;
		if (!st_Input.IsOrthorhombic()) {
			std::cout << "  - Warning: The input atomic structure is not orthorhombic and should be modified before it can be used." << std::endl;
		}
		std::cout << "  - composition: " << st_Input.GetCompositionString() << std::endl;
		std::cout << "  - number of atom types: " << st_Input.m_vTypes.size() << std::endl;
		std::cout << "  - number of atoms: " << st_Input.m_vAtoms.size() << std::endl;
	}
	//
_exit:
	return nerr;
}


int prm_sample::setup_slices(void)
{
	int nerr = 0, ierr = 0;
	int i = 0, ichk = 0, nz = 1, nzsugg = 1;
	float dz = 0.f, fz = 0.f, f0 = 0.f, f1 = 0.f;
	std::string stmp = "", sprm = "";
	prm_slice slc_tmp;
	CStructure st_slc;

	slc_tmp.set_ctrl(*this);

	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  The atomic structure model has a size of " << format("%9.5f", st_Used.m_vCellSize.z) << " nm along the beam direction." << std::endl;
		nzsugg = st_Used.SuggestEquidistantZSlices(0.05f);
		if (nzsugg > 1) {
			std::cout << "  - suggested number of equidistant slices: " << nzsugg << " (" << format("%9.6f", st_Used.m_vCellSize.z / (float)nzsugg) << " nm/slice)" << std::endl;
		}
		std::cout << std::endl;
		std::cout << "  Enter the number of slices: ";
		std::cin >> ichk;
		nz = std::max(1, ichk);
		// store commands
		v_str_ctrl.push_back("set_structure_slices_equidistant");
		v_str_ctrl.push_back(format("%d", nz));
	}
	else { // get data from control file input
		i = ctrl_find_param("set_structure_slices_equidistant", &stmp);
		if (i >= 0) {
			read_param(0, &stmp, &sprm);
			nz = std::max(1, to_int(sprm)); // structure file format pre-load
		}
	}
	// nz is now used to generate slices and store the structure of each slice
	fz = 1.f / (float)nz; // fractional slice thickness
	dz = st_Used.m_vCellSize.z * fz; // absolute slice thickness (nm)
	v_slc.clear(); // erase previous slice data
	v_slc_det.clear();
	v_slc_obj.clear();
	if (btalk || binteractive) {
		std::cout << std::endl;
	}
	for (i = 0; i < nz; i++) {
		f0 = fz * i;
		f1 = fz * (i + 1);
		ierr = st_Used.GetZSliceStructure(&st_slc, f0, f1);
		if (0 < ierr) {
			std::cout << "Error (setup_slices): failed to create structure data of slice #" << i << " of " << nz << "." << std::endl;
			nerr = 1; goto _exit;
		}
		slc_tmp.st = st_slc;
		slc_tmp.header.sx = st_slc.m_vCellSize.x;
		slc_tmp.header.sy = st_slc.m_vCellSize.y;
		slc_tmp.header.sz = st_slc.m_vCellSize.z;
		stmp = st_slc.GetCompositionString();
		memset(slc_tmp.header.s_slc_name, 0, sizeof(char) * SLICENAME_MAX);
		snprintf(slc_tmp.header.s_slc_name, sizeof(char) * std::min((size_t)SLICENAME_MAX, stmp.size()), "%s", stmp.data());
		if (btalk || binteractive) {
			std::cout << "  - slice #" << i + 1 << ": " << stmp << std::endl;
		}
		v_slc.push_back(slc_tmp);
	}

_exit:
	return nerr;
}


int prm_sample::setup_grid(float k_max)
{
	int nerr = 0;
	int i = 0, j = 0;
	std::string stmp = "", sprm = "";
	bool bsuggest = false;
	float k_max_tmp = k_max;
	size_t npx = 0, npy = 0;
	if (k_max > 0.f) bsuggest = true;
	grid_nx = 512;
	grid_ny = 512;
	grid_nz = 1;

	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Setup the number of samples (pixels) used to calculate object tansmission functions." << std::endl;
		if (bsuggest) { // suggest number of samples based on cell size and max. scattering vector
			k_max_tmp = std::max(k_max, 60.f); // use at least 60 1/nm (This is about where a Debye-Waller factor with B = 0.5 A^2 drops to 0.01.)
			// calculate the next integer sampling which would include k_max_tmp, considering a band-width limit at 2/3 Nyquist
			npx = (size_t)std::ceil(3.01f * st_Used.m_vCellSize.x * k_max_tmp);
			npy = (size_t)std::ceil(3.01f * st_Used.m_vCellSize.y * k_max_tmp);
			// update the number to the next number with low prime factors (2, 3, 5)
			npx = next_low_prime(npx, 5);
			npy = next_low_prime(npy, 5);
			//
			std::cout << "  - suggested number of samples for scattering up to k_max = " << k_max_tmp << " 1/nm: " << npx << " x " << npy << " pixels." << std::endl;
		}
		std::cout << std::endl;
		std::cout << "  Set number of samples along x: ";
		std::cin >> grid_nx;
		std::cout << std::endl;
		std::cout << "  Set number of samples along y: ";
		std::cin >> grid_ny;

		// store commands
		v_str_ctrl.push_back("set_slc_grid_samples");
		v_str_ctrl.push_back(format("%d %d", grid_nx, grid_ny));
	}
	else {
		i = ctrl_find_param("set_slc_grid_samples", &stmp);
		if (i >= 0) {
			j = read_param(0, &stmp, &sprm);
			if (0 == sprm.size()) { nerr = 1; goto _exit; }
			grid_nx = (unsigned int)to_int(sprm);
			j = read_param(j, &stmp, &sprm);
			if (0 == sprm.size()) { nerr = 2; goto _exit; }
			grid_ny = (unsigned int)to_int(sprm);
		}
	}

	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  - grid used to sample along x and y in the multislice calculation: " << grid_nx << " x " << grid_ny << " pixels." << std::endl;
		std::cout << "  - real space sampling rates (x,y) (nm/pix):     " << format("(%12.5E, %12.5E)", grid_a0 / grid_nx, grid_a1 / grid_ny) << std::endl;
		std::cout << "  - Fourier sampling rates (x,y) (1/nm/pix):      " << format("(%12.5E, %12.5E)", 1.0f/grid_a0, 1.0f/grid_a1) << std::endl;
		std::cout << "  - Frequency band-width limitation (x,y) (1/nm): " << format("(%12.5E, %12.5E)", 0.3333f * grid_nx / grid_a0, 0.3333f * grid_ny / grid_a1) << std::endl;
	}

_exit:
	return nerr;
}


int prm_sample::setup_slc_calc(void)
{
	int nerr = 0, ichk = 0;
	int i = 0, j = 0;
	std::string stmp = "", sprm = "";
	//
	// preset defaults
	mod_thermal = 0;
	mod_absorb = 0;
	num_var = 1;
	vabf = 0.f;
	save_sli = 0;
	//
	if (binteractive) {
		//
		// * TODO * let the user choose which form factors to use
		//
	_input_mod_thermal: // model the effect of thermal motion
		std::cout << std::endl;
		std::cout << "  Choose a model to simulate the thermal motion of atoms in the sample:" << std::endl;
		std::cout << "  <0> none" << std::endl;
		std::cout << "  <1> damping of diffracted beams by Debye-Waller factors (fast but no TDS)" << std::endl;
		std::cout << "  <2> quantum excitation of phonons (Einstein model, slow but more accurate)" << std::endl;
		std::cout << "  Enter a number from the above list: ";
		std::cin >> ichk;
		if (ichk < 0 || ichk > 2) goto _input_mod_thermal; // out of range repeat input
		mod_thermal = (unsigned int)ichk;
		if (mod_thermal == 2) {
			std::cout << std::endl;
			std::cout << "  Enter the number of thermal states to calculate: ";
			std::cin >> num_var;
		}
		else {
			num_var = 1;
		}
	_input_mod_absorb: // model absorptive form fators
		std::cout << std::endl;
		std::cout << "  Choose how absorptive form factors are calculated:" << std::endl;
		std::cout << "  <0> none" << std::endl;
		std::cout << "  <1> fix fraction of elastic form factors (Hashimoto et. al.)" << std::endl;
		if (mod_thermal == 1) {
			std::cout << "  <2> loss in the elastic channel due to TDS (Hall & Hirsch)" << std::endl;
		}
		if (mod_thermal == 0 || mod_thermal == 2) {
			std::cout << "  <2> loss due to scattering beyond the calculation band-width limit (Hall & Hirsch)" << std::endl;
		}
		std::cout << "  Enter a number from the above list: ";
		std::cin >> ichk;
		if (ichk < 0 || ichk > 2) goto _input_mod_absorb; // out of range repeat input
		mod_absorb = (unsigned int)ichk;
		if (mod_absorb == 1) {
			std::cout << std::endl;
			std::cout << "  Enter the strength of absorptive form factors, relative to elastic form factors: ";
			std::cin >> vabf;
			vabf = std::abs(vabf);
		}
		//
		std::cout << std::endl;
		std::cout << "  Do you want to save transmission functions to files?  <0> No  <1> Yes : ";
		std::cin >> ichk;
		if (ichk == 1) {
			std::cout << std::endl;
			std::cout << "  Enter a file name (prefix): ";
			std::cin >> str_slc_file_pre;
			str_slc_file_suf = ".sli";
			if (str_slc_file_pre.size() == 0) str_slc_file_pre = "slice"; // give a default name 
			save_sli = 1;
		}
		//
		// store commands
		v_str_ctrl.push_back("set_slc_mod_thermal");
		v_str_ctrl.push_back(format("%d", mod_thermal));
		if (mod_thermal == 2) {
			v_str_ctrl.push_back("set_slc_var_num");
			v_str_ctrl.push_back(format("%d", num_var));
		}
		v_str_ctrl.push_back("set_slc_mod_absorb");
		v_str_ctrl.push_back(format("%d", mod_absorb));
		if (mod_absorb == 1) {
			v_str_ctrl.push_back("set_slc_abs_fac");
			v_str_ctrl.push_back(format("%9.6f", vabf));
		}
		v_str_ctrl.push_back("set_slc_save");
		v_str_ctrl.push_back(format("%d", save_sli));
		if (save_sli == 1) {
			v_str_ctrl.push_back("set_slc_file_name");
			v_str_ctrl.push_back(str_slc_file_pre);
		}
	}
	else {
		// mod_thermal
		i = ctrl_find_param("set_slc_mod_thermal", &stmp);
		if (i >= 0) {
			j = read_param(0, &stmp, &sprm);
			if (0 == sprm.size()) {
				std::cerr << "Error (setup_slc_calc): missing parameter on command 'set_slc_mod_thermal'." << std::endl;
				nerr = 1; goto _exit;
			}
			mod_thermal = (unsigned int)to_int(sprm);
			if (mod_thermal > 2) {
				std::cerr << "Error (setup_slc_calc): parameter on command 'set_slc_mod_thermal' out of range (0 - 2)." << std::endl;
				nerr = 1001; goto _exit;
			}
		}
		if (mod_thermal == 2) { // var_num
			i = ctrl_find_param("set_slc_var_num", &stmp);
			if (i >= 0) {
				j = read_param(0, &stmp, &sprm);
				if (0 == sprm.size()) {
					std::cerr << "Error (setup_slc_calc): missing parameter on command 'set_slc_var_num'." << std::endl;
					nerr = 2; goto _exit;
				}
				num_var = (unsigned int)to_int(sprm);
			}
		}
		// mod_absorb
		i = ctrl_find_param("set_slc_mod_absorb", &stmp);
		if (i >= 0) {
			j = read_param(0, &stmp, &sprm);
			if (0 == sprm.size()) {
				std::cerr << "Error (setup_slc_calc): missing parameter on command 'set_slc_mod_absorb'." << std::endl;
				nerr = 3; goto _exit;
			}
			mod_absorb = (unsigned int)to_int(sprm);
			if (mod_absorb > 2) {
				std::cerr << "Error (setup_slc_calc): parameter on command 'set_slc_mod_absorb' out of range (0 - 2)." << std::endl;
				nerr = 1003; goto _exit;
			}
		}
		if (mod_absorb == 1) { // vabf
			i = ctrl_find_param("set_slc_abs_fac", &stmp);
			if (i >= 0) {
				j = read_param(0, &stmp, &sprm);
				if (0 == sprm.size()) {
					std::cerr << "Error (setup_slc_calc): missing parameter on command 'set_slc_abs_fac'." << std::endl;
					nerr = 4; goto _exit;
				}
				vabf = to_float(sprm);
			}
		}
		// save_sli
		i = ctrl_find_param("set_slc_save", &stmp);
		if (i >= 0) {
			j = read_param(0, &stmp, &sprm);
			if (0 == sprm.size()) {
				std::cerr << "Error (setup_slc_calc): missing parameter on command 'set_slc_save'." << std::endl;
				nerr = 5; goto _exit;
			}
			save_sli = (unsigned int)to_int(sprm);
			if (save_sli > 1) {
				std::cerr << "Error (setup_slc_calc): parameter on command 'set_slc_save' out of range (0 - 1)." << std::endl;
				nerr = 1005; goto _exit;
			}
		}
		if (save_sli > 0) {
			i = ctrl_find_param("set_slc_file_name", &stmp);
			if (i >= 0) {
				j = read_param(0, &stmp, &sprm);
				if (0 == sprm.size()) {
					std::cerr << "Error (setup_slc_calc): missing parameter on command 'set_slc_file_name'." << std::endl;
					nerr = 6; goto _exit;
				}
				str_slc_file_pre = to_string(sprm);
				str_slc_file_suf = ".sli"; // only sli format supported yet
			}
		}
	}

	if (btalk || binteractive) {
		std::cout << std::endl;
		switch (mod_thermal) {
		case 1:
			std::cout << "  - applying Debye-Waller factors to the potentials." << std::endl;
			break;
		case 2:
			std::cout << "  - calculating " << num_var << " thermal states for calculations in the QEP model." << std::endl;
			break;
		}
		switch (mod_absorb) {
		case 1:
			std::cout << "  - calculating absorptive form factors with relative strength: " << format("%8.4f", vabf) << std::endl;
			break;
		case 2:
			if (mod_thermal == 1) {
				std::cout << "  - calculating absorptive form factors due to thermal diffuse scattering." << std::endl;
			}
			if (mod_thermal == 0 || mod_thermal == 2) {
				std::cout << "  - calculating absorptive form factors due to scattering beyond the band-width limit." << std::endl;
			}
			break;
		}
		if (save_sli > 0) {
			std::cout << "  - calculated object transmission functions will be saved to files:" << std::endl;
			std::cout << "    " << str_slc_file_pre << "_###" << str_slc_file_suf << std::endl;
		}
	}

_exit:
	return nerr;
}

unsigned int prm_sample::get_num_slc_det(void)
{
	unsigned int num_slc_det = 0;
	unsigned int num_slc_lst = (unsigned int)v_slc_det.size();
	if (num_slc_lst > 0) {
		for (unsigned int i = 0; i < num_slc_lst; i++) {
			if (v_slc_det[i] >= 0) {
				num_slc_det++;
			}
		}
	}
	return num_slc_det;
}

float prm_sample::get_slc_obj_thickness(unsigned int idx)
{
	float f_dz = 0.f;
	unsigned int num_slc_obj = (unsigned int)v_slc_obj.size();
	unsigned int num_slc = (unsigned int)v_slc.size();
	if (num_slc > 0 && num_slc_obj > 0 && idx < num_slc_obj) {
		unsigned int idx_slc = v_slc_obj[idx];
		if (idx_slc < num_slc) {
			f_dz = v_slc[idx_slc].header.sz;
		}
	}
	return f_dz;
}

float prm_sample::get_thickness_max()
{
	float thick = 0.f;
	size_t num_slc_obj = v_slc_obj.size();
	size_t num_slc = v_slc.size();
	size_t islc_obj = 0, islc = 0;
	if (num_slc > 0 && num_slc_obj > 0) {
		for (islc_obj = 0; islc_obj < num_slc_obj; islc_obj++) {
			islc = (size_t)v_slc_obj[islc_obj];
			thick += v_slc[islc].get_thickness();
		}
	}
	return thick;
}

float prm_sample::get_thickness_at(unsigned int idx)
{
	float thick = 0.f;
	size_t num_slc_obj = v_slc_obj.size();
	size_t num_slc = v_slc.size();
	size_t islc_obj = 0, islc = 0, nobjslc = std::min((size_t)idx, num_slc_obj);
	if (num_slc > 0 && nobjslc > 0) {
		for (islc_obj = 0; islc_obj < nobjslc; islc_obj++) {
			islc = (size_t)v_slc_obj[islc_obj];
			thick += v_slc[islc].get_thickness();
		}
	}
	return thick;
}


int prm_sample::get_variant_num_max(void)
{
	int n = (int)v_slc.size();
	int m = 0;
	if (n > 0) {
		for (int i = 0; i < n; i++) {
			m = std::max(m, (int)v_slc[i].header.num_var);
		}
	}
	return m;
}

int prm_sample::setup_pgr_jms(CJMultiSlice *pjms, int whichcode)
{
	int nerr = 0, ierr = 0;
	int nslc = 0, islc = 0;
	int nmaxvar = 0;
	int *pnslcvar = NULL;
	fcmplx *pgr = NULL;

	if (NULL == pjms) {
		nerr = 1;
		std::cerr << "Error: (setup_pgr_jms) invalid pointer to CJMultiSlice object." << std::endl;
		goto _exit;
	}

	nslc = (int)v_slc.size();

	// prepare variant index list and max variant count
	if (nslc > 0) {
		nmaxvar = get_variant_num_max();
		if (nmaxvar > 0) {
			pnslcvar = (int*)malloc(sizeof(int)*nslc);
			if (NULL == pnslcvar) {
				nerr = 2;
				std::cerr << "Error: (setup_pgr_jms) failed to allocate slie variant list." << std::endl;
				goto _exit;
			}
			memset(pnslcvar, 0, sizeof(int)*nslc);
			for (islc = 0; islc < nslc; islc++) {
				pnslcvar[islc] = (int)v_slc[islc].header.num_var;
			}
		}
	}

	// setup phase grating memory in JMS
	ierr = pjms->PhaseGratingSetup(whichcode, grid_nx, grid_ny, nslc, nmaxvar, pnslcvar);
	if (ierr != 0) {
		nerr = 3;
		std::cerr << "Error: (setup_pgr_jms) phase grating setup failed in multislice module (" << ierr << ")." << std::endl;
		goto _exit;
	}

	// transfer the phase grating data
	if (nslc > 0) {
		for (islc = 0; islc < nslc; islc++) {
			pgr = v_slc[islc].get_pgr();
			if (NULL == pgr) {
				nerr = 4;
				std::cerr << "Error: (setup_pgr_jms) no grating data for slice #" << islc << "." << std::endl;
				goto _exit;
			}
			ierr = pjms->SetPhaseGratingData(whichcode, islc, pnslcvar[islc], pgr);
			if (ierr != 0) {
				nerr = 5;
				std::cerr << "Error: (setup_pgr_jms) failed to transfer phase grating data of slice #" << islc << "to multislice module (" << ierr << ")." << std::endl;
				goto _exit;
			}
			pjms->SetSliceThickness(islc, v_slc[islc].get_thickness());
		}
	}
	
_exit:
	if (NULL != pnslcvar) {	free(pnslcvar); }
	return nerr;
}

int prm_sample::setup_object_slice_seq_jms(CJMultiSlice *pjms)
{
	int nerr = 0;
	int njmserr = 0;
	int nslcobj = (int)v_slc_obj.size();
	int i = 0;
	int *pslcobj = NULL;

	if (NULL == pjms) {
		nerr = 1;
		std::cerr << "Error: (setup_object_slice_seq_jms) called with invalid interface to JMS." << std::endl;
		goto _exit;
	}

	if (nslcobj > 0) {
		pslcobj = (int*)malloc(sizeof(int)*nslcobj);
		for (i = 0; i < nslcobj; i++) {
			pslcobj[i] = v_slc_obj[i];
		}
	}
	njmserr = pjms->ObjectSliceSetup(nslcobj, pslcobj);
	if (0 < njmserr) {
		std::cerr << "Error: (setup_object_slice_seq_jms): failed to setup object slice sequence in multislice module (" << njmserr << ")." << std::endl;
		nerr = 2;
		goto _exit;
	}

_exit:
	if (NULL != pslcobj) { free(pslcobj); }
	return nerr;
}

int prm_sample::prepare_pro(CJMultiSlice *pjms, int whichcode)
{
	int nerr = 0, ierr = 0;
	int nslc = 0, islc = 0;
	int *proidx = NULL;

	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  Preparing propagator functions ..." << std::endl;
	}

	if (NULL == pjms) {
		nerr = 1;
		std::cerr << "Error: (prepare_pro) called with invalid interface to JMS." << std::endl;
		goto _exit;
	}

	nslc = (int)v_slc.size();

	// calculate propagator in object member pro
	if (nslc > 0) {
		// allocate propagator index list
		proidx = (int*)malloc(sizeof(int)*nslc);
		if (NULL == proidx) {
			std::cerr << "Error: (prepare_pro): failed to allocate propagator index list." << std::endl;
			nerr = 2;
			goto _exit;
		}
		// calculate propagators
		for (islc = 0; islc < nslc; islc++) {
			ierr = v_slc[islc].calculate_propagator(pjms, tilt_x, tilt_y);
			if (0 < ierr) {
				nerr = 3;
				std::cerr << "Error: (prepare_pro) failed to calculate propagator #" << islc << " (" << ierr << ")." << std::endl;
				goto _exit;
			}
			proidx[islc] = islc; // set index
		}
	}

	// setup propagator memory in JMS // also works in case of no slices with proidx == NULL
	ierr = pjms->PropagatorSetup(whichcode, nslc, proidx);
	if (0 < ierr) {
		std::cerr << "Error: (prepare_pro): failed to setup propagators in multislice module (" << ierr << ")." << std::endl;
		nerr = 4;
		goto _exit;
	}

	if (nslc > 0) {
		// transfer propagator data to JMS
		for (islc = 0; islc < nslc; islc++) {
			ierr = pjms->SetPropagatorData(whichcode, islc, v_slc[islc].get_pro());
			if (0 < ierr) {
				std::cerr << "Error: (prepare_slices): failed to transfer propagator data #" << islc << " to multislice module (" << ierr << ")." << std::endl;
				nerr = 5;
				goto _exit;
			}
		}
	}

_exit:
	if (NULL != proidx) { free(proidx); }
	return nerr;
}


int prm_sample::calculate_pgr(int cpu_num, int gpu_id) {
	int nerr = 0, ierr = 0;
	int ncpu = cpu_num, igpu = gpu_id;
	unsigned int nx = 0, ny = 0, nz = 0;
	size_t nslc = 0, islc = 0, nvar = 0, ivar = 0, npix = 0, idx = 0;
	size_t mem_tot = 0, mem_pgr = 0, mem_chk = 0;
	size_t natt = 0, iatt = 0, nat = 0;
	bool bdwf = (1 == mod_thermal);
	bool bsthread = false;
	CJEElaSca elsa;
	fcmplx mip(0.f, 0.f); // mean inner potential of the structure
	fcmplx* ppgr = NULL; // pointer to the current buffer receiving data
	long long cl_0 = 0, cl_1 = 0;
	std::string sfile = "";

	// clock_0
	cl_0 = clock.getmsec();

	if (btalk) {
		std::cout << std::endl;
		std::cout << "  Calculation of object transmission functions ..." << std::endl;
	}

	if (igpu >= 0) { // check for valid GPU ID
		if (igpu >= GetGPUNum()) {
			std::cout << std::endl;
			std::cout << "  Warning: Invalid GPU ID (" << igpu << ") turning GPU calculation off." << std::endl;
			igpu = -1;
		}
	}

	if (ncpu == 0 && igpu < 0) { // need a single CPU fallback?
		std::cout << std::endl;
		std::cout << "  Warning: no CPU or GPU selected, fallback to run on 1 CPU." << std::endl;
		ncpu = 1;
	}

	// get basic sampling data
	nslc = v_slc.size();
	nx = grid_nx;
	ny = grid_ny;
	nz = (unsigned int)nslc; // slice projection mode // needs to change if pot3d > 0 (not implemented yet)
	nvar = (size_t)num_var;
	npix = (size_t)nx * ny; // items per phase grating
	mem_pgr = sizeof(fcmplx) * npix; // bytes per single phasegrating
	mem_tot = nvar * nslc * mem_pgr; // total number of bytes
	ierr = elsa.SetStructure(st_Used);
	if (0 < ierr) {
		nerr = 1;
		std::cerr << "Error (prm_sample::calculate_pgr): Failed to setup structure in elastic scattering module." << std::endl;
		goto _exit;
	}
	ierr = elsa.SetGrid(nx, ny, nz, pot3d);
	if (0 < ierr) {
		nerr = 2;
		std::cerr << "Error (prm_sample::calculate_pgr): Failed to setup grid in elastic scattering module." << std::endl;
		goto _exit;
	}
	ierr = elsa.SetPhysics(grid_ekv, mod_atffacs, mod_absorb, vabf, mod_thermal);
	if (0 < ierr) {
		nerr = 3;
		std::cerr << "Error (prm_sample::calculate_pgr): Failed to setup physics switches in elastic scattering module." << std::endl;
		goto _exit;
	}
	if (nslc > 0) {
		for (islc = 0; islc < nslc; islc++) { // update atom site lists
			ierr = v_slc[islc].st.SetAtomSites(0.05f); // treat atoms closer than 0.5 A in a slice as atoms sharing one site
			if (0 < ierr) {
				nerr = 4;
				std::cerr << "Error (prm_sample::calculate_pgr): Failed to determine atomic sites in slice # " << islc+1 << "." << std::endl;
				goto _exit;
			}
		}
	}

	// allocate memory for phase gratings
	if (nslc > 0) {
		if (btalk) {
			std::cout << "  - allocating memory for object transmission functions (" << (mem_tot >> 20) << " MB)." << std::endl;
		}
		for (islc = 0; islc < nslc; islc++) {
			v_slc[islc].header.grid_x = nx;
			v_slc[islc].header.grid_y = ny;
			v_slc[islc].header.ekv = grid_ekv;
			v_slc[islc].header.num_var = num_var;
			mem_chk = v_slc[islc].alloc_pgr();
			if (mem_chk != mem_pgr * nvar) {
				nerr = 10;
				std::cerr << "Error (prm_sample::calculate_pgr): Failed to allocate memory for phase gratings of slice #" << islc+1 << "." << std::endl;
				goto _exit;
			}
		}
	}

	// Calculate form factors for all atom types
	ierr = elsa.CalculateFormFactors(btalk || binteractive);
	if (0 < ierr) {
		nerr = 20;
		std::cerr << "Error (prm_sample::calculate_pgr): Failed to calculate form factors (" << ierr <<  ")." << std::endl;
		goto _exit;
	}
	// Prepare helper data for the calculation of phase gratings
	ierr = elsa.PrepareScatteringFunctions(igpu, ncpu, btalk || binteractive);
	if (0 < ierr) {
		nerr = 30;
		std::cerr << "Error (prm_sample::calculate_pgr): Failed to prepare scattering functions (" << ierr << ")." << std::endl;
		goto _exit;
	}

	if (btalk || binteractive) {
		ierr = elsa.CalculateMeanInnerPotential(&st_Used, mip);
		if (0 < ierr) {
			nerr = 30;
			std::cerr << "Error (prm_sample::calculate_pgr): Failed to calculate mean inner potential (" << ierr << ")." << std::endl;
			goto _exit;
		}
		std::cout << "  - mean inner potential: " << format("%8.3f V  (absorptive: %8.3f V).", mip.real(), mip.imag()) << std::endl;
	}

	// Calculate
	bsthread = ((igpu >= 0 && ncpu == 0) || (igpu < 0 && ncpu == 1));

	if (bsthread) { // calculate in this thread using either this CPU process or the GPU
		
		for (islc = 0; islc < nslc; islc++) { // loop over all slices
			for (ivar = 0; ivar < nvar; ivar++) { // loop over all variants
				if (btalk) {
					std::cout << "  - calculating slice " << islc + 1 << " of " << nslc << ", variant " << ivar + 1 << " of " << nvar << "     \r";
				}
				ppgr = v_slc[islc].get_pgr(ivar);
				if (ncpu == 1) {
					ierr = elsa.CalculatePhasegratingCPU(&v_slc[islc].st, ppgr, 0);
				}
				if (igpu >= 0) {
					ierr = elsa.CalculatePhasegratingGPU(&v_slc[islc].st, ppgr);
				}
				if (0 < ierr) {
					nerr = 50;
					std::cerr << "Error (prm_sample::calculate_pgr): Calculation of object transmission function failed (slice #" << islc+1 << ", variant #" << ivar+1 << ") (" << ierr << ")." << std::endl;
					goto _exit;
				}
			}
		}
		if (btalk) {
			std::cout << std::endl;
		}
	}
	else { // generate a task queue and let other threads do the work
		
		prm_sample_calc_queue calc_queue;
		int uoff = -1;
		size_t num_threads = cpu_num;
		size_t num_tasks = 0, num_finished = 0;
		float perc_solved = 0.f, perc_solved_prev = 0.f, perc_solved_delta = 0.5f;
		std::thread* pthread = NULL;
		std::vector<std::thread*> v_pthreads;

		ierr = calc_queue.init(this, &elsa);
		if (0 < ierr || 0 == calc_queue.total()) {
			nerr = 60;
			std::cerr << "Error (prm_sample::calculate_pgr): Failed to initialize multi-thread task queue (" << ierr << ")." << std::endl;
			goto _exit;
		}
		num_tasks = calc_queue.total();
		num_finished = calc_queue.solved() + calc_queue.failed();
		if (num_tasks > 0) {
			if (gpu_id >= 0) { // start a gpu thread
				pthread = new std::thread(WorkerCalculatePGR, gpu_id, uoff, &calc_queue);
				v_pthreads.push_back(pthread);
			}
			if (cpu_num > 0) { // start cpu threads
				for (int ithread = 0; ithread < cpu_num; ithread++) {
					pthread = new std::thread(WorkerCalculatePGR, uoff, ithread, &calc_queue);
					v_pthreads.push_back(pthread);
				}
			}
			if (btalk) {
				std::cout << "  - calculating object transmission functions: " << format("%5.1f", 0.0) << "%     \r";
			}
			num_finished = calc_queue.solved() + calc_queue.failed();
			while (num_finished < num_tasks) { // still calculations running, threads running asynchonously
				if (btalk) {
					perc_solved = 100.f * (float)num_finished / num_tasks;
					if (perc_solved - perc_solved_prev >= perc_solved_delta) {
						std::cout << "  - calculating object transmission functions: " << format("%5.1f", perc_solved) << "%     \r";
						perc_solved_prev = perc_solved;
					}
				}
				std::this_thread::sleep_for(std::chrono::milliseconds(5));
				num_finished = calc_queue.solved() + calc_queue.failed();
			}
			if (btalk) {
				std::cout << "  - calculating object transmission functions: " << format("%5.1f", 100.0) << "%       " << std::endl;
			}
			if (calc_queue.failed() > 0) {
				nerr = 61;
				std::cerr << "Error (prm_sample::calculate_pgr): " << calc_queue.failed() << " of " << calc_queue.total() << " calculation tasks failed." << std::endl;
				goto _exit;
			}
		}
		
	}

	// clock_1
	cl_1 = clock.getmsec();
	if (btalk) {
		std::cout << std::endl;
		std::cout << "  Object transmission functions calculated in " << std::setprecision(6) << 0.001f * (cl_1 - cl_0) << " s." << std::endl;
	}

	if (save_sli == 1 && nslc > 0) { // store the phase grating data in sli files.
		int ndig = std::max((int)3, (int)std::ceil(std::log10((double)nslc))); // number of digits are at least 3
		for (islc = 0; islc < nslc; islc++) {
			sfile = generate_ser_file_name(str_slc_file_pre + "_", (unsigned int)islc + 1, (unsigned int)ndig, str_slc_file_suf);
			ierr = v_slc[islc].save_ems(sfile);
			if (0 < ierr) {
				nerr = 90;
				std::cerr << "Error (prm_sample::calculate_pgr): Failed to save object transmission functions of slice #" << islc+1 <<" (" << ierr << ")." << std::endl;
				goto _exit;
			}
		}
	}

_exit:
	return nerr;
}




// ----------------------------------------------------------------------------
//
// class prm_sample_calc_task
//
// ----------------------------------------------------------------------------

prm_sample_calc_task::prm_sample_calc_task()
{
	state = 0;
	pelsa = NULL;
	pst = NULL;
	ppgr = NULL;
}

prm_sample_calc_task::~prm_sample_calc_task()
{

}




// ----------------------------------------------------------------------------
//
// class prm_sample_calc_queue
//
// ----------------------------------------------------------------------------

prm_sample_calc_queue::prm_sample_calc_queue()
{
	m_state = 0;
	m_q = 0;
	m_len_q = 0;
	m_num_tasks_open = 0;
	m_num_tasks_solved = 0;
	m_num_tasks_failed = 0;
}

prm_sample_calc_queue::~prm_sample_calc_queue()
{
	if (NULL != m_q) { delete[] m_q; }
}

int prm_sample_calc_queue::init(prm_sample* psample, CJEElaSca* pelsa)
{
	int nerr = 0;
	size_t num_slc = 0;
	size_t num_tasks = 0;
	size_t itask = 0;

	if (NULL == psample || NULL == pelsa) {
		nerr = 1;
		std::cerr << "Error (prm_sample_calc_queue::init): invalid parameter." << std::endl;
		goto _exit;
	}

	if (NULL != m_q) {
		delete[] m_q;
		m_q = NULL;
		m_state = 0;
		m_len_q = 0;
		m_num_tasks_open = 0;
		m_num_tasks_solved = 0;
		m_num_tasks_failed = 0;
	}

	num_slc = psample->v_slc.size();
	num_tasks = num_slc * psample->num_var;

	if (num_tasks == 0) {
		// finish without tasks
		goto _exit;
	}

	m_q = new prm_sample_calc_task[num_tasks];
	itask = 0;
	for (size_t islc = 0; islc < num_slc; islc++) {
		for (size_t ivar = 0; ivar < psample->num_var; ivar++) {
			m_q[itask].pelsa = pelsa;
			m_q[itask].pst = &psample->v_slc[islc].st;
			m_q[itask].ppgr = psample->v_slc[islc].get_pgr(ivar);
			m_q[itask].state = 0;
			itask++;
		}
	}
	m_num_tasks_solved = 0;
	m_num_tasks_failed = 0;
	m_num_tasks_open = num_tasks;
	m_len_q = num_tasks;
	m_state = 1;

_exit:
	return nerr;
}

bool prm_sample_calc_queue::empty(void)
{
	return (0 == open());
}

size_t prm_sample_calc_queue::open(void)
{
	/*std::lock_guard<std::mutex> lock(guard);
	size_t new_num_tasks_open = 0;
	if (m_len_q > 0 && m_state == 1) {
		for (size_t itask = 0; itask < m_len_q; itask++) {
			if (m_q[itask].state == 0) {
				new_num_tasks_open++;
			}
		}
	}
	m_num_tasks_open = new_num_tasks_open;*/
	return m_num_tasks_open;
}

size_t prm_sample_calc_queue::solved(void)
{
	/*std::lock_guard<std::mutex> lock(guard);
	size_t new_num_tasks_solved = 0;
	if (m_len_q > 0 && m_state == 1) {
		for (size_t itask = 0; itask < m_len_q; itask++) {
			if (m_q[itask].state == 3) {
				new_num_tasks_solved++;
			}
		}
	}
	m_num_tasks_solved = new_num_tasks_solved;*/
	return m_num_tasks_solved;
}

size_t prm_sample_calc_queue::failed(void)
{
	/*std::lock_guard<std::mutex> lock(guard);
	size_t new_num_tasks_failed = 0;
	if (m_len_q > 0 && m_state == 1) {
		for (size_t itask = 0; itask < m_len_q; itask++) {
			if (m_q[itask].state == 5) {
				new_num_tasks_failed++;
			}
		}
	}
	m_num_tasks_failed = new_num_tasks_solved;*/
	return m_num_tasks_failed;
}

size_t prm_sample_calc_queue::total(void)
{
	return m_len_q;
}

size_t prm_sample_calc_queue::request(void)
{
	size_t itask = 0;
	if (!empty()) {
		std::lock_guard<std::mutex> lock(guard);
		for (itask = 0; itask < m_len_q; itask++) {
			if (m_q[itask].state == 0) {
				m_q[itask].state = 1;
				m_num_tasks_open--;
				break;
			}
		}
		return itask;
	}
	return NULL;
}

// get specific task parameters
int prm_sample_calc_queue::get_task(size_t idx, prm_sample_calc_task& task_copy)
{
	if (idx < m_len_q) {
		if (m_q[idx].state == 1) {
			/*task_copy.pelsa = m_q[idx].pelsa;
			task_copy.ppgr = m_q[idx].ppgr;
			task_copy.pst = m_q[idx].pst;*/
			task_copy = m_q[idx];
			return 0;
		}
	}
	return 1;
}

// set task calculation state
int prm_sample_calc_queue::set_task_calculate(size_t idx)
{
	if (idx < m_len_q) {
		if (m_q[idx].state == 1) {
			m_q[idx].state = 2;
			return 0;
		}
	}
	return 1;
}

// set task solved state
int prm_sample_calc_queue::set_task_solved(size_t idx)
{
	if (idx < m_len_q) {
		if (m_q[idx].state == 2) {
			m_q[idx].state = 3;
			m_num_tasks_solved++;
			return 0;
		}
	}
	return 1;
}

// set task cancel state
int prm_sample_calc_queue::set_task_cancel(size_t idx)
{
	if (idx < m_len_q) {
		if (m_q[idx].state == 2) {
			m_q[idx].state = 4;
			return 0;
		}
	}
	return 1;
}

// set task error state
int prm_sample_calc_queue::set_task_error(size_t idx)
{
	if (idx < m_len_q) {
		if (m_q[idx].state == 2) {
			m_q[idx].state = 5;
			m_num_tasks_failed++;
			return 0;
		}
	}
	return 1;
}



// ----------------------------------------------------------------------------
//
// thread worker functions
//
// ----------------------------------------------------------------------------

void WorkerCalculatePGR(int gpu_id, int cpu_thread_id, prm_sample_calc_queue* pqueue)
{
	bool bgpu = (gpu_id >= 0);
	bool bcpu = (cpu_thread_id >= 0);
	prm_sample_calc_task my_task;
	int ierr = 0;
	size_t itask = 0;
	if (bcpu && bgpu) return; // exit due to invalid code choice
	if (!bcpu && !bgpu) return; // exit due to invalid code choice
	if (NULL == pqueue) return; // exit due to invalid interface
	while (!pqueue->empty()) { // keep on solving tasks until there are no more
		itask = pqueue->request();
		ierr = pqueue->get_task(itask, my_task);
		if (0 == ierr) {
			pqueue->set_task_calculate(itask);
			if (bgpu) {
				ierr = my_task.pelsa->CalculatePhasegratingGPU(my_task.pst, my_task.ppgr);
			}
			else {
				ierr = my_task.pelsa->CalculatePhasegratingCPU(my_task.pst, my_task.ppgr, cpu_thread_id);
			}
			if (0 < ierr) {
				pqueue->set_task_error(itask); // error
			}
			else {
				pqueue->set_task_solved(itask); // solved
			}
		}
	}
}