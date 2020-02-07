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


prm_sample::prm_sample()
{
	input_form = 0;
	str_slc_file_pre = "slice";
	str_slc_file_suf = ".sli";
	grid_nx = 0;
	grid_ny = 0;
	grid_a0 = 0.0f;
	grid_a1 = 0.0f;
	grid_ekv = 0.0f;
	tilt_x = 0.0f;
	tilt_y = 0.0f;
}


prm_sample::~prm_sample()
{
}


void prm_sample::copy_data_from(prm_sample *psrc)
{
	if (psrc) {
		input_form = psrc->input_form;
		str_slc_file_pre = psrc->str_slc_file_pre;
		str_slc_file_suf = psrc->str_slc_file_suf;
		grid_a0 = psrc->grid_a0;
		grid_a1 = psrc->grid_a1;
		grid_ekv = psrc->grid_ekv;
		grid_nx = psrc->grid_nx;
		grid_ny = psrc->grid_ny;
		tilt_x = psrc->tilt_x;
		tilt_y = psrc->tilt_y;

		v_slc = psrc->v_slc;
		
		v_slc_obj = psrc->v_slc_obj;
		
		slc_det_per = 0;
		v_slc_det = psrc->v_slc_det;
	}
}


void prm_sample::copy_setup_from(prm_sample *psrc)
{
	if (psrc) {
		input_form = psrc->input_form;
		str_slc_file_pre = psrc->str_slc_file_pre;
		str_slc_file_suf = psrc->str_slc_file_suf;
		grid_a0 = psrc->grid_a0;
		grid_a1 = psrc->grid_a1;
		grid_ekv = psrc->grid_ekv;
		grid_nx = psrc->grid_nx;
		grid_ny = psrc->grid_ny;
		tilt_x = psrc->tilt_x;
		tilt_y = psrc->tilt_y;

		v_slc.clear();
		int nslc = (int)psrc->v_slc.size();
		if (nslc > 0) {
			prm_slice slc_tmp;
			for (int i = 0; i <= nslc; i++) {
				slc_tmp.header = psrc->v_slc[i].header;
				v_slc.push_back(slc_tmp);
			}
		}

		v_slc_obj = psrc->v_slc_obj;

		slc_det_per = 0;
		v_slc_det = psrc->v_slc_det;
	}
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
				v_slc_det[i+1] = (int)slc_num_det;
				slc_num_det++;
			}
			else {
				v_slc_det[i+1] = -1;
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


unsigned int prm_sample::setup_tilt()
{
	int ierr = 0, i = 0, j = 0;
	std::string  stmp, sprm;

	if (binteractive) {
	_repeat_input:
		std::cout << std::endl;
		std::cout << "  Set the sampel tilt to the optics axis along x in mrad: ";
		std::cin >> tilt_x;
		std::cout << std::endl;
		std::cout << "  Set the sampel tilt to the optics axis along y in mrad: ";
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

unsigned int prm_sample::get_sli_file_num(std::string sfile_name_prefix)
{
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
		for (unsigned int i = 0; i < nslc; i++) {
			sfile = generate_ser_file_name(str_slc_file_pre + "_", i+1, 3, str_slc_file_suf);
			ierr = slc.load_ems_header(sfile);
			if (ierr > 0) {
				ierr = 100 + ierr;
				std::cerr << "Failed to read ems header from file " << sfile << std::endl;
				goto exit;
			}
			v_slc_tmp.push_back(slc);
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
	}
exit:
	return ierr;
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


int prm_sample::get_variant_num_max(void)
{
	int n = (int)v_slc.size();
	int m = 0;
	if (n > 0) {
		for (int i = 0; i < n; i++) {
			m = __max(m, (int)v_slc[i].header.num_var);
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
			memset(pnslcvar, 0, sizeof(int)*nslc);
			for (islc = 0; islc < nslc; islc++) {
				pnslcvar[islc] = (int)v_slc[islc].header.num_var;
			}
		}
	}

	// setup phase grating memory in JMS
	ierr = pjms->PhaseGratingSetup(whichcode, grid_nx, grid_ny, nslc, nmaxvar, pnslcvar);
	if (ierr != 0) {
		nerr = 2;
		std::cerr << "Error: (setup_pgr_jms) phase grating setup faileed in multislice module (" << ierr << ")." << std::endl;
		goto _exit;
	}

	// transfer the phase grating data
	if (nslc > 0) {
		for (islc = 0; islc < nslc; islc++) {
			pgr = v_slc[islc].get_pgr();
			if (NULL == pgr) {
				nerr = 3;
				std::cerr << "Error: (setup_pgr_jms) no grating data for slice #" << islc << "." << std::endl;
				goto _exit;
			}
			ierr = pjms->SetPhaseGratingData(whichcode, islc, pnslcvar[islc], pgr);
			if (ierr != 0) {
				nerr = 4;
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