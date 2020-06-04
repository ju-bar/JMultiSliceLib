// file: 'prm_main.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementations of the class prm_main.
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

#include "prm_main.h"
#include <thread>

prm_main::prm_main()
{
	gpu_id = -1;
	cpu_num = 0;
	cpu_num_max = 0;
	str_out_file = "";
	tem_mode = (int)(_PRM_TEM_MODE_STEM);
}


prm_main::~prm_main()
{
}


void prm_main::print_cuda_device_info(int id)
{
	int nerr = 0;
	if (id >= 0) {
		char * smsg = new char[1024];
		std::cout << std::endl;
		std::cout << "  - using CUDA device #" << id << std::endl;
		GetGPUName(id, smsg);
		std::cout << "    " << smsg << std::endl;
		int iCMajor = 0, iCMinor = 0, iMaxThread = 0;
		int nMultiProc = 0, nCores = 0, nThrPProc = 0;
		int64_t CUDAmemtotal = 0, CUDAmemfree = 0;
		nerr = GetGPUStats(id, iCMajor, iCMinor, iMaxThread, CUDAmemtotal, CUDAmemfree);
		nerr += GetGPUCores(id, nMultiProc, nCores, nThrPProc);
		if (0 == nerr) {
			std::cout << "    Compute class         : " << iCMajor << "." << iCMinor << std::endl;
			std::cout << "    Max. threads per block: " << iMaxThread << std::endl;
			std::cout << "    Number of Multiprocs  : " << nMultiProc << std::endl;
			std::cout << "    Max. threads per proc.: " << nThrPProc << std::endl;
			std::cout << "    Number of cores       : " << nCores << std::endl;
			std::cout << "    Memory total      [MB]: " << CUDAmemtotal / 1048576 << std::endl;
			std::cout << "    Memory available  [MB]: " << CUDAmemfree / 1048576 << std::endl;
		}
		else {
			std::cerr << "Error: Failed to retrieve device information." << std::endl;
		}
		delete[] smsg;
	}
	else {
		std::cout << std::endl;
		std::cout << "    No GPU device." << std::endl;
	}
	return;
}


int prm_main::setup_cuda_device()
{
	int idev = -1;
	int ndev = 0;
	int i = 0;
	int nerr = 0;
	std::string stmp, sprm;

	idev = gpu_id; // preset result with current gpu id

	ndev = GetGPUNum();
	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  Found " << GetGPUNum() << " CUDA device(s)" << std::endl;
		if (ndev > 0) {
			char gpuname[256];
			int iCMajor = 0, iCMinor = 0, iMaxThread = 0;
			int64_t CUDAmemtotal = 0, CUDAmemfree = 0;
			for (i = 0; i < ndev; i++) {
				nerr = GetGPUName(i, gpuname);
				nerr += GetGPUStats(i, iCMajor, iCMinor, iMaxThread, CUDAmemtotal, CUDAmemfree);
				if (0 == nerr) {
					std::cout << "  <" << i << "> " << gpuname << " (" << (CUDAmemtotal >> 20) << " MB)" << std::endl;
				}
				else {
					std::cerr << "  <" << i << "> ** failure ** " << std::endl;
				}
			}
			std::cout << "  <" << ndev << "> Do not use CUDA." << std::endl;
		}
	}
	if (ndev > 0) {
		if (binteractive) { // interactive user input
			std::cout << std::endl;
			std::cout << "  Select CUDA device (0 ... " << ndev << "): ";
			std::cin >> idev;
			// handle deselect of CUDA (out of range)
			if (idev < 0 || idev >= ndev) {
				idev = -1;
			}
			// insert control command
			v_str_ctrl.push_back("set_gpu_id");
			v_str_ctrl.push_back(format("%d", idev));
		}
		else { // get from input command strings
			i = ctrl_find_param("set_gpu_id", &stmp);
			if (i >= 0) {
				read_param(0, &stmp, &sprm);
				gpu_id = to_int(sprm);
				if (gpu_id < 0) gpu_id = -1;
				if (gpu_id >= ndev) {
					std::cerr << "Error: Invalid GPU ID (" << gpu_id << "), valid range: 0 ... " << ndev - 1 << std::endl;
					std::cout << "  - CUDA deactivated." << std::endl;
					gpu_id = -1;
				}
			}
			idev = gpu_id;
		}
		if (idev >= 0 && idev < ndev) {
			if (0 != SetCurrentGPU(idev)) {
				std::cerr << "Error: Failed to set current CUDA device number to " << idev << std::endl;
				idev = -1;
			}
		}
	}
	else {
		idev = -1; // invalidate CUDA, no device
	}
	
	return idev;
}


int prm_main::setup_cpu(void)
{
	int ncpu = cpu_num; // preset result with current member value
	int i = 0;
	std::string stmp, sprm;

	cpu_num_max = std::thread::hardware_concurrency(); // get max number of cpu threads
	
	if (cpu_num_max > 1) {
		if (btalk || binteractive) {
			std::cout << std::endl;
			std::cout << "  Multi-thread CPU calculation setup ... " << std::endl;
			std::cout << "  - detected number of logical processors (CPUs): " << cpu_num_max << std::endl;
			//std::cout << "  - current number of CPU calculation threads: " << cpu_num << std::endl;
		}
		if (binteractive) { // invalid ncpu preset, request new value
			std::cout << std::endl;
			std::cout << "  Select the number of CPU calculation threads (0 ... " << cpu_num_max - 1 << "): ";
			std::cin >> ncpu;
			// insert control command
			v_str_ctrl.push_back("set_cpu_num");
			v_str_ctrl.push_back(format("%i", ncpu));
		}
		else {
			i = ctrl_find_param("set_cpu_num", &stmp);
			if (i >= 0) {
				read_param(0, &stmp, &sprm);
				cpu_num = to_int(sprm);
			}
			ncpu = cpu_num;
		}
		// limit to reasonable range
		ncpu = std::min((int)(cpu_num_max - 1), std::max(0, ncpu));
	}
	else { // fallback to default (0 or 1)
		// limit to reasonable range
		ncpu = std::min((int)cpu_num_max, std::max(0, ncpu));
	}
	if (btalk) {
		if (ncpu <= 0) {
			std::cout << "  - no CPU calculation possible." << std::endl;
		}
		if (ncpu == 1) {
			std::cout << "  - CPU calculation restricted one thread." << std::endl;
		}
		if (ncpu > 1) {
			std::cout << "  - multi-threaded CPU calculation with " << ncpu << " threads." << std::endl;
		}
	}
	return ncpu;
}


int prm_main::setup_file_output(void)
{
	int ierr = 0, i = 0;
	std::string stmp = "", sprm = "";
	
	if (binteractive && str_out_file.length() == 0) { // interactive & no file name set yet
		std::cout << std::endl;
		std::cout << "  Enter the output file name prefix: ";
		std::cin >> str_out_file;
	}
	if (binteractive) { // store output file name command whether from current input or previous input (command line option -o)
		v_str_ctrl.push_back("set_output_file");
		v_str_ctrl.push_back("'" + str_out_file + "'");
	}
	else { // read from control file
		i = ctrl_find_param("set_output_file", &stmp);
		if (i >= 0) {
			read_param(0, &stmp, &sprm);
			str_out_file = to_string(sprm);
		}
	}

	return ierr;
}


int prm_main::setup_sample(void)
{
	int nerr = 0, ierr = 0;
	float k_max = 0.f;
	sample.set_ctrl(*this);
	
	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  Sample parameter setup ..." << std::endl;
	}

	sample.input_form = sample.setup_input_form();
	if (sample.input_form < SAMPLE_INF_MIN || sample.input_form > SAMPLE_INF_MAX) {
		ierr = 1;
		goto exit;
	}

	// switch input forms
	switch (sample.input_form)
	{
	//
	case (SAMPLE_INF_CEL): // input by atomic structure data files
		// load atomic structure data
		ierr = sample.load_structure();
		if (ierr != 0) {
			nerr = SAMPLE_INF_CEL * 10 + 1;
			goto exit;
		}
		// setup orientation & final cell size
		// * to be implemented later ... *
		sample.st_Used = sample.st_Input; // use the input cell now
		// write some data needed from other parameter objects to the sample object
		sample.grid_a0 = sample.st_Used.m_vCellSize.x;
		sample.grid_a1 = sample.st_Used.m_vCellSize.y;
		sample.grid_a2 = sample.st_Used.m_vCellSize.z;
		sample.grid_ekv = probe.ekv;
		// setup grid (x,y)
		k_max = 0.001f * detector.get_annular_betamax() / probe.get_wl(); // angle (mrad) -> to k (1/nm) in small angle approximation
		ierr = sample.setup_grid(k_max);
		if (ierr != 0) {
			nerr = SAMPLE_INF_CEL * 10 + 4;
			goto exit;
		}
		// setup grid (z) (slices)
		ierr = sample.setup_slices();
		if (ierr != 0) {
			nerr = SAMPLE_INF_CEL * 10 + 5;
			goto exit;
		}
		// setup parameters and options for calculation of potentials and transmission functions
		ierr = sample.setup_slc_calc();
		if (ierr != 0) {
			nerr = SAMPLE_INF_CEL * 10 + 6;
			goto exit;
		}
		//
		break;
	//
	case (SAMPLE_INF_SLI): // input by object transmission functions / slice files
		// get number of slice files -> number of slices
		if (0 == sample.find_sli_files()) {
			nerr = SAMPLE_INF_SLI * 10 + 1;
			goto exit;
		}
		// load slice header data from slices and store data of the first slice
		ierr = sample.load_sli_file_headers();
		if (ierr != 0) {
			nerr = SAMPLE_INF_SLI * 10 + 2;
			goto exit;
		}
		break;
	//
	default:
		nerr = 99;
		goto exit;
	}

	// set sample thickness / steps
	ierr = sample.setup_thickness();
	if (0 < ierr) {
		nerr = 100 + ierr;
		goto exit;
	}

	// set sample tilt
	ierr = sample.setup_tilt();
	if (0 < ierr) {
		nerr = 200 + ierr;
		goto exit;
	}

exit:
	if (nerr == 0 && binteractive && sample.v_str_ctrl.size()>0) { // append sample input control to main control
		v_str_ctrl = sample.v_str_ctrl;
	}
	if (!binteractive) { // erase control input copy
		sample.v_str_ctrl.clear();
	}
	return nerr;
}


int prm_main::setup_probe(void)
{
	int ierr = 0, nerr = 0;

	probe.set_ctrl(*this);

	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  Probe parameter setup ..." << std::endl;
	}

	ierr = probe.setup_ht();
	if (0 < ierr) {
		nerr = 100 + ierr;
		goto exit;
	}

	ierr = probe.setup_aperture();
	if (0 < ierr) {
		nerr = 200 + ierr;
		goto exit;
	}

	ierr = probe.setup_aberrations();
	if (0 < ierr) {
		nerr = 300 + ierr;
		goto exit;
	}

exit:
	if (nerr == 0 && binteractive && probe.v_str_ctrl.size() > 0) { // append probe input control to main control
		v_str_ctrl = probe.v_str_ctrl;
	}
	if (!binteractive) { // erase control input copy
		probe.v_str_ctrl.clear();
	}
	return nerr;
}



int prm_main::setup_detector(void)
{
	int ierr = 0, nerr = 0;

	detector.set_ctrl(*this);

	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  Detector setup ..." << std::endl;
	}

	// setup integrating annular detectors
	ierr = detector.setup_annular();
	if (0 < ierr) {
		nerr = 100 + ierr;
		goto exit;
	}

	// setup pixelated detectors
	ierr = detector.setup_pixelated();
	if (0 < ierr) {
		nerr = 200 + ierr;
		goto exit;
	}

	// * TODO * //
	// setup other detection schemes

exit:
	if (nerr == 0 && binteractive && detector.v_str_ctrl.size() > 0) { // append detector input control to main control
		v_str_ctrl = detector.v_str_ctrl;
	}
	if (!binteractive) { // erase control input copy
		detector.v_str_ctrl.clear();
	}
	return nerr;
}


int prm_main::setup_scan(void)
{
	int ierr = 0, nerr = 0;
	int ichk = 0;
	prm_scan scan_tmp;

	scan.set_ctrl(*this);
	

	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  Scanning setup ..." << std::endl;
	}

_repeat_input:
	
	scan_tmp = scan; // work on a temporary scan object
	
	ierr = scan_tmp.setup_beam_position();
	if (0 < ierr) {
		nerr = 100 + ierr;
		goto _exit;
	}
	ierr = scan_tmp.setup_frame_size();
	if (0 < ierr) {
		nerr = 200 + ierr;
		goto _exit;
	}
	ierr = scan_tmp.setup_frame_rotation();
	if (0 < ierr) {
		nerr = 300 + ierr;
		goto _exit;
	}
	ierr = scan_tmp.setup_sampling();
	if (0 < ierr) {
		nerr = 400 + ierr;
		goto _exit;
	}
	// set to default sub-frame
	scan_tmp.sub_ix1 = scan_tmp.nx - 1;
	scan_tmp.sub_iy1 = scan_tmp.ny - 1;
	ierr = scan_tmp.setup_subframe();
	if (0 < ierr) {
		nerr = 500 + ierr;
		goto _exit;
	}
	ierr = scan_tmp.setup_repeats();
	if (0 < ierr) {
		nerr = 600 + ierr;
		goto _exit;
	}

	if (btalk) {
		scan_tmp.print_setup();
	}

	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  <0> Accept  <1> Change  the current scan setup ? ";
		std::cin >> ichk;
		if (ichk == 1) goto _repeat_input;
	}

	scan = scan_tmp; // transfer data of accepted scan setup to parameters

_exit:
	if (nerr == 0 && binteractive && scan.v_str_ctrl.size() > 0) { // append scan input control to main control
		v_str_ctrl = scan.v_str_ctrl;
	}
	if (!binteractive) { // erase control input copy
		scan.v_str_ctrl.clear();
	}
	return nerr;
}


int prm_main::setup_qepopt(void)
{
	int nerr = 0;
	if (scan.num_repeat < 100) {
		std::cout << std::endl;
		std::cout << "Warning: You are about to setup additional options for calculations" << std::endl;
		std::cout << "         within the QEP model with a low number of statistical" << std::endl;
		std::cout << "         samples (" << scan.num_repeat << ")." << std::endl;
		std::cout << "         QEP = Quantum Excitations of Phonons" << std::endl;
		std::cout << "         Forbes et al., Phys.Rev.B 82, 104103 (2010)" << std::endl;
	}
	//sample.set_ctrl(*this);
	//// do something here
	//if (nerr == 0 && binteractive && sample.v_str_ctrl.size() > 0) { // append scan input control to main control
	//	v_str_ctrl = sample.v_str_ctrl;
	//}
	//if (!binteractive) { // erase control input copy
	//	sample.v_str_ctrl.clear();
	//}
	return nerr;
}


int prm_main::setup_lis(void)
{
	int nerr = 0;
	if (scan.num_repeat < 100) {
		std::cout << std::endl;
		std::cout << "Warning: You are about to setup the low-loss inelastic scattering" << std::endl;
		std::cout << "         Monte-Carlo procedure with a low number of statistical" << std::endl;
		std::cout << "         samples per scan position (" << scan.num_repeat << ")." << std::endl;
		std::cout << "         Barthel et al., Phys. Rev. B 101, 184109 (2020)." << std::endl;
	}
	sample.set_ctrl(*this);
	nerr = sample.setup_lis();
	if (nerr == 0 && binteractive && sample.v_str_ctrl.size() > 0) { // append scan input control to main control
		v_str_ctrl = sample.v_str_ctrl;
	}
	if (!binteractive) { // erase control input copy
		sample.v_str_ctrl.clear();
	}
	return nerr;
}


int prm_main::prepare_sample_pgr(void)
{
	int nerr = 0, ierr = 0;

	// Create or load phase gratings from information in member variable sample
	switch (sample.input_form) {
	case SAMPLE_INF_CEL:
		if (btalk || binteractive) {
			std::cout << std::endl;
			std::cout << "  Calculating object transmission functions ..." << std::endl;
		}
		ierr = sample.calculate_pgr(cpu_num, gpu_id);
		if (0 < ierr) {
			nerr = 2;
			std::cerr << "Error: (prepare_sample_pgr) failed to calculate object transmission functions (" << ierr << ")." << std::endl;
			goto _exit;
		}
		break;
	case SAMPLE_INF_SLI:
		if (btalk || binteractive) {
			std::cout << std::endl;
			std::cout << "  Loading object transmission functions ..." << std::endl;
		}
		ierr = sample.load_sli_file_data();
		if (0 < ierr) {
			nerr = 10;
			std::cerr << "Error: (prepare_sample_pgr) failed to load object transmission functions from files (" << ierr << ")." << std::endl;
		}
		break;
	default:
		nerr = 1;
		std::cerr << "Error: (prepare_sample_pgr) invalid input option for object transmission functions (" << sample.input_form << ")." << std::endl;
		goto _exit;
	}

_exit:
	return nerr;
}



int prm_main::prepare_result_params()
{
	int nerr = 0, ierr = 0;
	unsigned int num_ann = (unsigned int)detector.v_annular.size();
	unsigned int num_det_pln = sample.get_num_slc_det();

	if (detector.b_annular && num_ann > 0 && num_det_pln > 0) {
		stem_images.set_ctrl(*this);
		stem_images.data_type = (unsigned int)_RESULT_DATA_TYPE_FLOAT;
		stem_images.det_type = (unsigned int)_JMS_DETECT_INTEGRATED;
		stem_images.v_dim.push_back(scan.nx); // length of scan rows
		stem_images.v_dim.push_back(scan.ny); // number of scan rows
		stem_images.v_dim.push_back(num_det_pln); // number of thickness samples
		stem_images.v_dim.push_back(num_ann); // number of detectors
		stem_images.probe.copy_data_from(&probe); // copy probe parameters
		stem_images.sample.copy_setup_from(&sample); // copy sample setup (not all data)
		stem_images.detector.copy_setup_from(&detector); // copy detector parameters
		stem_images.scan = scan; // copy scan setup
		stem_images.str_out_file = str_out_file;
		stem_images.f_calc_scale = 1.f;
		stem_images.f_calc_weight = 0.f;
		if (0 < stem_images.init_buffer()) {
			nerr = 1;
			std::cerr << "Error: (prepare_result_params) failed to allocate scan image buffer." << std::endl;
			goto _exit;
		}
	}

	if (detector.b_difpat) {
		stem_pix_dif.set_ctrl(*this);
		stem_pix_dif.data_type = (unsigned int)_RESULT_DATA_TYPE_FLOAT;
		stem_pix_dif.det_type = (unsigned int)_JMS_DETECT_DIFFRACTION;
		stem_pix_dif.v_dim.push_back(sample.grid_nx); // grid row length
		stem_pix_dif.v_dim.push_back(sample.grid_ny); // grid column length
		stem_pix_dif.v_dim.push_back(num_det_pln); // number of thickness samples
		stem_pix_dif.probe.copy_data_from(&probe); // copy probe parameters
		stem_pix_dif.sample.copy_setup_from(&sample); // copy sample setup (not all data)
		stem_pix_dif.detector.copy_setup_from(&detector); // copy detector parameters
		stem_pix_dif.scan = scan; // copy scan setup
		stem_pix_dif.str_out_file = str_out_file;
		stem_pix_dif.f_calc_scale = 1.f;
		stem_pix_dif.f_calc_weight = 0.f;
		if (0 < stem_pix_dif.init_buffer()) {
			nerr = 2;
			std::cerr << "Error: (prepare_result_params) failed to allocate scanned diffraction buffer." << std::endl;
			goto _exit;
		}
	}

	if (detector.b_difpat_avg) {
		stem_pix_padif.set_ctrl(*this);
		stem_pix_padif.data_type = (unsigned int)_RESULT_DATA_TYPE_FLOAT;
		stem_pix_padif.det_type = (unsigned int)_JMS_DETECT_DIFFRACTION;
		stem_pix_padif.v_dim.push_back(sample.grid_nx); // grid row length
		stem_pix_padif.v_dim.push_back(sample.grid_ny); // grid column length
		stem_pix_padif.v_dim.push_back(num_det_pln); // number of thickness samples
		stem_pix_padif.probe.copy_data_from(&probe); // copy probe parameters
		stem_pix_padif.sample.copy_setup_from(&sample); // copy sample setup (not all data)
		stem_pix_padif.detector.copy_setup_from(&detector); // copy detector parameters
		stem_pix_padif.scan = scan; // copy scan setup
		stem_pix_padif.str_out_file = str_out_file;
		stem_pix_padif.f_calc_scale = 1.f;
		stem_pix_padif.f_calc_weight = 0.f;
		if (0 < stem_pix_padif.init_buffer()) {
			nerr = 2;
			std::cerr << "Error: (prepare_result_params) failed to allocate average diffraction buffer." << std::endl;
			goto _exit;
		}
	}

	if (detector.b_image) {
		stem_pix_img.set_ctrl(*this);
		stem_pix_img.data_type = (unsigned int)_RESULT_DATA_TYPE_FLOAT;
		stem_pix_img.det_type = (unsigned int)_JMS_DETECT_IMAGE;
		stem_pix_img.v_dim.push_back(sample.grid_nx); // grid row length
		stem_pix_img.v_dim.push_back(sample.grid_ny); // grid column length
		stem_pix_img.v_dim.push_back(num_det_pln); // number of thickness samples
		stem_pix_img.probe.copy_data_from(&probe); // copy probe parameters
		stem_pix_img.sample.copy_setup_from(&sample); // copy sample setup (not all data)
		stem_pix_img.detector.copy_setup_from(&detector); // copy detector parameters
		stem_pix_img.scan = scan; // copy scan setup
		stem_pix_img.str_out_file = str_out_file;
		stem_pix_img.f_calc_scale = 1.f / ((float)sample.grid_nx*sample.grid_ny);
		stem_pix_img.f_calc_weight = 0.f;
		if (0 < stem_pix_img.init_buffer()) {
			nerr = 2;
			std::cerr << "Error: (prepare_result_params) failed to allocate scanned image buffer." << std::endl;
			goto _exit;
		}
	}

	if (detector.b_image_avg) {
		stem_pix_paimg.set_ctrl(*this);
		stem_pix_paimg.data_type = (unsigned int)_RESULT_DATA_TYPE_FLOAT;
		stem_pix_paimg.det_type = (unsigned int)_JMS_DETECT_IMAGE;
		stem_pix_paimg.v_dim.push_back(sample.grid_nx); // grid row length
		stem_pix_paimg.v_dim.push_back(sample.grid_ny); // grid column length
		stem_pix_paimg.v_dim.push_back(num_det_pln); // number of thickness samples
		stem_pix_paimg.probe.copy_data_from(&probe); // copy probe parameters
		stem_pix_paimg.sample.copy_setup_from(&sample); // copy sample setup (not all data)
		stem_pix_paimg.detector.copy_setup_from(&detector); // copy detector parameters
		stem_pix_paimg.scan = scan; // copy scan setup
		stem_pix_paimg.str_out_file = str_out_file;
		stem_pix_paimg.f_calc_scale = 1.f / ((float)sample.grid_nx*sample.grid_ny);
		stem_pix_paimg.f_calc_weight = 0.f;
		if (0 < stem_pix_paimg.init_buffer()) {
			nerr = 2;
			std::cerr << "Error: (prepare_result_params) failed to allocate average image buffer." << std::endl;
			goto _exit;
		}
	}

_exit:
	return nerr;
}
