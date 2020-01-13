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
#include "string_format.h"
#include <thread>
#include "JMultiSliceLib.h"


prm_main::prm_main()
{
	gpu_id = -1;
	cpu_num = 0;
	cpu_num_max = 0;
	str_out_file = "output";
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
			std::stringstream ss; ss << ncpu;
			v_str_ctrl.push_back(ss.str());
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
		ncpu = __min((int)(cpu_num_max - 1), __max(0, ncpu));
	}
	else { // fallback to default (0 or 1)
		// limit to reasonable range
		ncpu = __min((int)cpu_num_max, __max(0, ncpu));
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


int prm_main::setup_sample(void)
{
	int nerr = 0, ierr = 0;

	sample.btalk = btalk;
	sample.binteractive = binteractive;
	sample.ndebug = ndebug;

	sample.v_str_ctrl = v_str_ctrl;

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
	case (SAMPLE_INF_SLI):
		// get number of slice files -> number of slices
		sample.find_sli_files();
		if (sample.get_num_slc() == 0) {
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

	probe.btalk = btalk;
	probe.binteractive = binteractive;
	probe.ndebug = ndebug;

	probe.v_str_ctrl = v_str_ctrl;

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

	detector.btalk = btalk;
	detector.binteractive = binteractive;
	detector.ndebug = ndebug;

	detector.v_str_ctrl = v_str_ctrl;

	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  Detector setup ..." << std::endl;
	}

	ierr = detector.setup_annular();
	if (0 < ierr) {
		nerr = 100 + ierr;
		goto exit;
	}

	// TODO: setup other detection schemes

exit:
	if (nerr == 0 && binteractive && detector.v_str_ctrl.size() > 0) { // append detector input control to main control
		v_str_ctrl = detector.v_str_ctrl;
	}
	if (!binteractive) { // erase control input copy
		detector.v_str_ctrl.clear();
	}
	return nerr;
}