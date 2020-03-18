// file: 'JMSBench1.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the 'main' function.
// Program execution begins and ends there.
//
// This program runs multislice calculations.
// Main purpose is to benchmark the JMS library.
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
#include "multislice.h"
#include "JMultiSliceLib.h"


int parseoptions(prm_main *pprm, int argc, char* argv[])
{
	if (NULL == pprm) {
		std::cerr << "Error (parseoptions): missing parameter interface." << std::endl;
		return 1;
	}
	pprm->btalk = true;
	pprm->ndebug = 0;
	if (argc > 1) {

		std::string cmd;
		// parse other arguments
		for (int iarg = 1; iarg < argc; iarg++) {
			cmd = argv[iarg];
			if (cmd == "/silent") {
				pprm->btalk = false;
				continue;
			}
			if (cmd == "/debug") { // direct debug switch
				pprm->ndebug = 1; // set debug level 1
				continue;
			}

			if (cmd == "-dbgl") { // leveled debug switch + level
				iarg++;
				if (iarg >= argc) {
					std::cerr << "Error: expecting a debug level number after option -dbgl.\n";
					return 1;
				}
				pprm->ndebug = atoi(argv[iarg]);
				if (pprm->ndebug < 0) pprm->ndebug = 0;
				if (pprm->ndebug > 5) pprm->ndebug = 5;
				continue;
			}

			if (cmd == "-o" || cmd == "--output") { // modified output name + name
				iarg++;
				if (iarg >= argc) {
					std::cerr << "Error: expecting a file name string after option --output (-o).\n";
					return 1;
				}
				pprm->str_out_file = argv[iarg];
				continue;
			}
			if (cmd == "-c" || cmd == "--control") { // modified name of control file
				iarg++;
				if (iarg >= argc) {
					std::cerr << "Error: expecting a file name string after option --control (-c).\n";
					return 1;
				}
				pprm->str_ctrl_file = argv[iarg];
				continue;
			}
		}
	}
	if (pprm->ndebug > 0) { // debug overrides
		pprm->btalk = true;
	}
	return 0;
}

int init(prm_main *pprm)
{
	int ierr = 0;
	int nerr = 0;
	int ichk = 0;
	int ichk2 = 0;

	if (NULL == pprm) {
		nerr = 1;
		std::cerr << "Error (init): missing parameter interface." << std::endl;
		goto exit;
	}

	pprm->bstore = false;
	pprm->binteractive = false; // activate interactive mode to record input
	pprm->v_str_ctrl.clear(); // clear input sequence

	pprm->clock.start();

	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << "  +------------------------------------------+  " << std::endl;
		std::cout << "  |                                          |  " << std::endl;
		std::cout << "  |   JMSBench1 - version 1.0 - 2020-02-19   |  " << std::endl;
		std::cout << "  |                                          |  " << std::endl;
		std::cout << "  |   by J. Barthel                          |  " << std::endl;
		std::cout << "  |      ju.barthel@fz-juelich.de            |  " << std::endl;
		std::cout << "  |      Forschungszentrum Juelich GmbH      |  " << std::endl;
		std::cout << "  |      52425 Juelich, Germany              |  " << std::endl;
		std::cout << "  |                                          |  " << std::endl;
		std::cout << "  +------------------------------------------+  " << std::endl;
		std::cout << std::endl;
	}

	ierr = pprm->ctrl_init();
	if (0 < ierr) {
		std::cerr << "Failed to initialize program control." << std::endl;
		nerr = ierr;
		goto exit;
	}
	if (pprm->binteractive) {
		std::cout << std::endl;
		std::cout << "  Do you want to store your input to a file?  <0> No  <1> Yes: ";
		std::cin >> ichk;
		if (ichk == 1) {
_repeat_input:
			std::cout << std::endl;
			std::cout << "  Specify a control file name: ";
			std::cin >> pprm->str_ctrl_file;
			if (pprm->file_exists(pprm->str_ctrl_file)) {
				std::cout << std::endl;
				std::cout << "  File (" << pprm->str_ctrl_file << ") already exists." << std::endl;
				std::cout << "  Do you want to override this file:  <0> No  <1> Yes ? ";
				std::cin >> ichk2;
				if (ichk2 == 1) {
					std::ofstream ofs;
					ofs.open(pprm->str_ctrl_file, std::ofstream::out | std::ofstream::trunc);
					ofs.close();
				}
				else {
					goto _repeat_input;
				}
			}
			pprm->bstore = true;

		}
	}

exit:
	return nerr;
}



int setup(prm_main *pprm)
{
	int ierr = 0;
	int nerr = 0;

	if (NULL == pprm) {
		nerr = 1;
		std::cerr << "Error (setup): missing parameter interface." << std::endl;
		goto exit;
	}

	if (pprm->btalk || pprm->binteractive) {
		std::cout << std::endl;
		std::cout << "  Running program setup routines ..." << std::endl;
	}

	// basic and general parameters
	// - CUDA device
	if (pprm->gpu_id < 0) {
		pprm->gpu_id = pprm->setup_cuda_device();
		if (pprm->btalk) {
			pprm->print_cuda_device_info(pprm->gpu_id);
		}
	}
	// - CPU threading
	if (pprm->cpu_num < 1) {
		pprm->cpu_num = pprm->setup_cpu();
	}
	if (pprm->gpu_id < 0 && pprm->cpu_num <= 0) {
		nerr = 1001;
		std::cerr << "Failed to setup processors. Unable to calculate." << std::endl;
		goto exit;
	}
	if (pprm->btalk) {
		std::cout << std::endl;
		if (pprm->gpu_id >= 0 && pprm->cpu_num > 1) {
			std::cout << "  Parallel calculations with GPU and CPU." << std::endl;
		}
		if (pprm->gpu_id < 0 && pprm->cpu_num > 1) {
			std::cout << "  Parallel calculations with CPU." << std::endl;
		}
		if (pprm->gpu_id >= 0 && pprm->cpu_num <= 1) {
			std::cout << "  Calculations with GPU." << std::endl;
		}
		if (pprm->gpu_id < 0 && pprm->cpu_num == 1) {
			std::cout << "  Single thread CPU calculations." << std::endl;
		}
	}

	// - Output file name
	ierr = pprm->setup_file_output();
	if (0 < ierr) {
		nerr = 1500 + ierr;
		std::cerr << "Failed to setup output file name. Falling back to default 'output'." << std::endl;
		pprm->str_out_file = "output";
	}

	// - Sample
	ierr = pprm->setup_sample();
	if (0 < ierr) {
		nerr = 2000 + ierr;
		std::cerr << "Failed to setup sample data. Unable to calculate." << std::endl;
		goto exit;
	}

	// - Simulation (TEM or STEM) *currently STEM only*

	// - Probe
	ierr = pprm->setup_probe();
	if (0 < ierr) {
		nerr = 3000 + ierr;
		std::cerr << "Failed to setup probe data. Unable to calculate." << std::endl;
		goto exit;
	}

	// - Detectors
	ierr = pprm->setup_detector();
	if (0 < ierr) {
		nerr = 4000 + ierr;
		std::cerr << "Failed to setup detectors. Unable to calculate." << std::endl;
		goto exit;
	}

	// - Scanning
	ierr = pprm->setup_scan();
	if (0 < ierr) {
		nerr = 5000 + ierr;
		std::cerr << "Failed to setup scan. Unable to calculate." << std::endl;
		goto exit;
	}
	
exit:
	return nerr;
}


int calculate(prm_main *pprm)
{
	int nerr = 0;

	if (NULL == pprm) {
		nerr = 1;
		std::cerr << "Error (calculate): missing parameter interface." << std::endl;
		goto exit;
	}

	// calculate phase gratings
	nerr = pprm->prepare_sample_pgr();
	if (0 < nerr) goto exit;


	// determine multislice calculation scheme
	if ((pprm->gpu_id >= 0 && pprm->cpu_num == 0) ||
		(pprm->gpu_id < 0 && pprm->cpu_num == 1)) {
		// pure gpu calculation (no multi-threading)
		nerr = singlethread_stem(pprm);
	}
	else {
		// multi-threading calculation
		nerr = multithread_stem(pprm);
	}

exit:
	return nerr;
}


int output(prm_main *pprm)
{
	int ierr = 0, nerr = 0;
	unsigned int det_num = 0, idet = 0;
	std::string sfile = "";
	float normfac = 1.f;

	if (NULL == pprm) {
		nerr = 1;
		std::cerr << "Error (output): missing parameter interface" << std::endl;
		goto exit;
	}
	
	if (pprm->stem_images.detector.b_annular && (0 < pprm->stem_images.get_data_bytes())) {
		// store stem images for each detector separately, but as series over thickness
		det_num = (unsigned int)pprm->stem_images.detector.v_annular.size(); // number of detectors
		if (det_num > 0) { 
			for (idet = 0; idet < det_num; idet++) { // loop over detectors
				// write STEM image data of this detector to file
				ierr = pprm->stem_images.save(idet, "_" + pprm->stem_images.detector.v_annular[idet].name, "bin");
				if (0 < ierr) {
					nerr = 10;
					std::cerr << "Error (output): failed to write " << pprm->stem_images.detector.v_annular[idet].name  <<
						" STEM images." << std::endl;
				}
			}
		}
	}

	if (pprm->stem_pix_padif.detector.b_difpat_avg && (0 < pprm->stem_pix_padif.get_data_bytes())) {
		// get origin shift
		int orgx = 1 + (pprm->stem_pix_padif.sample.grid_nx >> 1);
		int orgy = 1 + (pprm->stem_pix_padif.sample.grid_ny >> 1);
		// get normalization factor
		normfac = 1.f / (pprm->stem_pix_padif.f_calc_weight * pprm->stem_pix_padif.f_calc_scale);
		// correct origin shift such that the zero beam is at (orgx, orgy)
		ierr = pprm->stem_pix_padif.shift_org(orgx, orgy, normfac);
		if (0 < ierr) {
			nerr = 20;
			std::cerr << "Error (output): failed to normalize and shift origin of STEM averaged diffraction pattern." << std::endl;
		}
		// store averaged stem diffraction pattern as series over thickness
		ierr = pprm->stem_pix_padif.save(0, "_padif", "bin");
		if (0 < ierr) {
			nerr = 21;
			std::cerr << "Error (output): failed to write STEM averaged diffraction pattern." << std::endl;
		}
	}

	if (pprm->stem_pix_paimg.detector.b_image_avg && (0 < pprm->stem_pix_paimg.get_data_bytes())) {
		// normalize the images
		ierr = pprm->stem_pix_paimg.normalize(pprm->stem_pix_paimg.f_calc_weight / pprm->stem_pix_paimg.f_calc_scale);
		if (0 < ierr) {
			nerr = 30;
			std::cerr << "Error (output): failed to normalize STEM averaged probe image." << std::endl;
		}
		// store averaged stem diffraction pattern as series over thickness
		ierr = pprm->stem_pix_paimg.save(0, "_paimg", "bin");
		if (0 < ierr) {
			nerr = 31;
			std::cerr << "Error (output): failed to write STEM averaged probe image." << std::endl;
		}
	}

exit:
	return nerr;
}


int finish(prm_main *pprm)
{
	int ierr = 0, nerr = 0;
	Cleanup(); // JMS clean up

	if (pprm->binteractive && pprm->bstore) {
		std::cout << std::endl;
		std::cout << "  Writing control file: " << pprm->str_ctrl_file << std::endl;
		ierr = pprm->ctrl_write_file(&pprm->str_ctrl_file);
		if (0 < ierr) {
			std::cerr << "Failed to write control file (" << pprm->str_ctrl_file << ")." << std::endl;
			nerr = 3000 + ierr;
			goto exit;
		}
	}

	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << "  Total time elapsed: " << std::setprecision(6) << 0.001f * pprm->clock.getmsec(false) << " s" << std::endl;
		std::cout << std::endl;
	}

exit:
	return nerr;
}


int main(int argc, char* argv[])
{
	int nerr = 0, ierr = 0;
	
	prm_main prm;

	// process command line options
	parseoptions(&prm, argc, argv);

	// initialize the program control
	ierr = init(&prm);
	if (0 < ierr) {
		nerr = ierr;
		goto exit;
	}

	// run parameter setup routines
	ierr = setup(&prm);
	if (0 < ierr) {
		nerr = ierr;
		goto exit;
	}

	// run calculations
	ierr = calculate(&prm);
	if (0 < ierr) {
		nerr = ierr;
		goto exit;
	}

	// output results
	ierr = output(&prm);
	if (0 < ierr) {
		nerr = ierr;
		goto exit;
	}

exit:

	ierr = finish(&prm);
	if (0 < ierr) {
		nerr = ierr;
	}

	return nerr;
}

