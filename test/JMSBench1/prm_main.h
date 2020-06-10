// file: 'prm_main.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the declaration of the class prm_main.
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

#pragma once

#include "params.h"
#include "prm_result.h"
#include "JMultiSliceLib.h"

#ifndef __PRM_MAIN_H
#define __PRM_MAIN_H

#define _PRM_TEM_MODE_NUM		2 // number of tem application modes
#define _PRM_TEM_MODE_CTEM		0 // conventional TEM mode with parallel beam illumination
#define _PRM_TEM_MODE_STEM		1 // scanning TEM mode with convergent beam illumination

class prm_main :
	public params
{
public:
	prm_main();
	~prm_main();

	// Member variables
public:	
	int gpu_id; // GPU id, >=0 is valid
	int cpu_num; // CPU number of threads
	int tem_mode; // TEM application mode (0: CTEM, 1: STEM)
	unsigned int cpu_num_max; // maximum number of supported CPU threads
	std::string str_out_file; // output file name prefix
	prm_sample sample; // sample data
	prm_probe probe; // TEM probe data
	prm_detector detector; // detector data
	prm_scan scan; // scan data

	prm_result stem_images; // STEM imaging results
	prm_result stem_pix_dif; // STEM pixelated diffraction patterns
	prm_result stem_pix_padif; // STEM pixelated averaged diffraction patterns
	prm_result stem_pix_img; // STEM pixelated probe images
	prm_result stem_pix_paimg; // STEM pixelated averaged probe images
	prm_result stem_images_ela; // STEM imaging results (elastic channel)

	// Methods

	// user input and setup of CUDA GPU device
	int setup_cuda_device(void);

	// print CUDA GPU device info
	void print_cuda_device_info(int id);
	
	// user input and setup of multiple CPU threads
	int setup_cpu(void);

	// user input and setup of output file name
	int setup_file_output(void);

	// user input and setup of sample data
	int setup_sample(void);

	// user input and setup of probe data
	int setup_probe(void);

	// user input and setup of detector data
	int setup_detector(void);

	// user input and setup of probe scanning
	int setup_scan(void);

	// user input and setup of QEP options
	int setup_qepopt(void);

	// user input and setup of low-loss inelastic scattering (run after setup_sample and setup_scan)
	int setup_lis(void);

	// prepare transmission functions (phase gratings) from current sample setup
	int prepare_sample_pgr(void);

	// prepare result parameter members to recieve data and store backups of the calculation setup
	int prepare_result_params(void);
};

#endif // ifndef __PRM_MAIN_H