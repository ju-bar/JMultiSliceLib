// file: 'prm_sample.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the declaration of the class prm_sample.
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

#pragma once
#include "params.h"

// sample input form option range
#define SAMPLE_INF_MIN	2
#define SAMPLE_INF_MAX	2
// supported sample input forms
#define SAMPLE_INF_NONE	0 // no sample
#define SAMPLE_INF_CEL	1 // structure file
#define SAMPLE_INF_SLI	2 // sample data input from SLI files

class prm_sample :
	public params
{
public:
	prm_sample();
	~prm_sample();

	unsigned int input_form;
	std::string str_slc_file_pre; // file name prefix for slice transmission functions
	std::string str_slc_file_suf; // file name suffix for slice transmission functions
	unsigned int grid_nx; // potential samples along x
	unsigned int grid_ny; // potential samples along y
	float grid_a0; // potential grid size along x (nm)
	float grid_a1; // potential grid size along y (nm)
	float grid_ekv; // electron energy from sample data (keV)

protected:
	unsigned int slc_num; // number of slice transmission functions
	unsigned int slc_num_obj; // number of object slices
	unsigned int slc_det_per; // slice detection periodicity
	unsigned int* slc_obj; // slice stacking in object (length = slc_num_obj)
	int* slc_det; // detection slice flags (length = slc_num_obj)
	unsigned int slc_num_det; // number of detection slices
	float* slc_dz; // slice thickness in nm (length = slc_num)
	unsigned int* slc_var_num; // number of variants per slice (length = slc_num)
	unsigned __int64* slc_pgr_offset; // slice file data offset (bytes, length = slc_num)

public:

	// setup the sample input form
	unsigned int setup_input_form(void);

	// get number of slices of the periodic unit
	unsigned int get_num_slc(void);

	// get number of slices of the sample
	unsigned int get_num_slc_obj(void);
	
	// try to find slice files using current slice file name
	unsigned int find_sli_files(void);

	// load data from slice file headers and fill related data members
	int load_sli_file_headers(void);

	// setup the sample thickness / slice sequence
	unsigned int setup_thickness(void);

	// get number of registered detection planes
	unsigned int get_num_slc_det(void);

	// get object slice array
	unsigned int* get_slc_obj(void);

	// get detection slice array
	int* get_slc_det(void);

};

