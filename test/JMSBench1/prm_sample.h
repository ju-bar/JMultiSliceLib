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
#include "prm_slice.h"

// sample input form option range
#define SAMPLE_INF_MIN	2
#define SAMPLE_INF_MAX	2
// supported sample input forms
#define SAMPLE_INF_NONE	0 // no sample
#define SAMPLE_INF_CEL	1 // structure file
#define SAMPLE_INF_SLI	2 // sample data input from SLI files
// range of  input forms for Low-loss Inelastic Scattering
#define SAMPLE_LIS_INF_MIN	0 //
#define SAMPLE_LIS_INF_MAX	0 //
// supported input forms for Low-loss Inelastic Scattering
#define SAMPLE_LIS_INF_NONE	0 // none
#define SAMPLE_LIS_INF_PLA	1 // plasmonic
#define SAMPLE_LIS_INF_GEN	2 // generalized

class prm_sample :
	public params
{
public:
	
	// standard constructor
	prm_sample();

	// destructor
	~prm_sample();

	// member data

	unsigned int input_form; // switch input form of object data (structure files, phase grating files)
	std::string str_slc_file_pre; // file name prefix for slice transmission functions
	std::string str_slc_file_suf; // file name suffix for slice transmission functions
	unsigned int slc_det_per; // slice detection period (thickness steps)
	unsigned int grid_nx; // potential samples along x
	unsigned int grid_ny; // potential samples along y
	float grid_a0; // potential grid size along x (nm)
	float grid_a1; // potential grid size along y (nm)
	float grid_a2; // potential grid size along z (nm)
	float grid_ekv; // electron energy from sample data (keV)
	float tilt_x; // sample tilt along x (mrad)
	float tilt_y; // sample tilt along y (mrad)

	float lis_mfp; // mean free path for low-loss inelastic scattering (nm)
	float lis_qe; // characteristic angle for low-loss inelastic scattering (mrad)
	float lis_qc; // critical angle for low-loss inelastic scattering (mrad)
	unsigned int lis_exc_max; // max. number of excitations for los-loss inelastic scattering

	std::vector<unsigned int> v_slc_obj; // list of object slice indices stacked to max. thickness
	std::vector<int> v_slc_det; // list of object slice detection plane indices; 0: before object, 1: after first slice, ... to slice at max. thickness
	std::vector<prm_slice> v_slc; // list of slice data, length represents the number of potential samples along z

	// member functions

	// copies data from psrc to this object
	void copy_data_from(prm_sample *psrc);

	// copies setup data from psrc to this object
	// - excludes phase gratings and propagators
	void copy_setup_from(prm_sample *psrc);

	// setup the sample input form
	unsigned int setup_input_form(void);

	// returns the number of *.sli files found with the given file name prefix
	unsigned int get_sli_file_num(std::string sfile_name_prefix);

	// try to find slice files using current slice file name
	unsigned int find_sli_files(void);

	// load data from slice file headers and fill related data members
	int load_sli_file_headers(void);

	// prepare buffers and load daa from slice files
	// - assumes previous successful call of load_sli_file_headers
	// (actual loading may be postponed by load-on-demand switches, to be implemented later)
	int load_sli_file_data(void);

	// setup the sample thickness / slice sequence
	unsigned int setup_thickness(void);

	// setup sample tilt
	unsigned int setup_tilt(void);

	// setup low-loss inelastic scattering
	// /!\ Call this at a stage where the electron energy (grid_ekv) and the maximum sample thickness are known.
	unsigned int setup_lis(void);

	// get number of registered detection planes
	unsigned int get_num_slc_det(void);

	// returns the thickness of slice with index idx in the object slice stack
	float get_slc_obj_thickness(unsigned int idx);

	// returns the maximum sample thickness of the simulation
	float get_thickness_max(void);

	// returns the sample thickness upt to a given object slice index
	float get_thickness_at(unsigned int idx);

	// return the maximum number of variants over all slices
	int get_variant_num_max(void);

	// prepare propagator functions from current sample setup
	// The function requires previous setup of parameters in the JMultiSlice object addressed in the interface.
	// Propagator grid size is taken from JMS grid parameters due to possible periodic extension compared to sample data grids.
	int prepare_pro(CJMultiSlice *pjms, int whichcode);

	// setup the object slice sequence in JMS
	int setup_object_slice_seq_jms(CJMultiSlice *pjms);

	// setup phase gratings in JMS
	int setup_pgr_jms(CJMultiSlice *pjms, int whichcode);

};

