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
#include "Structure.h"
#include "JMultiSliceLib.h"
#include <mutex>

// sample input form option range
#define SAMPLE_INF_MIN	1
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
	std::string str_structure_file; // file name of the input atomic structure
	unsigned int structure_file_format; // index of the structure file format (see macros STRUCTURE_FILE_FORMAT_* in "Structure.h")
	std::string str_slc_file_pre; // file name prefix for slice transmission functions
	std::string str_slc_file_suf; // file name suffix for slice transmission functions
	unsigned int slc_det_per; // slice detection period (thickness steps)
	unsigned int grid_nx; // potential samples along x
	unsigned int grid_ny; // potential samples along y
	unsigned int grid_nz; // potential samples along z (for 3d potential mode only)
	float grid_a0; // potential grid size along x (nm)
	float grid_a1; // potential grid size along y (nm)
	float grid_a2; // potential grid size along z (nm)
	float grid_ekv; // electron energy from sample data (keV)
	float tilt_x; // sample tilt along x (mrad)
	float tilt_y; // sample tilt along y (mrad)
	unsigned int pot3d; // switch for projecting from a 3d potential
	unsigned int mod_atffacs; // switch the model for calculating form factors (0: Weickenmeier & Kohl)
	unsigned int mod_thermal; // switch the model of simulation for the effect of thermal motion of atoms (0: none, 1: Debye-Waller factors, 2: QEP)
	unsigned int mod_absorb; // switch hwo absorptive form factors are calculated (0: none, 1: Hashimoto, Howie & Welan, 2: Hall & Hirsch)
	float vabf; // fix absorption factor, the 
	unsigned int num_var; // number of frozen lattice variants generated per slice
	unsigned int save_sli; // switch for saving object transmission functions (>0 selects format: 1 = sli)

	float lis_mfp; // mean free path for low-loss inelastic scattering (nm)
	float lis_qe; // characteristic angle for low-loss inelastic scattering (mrad)
	float lis_qc; // critical angle for low-loss inelastic scattering (mrad)
	unsigned int lis_exc_max; // max. number of excitations for los-loss inelastic scattering

	std::vector<unsigned int> v_slc_obj; // list of object slice indices stacked to max. thickness
	std::vector<int> v_slc_det; // list of object slice detection plane indices; 0: before object, 1: after first slice, ... to slice at max. thickness
	std::vector<prm_slice> v_slc; // list of slice data, length represents the number of potential samples along z

	CStructure st_Input; // structure data read from input file
	CStructure st_Used; // structure data used for calculations

	// member functions

	// copies data from psrc to this object
	void copy_data_from(prm_sample *psrc);

	// copies setup data from psrc to this object
	// - excludes phase gratings and propagators
	void copy_setup_from(prm_sample *psrc);

	// returns the wave length for the current kinetic energy given by grid_ekv
	float get_wl(void);

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

	// runs a setup routine to load atomic structure data from files
	int load_structure(void);

	// calculates object transmission functions for slicsed of the atomic structure model 
	// - provide number of CPUs to use by cpu_num, set to 0 if no CPU should be used
	// - select the GPU ID by gpu_id, set to -1 if no GPU should be used
	int calculate_pgr(int cpu_num, int gpu_id);

	// runs a setup routine to slice the structure model stored in st_Used
	int setup_slices(void);

	// runs a setup routine to define the grid size of the calculation / phase gratings
	// a maximum spatial frequency k_max must be provided to use helper functions in interactive setup mode
	int setup_grid(float k_max = 0.f);

	// runs a setup routine to define the calculation options of potentials and object transmission functions
	int setup_slc_calc(void);

	// setup the sample thickness / slice sequence
	int setup_thickness(void);

	// setup sample tilt
	int setup_tilt(void);

	// setup low-loss inelastic scattering
	// /!\ Call this at a stage where the electron energy (grid_ekv) and the maximum sample thickness are known.
	// /!\ This activates a Monte-Carlo procedure. I should be used with a sufficient number of repeats per scan position.
	int setup_lis(void);

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

class prm_sample_calc_task
{
public:

	// standard constructor
	prm_sample_calc_task();

	// destructor
	~prm_sample_calc_task();

	// member variables

	CJEElaSca* pelsa; // address of an elastic scattering object
	CStructure* pst; // slice structure reference
	fcmplx* ppgr; // output buffer reference
	unsigned int state; // task status (0 = unsolved, 1 = reserved, 2 = calculation running, 3 = solved, 4 = cancelled, 5 = error)

};



class prm_sample_calc_queue
{
public:


	// standard constructor
	prm_sample_calc_queue();

	// destructor
	~prm_sample_calc_queue();

	// member variables

protected:

	unsigned int m_state; // queue state (0 = not ready, 1 = ready)
	prm_sample_calc_task* m_q; // task queue
	size_t m_len_q; // number of queue tasks
	size_t m_num_tasks_open; // number of open tasks (those with state = 0)
	size_t m_num_tasks_solved; // number of open tasks (those with state = 3)
	size_t m_num_tasks_failed; // number of failed tasks (those with state = 5)
	mutable std::mutex guard;

	// member functions

public:

	// initializes the task queue from a sample parameter object and an elastic scattering object
	int init(prm_sample* psample, CJEElaSca* pelsa);

	// returns true if there are no open tasks to distribute
	bool empty(void);

	// return sthe number of open tasks
	size_t open(void);

	// returns the number of solved tasks
	size_t solved(void);

	// returns the number of failed tasks
	size_t failed(void);

	// returns the total number of tasks
	size_t total(void);

	// returns a pointer to a task for a requesting thread
	size_t request(void);

	// get specific task parameters
	int get_task(size_t idx, prm_sample_calc_task &task_copy);

	// set task calculation state
	int set_task_calculate(size_t idx);

	// set task solved state
	int set_task_solved(size_t idx);

	// set task cancel state
	int set_task_cancel(size_t idx);

	// set task error state
	int set_task_error(size_t idx);
};



void WorkerCalculatePGR(int gpu_id, int cpu_thread_id, prm_sample_calc_queue* pqueue);