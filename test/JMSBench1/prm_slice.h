#pragma once
#include "fcomplex.h"
#include "params.h"
#include "prm_slice_file_info.h"
#include "JMultiSlice.h"
#include "Structure.h"

class prm_slice :
	public params
{
public:
	
	// standard constructor
	prm_slice();

	// copy constructor
	prm_slice(const prm_slice &src);

	// destructor
	~prm_slice();

	// member variables
public:

	prm_slice_file_info header; // header
	CStructure st; // structure data of this slice (atom types and positions)

protected:
	
	size_t pro_nx; // propagator samples along x
	size_t pro_ny; // propagator samples along y
	size_t sz_pgr; // size of buffer pgr
	size_t sz_pro; // size of buffer pro
	float* pgr; // pointer to the memory block holding all phase-grating data (variants in sequence)
	float* pro; // pointer to the memory block holding the propagator
	

	// member functions
public:

	// clear all data and free memory allocations
	void clear();

	// load ems slice file header and store information in members
	int load_ems_header(std::string sfile);

	// load ems slice data from file sfile (assuming header information has been loaded previously)
	int load_ems_data(std::string sfile);

	// load ems slice atomic structure table from an open input file stream
	int load_ems_itab(std::ifstream* pfs, unsigned int header_version);

	// saves ems slice information to a file (header, structure, and data)
	int save_ems(std::string sfile);

	// saves the structure table to an output file stream
	int save_ems_itab(std::ofstream* pfs);

	// data interface

	// allocates new memory for phase grating data based on information in header
	// returns number of allocated bytes
	size_t alloc_pgr(void);

	// allocates new memory for phase grating data based on information in header
	// returns number allocated bytes
	size_t alloc_pro(size_t nx, size_t ny);

	// returns the pointer to a phase grating (handle with care, never free or malloc this)
	fcmplx* get_pgr(size_t var = 0);

	// returns the pointer to the propagator data (handle with care, never free or malloc this)
	fcmplx* get_pro(void);

	// returns the 2d grid dimensions of the propagator
	void get_pro_dim(size_t &nx, size_t &ny);

	// calculates a propagator using current slice parameters
	// - provide the address of a CJMultiSlice object with energy and grid size set up
	// - tilt_x and tilt_y are sampel tilts to be simulated by tilted propagation (mrad)
	// - type is a flag switching: 0: geometric slice plane propagation, 1: Fresnel propagation
	int calculate_propagator(CJMultiSlice *pjms, float tilt_x, float tilt_y, int type = 0);
	
	// returns an value of the phase grating data at position pos of variant var
	fcmplx get_pgr_at(size_t pos, size_t var);

	// returns an value of the propagator at position pos
	fcmplx get_pro_at(size_t pos);

	// returns the 2d slice grid dimension
	void get_grid_dim(unsigned int* dim_x, unsigned int* dim_y);

	// returns the represented size of the grid in nm
	void get_grid_size(float* size_x, float* size_y);

	// returns the electron energy related to the slice data in keV
	float get_energy_kev(void);

	// returns the slice thickness
	float get_thickness(void);

	// return the number of variants of the slice data
	unsigned int get_variant_num(void);

	// returns the offset byte of slice data in the slice file
	unsigned long long get_file_data_offset(unsigned int header_version);

	// returns the offset byte of structure data in the slice file
	unsigned long long get_file_structure_offset(void);

protected:

	// returns the number of bytes required for the structure data in the slice file
	unsigned long long get_file_structure_size(unsigned int header_version);

};

