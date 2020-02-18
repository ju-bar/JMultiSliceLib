#include "prm_slice.h"
#include "NatureConstants.h"


prm_slice::prm_slice()
{
	pro_nx = 0;
	pro_ny = 0;
	sz_pgr = 0;
	sz_pro = 0;
	pgr = NULL; // pointer to the memory block holding phase grating data (variants in sequence)
	pro = NULL; // pointer to the memory block holding propagator data
}

prm_slice::prm_slice(const prm_slice &src)
{
	pro_nx = 0;
	pro_ny = 0;
	sz_pgr = 0;
	sz_pro = 0;
	pgr = NULL;
	pro = NULL;

	prm_slice *psrc = const_cast <prm_slice*>(&src);
	header = src.header;
	if (NULL != psrc->get_pgr()) {
		alloc_pgr();
		memcpy((void*)pgr, (void*)psrc->get_pgr(), sz_pgr);
	}
	if (NULL != psrc->get_pro()) {
		size_t nx = 0, ny = 0;
		psrc->get_pro_dim(nx, ny);
		alloc_pro(nx, ny);
		memcpy((void*)pro, (void*)psrc->get_pro(), sz_pro);
	}
}

prm_slice::~prm_slice()
{
	clear();
}

void prm_slice::clear()
{
	if (NULL != pro) {
		free(pro);
		pro = NULL;
		sz_pro = 0;
	}
	if (NULL != pgr) {
		free(pgr);
		pgr = NULL;
		sz_pgr = 0;
	}
}

int prm_slice::load_ems_header(std::string sfile)
{
	return header.load_ems_header(sfile);
}

int prm_slice::load_ems_data(std::string sfile)
{
	int ierr = 0;
	std::ifstream fs;
	size_t num_floats = (size_t)(2 * header.grid_x * header.grid_y * header.num_var);
	if (num_floats == 0 || header.data_offset < (unsigned long long)EMS_POS_DATA) {
		std::cerr << "Failed to load ems data due to invalid header presets." << std::endl;
		ierr = 1;
		goto exit;
	}
	if (!file_exists(sfile)) {
		std::cerr << "Failed to load ems data, file " << sfile << " not found." << std::endl;
		ierr = 2;
		goto exit;
	}
	// allocate destination buffer
	if (0 == alloc_pgr()) {
		ierr = 3;
		goto exit;
	}
	// open the file
	fs.open(sfile, std::ios::binary);
	if (fs.is_open()) {
		// goto data offset
		fs.seekg((std::streampos)header.data_offset);
		// load data from file
		if (!file_read_buf(&fs, (char*)pgr, sizeof(float), num_floats)) {
			std::cerr << "Failed to load ems data, file " << sfile << " not found." << std::endl;
			ierr = 5;
			goto exit;
		}
	}
	else {
		std::cerr << "Failed to open file " << sfile << "." << std::endl;
		ierr = 4;
		goto exit;
	}
exit:
	if (fs.is_open()) fs.close(); // close file
	return ierr;
}

size_t prm_slice::alloc_pgr()
{
	sz_pgr = sizeof(fcmplx) * header.grid_x * header.grid_y * header.num_var;
	if (NULL != pgr) {
		free(pgr);
		pgr = NULL;
	}
	if (sz_pgr > 0) {
		pgr = (float*)malloc(sz_pgr);
		if (NULL == pgr) {
			std::cerr << "Error: failed to allocate phase gratings memory." << std::endl;
			std::cerr << "- requested size [MB]: " << sz_pgr / 1048576 << std::endl;
			sz_pgr = 0;
			goto _exit;
		}
		memset(pgr, 0, sz_pgr); // preset with 0
	}
_exit:
	return sz_pgr;
}

size_t prm_slice::alloc_pro(size_t nx, size_t ny)
{
	if (NULL != pro) {
		free(pro);
		pro = NULL;
		sz_pro = 0;
		pro_nx = 0;
		pro_ny = 0;
	}
	sz_pro = sizeof(fcmplx) * nx * ny;
	if (sz_pro > 0) {
		pro = (float*)malloc(sz_pro);
		if (NULL == pro) {
			std::cerr << "Error: failed to allocate propagator memory." << std::endl;
			std::cerr << "- requested size [MB]: " << sz_pro / 1048576 << std::endl;
			sz_pro = 0;
			goto _exit;
		}
		memset(pro, 0, sz_pro); // preset with 0
		pro_nx = nx;
		pro_ny = ny;
	}
_exit:
	return sz_pro;
}

void prm_slice::get_pro_dim(size_t &nx, size_t &ny)
{
	nx = pro_nx;
	ny = pro_ny;
}

fcmplx* prm_slice::get_pgr(size_t var)
{
	size_t npixvar = (size_t)(header.grid_x * header.grid_y);
	return &(((fcmplx*)pgr)[var * npixvar]);
}

int prm_slice::calculate_propagator(CJMultiSlice *pjms, float tilt_x, float tilt_y, int type)
{
	int nerr = 0, ierr = 0;
	int nx = 0, ny = 0;
	float ttf = (float)_R2D * 0.001f; // scales from mrad to degrees
	float otx = ttf * tilt_x, oty = ttf * tilt_y; // tilts in degree

	if (NULL == pjms) {
		nerr = 1;
		std::cerr << "Error: (calculate_propagator) called with invalid interface to JMS." << std::endl;
		goto _exit;
	}

	pjms->GetGridSize(nx, ny);
	alloc_pro((size_t)nx, (size_t)ny);
	if (sz_pro == 0) {
		nerr = 3;
		std::cerr << "Error: (calculate_propagator) failed to allocate memory." << std::endl;
		goto _exit;
	}

	ierr = pjms->CalculatePropagator(header.sz, otx, oty, (fcmplx*)pro, type);
	if (0 < ierr) {
		nerr = 4;
		std::cerr << "Error: (calculate_propagator) failed to calculate propagator." << std::endl;
		goto _exit;
	}

_exit:
	return nerr;
}

fcmplx* prm_slice::get_pro(void)
{
	return (fcmplx*)pro;
}

fcmplx prm_slice::get_pgr_at(size_t pos, size_t var)
{
	size_t npixvar = (size_t)(header.grid_x * header.grid_y);
	return ((fcmplx*)pgr)[pos + var * npixvar];
}

fcmplx prm_slice::get_pro_at(size_t pos)
{
	return ((fcmplx*)pro)[pos];
}

void prm_slice::get_grid_dim(unsigned int* dim_x, unsigned int* dim_y)
{
	*dim_x = header.grid_x;
	*dim_y = header.grid_y;
}

void prm_slice::get_grid_size(float* size_x, float* size_y)
{
	*size_x = header.sx;
	*size_y = header.sy;
}


float prm_slice::get_energy_kev(void)
{
	return header.ekv;
}

float prm_slice::get_thickness(void)
{
	return header.sz;
}

unsigned int prm_slice::get_variant_num(void)
{
	return header.num_var;
}

unsigned long long prm_slice::get_file_data_offset(void)
{
	return header.data_offset;
}

unsigned long long prm_slice::get_file_structure_offset(void)
{
	return header.structure_offset;
}