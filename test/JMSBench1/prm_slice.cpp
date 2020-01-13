#include "prm_slice.h"



prm_slice::prm_slice()
{
	ekv = 0.f;
	sx = 0.f;
	sy = 0.f;
	sz = 0.f;
	grid_x = 0;
	grid_y = 0;
	data_type = 0;
	var_num = 0;
	data_offset = 0;
	structure_offset = 0;
	sz_slc = 0;
	sz_all = 0;
	memset(s_slc_name, 0, (size_t)SLICENAME_MAX);
	pdata = NULL; // pointer to the memory block holding all slice data (variants in sequence)
}


prm_slice::~prm_slice()
{
	clear();
}

void prm_slice::clear()
{
	if (NULL != pdata) free(pdata); pdata = NULL;
	in_byte_swap = false;
	ekv = 0.f;
	sx = 0.f;
	sy = 0.f;
	sz = 0.f;
	data_type = 0;
	grid_x = 0;
	grid_y = 0;
	var_num = 0;
	data_offset = 0;
	structure_offset = 0;
	sz_slc = 0;
	sz_all = 0;
	memset(s_slc_name, 0, (size_t)SLICENAME_MAX);
}

int prm_slice::load_ems_header(std::string sfile)
{
	int ierr = 0;
	unsigned int ext_hdr_ver = 0;
	unsigned int data_type = 0;
	unsigned int utmp = 0;
	data_offset = (unsigned long long)EMS_POS_DATA;
	structure_offset = (unsigned long long)EMS_POS_ITAB;
	std::ifstream fs;
	if (!file_exists(sfile)) {
		std::cerr << "Failed to load ems header, file " << sfile << " not found." << std::endl;
		ierr = 1;
		goto exit;
	}
	fs.open(sfile);
	if (fs.is_open()) {
		in_byte_swap = check_ems_endian(&fs);
		fs.seekg(EMS_POS_DATASIZE_X); file_read_buf(&fs, (char*)&grid_x, sizeof(unsigned int), 1);
		fs.seekg(EMS_POS_DATASIZE_Y); file_read_buf(&fs, (char*)&grid_y, sizeof(unsigned int), 1);
		fs.seekg(EMS_POS_CONTENT_TYPE); file_read_buf(&fs, (char*)&data_type, sizeof(unsigned int), 1);
		fs.seekg(EMS_POS_EXTHEADER_VER); file_read_buf(&fs, (char*)&ext_hdr_ver, sizeof(unsigned int), 1);
		switch (ext_hdr_ver) {
		case 2010112301:
			// load number of variants
			fs.seekg(EMS_POS_VARIANT_NUM); file_read_buf(&fs, (char*)&var_num, sizeof(unsigned int), 1);
			break;
		case 2012071101: //  added: nothing, just the version was increased
			// load number of variants
			fs.seekg(EMS_POS_VARIANT_NUM); file_read_buf(&fs, (char*)&var_num, sizeof(unsigned int), 1);
			// load alternative data offset in file
			fs.seekg(EMS_POS_ALT_OFFSET); file_read_buf(&fs, (char*)&utmp, sizeof(unsigned int), 1);
			data_offset = (unsigned long long)utmp;
			// load alternative table offset in file
			fs.seekg(EMS_POS_ITAB_OFFSET); file_read_buf(&fs, (char*)&utmp, sizeof(unsigned int), 1);
			structure_offset = (unsigned long long)utmp;
			// load the atom table
			// structure.load(&fs, structure_offset);
			break;
		default:
			// set default values for the extended header data in cases
			// where no version number match is given
			var_num = 1;
		}
		fs.seekg(EMS_POS_THICKNESS); file_read_buf(&fs, (char*)&sz, sizeof(float), 1);
		// structure.v_cell.z = sz;
		fs.seekg(EMS_POS_HIGHTENSION); file_read_buf(&fs, (char*)&ekv, sizeof(float), 1);
		fs.seekg(EMS_POS_PHYSSIZE_X); file_read_buf(&fs, (char*)&sx, sizeof(float), 1);
		// structure.v_cell.x = sx;
		fs.seekg(EMS_POS_PHYSSIZE_Y); file_read_buf(&fs, (char*)&sy, sizeof(float), 1);
		// structure.v_cell.y = sy;
	}
	else {
		std::cerr << "Failed to load ems header, file " << sfile << " could not be accessed." << std::endl;
		ierr = 2;
		goto exit;
	}
	fs.close();
exit:
	return ierr;
}

bool prm_slice::check_ems_endian(std::ifstream* pfs)
{
	bool swap = false;
	__int32 ichk = 0;
	if (pfs) {
		pfs->seekg(0);
		pfs->read((char*)&ichk, sizeof(__int32));
		if (ichk < 0 || ichk > 16384) swap = true;
	}
	return swap;
}

void prm_slice::get_grid_dim(unsigned int* dim_x, unsigned int* dim_y)
{
	*dim_x = grid_x;
	*dim_y = grid_y;
}

void prm_slice::get_grid_size(float* size_x, float* size_y)
{
	*size_x = sx;
	*size_y = sy;
}


float prm_slice::get_energy_kev(void)
{
	return ekv;
}


float prm_slice::get_thickness(void)
{
	return sz;
}

unsigned int prm_slice::get_variant_num(void)
{
	return var_num;
}

unsigned long long prm_slice::get_file_data_offset(void)
{
	return data_offset;
}