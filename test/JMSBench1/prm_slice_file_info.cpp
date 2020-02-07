// file: 'prm_slice_file_info.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementation of the class prm_slice_file_info.
//
// Class prm_slice_file_info is used to store information from slice file
// headers.
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

#include "prm_slice_file_info.h"



prm_slice_file_info::prm_slice_file_info()
{
	in_byte_swap = false;
	ekv = 0.f;
	sx = 0.f;
	sy = 0.f;
	sz = 0.f;
	data_type = 0;
	grid_x = 0;
	grid_y = 0;
	num_var = 0;
	data_offset = 0;
	structure_offset = 0;
	sz_slc = 0;
	sz_all = 0;
	str_file_name = "";
	memset(s_slc_name, 0, (size_t)SLICENAME_MAX);
}


prm_slice_file_info::~prm_slice_file_info()
{
}


void prm_slice_file_info::clear()
{
	in_byte_swap = false;
	ekv = 0.f;
	sx = 0.f;
	sy = 0.f;
	sz = 0.f;
	data_type = 0;
	grid_x = 0;
	grid_y = 0;
	num_var = 0;
	data_offset = 0;
	structure_offset = 0;
	sz_slc = 0;
	sz_all = 0;
	str_file_name = "";
	memset(s_slc_name, 0, (size_t)SLICENAME_MAX);
}

int prm_slice_file_info::load_ems_header(std::string sfile)
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
		fs.seekg(EMS_POS_TITLE); file_read_buf(&fs, s_slc_name, 1, SLICENAME_MAX);
		fs.seekg(EMS_POS_EXTHEADER_VER); file_read_buf(&fs, (char*)&ext_hdr_ver, sizeof(unsigned int), 1);
		fs.seekg(EMS_POS_CONTENT_TYPE); file_read_buf(&fs, (char*)&data_type, sizeof(unsigned int), 1);
		switch (ext_hdr_ver) {
		case 2010112301:
			// load number of variants
			fs.seekg(EMS_POS_VARIANT_NUM); file_read_buf(&fs, (char*)&num_var, sizeof(unsigned int), 1);
			break;
		case 2012071101: //  added: nothing, just the version was increased
			// load number of variants
			fs.seekg(EMS_POS_VARIANT_NUM); file_read_buf(&fs, (char*)&num_var, sizeof(unsigned int), 1);
			// load alternative data offset in file
			fs.seekg(EMS_POS_ALT_OFFSET); file_read_buf(&fs, (char*)&utmp, sizeof(unsigned int), 1);
			data_offset = (unsigned long long)utmp;
			// load alternative table offset in file
			fs.seekg(EMS_POS_ITAB_OFFSET); file_read_buf(&fs, (char*)&utmp, sizeof(unsigned int), 1);
			structure_offset = (unsigned long long)utmp;
			// load the atom table elsewhere using structure_offset
			break;
		default:
			// set default values for the extended header data in cases
			// where no version number match is given
			num_var = 1;
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
	if (ierr == 0) {
		str_file_name = sfile;
	}
	return ierr;
}

bool prm_slice_file_info::check_ems_endian(std::ifstream* pfs)
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
