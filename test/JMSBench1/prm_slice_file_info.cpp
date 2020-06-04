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
	ver = (unsigned int)(SLICEPARAMS_EXTHDRVER);
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
	memset(s_slc_name, 0, (size_t)(SLICENAME_MAX));
	data_offset = (unsigned long long)(EMS_POS_DATA);
	structure_offset = (unsigned long long)(EMS_POS_ITAB);
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
	data_offset = (unsigned long long)(EMS_POS_DATA);
	structure_offset = (unsigned long long)(EMS_POS_ITAB);
	sz_slc = 0;
	sz_all = 0;
	str_file_name = "";
	memset(s_slc_name, 0, (size_t)SLICENAME_MAX);
}

int prm_slice_file_info::load_ems_header(std::ifstream* pfs)
{
	int nerr = 0;
	unsigned int utmp = 0;

	if (NULL == pfs) {
		std::cerr << "Error (prm_slice_file_info::load_ems_header): Invalid parameter." << std::endl;
		nerr = 1;
		goto _exit;
	}

	if (!pfs->is_open()) {
		std::cerr << "Error (prm_slice_file_info::load_ems_header): Input file not opened." << std::endl;
		nerr = 2;
		goto _exit;
	}

	data_offset = (unsigned long long)EMS_POS_DATA;
	structure_offset = (unsigned long long)EMS_POS_ITAB;
	in_byte_swap = check_ems_endian(pfs);
	pfs->seekg(EMS_POS_DATASIZE_X); file_read_buf(pfs, (char*)&grid_x, sizeof(unsigned int), 1);
	pfs->seekg(EMS_POS_DATASIZE_Y); file_read_buf(pfs, (char*)&grid_y, sizeof(unsigned int), 1);
	pfs->seekg(EMS_POS_TITLE); file_read_buf(pfs, s_slc_name, 1, SLICENAME_MAX);
	pfs->seekg(EMS_POS_EXTHEADER_VER); file_read_buf(pfs, (char*)&ver, sizeof(unsigned int), 1);
	pfs->seekg(EMS_POS_CONTENT_TYPE); file_read_buf(pfs, (char*)&data_type, sizeof(unsigned int), 1);
	switch (ver) {
	case 2010112301:
		// load number of variants
		pfs->seekg(EMS_POS_VARIANT_NUM); file_read_buf(pfs, (char*)&num_var, sizeof(unsigned int), 1);
		break;
	case 2012071101: //  added: nothing, just the version was increased
		// load number of variants
		pfs->seekg(EMS_POS_VARIANT_NUM); file_read_buf(pfs, (char*)&num_var, sizeof(unsigned int), 1);
		// load alternative data offset in file
		pfs->seekg(EMS_POS_ALT_OFFSET); file_read_buf(pfs, (char*)&utmp, sizeof(unsigned int), 1);
		data_offset = (unsigned long long)utmp;
		// load alternative table offset in file
		pfs->seekg(EMS_POS_ITAB_OFFSET); file_read_buf(pfs, (char*)&utmp, sizeof(unsigned int), 1);
		structure_offset = (unsigned long long)utmp;
		// load the atom table elsewhere using structure_offset
		break;
	case 2020052801: //  added: nothing, the version was increased due to a change in atom table format
		// load number of variants
		pfs->seekg(EMS_POS_VARIANT_NUM); file_read_buf(pfs, (char*)&num_var, sizeof(unsigned int), 1);
		// load alternative data offset in file
		pfs->seekg(EMS_POS_ALT_OFFSET); file_read_buf(pfs, (char*)&utmp, sizeof(unsigned int), 1);
		data_offset = (unsigned long long)utmp;
		// load alternative table offset in file
		pfs->seekg(EMS_POS_ITAB_OFFSET); file_read_buf(pfs, (char*)&utmp, sizeof(unsigned int), 1);
		structure_offset = (unsigned long long)utmp;
		// load the atom table elsewhere using structure_offset
		break;
	default:
		// set default values for the extended header data in cases
		// where no version number match is given
		num_var = 1;
	}
	pfs->seekg(EMS_POS_THICKNESS); file_read_buf(pfs, (char*)&sz, sizeof(float), 1);
	// structure.v_cell.z = sz;
	pfs->seekg(EMS_POS_HIGHTENSION); file_read_buf(pfs, (char*)&ekv, sizeof(float), 1);
	pfs->seekg(EMS_POS_PHYSSIZE_X); file_read_buf(pfs, (char*)&sx, sizeof(float), 1);
	// structure.v_cell.x = sx;
	pfs->seekg(EMS_POS_PHYSSIZE_Y); file_read_buf(pfs, (char*)&sy, sizeof(float), 1);
	// structure.v_cell.y = sy;
	
_exit:
	return nerr;
}

int prm_slice_file_info::save_ems_header(std::ofstream* pfs)
{
	int nerr = 0, ierr = 0;
	unsigned int utmp = 0;
	char* buf = NULL;
	size_t szbuf = (size_t)(EMS_POS_ITAB); // the starndard header information is up to the atom and inelastic transition table
	if (NULL == pfs) { nerr = 1; goto _exit; }
	if (!pfs->is_open()) { nerr = 2; goto _exit; }
	buf = new char[szbuf];
	if (NULL == buf) { nerr = 3; goto _exit; }
	memset(buf, 0, szbuf); // preset zeroes
	// prepare data in buf
	memcpy(&buf[EMS_POS_DATASIZE_X], &grid_x, sizeof(unsigned int));
	memcpy(&buf[EMS_POS_DATASIZE_Y], &grid_y, sizeof(unsigned int));
	memcpy(&buf[EMS_POS_TITLE], s_slc_name, (size_t)(SLICENAME_MAX));
	utmp = (unsigned int)(SLICEPARAMS_EXTHDRVER);
	memcpy(&buf[EMS_POS_EXTHEADER_VER], &utmp, sizeof(unsigned int));
	memcpy(&buf[EMS_POS_VARIANT_NUM], &num_var, sizeof(unsigned int));
	memcpy(&buf[EMS_POS_CONTENT_TYPE], &data_type, sizeof(unsigned int));
	utmp = (unsigned int)data_offset;
	memcpy(&buf[EMS_POS_ALT_OFFSET], &utmp, sizeof(unsigned int));
	utmp = (unsigned int)structure_offset;
	memcpy(&buf[EMS_POS_ITAB_OFFSET], &utmp, sizeof(unsigned int));
	memcpy(&buf[EMS_POS_THICKNESS], &sz, sizeof(float));
	memcpy(&buf[EMS_POS_HIGHTENSION], &ekv, sizeof(float));
	memcpy(&buf[EMS_POS_PHYSSIZE_X], &sx, sizeof(float));
	memcpy(&buf[EMS_POS_PHYSSIZE_Y], &sy, sizeof(float));
	// write buf to file
	if (!file_write_buf(pfs, buf, szbuf, 1)) { nerr = 200; goto _exit; }

_exit:
	if (NULL != buf) { delete[] buf; }
	return nerr;
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
