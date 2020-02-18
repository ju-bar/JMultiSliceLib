// file: 'prm_result.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementation of the class prm_result.
//
// Class prm_result handles one result of an image simulation.
//   Member variables are copies of the relevant calculation parameter
//   sets required to reproduce the result. The base class member
//   v_str_ctrl contains the control command sequence that lead to the
//   calculation of the results.
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

#include "prm_result.h"



prm_result::prm_result()
{
	data_type = _RESULT_DATA_TYPE_UNKNOWN;
	det_type = _JMS_DETECT_NONE;
	sz_data = 0;
	pdata = NULL;
	str_out_file = "";
}


prm_result::prm_result(const prm_result &src)
{
	prm_result* psrc = const_cast <prm_result*>(&src);

	set_ctrl(src); // copy base class member

	data_type = src.data_type;
	det_type = src.det_type;
	v_dim = src.v_dim;
	str_out_file = src.str_out_file;
	
	pdata = NULL;
	sz_data = 0;
	psrc->copy_data_to((void**)pdata, sz_data);
	if (NULL == pdata) { // failed to allocate
		v_dim.clear();
		det_type = _JMS_DETECT_NONE;
		data_type = _RESULT_DATA_TYPE_UNKNOWN;
	}

	probe.copy_data_from(&psrc->probe);
	scan.copy_data_from(&psrc->scan);
	detector.copy_data_from(&psrc->detector);
	sample.copy_data_from(&psrc->sample);
}


prm_result::~prm_result()
{
	if (NULL != pdata) { 
		free(pdata);
		pdata = NULL;
		sz_data = 0;
	}
	v_dim.clear();
}


size_t prm_result::get_item_bytes()
{
	size_t szi = 0;
	switch (data_type) {
	case _RESULT_DATA_TYPE_FLOAT:
		szi = sizeof(float);
		break;
	case _RESULT_DATA_TYPE_COMPLEX:
		szi = sizeof(fcmplx);
		break;
	default:
		szi = 1;
	}
	return szi;
}


std::string prm_result::get_data_type_name(void)
{
	std::string sname = "";

	switch (data_type) {
	case _RESULT_DATA_TYPE_FLOAT:
		sname = "32-bit float";
		break;
	case _RESULT_DATA_TYPE_COMPLEX:
		sname = "64-bit complex";
		break;
	default:
		sname = "unknown";
	}

	return sname;
}


size_t prm_result::get_data_bytes()
{
	sz_data = 0;
	if (NULL == pdata || 0 == v_dim.size() || data_type == _RESULT_DATA_TYPE_UNKNOWN) {
		return sz_data;
	}
	sz_data = sizeof(float);
	for (int i = 0; i < v_dim.size(); i++) {
		sz_data *= (size_t)v_dim[i];
	}
	switch (data_type) {
	case _RESULT_DATA_TYPE_COMPLEX:
		sz_data *= 2;
		break;
	}
	return sz_data;
}


int prm_result::copy_data_to(void ** _data, size_t &bytes)
{
	sz_data = get_data_bytes();
	bytes = sz_data;
	if (sz_data > 0) {
		*_data = malloc(bytes);
		if (NULL == *_data) {
			std::cerr << "Error: failed to allocate memory of destination buffer for result data copy." << std::endl;
			std::cerr << "- requested size [MB]: " << sz_data / 1048576 << std::endl;
			return 1;
		}
		if (NULL == memcpy(*_data, pdata, sz_data)) {
			std::cerr << "Error: failed to copy result data to destination buffer." << std::endl;
			return 2;
		}
	}
	else {
		*_data = NULL;
	}
	return 0;
}

size_t prm_result::get_pos_index(std::vector<unsigned int> v_pos)
{
	size_t pos = 0;
	unsigned int num_dim = (unsigned int)v_pos.size();
	
	if (num_dim > 0) {
		unsigned int ndpln = 1;
		for (unsigned int idim = 0; idim < num_dim; idim++) {
			pos += v_pos[idim] * ndpln;
			ndpln *= v_dim[idim];
		}
	}

	return pos;
}

float prm_result::get_dataf(std::vector<unsigned int> v_pos)
{
	float f_val = 0.f;
	if (NULL != pdata) {
		f_val = ((float*)pdata)[get_pos_index(v_pos)];
	}
	return f_val;
}

fcmplx prm_result::get_datac(std::vector<unsigned int> v_pos)
{
	fcmplx c_val = fcmplx(0.f,0.f);
	if (NULL != pdata) {
		c_val = ((fcmplx*)pdata)[get_pos_index(v_pos)];
	}
	return c_val;
}

int prm_result::write_image_series(std::string sfilename, unsigned int idx_chan)
{
	int nerr = 0;
	size_t dim_num = v_dim.size(); // number of data dimensions
	size_t items_num = 0; // number of items per image (pixels)
	size_t sz_out = 0; // size of output in bytes
	size_t chan_pos = 0; // offset position (item) of the channel
	unsigned int img_num = 1; // number of images per channel
	unsigned int chan_num = 1; // number of channels
	float* pchan = (float*)pdata; // pointer to the channel in pdata
	std::ofstream ofs;
	
	if (dim_num < 2) {
		nerr = 1;
		std::cerr << "Error (write_image_series): data has insufficient number of dimensions to represent images." << std::endl;
		goto _exit;
	}

	if (NULL == pdata || 0 == sz_data) {
		nerr = 2;
		std::cerr << "Error (write_image_series): data buffer not allocated." << std::endl;
		goto _exit;
	}

	items_num = (size_t)v_dim[0] * v_dim[1];

	if (dim_num > 2) {
		img_num = v_dim[2];
	}
	if (dim_num > 3) {
		chan_num = v_dim[3];
	}

	if (idx_chan >= chan_num) {
		nerr = 3;
		std::cerr << "Error (write_image_series): index of selected output channel (" << 
			idx_chan << ") is too large for present number of channels (" << 
			chan_num << ")." << std::endl;
		goto _exit;
	}

	sz_out = get_item_bytes() * items_num * img_num; // size of the output in bytes
	chan_pos = items_num * img_num * idx_chan;
	pchan = &((float*)pdata)[chan_pos];

	ofs.open(sfilename, std::ios::out | std::ios::trunc | std::ios::binary);
	if (ofs.is_open()) {
		ofs.write((char*)pchan, sz_out);
		if (ofs.fail()) {
			nerr = 5;
			std::cerr << "Error (write_image_series): failed to write data to file [" << sfilename << "]." << std::endl;
			goto _exit;
		}
	}
	else {
		nerr = 4;
		std::cerr << "Error (write_image_series): failed to open file [" << sfilename << "]." << std::endl;
		goto _exit;
	}

	if (btalk) {
		std::cout << std::endl;
		std::cout << "  Image series of channel #" << idx_chan << " written to file [" << sfilename << "]" << std::endl;
		std::cout << "  - " << v_dim[0] << " x " << v_dim[1] << " x " << img_num << " x " << get_data_type_name() << std::endl;
	}


_exit:
	if (ofs.is_open()) { ofs.close(); }
	return nerr;
}