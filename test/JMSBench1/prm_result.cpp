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
	f_calc_scale = 1.f;
	f_calc_weight = 0.f;
	pshift_hash = NULL;
	sz_shift_hash = 0;
}


prm_result::prm_result(const prm_result &src)
{
	prm_result* psrc = const_cast <prm_result*>(&src);

	set_ctrl(src); // copy base class member

	data_type = src.data_type;
	det_type = src.det_type;
	v_dim = src.v_dim;
	str_out_file = src.str_out_file;
	f_calc_scale = src.f_calc_scale;
	f_calc_weight = src.f_calc_weight;
	
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
	if (NULL != pshift_hash) {
		free(pshift_hash);
		pshift_hash = NULL;
		sz_shift_hash = 0;
	}
	v_dim.clear();
}


int prm_result::init_buffer()
{
	size_t nd = v_dim.size();
	size_t sz = get_item_bytes();
	if (NULL != pdata) {
		free(pdata);
		pdata = NULL;
	}
	sz_data = 0;
	if (NULL != pshift_hash) {
		free(pshift_hash);
		pshift_hash = NULL;
	}
	sz_shift_hash = 0;
	
	if (nd > 0) {
		sz_data = sz;
		for (size_t i = 0; i < nd; i++) {
			sz_data *= ((size_t)v_dim[i]);
		}
		pdata = (float*)malloc(sz_data);
		if (NULL == pdata) {
			sz_data = 0;
			std::cerr << "Error: (init_buffer) failed to allocate result buffer." << std::endl;
			return 1;
		}
	}
	zero_buffer();
	return 0;
}


int prm_result::zero_buffer()
{
	if (NULL == pdata || 0 == sz_data) return 1;
	memset(pdata, 0, sz_data);
	f_calc_weight = 0.f;
	return 0;
}


int prm_result::add_buffer(float* src, float weight)
{
	if (NULL == pdata || 0 == sz_data) return 1;
	if (NULL == src) return 2;
	if (weight == 0.f) return 0;
	float* dst = (float*)pdata;
	size_t num_float = sz_data / sizeof(float);
	for (size_t i = 0; i < num_float; i++) {
		dst[i] = dst[i] + weight * src[i];
	}
	f_calc_weight += weight;
	return 0;
}

int prm_result::normalize(float weight)
{
	if (NULL == pdata || 0 == sz_data) return 1;
	if (weight == 0.f) return 2;
	float iwgt = 1.f / weight;
	float* pfdata = (float*)pdata;
	size_t num_float = sz_data / sizeof(float);
	for (size_t i = 0; i < num_float; i++) {
		pfdata[i] = pfdata[i] * iwgt;
	}
	return 0;
}

int prm_result::save(unsigned int idx_chan, std::string str_suffix, std::string str_format)
{
	int nerr = 0;
	int ifmt = 0;
	size_t dim_num = v_dim.size(); // number of data dimensions
	size_t items_num = 0; // number of items per image (pixels)
	size_t sz_out = 0; // size of output in bytes
	size_t chan_pos = 0; // offset position (item) of the channel
	unsigned int img_num = 1; // number of images per channel
	unsigned int chan_num = 1; // number of channels
	float* pchan = (float*)pdata; // pointer to the channel in pdata
	std::string str_file;
	std::ofstream ofs;

	// check for sufficient image dimensions
	if (dim_num < 2) {
		nerr = 1;
		std::cerr << "Error (save): data has insufficient number of dimensions to represent images." << std::endl;
		goto _exit;
	}

	// check state of data buffer
	if (NULL == pdata || 0 == sz_data) {
		nerr = 2;
		std::cerr << "Error (save): data buffer not allocated." << std::endl;
		goto _exit;
	}

	// calculate number of data items per image
	items_num = (size_t)v_dim[0] * v_dim[1];

	// get length of image series (thickness)
	if (dim_num > 2) {
		img_num = v_dim[2];
	}

	// get number of data channels (detectors)
	if (dim_num > 3) {
		chan_num = v_dim[3];
	}

	// check selected output channel
	if (idx_chan >= chan_num) {
		nerr = 3;
		std::cerr << "Error (save): index of selected output channel (" <<
			idx_chan << ") is too large for present number of channels (" <<
			chan_num << ")." << std::endl;
		goto _exit;
	}

	// calculate output bytes
	sz_out = get_item_bytes() * items_num * img_num; // size of the output in bytes
	chan_pos = items_num * img_num * idx_chan;
	pchan = &((float*)pdata)[chan_pos];

	// prepare output file name
	str_file = str_out_file + str_suffix;

	// make format specific settings
	if (str_format == "bin") {
		ifmt = 1;
		str_file = str_file + ".bin";
	}

	if (ifmt == 0) { // format not supported
		nerr = 999;
		std::cerr << "Error: failed to save result due to unsupported format (" << str_format << ")." << std::endl;
		goto _exit;
	}

	ofs.open(str_file, std::ios::out | std::ios::trunc | std::ios::binary);
	if (!ofs.is_open()) {
		nerr = 4;
		std::cerr << "Error (save): failed to open file [" << str_file << "]." << std::endl;
		goto _exit;
	}

	// save depending on format
	if (str_format == "bin") {
		ofs.write((char*)pchan, sz_out);
		if (ofs.fail()) {
			nerr = 5;
			std::cerr << "Error (save): failed to write data to file [" << str_file << "]." << std::endl;
			goto _exit;
		}
		if (btalk) {
			std::cout << std::endl;
			std::cout << "  Image series of channel #" << idx_chan << " written to file [" << str_file << "]" << std::endl;
			std::cout << "  - " << v_dim[0] << " x " << v_dim[1] << " x " << img_num << " x " << get_data_type_name() << std::endl;
		}
	}

_exit:
	if (ofs.is_open()) { ofs.close(); }
	return nerr;
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

int prm_result::shift_org(int orgx, int orgy, float fac)
{
	int nerr = 0;
	unsigned int dim_num = (unsigned int)v_dim.size(); // number of data dimensions
	unsigned int items_num = 0; // number of items per image (pixels)
	unsigned int float_num = 1; // number of floats per item
	unsigned int float_tot = 1; // total number of floats
	size_t sz_img = 0; // number of bytes per image
	unsigned int i = 0, h = 0, k = 0, l = 0, h0 = 0, k0 = 0, h1 = 0, k1 = 0;
	unsigned int idx = 0, idx1 = 0, idy = 0, idy1 = 0;
	unsigned int img_num = 1; // number of images
	float* pimg = (float*)pdata; // pointer to the current image
	float* ptmp = NULL; // temporary copy of an image
	unsigned int* pidx = NULL; // shift hash table (target index -> source index)

	// check for sufficient image dimensions
	if (dim_num < 2) {
		nerr = 1;
		std::cerr << "Error (shift_org): data has insufficient number of dimensions to represent images." << std::endl;
		goto _exit;
	}

	// check state of data buffer
	if (NULL == pdata || 0 == sz_data) {
		nerr = 2;
		std::cerr << "Error (shift_org): data buffer not allocated." << std::endl;
		goto _exit;
	}

	float_num = (unsigned int)(get_item_bytes() / sizeof(float)); // # floats per item
	items_num = v_dim[0] * v_dim[1]; // # items per image
	float_tot = float_num * items_num; // total # floats per image
	sz_img = sizeof(float) * (size_t)(float_tot); // # bytes per image

	// get location of the origin in the array
	h0 = (unsigned int)(orgx % (int)v_dim[0]);
	k0 = (unsigned int)(orgy % (int)v_dim[1]);
	if (h0 == 0 && k0 == 0) { // no shift needed
		goto _exit;
	}

	// allocate a local buffer to store a temporary image copy
	ptmp = (float*)malloc(sz_img);
	if (NULL == ptmp) {
		nerr = 3;
		std::cerr << "Error (shift_org): data buffer not allocated." << std::endl;
		goto _exit;
	}
	// allocate a local buffer for shift indices + setup the shifts
	if (NULL == pshift_hash) {
		sz_shift_hash = sizeof(unsigned int) * (size_t)(float_num * items_num);
		pshift_hash = (unsigned int*)malloc(sz_shift_hash);
		if (NULL == pshift_hash) {
			nerr = 4;
			std::cerr << "Error (shift_org): data buffer not allocated." << std::endl;
			goto _exit;
		}
		// the following calculation of shift assignments should happen only once per lifetime
		// of the data buffer
		for (k = 0; k < v_dim[1]; k++) { // loop over destination rows (y)
			idy = k * v_dim[0]; // target row offset
			k1 = (k + k0) % v_dim[1]; // get position + source origin wrapped to array bounds on dim 1
			idy1 = k1 * v_dim[0]; // source row offset
			for (h = 0; h < v_dim[0]; h++) { // loop through destination rows (x)
				idx = h + idy; // target column
				h1 = (h + h0) % v_dim[0]; // get position + source origin wrapped to array bounds on dim 0
				idx1 = h1 + idy1; // source column
				for (l = 0; l < float_num; l++) { // store assignment idx -> idx1 for all floats per item
					pshift_hash[l + idx] = l + idx1;
				}
			}
		}
	}

	// get length of image series (thickness)
	if (dim_num > 2) {
		img_num = v_dim[2];
	}

	// get number of data channels (detectors)
	if (dim_num > 3) {
		for (i = 3; i < dim_num; i++) {
			img_num = img_num * v_dim[i];
		}
	}

	// loop over images
	for (i = 0; i < img_num; i++) {
		pimg = &((float*)pdata)[i * float_tot]; // get pointer to current image
		memcpy(ptmp, pimg, sz_img); // copy the image
		// copy all floats shifted using the table pshift_hash
		if (fac == 1.f) { // shift without multiplication
			for (idx = 0; idx < float_tot; idx++) {
				idx1 = pshift_hash[idx];
				pimg[idx] = ptmp[idx1];
			}
		}
		else { // shift with multiplication
			for (idx = 0; idx < float_tot; idx++) {
				idx1 = pshift_hash[idx];
				pimg[idx] = ptmp[idx1] * fac;
			}
		}
	}

_exit:
	if (NULL != ptmp) { free(ptmp); }
	return nerr;
}