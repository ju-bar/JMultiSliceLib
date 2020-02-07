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
		f_val = pdata[get_pos_index(v_pos)];
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