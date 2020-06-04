// file: 'prm_result.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains:
// * the declaration of the class prm_result.
//   Class prm_result handles one result of an image simulation.
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

#pragma once

#include "params.h"
#include "prm_probe.h"
#include "prm_sample.h"
#include "prm_detector.h"
#include "prm_scan.h"
#include "JMultiSlice.h"

#define _RESULT_DATA_TYPE_UNKNOWN	0
#define _RESULT_DATA_TYPE_FLOAT		1
#define _RESULT_DATA_TYPE_COMPLEX	2

class prm_result :
	public params
{
public:
	
	// delauft constructor
	prm_result();

	// copy constructor
	prm_result(const prm_result &src);

	// destructor
	~prm_result();

public:
	// member data

	unsigned int data_type; // data type flag (1 = 32-bit float, 2 = 64-bit complex)
	unsigned int det_type; // copy of the _JMS_DETECT_* type used to calculate the result
	std::vector<unsigned int> v_dim; // data dimensions

	float f_calc_scale; // scaling factor of the data due to calculation properties
	float f_calc_weight; // weight of the data due to calculation repeats
	
	size_t sz_data; // size of the data buffer, modify with care
	void *pdata; // data buffer, modify with care

	unsigned int *pshift_hash; // origin shift hash table (allocated on demand and may be re-used with care)
	size_t sz_shift_hash; // size of origin shift hash table

	std::string str_out_file; // file name used to save the data

	prm_probe probe;
	prm_sample sample;
	prm_scan scan;
	prm_detector detector;

protected:

	mutable std::mutex guard; // threading guard used for some data accumulation functions

	// member functions
public:

	// Initializes the result data buffer from object parameters
	int init_buffer(void);

	// Resets the content of the result buffer to zero
	int zero_buffer(void);

	// Adds content of another buffer using an accumulation weight using float arithmetics.
	// Also increments the f_calc_weight member by weight.
	// Assumes identical size of the input buffer src and pdata
	int add_buffer(float* src, float weight = 1.f);

	// normalizes the data by dividing the data with a given weight
	int normalize(float weight);

	// Stores data of a selected channel to file.
	// Channels are assumed to be the fourth index in v_dim, if existing.
	// The file name is made from str_out_file and str_suffix is appended
	// The data is stored in a format selected by str_format
	// Currently supported formats:
	// - "bin" = raw binary file, data is dumped as is
	int save(unsigned int idx_chan, std::string str_suffix, std::string str_format);
	
	// Returns the number of bytes per data item as determined by data_type
	size_t get_item_bytes(void);

	// Returns a string label for teh current data_type
	std::string get_data_type_name(void);

	// Returns the number of bytes representing the data = buffer size (pdata)
	size_t get_data_bytes(void); 

	// generates a copy of the data in a new memory block
	// returns the address of the memory to _data and the size in bytes to bytes
	int copy_data_to(void ** _data, size_t &bytes);

	// returns the index for a n-dimensional position v_pos in the result data
	size_t get_pos_index(std::vector<unsigned int> v_pos);

	// returns a float element of result data
	// Warning: The routine works without controlling consistency with object data and assumes pdata is allocated.
	float get_dataf(std::vector<unsigned int> v_pos);

	// returns a complex element of result data
	// Warning: The routine works without controlling consistency with object data and assumes pdata is allocated.
	fcmplx get_datac(std::vector<unsigned int> v_pos);

	// Rigidly shifts content of all images (first two dimensions of pdata),
	// such that the pixel specified by (orgx, orgy) is on (0, 0).
	// Applies periodic boundary conditions.
	// This function may also be used to normalize the array with a factor fac.
	// Returns an error code in case of failure.
	int shift_org(int orgx, int orgy, float fac = 1.f);

};

