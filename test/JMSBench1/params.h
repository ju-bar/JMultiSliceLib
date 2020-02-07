// file: 'params.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the declaration of the class params.
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
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "perfclock.h"
#include "string_format.h"

class params
{
public:

	params(); // standard constructor
	params(const params &prm); // copy constructor
	~params(); // destructor

	// member parameters

	bool btalk; // flag talkative mode
	bool binteractive; // flag interactive control mode
	bool bstore; // flag for storing the command list
	bool in_byte_swap; // flag byte swap on input
	int ndebug; // debug level
	perfclock clock; // clock object for taking times

	// control interface members
	
	std::string str_ctrl_file; // control file name
	std::vector<std::string> v_str_ctrl; // list of commands read from control file

protected:

	int input_getline(std::istream * pfin, std::string * str);

	int file_getline(std::ifstream * pfin, std::string * str);

public:

	// set control interface parameters in base class members
	void set_ctrl(const params &ctrl);

	// checks whether a file specified by its file name (sfile) exists
	bool file_exists(std::string sfile);

	// reads input from file stream to a buffer with byte swap depending on member in_byte_swap
	// num_items determines the number of items of size sz_item contained by buf
	bool file_read_buf(std::ifstream *pfs, char * buf, size_t sz_item, size_t num_items);

	// generates a file name of a series with numerical index (idx)
	// concateneting prefix (pre), a string made of (digits) digits
	// represnting (idx), and a suffix (suf).
	std::string generate_ser_file_name(std::string pre, unsigned int idx, unsigned int digits, std::string suf);

	// separates a single parameter string (*prm)
	// from a list of parameters in a string (*pstr)
	// beginning with position (ipos) of (*pstr)
	// returns a position in pstr which might be the offset for the next parameter
	int read_param(int ipos, std::string * pstr, std::string * prm);

	// converts a string to int
	int to_int(std::string sprm);

	// converts a string to float
	float to_float(std::string sprm);

	// converts a string to double
	double to_double(std::string sprm);

	// converts an input string to a string value (removes quotation marks if present)
	std::string to_string(std::string sprm);

	// reads a new command or parameter line from std::cin or the control
	// file lines to (str) and increments (iline)
	// (inc) is a prefix shouted by the interactive input console
	int ctrl_getline(size_t & iline, std::string inc, std::string * str);

	// look up of stag in list of input control lines, beginning from index start
	// returns index of stag in v_str_ctrl and stores next line string to prm
	// returns a value < 0 in case of error or if stag is not found
	int ctrl_find_param(std::string stag, std::string * prm, int start = 0);

	int ctrl_init(void);

	int ctrl_read_file(std::string * str_ctrl = NULL);

	int ctrl_write_file(std::string * str_ctrl = NULL);
};

