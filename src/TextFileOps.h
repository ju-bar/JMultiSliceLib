// file: 'TextFileOps.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the declarations of some functions operating on
// text file I/O and parsing.
//
/* -----------------------------------------------------------------------

	This file is part of JMultiSlice.

	JMultiSlice is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	JMultiSlice is distributed in the hope that it will be useful,
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
#include "string_format.h"

// checks whether a file specified by its file name (sfile) exists
bool file_exists(std::string sfile);

// separates a single parameter string (*prm)
// from a list of parameters in a string (*pstr)
// beginning with position (ipos) of (*pstr)
// returns a position in pstr which might be the offset for the next parameter
int read_param(int ipos, std::string* pstr, std::string* prm);

// converts a string to int
int to_int(std::string sprm);

// converts a string to float
float to_float(std::string sprm);

// converts a string to double
double to_double(std::string sprm);

// converts an input string to a string value (removes quotation marks if present)
std::string to_string(std::string sprm);

// reads lines of text from file into a std::vector<std::string> object
int read_textfile_lines(std::string filename, std::vector<std::string>& v_lines);