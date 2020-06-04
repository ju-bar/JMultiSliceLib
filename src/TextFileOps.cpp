// file: 'TextFileOps.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the defenitions of some functions operating on
// text file I/O and parsing.
//
/* -----------------------------------------------------------------------

	This file is part of JMultiSlice.

	JMultiSlice is free software: you can redistribute it and/or modify
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

#include "TextFileOps.h"

bool file_exists(std::string sfile)
{
	struct stat statbuf;
	return (stat(sfile.c_str(), &statbuf) == 0);
}

int read_param(int ipos, std::string* pstr, std::string* prm)
{
	size_t lipos = (std::size_t)(ipos >= 0 ? ipos : 0);
	size_t lcpos = lipos;
	size_t lmpos = 0;
	size_t lsep = 0;
	std::string str_sep = ", ";
	if (NULL == pstr) {
		return -1;
	}
	if (NULL == prm) {
		return -2;
	}
	lmpos = pstr->size();
	if (lcpos >= lmpos) {
		*prm = "";
		return (int)lcpos;
	}
	// find first non-separator character
	lipos = pstr->find_first_not_of(str_sep, lipos);
	// find next separator ...
	lsep = pstr->find_first_of(str_sep, lipos);
	if (lsep == pstr->npos) {
		lcpos = lmpos;
	}
	else {
		lcpos = lsep;
	}
	// lcpos is now the position of the next separator or eos
	if (lcpos - lipos > 0) { // extract sub-string
		*prm = pstr->substr(lipos, lcpos - lipos);
	}
	// find next non-separator character or eos
	if (lcpos < lmpos) {
		lsep = pstr->find_first_not_of(str_sep, lcpos);
		if (lsep == pstr->npos) { // eos
			lcpos = lmpos;
		}
		else { // next param
			lcpos = lsep;
		}
	}
	return (int)lcpos; // return position which is the possible begin of a new parameter or eos
}


int to_int(std::string sprm)
{
	return atoi(sprm.c_str());
}

float to_float(std::string sprm)
{
	return (float)atof(sprm.c_str());
}

double to_double(std::string sprm)
{
	return atof(sprm.c_str());
}

std::string to_string(std::string sprm)
{
	size_t np = std::string::npos;
	size_t i1 = np, i2 = np;
	std::string sqm = "";
	std::string stmp = "";
	i1 = sprm.find_first_of("'");
	i2 = sprm.find_last_of("'");
	if (i1 == np || i2 == np) {
		i1 = sprm.find_first_of('"');
		i2 = sprm.find_last_of('"');
	}
	if (i1 == np || i2 == np) {
		stmp = sprm;
	}
	else {
		stmp = sprm.substr((size_t)(i1 + 1), (size_t)(i2 - i1 - 1));
	}
	return stmp;
}

int read_textfile_lines(std::string filename, std::vector<std::string>& v_lines)
{
	int nerr = 0;
	std::string sline = "";
	std::ifstream fcin; // file stream
	fcin.open(filename);  // open the file
	if (!fcin.is_open()) {
		nerr = 1;
		goto _exit;
	}
	// read all lines of the file
	v_lines.clear();
	while (fcin.good()) { // file state is still good and no end of structure character appeared
		getline(fcin, sline);
		v_lines.push_back(sline);
	}
	fcin.close(); // close the file
_exit:
	return nerr;
}