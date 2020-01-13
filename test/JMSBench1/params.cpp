// file: 'params.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementations of the class params.
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


#include "params.h"
#include <sstream>
#include <algorithm>
#include "string_format.h"

params::params()
{
	btalk = true;
	binteractive = false;
	bstore = false;
	in_byte_swap = false;
	ndebug = 0;
	str_ctrl_file = "control";
}


params::~params()
{
	v_str_ctrl.clear();
}

bool params::file_exists(std::string sfile)
{
	struct stat statbuf;
	return (stat(sfile.c_str(), &statbuf) == 0);
}

bool params::file_read_buf(std::ifstream *pfs, char * buf, size_t sz_item, size_t num_items)
{
	if (!pfs->is_open()) return false;
	pfs->read(buf, (std::streamsize)sz_item*num_items); // always read data from stream
	if (pfs->bad()) return false;
	if (in_byte_swap && sz_item > 1) { // support 2, 4, and 8 byte items
		// swap bytes in buf
		char tmp[32];
		char* ibuf = NULL;
		size_t i = 0, j = 0, k = 0, l = 0;
		for (i = 0; i < num_items; i++) { // loop over items
			j = i * sz_item;
			ibuf = &buf[j]; // pointer to item
			memcpy(tmp, ibuf, sz_item); // copy input bytes to tmp
			for (k = 0; k < sz_item; k++) { // swapping from tmp back to buf via ibuf
				l = sz_item - k - 1;
				ibuf[k] = tmp[l];
			}
		}
	}
	return true;
}

std::string params::generate_ser_file_name(std::string pre, unsigned int idx, unsigned int digits, std::string suf)
{
	std::stringstream ss;
	ss << std::setw(digits) << std::setfill('0') << idx;
	return pre + ss.str() + suf;
}


int params::read_param(int ipos, std::string * pstr, std::string * prm)
{
	std::size_t lipos = (std::size_t)(ipos >= 0 ? ipos : 0);
	std::size_t lcpos = lipos;
	std::size_t lmpos = 0;
	std::size_t lsep = 0;
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
	// find next separator ...
	lsep = pstr->find_first_of(str_sep, lcpos);
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
	return (int)lcpos; // return position which the possible begin of a new parameter or eos
}


int params::to_int(std::string sprm)
{
	return atoi(sprm.c_str());
}

float params::to_float(std::string sprm)
{
	return (float)atof(sprm.c_str());
}

double params::to_double(std::string sprm)
{
	return atof(sprm.c_str());
}

std::string params::to_string(std::string sprm)
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


int params::input_getline(std::istream * pfin, std::string * str)
{
	if (NULL == pfin) {
		return 1; // missing or invalid parameter 1
	}
	if (NULL == str) {
		return 2; // missing or invalid parameter 2
	}
	str->clear();
	if (!pfin->good()) {
		return 11; // input stream is not good
	}
	getline(*pfin, *str);
	if (!pfin->good()) {
		std::cerr << "Error reading command from input.\n";
		return 100;
	}
	return 0;
}

int params::file_getline(std::ifstream * pfin, std::string * str)
{
	if (NULL == pfin) {
		return 1; // missing or invalid parameter 1
	}
	if (NULL == str) {
		return 2; // missing or invalid parameter 2
	}
	str->clear();
	if (!pfin->good()) {
		return 11; // input stream is not good
	}
	getline(*pfin, *str);
	if (pfin->fail() && !pfin->eof()) {
		std::cerr << "Error reading command from input file.\n";
		return 100;
	}
	return 0;
}



int params::ctrl_getline(size_t & iline, std::string inc, std::string * str)
{
	bool bsuccess = false;
	int nerr = 0;
	if (NULL == str) {
		return 1;
	}
	if (binteractive) { // interactive line input
		std::cout << inc << " > "; // prefix
		nerr = input_getline(&(std::cin), str); // get input from std::cin to string
		if (nerr == 0) { // success
			bsuccess = true;
			v_str_ctrl.push_back(*str);
			iline++;
		}
	}
	else {
		if (iline >= 0 && iline < v_str_ctrl.size()) {
			*str = v_str_ctrl[iline];
			bsuccess = true;
			iline++;
			if (btalk) {
				std::cout << *str << std::endl;
			}
		}
	}
	if (!bsuccess && nerr == 0) {
		nerr = 300;
	}
	return nerr;
}


int params::ctrl_find_param(std::string stag, std::string * prm, int start) {
	int itag = -1, i0 = start;
	std::string lstag, lscmd;
	lstag = stag;
	std::transform(lstag.begin(), lstag.end(), lstag.begin(), ::tolower);
	int ncmd = (int)v_str_ctrl.size();
	if (ncmd > 0 && lstag.length() > 0) {
		if (i0 < 0) i0 = 0;
		if (i0 > ncmd) i0 = ncmd - 1;
		for (int i = i0; i < ncmd; i++) {
			lscmd = v_str_ctrl[i];
			std::transform(lscmd.begin(), lscmd.end(), lscmd.begin(), ::tolower);
			if (lstag == lscmd) {
				itag = i;
				break;
			}
		}
		if (itag >= 0 && itag < ncmd - 1 && prm != NULL) {
			*prm = v_str_ctrl[itag + 1];
		}
	}
	return itag;
}


int params::ctrl_init(void)
{
	int nerr = 0;
	binteractive = true;
	if (btalk) {
		std::cout << std::endl;
	}
	if (file_exists(str_ctrl_file)) {
		binteractive = false;
		if (btalk) {
			std::cout << "  Running control file: " << str_ctrl_file << std::endl;
		}
		nerr = ctrl_read_file(&str_ctrl_file);
		if (btalk) {
			if (nerr == 0) {
				std::cout << "  - number of lines: " << v_str_ctrl.size() << std::endl;
			}
			else {
				std::cerr << "Error: Failed to read lines from control file." << std::endl;
			}
		}
	}
	else {
		binteractive = true;
		if (btalk) {
			std::cout << "  Interactive control input" << std::endl;
		}
	}
	return nerr;
}

int params::ctrl_read_file(std::string * str_ctrl)
{
	int nerr = 0;
	std::ifstream fcin;
	std::string sline;
	std::string str_file = str_ctrl_file;
	if (str_ctrl) str_file = *str_ctrl; // use file name from routine parameter
	fcin.open(str_file);
	if (!fcin.is_open()) {
		std::cerr << "Error: failed to open control file " << str_file << " for reading." << std::endl;
		return 1;
	}
	v_str_ctrl.clear();
	while (nerr == 0) {
		nerr = file_getline(&fcin, &sline);
		if (nerr == 0) {
			v_str_ctrl.push_back(sline);
		}
		else {
			break;
		}
	}
	if (nerr != 0 && fcin.eof()) { // no reading error, just eof
		nerr = 0;
	}
	if (nerr != 0) {
		std::cerr << "Error reading command from file.\n ";
	}
	if (fcin.is_open()) fcin.close();
	return nerr;
}

int params::ctrl_write_file(std::string * str_ctrl)
{
	int nerr = 0;
	size_t nlines = 0, iline = 0;
	std::ofstream fcout;
	std::string sline;
	std::string str_file = str_ctrl_file;
	if (str_ctrl) str_file = *str_ctrl; // use file name from routine parameter
	fcout.open(str_file, std::ios::trunc);
	if (!fcout.is_open()) {
		std::cerr << "Error: failed to open control file " << str_file << " for writing." << std::endl;
		return 1;
	}
	nlines = v_str_ctrl.size();
	if (nlines > 0) {
		for (iline = 0; iline < nlines; iline++) {
			fcout << v_str_ctrl[iline] << std::endl;
		}
	}
	if (fcout.is_open()) fcout.close();
	if (nerr == 0 && btalk) {
		std::cout << "  - command list written to file: " << str_file << std::endl;
	}
	return nerr;
}