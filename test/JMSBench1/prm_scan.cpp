// file: 'prm_scan.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementation of the class prm_scan.
//
// Class prm_scan handles scanning parameters and setup routines for
// STEM image simulations.
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

#include "prm_scan.h"


prm_scan::prm_scan()
{
	nx = 1;
	ny = 1;
	num_repeat = 1;
	offset_x = 0.f;
	offset_y = 0.f;
	size_x = 0.f;
	size_y = 0.f;
	rotation = 0.f;
	sub_ix0 = 0;
	sub_iy0 = 0;
	sub_ix1 = 0;
	sub_iy1 = 0;
}


prm_scan::~prm_scan()
{
}

void prm_scan::copy_data_from(prm_scan *psrc)
{
	if (psrc) {
		nx = psrc->nx;
		ny = psrc->ny;
		num_repeat = psrc->num_repeat;
		offset_x = psrc->offset_x;
		offset_y = psrc->offset_y;
		size_x = psrc->size_x;
		size_y = psrc->size_y;
		rotation = psrc->rotation;
		sub_ix0 = psrc->sub_ix0;
		sub_iy0 = psrc->sub_iy0;
		sub_ix1 = psrc->sub_ix1;
		sub_iy1 = psrc->sub_iy1;
	}
}


float prm_scan::get_step_x(void)
{
	if (nx > 0) {
		return size_x / (float)nx;
	}
	return 0.f;
}


float prm_scan::get_step_y(void)
{
	if (nx > 0) {
		return size_y / (float)ny;
	}
	return 0.f;
}


void prm_scan::get_min_scan_samples(prm_probe *probe, unsigned int &min_x, unsigned int &min_y)
{
	float wl = probe->get_wl();
	min_x = (unsigned int)(0.004f * size_x * probe->alpha / wl + 0.5f);
	min_y = (unsigned int)(0.004f * size_y * probe->alpha / wl + 0.5f);
}


void prm_scan::print_setup()
{
	std::cout << std::endl;
	std::cout << "  Scan frame setup (position and size in nanometers):" << std::endl;
	std::cout << "      | position |     size |  samples |     step" << std::endl;
	if (nx > 0) {
		std::cout << "    x " << format("| %8.3f | %8.3f | %8d | %8.5f", offset_x, size_x, nx, size_x / (float)nx) << std::endl;
	}
	else {
		std::cout << "    x " << format("| %8.3f | no scan along x", offset_x) << std::endl;
	}
	if (ny > 0) {
		std::cout << "    y " << format("| %8.3f | %8.3f | %8d | %8.5f", offset_y, size_y, ny, size_y / (float)ny) << std::endl;
	}
	else {
		std::cout << "    y " << format("| %8.3f | no scan along y", offset_y) << std::endl;
	}
	std::cout << "    rotation towards sample box x-axis:  " << format("%8.2f", rotation) << " deg." << std::endl;
	std::cout << "    active scan range: " << format("((%d, %d), (%d, %d))", sub_ix0, sub_iy0, sub_ix1, sub_iy1) << std::endl;
	std::cout << "    number of repeats per scan position: " << format("%8d", num_repeat) << std::endl;
}

int prm_scan::setup_sampling()
{
	int ierr = 0, i = 0, j = 0;
	std::string stmp = "", sprm = "";
	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Enter the number of scan steps along the fast scan direction (x): ";
		std::cin >> nx;
		std::cout << std::endl;
		std::cout << "  Enter the number of scan steps along the fast scan direction (y): ";
		std::cin >> ny;
		v_str_ctrl.push_back("set_scan_sampling");
		v_str_ctrl.push_back(format("%i %i", nx, ny));
	}
	else {
		nx = 0; ny = 0;
		i = ctrl_find_param("set_scan_sampling", &stmp);
		if (i >= 0) {
			j = read_param(0, &stmp, &sprm);
			nx = to_int(sprm);
			if (j > 0 && j < stmp.length()) {
				read_param(j, &stmp, &sprm);
				ny = to_int(sprm);
			}
			else {
				ierr = 2;
				std::cerr << "Error: missing second parameter following control command set_scan_sampling." << std::endl;
				goto _exit;
			}
		}
		else {
			ierr = 1;
			std::cerr << "Error: missing control command set_scan_sampling." << std::endl;
			goto _exit;
		}
	}
_exit:
	return ierr;
}

int prm_scan::setup_subframe()
{
	int ierr = 0, i = 0, j = 0, ichk = 0;
	std::string stmp = "", sprm = "";
	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Do you want to define a special scan sub-frame,  <0> No  <1> Yes  ? ";
		std::cin >> ichk;
		if (1 == ichk) {
			std::cout << std::endl;
			std::cout << "  Enter the first scan column index (0 <= x0 <= " << nx - 1 << "): ";
			std::cin >> sub_ix0;
			sub_ix0 = std::min(sub_ix0, nx - 1);
			std::cout << std::endl;
			std::cout << "  Enter the first scan row index (0 <= y0 <= " << ny - 1 << "): ";
			std::cin >> sub_iy0;
			sub_iy0 = std::min(sub_iy0, ny - 1);
			std::cout << std::endl;
			std::cout << "  Enter the last scan column index (" << sub_ix0 << " <= x1 <= " << nx - 1 << "): ";
			std::cin >> sub_ix1;
			sub_ix1 = std::max(std::min(sub_ix1, nx - 1), sub_ix0);
			std::cout << std::endl;
			std::cout << "  Enter the last scan row index (" << sub_iy0 << " <= y1 <= " << ny - 1 << "): ";
			std::cin >> sub_iy1;
			sub_iy1 = std::max(std::min(sub_iy1, ny - 1), sub_iy0);
		}
		else {
			sub_ix0 = 0; sub_iy0 = 0;
			sub_ix1 = nx - 1; sub_iy1 = ny - 1;
		}
		v_str_ctrl.push_back("set_scan_sub_frame");
		v_str_ctrl.push_back(format("%i %i %i %i", sub_ix0, sub_iy0, sub_ix1, sub_iy1));
	}
	else {
		sub_ix0 = 0; sub_iy0 = 0;
		sub_ix1 = nx - 1; sub_iy1 = ny - 1;
		i = ctrl_find_param("set_scan_sub_frame", &stmp);
		if (i >= 0) { //
			// * TODO * //
			// implement a working error handling here, this is user input and may not work
			j = read_param(0, &stmp, &sprm);
			sub_ix0 = to_int(sprm);
			sub_ix0 = std::min(sub_ix0, nx - 1);
			if (j > 0 && j < stmp.length()) {
				j = read_param(j, &stmp, &sprm);
				sub_iy0 = to_int(sprm);
				sub_iy0 = std::min(sub_iy0, ny - 1);
			}
			if (j > 0 && j < stmp.length()) {
				j = read_param(j, &stmp, &sprm);
				sub_ix1 = to_int(sprm);
				sub_ix1 = std::max(std::min(sub_ix1, nx - 1), sub_ix0);
			}
			if (j > 0 && j < stmp.length()) {
				j = read_param(j, &stmp, &sprm);
				sub_iy1 = to_int(sprm);
				sub_iy1 = std::max(std::min(sub_iy1, ny - 1), sub_iy0);
			}
		}
	}
//_exit:
	return ierr;
}

int prm_scan::setup_beam_position()
{
	int ierr = 0, i = 0, j = 0;
	std::string stmp = "", sprm = "";

	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Set the beam x position in nm: ";
		std::cin >> offset_x;
		std::cout << std::endl;
		std::cout << "  Set the beam y position in nm: ";
		std::cin >> offset_y;
		//
		v_str_ctrl.push_back("set_scan_beam_pos");
		v_str_ctrl.push_back(format("%g %g", offset_x, offset_y));
	}
	else {
		offset_x = 0.f; offset_y = 0.f;
		i = ctrl_find_param("set_scan_beam_pos", &stmp);
		if (i >= 0) {
			j = read_param(0, &stmp, &sprm);
			offset_x = to_float(sprm);
			if (j > 0 && j < stmp.length()) {
				read_param(j, &stmp, &sprm);
				offset_y = to_float(sprm);
			}
			else {
				ierr = 2;
				std::cerr << "Error: missing second parameter following control command set_scan_beam_pos." << std::endl;
				goto _exit;
			}
		}
		else {
			ierr = 1;
			std::cerr << "Error: missing control command set_scan_beam_pos." << std::endl;
			goto _exit;
		}
	}

_exit:
	return ierr;
}

int prm_scan::setup_frame_size()
{
	int ierr = 0, i = 0, j = 0;
	std::string stmp = "", sprm = "";

	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Set the scan frame size along the fast scan axis (x) in nm: ";
		std::cin >> size_x;
		std::cout << std::endl;
		std::cout << "  Set the scan frame size along the slow scan axis (y) in nm: ";
		std::cin >> size_y;
		//
		v_str_ctrl.push_back("set_scan_frame_size");
		v_str_ctrl.push_back(format("%g %g", size_x, size_y));
	}
	else {
		size_x = 0.f; size_y = 0.f;
		i = ctrl_find_param("set_scan_frame_size", &stmp);
		if (i >= 0) {
			j = read_param(0, &stmp, &sprm);
			size_x = to_float(sprm);
			if (j > 0 && j < stmp.length()) {
				read_param(j, &stmp, &sprm);
				size_y = to_float(sprm);
			}
			else {
				ierr = 2;
				std::cerr << "Error: missing second parameter following control command set_scan_frame_size." << std::endl;
				goto _exit;
			}
		}
		else {
			ierr = 1;
			std::cerr << "Error: missing control command set_scan_frame_size." << std::endl;
			goto _exit;
		}
	}

_exit:
	return ierr;
}

int prm_scan::setup_frame_rotation()
{
	int ierr = 0, i = 0, j = 0;
	std::string stmp = "", sprm = "";

	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Set the scan frame rotation in degrees: ";
		std::cin >> rotation;
		//
		v_str_ctrl.push_back("set_scan_frame_rot");
		v_str_ctrl.push_back(format("%g", rotation));
	}
	else {
		rotation = 0.0f;
		i = ctrl_find_param("set_scan_frame_rot", &stmp);
		if (i >= 0) {
			j = read_param(0, &stmp, &sprm);
			rotation = to_float(sprm);
		}
	}

	return ierr;
}

int prm_scan::setup_repeats()
{
	int ierr = 0, i = 0, j = 0;
	std::string stmp = "", sprm = "";

	if (binteractive) {
		std::cout << std::endl;
		std::cout << "  Set number of repeats per scan position: ";
		std::cin >> num_repeat;
		//
		v_str_ctrl.push_back("set_num_repeat");
		v_str_ctrl.push_back(format("%d", num_repeat));
	}
	else {
		num_repeat = 1;
		i = ctrl_find_param("set_num_repeat", &stmp);
		if (i >= 0) {
			j = read_param(0, &stmp, &sprm);
			num_repeat = to_int(sprm);
			if (num_repeat < 1) num_repeat = 1;
		}
	}

	return ierr;
}
