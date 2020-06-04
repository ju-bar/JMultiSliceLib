// file: 'prm_probe.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementation of the class prm_probe.
//
// Class prm_probe handles parameters and setup routines around the
// TEM probe / beam and lens parameters.
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

#include "prm_probe.h"
#include "NatureConstants.h" // from JMS


prm_probe::prm_probe()
{
	ekv = 200.f;
	alpha = 25.f;
	alpha_edge_smooth = 0.03f;
	alpha_aniso_amp = 0.f;
	alpha_aniso_dir = 0.f;
	aberrations = NULL;
	focus_spread = 5.f;
	size = 0.05f;
	size_distribution = 1;
	num_aberrations = 0;
}

prm_probe::~prm_probe()
{
	if (NULL != aberrations) free(aberrations);
}


unsigned int prm_probe::prepare_aberrations(unsigned int num)
{
	unsigned int nready = 0;
	if (NULL != aberrations) {
		free(aberrations);
		aberrations = NULL;
		num_aberrations = 0;
	}
	if (num > 0) {
		size_t sz_list = sizeof(float) * num * 2;
		aberrations = (float*)malloc(sz_list);
		if (NULL != aberrations) {
			nready = num;
			num_aberrations = num;
			memset(aberrations, 0, sz_list);
		}
	}
	return nready;
}


void prm_probe::copy_data_from(prm_probe *psrc)
{
	if (psrc) {
		ekv = psrc->ekv;
		alpha = psrc->alpha;
		alpha_edge_smooth = psrc->alpha_edge_smooth;
		alpha_aniso_amp = psrc->alpha_aniso_amp;
		alpha_aniso_dir = psrc->alpha_aniso_dir;
		focus_spread = psrc->focus_spread;
		size = psrc->size;
		size_distribution = psrc->size_distribution;

		prepare_aberrations(psrc->get_num_aberrations());
		if (num_aberrations > 0 && NULL != aberrations) {
			psrc->copy_aberrations_to(aberrations);
		}
	}
}


int prm_probe::copy_aberrations_to(float *buf)
{
	int nerr = 0;
	size_t sz_list = sizeof(float) * 2 * num_aberrations;
	if (NULL == buf) {
		return 1;
	}
	if (sz_list > 0) {
		memcpy(buf, aberrations, sz_list);
	}
	return nerr;
}


unsigned int prm_probe::get_num_aberrations()
{
	return num_aberrations;
}


int prm_probe::get_aberration_coeff(unsigned int idx, float &x, float &y)
{
	if (idx < num_aberrations) {
		x = aberrations[2 * idx];
		y = aberrations[1 + 2 * idx];
	}
	else {
		return 1; // error: invalid index
	}
	return 0;
}


int prm_probe::set_aberration_coeff(unsigned int idx, float x, float y)
{
	if (idx < num_aberrations) {
		aberrations[2 * idx] = x;
		aberrations[1 + 2 * idx] = y;
	}
	else {
		return 1; // error: invalid index
	}
	return 0;
}


float prm_probe::get_wl()
{
	return (float)(_WLELKEV / sqrt(ekv*(ekv + 2.*_EEL0KEV)));
}


// user input and setup of electron energy (high tension)
int prm_probe::setup_ht(void)
{
	int ierr = 0, i = 0;
	std::string  stmp, sprm;
	float ftmp = 0.f;
	
	if (binteractive) {
_repeat_input:
		std::cout << std::endl;
		std::cout << "  Set the probe electron energy in keV: ";
		std::cin >> ekv;
		if (ekv < 10.f || ekv > 10000.f) {
			std::cout << std::endl;
			std::cerr << "Error: Invalid probe electron energy, accepted range: 10 keV - 10000 keV.";
			goto _repeat_input;
		}
		v_str_ctrl.push_back("set_probe_ekv");
		v_str_ctrl.push_back(format("%8.2f", ekv));
	}
	else {
		i = ctrl_find_param("set_probe_ekv", &stmp);
		if (i >= 0) {
			read_param(0, &stmp, &sprm);
			ekv = to_float(sprm);
		}
	}

	if (btalk || binteractive) {
		std::cout << std::endl;
		ftmp = 1000.f * get_wl();
		std::cout << "  - electron wavelength in vacuum: " << format("%6.4f", ftmp) << " pm" << std::endl;
	}

	return ierr;
}

// user input and setup of the probe-forming aperture
int prm_probe::setup_aperture(void)
{
	int ierr = 0, i = 0, j = 0, ichk = 0;
	std::string  stmp = "", sprm = "";
	float ftmp1 = 0.f;

	if (binteractive) {
	_repeat_input:
		std::cout << std::endl;
		std::cout << "  Set the convergence semi-angle of the electron probe in mrad: ";
		std::cin >> alpha;
		if (alpha < 0.001f || alpha > 200.f) {
			std::cout << std::endl;
			std::cerr << "Error: Invalid probe convergence, accepted range: 0.001 mrad - 200 mrad.";
			goto _repeat_input;
		}
		//
		std::cout << std::endl;
		std::cout << "  Specify more characteristics of the effective probe forming aperture? " << std::endl;
		std::cout << "  <0> Proceed with default  <1> Modify effective aperture shape: ";
		std::cin >> ichk;
		if (ichk == 1) {
			//
		_repeat_input_smooth:
			std::cout << std::endl;
			std::cout << "  Set the smoothness of the aperture edge relative to the convergence semi-angle: ";
			std::cin >> alpha_edge_smooth;
			if (alpha_edge_smooth < 0.0f || alpha_edge_smooth >= 1.f) {
				std::cout << std::endl;
				std::cerr << "Error: Invalid edge smoothness, accepted range: 0.0 - 1.0";
				goto _repeat_input_smooth;
			}
			//
		_repeat_input_aniso_amp:
			std::cout << std::endl;
			std::cout << "  Set the anisotropy of the aperture relative to the convergence semi-angle: ";
			std::cin >> alpha_aniso_amp;
			if (alpha_aniso_amp < 0.0f || alpha_aniso_amp >= 0.5f) {
				std::cout << std::endl;
				std::cerr << "Error: Invalid anisotropy, accepted range: 0.0 - 0.5";
				goto _repeat_input_aniso_amp;
			}
			//
			std::cout << std::endl;
			std::cout << "  Set the direction of aperture anisotropy (large axis) to the grid x-axis in degrees: ";
			std::cin >> alpha_aniso_dir;
		}
		v_str_ctrl.push_back("set_probe_convergence");
		v_str_ctrl.push_back(format("%8.2f %8.5f %8.5f %8.2f", alpha, alpha_edge_smooth, alpha_aniso_amp, alpha_aniso_dir));
	}
	else {
		i = ctrl_find_param("set_probe_convergence", &stmp);
		if (i >= 0) {
			j = read_param(0, &stmp, &sprm);
			alpha = to_float(sprm);
			if (j > 0 && j < (int)(stmp.length()-1)) {
				j = read_param(j, &stmp, &sprm);
				alpha_edge_smooth = to_float(sprm);
			}
			if (j > 0 && j < (int)(stmp.length() - 1)) {
				j = read_param(j, &stmp, &sprm);
				alpha_aniso_amp = to_float(sprm);
			}
			if (j > 0 && j < (int)(stmp.length() - 1)) {
				j = read_param(j, &stmp, &sprm);
				alpha_aniso_dir = to_float(sprm);
			}
		}
		
	}

	if (btalk || binteractive) {
		std::cout << std::endl;
		ftmp1 = 0.001f * alpha / get_wl();
		std::cout << "  - probe max. diffraction vector: " << format("%6.4f", ftmp1 ) << " nm-1" << std::endl;
		std::cout << "  - probe-forming aperture relative edge smoothness: " << format("%6.4f", alpha_edge_smooth) << std::endl;
		std::cout << "  - probe-forming aperture relative anisotropy: " << format("%6.4f   at %7.2f deg", alpha_aniso_amp, alpha_aniso_dir) << std::endl;
	}

	return ierr;
}

void prm_probe::print_setup_aberrations(int notation, CJProbeGen* ppg)
{
	if (NULL == ppg) return;
	int i = 0, n = 0, l = 0;
	int inot = notation;
	float rad2deg = 45.f / atanf(1.f);
	if (inot != 1) inot = 0;
	std::string ssym = "";
	std::cout << std::endl;
	if (inot == 0) {
		std::cout << "    Index | Symbol |  x (nm)      |  y (nm)" << std::endl;
	}
	else {
		std::cout << "    Index | Symbol |  abs. (nm)   |  angle (deg)" << std::endl;
	}
	for (i = 0; i < (int)num_aberrations; i++) {
		ppg->GetAberrationSymbol(i, &ssym);
		n = ppg->GetAberrationTermOrder(i);
		l = ppg->GetAberrationTermSymmetry(i);
		if (inot == 0) { // cartesian
			if (l == 0) { // rotational symmetric aberration (only 1 coefficient)
				std::cout << format("    <%2d>  | ", i + 1) << ssym << format("    | %12.5E |", aberrations[2 * i]) << std::endl;
			} 
			else { // lower symmetry (2 coefficients)
				std::cout << format("    <%2d>  | ", i + 1) << ssym << format("    | %12.5E | %12.5E", aberrations[2 * i], aberrations[1 + 2 * i]) << std::endl;
			}
		}
		else { // polar
			float a = 0.f, p = 0.f;
			if (l == 0) { // rotational symmetric aberration (only 1 coefficient, no polar form)
				std::cout << format("    <%2d>  | ", i + 1) << ssym << format("    | %12.5E |", aberrations[2 * i]) << std::endl;
			}
			else { // lower symmetry (2 coefficients)
				a = sqrtf(aberrations[2 * i] * aberrations[2 * i] + aberrations[1 + 2 * i] * aberrations[1 + 2 * i]);
				p = rad2deg * atan2f(aberrations[1 + 2 * i], aberrations[2 * i]);
				std::cout << format("    <%2d>  | ", i + 1) << ssym << format("    | %12.5E | %8.2f", a, p) << std::endl;
			}

		}
	}
}

int prm_probe::input_aberration_coefficient(int idx, int notation, CJProbeGen* ppg)
{
	int nerr = 0;
	int inot = notation, l = 0;
	float ax = 0.f, ay = 0.f, aa = 0.f, ap = 0.f;
	float ox = 0.f, oy = 0.f, oa = 0.f, op = 0.f;
	float deg2rad = atanf(1.f) / 45.f;
	std::string sabr = "", snum = "";
	if (inot != 1) inot = 0;
	if (idx < 0 || idx >= (int)num_aberrations) {
		nerr = 1;
		std::cerr << "Error (input_aberration_coefficient): invalid aberration index (" << idx + 1 << ")." << std::endl;
		goto _exit;
	}
	if (NULL == ppg) {
		nerr = 2;
		std::cerr << "Error (input_aberration_coefficient): invalid pointer to CJProbeGen object." << std::endl;
		goto _exit;
	}
	l = ppg->GetAberrationTermSymmetry(idx);
	ppg->GetAberrationName(idx, &sabr);
	ox = aberrations[2 * idx]; oy = aberrations[1 + 2 * idx];
	oa = sqrtf(ox * ox + oy * oy);
	if (l > 0) op = atan2(oy, ox) / deg2rad;
	std::cout << std::endl;
	if (inot == 0) { // cartesian
		if (l > 0) {
			snum = format("(%12.5E, %12.5E) nm", ox, oy);
		}
		else {
			snum = format("%12.5E nm", ox);
		}
		std::cout << "    - Current " << sabr << ": " << snum << std::endl;
		std::cout << "      Enter new x component (nm): ";
		std::cin >> ax;
		aberrations[2 * idx] = ax;
		if (l > 0) {
			std::cout << "      Enter new y component (nm): ";
			std::cin >> ay;
			aberrations[1 + 2 * idx] = ay;
		}
		else {
			aberrations[1 + 2 * idx] = 0.f;
		}
		
	}
	else {
		if (l > 0) {
			snum = format("%12.5E nm at %8.2f deg", oa, op);
		}
		else {
			snum = format("%12.5E nm", ox);
		}
		std::cout << "    - Current " << sabr << ": " << snum << std::endl;
		std::cout << "      Enter new abs. value (nm): ";
		std::cin >> aa;
		if (l > 0) {
			std::cout << "      Enter new angle (deg): ";
			std::cin >> ap;
			ax = fabs(aa) * cosf(ap * deg2rad);
			ay = fabs(aa) * sinf(ap * deg2rad);
			aberrations[2 * idx] = ax;
			aberrations[1 + 2 * idx] = ay;
		}
		else {
			ax = aa;
			aberrations[2 * idx] = ax;
			aberrations[1 + 2 * idx] = 0.f;
		}
		
	}
_exit:
	return nerr;
}

// user input and setup of aberrations
int prm_probe::setup_aberrations(void)
{
	int ierr = 0, nerr = 0;
	int ichk = 0, inot = 0, idx = 0, i = 0, j = 0, l = 0;
	std::string  stmp, sprm;
	
	// use a CJProbeGen object for information around aberrations
	CJProbeGen pg;
	prepare_aberrations((unsigned int)pg.GetAberrationNum());

	if (btalk || binteractive) {
		std::cout << std::endl;
		std::cout << "  - probe aberration setup ... " << std::endl;
	}

	if (binteractive) { // interactive aberration setup
		std::cout << std::endl;
		std::cout << "    Choose probe aberration notation:  <0> cartesian (x,y),  <1> polar (abs, angle): ";
		std::cin >> inot;
		if (inot != 1) inot = 0;
	_repeat_input: // start over from here to make more changes to the aberration setup
		print_setup_aberrations(inot, &pg); // print the loaded aberration setup
		// allow to modify the setup
		std::cout << std::endl;
		std::cout << "    <0> Accept, or choose aberration from list: ";
		std::cin >> ichk;
		if (ichk > 0) {
			idx = ichk - 1;
			ierr = input_aberration_coefficient(idx, inot, &pg);
			if (0 < ierr) {
				std::cerr << "Error (setup_aberrations): aberration input failed." << std::endl;
			}
			goto _repeat_input;
		}
		if (num_aberrations > 0) { // write the setup to the list of control lines
			for (i = 0; i < (int)num_aberrations; i++) {
				v_str_ctrl.push_back("set_probe_aberration");
				pg.GetAberrationName(idx, &stmp);
				v_str_ctrl.push_back(format("%d %13.5E %13.5E", ichk, aberrations[2*idx], aberrations[1+2*idx]) + " (" + stmp + ")");
			}
		}
	}
	else { // read probe aberration setup form
		i = 0; j = 0;
		while (i >= 0) {
			i = ctrl_find_param("set_probe_aberration", &stmp, j);
			if (i >= 0) {
				int ipos = 0;
				// get aberration index (+1)
				ipos = read_param(ipos, &stmp, &sprm);
				if (ipos < 0) { // parsing error
					std::cerr << "Error (setup_aberrations): Failed to parse probe aberration index from line #" << i + 1 << "(" << stmp << ")" << std::endl;
					j = i + 2;
					continue;
				} 
				ichk = to_int(sprm);
				idx = ichk - 1;
				if (idx < 0 || idx > (int)num_aberrations) { // parsing error
					std::cerr << "Error (setup_aberrations): Invalid probe aberration index in line #" << i + 1 << "(" << stmp << ")" << std::endl;
					j = i + 2;
					continue;
				}
				// get aberration x component (nm)
				ipos = read_param(ipos, &stmp, &sprm);
				if (ipos < 0) { // parsing error
					std::cerr << "Error (setup_aberrations): Failed to parse probe aberration x component from line #" << i + 1 << "(" << stmp << ")" << std::endl;
					j = i + 2;
					continue;
				}
				aberrations[2 * idx] = to_float(sprm);
				l = pg.GetAberrationTermSymmetry(idx); // get aberration term rotational symmetry index
				if (l > 0) {
					// get aberration y component (nm)
					ipos = read_param(ipos, &stmp, &sprm);
					if (ipos < 0) { // parsing error
						std::cerr << "Error (setup_aberrations): Failed to parse probe aberration y component from line #" << i + 1 << "(" << stmp << ")" << std::endl;
						j = i + 2;
						continue;
					}
					aberrations[1 + 2 * idx] = to_float(sprm);
				}
				//
				j = i + 2;
			}
		}
		if (btalk) {
			print_setup_aberrations(0, &pg); // print the loaded aberration setup
		}
	}
	
	return nerr;
}

int prm_probe::get_jpg_params(CJProbeParams *pjpp)
{
	int ierr = 0, nab = 0;
	if (NULL == pjpp) {
		ierr = 1;
		std::cerr << "Error: (set_jpg_params) invalid argument." << std::endl;
		goto _exit;
	}
	pjpp->m_wl = get_wl();
	pjpp->m_alpha = alpha;
	pjpp->m_alpha_adir = alpha_aniso_dir;
	pjpp->m_alpha_asym = alpha_aniso_amp;
	pjpp->m_alpha_rs = alpha_edge_smooth;
	pjpp->m_alpha_x0 = 0.f;
	pjpp->m_alpha_y0 = 0.f;
	pjpp->m_btx = 0.f;
	pjpp->m_bty = 0.f;
	pjpp->m_fspread_kernel_samples = 7;
	pjpp->m_fspread_kernel_width = 2.f;
	pjpp->m_fspread_width = focus_spread;
	pjpp->m_source_shape = (int)size_distribution;
	pjpp->m_source_width = size;
	nab = std::min((int)num_aberrations, pjpp->GetAberrationNum());
	if (nab > 0) {
		for (int i = 0; i < nab; i++) {
			pjpp->m_abrr_coeff[2 * i] = aberrations[2 * i];
			pjpp->m_abrr_coeff[1 + 2 * i] = aberrations[1 + 2 * i];
		}
	}

_exit:
	return ierr;
}