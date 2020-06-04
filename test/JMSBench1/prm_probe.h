// file: 'prm_probe.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the declaration of the class prm_probe.
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

#pragma once

#include "params.h"
#include "JProbeGen.h"

class prm_probe :
	public params
{
public:

	// standard constructor
	prm_probe();

	// destructor
	~prm_probe();

public:
	
	float ekv; // probe electron eneryg in keV

	float alpha; // probe convergence angle in mrad

	float alpha_edge_smooth; // probe forming aperture edge smoothness relative to alpha

	float alpha_aniso_amp; // probe forming aperture radius anisotropy relative to alpha

	float alpha_aniso_dir; // probe forming aperture radius anisotropy orientation (of large axis) relative to grid x axis in degrees

	float size; // effective probe size in nm (HWHM)

	float focus_spread; // effective focus spread in nm (1/e half width)

	unsigned int size_distribution; // source size distribution (0: point, 1: Gaussian, 2: Lorentz, 3: Disk)

protected:

	float* aberrations; // probe aberration coefficient list, coefficients (x,y)

	unsigned int num_aberrations; // number of aberrations, half length of aberrations

public:

	// Returns number of abarration coefficients described by the list aberrations in (x,y) sequence.
	// This should be at maximum half of the number of float items in the aberrations member.
	unsigned int get_num_aberrations(void);

	// calculates the electron wavelength in nm from the ekv member
	float get_wl(void);

	// returns the two coefficients (x, y) of an aberration determined by its index
	// returns an error code (0 = success)
	int get_aberration_coeff(unsigned int idx, float &x, float &y);

	// sets the two coefficients (x, y) of an aberration determined by its index
	// returns an error code (0 = success)
	int set_aberration_coeff(unsigned int idx, float x, float y);

	// copies all data from psrc to this object
	void copy_data_from(prm_probe *psrc);

	// copies the list of aberration coefficients to buffer buf
	// returns an error code (0 = success)
	// /!\ Take care to provide buf with sufficient memory to hold sizeof(float) * 2 * num_aberrations bytes.
	int copy_aberrations_to(float *buf);

	// user input and setup of electron energy (high tension)
	int setup_ht(void);

	// user input and setup of the probe-forming aperture
	int setup_aperture(void);

	// prepares the aberration list for num pair of aberration coefficients (x,y)
	// returns the number of aberrations ready to store
	unsigned int prepare_aberrations(unsigned int num);

	// prints the current aberration setup
	void print_setup_aberrations(int notation, CJProbeGen* ppg);

	// gets aberration input from user, apply either cartesian (notation==0) or polar (notation==1) form
	int input_aberration_coefficient(int idx, int notation, CJProbeGen* ppg);

	// user input and setup of aberrations
	int setup_aberrations(void);


	// Transfers values to a CJProbeParams object
	int get_jpg_params(CJProbeParams *pjpp);
};

