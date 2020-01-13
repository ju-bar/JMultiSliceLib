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

	float size; // effective probe size in nm (HWHM)

	float focus_spread; // effective focus spread in nm (1/e half width)

	unsigned int size_distribution; // source size distribution (0: point, 1: Gaussian, 2: Lorentz, 3: Disk)

	float* aberrations; // probe aberration coefficient list, coefficients (x,y)

protected:

	unsigned int num_aberrations; // number of aberrations, half length of aberrations

public:

	// Returns number of abarration coefficients described by the list aberrations in (x,y) sequence.
	// This should be at maximum half of the number of float items in the aberrations member.
	unsigned int get_num_aberrations(void);

	// calculates the electron wavelength in nm from the ekv member
	float get_wl(void);

	// user input and setup of electron energy (high tension)
	int setup_ht(void);

	// user input and setup of the probe-forming aperture
	int setup_aperture(void);

	// user input and setup of aberrations
	int setup_aberrations(void);
};

