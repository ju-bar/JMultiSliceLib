// file: 'prm_scan.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains:
// * the declaration of the class prm_scan.
//   Class prm_scan handles scanning parameters and setup routines for
//   STEM image simulations.
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

class prm_scan :
	public params
{
public:

	// standard constructor
	prm_scan();

	// destructor
	~prm_scan();

public:

	float size_x; // size of the scan frame in nm, along fast scan direction
	float size_y; // size of the scan frame in nm, along slow scan direction
	float offset_x; // horizontal offset position of the scan frame in the calculation box in nm
	float offset_y; // vertical offset position of the scan frame in the calculation box in nm
	float rotation; // rotation of the scan frame fast axis to the calculation box x axis in rad
	unsigned int nx; // number of samples taken along the fast scan direction
	unsigned int ny; // number of samples taken along the slow scan direction
	unsigned int num_repeat; // number of calculation repeats per scan position

public:

	// copies all data from psrc to this object
	void copy_data_from(prm_scan *psrc);

	// returns the scan step along the fast scan direction in nm
	float get_step_x(void);

	// returns the scan step along the slow scan direction in nm
	float get_step_y(void);

	// returns the minimum scan samples for a given probe setup
	// - assumes periodic boundary conditions for the scan frame
	void get_min_scan_samples(prm_probe *probe, unsigned int &min_x, unsigned int &min_y);

	// user input and setup of the scan sampling
	int setup_sampling(void);

	// user input and setup of the beam position (scan frame offset)
	int setup_beam_position(void);
	
	// user input and setup of the scan frame size
	int setup_frame_size(void);

	// user input and setup of the scan frame rotation
	int setup_frame_rotation(void);

	// user input and setup of the scan frame rotation
	int setup_repeats(void);

	// prints the current scan frame setup
	void print_setup(void);

};

