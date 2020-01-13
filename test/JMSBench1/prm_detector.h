// file: 'prm_detector.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains:
// * the declaration of the class prm_detector.
//   Class prm_detector handles parameters and setup routines around the
//   detection of signal for TEM image simulations.
// * the declaration of the class prm_annular.
//   Class prm_annular handles parameters and setup routines for annular
//   STEM detectors.
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

#ifndef _ANNULAR_TYPES
#define _ANNULAR_TYPES
#define ANNULAR_RING		0
#define ANNULAR_SEGMENTED	1
#endif // _ANNULAR_TYPES

class prm_annular :
	public params
{
public:

	// standard constructor
	prm_annular();

	// destructor
	~prm_annular();

public:
	std::string name; // detector name
	float beta_inner; // inner detection angle (mrad)
	float beta_outer; // outer detection angle (mrad)
	unsigned int type; // type flags
	float phi_begin; // segement azimuth begin (deg) with ANNULAR_SEGMENTED
	float phi_end; // segement azimuth end (deg) with ANNULAR_SEGMENTED

	// resets the data to default
	void reset(void);

	// returns a formatted string with annular detector information
	// call with argument true to get a table header string
	std::string get_str_info(bool header=false);

	// returns a formatted string with annular detector information
	// for control string output
	std::string get_str_params();

	// user input and setup of annular detector data
	// returns an error code (0: success)
	int setup_data(std::string str_ctrl = "");
};


class prm_detector :
	public params
{
public:

	// standard constructor
	prm_detector();

	// destructor
	~prm_detector();

public:
	
	bool b_annular; // record signal of annular detectors
	bool b_difpat; // record diffraction pattern (CBED)
	bool b_difpat_avg; // record average diffraction pattern (PACBED)
	bool b_image; // record image pattern (probe image)
	bool b_wave; // record wave function
	bool b_waveft; // record wave function in reciprocal space
	bool b_wave_avg; // record average wave function (elastic channel)
	bool b_waveft_avg; // record average wave function in reciprocal space (elastic channel)

	std::vector<prm_annular> v_annular; // list of registered annular detectors

protected:

	unsigned int num_det_annular; // number of registered annular detectors

public:

	// Returns number of registered annular detectors.
	unsigned int get_num_det_annular(void);

	// Prints an annular detector setup.
	// The setup may be provided as pointer to a std::vector<prm_annular> object
	// The current object member v_annular is printed, if the argument is NULL (default).
	void print_setup_annular(std::vector<prm_annular> *pv_ad = NULL);

	// user input and setup of annular detectors
	int setup_annular(void);

};

