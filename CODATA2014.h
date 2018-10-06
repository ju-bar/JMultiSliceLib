//
// C++ header file: CODATA2014.h
// Declaration of physical constants based on published CODATA 2014
// Source: <https://doi.org/10.5281/zenodo.22826>
//
// Copyright (C) 2018 Juri Barthel
// Email: juribarthel@gmail.com
//
/*
	This program is free software : you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.If not, see <https://www.gnu.org/licenses/>
*/
//
//
//
#pragma once
//
#ifndef __CODATA__
#define __CODATA__
#define __CODATA_YEAR__			2014
//
// Fundamental constants in physics and chemistry
// Exact constants

#define _C0						299792458. // vacuum speed of light [ m s^(-1) ]
#define _MU0					1.2566370614359172953850573533118E-6 // magnetic constant = 4*Pi * 10^(-7) [ N A^(-1) ]
#define _EPS0					8.8541878176203898505365630317108E-12 // electric constant = 1/(_MU0 _C0^2) [ F m^(-1) ]

// Derived constants

#define _HPL					6.626070040E-34	// Planck's constant, error (81) in last digits [ J s ]
#define _HPLEV					4.135667662E-15	// Planck's constant, error (25) in last digits [ eV s ]
#define _HBAR					1.054571800E-34	// reduced Planck's constant, error (13) in last digits [ J s ]
#define _HBAREV					6.582119514E-16	// reduced Planck's constant, error (40) in last digits [ eV s ]

#define _QEL					1.6021766208E-19 // elementary charge, error (98) in last digits [ C ]

#define _ME0					9.10938356E-31 // electron rest mass, error (11) in last digits [ kg ]

//
#endif
