// *********************************************************************
//
// file "wekoscat.h"
//
// *********************************************************************
//
// Implementation of Weickenmeier & Kohl electron atomic form factors
//
// Acta Cryst.A 47 (1991) 590 - 597
// by J.Barthel, Forschungszentrum Jülich GmbH, Jülich, Germany
// 2018 - Dec - 23
//
// *********************************************************************
//
// Original implementation by A.Weickenmeier(wscatt.f)
// Modified:
// A.W.	07.08.90
// P.S.	01.14.91
// J.B.    23.11.2010
// J.B.    03.03.2011 (added parameters for H)
// J.B.    14.11.2014 (added function "dwfjbr" returning a dwf)
// J.B.    30.06.2017 (added parameters for L, z = 0, vacancy)
// J.B.    29.05.2018 (added f ^ 2(theta)integrals, wekof2)
// J.B.    23.12.2018 (translated to C code)
// 
// *********************************************************************
//
// ---------------------------------------------------------------------
//
// This program is free software : you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.If not, see < http://www.gnu.org/licenses/>. 
//
// ---------------------------------------------------------------------

#pragma once

#include "fcomplex.h"
#include "NatureConstants.h"

#ifndef __WEKOSCAT_H
#define __WEKOSCAT_H

#define _WEKOSCAT_NPRM			6 // number of table parameters
								  // B1 ... B6 as function of Z
#define _WEKOSCAT_SYMLEN		2 // length of the symbol strings
#define _WEKOSCAT_MINZ			0 // smallest Z supported
#define _WEKOSCAT_MAXZ			98 // largest Z supported
#define _WEKOSCAT_NPRMF			21 // number of extra parameters for rih2
#define _WEKOSCAT_NPRMEI		4 // number of parameters for ei
#define _WEKOSCAT_R8PI2A		(100. / (8.*_PI*_PI))  // reciprocal of 8*Pi^2 times 100
#define _WEKOSCAT_K0PA			(_TPI * _QEL / (_HPL*_C0) * 1.E-7) // prefactor for wavenumber [A^-1 * kV^-1]



class CWEKOScat {
public:

	// Constructor
	CWEKOScat();

	// Destructor
	~CWEKOScat();

	// Data Member
protected:
	double m_a[2]; // current a parameters (backup)
	double m_b[_WEKOSCAT_NPRM]; // current b parameters (backup)
	double m_dv[_WEKOSCAT_MAXZ + 1]; // V parameters
	double m_dprm[_WEKOSCAT_MAXZ + 1][_WEKOSCAT_NPRM]; // B parameters
	double m_dprmf[_WEKOSCAT_NPRMF]; // extra parameters for rih2
	double m_dprmeia[_WEKOSCAT_NPRMEI]; // a parameters of ei
	double m_dprmeib[_WEKOSCAT_NPRMEI]; // b parameters of ei
	char m_ssym[_WEKOSCAT_MAXZ + 1][_WEKOSCAT_SYMLEN + 1]; // Atom symbols

	// Member Functions
protected:

	// retrieve A and B parameters for atomic number Z
	// returns 0 if successful, otherwise an error code
	// - z : atmomic number
	// - a : double[2] : receives computed A parameters
	// - b : double[6] : receives copies of B parameters
	int getweko(int z, double* a, double* b);

	// returns elastic atomic form factor [A^2] from backup parameters
	// - s : scattering angle [1/A]
	double atffr1(double s);

	// calculates the exponential integral (tested for -60 < x < 60)
	double ei(double x);
	
	// calculates exp(-x1) * ( ei(x2)-ei(x3) )
	double rih1(double x1, double x2, double x3);
	
	// calculates x*exp(-x)*ei(x) for large x
	double rih2(double x);

	// form factor integral #1
	double ri1(double bi, double bj, double g);

	// form factor integral #2
	double ri2(double bi, double bj, double g, double u);

	// returns elastic atomic form factor [A^2]
	// - a : double[2] : A parameters
	// - b : double[6] : B parameters
	// - s : scattering angle [1/A]
	double weko(double* a, double* b, double s);

	// returns absorptive form factor [A]
	// - g : scattering vector [1/A]
	// - ul : mean square displacement of vibration [A^2]
	// - a : double[2] : A parameters
	// - b : double[6] : B parameters
	double wekoimag(double g, double ul, double* a, double *b);

public:

	// retrieve atom symbol for atomic number Z
	// returns 0 if successful, otherwise an error code
	// - z : atmomic number
	// - sym : receives atomic symbol
	int getwekosym(int z, char* sym);

	// calculate the Debye-Waller factor
	// - g : scattering vector [1/nm] ( = 2*s )
	// - dw : Debye-Waller parameter [nm^2] ( = 8*Pi^2 * <u^2> )
	// - dwflg : Flag for using the Debye Waller factor
	double getdwf(double g, double dw, bool dwflg);

	// calculate the elastic scattering factor [nm] with relativistic correction
	// - g : scattering vector [1/nm] ( = 2*s )
	// - dw : Debye-Waller parameter [nm^2] ( = 8*Pi^2 * <u^2> )
	// - z : atmomic number
	// - akv : accelerating voltage [kV]
	// - dwflg : Flag for using the Debye Waller factor
	double wekoscar(double g, double dw, int z, double akv, bool dwflg);
	
	// calculate the scattering factor [nm] with relativistic correction with absorptive component
	// - g : scattering vector [1/nm] ( = 2*s )
	// - dw : Debye-Waller parameter [nm^2] ( = 8*Pi^2 * <u^2> )
	// - z : atmomic number
	// - akv : accelerating voltage [kV]
	// - dwflg : Flag for using the Debye Waller factor
	// - absflg : Flag for calculating the absorptive form factor
	dcmplx wekosca(double g, double dw, int z, double akv, bool dwflg, bool absflg);
};

#endif
