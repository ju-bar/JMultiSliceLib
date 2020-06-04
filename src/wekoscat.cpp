// *********************************************************************
//
// file "wekoscat.cpp"
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

//#include "stdafx.h"
#include "wekoscat.h"
#include <stdlib.h>
#include <malloc.h>
#include <memory>
#include "integration.h"


CWEKOScat::CWEKOScat() :
	m_ssym{  "L",  "H", "He", "Li", "Be",  "B",  "C",  "N",  "O",  "F", "Ne",
			"Na", "Mg", "Al", "Si",  "P",  "S", "Cl", "Ar",	 "K", "Ca", "Sc",
			"Ti",  "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
			"As", "Se",  "B", "Kr",	"Rb", "Sr",  "Y", "Zr", "Nb", "Mo", "Tc",
			"Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",  "I", "Xe",
			"Cs", "Ba", "La", "Ce",	"Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
			"Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",  "W", "Re", "Os",
			"Ir", "Pt", "Au", "Hg", "Tl", "Pb",	"Bi", "Po", "At", "Rn",	"Fr",
			"Ra", "Ac", "Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf" },
	m_dv{	 0.0,  0.5,  0.5,  0.5,  0.3,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5, 
			 0.5,  0.5,  0.4,  0.5,  0.5,  0.5,  0.5,  0.5,  0.2,  0.3,  0.5,
			 0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,
	         0.5,  0.5,  0.5,  0.5,  0.2,  0.3,  0.5,  0.5,  0.5,  0.5,  0.5,
			 0.4,  0.5,  0.5,  0.5,  0.3,  0.4,  0.6,  0.6,  0.6,  0.4,  0.4,
			 0.1,  0.1,  0.3,  0.3,  0.2,  0.2,  0.2,  0.2,  0.1,  0.2,  0.1,
			 0.2,  0.1,  0.2,  0.1,  0.1,  0.1,  0.1,  0.4,  0.2,  0.5,  0.4,
			 0.5,  0.5,  0.4,  0.4,  0.4,  0.3,  0.4,  0.4,  0.4,  0.4,  0.1,
			 0.2,  0.2,  0.3,  0.2,  0.2,  0.2,  0.2,  0.2,  0.3,  0.2,  0.3 },
	m_dprm{ { 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000}, // 0
			{48.75740,  4.96588, 18.24440, 18.24440, 18.24440, 18.24440}, // 1
			{ 2.54216,  8.74302, 12.69098,  0.43711,  5.29446, 28.25045},
			{ 0.68454,  3.06497,  6.23974,126.17816,131.20160,131.76538},
			{ 0.53996,  3.38752, 55.62340, 50.78098, 67.00502, 96.36635},
			{ 0.33138,  2.97485, 34.01118, 35.98365, 36.68364, 60.80991},
			{ 0.29458,  3.93381, 24.97836, 25.27916, 25.46696, 46.70328},
			{ 0.23925,  4.93515, 18.11895, 15.69698, 15.81922, 40.24150},
			{ 6.37582,  8.03744, 27.20649,  0.11157,  0.38686, 10.89944},
			{ 0.21800,  6.76987,  7.05056,  6.67484, 12.38148, 28.08398},
			{ 0.20055,  5.49814,  6.28052,  7.19211,  7.54763, 23.26388}, // 10
			{ 0.21902,  5.30022,  5.31938,  5.28281,  5.28546,128.18391},
			{ 1.97633,  2.80902, 16.39184,  0.05494,  2.06121,121.70512},
			{ 2.29692,  2.35822, 24.98576,  0.07462,  0.55953,128.50104},
			{ 1.73656,  3.04329, 30.57191,  0.05070,  0.99181, 86.18340},
			{ 0.17949,  2.63250,  2.67559, 34.57098, 36.77888, 54.06180},
			{ 1.00609,  4.90414, 31.34909,  0.03699,  0.98700, 44.94354},
			{ 0.18464,  1.47963,  5.20989, 24.79470, 32.06184, 39.09933},
			{ 0.20060,  6.53262, 22.72092,  1.20022,  1.27398, 36.25907},
			{ 0.44424,  3.36735, 19.63031,  0.01824, 23.51332,212.86819},
			{ 0.18274,  2.06638, 16.99062, 11.57795, 13.97594,186.10446}, // 20
			{ 0.14245,  1.46588, 15.46955,  4.24287,  9.80399,121.46864},
			{ 0.12782,  1.45591, 12.09738,  4.61747, 11.96791,105.00546},
			{ 0.13126,  1.39923,  8.00762,  7.98129, 13.41408, 95.30811},
			{ 0.12311,  2.38386,  9.92149,  1.64793, 11.00035, 68.45583},
			{ 0.48173,  3.78306,  8.47337,  0.04690,  8.74544, 77.44405},
			{ 0.44704,  6.89364,  6.90335,  0.05691,  3.02647, 70.86599},
			{ 0.10705,  3.63573,  7.55825,  1.27986,  5.14045, 67.16051},
			{ 0.11069,  1.61889,  6.00325,  5.97496,  6.06049, 59.41419},
			{ 0.11293,  1.89077,  5.08503,  5.07335,  5.09928, 46.38955},
			{ 0.10209,  1.73365,  4.78298,  4.80706,  5.64485, 51.21828}, // 30
			{ 0.10642,  1.53735,  5.13798,  4.74298,  4.99974, 61.42872},
			{ 0.09583,  1.67715,  4.70275,  2.91198,  7.87009, 64.93623},
			{ 0.09428,  2.21409,  3.95060,  1.52064, 15.81446, 52.41380},
			{ 0.09252,  1.60168,  3.04917,  3.18476, 18.93890, 47.62742},
			{ 0.09246,  1.77298,  3.48134,  1.88354, 22.68630, 40.69434},
			{ 0.49321,  2.08254, 11.41282,  0.03333,  2.09673, 42.38068},
			{ 0.15796,  1.71505,  9.39164,  1.67464, 23.58663,152.53635},
			{ 0.36052,  2.12757, 12.45815,  0.01526,  2.10824,133.17088},
			{ 0.09003,  1.41396,  2.05348, 10.25766, 10.74831, 90.63555},
			{ 0.10094,  1.15419,  2.34669, 10.58145, 10.94962, 82.82259}, // 40
			{ 0.09243,  1.16977,  5.93969,  1.30554, 13.43475, 66.37486},
			{ 0.43543,  1.24830,  7.45369,  0.03543,  9.91366, 61.72203},
			{ 0.45943,  1.18155,  8.31728,  0.03226,  8.32296, 64.97874},
			{ 0.08603,  1.39552, 11.69728,  1.39552,  3.45200, 55.55519},
			{ 0.09214,  1.11341,  7.65767,  1.12566,  8.32517, 48.38017},
			{ 0.09005,  1.12460,  9.69801,  1.08539,  5.70912, 33.48585},
			{ 0.08938,  3.19060,  9.10000,  0.80898,  0.81439, 41.34453},
			{ 0.28851,  1.61312,  8.99691,  0.01711,  9.46666, 58.13256},
			{ 0.08948,  1.23258,  8.23129,  1.22390,  7.06201, 59.69622},
			{ 0.07124,  0.85532,  6.40081,  1.33637,  6.38240, 50.92361}, // 50
			{ 0.35749,  1.32481,  6.51696,  0.03550,  6.51913, 50.80984},
			{ 0.50089,  3.95301,  7.62830,  0.03005,  0.50737, 49.62628},
			{ 0.08429,  1.12959,  8.86209,  1.12981,  9.13243, 56.01965},
			{ 0.27796,  1.62147, 11.45200,  0.02032,  3.27497, 51.44078},
			{ 0.12045,  1.53654,  9.81569, 41.21656, 42.62216,224.34816},
			{ 0.12230,  1.44909,  9.50159, 49.40860, 74.94942,217.04485},
			{ 0.08930,  1.26225,  8.09703,  1.20293, 17.65554,116.61481},
			{ 0.08504,  1.28286, 11.22123,  1.32741,  4.61040,112.19678},
			{ 0.09805,  1.52628,  8.58953,  1.23893, 22.49126,140.02856},
			{ 0.09413,  1.26616,  5.98844, 17.78775, 18.14397,132.59305}, // 60
			{ 0.09447,  1.25111,  5.91205, 16.28675, 16.73089,127.90916},
			{ 0.09061,  1.59281, 10.64077,  1.78861,  2.22148,124.56328},
			{ 0.10485,  1.54396,  8.65223,  7.09290, 53.36537,183.69014},
			{ 0.09338,  1.38681,  7.35883,  1.55122, 20.81916,111.03201},
			{ 0.10190,  1.52368,  7.16923, 20.86269, 49.29465,166.09206},
			{ 0.08402,  1.40890,  7.14042,  1.34848, 11.42203,108.01204},
			{ 0.09441,  1.61807,  6.27142, 40.34946, 42.82722,130.59616},
			{ 0.08211,  1.25106,  4.81241, 10.84493, 10.90164,100.07855},
			{ 0.09662,  1.60236,  5.67480, 30.59014, 31.12732,138.69682},
			{ 0.09493,  1.60220,  5.43916, 28.31076, 29.27660,138.08665}, // 70
			{ 0.09658,  1.56751,  5.32170, 34.18217, 35.25187,121.42893},
			{ 0.09294,  1.55499,  5.25121, 37.51883, 38.88302,105.16978},
			{ 0.06298,  0.81950,  2.89124,  5.54290,  5.98101, 54.42459},
			{ 0.07902,  1.37096,  8.23364,  1.38300,  1.39219, 77.11813},
			{ 0.05266,  0.90718,  4.43830,  0.94590,  4.37477, 43.97909},
			{ 0.22700,  1.56975,  6.34451,  0.01564,  1.61769, 46.15815},
			{ 0.05055,  0.86775,  5.09325,  0.88123,  3.56919, 39.77390},
			{ 0.05253,  0.83773,  3.95899,  0.81515,  6.44217, 34.21146},
			{ 0.54927,  1.72752,  6.71952,  0.02637,  0.07253, 35.45745},
			{ 0.21941,  1.41611,  6.68241,  0.01472,  1.57578, 37.15826}, // 80
			{ 0.22459,  1.12822,  4.30289,  0.01485,  7.15607, 43.08737},
			{ 0.06432,  1.19406,  7.39342,  1.14160,  1.28905, 51.13401},
			{ 0.05380,  0.86719,  1.87540,  7.64796,  7.86794, 45.63897},
			{ 0.50112,  1.63784,  6.78551,  0.02187,  0.08602, 46.72951},
			{ 0.22321,  1.10827,  3.59116,  0.01011, 11.63732, 45.06839},
			{ 0.21152,  1.14015,  3.41473,  0.01188, 13.41211, 43.11389},
			{ 0.09435,  1.02649,  6.25480, 32.51444, 36.29119,149.11722},
			{ 0.07300,  1.01825,  5.89629,  1.03089, 20.37389,115.34722},
			{ 0.07515,  0.94941,  3.72527, 17.58346, 19.75388,109.12856},
			{ 0.06385,  0.90194,  4.65715,  0.90253, 15.70771, 83.69695}, // 90
			{ 0.07557,  0.84920,  4.00991, 16.95003, 17.78767,100.20415},
			{ 0.07142,  1.14907,  9.21231,  0.95923,  1.20275,104.32746},
			{ 0.06918,  0.98102,  5.95437,  0.99086, 22.06437, 90.98156},
			{ 0.07136,  0.95772,  6.13183,  0.97438, 15.67499, 89.86625},
			{ 0.07301,  0.93267,  6.34836,  0.91032, 13.26179, 86.85986},
			{ 0.05778,  0.72273,  3.01146,  9.21882,  9.53410, 65.86810},
			{ 0.07088,  0.77587,  6.14295,  1.79036, 15.12379, 83.56983},
			{ 0.06164,  0.81363,  6.56165,  0.83805,  4.18914, 61.41408} // 98
			 },
	m_dprmf{  1.000000, 1.005051, 1.010206, 1.015472, 1.020852,
			  1.026355, 1.031985, 1.037751, 1.043662, 1.049726,
			  1.055956, 1.062364, 1.068965, 1.075780, 1.082830,
			  1.090140, 1.097737, 1.105647, 1.113894, 1.122497, 1.131470 },
	m_dprmeia{8.57332, 18.05901,  8.63476,  0.26777 },
	m_dprmeib{9.57332, 25.63295, 21.09965,  3.95849 }
{
	m_g = 0.;
	m_k0 = 0.;
	m_dw = 0.;
}

CWEKOScat::~CWEKOScat()
{
}

int CWEKOScat::getweko(int z, double *a, double* b)
{
	if (z<_WEKOSCAT_MINZ || z>_WEKOSCAT_MAXZ) { return 1; }
	if (NULL == a || NULL == b) { return 2; }
	a[0] = _WKS_CFFA * (double)z / (3.0 * (1.0 + m_dv[z]));
	a[1] = m_dv[z] * a[0];
	if (b != memcpy(b, m_dprm[z], sizeof(double)*_WEKOSCAT_NPRM)) {
		return 3;
	}
	return 0;
}

int CWEKOScat::getwekosym(int z, char* sym)
{
	if (z<_WEKOSCAT_MINZ || z>_WEKOSCAT_MAXZ) { return 1; }
	if (NULL == sym) { return 2; }
	sym = m_ssym[z];
	return 0;
}

double CWEKOScat::atffr1(double s)
{
	return weko(m_a, m_b, s);
}


double CWEKOScat::weko(double* a, double* b, double s)
{
	double ff = 0.0;
	double s2 = 0.0;
	double argu = 0.0;
	int i = 0, j = 0;
	if (s > 0.0) {
		s2 = 1. / (s*s);
	}
	for (i = 0; i < _WEKOSCAT_NPRM; i++) {
		j = (int)((double)i / 3.);
		argu = b[i] * s * s;
		if (argu < 0.1) {
			ff += (a[j] * b[i] * (1.0 - 0.5 * argu));
		}
		else if (argu > 20.0) {
			ff += (a[j] * s2);
		}
		else {
			ff += (a[j] * (1.0 - exp(-argu)) * s2);
		}
	}
	return ff;
}


double CWEKOScat::ei(double x)
{
	double r = 0.0;
	if (x > 60.) { // error >>> 'ei fuer x= ', x, ' nicht getestet <<<'
		return r;
	}
	if (x < -60.) {
		return r;
	}
	if (x < -1.) { // abramowitz(5.1.56)
		double xp = abs(x);
		r = -(m_dprmeia[3] + xp * (m_dprmeia[2] + xp * (m_dprmeia[1] + xp * (m_dprmeia[0] + xp)))) /
			(m_dprmeib[3] + xp * (m_dprmeib[2] + xp * (m_dprmeib[1] + xp * (m_dprmeib[0] + xp)))) * exp(-xp) / xp;
		return r;
	}
	r = .577216 + log(abs(x));
	int i = 0;
	double si = x;
	double summ = si;
	double ri = 0., ri1 = 0.;
	while(abs(si / x) > 1.e-6) {
		i++;
		ri = (double)i;
		ri1 = ri + 1.;
		si *= (x * ri / (ri1*ri1));
		summ += si;
	}
	r += summ;
	return r;
}


double CWEKOScat::rih1(double x1, double x2, double x3)
{
	double r = 0.0;
	if (x2 <= 20.  && x3 <= 20.) {
		r = exp(-x1) * (ei(x2) - ei(x3));
		return r;
	}
	if (x2 > 20.) {
		r = exp(x2 - x1)*rih2(x2) / x2;
	}
	else {
		r = exp(-x1)*ei(x2);
	}
	if (x3 > 20.) {
		r -= exp(x3 - x1)*rih2(x3) / x3;
	}
	else {
		r -= exp(-x1)*ei(x3);
	}
	return r;
}


double CWEKOScat::rih2(double x)
{
	double r = 0.0;
	double x1 = 1. / x;
	int i = (int)(200.*x1);
	int	i1 = i + 1;
	r = m_dprmf[i] + 200.*(m_dprmf[i1] - m_dprmf[i]) * (x1 - 0.005*(double)i);
	return r;
}


double CWEKOScat::ri1(double bi, double bj, double g)
{
	double r = 0.0;
	double c = 0.5772157;
	if (g == 0.0) {
		r = bi * log((bi + bj) / bi) + bj * log((bi + bj) / bj);
		r *= _WKS_PI;
		return r;
	}
	double g2 = g * g;
	double big2 = bi * g2;
	double bjg2 = bj * g2;
	r = 2.*c + log(big2) + log(bjg2) - 2.*ei(-bi * bj*g2 / (bi + bj));
	double x1 = big2;
	double x2 = big2 * bi / (bi + bj);
	double x3 = big2;
	r += rih1(x1, x2, x3);
	x1 = bjg2;
	x2 = bjg2 * bj / (bi + bj);
	x3 = bjg2;
	r += rih1(x1, x2, x3);
	r *= (_WKS_PI / g2);
	return r;
}

double CWEKOScat::ri2(double bi, double bj, double g, double u)
{
	double r = 0.0;
	double u2 = u * u;
	if (g == 0.) {
		r = (bi + u2) * log((bi + bj + u2) / (bi + u2));
		r += bj * log((bi + bj + u2) / (bj + u2));
		r += u2 * log(u2 / (bj + u2));
		r *= _WKS_PI;
		return r;
	}
	double g2 = g * g;
	double biuh = bi + 0.5 * u2;
	double bjuh = bj + 0.5 * u2;
	double	biu = bi + u2;
	double	bju = bj + u2;
	r = ei(-0.5*u2*g2*biuh / biu) + ei(-0.5*u2*g2*bjuh / bju);
	r -= (ei(-biuh * bjuh*g2 / (biuh + bjuh)) + ei(-.25*u2*g2));
	r *= 2.0;
	double x1 = 0.5*u2*g2;
	double x2 = 0.25*u2*g2;
	double x3 = 0.25*u2*u2*g2 / biu;
	r += rih1(x1, x2, x3);
	x1 = 0.5*u2*g2;
	x2 = 0.25*u2*g2;
	x3 = 0.25*u2*u2*g2 / bju;
	r += rih1(x1, x2, x3);
	x1 = biuh * g2;
	x2 = biuh * biuh*g2 / (biuh + bjuh);
	x3 = biuh * biuh*g2 / biu;
	r += rih1(x1, x2, x3);
	x1 = bjuh * g2;
	x2 = bjuh * bjuh*g2 / (biuh + bjuh);
	x3 = bjuh * bjuh*g2 / bju;
	r += rih1(x1, x2, x3);
	r *= (_WKS_PI / g2);
	return r;
}


double CWEKOScat::wekoimag(double g, double ul, double* a, double *b)
{
	int i = 0, j = 0, ii = 0, jj = 0;
	double fi = 0.0;
	double fp2 = _WKS_FPI * _WKS_FPI;
	double a1[2], b1[_WEKOSCAT_NPRM];
	double u2 = ul * ul;
	for (i = 0; i < 2; i++) { a1[i] = a[i] * fp2; }
	for (i = 0; i < (int)_WEKOSCAT_NPRM; i++) { b1[i] = b[i] / fp2; }
	double g2 = g * g;
	double dewa = exp(-0.5*ul*g2);
	for (j = 0; j < (int)_WEKOSCAT_NPRM; j++) {
		jj = (int)((double)j / 3.);
		for (i = 0; i < (int)_WEKOSCAT_NPRM; i++) {
			ii = (int)((double)i / 3.);
			fi += ( a1[jj] * a1[ii]
					*( dewa * ri1(b1[i],b1[j],g) 
						- ri2(b1[i],b1[j],g,ul) ) );
		}
	}
	return fi;
}

double CWEKOScat::wekoscar1(double s)
{
	return weko(m_a, m_b, s);
}

double CWEKOScat::wekomug(double theta, double phi)
{
	double vap = 0.;
	double ct = std::cos(theta);
	double st = std::sin(theta);
	double cp = std::cos(phi);
	// q = scattering vector length of Q in the ewald sphere
	// Q = (qx, qy, qz) = K' - K
	// qx = k * sin(theta) * cos(phi)
	// qy = k * sin(theta) * sin(phi)
	// qz = k * (cos(theta) - 1)
	// k = | K | = wekok
	double k = m_k0;
	double twok = k + k;
	// q = Sqrt[2 k ^ 2 (1 - Cos[q])]
	double q = k * std::sqrt(2. - 2. * ct);
	// g = length of some reciprocal space vector G = (gx, gy, gz)
	// with gx = g, gy = 0, gz = 0 (obda)
	double g = m_g;
	// qg = length of the difference vector Q - G
	double qg = std::sqrt(g * g + twok * k - twok * (k * ct + g * cp * st));
	double biso = m_dw; // biso value
	double sg = 0.5 * g;
	double sq = 0.5 * q;
	double sqg = 0.5 * qg;
	double fq = wekoscar1(sq);
	double fqg = wekoscar1(sqg);
	// Debye-Waller factors values
	double ag = std::exp(-biso * sg * sg);
	double aq = std::exp(-biso * sq * sq);
	double aqg = std::exp(-biso * sqg * sqg);
	vap = fq * fqg * (ag - aq * aqg) * st;
	return vap;
}

double CWEKOScat::wekomugap(double theta, double phi)
{
	double vap = 0.;
	double ct = std::cos(theta);
	double st = std::sin(theta);
	double cp = std::cos(phi);
	// q = scattering vector length of Q in the ewald sphere
	// Q = (qx, qy, qz) = K' - K
	// qx = k * sin(theta) * cos(phi)
	// qy = k * sin(theta) * sin(phi)
	// qz = k * (cos(theta) - 1)
	// k = | K | = wekok
	double k = m_k0;
	double twok = k + k;
	// q = Sqrt[2 k ^ 2 (1 - Cos[q])]
	double q = k * std::sqrt(2. - 2. * ct);
	// g = length of some reciprocal space vector G = (gx, gy, gz)
	// with gx = g, gy = 0, gz = 0 (obda)
	double g = m_g;
	// qg = length of the difference vector Q - G
	double qg = std::sqrt(g * g + twok * k - twok * (k * ct + g * cp * st));
	double sgmax = m_dw * 0.5; // aperture radius on the scale of s = g / 2
	double sg = 0.5 * g;
	double sq = 0.5 * q;
	double sqg = 0.5 * qg;
	double fq = wekoscar1(sq);
	double fqg = wekoscar1(sqg);
	// aperture values
	double ag = 1.;
	double aq = 1.;
	double aqg = 1.;
	if (sg >= sgmax) ag = 0.;
	if (sq >= sgmax) aq = 0.;
	if (sqg >= sgmax) aqg = 0.;
	vap = fq * fqg * (ag - aq * aqg) * st;
	return vap;
}

double CWEKOScat::wekoabs(double g, double dw, double* a, double* b, double k0)
{
	double aff = 0.;
	double si0 = 0., si1 = 0., pi0 = 0., pi1 = 0.;
	// copy parameters to module members to be used by the integrator
	memcpy(m_a, a, sizeof(double) * 2);
	memcpy(m_b, b, sizeof(double) * _WEKOSCAT_NPRM);
	m_g = g;
	m_dw = dw;
	m_k0 = k0;
	// set the integration range
	si0 = 0.; si1 = _WKS_PI; // theta (0 ... Pi)
	pi0 = 0.; pi1 = _WKS_PI; // phi (0 ... Pi)
	// call the grid integrator on 128 x 64 pixels
	// !  theta with quadratic sampling, phi with linear sampling
	// !  phi on the half side only (symmteric) therefore:    * 2
	aff = dsgrid2d<CWEKOScat>(&CWEKOScat::wekomug, si0, si1, pi0, pi1, 2.0, 1.0, 128, 64, *this) * 2.0;
	return aff;
}

double CWEKOScat::wekoabsap(double g, double ap, double* a, double* b, double k0)
{
	double aff = 0.;
	double si0 = 0., si1 = 0., pi0 = 0., pi1 = 0.;
	if (g < ap) {
		// copy parameters to module members to be used by the integrator
		memcpy(m_a, a, sizeof(double) * 2);
		memcpy(m_b, b, sizeof(double) * _WEKOSCAT_NPRM);
		m_g = g;
		m_dw = ap;
		m_k0 = k0;
		// set the integration range
		si0 = 0.; si1 = _WKS_PI; // theta (0 ... Pi)
		pi0 = 0.; pi1 = _WKS_PI; // phi (0 ... Pi)
		// call the grid integrator on 128 x 64 pixels
		// !  theta with quadratic sampling, phi with linear sampling
		// !  phi on the half side only (symmteric) therefore:    * 2
		aff = dsgrid2d<CWEKOScat>(&CWEKOScat::wekomugap, si0, si1, pi0, pi1, 2.0, 1.0, 128, 64, *this) * 2.0;
	}
	return aff;
}

double CWEKOScat::relcorr(double ekv)
{
	return (_WKS_E0KV + ekv) / _WKS_E0KV;
}

double CWEKOScat::getdwf(double g, double dw, bool dwflg)
{
	if (!dwflg) return 1.0;
	return exp(-0.25 * dw * g * g);
}


double CWEKOScat::wekoscar(double g, double dw, int z, double akv, bool dwflg)
{
	if (z<_WEKOSCAT_MINZ || z>_WEKOSCAT_MAXZ) { return 0.0; }
	double ff = 0.0, fa = 0.0;
	double rc = relcorr(akv);
	double dwf = getdwf(g, dw, dwflg);
	double sa = 0.050 * g;  // s = g / 2 [A]
	double a[2], b[_WEKOSCAT_NPRM];
	if (0 == getweko(z, a, b)) {
		fa = weko(a, b, sa); // ff in [A]
		ff = 0.1 * fa * rc * _WKS_FPI * dwf; // to [nm] with rel. corr. DWF and factor 4*Pi (used by EMS)
	}
	return ff;
}


dcmplx CWEKOScat::wekosca(double g, double dw, int z, double akv, bool dwflg, bool absflg)
{
	dcmplx dff = dcmplx(0, 0);
	if (z<_WEKOSCAT_MINZ || z>_WEKOSCAT_MAXZ) { return dff; }
	double ff = 0.0, fa = 0.0;
	double rc = relcorr(akv);
	double dwf = getdwf(g, dw, dwflg);
	double sa = 0.050 * g;  // s = g / 2 [A]
	double a[2], b[_WEKOSCAT_NPRM];
	if (0 == getweko(z, a, b)) {
		fa = weko(a, b, sa); // ff in [A]
		ff = 0.1 * fa * rc * _WKS_FPI * dwf; // to [nm] with rel. corr. DWF and factor 4*Pi (used by EMS)
		dff.real(ff);
	}
	if (absflg && dwflg && dw > _WKS_EPS && ff > _WKS_EPS) { // calculate absorptive form factor for the elastic channel due to TDS
		double fi = 0.0;
		// wave number in [A^-1] 0.506774 [A^-1 * kV^-1]
		double k0 = _WKS_K0PA * sqrt((2. * _WKS_E0KV + akv) * akv); // 2 * Pi*k !andere k - Notation hier : Aufpassen!
		double ga = 0.1 * g * _WKS_TPI;
		double ua = sqrt(dw * _WKS_R8PI2A);
		fi = 0.1 * rc * rc * wekoimag(ga, ua, a, b) / k0;
		dff.imag(fi);
	}
	if (absflg && !dwflg && g < dw) { // calculate absorptive form factor due to a band-width limiting aperture
		double fi = 0.0;
		// wave number in [A^-1] 0.506774 [A^-1 * kV^-1]
		double k0 = _WKS_K0PA * sqrt((2. * _WKS_E0KV + akv) * akv); // 2 * Pi*k !andere k - Notation hier : Aufpassen!
		fi = 0.1 * rc * rc * wekoabsap(0.1 * g, 0.1 * dw, a, b, k0 / _WKS_TPI) * k0;
		dff.imag(fi);
	}
	return dff;
}