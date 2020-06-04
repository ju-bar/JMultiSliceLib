// file: 'integration.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the declaration of integration routines.
//
/* -----------------------------------------------------------------------

	This file is part of JMultiSlice.

	JMultiSlice is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	JMultiSlice is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with JMultiSlice.  If not, see <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------- */

#pragma once

#include <cmath>
#include <algorithm>

constexpr double _INT_GRID_PMIN = 0.1; // minimum grid sampling power
constexpr double _INT_GRID_PMAX = 10.; // maximum grid sampling power
constexpr unsigned int _INT_GRID_NMIN = 8; // minimum grid size
constexpr unsigned int _INT_GRID_NMAX = 2048; // maximum grid size

// numerical integration of a function foo on a 2d grid of (n1, n2) samples
// over the range (a1, b1)_n1 and (a2, b2)_n2
// using a grid sampling with powers p1 on dimension n1 and p2 on dimension n2
// can be called with func being a member of a template class C
template<class C>
double dsgrid2d(double (C::*func)(double, double), double a1, double b1, double a2, double b2, double p1, double p2, unsigned int n1, unsigned int n2, C& c)
{
	double dresult = 0.;
	// (this may not be thread safe)
	static double* _x1 = NULL; // grid x values
	static double* _x2 = NULL; // grid y values
	static double* _s1 = NULL; // intermediate sums
	static double* _y = NULL; // data
	static unsigned int _n1 = 0; // pre-allocated dimensions of the static grid on dimension 1
	static unsigned int _n2 = 0; // pre-allocated dimensions of the static grid on dimension 2
	//
	unsigned int ni1 = std::min(std::max(n1, _INT_GRID_NMIN), _INT_GRID_NMAX);
	unsigned int ni2 = std::min(std::max(n2, _INT_GRID_NMIN), _INT_GRID_NMAX);
	double pi1 = std::min(std::max(p1, _INT_GRID_PMIN), _INT_GRID_PMAX);
	double pi2 = std::min(std::max(p2, _INT_GRID_PMIN), _INT_GRID_PMAX);
	unsigned int i = 0, j = 0, idx = 0, ix = 0;
	double cx1 = 0., cx2 = 0., cy = 0.;
	//
	// check static grid data
	if (ni1 > _n1 || ni2 > _n2) { // increased size, update allocations
		if (_y != NULL) { free(_y); _y = NULL; }
		if (_x1 != NULL) { free(_x1); _x1 = NULL; }
		if (_x2 != NULL) { free(_x2); _x2 = NULL; }
		if (_s1 != NULL) { free(_s1); _s1 = NULL; }
		if (NULL == _y) _y = (double*)malloc(sizeof(double) * ni1 * ni2);
		if (NULL == _x1) _x1 = (double*)malloc(sizeof(double) * ni1);
		if (NULL == _s1) _s1 = (double*)malloc(sizeof(double) * ni1);
		if (NULL == _x2) _x2 = (double*)malloc(sizeof(double) * ni2);
		_n1 = ni1; _n2 = ni2;
	}
	// fill data
	cx1 = 1.0 / (double)(ni1 - 1);
	for (i = 0; i < ni1; i++) {
		_x1[i] = a1 + (b1 - a1) * std::pow(cx1 * (double)i, pi1);
		_s1[i] = 0.;
	}
	cx2 = 1.0 / (double)(ni2 - 1);
	for (j = 0; j < ni2; j++) {
		_x2[j] = a2 + (b2 - a2) * std::pow(cx2 * (double)j, pi2);
		ix = j * ni1;
		for (i = 0; i < ni1; i++) {
			idx = i + ix;
			_y[idx] = (c.*func)(_x1[i], _x2[j]);
		}
	}
	// integrate with trapezoidal summation
	for (j = 0; j < ni2; j++) { // - inner sums
		ix = j * ni1;
		for (i = 0; i < ni1 - 1; i++) {
			idx = i + ix;
			_s1[j] += (_y[idx + 1] + _y[idx]) * (_x1[i + 1] - _x1[i]);
		}
		_s1[j] *= 0.5; // half sum (trapez)
	}
	for (j = 0; j < ni2 - 1; j++) { // - outer sums
		cy += (_s1[j + 1] + _s1[j]) * (_x2[j + 1] - _x2[j]);
	}
	dresult = 0.5 * cy; // half sum (trapez)
	return dresult;
}

// numerical integration of a function foo on a 2d grid of (n1, n2) samples
// over the range (a1, b1)_n1 and (a2, b2)_n2
// using a grid sampling with powers p1 on dimension n1 and p2 on dimension n2
double dsgrid2d(double (*func)(double, double), double a1, double b1, double a2, double b2, double p1, double p2, unsigned int n1, unsigned int n2)
{
	double dresult = 0.;
	// (this may not be thread safe)
	static double* _x1 = NULL; // grid x values
	static double* _x2 = NULL; // grid y values
	static double* _s1 = NULL; // intermediate sums
	static double* _y = NULL; // data
	static unsigned int _n1 = 0; // pre-allocated dimensions of the static grid on dimension 1
	static unsigned int _n2 = 0; // pre-allocated dimensions of the static grid on dimension 2
	//
	unsigned int ni1 = std::min(std::max(n1, _INT_GRID_NMIN), _INT_GRID_NMAX);
	unsigned int ni2 = std::min(std::max(n2, _INT_GRID_NMIN), _INT_GRID_NMAX);
	double pi1 = std::min(std::max(p1, _INT_GRID_PMIN), _INT_GRID_PMAX);
	double pi2 = std::min(std::max(p2, _INT_GRID_PMIN), _INT_GRID_PMAX);
	unsigned int i = 0, j = 0, idx = 0, ix = 0;
	double cx1 = 0., cx2 = 0., cy = 0.;
	//
	// check static grid data
	if (ni1 > _n1 || ni2 > _n2) { // increased size, update allocations
		if (_y != NULL) { free(_y); _y = NULL; }
		if (_x1 != NULL) { free(_x1); _x1 = NULL; }
		if (_x2 != NULL) { free(_x2); _x2 = NULL; }
		if (_s1 != NULL) { free(_s1); _s1 = NULL; }
		if (NULL == _y) _y = (double*)malloc(sizeof(double) * ni1 * ni2);
		if (NULL == _x1) _x1 = (double*)malloc(sizeof(double) * ni1);
		if (NULL == _s1) _s1 = (double*)malloc(sizeof(double) * ni1);
		if (NULL == _x2) _x2 = (double*)malloc(sizeof(double) * ni2);
		_n1 = ni1; _n2 = ni2;
	}
	// fill data
	cx1 = 1.0 / (double)(ni1 - 1);
	for (i = 0; i < ni1; i++) {
		_x1[i] = a1 + (b1 - a1) * std::pow(cx1 * (double)i, pi1);
		_s1[i] = 0.;
	}
	cx2 = 1.0 / (double)(ni2 - 1);
	for (j = 0; j < ni2; j++) {
		_x2[j] = a2 + (b2 - a2) * std::pow(cx2 * (double)j, pi2);
		ix = j * ni1;
		for (i = 0; i < ni1; i++) {
			idx = i + ix;
			_y[idx] = func(_x1[i], _x2[j]);
		}
	}
	// integrate with trapezoidal summation
	for (j = 0; j < ni2; j++) { // - inner sums
		ix = j * ni1;
		for (i = 0; i < ni1 - 1; i++) {
			idx = i + ix;
			_s1[j] += (_y[idx + 1] + _y[idx]) * (_x1[i + 1] - _x1[i]);
		}
		_s1[j] *= 0.5; // half sum (trapez)
	}
	for (j = 0; j < ni2 - 1; j++) { // - outer sums
		cy += (_s1[j + 1] + _s1[j]) * (_x2[j + 1] - _x2[j]);
	}
	dresult = 0.5 * cy; // half sum (trapez)
	return dresult;
}