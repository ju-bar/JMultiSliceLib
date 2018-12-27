//
// C++ header file: fcomplex.h
// declaration of standard float complex type
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
#pragma once
#include <iostream>
#include <math.h>
#include <complex>
//
// declare basic complex data type
#ifndef FCOMPLEX_H_
#define FCOMPLEX_H_
typedef std::complex<float> fcmplx;
#define FJ fcmplx(0.F,1.F)
#endif /* FCOMPLEX_H_ */
//
#ifndef DCOMPLEX_H_
#define DCOMPLEX_H_
typedef std::complex<double> dcmplx;
#define DJ dcmplx(0.,1.)
#endif /* DCOMPLEX_H_ */

//

