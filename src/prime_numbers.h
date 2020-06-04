// file: 'prime_numbers.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the declaration of routines on prime numbers.
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

// returns the largest prime factor of n (recursive)
size_t largest_prime_factor(size_t n);

// returns the next number following n which has only prime factors smaller or equal to pfmax
size_t next_low_prime(size_t n, size_t pfmax);