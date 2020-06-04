// file: 'prime_numbers.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementation of routine on prime numbers.
//
/* -----------------------------------------------------------------------

	This file is part of JMSBench1.

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

#include "prime_numbers.h"
#include <math.h>
#include <algorithm>

size_t largest_prime_factor(size_t n)
{
	if (n == 0) return 0;
	if (n == 1) return 1;
	if (n == 2) return 2;
	if (n == 3) return 3;
	if (n == 4) return 2;
	if (n == 5) return 5;
	size_t ntest = (size_t)sqrt((double)n);
	size_t i = ntest, j = 1, n1 = 1, n2 = 1;
	for (i = ntest; i > 1; i--) {
		if (0 == n % i) {
			j = (size_t)n / i;
			n1 = largest_prime_factor(i);
			n2 = largest_prime_factor(j);
			j = std::max(n1, n2);
			return j;
		}
	}
	return n;
}

size_t next_low_prime(size_t n, size_t pfmax)
{
	size_t n2 = (size_t)pow(2.0, 1.0 + floor(log((double)n) / log(2.0))); // next power of two number, limit search to there
	size_t nn = n2;
	size_t n1 = n;
	size_t nlp = 0;
	if (n > 6 && n < n2) {
		for (size_t i = n1; i < n2; i++) {
			nlp = largest_prime_factor(i);
			if (nlp <= pfmax) {
				nn = i;
				break;
			}
		}
	}
	return nn;
}