// file: 'rng.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the declaration of the class CRng to be used as
// pseudo random number generator based on the std <random> library.
//
/* -----------------------------------------------------------------------

	This file is distributed as part of the JMultiSlice library.

	JMultiSlice is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Merlinio is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with Merlinio.  If not, see <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------- */

#pragma once
#include <random>
#define RAND_MAXROLL	0x7FFFFFFF
class CRng
{
private:
	std::random_device rd;
	typedef std::mt19937 urng_t;
	urng_t urng;
	std::uniform_int_distribution<int> uni_dist;
	std::normal_distribution<double> norm_dist;
	int m_seed;
public:
	CRng();
	CRng(int seed);
	~CRng();

	// re-seed the random number generator
	void seed(int seed = 0);

	// returns the last seed number used
	int get_seed(void);
	
	// get integer pseudo random number from 0 to RAND_MAXROLL
	int rand();

	// get uniform pseudo random number in closed range [rmin, rmax]
	double unirand(double rmin=0.0, double rmax=1.0);

	// get uniform pseudo random number in open range [rmin, rmax)
	double unirand_oub(double rmin = 0.0, double rmax = 1.0);

	// get normal pseudo random number with mean mu and std. deviation sdev
	double normrand(double mu = 0.0, double sdev = 1.0);

	// get poisson pseudo random number with mean mu
	int poisrand(double mu = 1.0);
};

