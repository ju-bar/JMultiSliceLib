// file: 'rng.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementation of the class CRng to be used as
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

#include "rng.h"
#include <time.h>

CRng::CRng() :
	m_seed((int)::time(NULL)), 
	urng(rd()), 
	uni_dist((std::uniform_int_distribution<>(0, RAND_MAXROLL))),
	norm_dist((std::normal_distribution<>(0., 1.)))
{
	urng.seed(m_seed);
}


CRng::CRng(int seed) : urng(rd()), uni_dist((std::uniform_int_distribution<>(0, RAND_MAXROLL)))
{
	m_seed = seed;
	urng.seed(m_seed);
}

CRng::~CRng()
{
}

void CRng::seed(int seed)
{
	if (seed == 0) {
		m_seed = (int)::time(NULL);
	}
	else {
		m_seed = seed;
	}
	urng.seed(m_seed);
}

int CRng::get_seed()
{
	return m_seed;
}

int CRng::rand()
{
	return uni_dist(urng);
}

double CRng::unirand(double rmin, double rmax)
{
	return rmin + ((double)(rand()) / (double)((unsigned long)RAND_MAXROLL)) * (rmax - rmin);
}

double CRng::unirand_oub(double rmin, double rmax)
{
	return rmin + ((double)(rand()) / (double)((unsigned long)RAND_MAXROLL+1)) * (rmax - rmin);
}

double CRng::normrand(double mu, double sdev)
{
	return mu + sdev * norm_dist(urng);
}

int CRng::poisrand(double mu)
{
	std::poisson_distribution<int> pd(mu);
	return pd(urng);
}