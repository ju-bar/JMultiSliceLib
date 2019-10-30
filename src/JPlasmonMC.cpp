//
// C++ source file: JPlasmonMC.cpp
// implementation for library JMultislice.lib (declarations see JPlasmonMC.h)
//
//
// Copyright (C) 2018, 2019 - Juri Barthel (juribarthel@gmail.com)
// Copyright (C) 2018, 2019 - Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
//
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

#include "stdafx.h"
#include "JPlasmonMC.h"
#include "NatureConstants.h"
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

using namespace std;


// ****************************************************************************
//
// Definition of class CJProbeParams
//

CJPlasmonMC::CJPlasmonMC()
{
	m_sca_num_max = 0;
	m_q_e = 0.f;
	m_q_c = 0.f;
	m_meanfreepath = 0.f;
	m_sca_num = 0;
	m_sca_tot_qx = 0.f;
	m_sca_tot_qy = 0.f;
}

CJPlasmonMC::CJPlasmonMC(const CJPlasmonMC & src)
{
	m_sca_num_max = src.m_sca_num_max;
	m_q_e = src.m_q_e;
	m_q_c = src.m_q_c;
	m_meanfreepath = src.m_meanfreepath;
	m_sca_num = src.m_sca_num;
	m_sca_tot_qx = src.m_sca_tot_qx;
	m_sca_tot_qy = src.m_sca_tot_qy;
}

CJPlasmonMC::~CJPlasmonMC()
{
	
}

void CJPlasmonMC::operator=(const CJPlasmonMC &src)
{
	m_sca_num_max = src.m_sca_num_max;
	m_q_e = src.m_q_e;
	m_q_c = src.m_q_c;
	m_meanfreepath = src.m_meanfreepath;
	m_sca_num = src.m_sca_num;
	m_sca_tot_qx = src.m_sca_tot_qx;
	m_sca_tot_qy = src.m_sca_tot_qy;
}

bool CJPlasmonMC::operator==(const CJPlasmonMC &src) const
{
	bool bResult = true;
	bResult &= (m_sca_num_max == src.m_sca_num_max);
	bResult &= (m_q_e == src.m_q_e);
	bResult &= (m_q_c == src.m_q_c);
	bResult &= (m_meanfreepath == src.m_meanfreepath);
	bResult &= (m_sca_num == src.m_sca_num);
	bResult &= (m_sca_tot_qx == src.m_sca_tot_qx);
	bResult &= (m_sca_tot_qy == src.m_sca_tot_qy);
	return bResult;
}

float CJPlasmonMC::UniRand(void)
{
	float rnd = (float)rand() / (float)(RAND_MAX + 1);
	return rnd;
}

UINT CJPlasmonMC::PoissonRand(float m)
{
	UINT u = 0;
	double dm, dl, dp;
	UINT k = 0;
	dm = (double)abs(m);
	// D. Knuth's algorithm is fast for small mean values < 30
	dl = exp(-dm);
	dp = 1.;
	while (dp > dl) {
		k++;
		dp = dp * (double)UniRand();
	}
	if (k > 0) {
		u = k - 1;
	}
	return u;
}

void CJPlasmonMC::Init(void)
{
	m_q_e2 = m_q_e * m_q_e;
	m_q_c2 = m_q_c * m_q_c;
	m_meanfreepath = abs(m_meanfreepath);
	ResetMC();
}

void CJPlasmonMC::ResetMC()
{
	m_sca_num = 0;
	m_sca_tot_qx = 0.f;
	m_sca_tot_qy = 0.f;
}

void CJPlasmonMC::ScatteringMC(float dt, UINT num_sca_max, UINT & num_sca, float & sca_qx, float & sca_qy)
{
	num_sca = 0;
	sca_qx = 0.f;
	sca_qy = 0.f;
	if (m_meanfreepath <= 0.f || m_q_e <= 0.f || m_q_c < m_q_e) return;
	// calculate a poissonian random number(limited to max.allowed excitation level)
	num_sca = min(PoissonRand(dt / m_meanfreepath), num_sca_max);
	if (num_sca > 0) {
		float pmc = 0.f;
		float qmc = 0.f;
		float emc = 0.f;
		for (UINT i = 0; i < num_sca; i++) {
			pmc = (float)(_TPI *UniRand());
			emc = UniRand();
			qmc = sqrt(m_q_e2 * pow(m_q_c2 / m_q_e2 + 1.f, emc) - m_q_e2);
			sca_qx += (qmc * cos(pmc));
			sca_qy += (qmc * sin(pmc));
		}
	}
}

int CJPlasmonMC::ScatGridMC(float dt, float a_x, float a_y, int *di, int *dj)
{
	float dqx = 0.f, dqy = 0.f;
	UINT nexmax = m_sca_num_max;
	UINT nex = m_sca_num;
	UINT nexadd = 0;
	*di = 0;
	*dj = 0;
	// get new scattering angle due to additional excitations
	ScatteringMC(dt, nexmax - nex, nexadd, dqx, dqy);
	if (nexadd > 0) { // excitation occurred
		int ntqi = (int)round(m_sca_tot_qx * a_x); // current hor. pixel shift
		int ntqj = (int)round(m_sca_tot_qy * a_y); // current ver. pixel shift
		m_sca_num += nexadd; // new excitation level
		m_sca_tot_qx += dqx; // new total scattering vector, x-component [1/nm]
		m_sca_tot_qy += dqy; // new total scattering vector, y-component [1/nm]
		*di = (int)round(m_sca_tot_qx * a_x) - ntqi; // change of hor. pixel shift
		*dj = (int)round(m_sca_tot_qy * a_y) - ntqj; // change of ver. pixel shift
	}
	return nexadd;
}