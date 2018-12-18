//
// JMultiSliceLib.cpp
// compile with: cl /c /EHsc JMultiSliceLib.cpp
// post-build command: lib JMultiSliceLib.obj
//
// This defines library wrapper function around class CJMultiSlice
//
// Copyright (C) 2018 - Juri Barthel (juribarthel@gmail.com)
// Copyright (C) 2018 - Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
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
#include <stdexcept>
#include "JMultiSliceLib.h"


using namespace std;

float __stdcall GetJMSVersion(void)
{
	return (float)__JMS_VERSION__ + 0.1f*(float)__JMS_VERSION_SUB__ + 0.01f*(float)__JMS_VERSION_SUB_SUB__;
}

int __stdcall SetJMSDebugLevel(int debuglevel)
{
	return JMS.SetDebugLevel(debuglevel);
}

int __stdcall SetJMSRNGSeed(int rngseed)
{
	return JMS.SetRNGSeedEx(rngseed);
}

float __stdcall SetHighTension(float ht)
{
	return JMS.SetHighTension(ht);
}

float __stdcall GetHighTension(void)
{
	return JMS.GetHighTension();
}

float __stdcall GetWaveLength(void)
{
	return JMS.GetWaveLength();
}

int __stdcall GetDetetionSliceNum(void)
{
	return JMS.GetDetetionSliceNum();
}

int __stdcall GetDetNum(void)
{
	return JMS.GetDetNum();
}
void __stdcall SetGridSize(int nx, int ny)
{
	JMS.SetGridSize(nx, ny);
}

void __stdcall GetGridSize(int &nx, int &ny)
{
	JMS.GetGridSize(nx, ny);
}

void __stdcall SetSupercellSize(float *a0)
{
	JMS.SetSupercellSize(a0);
}

void __stdcall SetSupercellSizeABC(float a, float b, float c)
{
	JMS.SetSupercellSize(a, b, c);
}

void __stdcall SetSliceThickness(int islc, float fthickness)
{
	JMS.SetSliceThickness(islc, fthickness);
}

void __stdcall DiffractionDescan(bool bActivate)
{
	JMS.DiffractionDescan(bActivate);
}

void __stdcall SetDiffractionDescanN(int whichcode, int ndescanx, int ndescany, int iThread)
{
	JMS.SetDiffractionDescanN(whichcode, ndescanx, ndescany, iThread);
}

void __stdcall SetDiffractionDescan(int whichcode, float descanx, float descany, int iThread)
{
	JMS.SetDiffractionDescan(whichcode, descanx, descany, iThread);
}

void __stdcall SetDiffractionDescanMRad(int whichcode, float descanx, float descany, int iThread)
{
	JMS.SetDiffractionDescanMRad(whichcode, descanx, descany, iThread);
}

int __stdcall PhaseGratingSetup(int whichcode, int nx, int ny, int nslc, int nvarmax, int* nslcvar)
{
	return JMS.PhaseGratingSetup(whichcode, nx, ny, nslc, nvarmax, nslcvar);
}

int __stdcall ObjectSliceSetup(int nobjslc, int* objslc)
{
	return JMS.ObjectSliceSetup(nobjslc, objslc);
}

int __stdcall PropagatorSetup(int whichcode, int npro, int *proidx)
{
	return JMS.PropagatorSetup(whichcode, npro, proidx);
}

int __stdcall DetectorSetup(int whichcode, int ndetper, int ndetint, int imagedet, int nthreads_CPU)
{
	return JMS.DetectorSetup(whichcode, ndetper, ndetint, imagedet, nthreads_CPU);
}

int __stdcall SetSlicePhaseGratings(int whichcode, int islc, int nvar, fcmplx* pgr)
{
	return JMS.SetPhaseGratingData(whichcode, islc, nvar, pgr);
}

int __stdcall SetPropagatorData(int whichcode, int ipro, fcmplx* pro)
{
	return JMS.SetPropagatorData(whichcode, ipro, pro);
}

int __stdcall SetDetectorData(int whichcode, int idet, float* det, int msklen, int* msk)
{
	return JMS.SetDetectorData(whichcode, idet, det, msklen, msk);
}

int __stdcall InitCore(int whichcode, int nCPUthreads)
{
	return JMS.InitCore(whichcode, nCPUthreads);
}

int __stdcall SetIncidentWave(int whichcode, fcmplx* wav)
{
	return JMS.SetIncidentWave(whichcode, wav);
}

int __stdcall GetUnscrambleHash(UINT* phash)
{
	return JMS.GetUnscrambleHash(phash);
}
int __stdcall Cleanup(void)
{
	return JMS.Cleanup();
}

void __stdcall CleanFFTW(void)
{
	JMS.CleanFFTW();
	return;
}


////////////////////////////////////////////////////////////////////////////////
//
// CALCULATIONS
//
////////////////////////////////////////////////////////////////////////////////

int __stdcall OffsetIncomingWave(int whichcode, float dx, float dy, float dz, int iThread)
{
	return JMS.OffsetIncomingWave(whichcode, dx, dy, dz, iThread);
}

int __stdcall CalculateProbeWaveFourier(CJProbeParams* prm, fcmplx *wav)
{
	return JMS.CalculateProbeWaveFourier(prm, wav);
}

int __stdcall CalculatePropagator(float fthick, float otx, float oty, fcmplx *pro, int ntype)
{
	return JMS.CalculatePropagator(fthick, otx, oty, pro, ntype);
}

int __stdcall CalculateRingDetector(float beta0, float beta1, float phi0, float phi1, float theta0x, float theta0y, std::string sdsprofile, float *det, int &msklen, int *msk)
{
	return JMS.CalculateRingDetector(beta0, beta1, phi0, phi1, theta0x, theta0x, sdsprofile, det, msklen, msk);
}

int __stdcall GetResult(int whichcode, int whichresult, float *dst, int iThread)
{
	return JMS.GetResult(whichcode, whichresult, dst, iThread);
}

////////////////////////////////////////////////////////////////////////////////
//
// CPU CODE
//
////////////////////////////////////////////////////////////////////////////////

int __stdcall GetCPUNum(void)
{
	return JMS.GetCPUNum();
}

int __stdcall SetIncomingWaveCPU(fcmplx* wav, bool bTranspose, int iThread)
{
	return JMS.SetIncomingWaveCPU(wav, bTranspose, iThread);
}

int __stdcall CPUMultislice(int islc0, int accmode, float weight, int iThread)
{
	return JMS.CPUMultislice(islc0, accmode, weight, iThread);
}

int __stdcall ClearDetMem_h(int iThread)
{
	return JMS.ClearDetMem_h(iThread);
}

////////////////////////////////////////////////////////////////////////////////
//
// GPU CODE
//
////////////////////////////////////////////////////////////////////////////////

int __stdcall GetGPUNum(void)
{
	return JMS.GetGPUNum();
}

int __stdcall GetGPUName(int idev, char* name)
{
	return JMS.GetGPUName(idev, name);
}

int __stdcall GetCurrentGPU(void)
{
	return JMS.GetCurrentGPU();
}

int __stdcall SetCurrentGPU(int idev)
{
	return JMS.SetCurrentGPU(idev);
}

int __stdcall GetGPUStats(int idev, int &iCMajor, int &iCMinor, int &iMaxThread, int64_t &CUDAmemtotal, int64_t &CUDAmemfree)
{
	return JMS.GetGPUStats(idev, iCMajor, iCMinor, iMaxThread, CUDAmemtotal, CUDAmemfree);
}

int __stdcall GetGPUCores(int idev, int &nMultiProc, int &nCores, int& nMaxThreadPerProc)
{
	return JMS.GetGPUCores(idev, nMultiProc, nCores, nMaxThreadPerProc);
}

int __stdcall SetIncomingWaveGPU(fcmplx* wav, bool bTranspose)
{
	return JMS.SetIncomingWaveGPU(wav, bTranspose);
}

int __stdcall GPUMultislice(int islc0, int accmode, float weight)
{
	return JMS.GPUMultislice(islc0, accmode, weight);
}

int __stdcall ClearDetMem_d(void)
{
	return JMS.ClearDetMem_d();
}