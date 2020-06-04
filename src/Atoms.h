// ************************************************************************** //
//
// FILE "Atoms.h"
// HEADERFILE for "Atoms.cpp"
// AUTHOR Juri Bartthel (ju.barthel@fz-juelich.de)
// DATE 2015/02/19
// Ernst Ruska-Centre for Microscopy and Spectroscopy with Electrons
//   Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
//
// LINK with "Params.h" and "V3D.h"
//
// DECLARES classes CAtomType, CAtom
//
// ************************************************************************** //
//
#pragma once

#include <string>
#include "fcomplex.h"
#include "V3D.h"

//
// DECLARATION of class CAtomType
//
#ifndef ATOMTYPE_H
#define ATOMTYPE_H

// Define atom and atom type parameter difference thresholds.
// Differences smaller than those thresholds are not recognized.
// In other words, atoms and types with parameter differences
// below these (default) thresholds will be treated as equal.
#define ATOMTYPE_CHRGDIF_THR			0.001
#define ATOMTYPE_BISODIF_THR			0.00001
#define ATOMTYPE_OCCUDIF_THR			0.001
#define ATOMTYPE_ZNUM_MAX				128
#define ATOMTYPE_ZSYM_LEN				4
#define ATOM_POSDIF_THR					0.00001

#ifndef ATOMTYPE_SYMBOLS
#define ATOMTYPE SYMBOLS
static char sAtomSymbols[ATOMTYPE_ZNUM_MAX+1][ATOMTYPE_ZSYM_LEN+1] = {
	"Void",														// 0
	"H" , "He", "Li", "Be", "B",  "C",  "N",  "O", "F",  "Ne",	// 10
	"Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", // 20
	"Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", //
	"Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", //
	"Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", // 50
	"Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", //
	"Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", //
	"Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", //
	"Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", //
	"Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", // 100 // scattering tables go up to Cf (Z=98)
	"Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", //
	"Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo", "Uun", "Udc", // 120
	"Udu", "Udd", "Udt", "Udq", "Udp", "Udh", "Uds", "Udo"				// 128
};
#endif // ATOMTYPE_SYMBOLS

class CAtomType
{
public:
	// constructors
	CAtomType(void); // standard
	CAtomType(const CAtomType &attSrc); // copy constructor
	// deconstructor
	~CAtomType(void); // standard

public:
	// data members
	int m_nAtomicNumber; // atomic number, Z
	float m_fIonicCharge; // ionic charge (e)
	float m_fBiso; // Biso (nm^2)
	float m_fOccupancy; // occupancy (sum of atom occupancies with this type)
	

protected:
	size_t nk_scaf; // number of sampled for prepared scattering factors: == 0 if not prepared
	float sk_scaf; // spatial frequency sampling rate of the scattering factor table (1/nm/sample)
	fcmplx* scaf; // array holding prepared scattering factors


public:
	// operators
	void operator=(const CAtomType &other); // copies data
	const bool operator==(const CAtomType &other) const; // compares data

public:
	// member functions

	// Resets data to default as with standard constructor
	void ResetData();

	// Returns the atomic number for an atom symbol given as string ("H" -> 1, "Ti" -> 22)
	// Returns -1 on failure
	int GetAtomicNumber(std::string symbol);

	// Returns the ionic charge from an atom symbol ("Sr2+" -> 2.0f, "O2-" -> -2.0f)
	float GetIonicCharge(std::string symbol);

	// Checks if the two atom types are of the same species (atomic number and ionic charge).
	bool SameSpecies(CAtomType *pAtt);

	// Checks if the two atom types are of the same species and Biso parameter (atomic number and ionic charge).
	bool SameSpeciesBiso(CAtomType *pAtt);
	
	// Returns the name string of the species (symbol + charge, if charged)
	std::string SpeciesName(void);

	// Returns species name + occupancy sub-suffix
	std::string SpeciesComposition(void);

	// Calculates radial form factors for this atom type
	// - with nk samples and k-space sampling rate sk (1/nm/sample)
	// - for band-width limit bwlk (1/nm)
	// - for kinetic energy of electrons ekv (keV)
	// - using Debye-Waller factors if dwf == true
	// - with absorption model switch absorb and fix absorption factor abf if absorb == 1
	// - atomic form factors from table aaffacs
	// - bverbose switches processing output to std::cout
	int CalculateRadialFormFactors(size_t nk, float sk, float bwlk, float ekv, bool dwf, unsigned int absorb, float abf, unsigned int atffacs, bool bverbose = false);

	// Copies the form fractor buffer to the buffer adressed by ptr_buff and returns the size in fcmplx itmes.
	// The buffer pointed to by ptr_buf will be allocated on heap to the same size as the scaf member if this object.
	// Also copies the sampling rate to the pointed adress sk
	size_t CopyFormFactors(void** ptr_buf, float* psk) const;

	// Returns the radial form factor value at a given spatial frequency q (1/nm) from the precalulated table
	// - call CalculateRadialFormFactors before using this function
	// - ipo : order of polynomial interpolation
	fcmplx GetFormFactor(float q, unsigned int ipo = 1);

private:

	// Internal plynomial interpolation of scattering factors
	fcmplx InterpolateScaf(float q, unsigned int ipo = 3);
};

#endif // ATOMTYPE_H


//
// DECLARATION of class CAtom
//
#ifndef ATOM_H
#define ATOM_H

class CAtom
{
	public:
	// constructors
	CAtom(void); // standard
	CAtom(const CAtom &atSrc); // copy constructor
	// deconstructor
	~CAtom(void); // standard

public:
	// data members
	int m_nAtomTypeIndex; // links to a member in a list of CAtomType
	float m_fOccupancy; // site occupancy
	CV3D m_vPos; // fractional cell position

public:
	// operators
	void operator=(const CAtom &other); // copies data
	const bool operator==(const CAtom &other) const; // compares data

public:
	// member functions
	void ResetData(); // Resets data to default as with standard constructor
};

#endif // ATOM_H
//
// ************************************************************************** //