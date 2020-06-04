// ************************************************************************** //
//
// FILE "AtomSite.h"
// HEADERFILE for "AtomSite.cpp"
// AUTHOR Juri Barthel (ju.barthel@fz-juelich.de)
// DATE 2015/01/29
// Ernst Ruska-Centre for Microscopy and Spectroscopy with Electrons
//   Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
//
// LINK with "Atoms.h"
//
// DECLARES class CAtomSite
//
// ************************************************************************** //
//
// 2020/05/19 (JB)
// - added a member variable to hold an effective biso
//
// ************************************************************************** //

#include <vector>
#include "Atoms.h"


//
// DECLARATION of class CAtomSite
//
#ifndef ATOMSITE_H
#define ATOMSITE_H

#define ATOMSITEDISTANCE				0.02 // minimum distance to distinguish atomic sites [nm]

class CAtomSite
{

	// constructors
public:

	CAtomSite(void); // standard
	CAtomSite(const CAtomSite &atsSrc); // copy constructor
	// deconstructor
	~CAtomSite(void); // standard

	// data members
public:

	std::vector<int> m_vatocc; // list of atoms occupying this site (managed)
	float m_fOccupancy; // total site occupancy by atoms
	float m_fBiso; // effective Biso parameter
	CV3D m_vPos; // average position of the atomic site

protected:

	// Selection flag
	bool m_bSelected;

	// operators
public:

	void operator=(const CAtomSite &other); // copies data
	const bool operator==(const CAtomSite &other) const; // compares data
	
	// member functions
public:
	
	void ResetData(void); // Resets data to default as with standard constructor
	bool EqualOccupationList(CAtomSite *pOther); // Returns TRUE if the occupation list is identical

	bool IsDefined(void); // Returns TRUE if the atom site is defined by containing atoms.
	
	bool IsSelected(void); // returns the selection state of the atom site.
	bool Select(void); // Sets the selection state of the atom site to selected. Returns TRUE if the state was changed.
	bool Deselect(void); // Sets the selection state of the atom site to deselected. Returns TRUE if the state was changed.
	bool SwitchSelection(void); // Switches the selection state from TRUE to FALSE or from FALSE to TRUE. Returns the new selection state of the atom site object.

};

#endif // ATOMSITE_H

// ************************************************************************** //
