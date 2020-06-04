// ************************************************************************** //
//
// FILE "StructureParams.h"
// HEADERFILE for "StructureParams.cpp"
// AUTHOR Juri Barthel (ju.barthel@fz-juelich.de)
// DATE 2020/04/24
// Ernst Ruska-Centre for Microscopy and Spectroscopy with Electrons
//   Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
//
// LINK with "Atoms.h", "AtomSite.h", "M3D.h"
//
// DECLARES class CStructure
//
// ************************************************************************** //
//


#include "Atoms.h"
#include "AtomSite.h"
#include "M3D.h"

//
// DECLARATION of class CStructure
//
#ifndef STRUCTURE_H
#define STRUCTURE_H

#define STRUCTURE_TWOPI			6.2831853072
#define STRUCTURE_PIH			1.5707963268

#define STRUCTURE_FILE_FORMAT_UNKNOWN	0
#define STRUCTURE_FILE_FORMAT_CEL		1
#define STRUCTURE_FILE_FORMAT_CIF		2
#define STRUCTURE_FILE_FORMAT_XTL		3
#define STRUCTURE_FILE_FORMAT_MAX		1

class CStructure
{
	// constructors
public:
	CStructure(void);
	CStructure(const CStructure &other); // copy
	// destructor
	~CStructure(void);

	// member variables

public:
	CV3D m_vCellSize; // (a,b,c) size of the structure super-cell [nm]
	CV3D m_vCellAngles; // (alpha,beta,gamma) angles of the super-cell [deg]
	CM3D m_mCellBasis; // Super-cell basis vectors
	
	std::vector<CAtomType> m_vTypes; // list of atom types as received, always add new atoms at the end
	std::vector<CAtom> m_vAtoms; // list of atoms as received, always add new atoms at the end
	std::vector<CAtomSite> m_vSites; // list of atom sites, linking to m_vAtoms (call SetAtomSites() to refresh from current m_vAtoms)
	

protected:
	std::string m_str_last_error; // contains information on the last error message

	// operators

public:
	void operator=(CStructure other); // copies data
	const bool operator==(const CStructure other) const; // compares data

	// member functions

public:
	void DeleteStructureData(void); // Deletes all structure data content.
	void AddAtom(int nAtomicNumber, float fCharge, float fBiso, float fOccupancy, float fPosX, float fPosY, float fPosZ); // Adds an atom to the lists
	bool GetAtomData(int iAtom, int &nAtomicNumber, float &fCharge, float &fBiso, float &fOccupancy, float &fPosX, float &fPosY, float &fPosZ); // retrieves data about an atom in the atom list
	size_t GetAtomNumber(int nType); // Returns the number of stored atoms of the given type index.
	size_t GetAtomNumber(CAtomType att); // Returns the number of stored atoms of the given type.
	size_t GetMaxAtomNumber(void); // returns the number of atoms of the atom type which has a maximum of atoms in the structure
	float GetVolume(void); // determines the volume of the structure
	float GetDistance(CV3D p1, CV3D p2); // returns the distance between the two fractional positions p1 and p2.
	std::string GetCompositionString(void); // returns a string with the current composition in the form "At1_Occ1 At2_Occ2 ... "
	CM3D CalculateBasis(void); // calculates the cell basis vectors from the current cell paramater m_vCellSize and m_vCellAngles
	bool IsOrthorhombic(void); // returns true if the cell has only 90 degree angles
	int SuggestEquidistantZSlices(float dmin); // determines an appropriate number of equidistant slices along z from the structure model
	int GetZSliceStructure(CStructure* pst, float fz_0, float fz_1); // fills pst with structure data from a z-slice of the current structure object with fractional coordinates in [fz_0, fz_1[ (excludes fz_1).
	int SetAtomSites(float thr_close = 0.05); // Refreshes the list of atom sites (m_vSites) from the current list of atoms (m_vAtoms). Close atoms are related to one site. Close means having a distance below thr_close in nm.
	
	int LoadCEL(std::string filename); // loads structure data from a file in CEL format
	int SaveCEL(std::string filename); // saves structure data to a file in CEL format
	
	std::string GetLastErrorMessage(); // returns the last error message

protected:

};
//
#endif // STRUCTURE_H