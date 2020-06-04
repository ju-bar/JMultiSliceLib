// ************************************************************************** //
//
// FILE "AtomSite.cpp"
// SOURCEFILE
// AUTHOR Juri Barthel (ju.barthel@fz-juelich.de)
// DATE 2015/01/29
// Ernst Ruska-Centre for Microscopy and Spectroscopy with Electrons
//   Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
//
// LINK with "StructureParams.h"
//
// IMPLEMENTS class AtomSite and CStructure
//
// ************************************************************************** //
//

#include "AtomSite.h"

// ************************************************************************** //
//
// IMPLEMENTATION OF CLASS CATOMSITE
//
// * CONSTRUCTORS *********************************************************** //

CAtomSite::CAtomSite(void)
{
	// standard constructor, preset with defaults
	ResetData();
}

CAtomSite::CAtomSite(const CAtomSite &atsSrc)
{
	ResetData();
	//
	CAtomSite *pats = const_cast <CAtomSite*>(&atsSrc);
	m_vatocc = pats->m_vatocc;
	m_vPos = pats->m_vPos;
	m_fOccupancy = pats->m_fOccupancy;
	m_fBiso = pats->m_fBiso;
	// The selection state is not copied.
}

// * DECONSTRUCTORS ********************************************************* //

CAtomSite::~CAtomSite(void)
{
	ResetData();
}

// * OPERATORS ************************************************************** //

void CAtomSite::operator=(const CAtomSite &other)
{
	ResetData();
	//
	CAtomSite *pats = const_cast <CAtomSite*>(&other);
	m_vatocc = pats->m_vatocc;
	m_vPos = pats->m_vPos;
	m_fOccupancy = pats->m_fOccupancy;
	m_fBiso = pats->m_fBiso;
	// The selection state is not copied.
}

const bool CAtomSite::operator==(const CAtomSite &other) const
{
	bool bResult = true;
	CAtomSite *pats = const_cast <CAtomSite*>(&other);
	bResult &= pats->EqualOccupationList(pats);
	CV3D vDif = m_vPos - other.m_vPos;
	bResult &= ( ATOM_POSDIF_THR > vDif.Length() );
	// The selection state is not compared.

	return bResult;
}

// * MEMBER FUNCTIONS ******************************************************* //

void CAtomSite::ResetData(void)
{
	m_vatocc.clear();
	m_fOccupancy = 0.0f;
	m_fBiso = 0.0f;
	m_vPos = CV3D();
	m_bSelected = false; // The selection state is reset as well.
}



bool CAtomSite::EqualOccupationList(CAtomSite *pOther)
{
	bool bResult = true;
	size_t n1 = m_vatocc.size();
	size_t n2 = pOther->m_vatocc.size();
	bResult &= (n1==n2); // check same list size
	if (bResult &&  n1>0) {
		for (size_t i = 0; i < n1; i++) { // check same list
			bResult &= (pOther->m_vatocc[i] == m_vatocc[i]);
		}
	}
	return bResult;
}


bool CAtomSite::IsDefined(void)
{
	if (0 == m_vatocc.size())
		return false;
	return true;
}


bool CAtomSite::IsSelected(void)
{
	return m_bSelected;
}

bool CAtomSite::Select(void)
{
	bool bIS = m_bSelected;
	m_bSelected = true;
	return (!bIS);
}

bool CAtomSite::Deselect(void)
{
	bool bIS = m_bSelected;
	m_bSelected = false;
	return bIS;
}

bool CAtomSite::SwitchSelection(void)
{
	m_bSelected = !m_bSelected;
	return m_bSelected;
}


// ************************************************************************** //



