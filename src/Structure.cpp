// ************************************************************************** //
//
// FILE "Structure.cpp"
// SOURCEFILE
// AUTHOR Juri Barthel (ju.barthel@fz-juelich.de)
// DATE 2020/04/24
// Ernst Ruska-Centre for Microscopy and Spectroscopy with Electrons
//   Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
//
// LINK with "Atoms.h", "AtomSite.h"
//
// IMPLEMENTS class CStructure
//
// ************************************************************************** //
//

#include <fstream>
#include <algorithm>
#include "Structure.h"
#include "TextFileOps.h"
#include "JFFTMKLcore.h"



//
// ************************************************************************** //
//
// IMPLEMENTATION OF CLASS CStructure
//
// * CONSTRUCTORS *********************************************************** //

CStructure::CStructure(void)
{
	m_vCellSize = CV3D();
	m_vCellAngles = CV3D(90.0f, 90.0f, 90.0f);
	m_str_last_error = "";

}

CStructure::CStructure(const CStructure &other)
{
	m_vCellSize = CV3D();
	m_vCellAngles = CV3D(90.0f, 90.0f, 90.0f);
	m_str_last_error = "";
	
	CStructure *psrc = const_cast <CStructure*>(&other);
	m_vCellSize = psrc->m_vCellSize;
	m_vCellAngles = psrc->m_vCellAngles;
	m_vTypes = psrc->m_vTypes;
	m_vAtoms = psrc->m_vAtoms;
	m_vSites = psrc->m_vSites;
}

CStructure::~CStructure(void)
{
	m_vAtoms.clear();
	m_vTypes.clear();
	m_vSites.clear();
}


// * OPERATORS ************************************************************** //

void CStructure::operator=(CStructure other)
{
	DeleteStructureData();
	//
	// copy structure data
	m_vCellSize = other.m_vCellSize;
	m_vCellAngles = other.m_vCellAngles;
	m_vTypes = other.m_vTypes;
	m_vAtoms = other.m_vAtoms;
	m_vSites = other.m_vSites;
}


const bool CStructure::operator==(const CStructure other) const
{
	bool bResult = true;
	size_t i = 0, n = 0;
	// structure data
	bResult &= (m_vCellSize == other.m_vCellSize);
	bResult &= (m_vCellAngles == other.m_vCellAngles);
	bResult &= (n = m_vAtoms.size() == other.m_vAtoms.size());
	if (bResult && n > 0 ) { for (i=0; i<n; i++) { bResult &= (m_vAtoms[i]==other.m_vAtoms[i]); } }
	bResult &= (n = m_vTypes.size() == other.m_vTypes.size());
	if (bResult && n > 0 ) { for (i=0; i<n; i++) { bResult &= (m_vTypes[i]==other.m_vTypes[i]); } }
	// no need to compare m_vSites, since the information is redunant with m_vAtoms
	return bResult;
}


// * MEMBER FUNCTIONS ******************************************************* //

std::string CStructure::GetLastErrorMessage()
{
	return m_str_last_error;
}

void CStructure::DeleteStructureData(void)
{
	// reset the structure data
	m_vAtoms.clear();
	m_vTypes.clear();
	m_vSites.clear();
	// also reset error message
	m_str_last_error = "";
}

size_t CStructure::GetAtomNumber(int nType)
{
	size_t nResult = 0;
	size_t iType = (size_t)nType;
	if (iType >=0 && iType < m_vTypes.size()) {
		CAtomType att = m_vTypes.at(iType);
		nResult = GetAtomNumber(att);
	}
	return nResult;
}

size_t CStructure::GetAtomNumber(CAtomType att)
{
	size_t nResult = 0;
	size_t natty = m_vTypes.size();
	size_t nat = m_vAtoms.size();
	if (nat > 0 && natty > 0) {
		size_t nType = 0;
		for (size_t i = 0; i < nat; i++) {
			nType = m_vAtoms[i].m_nAtomTypeIndex;
			if (nType >= 0 && nType < natty) {
				if ( att == m_vTypes[nType] )
					nResult++;
			}
		}
	}
	return nResult;
}


size_t CStructure::GetMaxAtomNumber()
{
	size_t nResult = 0;
	size_t natty = m_vTypes.size();
	if (natty>0) {
		for (size_t i = 0; i < natty; i++) {
			nResult = std::max(nResult, GetAtomNumber((int)i));
		}
	}
	return nResult;
}

bool CStructure::GetAtomData(int iAtom, int &nAtomicNumber, float &fCharge, float &fBiso, float &fOccupancy, float &fPosX, float &fPosY, float &fPosZ)
{
	size_t nat = m_vAtoms.size();
	size_t natt = m_vTypes.size();
	if ( nat <= 0 || nat <= iAtom || iAtom < 0 ) return false; // evade empty atom list and invalid atom index
	CAtom *pat = &m_vAtoms[iAtom];
	if (NULL == pat) return false; // evade invalid atom data
	int iAtomType = pat->m_nAtomTypeIndex; // get atom type index
	if ( natt <= 0 || natt <= iAtomType || iAtomType < 0 ) return false; // evade empty atom type list and invalid atom type index
	CAtomType *patt = &m_vTypes[iAtomType];
	if (NULL == patt) return false; // evade invalid atom type data
	// transfer data
	nAtomicNumber = patt->m_nAtomicNumber;
	fCharge = patt->m_fIonicCharge;
	fBiso = patt->m_fBiso;
	fOccupancy = pat->m_fOccupancy;
	fPosX = pat->m_vPos.x;
	fPosY = pat->m_vPos.y;
	fPosZ = pat->m_vPos.z;
	return true;
}

float CStructure::GetVolume(void)
{
	float fVol = 0.0f;
	float fd2r = .01745329252f;
	float ca = cosf(m_vCellAngles.x*fd2r);
	float cb = cosf(m_vCellAngles.y*fd2r);
	float cc = cosf(m_vCellAngles.z*fd2r);
	float da = sqrtf(fabsf(1.0f-ca*ca-cb*cb-cc*cc+2.0f*ca*cb*cc));
	fVol = m_vCellSize.x * m_vCellSize.y * m_vCellSize.z * da;
	return fVol;
}

bool CStructure::IsOrthorhombic(void)
{
	bool bresult = false;
	bresult = (fabs(m_vCellAngles.x - 90.0f) < 0.5f && fabs(m_vCellAngles.y - 90.0f) < 0.5f && fabs(m_vCellAngles.z - 90.0f) < 0.5f);
	return bresult;
}


int CStructure::SuggestEquidistantZSlices(float dmin)
{
	int nslc = 1;
	int ndimh = 0, nnyq = 0;
	int i = 0;
	int pdims[1];
	float dstep = std::max(0.0001f, dmin);
	float hstep = 0.f, zcur = 0.f, dz = 0.f, zpos = 0.f, zamp = 0.f, rfac = -1.f / (dstep * dstep);
	float amp_max = 0.f;
	float* phist = NULL;
	float* pfreq = NULL;
	CJFFTMKLcore jco;
	size_t idx = 0, iat = 0, iatt = 0, nat = m_vAtoms.size();
	if (2 > nat) {
		goto _exit; // zero or one atom -> 1 slice
	}
	if (m_vCellSize.z <= 0.f) { // no cell
		nslc = 0;
		goto _exit;
	}
	// prepare histogram
	i = 8 * (1 + (int)(m_vCellSize.z / dstep)); // get minimum histogram size
	ndimh = std::max( 64, (int)std::pow(2.f, std::ceil( log((float)i) / log(2.0f) ) ) ); // extend to next power of 2 number
	nnyq = ndimh >> 1; // Nyquist frequency number
	pdims[0] = ndimh; // store size in array
	hstep = m_vCellSize.z / (float)ndimh; // set histogram sampling rate
	phist = (float*)malloc(sizeof(float) * ndimh); // allocate histogram
	memset(phist, 0, sizeof(float) * ndimh); // initialize histogram
	for (iat = 0; iat < nat; iat++) { // loop over all atoms
		iatt = m_vAtoms[iat].m_nAtomTypeIndex;
		zamp = (float)m_vTypes[iatt].m_nAtomicNumber; // atom core charge
		zpos = m_vAtoms[iat].m_vPos.z * m_vCellSize.z; // atom z position in nm
		// add an atom function
		for (i = 0; i < ndimh; i++) {
			idx = (size_t)i;
			zcur = (float)i * hstep; // current z position in nm
			dz = zcur - zpos; // distance to atom in nm
			phist[idx] += zamp * expf(rfac * dz * dz); // add a surrogate atomic function as Gaussian with amplitude of core charge and width of requested resolution
		}
	}
	jco.Init(1, pdims); // init the fft of size ndimh
	jco.SetDataRe(phist); // set the histogram to the real part
	jco.FT(); // transform to frequency space
	// get frequency data
	pfreq = (float*)malloc(sizeof(float) * ndimh); // allocate histogram
	memset(pfreq, 0, sizeof(float) * ndimh); // initialize histogram
	jco.GetDataAbs(pfreq);
	// analyze the amplitudes from frequency 2 on
	nslc = 1;
	amp_max = pfreq[1]; // preset with 1 slice
	for (i = 2; i < nnyq; i++) {
		idx = (size_t)i;
		if (pfreq[idx] > amp_max) { // update suggested period
			amp_max = pfreq[idx];
			nslc = i;
		}
	}
_exit:
	jco.Deinit();
	if (NULL != phist) { free(phist); }
	if (NULL != pfreq) { free(pfreq); }
	return nslc;
}

int CStructure::GetZSliceStructure(CStructure* pst, float fz_0, float fz_1)
{
	int nerr = 0;
	float fz = fz_1 - fz_0;
	float fz_cur = 0.f, fz_new = 0.f;
	size_t iat = 0, nat = m_vAtoms.size();
	CAtom at;
	if (NULL == pst) {
		nerr = 1; goto _exit;
	}
	if (fz <= 0.f) {
		nerr = 2; goto _exit;
	}
	pst->m_vCellSize = m_vCellSize; // copy x, y 
	pst->m_vCellSize.z = fz * m_vCellSize.z; // new z
	pst->m_vCellAngles = m_vCellAngles; // copy angles
	pst->m_vTypes = m_vTypes; // copy all atom types, even if not needed (this keeps up the atom type index in m_vAtoms)
	pst->m_vAtoms.clear();
	if (nat > 0) { // copy atoms if on range [fz_0, fz_1[
		for (iat = 0; iat < nat; iat++) {
			fz_cur = m_vAtoms[iat].m_vPos.z; // fractional z coordinate in current cell
			fz_cur = fmodf(fz_cur, 1.f); // get fractional
			if (fz_cur < 0.f) fz_cur += 1.f; // bring to positive interval 
			if (fz_cur >= fz_0 && fz_cur < fz_1) { // atom is in the slice range
				fz_new = (fz_cur - fz_0) / fz; // get atom fractional z position in the slice
				at = m_vAtoms[iat];
				at.m_vPos.z = fz_new;
				pst->m_vAtoms.push_back(at);
			}
		}
	}
_exit:
	return nerr;
}

std::string CStructure::GetCompositionString(void)
{
	std::string sResult = "";

	size_t natt = m_vTypes.size();
	size_t nats = m_vAtoms.size();

	std::vector<CAtomType> attComposition;
	size_t natcomp = 0;
	size_t iatcomp = 0;
	size_t fatcomp = 0;
	size_t iat = 0;
	size_t iatt = 0;
	bool bfound = false;
	bool binsert = false;

	CAtom atx;
	CAtomType attx;
	
	if (nats > 0 && natt > 0) { // structure contains atoms
		//
		// accumulate atom types in attComposition
		//
		for(iat = 0; iat < nats; iat++) { // loop over atoms
			atx = m_vAtoms[iat]; // get current atom
			iatt = atx.m_nAtomTypeIndex; // get atom type index of this atom
			attx = m_vTypes[iatt]; // get a copy of the atom type object
			//
			// check if this atom type is already in the composition list
			fatcomp = 0;
			bfound = false; // preset with "NO"
			natcomp = attComposition.size();
			if (natcomp>0) { // there are atom types in the composition
				for ( iatcomp=0; iatcomp<natcomp; iatcomp++ ) { // loop through ...
					if (attx.m_nAtomicNumber == attComposition[iatcomp].m_nAtomicNumber) { // found type with same atomic number
						fatcomp = iatcomp; // store index
						bfound = true; // mark found
						break; // stop searching
					}
				}
			}
			if (bfound) { // found atom type in the current composition list
				// update the type occupancy
				attComposition[iatcomp].m_fOccupancy += atx.m_fOccupancy;
			}
			else { // not found -> new atom type
				attx.m_fOccupancy = atx.m_fOccupancy; // set first occupancy
				// add to type composition list
				binsert = false; // preset insert index with "end"
				fatcomp = 0;
				if (natcomp>0) { // there are atom types in the composition
					for ( iatcomp=0; iatcomp<natcomp; iatcomp++ ) { // loop through ...
						// insert heavy atoms first
						if (attx.m_nAtomicNumber > attComposition[iatcomp].m_nAtomicNumber) { // found type with smaller atomic number
							fatcomp = iatcomp; // store index
							binsert = true;
							break; // stop searching
						}
					}
				}
				if (binsert) { // insert position found
					auto it = attComposition.begin();
					attComposition.insert(it+fatcomp, attx);
				}
				else { // no insert position found
					// just add
					attComposition.push_back(attx);
				}
			}
		}
	}
	//
	// create string from attComposition
	//
	natcomp = attComposition.size();
	if ( natcomp>0 ) { // composition present
		for ( iatcomp=0; iatcomp<natcomp; iatcomp++) {
			if ( iatcomp>0 ) sResult += " ";
			sResult += attComposition[iatcomp].SpeciesComposition();
		}
	}
	else { // empty structure
		sResult = "vacuum";
	}

	return sResult;
}

CM3D CStructure::CalculateBasis(void)
{
	CM3D basis;
	CV3D vAngs;
	float fd2r = .01745329252f;
	vAngs = m_vCellAngles*fd2r;
	basis.x.x = m_vCellSize.x;
	basis.x.y = m_vCellSize.y*cosf(vAngs.z);
	basis.y.y = m_vCellSize.y*sinf(vAngs.z);
	basis.x.z = m_vCellSize.z*cosf(vAngs.y);
	basis.y.z = m_vCellSize.z*( cosf(vAngs.x) - cosf(vAngs.y)*cosf(vAngs.z) )/sinf(vAngs.z);
	basis.z.z = m_vCellSize.z*sqrtf( sinf(vAngs.z)*sinf(vAngs.z) - cosf(vAngs.y)*cosf(vAngs.y) - cosf(vAngs.x)*cosf(vAngs.x)
            + 2.f*cosf(vAngs.x)*cosf(vAngs.y)*cosf(vAngs.z) ) /sinf(vAngs.z);
	//basis.Chop();
	return basis;
}


void CStructure::AddAtom(int nAtomicNumber, float fCharge, float fBiso, float fOccupancy, float fPosX, float fPosY, float fPosZ)
{
	CAtomType att;
	CAtom at;
	att.m_nAtomicNumber = nAtomicNumber;
	att.m_fIonicCharge = fCharge;
	att.m_fBiso = fBiso;
	at.m_fOccupancy = fOccupancy;
	at.m_vPos.x = fPosX;
	at.m_vPos.y = fPosY;
	at.m_vPos.z = fPosZ;
	//
	size_t natty = m_vTypes.size();
	bool bfound = false;
	size_t idxatty = 0;
	//
	if (natty>0) {
		for(size_t i=0; i<natty; i++) {
			if (m_vTypes[i] == att) {
				idxatty = i;
				bfound = true;
				break;
			}
		}
	}
	if (!bfound) {
		m_vTypes.push_back(att);
		idxatty = natty;
		natty++;
	}
	at.m_nAtomTypeIndex = (int)idxatty;
	m_vAtoms.push_back(at);
	m_vTypes[idxatty].m_fOccupancy += at.m_fOccupancy;
}


float CStructure::GetDistance(CV3D p1, CV3D p2)
{
	CV3D v = m_mCellBasis.VProduct(p2 - p1);
	return v.Length();
}

int CStructure::SetAtomSites(float thr_close)
{
	int nerr = 0;
	int iat_link = 0;
	size_t nat = m_vAtoms.size(), iat = 0, jat = 0, jatt = 0;
	size_t nsite = 0, isite = 0, fsite = 0;
	size_t nlink = 0, ilink = 0;
	bool bnew_site = true;
	float ftmp = 0.f, focc = 0.f;
	CV3D ptmp(0.f, 0.f, 0.f);
	CV3D pat(0.f, 0.f, 0.f);
	CAtomSite site; // temp. site
	m_vSites.clear(); // reset
	CalculateBasis(); // update the basis matrix from cell parameters
	if (nat > 0) {
		for (iat = 0; iat < nat; iat++) {
			bnew_site = true;
			pat = m_vAtoms[iat].m_vPos;
			nsite = m_vSites.size();
			if (nsite > 0) {
				for (isite = 0; isite < nsite; isite++) {
					if (GetDistance(pat, m_vSites[isite].m_vPos) < thr_close) {
						bnew_site = false;
						fsite = isite;
						break;
					}
				}
			}
			if (bnew_site) { // this atom opens a new site
				site.ResetData();
				iat_link = (int)iat;
				site.m_vatocc.push_back(iat_link);
				site.m_fOccupancy = m_vAtoms[iat].m_fOccupancy;
				site.m_vPos = pat;
				jatt = (size_t)m_vAtoms[iat].m_nAtomTypeIndex;
				site.m_fBiso = m_vTypes[jatt].m_fBiso;
				m_vSites.push_back(site);
			}
			else { // atom links to an existing site
				iat_link = (int)iat;
				m_vSites[fsite].m_vatocc.push_back(iat_link);
				// update effective position and biso of the site by weighted occupancy average
				ftmp = 0.f; focc = 0.f; ptmp = CV3D(0.f, 0.f, 0.f);
				nlink = m_vSites[fsite].m_vatocc.size();
				for (ilink = 0; ilink < nlink; ilink++) {
					jat = m_vSites[fsite].m_vatocc[ilink];
					jatt = (size_t)m_vAtoms[jat].m_nAtomTypeIndex;
					focc = focc + m_vAtoms[jat].m_fOccupancy; // accumulate occupancy
					ptmp = ptmp + m_vAtoms[jat].m_vPos * m_vAtoms[jat].m_fOccupancy; // accumulate effective position
					ftmp = ftmp + m_vTypes[jatt].m_fBiso * m_vAtoms[jat].m_fOccupancy; // accumulate effective biso
				}
				m_vSites[fsite].m_fOccupancy = focc; // set total occupancy
				m_vSites[fsite].m_vPos = ptmp / focc; // normalize to effective position
				m_vSites[fsite].m_fBiso = ftmp / focc; // normalize to effective biso
			}
		}
	}
	return nerr;
}

int CStructure::LoadCEL(std::string filename)
{
	int nerr = 0, ierr = 0;
	int i = 0, j = 0;
	CAtomType att;
	CAtom at;
	CStructure stTmp; // temp structure object for loading, will only replace current data if loading is successful
	std::string sprm = ""; // parameter as string
	std::string ssymb = "";
	std::string sline = "";
	std::vector<std::string> v_lines; // list of all lines read from the file
	size_t nlines = 0, natoms = 0;
	float vol = 0.f;
	ierr = read_textfile_lines(filename, v_lines);
	if (0 < ierr) {
		m_str_last_error = "Error: failed to open file " + filename + " for reading.";
		nerr = 1;
		goto _exit;
	}
	nlines = v_lines.size();
	if (nlines < 3) { // check that at least 3 lines are present (title, cell data, end of structure signal)
		m_str_last_error = "Error: failed to read CEL structure from file " + filename + ".";
		nerr = 2;
		goto _exit;
	}
	// parse the cell parameters from the second line of the file
	sline = v_lines[1];
	i = read_param(0, &sline, &sprm); // there should be a '0' leading input
	if (i > 0 && i < (int)sline.size()) { // read cell / lattice length constants a, b, c in nm
		i = read_param(i, &sline, &sprm);
		stTmp.m_vCellSize.x = to_float(sprm);
	}
	if (i > 0 && i < (int)sline.size()) {
		i = read_param(i, &sline, &sprm);
		stTmp.m_vCellSize.y = to_float(sprm);
	}
	if (i > 0 && i < (int)sline.size()) {
		i = read_param(i, &sline, &sprm);
		stTmp.m_vCellSize.z = to_float(sprm);
	}
	if (i > 0 && i < (int)sline.size()) { // read cell angles in deg
		i = read_param(i, &sline, &sprm);
		stTmp.m_vCellAngles.x = to_float(sprm);
	}
	if (i > 0 && i < (int)sline.size()) {
		i = read_param(i, &sline, &sprm);
		stTmp.m_vCellAngles.y = to_float(sprm);
	}
	if (i > 0 && i < (int)sline.size()) {
		i = read_param(i, &sline, &sprm);
		stTmp.m_vCellAngles.z = to_float(sprm);
	}
	stTmp.m_mCellBasis = CalculateBasis();
	vol = stTmp.GetVolume();
	if (vol <= 0.f) {
		m_str_last_error = "Error: the structure in file " + filename + " defines a cell of zero volume.";
		nerr = 4;
		goto _exit;
	}

	// parse all atom lines
	for (size_t iline = 2; iline < nlines - 2; iline++) {
		sline = v_lines[iline];
		if (std::string::npos != sline.find("*", 0)) break; // reached end of structure signal
		j = 0; // counts number of parameters successfully read
		i = read_param(0, &sline, &sprm);
		if (i >= 0 && i < (int)sline.size()) {
			att.m_nAtomicNumber = att.GetAtomicNumber(sprm);
			if (att.m_nAtomicNumber < 0) {
				m_str_last_error = "Error: failed to identify atom symbol '" + sprm + "' in line " + format("%d", iline) + ".";
				nerr = 5;
				goto _exit;
			}
			else {
				att.m_fIonicCharge = att.GetIonicCharge(sprm);
				j++;
			}
		}
		if (i >= 0 && i < (int)sline.size()) {
			i = read_param(i, &sline, &sprm);
			at.m_vPos.x = to_float(sprm);
			j++;
		}
		if (i >= 0 && i < (int)sline.size()) {
			i = read_param(i, &sline, &sprm);
			at.m_vPos.y = to_float(sprm);
			j++;
		}
		if (i >= 0 && i < (int)sline.size()) {
			i = read_param(i, &sline, &sprm);
			at.m_vPos.z = to_float(sprm);
			j++;
		}
		if (i >= 0 && i < (int)sline.size()) {
			i = read_param(i, &sline, &sprm);
			at.m_fOccupancy = to_float(sprm);
			j++;
		}
		if (i >= 0 && i < (int)sline.size()) {
			i = read_param(i, &sline, &sprm);
			att.m_fBiso = to_float(sprm); // nm
			j++;
		}
		if (j != 6) {
			m_str_last_error = "Error: failed to parse line " + format("%d", iline) + " for atomic data.";
			nerr = 6;
			goto _exit;
		}
		stTmp.AddAtom(att.m_nAtomicNumber, att.m_fIonicCharge, att.m_fBiso, at.m_fOccupancy, at.m_vPos.x, at.m_vPos.y, at.m_vPos.z);
	}

_exit:
	if (nerr == 0) { // loading successful, copy data to this object, replaces current data 
		*this = stTmp;
	}
	return nerr;
}


int CStructure::SaveCEL(std::string filename)
{
	int nerr = 0;
	std::string stmp = "";
	std::string sline = "";
	std::ofstream fcout; // file stream
	fcout.open(filename, std::ofstream::trunc);  // open the file
	size_t natoms = m_vAtoms.size();
	size_t iat = 0, iatt = 0;
	if (!fcout.is_open()) {
		nerr = 1;
		goto _exit;
	}
	// header info line (comment)
	stmp = GetCompositionString();
	sline = "# SUPERCELL of " + stmp + " (file generated by JMultislice/CStructure)";
	fcout << sline << std::endl;
	// supercell size and angles
	sline = format(" 0 %9.5f %9.5f %9.5f %8.4f %8.4f %8.4f", m_vCellSize.x, m_vCellSize.y, m_vCellSize.z, m_vCellAngles.x, m_vCellAngles.y, m_vCellAngles.z);
	fcout << sline << std::endl;
	if (0 < natoms) { // write atom data
		for (iat = 0; iat < natoms; iat++) {
			iatt = m_vAtoms[iat].m_nAtomTypeIndex;
			stmp = m_vTypes[iatt].SpeciesName();
			sline = " " + stmp + "  " + format("%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f",
				m_vAtoms[iat].m_vPos.x, m_vAtoms[iat].m_vPos.y, m_vAtoms[iat].m_vPos.z,
				m_vAtoms[iat].m_fOccupancy, m_vTypes[iatt].m_fBiso, 0.f, 0.f, 0.f);
			fcout << sline << std::endl;
		}
	}
	// write end of structure signal
	fcout << "*" << std::endl;
	// close file
	fcout.close();
_exit:
	return nerr;
}

// ************************************************************************** //