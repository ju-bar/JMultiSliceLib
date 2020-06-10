// ************************************************************************** //
//
// FILE "Atoms.cpp"
// SOURCEFILE
// AUTHOR Juri Bartthel (ju.barthel@fz-juelich.de)
// DATE 2015/02/19
// Ernst Ruska-Centre for Microscopy and Spectroscopy with Electrons
//   Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
//
// IMPLEMENTS classes CAtomType and CAtom
//
// ************************************************************************** //
//

#include "Atoms.h"
#include "TextFileOps.h"
#include "wekoscat.h"
#include <algorithm>

//
// ************************************************************************** //
//
// IMPLEMENTATION OF CLASS CATOMTYPE
//
// * CONSTRUCTORS *********************************************************** //

CAtomType::CAtomType(void)
{
	// standard constructor, preset with defaults
	//
	scaf = NULL;
	nk_scaf = 0;
	sk_scaf = 0.f;
	ResetData();
}

CAtomType::CAtomType(const CAtomType &attSrc)
{
	scaf = NULL;
	nk_scaf = 0;
	sk_scaf = 0.f;
	ResetData();
	//
	m_nAtomicNumber = attSrc.m_nAtomicNumber;
	m_fIonicCharge = attSrc.m_fIonicCharge;
	m_fBiso = attSrc.m_fBiso;
	m_fOccupancy = attSrc.m_fOccupancy;
	//
	nk_scaf = attSrc.CopyFormFactors((void**)&scaf, &sk_scaf);
}

// * DECONSTRUCTORS ********************************************************* //

CAtomType::~CAtomType(void)
{
	if (scaf != NULL) {
		free(scaf);
		scaf = NULL;
	}
	nk_scaf = 0;
}

// * OPERATORS ************************************************************** //

void CAtomType::operator=(const CAtomType &other)
{
	//
	m_nAtomicNumber = other.m_nAtomicNumber;
	m_fIonicCharge = other.m_fIonicCharge;
	m_fBiso = other.m_fBiso;
	m_fOccupancy = other.m_fOccupancy;
	//
	nk_scaf = other.CopyFormFactors((void**)&scaf, &sk_scaf);
}

const bool CAtomType::operator==(const CAtomType &other) const
{
	bool bResult = true;

	bResult &= (m_nAtomicNumber == other.m_nAtomicNumber);
	bResult &= ((float)ATOMTYPE_CHRGDIF_THR > fabsf( m_fIonicCharge - other.m_fIonicCharge) );
	bResult &= ((float)ATOMTYPE_BISODIF_THR > fabsf( m_fBiso - other.m_fBiso) );
	//bResult &= ((float)ATOMTYPE_OCCUDIF_THR > fabsf( m_fOccupancy - other.m_fOccupancy ) );

	return bResult;
}

// * MEMBER FUNCTIONS ******************************************************* //

void CAtomType::ResetData(void)
{
	m_nAtomicNumber = 0;
	m_fIonicCharge = 0.0f;
	m_fBiso = 0.0f;
	m_fOccupancy = 0.0f;
	// reset protected buffer for form factors
	if (scaf != NULL) {
		free(scaf);
		scaf = NULL;
	}
	nk_scaf = 0;
	sk_scaf = 0.f;
}


size_t CAtomType::CopyFormFactors(void** ptr_buf, float* psk) const
{
	fcmplx* dst = NULL;
	size_t mem = 0;
	if (NULL == ptr_buf || NULL == psk) return 0; // theres no adress given to the destination buffer
	dst = (fcmplx*)*ptr_buf; 
	if (NULL != dst) { // free existing memory allocation for destination buffer
		free(dst);
		dst = NULL;
	}
	if (0 == nk_scaf || NULL == scaf) { // source buffer is not allocated
		dst = NULL; // return with no allocated destination
		*psk = 0.f;
		return 0;
	}
	// we get here if the source buffer is allocated and nk_scaf > 0
	mem = sizeof(fcmplx) * nk_scaf; // get buffer size
	dst = (fcmplx*)malloc(mem); // allocate new memory for destination buffer
	memcpy(dst, scaf, mem); // copy data
	*psk = sk_scaf; // copy sampling rate
	return nk_scaf;
}

int CAtomType::GetAtomicNumber(std::string symbol)
{
	int z = -1, iz = 0;
	std::string scmp = "";
	std::string snum = "012345689";
	size_t pnum = symbol.find_first_of(snum, 0); // find first numerical character
	if (pnum == symbol.npos) pnum = symbol.size();
	bool bfound = false;
	while (!bfound && iz < (int)ATOMTYPE_ZNUM_MAX) {
		scmp = sAtomSymbols[iz];
		if (std::string::npos != symbol.find(scmp, 0) && pnum == scmp.size()) { // symbol found
			bfound = true;
			z = iz;
			break;
		}
		iz++; // symbol not found, try next
	}
	return z;
}

float CAtomType::GetIonicCharge(std::string symbol)
{
	float crg = 0.f;
	int iz = 0;
	std::string scmp = "";
	std::string snum = "012345689";
	size_t pnum = symbol.find_first_of(snum, 0); // find first numerical character
	if (pnum == symbol.npos) pnum = symbol.size();
	size_t p1 = 0, p2 = 0, nc = 0;
	bool bfound = false;
	while (!bfound && iz < (int)ATOMTYPE_ZNUM_MAX) {
		scmp = sAtomSymbols[iz];
		p1 = symbol.find(scmp, 0);
		if (std::string::npos != p1 && pnum == scmp.size()) { // symbol found
			bfound = true;
			break;
		}
		iz++; // symbol not found, try next
	}
	p2 = symbol.find_first_of("+-", 0);
	if (std::string::npos != p2) { // charge sign found
		nc = p2 - p1 - scmp.size();
		if (nc > 0) {
			crg = to_float(symbol.substr(p1 + scmp.size(), nc));
			if (symbol[p2] == '-') crg = -crg;
		}
	}
	return crg;
}

bool CAtomType::SameSpecies(CAtomType *pAtt)
{
	if (NULL==pAtt)
		return false;
	return (	m_nAtomicNumber==pAtt->m_nAtomicNumber && 
				fabsf(m_fIonicCharge-pAtt->m_fIonicCharge)<(float)ATOMTYPE_CHRGDIF_THR );
}


bool CAtomType::SameSpeciesBiso(CAtomType *pAtt)
{
	if (NULL==pAtt)
		return false;
	return (	m_nAtomicNumber==pAtt->m_nAtomicNumber && 
				fabsf(m_fIonicCharge-pAtt->m_fIonicCharge)<(float)ATOMTYPE_CHRGDIF_THR && 
				fabsf(m_fBiso-pAtt->m_fBiso)<(float)ATOMTYPE_BISODIF_THR );
}

std::string CAtomType::SpeciesName(void)
{
	std::string sName = "";
	std::string sSymb = sAtomSymbols[m_nAtomicNumber];
	float fchrg = fabsf(m_fIonicCharge);
	if (fchrg>=0.1f) { // name string with charge
		size_t lenstr = 4;
		std::string sSign = "+";
		if (m_fIonicCharge<0.0f)
			sSign = "-";
		std::string sChrg = format("%.2f",fchrg);
		if (sChrg.length() > lenstr) {
			sChrg = sChrg.substr(0, lenstr);
		}
		else {
			lenstr = sChrg.length();
		}
		while (lenstr > 1 && sChrg[lenstr - 1] == '0') { // rem trailing zeros
			lenstr--;
		}
		if (lenstr > 1 && sChrg[lenstr - 1] == '.') { // rem trailing point
			lenstr--;
		}
		sChrg = sChrg.substr(0, lenstr);
		sName = sSymb + sChrg + sSign;
	}
	else {
		sName = sSymb;
	}
	return sName;
}

std::string CAtomType::SpeciesComposition(void)
{
	std::string sName = SpeciesName();
	std::string sOcc = "";
	std::string sComp = "";
	size_t lenstr = 4;
	if (m_fOccupancy>=1000.0f) {
		int nOcc = (int)(m_fOccupancy+0.5f);
		sOcc = format("%d",nOcc);
	}
	else {
		sOcc = format("%f",m_fOccupancy);
		if (sOcc.length() > lenstr) {
			sOcc = sOcc.substr(0, lenstr);
		}
		else {
			lenstr = sOcc.length();
		}
		while (lenstr > 1 && sOcc[lenstr - 1] == '0') { // rem trailing zeros
			lenstr--;
		}
		if (lenstr > 1 && sOcc[lenstr - 1] == '.') { // rem trailing point
			lenstr--;
		}
		sOcc = sOcc.substr(0, lenstr);
	}
	sComp = sName + "_" + sOcc;
	return sComp;
}


int CAtomType::CalculateRadialFormFactors(size_t nk, float sk, float bwlk, float ekv, bool dwf, unsigned int absorb, float abf, unsigned int atffacs, bool bverbose)
{
	int nerr = 0;
	size_t ik = 0;
	double q = 0., sq = (double)sk, ht = (double)ekv;
	if (nk != nk_scaf || NULL == scaf) {
		if (NULL != scaf) {
			free(scaf);
			scaf = NULL;
			nk_scaf = 0;
		}
		scaf = (fcmplx*)malloc(sizeof(fcmplx) * nk);
		if (NULL == scaf) {
			nerr = 1;
			std::cerr << "Error (CAtomType::CalculateRadialFormFactors): failed to allocate memory for radial form factors of " << SpeciesName() << "." << std::endl;
			goto _exit;
		}
		nk_scaf = nk;
		memset(scaf, 0, sizeof(fcmplx) * nk);
	}
	if (bverbose) {
		std::cout << "  - calculating radial form factors of " << SpeciesName() << std::endl;
	}

	// set smapling rate for keepting track
	sk_scaf = sk;

	if (atffacs == 0) { // Weickenmeier & Kohl tables
		CWEKOScat weko;
		dcmplx dtmp(0.,0.), dtmp0(0.,0.);
		double biso = 0.;
		if (dwf) {
			biso = (double)m_fBiso;
		}
		else {
			biso = (double)bwlk;
		}
		for (ik = 0; ik < nk_scaf; ik++) { // loop over spatial frequencies
			q = sq * (double)ik; // get current spatial frequency
			dtmp = dtmp0; // preset result with zero
			if (absorb == 2) { // Hall & Hirsch -> absorptive form factors
				dtmp = weko.wekosca(q, biso, m_nAtomicNumber, ht, dwf, true);
			}
			else { // no integration, get elastic ff
				dtmp.real((float)weko.wekoscar(q, biso, m_nAtomicNumber, ht, dwf));
				if (absorb == 1) { // Hashimoto et al. -> absorptive form factor as fraction of the elastic form factor
					dtmp.imag(dtmp.real() * (double)abf);
				}
			}
			scaf[ik] = fcmplx((float)dtmp.real(), (float)dtmp.imag()); // set the value (to float)
		}
	}

	//
	// * TODO *
	//
	// * implement the use of other tables for scattering factors here
	//

_exit:
	return nerr;
}


fcmplx CAtomType::InterpolateScaf(float q, unsigned int ipo)
{
	// static data only changing if ipo changes
	static size_t n_l = 0;
	static size_t n_l2 = 0;
	static float* lq = NULL; // list of q values
	static fcmplx* ly = NULL; // list of values
	static fcmplx* lyc = NULL; // helper list
	static fcmplx* lyd = NULL; // helper list
	// variable data changing if q changes
	long ix = 0, i1 = 0, i2 = 0, i = 0, j = 0, j1 = 0;
	long ns = 0;
	float fp = 0.f, dif = 0.f, difl = 0.f, dift = 0.f, ho = 0.f, hp = 0.f, den = 0.f;
	fcmplx y(0.f, 0.f);
	fcmplx dy(0.f, 0.f);
	fcmplx w(0.f, 0.f);
	fcmplx cden(0.f, 0.f);
	if (nk_scaf == 0 || NULL == scaf) goto _exit; // skip if table not prepared
	// setup of static helper arrays for interpolation
	if (n_l < (size_t)ipo + 1) { // n_l decides whether we need a re-allocation of static helper lists
		if (NULL != lq) {
			free(lq); lq = NULL;
		}
		if (NULL != ly) {
			free(ly); ly = NULL;
		}
		n_l = (size_t)ipo + 1;
		n_l2 = (n_l - (n_l % 2)) >> 1;
		lq = (float*)malloc(sizeof(float) * n_l);
		ly = (fcmplx*)malloc(sizeof(fcmplx) * n_l);
		lyc = (fcmplx*)malloc(sizeof(fcmplx) * n_l);
		lyd = (fcmplx*)malloc(sizeof(fcmplx) * n_l);
	}
	// setup of variable helper data for interpolation
	fp = fabsf(q) / sk_scaf;
	ix = (long)std::min(nk_scaf, (size_t)(0.5 + fp));
	i1 = ix - (long)n_l2;
	i2 = i1 + (long)n_l - 1;
	for (i = i1; i <= i2; i++) {
		j1 = i - i1; // index in the helper array
		j = i;
		if (j < 0) j = 0 - i; // reversal of data at q = 0
		j = std::min((long)nk_scaf - 1, j); // limit to tabulated range
		lq[j1] = sk_scaf * (float)i;
		ly[j1] = scaf[j];
	}
	// run the interpolation
	ns = 0;
	y = ly[0];
	dif = fabsf(q - lq[0]);
	for (i = 0; i < (long)n_l; i++) { // find nearest neighbor solution
		dift = fabsf(q - lq[i]);
		if (dift < dif) {
			ns = i;
			dif = dift;
		}
		lyc[i] = ly[i];
		lyd[i] = ly[i];
	}
	y = ly[ns]; // set nearest neighbor solution
	dif = fabsf(q - lq[ns]);
	if (dif / sk_scaf < 0.00001f) goto _exit; // the point is very close to the next sample, use the sample
	ns = std::max((long)0, ns - 1);
	for (j = 0; j < (long)n_l - 1; j++) { // for each column of the table
		for (i = 0; i < (long)n_l - j - 1; i++) { // ... update c and d
			ho = lq[i] - q;
			hp = lq[i + j + 1] - q;
			w = lyc[i + 1] - lyd[i];
			den = ho - hp; // in any case den should never be 0.f of sk_scaf > 0.f
			cden = w / den;
			lyd[i] = cden * hp;
			lyc[i] = cden * ho;
		}
		//
		// Decide now, which correction must be applied to y (c or d)
		//
		if (2 * ns < (long)n_l - j) {
			dy = lyc[ns + 1];
		}
		else {
			dy = lyd[ns];
			ns--;
		}
		y = y + dy;
	}
_exit:
	return y;
}


fcmplx CAtomType::GetFormFactor(float q, unsigned int ipo)
{
	fcmplx ff(0.f, 0.f);
	size_t ix = 0, i0 = 0, i1 = 0;
	float fp = 0.f, fi = 0.f;

	if (nk_scaf == 0 || NULL == scaf) goto _exit; // skip if table not prepared

	fp = fabsf(q) / sk_scaf;
	ix = std::min(nk_scaf - 1, (size_t)( 0.5 + fp )); // nearest neighbor interpolation (ipo == 0)
	ff = scaf[ix];

	switch (ipo) { // higher order interpolations
	case 1: // linear interpolation
		i0 = (size_t)fp; i1 = i0 + 1;
		fi = fp - (float)i0;
		ff = scaf[i0] * (1.0f - fi) + scaf[i1] * fi;
		break;
	default: // polynomial interpolation
		ff = InterpolateScaf(q, ipo);
	}

_exit:
	return ff;
}

//
// ************************************************************************** //







// ************************************************************************** //
//
// IMPLEMENTATION OF CLASS CATOM
//
// * CONSTRUCTORS *********************************************************** //

CAtom::CAtom(void)
{
	// standard constructor, preset with defaults
	//
	ResetData();
}

CAtom::CAtom(const CAtom &atSrc)
{
	ResetData();
	//
	m_nAtomTypeIndex = atSrc.m_nAtomTypeIndex;
	m_vPos = atSrc.m_vPos;
	m_fOccupancy = atSrc.m_fOccupancy;
}

// * DECONSTRUCTORS ********************************************************* //

CAtom::~CAtom(void)
{
}

// * OPERATORS ************************************************************** //

void CAtom::operator=(const CAtom &other)
{
	m_nAtomTypeIndex = other.m_nAtomTypeIndex;
	m_vPos = other.m_vPos;
	m_fOccupancy = other.m_fOccupancy;
}

const bool CAtom::operator==(const CAtom &other) const
{
	bool bResult = true;

	bResult &= ( m_nAtomTypeIndex == other.m_nAtomTypeIndex );
	CV3D vDif = m_vPos - other.m_vPos;
	bResult &= ((float)ATOMTYPE_OCCUDIF_THR > fabsf(m_fOccupancy - other.m_fOccupancy));
	bResult &= ( (float)ATOM_POSDIF_THR > vDif.Length() );

	return bResult;
}

// * MEMBER FUNCTIONS ******************************************************* //

void CAtom::ResetData(void)
{
	m_nAtomTypeIndex = -1;
	m_vPos = CV3D();
	m_fOccupancy = 1.0f;
}

// ************************************************************************** //

