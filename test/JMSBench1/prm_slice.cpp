#include "prm_slice.h"
#include "NatureConstants.h"


prm_slice::prm_slice()
{
	pro_nx = 0;
	pro_ny = 0;
	sz_pgr = 0;
	sz_pro = 0;
	pgr = NULL; // pointer to the memory block holding phase grating data (variants in sequence)
	pro = NULL; // pointer to the memory block holding propagator data
}

prm_slice::prm_slice(const prm_slice &src) : params(src)
{
	pro_nx = 0;
	pro_ny = 0;
	sz_pgr = 0;
	sz_pro = 0;
	pgr = NULL;
	pro = NULL;

	prm_slice *psrc = const_cast <prm_slice*>(&src);
	header = src.header;
	st = src.st;
	if (NULL != psrc->get_pgr()) {
		alloc_pgr();
		memcpy((void*)pgr, (void*)psrc->get_pgr(), sz_pgr);
	}
	if (NULL != psrc->get_pro()) {
		size_t nx = 0, ny = 0;
		psrc->get_pro_dim(nx, ny);
		alloc_pro(nx, ny);
		memcpy((void*)pro, (void*)psrc->get_pro(), sz_pro);
	}
}

prm_slice::~prm_slice()
{
	clear();
}

void prm_slice::clear()
{
	st.DeleteStructureData();
	if (NULL != pro) {
		free(pro);
		pro = NULL;
		sz_pro = 0;
	}
	if (NULL != pgr) {
		free(pgr);
		pgr = NULL;
		sz_pgr = 0;
	}
}

int prm_slice::load_ems_header(std::string sfile)
{
	int nerr = 0, ierr = 0;
	std::ifstream fs;

	if (!file_exists(sfile)) {
		std::cerr << "Error (prm_slice::load_ems_data): File " << sfile << " not found." << std::endl;
		nerr = 1;
		goto _exit;
	}

	// open the file
	fs.open(sfile, std::ios::binary);
	if (!fs.is_open()) {
		std::cerr << "Error (prm_slice::load_ems_data): Failed to open file " << sfile << "." << std::endl;
		nerr = 2;
		goto _exit;
	}

	// load the header
	ierr = header.load_ems_header(&fs);
	if (0 < ierr) {
		std::cerr << "Error (prm_slice::load_ems_data): Failed to open file " << sfile << "." << std::endl;
		nerr = 3;
		goto _exit;
	}
	header.str_file_name = sfile;

	ierr = load_ems_itab(&fs, header.ver);
	if (0 < ierr) {
		std::cerr << "Error (prm_slice::load_ems_data): Failed read structure data from file " << sfile << "." << std::endl;
		nerr = 4;
		goto _exit;
	}

_exit:
	if (fs.is_open()) fs.close(); // close file
	return nerr;
}

int prm_slice::load_ems_itab(std::ifstream* pfs, unsigned int header_version)
{
	int nerr = 0;

	unsigned __int32 natt = 0;
	unsigned __int32 natp_max = 0, natp = 0;
	unsigned __int32 natr_max = 0, natr = 0;
	unsigned __int32 iatt = 0, iat = 0;
	__int32 itmp = 0;
	float ftmp = 0.f;
	size_t szi4 = sizeof(__int32);
	size_t szf4 = sizeof(float);
	CAtomType att;
	CAtom at;
	bool readok = true;

	if (NULL == pfs) {
		std::cerr << "Error (prm_slice::load_ems_itab): Invalid parameter." << std::endl;
		nerr = 1;
		goto _exit;
	}

	if (!pfs->is_open()) {
		std::cerr << "Error (prm_slice::load_ems_itab): Input file stream is not opened." << std::endl;
		nerr = 2;
		goto _exit;
	}

	if (!pfs->good()) {
		std::cerr << "Error (prm_slice::load_ems_itab): Input file stream is not good for reading." << std::endl;
		nerr = 3;
		goto _exit;
	}

	if (header.ver != header_version) {
		std::cerr << "Error (prm_slice::load_ems_itab): Version parameter (" << header_version  << ") is inconsistent with current header data (" << header.ver  << ")." << std::endl;
		nerr = 4;
		goto _exit;
	}

	st.m_vTypes.clear();
	st.m_vAtoms.clear();
	st.m_vSites.clear();

	pfs->seekg((std::streampos)header.structure_offset);
	if (readok) readok &= file_read_buf(pfs, (char*)&natt, szi4, 1);
	if (readok) readok &= file_read_buf(pfs, (char*)&natp_max, szi4, 1);
	if (readok) readok &= file_read_buf(pfs, (char*)&natr_max, szi4, 1);

	if (!readok) {
		std::cerr << "Error (prm_slice::load_ems_itab): Failed reading table header." << std::endl;
		nerr = 5;
		goto _exit;
	}

	switch (header_version) {

	case 2012071101: // old slice data version
		if (natt > 0) { // there is atom type data
			for (iatt = 0; iatt < natt; iatt++) { // loop iatt over atom types
				att.ResetData();
				if (readok) readok &= file_read_buf(pfs, (char*)&att.m_nAtomicNumber, szi4, 1);
				if (readok) readok &= file_read_buf(pfs, (char*)&att.m_fIonicCharge, szf4, 1);
				if (readok) readok &= file_read_buf(pfs, (char*)&att.m_fBiso, szf4, 1);
				if (readok) readok &= file_read_buf(pfs, (char*)&att.m_fOccupancy, szf4, 1);
				if (readok) readok &= file_read_buf(pfs, (char*)&natp, szi4, 1); // number of atom data following
				if (!readok) {
					std::cerr << "Error (prm_slice::load_ems_itab): Failed reading data of atom type " << iatt << " of " << iatt << "." << std::endl;
					nerr = 5;
					goto _exit;
				}
				if (natp > 0) { // there are atom positions to read
					for (iat = 0; iat < natp; iat++) { // loop iat over atom positions of the current type
						at.ResetData();
						if (readok) readok &= file_read_buf(pfs, (char*)&at.m_vPos.x, szf4, 1);
						if (readok) readok &= file_read_buf(pfs, (char*)&at.m_vPos.y, szf4, 1);
						if (readok) readok &= file_read_buf(pfs, (char*)&at.m_vPos.z, szf4, 1);
						if (!readok) {
							std::cerr << "Error (prm_slice::load_ems_itab): Failed reading data of atom " << iat << " of " << natp << " of type " << iatt << " of " << iatt << "." << std::endl;
							nerr = 6;
							goto _exit;
						}
						// add the atom to the structure
						st.AddAtom(att.m_nAtomicNumber, att.m_fIonicCharge, att.m_fBiso, att.m_fOccupancy, at.m_vPos.x, at.m_vPos.y, at.m_vPos.z);
					} // loop iat over atom positions of the current type
				}
				if (readok) readok &= file_read_buf(pfs, (char*)&natr, szi4, 1); // number of transition data following
				if (natr > 0) { // there are transition codes to load
					for (iat = 0; iat < natr; iat++) { // loop iat over atom transitions of the current type
						if (readok) readok &= file_read_buf(pfs, (char*)&itmp, szi4, 1); // transition code
						if (readok) readok &= file_read_buf(pfs, (char*)&ftmp, szf4, 1); // transition energy
					}
				}
			} // loop iatt over atom types
		}
		st.SetAtomSites(0.05f); // fill atom site list
		break;

	case 2020052801: // new slice data version
		if (natt > 0) { // there is atom type data
			for (iatt = 0; iatt < natt; iatt++) { // loop iatt over atom types
				att.ResetData();
				if (readok) readok &= file_read_buf(pfs, (char*)&att.m_nAtomicNumber, szi4, 1);
				if (readok) readok &= file_read_buf(pfs, (char*)&att.m_fIonicCharge, szf4, 1);
				if (readok) readok &= file_read_buf(pfs, (char*)&att.m_fBiso, szf4, 1);
				if (readok) readok &= file_read_buf(pfs, (char*)&ftmp, szf4, 1); // do not read the type Occupancy, we build this up by adding the atoms to the structure
				if (readok) readok &= file_read_buf(pfs, (char*)&natp, szi4, 1); // number of atom data following this type
				if (!readok) {
					std::cerr << "Error (prm_slice::load_ems_itab): Failed reading data of atom type " << iatt << " of " << iatt << "." << std::endl;
					nerr = 5;
					goto _exit;
				}
				if (natp > 0) { // there are atom positions to read
					for (iat = 0; iat < natp; iat++) { // loop iat over atom positions of the current type
						at.ResetData();
						if (readok) readok &= file_read_buf(pfs, (char*)&at.m_vPos.x, szf4, 1);
						if (readok) readok &= file_read_buf(pfs, (char*)&at.m_vPos.y, szf4, 1);
						if (readok) readok &= file_read_buf(pfs, (char*)&at.m_vPos.z, szf4, 1);
						if (readok) readok &= file_read_buf(pfs, (char*)&at.m_fOccupancy, szf4, 1);
						if (!readok) {
							std::cerr << "Error (prm_slice::load_ems_itab): Failed reading data of atom " << iat << " of " << natp << " of type " << iatt << " of " << iatt << "." << std::endl;
							nerr = 6;
							goto _exit;
						}
						// add the atom to the structure
						st.AddAtom(att.m_nAtomicNumber, att.m_fIonicCharge, att.m_fBiso, at.m_fOccupancy, at.m_vPos.x, at.m_vPos.y, at.m_vPos.z);
					} // loop iat over atom positions of the current type
				}
				if (readok) readok &= file_read_buf(pfs, (char*)&natr, szi4, 1); // number of transition data following
				if (natr > 0) { // there are transition codes to load
					for (iat = 0; iat < natr; iat++) { // loop iat over atom transitions of the current type
						if (readok) readok &= file_read_buf(pfs, (char*)&itmp, szi4, 1); // transition code
						if (readok) readok &= file_read_buf(pfs, (char*)&ftmp, szf4, 1); // transition energy
					}
				}
			} // loop iatt over atom types
		}
		st.SetAtomSites(0.05f); // fill atom site list
		break;
	default: // unknown version
		std::cerr << "Error (prm_slice::load_ems_itab): Unsupported header version (" << header_version << "). Atom table not loaded." << std::endl;
		nerr = 10;
		goto _exit;
	}

	if (btalk) {
		std::cout << "  - loaded slice structure with " << st.m_vTypes.size() << " atom types and " << st.m_vAtoms.size() << " atoms." << std::endl;
		std::cout << "    slice structure composition: " << st.GetCompositionString() << std::endl;
	}

_exit:
	return nerr;
}

int prm_slice::load_ems_data(std::string sfile)
{
	int ierr = 0;
	std::ifstream fs;
	size_t num_floats = (size_t)2 * header.grid_x * header.grid_y * header.num_var;

	if (num_floats == 0 || header.data_offset < (unsigned long long)EMS_POS_DATA) {
		std::cerr << "Error (prm_slice::load_ems_data): Failed to load ems data due to invalid header presets." << std::endl;
		ierr = 1;
		goto _exit;
	}

	if (!file_exists(sfile)) {
		std::cerr << "Error (prm_slice::load_ems_data): Failed to load ems data, file " << sfile << " not found." << std::endl;
		ierr = 2;
		goto _exit;
	}

	if (btalk) {
		std::cout << "  - loading data from file " << sfile << "." << std::endl;
	}

	// allocate destination buffer
	if (0 == alloc_pgr()) {
		ierr = 3;
		goto _exit;
	}

	// open the file
	fs.open(sfile, std::ios::binary);
	if (fs.is_open()) {
		// goto data offset
		fs.seekg((std::streampos)header.data_offset);
		// load data from file
		if (!file_read_buf(&fs, (char*)pgr, 1, sz_pgr)) {
			std::cerr << "Error (prm_slice::load_ems_data): Failed to load ems data, file " << sfile << " not found." << std::endl;
			ierr = 5;
			goto _exit;
		}
	}
	else {
		std::cerr << "Error (prm_slice::load_ems_data): Failed to open file " << sfile << "." << std::endl;
		ierr = 4;
		goto _exit;
	}

_exit:
	if (fs.is_open()) fs.close(); // close file
	return ierr;
}

int prm_slice::save_ems_itab(std::ofstream* pfs)
{
	int nerr = 0, ierr = 0;
	size_t sz_itab = header.data_offset - header.structure_offset;
	size_t buf_pos = 0; // buffer output position
	unsigned __int32 natt = (unsigned __int32)st.m_vTypes.size(); // number of types
	unsigned __int32 nat = (unsigned __int32)st.m_vAtoms.size(); // number of atoms
	unsigned __int32 natp = 0; // number of atoms of a type
	unsigned __int32 natr = 0; // number of transistions of a type
	unsigned __int32 natp_max = (unsigned __int32)st.GetMaxAtomNumber(); // max. number of atoms per type
	unsigned __int32 natr_max = 0; // max. number of inealstic transitions per type
	unsigned __int32 iatt = 0, iat = 0;
	char* buf = NULL;

	if (NULL == pfs) {
		nerr = 1;
		std::cerr << "Error (prm_slice::save_ems_itab): invalid parameter." << std::endl;
		goto _exit;
	}

	if (!pfs->is_open()) {
		nerr = 2;
		std::cerr << "Error (prm_slice::save_ems_itab): output file stream is not open." << std::endl;
		goto _exit;
	}

	buf = new char[sz_itab];
	if (NULL == buf) {
		nerr = 4;
		std::cerr << "Error (prm_slice::save_ems_itab): output buffer allocation failed." << std::endl;
		goto _exit;
	}

	memcpy(&buf[buf_pos], &natt, sizeof(unsigned __int32)); buf_pos += sizeof(unsigned __int32); // number of atom types
	memcpy(&buf[buf_pos], &natp_max, sizeof(unsigned __int32)); buf_pos += sizeof(unsigned __int32); // max. number of atoms per types
	memcpy(&buf[buf_pos], &natr_max, sizeof(unsigned __int32)); buf_pos += sizeof(unsigned __int32); // max. number of transitions per types

	if (natt > 0) { // there are atom types to store
		for (iatt = 0; iatt < natt; iatt++) { // loop over atom types
			natp = (unsigned __int32)st.GetAtomNumber((int)iatt); // number of atoms of this type
			memcpy(&buf[buf_pos], &st.m_vTypes[iatt].m_nAtomicNumber, sizeof(int)); buf_pos += sizeof(int); // atomic number
			memcpy(&buf[buf_pos], &st.m_vTypes[iatt].m_fIonicCharge, sizeof(float)); buf_pos += sizeof(float); // ionic charge
			memcpy(&buf[buf_pos], &st.m_vTypes[iatt].m_fBiso, sizeof(float)); buf_pos += sizeof(float); // Biso
			memcpy(&buf[buf_pos], &st.m_vTypes[iatt].m_fOccupancy, sizeof(float)); buf_pos += sizeof(float); // type occupancy
			memcpy(&buf[buf_pos], &natp, sizeof(unsigned __int32)); buf_pos += sizeof(unsigned __int32); // number of positions with the current atom type
			if (natp > 0) { // store related atom positions
				for (iat = 0; iat < nat; iat++) { // loop over all atoms
					if (st.m_vAtoms[iat].m_nAtomTypeIndex == (int)iatt) { // this atom is of the current type
						memcpy(&buf[buf_pos], &st.m_vAtoms[iat].m_vPos.x, sizeof(float)); buf_pos += sizeof(float); // fractional x position
						memcpy(&buf[buf_pos], &st.m_vAtoms[iat].m_vPos.y, sizeof(float)); buf_pos += sizeof(float); // fractional y position
						memcpy(&buf[buf_pos], &st.m_vAtoms[iat].m_vPos.z, sizeof(float)); buf_pos += sizeof(float); // fractional z position
						memcpy(&buf[buf_pos], &st.m_vAtoms[iat].m_fOccupancy, sizeof(float)); buf_pos += sizeof(float); // site occupancy
					}
				} // loop iat over atoms
			}
			memcpy(&buf[buf_pos], &natr, sizeof(unsigned __int32)); buf_pos += sizeof(unsigned __int32); // number of transitions with the current atom type
			if (natr > 0) { // store related transition code
				// 2020-05-28 - This is not implemented.
			}
		} // loop iatt over atom types
	}

	if (!file_write_buf(pfs, buf, sz_itab, 1)) {
		nerr = 100;
		std::cerr << "Error (prm_slice::save_ems_itab): failed to write atom table to output file stream (size: " << (sz_itab >> 10) << " kB)." << std::endl;
		goto _exit;
	}

_exit:
	if (NULL != buf) { delete[] buf; }
	return nerr;
}

int prm_slice::save_ems(std::string sfile)
{
	int nerr = 0, ierr = 0;
	int itmp = 0;
	float ftmp = 0.f;
	size_t num_floats = 0;
	std::ofstream fs;
	
	// check reasonable slice data
	if (NULL == pgr) {
		nerr = 1;
		std::cerr << "Error (prm_slice::save_ems): No data to save." << std::endl;
		goto _exit;
	}

	num_floats = 

	header.structure_offset = get_file_structure_offset();
	header.data_offset = get_file_data_offset((unsigned int)(SLICEPARAMS_EXTHDRVER));
	header.data_type = (unsigned int)(SLICEPARAMS_IFORM);

	fs.open(sfile, std::ios::binary | std::ios::trunc);

	if (fs.is_open()) {
		// basic header
		ierr = header.save_ems_header(&fs);
		if (0 < ierr) {
			nerr = 11;
			std::cerr << "Error (prm_slice::save_ems): Failed to write slice file header (" << ierr << ")." << std::endl;
			goto _exit;
		}
		// atom table
		ierr = save_ems_itab(&fs);
		if (0 < ierr) {
			nerr = 12;
			std::cerr << "Error (prm_slice::save_ems): Failed to write slice structure table (" << ierr << ")." << std::endl;
			goto _exit;
		}
		// object functions
		if (!file_write_buf(&fs, (char*)pgr, 1, sz_pgr)) {
			std::cerr << "Error (prm_slice::save_ems): Failed to save transmission functions (size: " << (sz_pgr >> 20) << " MB)." << std::endl;
			nerr = 13;
			goto _exit;
		}
	}
	else {
		std::cerr << "Error (prm_slice::save_ems): Failed to open file " << sfile << "." << std::endl;
		nerr = 2;
		goto _exit;
	}

_exit:
	if (fs.is_open()) fs.close(); // close file
	if (nerr == 0 && btalk) {
		std::cout << "  - saved object transmission functions to file " << sfile << "." << std::endl;
	}
	return nerr;
}

size_t prm_slice::alloc_pgr()
{
	sz_pgr = sizeof(fcmplx) * header.grid_x * header.grid_y * header.num_var;
	if (NULL != pgr) {
		free(pgr);
		pgr = NULL;
	}
	if (sz_pgr > 0) {
		pgr = (float*)malloc(sz_pgr);
		if (NULL == pgr) {
			std::cerr << "Error (prm_slice::alloc_pgr): failed to allocate phase gratings memory." << std::endl;
			std::cerr << "- requested size [MB]: " << (sz_pgr >> 20) << std::endl;
			sz_pgr = 0;
			goto _exit;
		}
		memset(pgr, 0, sz_pgr); // preset with 0
	}
_exit:
	return sz_pgr;
}

size_t prm_slice::alloc_pro(size_t nx, size_t ny)
{
	if (NULL != pro) {
		free(pro);
		pro = NULL;
		sz_pro = 0;
		pro_nx = 0;
		pro_ny = 0;
	}
	sz_pro = sizeof(fcmplx) * nx * ny;
	if (sz_pro > 0) {
		pro = (float*)malloc(sz_pro);
		if (NULL == pro) {
			std::cerr << "Error (prm_slice::alloc_pro): failed to allocate propagator memory." << std::endl;
			std::cerr << "- requested size [MB]: " << sz_pro / 1048576 << std::endl;
			sz_pro = 0;
			goto _exit;
		}
		memset(pro, 0, sz_pro); // preset with 0
		pro_nx = nx;
		pro_ny = ny;
	}
_exit:
	return sz_pro;
}

void prm_slice::get_pro_dim(size_t &nx, size_t &ny)
{
	nx = pro_nx;
	ny = pro_ny;
}

fcmplx* prm_slice::get_pgr(size_t var)
{
	size_t npixvar = (size_t)header.grid_x * header.grid_y;
	return &(((fcmplx*)pgr)[var * npixvar]);
}

int prm_slice::calculate_propagator(CJMultiSlice *pjms, float tilt_x, float tilt_y, int type)
{
	int nerr = 0, ierr = 0;
	int nx = 0, ny = 0;
	float ttf = (float)_R2D * 0.001f; // scales from mrad to degrees
	float otx = ttf * tilt_x, oty = ttf * tilt_y; // tilts in degree

	if (NULL == pjms) {
		nerr = 1;
		std::cerr << "Error (prm_slice::calculate_propagator) called with invalid interface to JMS." << std::endl;
		goto _exit;
	}

	pjms->GetGridSize(nx, ny);
	alloc_pro((size_t)nx, (size_t)ny);
	if (sz_pro == 0) {
		nerr = 3;
		std::cerr << "Error (prm_slice::calculate_propagator) failed to allocate memory." << std::endl;
		goto _exit;
	}

	ierr = pjms->CalculatePropagator(header.sz, otx, oty, (fcmplx*)pro, type);
	if (0 < ierr) {
		nerr = 4;
		std::cerr << "Error (prm_slice::calculate_propagator) failed to calculate propagator." << std::endl;
		goto _exit;
	}

_exit:
	return nerr;
}

fcmplx* prm_slice::get_pro(void)
{
	return (fcmplx*)pro;
}

fcmplx prm_slice::get_pgr_at(size_t pos, size_t var)
{
	size_t npixvar = (size_t)header.grid_x * header.grid_y;
	return ((fcmplx*)pgr)[pos + var * npixvar];
}

fcmplx prm_slice::get_pro_at(size_t pos)
{
	return ((fcmplx*)pro)[pos];
}

void prm_slice::get_grid_dim(unsigned int* dim_x, unsigned int* dim_y)
{
	*dim_x = header.grid_x;
	*dim_y = header.grid_y;
}

void prm_slice::get_grid_size(float* size_x, float* size_y)
{
	*size_x = header.sx;
	*size_y = header.sy;
}


float prm_slice::get_energy_kev(void)
{
	return header.ekv;
}

float prm_slice::get_thickness(void)
{
	return header.sz;
}

unsigned int prm_slice::get_variant_num(void)
{
	return header.num_var;
}

unsigned long long prm_slice::get_file_data_offset(unsigned int header_version)
{
	unsigned long long offset = (unsigned long long)(EMS_POS_DATA);
	unsigned long long offset2 = get_file_structure_offset() + get_file_structure_size(header_version) + 16;
	offset = std::max(offset, offset2);
	return offset;
}

unsigned long long prm_slice::get_file_structure_offset(void)
{
	return (unsigned long long)(EMS_POS_ITAB);
}

unsigned long long prm_slice::get_file_structure_size(unsigned int header_version)
{
	unsigned long long itab_size = 0;
	unsigned long long natt = st.m_vTypes.size(); // number of atom types
	unsigned long long nadat = 3; // number of data per atom type (charge, biso, occupancy)
	unsigned long long niatr_max = 0; // max. number of inelastic transition data per atom type (placeholder = 0, not used)
	unsigned long long natp_max = 0; // max. number of position data per atom type
	unsigned long long nat = st.m_vAtoms.size(); // number of atoms in the structure
	if (natt > 0 && nat > 0) {
		natp_max = st.GetMaxAtomNumber();
		switch (header_version) {
		case 2012071101: // older code
			itab_size = 4 * (3 + natt * (3 + nadat + 3 * natp_max + 2 * niatr_max));
			break;
		case 2020052801: // saved position with occupancy, type occupancy is the sum over related position occupancies
			itab_size = 4 * (3 + natt * (3 + nadat + 4 * natp_max + 2 * niatr_max));
			break;
		default:
			itab_size = 4 * 3;
		}
	}
	else {
		itab_size = 4 * 3;
	}	
	return itab_size;
}