#pragma once
#include "params.h"
// #include "prm_structure.h"

#define SLICENAME_MAX				40
#define SLICEPARAMS_EXTHDRVER		2012071101 // 2010112301
#define SLICEPARAMS_IFORM			1 // = COMPLEX REAL SPACE DATA
#define ENDIANTEST_MAX				0x00002000

#ifndef EMSPARAMS
#define EMSPARAMS
#define EMS_POS_DATASIZE_X			0x00000000
#define EMS_POS_DATASIZE_Y			0x00000004
#define EMS_POS_IFORM				0x00000010
#define EMS_POS_TITLE				0x0000004C
#define EMS_POS_EXTHEADER_VER		0x00000074
#define EMS_POS_VARIANT_NUM			0x00000078
#define EMS_POS_CONTENT_TYPE		0x00000080
#define EMS_POS_SUBSLICE_NUM		0x00000084
#define EMS_POS_ALT_OFFSET			0x00000088
#define EMS_POS_ITAB_OFFSET			0x0000008c
#define EMS_POS_THICKNESS			0x00001034
#define EMS_POS_HIGHTENSION			0x00001040
#define EMS_POS_PHYSSIZE_X			0x00001044
#define EMS_POS_PHYSSIZE_Y			0x00001048
#define EMS_POS_ITAB				0x00001100
#define EMS_POS_DATA				0x00002000
#endif // EMSPARAMS

/*
!**********************************************************************!
!
! Information on the EMS image and slice files
! ----------------------------------------
!
! Structure information was interpreted from the fortran source code
! of EMS by P. Stadelmann, files 'openfi.f' and 'openfi.com'.
!
! The file consists of a header block and a data block.
! The header block is 8192 bytes long.
! The data block consists of nsam*nrow data items, where the
! item size depends on the data type.
!
! Header structure:
!    OFFSET(DEC,HEX) TYPE            SIZE        INFO
!    0      0000     integer*4       4           nsam, image dims, number of samples
!    4      0004     integer*4       4           nrow, image dims, number of rows
!    8      0008     integer*4       4           nsaeff, effective x dim of the image
!    12     000C     integer*4       4           nroeff, effective y dim of the image
!    16     0010     integer*4       4           iform, data format, -1=='F' FT=complex_64, 0=='R' real_32, 1=='C' complex_64, 2=='I' int_8
!    20     0014     integer*4       4           iprot, ? flag for title interpretation ?
!    24     0018     integer*4       4           imami, flags that min. and max. values have been calculated already (0=no, 1=yes)
!    28     001C     integer*4       4           igeom, flags that geometry correction was done (0=no, 1=yes)
!    32     0020     integer*4       4           nslice, number of multislice iteration
!    36     0024     integer*4       4           iu, first zone index
!    40     0028     integer*4       4           iv, second zone index
!    44     002C     integer*4       4           iw, third zone index
!    48     0030     integer*4       4           ispcel, super-cell flag (0: unit cell, 1: super-cell)
!    52     0034     integer*4       24          itr(3,2), 2-d to 3-d indices transform matrix
!    56     0038     integer*4       ||          itr(3,2) (2,1)
!    60     003C     integer*4       ||          itr(3,2) (3,1)
!    64     0040     integer*4       ||          itr(3,2) (1,2)
!    68     0044     integer*4       ||          itr(3,2) (2,2)
!    72     0048     integer*4       ||          itr(3,2) (3,2)
!    76     004C     character*40    40          imatit, title string of image / slice
!    ...    ...      ...             ||
!*** 116    0074     integer*4       4           extended header version (added by J.B. 101202)
!*** 120	0078     integer*4       4           number of data variants
!*** 128    0080     integer*4       4           content type (0=phase grating, 1=potential) (added by J.B. 110407)
!*** 132    0084     integer*4       4           number of sub-slices (>1 if sub-slicing is done) (added by J.B. 110407)
!*** 136    0088     integer*4       4           alternative data offset, overrides default data offset (added by J.B. 120711)
!*** 140    008C     integer*4       4           offset of the inelastic data table [see below for table structure definition] (added by J.B. 120711)
!    ...
!    4096   1000     real*4          4           beamh, 1st beam trace index (3-d)
!    4100   1004     real*4          4           beamk, 2nd beam trace index (3-d)
!    4104   1008     real*4          4           beaml, 3rd beam trace index (3-d)
!    4108   100C     real*4          4           picmin, min. of the image (real part)
!    4112   1010     real*4          4           picmax, max. of the image (real part)
!    4116   1014     real*4          4           commin, min. of the image (imag part)
!    4120   1018     real*4          4           commax, max. of the image (imag part)
!    4124   101C     real*4          4           a, 1st lattice side [nm]
!    4128   1020     real*4          4           b, 2nd lattice side [nm]
!    4132   1024     real*4          4           c, 3rd lattice side [nm]
!    4136   1028     real*4          4           alfa, (b,c) lattice angle [deg]
!    4140   102C     real*4          4           beta, (a,c) lattice angle [deg]
!    4144   1030     real*4          4           gama, (a,b) lattice angle [deg]
!    4148   1034     real*4          4           thick, slice thickness [nm]
!    4152   1038     real*4          4           absor, absorption coefficiant [a.u.]
!    4156   103C     real*4          4           sigma, interaction constant
!    4160   1040     real*4          4           voltag, accelerating voltage  [kv]
!    4164   1044     real*4          4           ca, image geometry in real*4 space [a.u.]
!    4168   1048     real*4          4           cb, image geometry in real*4 space [a.u.]
!    4172   104C     real*4          4           angle, image geometry in real*4 space [a.u.]
!    4176   1050     real*4          4           b2dh, 1st beam trace index (2-d)
!    4180   1054     real*4          4           b2dk, 2nd beam trace index (2-d)
!    4184   1056
!    4188   105C
!    4192   1060     real*4          4           sub-slice thickness [nm] (added by J.B. 110407)
!
!    !!!! Atom & Inelastic potential table (max size to default data offset:
!*** 4352   1100     integer*4       4           naty, number of atom types (added by J.B.) 120711)
!*** 4360   1108     integer*4       4           natpm, max. number of atom positions (for pre-allocation of arrays)
!*** 4356   1104     integer*4       4           natrm, max. number of atom transitions (for pre-allocation of arrays)
!*** 4364   110C     integer*4       4    i-+    attyz(i), atomic number of the current atom type
!*** 4368   1110     3*real*4        12     |    atdat(1:3,i), data of the current atom type
!*** 4380   1118     integer*4       4      |    niatp(i), max. number of atom positions for each atom type in the slice
!*** 4384   111C     3*real*4        12     +-k  iaatpo(1:3,k,i), fractional atom coordinates x, y, z [ca, cb, cz]
!*** 4384+12*niatp   integer*4       4      |    niatr(i), number of atom transitions for each atom type [= number of potentials for the atom type]
!*** 4388+12*niatp   integer*4       4      +-j  iatrcd(j,i), inelastic transistion codes [= 100*shell index + l-index, K -> 100, L2,3-> 223, M4,5 -> 345]
!*** 4388+4*niatr+12*niatp  real*4   4      +-j  iatren(j,i), energy losses [eV]
!*** 4392+8*niatr+12*niatp ... next attyz(i+1)
!
!    TABLE STRUCTURE IN CLEAR WRITING:
!
!    1) number of atom types
!    2) max. number of positions
!    3) max. number of transitions
!    4) LOOP atom types
!    4.1) atomic number Z
!    4.2) atom data (ionic charge, Biso, Occupancy)
!    4.3) number of atom positions
!    4.4) LOOP atom positions
!    4.4.1) atom position (ca, cb, cdz)
!    4.5) number of atom state transitions
!    4.5.1) atom transition code
!    4.5.2) atom transition energy
!
!
! Data offset (default, may change):
!    8192   2000
!
!**********************************************************************!
*/

class prm_slice :
	public params
{
public:
	prm_slice();
	~prm_slice();

protected:
	float ekv; // slice electron energy (keV)
	float sx; // size of the slice along x (nm)
	float sy; // size of the slice along y (nm)
	float sz; // slice thickness (nm)
	unsigned int grid_x; // number of grid points along x
	unsigned int grid_y; // number of grid points along y
	unsigned int data_type; // input data type
	unsigned int var_num; // number of variants
	unsigned long long data_offset; // slice data offset byte
	unsigned long long structure_offset; // structure data offset byte
	size_t sz_slc; // size of one slice data set in bytes
	size_t sz_all; // size of all slice data sets

	char s_slc_name[SLICENAME_MAX]; // name for the Slice

	// prm_structure structure; // structure data object for this slice

	float* pdata; // pointer to the memory block holding all slice data (variants in sequence)

public:

	// clear all data and free memory allocations
	void clear();

	// load ems slice file header and store information in members
	int load_ems_header(std::string sfile);

	// data interface

	void get_grid_dim(unsigned int* dim_x, unsigned int* dim_y);

	void get_grid_size(float* size_x, float* size_y);

	float get_energy_kev(void);

	float get_thickness(void);
	
	unsigned int get_variant_num(void);

	unsigned long long get_file_data_offset(void);

protected:

	// checks whether endian needs to be swapped on input
	bool check_ems_endian(std::ifstream* pfs);
};

