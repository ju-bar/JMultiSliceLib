========================================================================
	STATIC LIBRARY : JMultiSliceLib Project Overview
========================================================================

	Copyright (C) 2018 - Juri Barthel
	Copyright (C) 2018 - Forschungszentrum Juelich GmbH, Germany
	Email: juribarthel@gmail.com

	Last Version: 0.1.5 - 2018-11-13

========================================================================


Content:
	1. Purpose
	2. Structure
	3. Dependencies
	4. Statement of license
	X. Dependencies statements of licenses
	X.1 FFTW statement of license
	X.2 CUDA statement of license



1) Purpose:

JMultiSliceLib implements routines to perform electron diffraction
calculations with the multislice algorithm. The calculations can run
on multiple CPU threads and/or a single GPU thread.



2) Structure:

Routine interfaces are declared in the header file JMultiSliceLib.h
You should include this header in your code when using the library.
Interfaces are available to 
* setup and check code parameters
* prepare library internal memory management
* setup input data
* perform multislice runs
* check setup status and calculation status
* retrieve calculation results
* clean up library internal memory



3) Dependencies:

The code of JMultiSliceLib links to

* libfftw3f-3.lib (see FFTW statement of license X.1)
  JMultiSliceLib uses data output by FFTW and is in no form based on
  work represented by FFTW. Source code and library binary code of
  FFTW are available from 
  http://www.fftw.org/ (accessed April 2018)

* cudart_static.lib and cufft.lib (see CUDA statement of license X.2)
  JMultiSliceLib uses data output and functionality probided by these
  two libraries from the CUDA Toolkit 9.0 as redistributable software.
  The CUDA Toolkit is available from
  https://developer.nvidia.com/cuda-toolkit (accessed April 2018)
  In order to be able to use the GPU routnines of JMultiSliceLib, keep
  your CUDA device drivers updated.

You should link these libraries in the calling code and take care to
comply with the FFTW and CUDA software licenses.


4) Statement of license

	1)	JMultiSliceLib is provided as freeware, but only for private, 
		educational, academic, non-commercial use (that means at home
		or in educational or academic institutions).
	1a)	JMultiSliceLib is free for educational use (schools,
		universities, museums and libraries) and for use in charity
		or humanitarian organisations.
	1b)	If you intend to use JMultiSliceLib at your place of business
		or for commercial purposes, please register and purchase it.
		Commercial users: please contact the author by E-Mail for 
		prices, discounts and payment methods.
	2)	JMultiSliceLib is owned by the RWTH Aachen University and is 
		protected by copyright laws and international treaty provisions.
		Therefore, you must treat this Software library like any other 
		copyrighted material.
	3)	You may not distribute, rent, sub-license or otherwise make
		available to others the Software or documentation or copies
		thereof, except as expressly permitted in this License without
		prior written consent from the RWTH Aachen University. In the
		case of an authorized transfer, the transferee must agree to be
		bound by the terms and conditions of this License Agreement.
	4)	You may not remove any proprietary notices, labels, trademarks
		on the Software or documentation. You may not modify, de-compile,
		disassemble or reverse engineer the Software.
	5)	Limited warranty: JMultiSliceLib, related programs in the
		package, and documentation are "as is" without any warranty
		as to their performance, merchantability or fitness for any
		particular purpose. The licensee assumes the entire risk as to
		the quality and performance of the software. In no event shall
		the RWTH Aachen University or anyone else who has been involved
		in the creation, development, production, or delivery of this
		software be liable for any direct, incidental or consequential
		damages, such as, but not limited to, loss of anticipated
		profits, benefits, use, or data resulting from the use of this
		software, or arising out of any breach of warranty.
	6)	Dependency licenses remain active and must be treated separately
		to this license.



X) Dependencies statement of licenses

X.1) FFTW statement of license

	Source: http://www.fftw.org/doc/License-and-Copyright.html

	FFTW is Copyright © 2003, 2007-11 Matteo Frigo,
	Copyright © 2003, 2007-11 Massachusetts Institute of Technology.
	FFTW is free software; you can redistribute it and/or modify it
	under the terms of the GNU General Public License as published
	by the Free	Software Foundation; either version 2 of the License,
	or (at your	option) any later version.
	This program is distributed in the hope that it will be useful,
	but	WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
	General Public License for more details.
	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
	MA 02110-1301 USA
	You can also find the GPL on the GNU web site:
	http://www.gnu.org/copyleft/gpl.html
	In addition, we	kindly ask you to acknowledge FFTW and its authors
	in any program or publication in which you use FFTW. (You are not
	required to do so; it is up to your common sense to decide whether
	you want to	comply with this request or not.) For general
	publications, we suggest referencing:
	Matteo Frigo and Steven G. Johnson, “The design and implementation
	of FFTW3,” Proc. IEEE 93 (2), 216–231 (2005).
	Non-free versions of FFTW are available under terms different from
	those of the General Public License. (e.g. they do not require you
	to accompany any object code using FFTW with the corresponding
	source code.) For these alternative terms you must purchase a
	license from MIT’s Technology Licensing Office. Users interested
	in such a license should contact us (fftw@fftw.org) for more
	information.

	A version of the GPL license is added as "gpl-3.0.txt" to the
	source tree of JMultiSliceLib.
	
X.2) CUDA statement of license

	Source: https://docs.nvidia.com/cuda/eula/index.html 

	CUDA is used according to the NVIDIA Software License Agreement
	linked above and made available with JMultiSliceLib as
	"CUDA EULA.pdf". The two CUDA libraries used by JMultiSliceLib
	are explicitly listed under 2.7. Attachment A of the NVIDIA Software
	License Agreement.

========================================================================