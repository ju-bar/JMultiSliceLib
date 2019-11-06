========================================================================
	STATIC LIBRARY : JMultiSliceLib Project Overview
========================================================================

	Copyright (C) 2018, 2019 - Juri Barthel
	Copyright (C) 2018, 2019 - Forschungszentrum Juelich GmbH, Germany
	Email: juribarthel@gmail.com

	Last Version: 0.33 - 2019-10-30

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

* mkl_intel_lp64.lib, mkl_sequential.lib, and mkl_core.lib
  (see MKL statement of license X.2)
  JMultiSliceLib uses data output by the FFT code of the Intel Math
  Kernel Library (MKL) and is in no form based on work represented by
  the Intel MKL. Intel MKL binary code is available from 
  https://software.intel.com/en-us/mkl (accessed Nov 2019)

* cudart_static.lib and cufft.lib (see CUDA statement of license X.3)
  JMultiSliceLib uses data output and functionality probided by these
  two libraries from the CUDA Toolkit 9.0 as redistributable software.
  The CUDA Toolkit is available from
  https://developer.nvidia.com/cuda-toolkit (accessed April 2018)
  In order to be able to use the GPU routnines of JMultiSliceLib, keep
  your CUDA device drivers updated.

You should link these libraries in the calling code and take care to
comply with the FFTW and CUDA software licenses.


4) Statement of license

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
	
X.2) Intel MKL statement of license

    Source: https://software.intel.com/en-us/mkl/license-faq

	The Intel MKL is used according to the Intel Simplified Software
	License. A copy of the license from April 2018 is made available
	with the distribution of JMultiSliceLib as
	"Intel-Simplified-Software-License.pdf".

X.3) CUDA statement of license

	Source: https://docs.nvidia.com/cuda/eula/index.html 

	CUDA is used according to the NVIDIA Software License Agreement
	linked above and made available with JMultiSliceLib as
	"CUDA EULA.pdf". The two CUDA libraries used by JMultiSliceLib
	are explicitly listed under 2.7. Attachment A of the NVIDIA Software
	License Agreement.

========================================================================