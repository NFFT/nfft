Master: [![Build Status](https://travis-ci.org/NFFT/nfft_new.svg?branch=master)](https://travis-ci.org/NFFT/nfft_new)
Develop: [![Build Status](https://travis-ci.org/NFFT/nfft_new.svg?branch=develop)](https://travis-ci.org/NFFT/nfft_new)

NFFT - Nonequispaced FFT
=========================

Overview
--------
NFFT is a software library, written in C, for computing non-equispaced fast
Fourier transforms and related variations. It implements the following
transforms:

1. Non-equispaced fast Fourier transform (NFFT)
    - forward transform *(NFFT)*, i.e. frequency to time/space domain
    - adjoint transform *(adjoint NFFT)*, i.e. time/space to frequency domain

2. Generalisations
    - to arbitrary nodes in time *and* frequency domain *(NNFFT)*
    - to real-valued data, i.e. (co)sine transforms, *(NFCT, NFST)*
    - to the sphere S^2 *(NFSFT)*
    - to the rotation group *(NFSOFT)*
    - to the hyperbolic cross *(NSFFT)*

3. Generalised inverse transformations based on iterative methods, e.g. CGNR/CGNE

Some examples for application of these transforms are provided:

1. Medical imaging
    - magnetic resonance imaging
    - computerised tomography

2. Summation schemes
    - fast Gauss transform (FGT)
    - singular kernels
    - zonal kernels

3. polar FFT, discrete Radon transform, ridgelet transform

Detailed API documentation in HTML format can be found in
`doc/api/html/index.html`, if you are working from a release tarball.
When working from a source repository, the documentation can be
generated with Doxygen.
```
make doc
```

Building
--------
When working from a source repository, your need to run libtoolize and autoreconf first. A bash script to do this is provided.
```
./bootstrap.sh
```

The rest of the build process is standard.
```
./configure (add options as necessary)
make
make install
```

Optionally, unit tests may be run.
```
make check
```

Citing
------
The most current general paper, the one that we recommend if you wish to cite NFFT, is *Keiner, J., Kunis, S., and Potts, D.
''Using NFFT 3 - a software library for various nonequispaced fast Fourier transforms''
ACM Trans. Math. Software,36, Article 19, 1-30, 2009*.

Directory structure
-------------------

File/Folder        | Purpose
------------------:| ------------------------------------------------------
3rdparty (dir)	   | Third-party source code
aclocal.m4		   | Macros for configure script
applications (dir) | Application programs (see 4) above)
AUTHORS			   | Information about the authors of NFFT
bootstrap.sh       | Bootstrap shell script that call Autoconf and friends
ChangeLog          | A short version history
config.guess       | Used by configure script
config.sub         | Used by configure script
configure          | Configure script (created by calling ./bootstrap.sh)
configure.in       | Autoconf configure script template
CONVENTIONS        | Internal coding conventions
COPYING            | Information about redistributing NFFT
depcomp            | Used by configure script
doc (dir)          | User and developer documentation
examples (dir)     | Simple examples for using NFFT routines
include (dir)      | Header files
INSTALL            | Installation instructions
install-sh         | Used by configure script
kernel (dir)       | Source code for core library routines
ltmain.sh          | Used by configure script
Makefile.am        | Automake Makefile template
Makefile.in        | Makefile template generated from Makefile.am, processed by configure script
matlab (dir)       | Matlab MEX interfaces for nfft, nfsft, nfsoft, nfft
missing            | Used by configure script
NEWS               | New and noteworthy
README             | This file
tests (dir)        | CUnit tests
TODO               | Work to be done

Feedback
--------
Your comments are welcome! This is the third version of the library and may
not be as robust or well documented as it should be. Please keep track of bugs
or missing/confusing instructions and report them to

  Prof. Dr. Daniel Potts <potts@mathematik.tu-chemnitz.de>
  TU Chemnitz, Fakultaet fuer Mathematik
  Reichenhainer Str. 39
  09107 Chemnitz, GERMANY

or

  Stefan Kunis           <stefan.kunis@math.uos.de>
  Jens Keiner            <jens@nfft.org>>

If you find NFFT useful, we would be delighted to hear about what application
you are using NFFT for!

Legal Information & Credits
---------------------------
Copyright (c) 2002, 2016 Jens Keiner, Stefan Kunis, Daniel Potts

This software was written by Jens Keiner, Stefan Kunis and Daniel Potts.
It was developed at the Mathematical Institute, University of
Luebeck, and at the Faculty of Mathematics, Chemnitz University of Technology.

NFFT3 is free software. You can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. If not stated otherwise, this applies to all files contained in this
package and its sub-directories.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
