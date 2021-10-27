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
    - magnetic resonance imaging (mri)
    - computerised tomography (radon)

2. Summation schemes
    - fast summation (fastsum)
    - fast Gauss transform (FGT)
    - singular kernels
    - zonal kernels

3. polar FFT, discrete Radon transform, ridgelet transform

Detailed API documentation in HTML format can be found in
`doc/html/index.html`, if you are working from a release tarball.
When working from a source repository, the documentation can be
generated with Doxygen.
```
make doc
```

Building
--------
The NFFT depends on the [FFTW](https://fftw.org) library, which is available for many Linux distros, Homebrew on macOS and MSYS2 on Windows. If you compile the FFTW yourself, it should be configured with the flag `--enable-shared`.

When working from a source repository, you need to run libtoolize and autoreconf first. A bash script to do this is provided.
```
./bootstrap.sh
```

The rest of the build process is standard.
```
./configure --enable-all --enable-openmp [add options as necessary, see below]
```

Alternatively, you might run the configure script for Matlab.
```
./configure --enable-all --enable-openmp --with-matlab=/path/to/matlab
```

Here are some useful optional flags for `./configure`:
* `--enable-all` specifies that all modules should be compiled,
* `--enable-openmp` enables the multicore support and
* `--enable-julia` specifies that the julia interface will be compiled.
* `--with-matlab=/path/to/matlab` specifies a path of Matlab, and
* `--with-octave=/path/to/octave` does the same for GNU Octave.
* For a list of all available options, run `./configure --help`.

Build the software.
```
make
```

Optionally, unit tests may be run. Some of the unit tests require an installation of [cunit](http://cunit.sourceforge.net).
```
make check
```

Optionally, install NFFT on your system.
```
make install
```

Citing
------
The most current general paper, the one that we recommend if you wish to cite NFFT, is *Keiner, J., Kunis, S., and Potts, D.
''Using NFFT 3 - a software library for various nonequispaced fast Fourier transforms''
ACM Trans. Math. Software,36, Article 19, 1-30, 2009*.

Feedback
--------
Your comments are welcome! This is the third version of the library and may
not be as robust or well documented as it should be. Please keep track of bugs
or missing/confusing instructions and report them to
[Daniel Potts](mailto:potts@mathematik.tu-chemnitz.de).
The postal address is

```
  Prof. Dr. Daniel Potts
  TU Chemnitz, Fakultaet fuer Mathematik
  Reichenhainer Str. 39
  09107 Chemnitz
  GERMANY
```

Alternatively, you might contact
[Stefan Kunis](mailto:stefan.kunis@math.uos.de)
or
[Jens Keiner](mailto:jens@nfft.org).

If you find NFFT useful, we would be delighted to hear about what application
you are using NFFT for!

Legal Information & Credits
---------------------------
Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts

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

Directory structure
-------------------

File/Folder        | Purpose
------------------:| ------------------------------------------------------
3rdparty (dir) 	   | Third-party source code
aclocal.m4		   | Macros for configure script
applications (dir) | Application programs (see 4) above)
AUTHORS			   | Information about the authors of NFFT
bootstrap.sh       | Bootstrap shell script that call Autoconf and friends
ChangeLog          | A short version history
config (dir)       | Used by configure script
configure          | Configure script (created by calling ./bootstrap.sh)
configure.ac       | Autoconf configure script template
CONVENTIONS        | Internal coding conventions
COPYING            | Information about redistributing NFFT
doc (dir)          | User and developer documentation
examples (dir)     | Simple examples for using NFFT routines
include (dir)      | Header files
INSTALL            | Installation instructions
julia (dir)        | Julia interface for nfft
kernel (dir)       | Source code for core library routines
Makefile.am        | Automake Makefile template
Makefile.in        | Makefile template generated from Makefile.am, processed by configure script
matlab (dir)       | Matlab MEX interfaces for nfft, nfsft, nfsoft, nfft
NEWS               | New and noteworthy
README             | This file
README.md          | This file
tests (dir)        | CUnit tests
