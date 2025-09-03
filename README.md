# NFFT - Nonequispaced Fast Fourier Transforms

## Overview

NFFT is a software library, written in C, for computing non-equispaced fast
Fourier transforms and related variations. It implements the following
transforms:

- **Non-equispaced fast Fourier transform (NFFT)**: A generalization of the 
  [Fast Fourier transform (FFT)](https://en.wikipedia.org/wiki/Fast_Fourier_transform) 
  to non-equispaced nodes in time/space domain. This includes
    - the forward transform, *NFFT* or *forward NFFT*, from frequency into the time/space domain,
    - the adjoint transform, *adjoint NFFT*, from time/space into the frequency domain.
  
  Note that in contrast to the regular FFT, the adjoint transform is generally **not** the inverse of the forward transform, nor does such inverse generally exist. Depending on the number of time/space nodes and the number of frequency nodes, generalized inverse transforms can be defined through solving either an under- or overdetermined system of linear equations. Solving these numerically typically requires iterative application of both the forward and the adjoint transforms.

- **NFCT & NFST**: a specialization of the NFFT to real-valued data, i.e. (co)sine transforms.

- **NSFFT**: sparse NFFT on a subset of viable frequencies, e.g. constrained by a hyperbolic cross.

- **NNFFT**: a further generalization of the NFFT to arbitrary nodes in time/space **and** frequency domains at the same time.

- **NFSFT**: a generalization of the NFFT to the sphere S^2.

- **NFSOFT**: a generalization of the NFFT to the rotation group SO(3).

- **Inverse transforms**: numerical inversion of forward transforms using e.g. CGNR/CGNE algorithms. 

## Applications
Example code for some applications of these transforms is also provided. This includes

- Medical imaging
    - magnetic resonance imaging (mri)
    - computerised tomography (radon)

- Summation schemes
    - fast summation (fastsum)
    - fast Gauss transform (FGT)
    - singular kernels
    - zonal kernels

- polar FFT, discrete Radon transform, ridgelet transform

Detailed API documentation in HTML format can be found in
`doc/html/index.html`, if you are working from a release tarball.
When working from a source repository, the documentation can be
generated with Doxygen (which requires the `doxygen-latex` and `perl` packages):
```shell
make doc
```

## Building

The NFFT depends on the [FFTW](https://fftw.org) library, which is available for many Linux distros, Homebrew on macOS and MSYS2 on Windows. If you compile the FFTW yourself, it should be configured with the flag `--enable-shared` (and `--enable-threads` for the multi-threaded version). Building the NFFT requires `make` and a C compiler such as `gcc`.

When working from a source repository, you need to run libtoolize and autoreconf first. A bash script to do this is provided. This step requries the tools `autoconf`, `automake` and `libtool`.
```shell
./bootstrap.sh
```

The rest of the build process is standard.
```shell
./configure --enable-all --enable-openmp [add options as necessary, see below]
```

Alternatively, you might run the configure script for Matlab.
```shell
./configure --enable-all --enable-openmp --with-matlab=/path/to/matlab
```

Here are some useful optional flags for `./configure`:
* `--enable-all` specifies that all modules should be compiled,
* `--enable-openmp` enables the multicore support and
* `--enable-julia` specifies that the julia interface will be compiled.
* `--with-matlab=/path/to/matlab` specifies the path of a Matlab installation, and
* `--with-octave=/path/to/octave` does the same for GNU Octave.
* For a list of all available options, run `./configure --help`.

Build the software.
```shell
make
```

Optionally, unit tests may be run. Some of the unit tests require an installation of [cunit](http://cunit.sourceforge.net).
```shell
make check
```

Optionally, install NFFT on your system.
```shell
make install
```

## CodSpeed benchmarks

Optionally, NFFT can build benchmark programs using the [CodSpeed](https://codspeed.io) C++ integration library.
[CodSpeed](https://codspeed.io) is a continuous benchmarking service that can help with tracking performance 
regressions and improvements.

While these benchmarks can run locally, they are actually only intended to run in CI. They also require a build of the CodSpeed integration library from [source](https://github.com/CodSpeedHQ/codspeed-cpp).

To enable building the benchmarks with the CodSpeed integration library, run 
```shell
./configure --enable-benchmarks --with-codspeed=<path to codspeed-cpp source directory>`
```

After the build, the benchmrks can be found in the `benchmarks` directory.

## Citing

The current general paper, the one that we recommend if you wish to cite NFFT, is *Keiner, J., Kunis, S., and Potts, D.
''Using NFFT 3 - a software library for various nonequispaced fast Fourier transforms''
ACM Trans. Math. Software 36, Article 19, 1-30, 2009*. BibTeX entry:
```bibtex
@article{KeKuPo09,
 author = {Jens Keiner and Stefan Kunis and Daniel Potts},
 title = {Using {NFFT3} - a Software Library for Various Nonequispaced Fast {Fourier} Transforms},
 journal = {{ACM} Trans. Math. Software},
 year = {2009},
 volume = {36},
 pages = {Article 19, 1--30},
 doi = {10.1145/1555386.1555388}}
```

## Feedback

Your comments are welcome! This is the third version of the library and may
not be as robust or well documented as it should be. Please keep track of bugs
or missing/confusing instructions and report them in our issue tracker or directly to
[Daniel Potts](mailto:potts@mathematik.tu-chemnitz.de).
The postal address is

```text
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

## Legal Information & Credits

Copyright (c) 2002, 2025 Jens Keiner, Stefan Kunis, Daniel Potts

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

## Directory structure

```text
nfft/
├── .github/              # GitHub-related configuration
├── 3rdparty/             # Bundled third-party source code
├── applications/         # Programs to illustrate demonstrate transform applications
├── benchmarks/           # Benchmark programs, intended to be run only in CI
├── doc/                  # User and developer documentation
├── doxygen/              # Doxygen configuration 
├── examples/             # Simple examples demonstrating the use of the NFFT library
├── include/              # Header files
├── julia/                # Julia interface for NFFT
├── kernel/               # Core NFFT library implementation
├── m4                    # Autoconf macros
├── matlab/               # Matlab MEX interfaces for NFFT, NFSFT,, NFSOFT, NNFFT
├── support/              # Support files, for maintainers only
├── tests/                # Unit tests
├── .gitignore            # Git ignore file
├── .aminclude.am         # Automake include file
├── AUTHORS               # Information about the authors of the NFFT library
├── bootstrap.sh          # Bootstrap shell script to call Autoconf and friends
├── ChangeLog             # A short version history
├── configure.ac          # Autoconf configure script template
├── CONVENTIONS           # Internal coding conventions
├── COPYING               # Information about redistributing NFFT
├── doxygen.dox           # Doxygen configuration
├── linux-build-mex.sh    # Script to build MEX files for Linux
├── macos-build-mex.sh    # Script to build MEX files for macOS
├── Makefile.am           # Automake Makefile template
├── NEWS                  # News file (empty, but required by automake)
├── nfft3.pc.in           # pkg-config file template
├── README                # Softlink to README.md
├── README.md             # This file
├── windows-build-dll.sh  # Script to build DLL files for Windows
```
