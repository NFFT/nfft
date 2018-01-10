#!/bin/sh

# This script builds statically linked DLL and Octave / Matlab interfaces for Windows.
# It is meant to be called with MSYS2-MinGW64.
# A Matlab installation must be specified in order to build the Matlab interface.
# The Matlab path should not contain spaces!
# 
# Example call:
# ./nfft-build-dll.sh --fftw=3.3.7 --octave=4.2.1 --matlab=/c/path/to/matlab
# 
# WARNING: This script downloads and compiles FFTW and downloads GCC and Octave (requires ~ 2GB).

# Any subsequent commands which fail will cause the shell script to exit immediately
set -ex

# default values (to be overwritten if respective parameters are set)
FFTWVERSION="3.3.7"
OCTAVEVERSION='4.2.1'

# read the options
TEMP=`getopt -o o:m:f: --long octave:,matlab:,fftw: -n 'nfft-build-dll.sh' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -o|--octave)
            case "$2" in
                "")  shift 2 ;;
                *) OCTAVEVERSION=$2; shift 2 ;;
            esac ;;
        -m|--matlab)
            case "$2" in
                "")  shift 2 ;;
                *) MATLABDIR=$2; shift 2 ;;
            esac ;;
        -f|--fftw)
            case "$2" in
                "")  shift 2 ;;
                *) FFTWVERSION=$2; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

# Install required packages
pacman -S --needed autoconf perl libtool automake mingw-w64-x86_64-gcc make tar zip unzip wget dos2unix

HOMEDIR=$(pwd)
FFTWDIR="$HOMEDIR"/fftw-$FFTWVERSION
GCCVERSION=$(gcc -dumpversion)
READMECONTENT='
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

Citing
------
The most current general paper, the one that we recommend if you wish to cite NFFT, is *Keiner, J., Kunis, S., and Potts, D.
Using NFFT 3 - a software library for various nonequispaced fast Fourier transforms
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

FFTW
----
The compiled NFFT files contain parts of the FFTW library (http://www.fftw.org)
Copyright (c) 2003, 2007-14 Matteo Frigo
Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology'

# Build FFTW
if [ ! -f "$FFTWDIR/build-success" ]; then
  rm -rf "$FFTWDIR"
  wget http://fftw.org/fftw-$FFTWVERSION.tar.gz
  tar -zxf fftw-$FFTWVERSION.tar.gz
  rm fftw-$FFTWVERSION.tar.gz
  cd "$FFTWDIR"

  mkdir build build-static build-threads build-threads-static
  FFTWFLAGS="--with-pic --with-our-malloc16 --enable-sse2 --enable-avx --with-incoming-stack-boundary=2 --disable-fortran"

  cd build-threads
  ../configure $FFTWFLAGS --with-windows-f77-mangling --enable-shared --enable-threads --with-combined-threads
  make -j4
  cd ../build-threads-static
  ../configure $FFTWFLAGS --host=x86_64-w64-mingw32.static --disable-fortran --enable-static --disable-shared --enable-threads --with-combined-threads
  make -j4
  cd ../build
  ../configure $FFTWFLAGS --with-windows-f77-mangling --enable-shared
  make -j4
  cd ../build-static
  ../configure $FFTWFLAGS --host=x86_64-w64-mingw32.static --enable-static --disable-shared
  make -j4
  touch "$FFTWDIR/build-success"
fi


# Get Octave
cd "$HOMEDIR"
OCTAVEDIR="$HOMEDIR/octave-$OCTAVEVERSION"
if [ ! -d "$OCTAVEDIR" ]; then
  rm -f "octave-$OCTAVEVERSION-w64.zip"
  wget https://ftp.gnu.org/gnu/octave/windows/octave-$OCTAVEVERSION-w64.zip
  unzip -q octave-$OCTAVEVERSION-w64.zip
  rm "octave-$OCTAVEVERSION-w64.zip"
fi
OCTLIBDIR=$("$OCTAVEDIR"/bin/octave-config -p OCTLIBDIR)
rm -f "$OCTLIBDIR"/liboctave.la "$OCTLIBDIR"/liboctinterp.la

# Build NFFT
NFFTDIR="$HOMEDIR"
./bootstrap.sh
make distclean || true

for OMPYN in 0 1
do
if [ $OMPYN = 1 ]; then
  FFTWBUILDDIR="$FFTWDIR/build-threads"
  BUILDDIR="$NFFTDIR/build-openmp"
  OMPFLAG="--enable-openmp"
  OMPLIBS="-fopenmp -static-libgcc"
  THREADSSUFFIX="_threads"
  OMPSUFFIX="-openmp"
else
  FFTWBUILDDIR="$FFTWDIR/build"
  BUILDDIR="$NFFTDIR/build"
  OMPFLAG=""
  OMPLIBS=""
  THREADSSUFFIX=""
  OMPSUFFIX=""
fi

rm -f -r "$BUILDDIR"
mkdir "$BUILDDIR"
cd "$BUILDDIR"

../configure --enable-all $OMPFLAG --with-fftw3-libdir="$FFTWBUILDDIR"/.libs --with-fftw3-includedir="$FFTWDIR"/api --with-octave="$OCTAVEDIR" --with-gcc-arch=core2 --disable-static --enable-shared
make


# Create DLL release
NFFTVERSION=$( grep 'Version: ' nfft3.pc | cut -c10-)
DIR=nfft-$NFFTVERSION-core2$OMPSUFFIX

gcc -shared  -Wl,--whole-archive 3rdparty/.libs/lib3rdparty.a kernel/.libs/libkernel$THREADSSUFFIX.a -lfftw3 -lm -L"$FFTWBUILDDIR-static/.libs" -Wl,--no-whole-archive -O3 -malign-double   -o .libs/libnfft3$THREADSSUFFIX-2.dll -Wl,-Bstatic -lwinpthread $OMPLIBS

mkdir $DIR-dll64
cp .libs/libnfft3$THREADSSUFFIX-2.dll $DIR-dll64/libnfft3$THREADSSUFFIX-2.dll
cp "$NFFTDIR"/include/nfft3.h $DIR-dll64/nfft3.h
cp "$NFFTDIR"/include/nfft3mp.h $DIR-dll64/nfft3mp.h
cp examples/nfft/simple_test.c $DIR-dll64/simple_test.c
echo 'NFFT - Nonequispaced FFT
========================

This archive contains the NFFT' $NFFTVERSION 'library and the associated header files.
The NFFT library was compiled with double precision support for 64-bit Windows
using GCC' $GCCVERSION 'x86_64-w64-mingw32 and FFTW' $FFTWVERSION '.

As a small example, you can compile the NFFT simple test with the following command

	gcc -O3 simple_test.c -o simple_test.exe -L. libnfft3'$THREADSSUFFIX'-2.dll
' "$READMECONTENT" > $DIR-dll64/readme-windows.txt
unix2dos $DIR-dll64/readme-windows.txt
cp "$NFFTDIR"/COPYING $DIR-dll64/COPYING
zip -r -q $DIR-dll64.zip $DIR-dll64


# Create Octave release
cp --recursive matlab octave-$DIR
cp --recursive --update ../matlab/* octave-$DIR

for SUBDIR in . nfft nfsft nfsft/@f_hat nfsoft nnfft fastsum nfct nfst
do
  cd octave-$DIR/$SUBDIR
  rm -f -r .deps .libs
  rm -f Makefile* README *.lo *.la *.c *.h
  cd "$BUILDDIR"
done

cp "$NFFTDIR"/COPYING octave-$DIR/COPYING
echo 'NFFT - Nonequispaced FFT
========================

This archive contains the Octave interface of NFFT' $NFFTVERSION 'compiled for 64-bit Windows
using GCC' $GCCVERSION 'x86_64-w64-mingw32 and Octave '$OCTAVEVERSION'.
' "$READMECONTENT" > octave-$DIR/readme-octave.txt
unix2dos octave-$DIR/readme-octave.txt
zip -r -q octave-$DIR.zip octave-$DIR


# Create Matlab release
if [ -n "$MATLABDIR" ]; then
  ../configure --enable-all $OMPFLAG --with-fftw3-libdir="$FFTWBUILDDIR"/.libs --with-fftw3-includedir="$FFTWDIR"/api --with-matlab="$MATLABDIR" --with-gcc-arch=core2 --disable-static --enable-shared
  make
  for LIB in nfft nfsft nfsoft nnfft fastsum nfct nfst
  do
    cd matlab/"$LIB"
    gcc -shared  .libs/lib"$LIB"_la-"$LIB"mex.o  -Wl,--whole-archive ../../.libs/libnfft3_matlab.a ../../matlab/.libs/libmatlab.a -Wl,--no-whole-archive  -L"$MATLABDIR"/bin/win64 -lmwfftw3 -lmx -lmex -lmat -O3 -malign-double -march=core2 -o .libs/lib"$LIB".mexw64 -Wl,--enable-auto-image-base -Xlinker --out-implib -Xlinker .libs/lib"$LIB".dll.a -static-libgcc -Wl,-Bstatic -lwinpthread $OMPLIBS
    cp .libs/lib"$LIB".mexw64 "$LIB"mex.mexw64
    cd ../..
  done

  cp --recursive matlab matlab-$DIR
  cp --recursive --update ../matlab/* matlab-$DIR

  for SUBDIR in . nfft nfsft nfsft/@f_hat nfsoft nnfft fastsum nfct nfst
  do
    cd matlab-$DIR/$SUBDIR
    rm -f -r .deps .libs
    rm -f Makefile* README *.lo *.la *.c *.h
    cd "$BUILDDIR"
  done
fi

cp "$NFFTDIR"/COPYING matlab-$DIR/COPYING
echo 'NFFT - Nonequispaced FFT
========================

This archive contains the Matlab interface of NFFT' $NFFTVERSION 'compiled for 64-bit Windows
using GCC' $GCCVERSION 'x86_64-w64-mingw32 and FFTW' $FFTWVERSION'.
' "$READMECONTENT" > matlab-$DIR/readme-matlab.txt
unix2dos matlab-$DIR/readme-matlab.txt
zip -r -q matlab-$DIR.zip matlab-$DIR

done