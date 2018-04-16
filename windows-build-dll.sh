#!/bin/bash

# This script builds statically linked DLL and Octave / Matlab interfaces for Windows.
# It is meant to be called with MSYS2-MinGW64.
# A Matlab installation must be specified in order to build the Matlab interface.
# The Matlab path should not contain spaces!
# 
# Example call:
# ./nfft-build-dll.sh --fftw=3.3.7 --octave=4.2.2 --matlab=/c/path/to/matlab
# 
# WARNING: This script downloads and compiles FFTW and downloads GCC and Octave (requires ~ 2GB).
# 
# Optional flags: --arch=64 (default) or 32 (swich wheater to build 64 or 32 bit binaries)
# To build for 32 bit you should use MinGW-w64 Win32 Shell.
# --gcc-arch (Flag for gcc --march=)

# Any subsequent commands which fail will cause the shell script to exit immediately
set -ex

# default values (to be overwritten if respective parameters are set)
FFTWVERSION=3.3.7
OCTAVEVERSION=4.2.2
ARCH=64
GCCARCH=""

# read the options
TEMP=`getopt -o o:m:f:a:g: --long octave:,matlab:,fftw:,arch:,gcc-arch: -n 'nfft-build-dll.sh' -- "$@"`
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
        -a|--arch)
            case "$2" in
                "")  shift 2 ;;
                *) ARCH=$2; shift 2 ;;
            esac ;;
        -g|--gcc-arch)
            case "$2" in
                "")  shift 2 ;;
                *) GCCARCH=$2; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

if [ "$ARCH" == "32" ]; then
  ARCHNAME=i686
  MATLABARCHFLAG="--with-matlab-arch=win32"
  if [ "$GCCARCH" == "" ]; then
    GCCARCH=pentium4;
  fi
elif [ "$ARCH" == "64" ]; then
  ARCHNAME=x86_64
  MATLABARCHFLAG=""
  if [ "$GCCARCH" == "" ]; then
    GCCARCH=core2;
  fi
else
  echo "Unknown architecture!" ; exit 1 ;
fi

# Install required packages
pacman -S --needed autoconf perl libtool automake mingw-w64-$ARCHNAME-gcc make mingw-w64-$ARCHNAME-cunit mingw-w64-$ARCHNAME-ncurses tar zip unzip wget dos2unix

#NFFTDIR=$(pwd)
NFFTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
HOMEDIR="$NFFTDIR"/windows-build-dll-$ARCHNAME
mkdir -p "$HOMEDIR"
cd "$HOMEDIR"
FFTWDIR=$HOMEDIR/fftw-$FFTWVERSION
GCCVERSION=$(gcc -dumpversion)

# Build FFTW
if [ ! -f "$FFTWDIR/build-success" ]; then
  rm -rf "$FFTWDIR"
  wget "http://fftw.org/fftw-$FFTWVERSION.tar.gz"
  tar -zxf "fftw-$FFTWVERSION.tar.gz"
  rm "fftw-$FFTWVERSION.tar.gz"
  cd "$FFTWDIR"

  mkdir build build-static build-threads build-threads-static
  FFTWFLAGS="--with-pic --with-our-malloc16 --enable-sse2 --enable-avx --with-incoming-stack-boundary=2 --disable-fortran"

  cd build-threads
  ../configure $FFTWFLAGS --with-windows-f77-mangling --enable-shared --enable-threads --with-combined-threads
  make -j4
  cd ../build-threads-static
  ../configure $FFTWFLAGS --host=$ARCHNAME-w64-mingw32.static --disable-fortran --enable-static --disable-shared --enable-threads --with-combined-threads
  make -j4
  cd ../build
  ../configure $FFTWFLAGS --with-windows-f77-mangling --enable-shared
  make -j4
  cd ../build-static
  ../configure $FFTWFLAGS --host=$ARCHNAME-w64-mingw32.static --enable-static --disable-shared
  make -j4
  touch "$FFTWDIR/build-success"
fi


# Get Octave
cd "$HOMEDIR"
OCTAVEDIR="$HOMEDIR/octave-$OCTAVEVERSION"
if [ ! -d "$OCTAVEDIR" ]; then
  rm -f "octave-$OCTAVEVERSION-w$ARCH.zip"
  wget "https://ftp.gnu.org/gnu/octave/windows/octave-$OCTAVEVERSION-w$ARCH.zip"
  unzip -q "octave-$OCTAVEVERSION-w$ARCH.zip"
  rm "octave-$OCTAVEVERSION-w$ARCH.zip"
fi
#OCTLIBDIR=$("$OCTAVEDIR"/bin/octave-config -p OCTLIBDIR)
OCTLIBDIR="$OCTAVEDIR"/lib/octave/"$OCTAVEVERSION"
rm -f "$OCTLIBDIR"/liboctave.la "$OCTLIBDIR"/liboctinterp.la

# Build NFFT
READMECONTENT="
$(sed '/Directory structure/Q' $NFFTDIR/README) 
"
cd "$NFFTDIR"
./bootstrap.sh
make distclean || true

for OMPYN in 0 1
do
if [ $OMPYN = 1 ]; then
  FFTWBUILDDIR="$FFTWDIR/build-threads"
  NFFTBUILDDIR="$HOMEDIR/build-openmp"
  OMPFLAG="--enable-openmp"
  OMPLIBS="-fopenmp -static-libgcc"
  THREADSSUFFIX="_threads"
  OMPSUFFIX="-openmp"
else
  FFTWBUILDDIR="$FFTWDIR/build"
  NFFTBUILDDIR="$HOMEDIR/build"
  OMPFLAG=""
  OMPLIBS=""
  THREADSSUFFIX=""
  OMPSUFFIX=""
fi

rm -f -r "$NFFTBUILDDIR"
mkdir "$NFFTBUILDDIR"
cd "$NFFTBUILDDIR"

"$NFFTDIR/configure" --enable-all $OMPFLAG --with-fftw3-libdir="$FFTWBUILDDIR"/.libs --with-fftw3-includedir="$FFTWDIR"/api --with-octave="$OCTAVEDIR" --with-gcc-arch=$GCCARCH --disable-static --enable-shared
make
make check

# Create DLL release
NFFTVERSION=$( grep 'Version: ' nfft3.pc | cut -c10-)
DLLDIR=nfft-"$NFFTVERSION"-dll$ARCH$OMPSUFFIX

gcc -shared  -Wl,--whole-archive 3rdparty/.libs/lib3rdparty.a kernel/.libs/libkernel$THREADSSUFFIX.a -lfftw3 -lm -L"$FFTWBUILDDIR-static/.libs" -Wl,--no-whole-archive -O3 -malign-double   -o .libs/libnfft3$THREADSSUFFIX-2.dll -Wl,-Bstatic -lwinpthread $OMPLIBS

mkdir "$DLLDIR"
cp ".libs/libnfft3$THREADSSUFFIX-2.dll" "$DLLDIR/libnfft3$THREADSSUFFIX-2.dll"
cp "$NFFTDIR"/include/nfft3.h "$DLLDIR"/nfft3.h
cp "$NFFTDIR"/include/nfft3mp.h "$DLLDIR"/nfft3mp.h
cp examples/nfft/simple_test.c "$DLLDIR"/simple_test.c
echo 'This archive contains the NFFT' $NFFTVERSION 'library and the associated header files.
The NFFT library was compiled with double precision support for' $ARCH'-bit Windows
using GCC' $GCCVERSION $ARCHNAME'-w64-mingw32 with march='$GCCARCH 'and FFTW' $FFTWVERSION'.

As a small example, you can compile the NFFT simple test with the following command

	gcc -O3 simple_test.c -o simple_test.exe -L. libnfft3'$THREADSSUFFIX'-2.dll
' "$READMECONTENT" "
FFTW
----
The compiled NFFT files contain parts of the FFTW library (http://www.fftw.org)
Copyright (c) 2003, 2007-14 Matteo Frigo
Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology" > "$DLLDIR"/readme-windows.txt
unix2dos "$DLLDIR"/readme-windows.txt
cp "$NFFTDIR"/COPYING "$DLLDIR"/COPYING
rm -f "$HOMEDIR/$DLLDIR".zip
zip -r -q -9 "$HOMEDIR/$DLLDIR".zip "$DLLDIR"


# Compile with Matlab
if [ -n "$MATLABDIR" ]; then
  "$MATLABDIR"/bin/matlab -wait -nodesktop -nosplash -nodisplay -r "fid=fopen('matlab_version.txt','wt'); fprintf(fid,'MATLAB_VERSION=%s\n', version); exit;" 
  MATLABSTRING=" and Matlab `grep MATLAB_VERSION matlab_version.txt | sed 's/.*(//' | sed 's/)//'`"
  cd "$NFFTBUILDDIR"
  "$NFFTDIR/configure" --enable-all $OMPFLAG --with-fftw3-libdir="$FFTWBUILDDIR"/.libs --with-fftw3-includedir="$FFTWDIR"/api --with-matlab="$MATLABDIR" "$MATLABARCHFLAG" --with-gcc-arch=$GCCARCH --disable-static --enable-shared
  make
  if [ "$ARCH" = "64" ]; then
    make check
  fi
  for LIB in nfft nfsft nfsoft nnfft fastsum nfct nfst
  do
    cd matlab/"$LIB"
    gcc -shared  .libs/lib"$LIB"_la-"$LIB"mex.o  -Wl,--whole-archive ../../.libs/libnfft3_matlab.a ../../matlab/.libs/libmatlab.a -Wl,--no-whole-archive  -L"$MATLABDIR"/bin/win$ARCH -lmwfftw3 -lmx -lmex -lmat -O3 -malign-double -march=$GCCARCH -o .libs/lib"$LIB".mexw$ARCH -Wl,--enable-auto-image-base -Xlinker --out-implib -Xlinker .libs/lib"$LIB".dll.a -static-libgcc -Wl,-Bstatic -lwinpthread $OMPLIBS
    cp .libs/lib"$LIB".mexw$ARCH "$LIB"mex.mexw$ARCH
    cd ../..
  done
fi

# Create Matlab/Octave release
MEXDIR=nfft-"$NFFTVERSION"-mexw$ARCH$OMPSUFFIX
for SUBDIR in nfft nfsft nfsft/@f_hat nfsoft nnfft fastsum nfct nfst infft1d
  do
  mkdir -p "$MEXDIR"/$SUBDIR
  cp -f -r matlab/$SUBDIR/*.mex* "$MEXDIR"/$SUBDIR/ || true
  cp -f -r "$NFFTDIR"/matlab/$SUBDIR/*.m "$MEXDIR"/$SUBDIR/
done
cp "$NFFTDIR"/matlab/infft1d/README "$MEXDIR"/infft1d/

cd "$NFFTBUILDDIR"
cp "$NFFTDIR"/COPYING "$MEXDIR"/COPYING
echo 'This archive contains the Matlab and Octave interface of NFFT' $NFFTVERSION 'compiled for
'$ARCH'-bit Windows using GCC' $GCCVERSION $ARCHNAME'-w64-mingw32 with march='$GCCARCH', FFTW' $FFTWVERSION$MATLABSTRING' and Octave '$OCTAVEVERSION'.
' "$READMECONTENT" > "$MEXDIR"/readme-matlab.txt
unix2dos "$MEXDIR"/readme-matlab.txt
rm -f "$HOMEDIR/$MEXDIR".zip
zip -r -q -9 "$HOMEDIR/$MEXDIR".zip "$MEXDIR"

done
