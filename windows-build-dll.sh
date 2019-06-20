#!/bin/bash

# This script builds statically linked DLL and Octave / Matlab interfaces for Windows.
# It is meant to be called with MSYS2-MinGW64.
# A Matlab installation must be specified in order to build the Matlab interface.
# The Matlab path should not contain spaces!
# 
# Example call:
# ./nfft-build-dll.sh --fftw=3.3.8 --octave=4.4.1 --matlab=/c/path/to/matlab
# 
# WARNING: This script downloads and compiles FFTW and downloads GCC and Octave (requires ~ 2GB).
# 
# Optional flags: --arch=64 (default) or 32 (swich wheater to build 64 or 32 bit binaries)
# To build for 32 bit you should use MinGW-w64 Win32 Shell.
# --gcc-arch (Flag for gcc --march=)

# Any subsequent commands which fail will cause the shell script to exit immediately
set -ex

# Write log file
# exec > >(tee -a windows-build-dll.log)
# exec 2>&1


# default values (to be overwritten if respective parameters are set)
FFTWVERSION=3.3.8
OCTAVEVERSION=4.4.1
MATLABVERSION=""
ARCH=64
GCCARCH=""
SOVERSION=4

# read the options
TEMP=`getopt -o o:m:f:a:g:v: --long octave:,matlab:,fftw:,arch:,gcc-arch:,matlab-version: -n 'nfft-build-dll.sh' -- "$@"`
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
        -v|--matlab-version)
            case "$2" in
                "")  shift 2 ;;
                *) MATLABVERSION="$2"; shift 2 ;;
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
pacman -S --needed autoconf perl libtool automake mingw-w64-$ARCHNAME-gcc make mingw-w64-$ARCHNAME-cunit tar zip unzip wget dos2unix rsync p7zip

#NFFTDIR=$(pwd)
NFFTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
HOMEDIR="$(pwd)"/windows-build-dll-$ARCHNAME
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
  FFTWFLAGS="--with-pic --with-our-malloc16 --enable-sse2 --enable-avx --enable-avx2 --with-incoming-stack-boundary=2 --disable-fortran"

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
  rm -f -r "octave-$OCTAVEVERSION-w$ARCH"
  wget "https://ftp.gnu.org/gnu/octave/windows/octave-$OCTAVEVERSION-w$ARCH.zip"
  unzip -q -o "octave-$OCTAVEVERSION-w$ARCH.zip"
  rm "octave-$OCTAVEVERSION-w$ARCH.zip"
  mv "$OCTAVEDIR-w$ARCH" "$OCTAVEDIR" || true	# Folder name has suffix -w64 for Octave >=4.4
fi
#OCTLIBDIR=$("$OCTAVEDIR"/bin/octave-config -p OCTLIBDIR)
# check which folder contains the Octave binaries
if [ -f "$OCTAVEDIR"/mingw$ARCH/bin/octave-cli.exe ]; then
  OCTAVEDIR="$OCTAVEDIR"/mingw$ARCH
fi
# remove files that prevent compilation
OCTLIBDIR="$OCTAVEDIR"/lib/octave/"$OCTAVEVERSION"
rm -f "$OCTLIBDIR"/liboctave.la "$OCTLIBDIR"/liboctinterp.la

# Build NFFT
READMECONTENT="
$(sed '/Directory structure/Q' $NFFTDIR/README) 
"
FFTWREADME='
FFTW
----
The compiled NFFT files contain parts of the FFTW library (http://www.fftw.org)
Copyright (c) 2003, 2007-14 Matteo Frigo
Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology'
cd "$NFFTDIR"
#./bootstrap.sh
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


# Compile with Octave
"$NFFTDIR/configure" --enable-all $OMPFLAG --with-fftw3-libdir="$FFTWBUILDDIR"/.libs --with-fftw3-includedir="$FFTWDIR"/api --with-octave="$OCTAVEDIR" --with-gcc-arch=$GCCARCH --disable-static --enable-shared --disable-applications --disable-examples
make
make check

# Create DLL release
NFFTVERSION=$( grep 'Version: ' nfft3.pc | cut -c10-)
DLLDIR=nfft-"$NFFTVERSION"-dll$ARCH$OMPSUFFIX

gcc -shared  -Wl,--whole-archive 3rdparty/.libs/lib3rdparty.a kernel/.libs/libkernel$THREADSSUFFIX.a -lfftw3 -lm -L"$FFTWBUILDDIR-static/.libs" -Wl,--no-whole-archive -O3 -malign-double   -o .libs/libnfft3$THREADSSUFFIX-$SOVERSION.dll -Wl,-Bstatic -lwinpthread $OMPLIBS -Wl,--output-def,.libs/libnfft3$THREADSSUFFIX-$SOVERSION.def

mkdir "$DLLDIR"
cp ".libs/libnfft3$THREADSSUFFIX-$SOVERSION.dll" "$DLLDIR/libnfft3$THREADSSUFFIX-$SOVERSION.dll"
cp ".libs/libnfft3$THREADSSUFFIX-$SOVERSION.def" "$DLLDIR/libnfft3$THREADSSUFFIX-$SOVERSION.def"
cp "$NFFTDIR"/include/nfft3.h "$DLLDIR"/nfft3.h
cp "$NFFTDIR"/include/nfft3mp.h "$DLLDIR"/nfft3mp.h
cp "$FFTWDIR"/api/fftw3.h "$DLLDIR"/fftw3.h
cp examples/nfft/simple_test.c "$DLLDIR"/simple_test.c
echo 'This archive contains the NFFT' $NFFTVERSION 'library and the associated header files.
The NFFT library was compiled with double precision support for' $ARCH'-bit Windows
using GCC' $GCCVERSION $ARCHNAME'-w64-mingw32 with march='$GCCARCH 'and FFTW' $FFTWVERSION'.

In order to link the .dll file from Visual C++, you should run
    lib /def:libnfft3'$THREADSSUFFIX'-'$SOVERSION'.def

As a small example, you can compile the NFFT simple test with the following command

	gcc -O3 simple_test.c -o simple_test.exe -I. -L. libnfft3'$THREADSSUFFIX'-'$SOVERSION'.dll
' "$READMECONTENT" "$FFTWREADME" > "$DLLDIR"/readme-windows.txt
unix2dos "$DLLDIR"/readme-windows.txt
cp "$NFFTDIR"/COPYING "$DLLDIR"/COPYING
rm -f "$HOMEDIR/$DLLDIR".zip
7z a -r "$HOMEDIR/$DLLDIR".zip "$DLLDIR"


# Compile with Matlab
if [ -n "$MATLABDIR" ]; then
  if [ "$MATLABVERSION" == "" ]; then
    "$MATLABDIR"/bin/matlab -wait -nodesktop -nosplash -nodisplay -r "fid=fopen('matlab_version.txt','wt'); fprintf(fid,'MATLAB_VERSION=%s\n', version); exit;" 
    MATLABVERSION=" and Matlab `grep MATLAB_VERSION matlab_version.txt | sed 's/.*(//' | sed 's/)//'`"
  fi
  MATLABSTRING=" and Matlab $MATLABVERSION"
  cd "$NFFTBUILDDIR"
  "$NFFTDIR/configure" --enable-all $OMPFLAG --with-fftw3-libdir="$FFTWBUILDDIR"/.libs --with-fftw3-includedir="$FFTWDIR"/api --with-matlab="$MATLABDIR" "$MATLABARCHFLAG" --with-gcc-arch=$GCCARCH --disable-static --enable-shared --disable-applications --disable-examples
  make
  if [ -f "$MATLABDIR"/bin/matlab.exe ]; then
    make check
  fi
  for LIB in nfft nfsft nfsoft nnfft fastsum nfct nfst fpt
  do
    cd matlab/"$LIB"
    gcc -shared  .libs/lib"$LIB"_la-"$LIB"mex.o  -Wl,--whole-archive ../../.libs/libnfft3_matlab.a ../../matlab/.libs/libmatlab.a -Wl,--no-whole-archive  -L"$MATLABDIR"/bin/win$ARCH -lmwfftw3 -lmx -lmex -lmat -O3 -malign-double -march=$GCCARCH -o .libs/lib"$LIB".mexw$ARCH -Wl,--enable-auto-image-base -Xlinker --out-implib -Xlinker .libs/lib"$LIB".dll.a -static-libgcc -Wl,-Bstatic -lwinpthread $OMPLIBS
    cp .libs/lib"$LIB".mexw$ARCH "$LIB"mex.mexw$ARCH
    cd ../..
  done
fi

# Create Matlab/Octave release
MEXDIR=nfft-"$NFFTVERSION"-mexw$ARCH$OMPSUFFIX
for SUBDIR in nfft nfsft/@f_hat nfsft nfsoft nnfft fastsum nfct nfst infft1d fpt
  do
  mkdir -p "$MEXDIR"/$SUBDIR
  cp -f -r matlab/$SUBDIR/*.mex* "$MEXDIR"/$SUBDIR/ || true
  cp -f -r "$NFFTDIR"/matlab/$SUBDIR/README "$MEXDIR"/$SUBDIR/ || true
  cp -f -r "$NFFTDIR"/matlab/$SUBDIR/*.m "$MEXDIR"/$SUBDIR/
  "$OCTAVEDIR"/bin/octave-cli --no-window-system --eval="cd $MEXDIR/$SUBDIR; if exist('simple_test')==2; simple_test; end; if exist('test_$SUBDIR')==2; test_$SUBDIR; end"
done

cd "$NFFTBUILDDIR"
cp "$NFFTDIR"/COPYING "$MEXDIR"/COPYING
echo 'This archive contains the Matlab and Octave interface of NFFT' $NFFTVERSION 'compiled for
'$ARCH'-bit Windows using GCC' $GCCVERSION $ARCHNAME'-w64-mingw32 with march='$GCCARCH', FFTW' $FFTWVERSION$MATLABSTRING' and Octave '$OCTAVEVERSION'.
' "$READMECONTENT" > "$MEXDIR"/readme-matlab.txt
unix2dos "$MEXDIR"/readme-matlab.txt
rm -f "$HOMEDIR/$MEXDIR".zip
7z a -r "$HOMEDIR/$MEXDIR".zip "$MEXDIR"


# Build applications
cd "$NFFTBUILDDIR"
CPPFLAGS="--static" LDFLAGS="--static" "$NFFTDIR/configure" --enable-all $OMPFLAG --with-fftw3-libdir="$FFTWBUILDDIR-static"/.libs --with-fftw3-includedir="$FFTWDIR"/api --with-gcc-arch=$GCCARCH --enable-static --disable-shared
make all check

APPSDIR=nfft-"$NFFTVERSION"-applications-w$ARCH$OMPSUFFIX
mkdir "$APPSDIR"
rsync -a --exclude='Makefile*' --exclude='.deps' --exclude='.libs' --exclude='*.la' --exclude='*.lo' --exclude='*.o' 'applications/' "$APPSDIR"
rsync -a --exclude='Makefile*' --exclude='doxygen*' --exclude='*.c.in' "$NFFTDIR/applications/" "$APPSDIR"

echo 'This archive contains the NFFT' $NFFTVERSION 'applications.
The NFFT library was compiled with double precision support for' $ARCH'-bit Windows
using GCC' $GCCVERSION $ARCHNAME'-w64-mingw32 with march='$GCCARCH 'and FFTW' $FFTWVERSION'.
' "$READMECONTENT" "$FFTWREADME" > "$APPSDIR"/readme.txt
unix2dos "$APPSDIR"/readme.txt
rm -f "$HOMEDIR/$APPSDIR".zip
7z a -r "$HOMEDIR/$APPSDIR".zip "$APPSDIR"

done
