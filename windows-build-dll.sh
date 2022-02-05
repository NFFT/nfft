#!/bin/bash

# This script builds statically linked DLL and Octave / Matlab interfaces for Windows.
# It is meant to be called with MSYS2-MinGW64.
# A Matlab installation must be specified in order to build the Matlab interface.
# The Matlab path should not contain spaces!
# 
# Example call:
# ./nfft-build-dll.sh --fftw=3.3.10 --octave=6.4.0 --matlab=/c/path/to/matlab
# 
# WARNING: This script downloads and compiles FFTW and downloads GCC, Julia and Octave (requires ~ 3GB).
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
FFTWVERSION=3.3.10
OCTAVEVERSION=6.4.0
MATLABVERSION=""
ARCH=64
GCCARCH=""
BINARIES_ARCH_README='
Please note that since the binaries were compiled with gcc flag -march=haswell,
they may not work on older CPUs (below Intel i3/i5/i7-4xxx or
AMD Excavator/4th gen Bulldozer) as well as on some Intel Atom/Pentium CPUs.
'
SOVERSION=4
OMP_ONLY=

# read the options
TEMP=`getopt -o o:m:f:a:g:v:s: --long octave:,matlab:,fftw:,arch:,gcc-arch:,matlab-version:,without-threads: -n 'nfft-build-dll.sh' -- "$@"`
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
        -s|--without-threads)
            OMP_ONLY=0; shift 2 ;;
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
    GCCARCH=haswell;
  fi
else
  echo "Unknown architecture!" ; exit 1 ;
fi

# Install required packages
pacman -S --needed autoconf perl libtool automake mingw-w64-$ARCHNAME-gcc make mingw-w64-$ARCHNAME-cunit tar zip unzip wget dos2unix rsync p7zip

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
  cd "$FFTWDIR"/build-threads-static
  ../configure $FFTWFLAGS --host=$ARCHNAME-w64-mingw32.static --disable-fortran --enable-static --disable-shared --enable-threads --with-combined-threads
  make -j4
  touch "$FFTWDIR/build-success"
fi
if [ "$OMP_ONLY" = "0" ] && [ ! -f "$FFTWDIR/build/build-success" ]; then
  cd "$FFTWDIR"/build
  "$FFTWDIR"/configure $FFTWFLAGS --with-windows-f77-mangling --enable-shared
  make -j4
  cd ../build-static
  "$FFTWDIR"/configure $FFTWFLAGS --host=$ARCHNAME-w64-mingw32.static --enable-static --disable-shared
  make -j4
  touch "$FFTWDIR/build/build-success"
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

# Get Julia
if [ ! -f julia/bin/julia.exe ]; then
  rm -f -r julia
  mkdir julia
  cd julia
  wget https://julialang-s3.julialang.org/bin/winnt/x$ARCH/1.3/julia-1.3.1-win$ARCH.exe
  7z x julia-*.exe
  7z x julia-installer.exe
  rm julia-*.exe
fi

# Build NFFT
READMECONTENT="
$(sed -e '/^\[!/d' -e '/Directory structure/Q' $NFFTDIR/README) 
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

for OMPYN in $OMP_ONLY  1
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
"$NFFTDIR/configure" --enable-all $OMPFLAG --with-fftw3-libdir="$FFTWBUILDDIR"/.libs --with-fftw3-includedir="$FFTWDIR"/api --with-octave="$OCTAVEDIR" --with-gcc-arch=$GCCARCH --disable-static --enable-shared --disable-examples --enable-exhaustive-unit-tests
make
make check

# Create DLL release
NFFTVERSION=$( grep 'Version: ' nfft3.pc | cut -c10-)
DLLDIR=nfft-"$NFFTVERSION"-dll$ARCH$OMPSUFFIX

gcc -shared  -Wl,--whole-archive kernel/.libs/libkernel$THREADSSUFFIX.a -lfftw3 -lm -L"$FFTWBUILDDIR-static/.libs" -Wl,--no-whole-archive -O3 -malign-double   -o .libs/libnfft3$THREADSSUFFIX-$SOVERSION.dll -Wl,-Bstatic -lwinpthread $OMPLIBS -Wl,--output-def,.libs/libnfft3$THREADSSUFFIX-$SOVERSION.def

mkdir "$DLLDIR"
cp ".libs/libnfft3$THREADSSUFFIX-$SOVERSION.dll" "$DLLDIR/libnfft3$THREADSSUFFIX-$SOVERSION.dll"
cp ".libs/libnfft3$THREADSSUFFIX-$SOVERSION.def" "$DLLDIR/libnfft3$THREADSSUFFIX-$SOVERSION.def"
cp "$NFFTDIR"/include/nfft3.h "$DLLDIR"/nfft3.h
cp "$NFFTDIR"/include/nfft3mp.h "$DLLDIR"/nfft3mp.h
cp "$FFTWDIR"/api/fftw3.h "$DLLDIR"/fftw3.h
cp examples/nfft/simple_test.c "$DLLDIR"/simple_test.c
echo 'This archive contains the NFFT' $NFFTVERSION 'library and the associated header files.
The NFFT library was compiled with double precision support for' $ARCH'-bit Windows
using GCC' $GCCVERSION $ARCHNAME'-w64-mingw32 with -march='$GCCARCH 'and FFTW' $FFTWVERSION'.

In order to link the .dll file from Visual C++, you should run
    lib /def:libnfft3'$THREADSSUFFIX'-'$SOVERSION'.def

As a small example, you can compile the NFFT simple test with the following command

	gcc -O3 simple_test.c -o simple_test.exe -I. -L. libnfft3'$THREADSSUFFIX'-'$SOVERSION'.dll
'"$BINARIES_ARCH_README""$READMECONTENT""$FFTWREADME" > "$DLLDIR"/readme-windows.txt
unix2dos "$DLLDIR"/readme-windows.txt
cp "$NFFTDIR"/COPYING "$DLLDIR"/COPYING
rm -f "$HOMEDIR/$DLLDIR".zip
7z a -r "$HOMEDIR/$DLLDIR".zip "$DLLDIR"


# Julia interface
cd julia
for LIB in nf*t fastsum
do
  cd "$LIB"
  if [ "$LIB" = "fastsum" ]; then FASTSUM_LIBS="../../applications/fastsum/.libs/libfastsum$THREADSSUFFIX.a ../../applications/fastsum/.libs/libkernels.a"; else FASTSUM_LIBS=""; fi
  gcc -shared  .libs/lib"$LIB"julia.o  -Wl,--whole-archive ../../.libs/libnfft3_julia.a $FASTSUM_LIBS -Wl,--no-whole-archive -L"$FFTWBUILDDIR-static/.libs"  -O3 -malign-double -ffast-math -march=$GCCARCH  -o .libs/lib"$LIB"julia.dll -Wl,--enable-auto-image-base -Xlinker --out-implib -Xlinker .libs/lib"$LIB"julia.dll.a -static-libgcc -Wl,-Bstatic -lwinpthread -lfftw3 $OMPLIBS
  cp .libs/lib"$LIB"julia.dll lib"$LIB"julia.dll
  cd ..
done

cd "$NFFTBUILDDIR"
JULIADIR=nfft-"$NFFTVERSION"-julia-w$ARCH$OMPSUFFIX
mkdir "$JULIADIR"
cp "$NFFTDIR"/COPYING "$JULIADIR"/COPYING
rsync -rLt --exclude='Makefile*' --exclude='doxygen*' --exclude='*.c.in' --exclude='*.c' --exclude='*.h' "$NFFTDIR/julia/" "$JULIADIR"
rsync -rLt --exclude='Makefile*' --exclude='.deps' --exclude='.libs' --exclude='*.la' --exclude='*.lo' --exclude='*.o' --exclude='*.c' 'julia/' "$JULIADIR"
for DIR in $JULIADIR/nf*t $JULIADIR/fastsum; do cd $DIR; for NAME in simple_test*.jl; do PATH=/c/Windows/System32 $HOMEDIR/julia/bin/julia.exe "$NAME"; done; cd "$NFFTBUILDDIR"; done;

echo 'This archive contains the NFFT' $NFFTVERSION 'Julia interface.
The NFFT library was compiled with double precision support for' $ARCH'-bit Windows
using GCC' $GCCVERSION $ARCHNAME'-w64-mingw32 with -march='$GCCARCH 'and FFTW' $FFTWVERSION'.
'"$BINARIES_ARCH_README""$READMECONTENT""$FFTWREADME" > "$JULIADIR"/readme-julia.txt
unix2dos "$JULIADIR"/readme-julia.txt
rm -f "$HOMEDIR/$JULIADIR".zip
7z a -r "$HOMEDIR/$JULIADIR".zip "$JULIADIR"


# Compile with Matlab
if [ -n "$MATLABDIR" ]; then
  if [ "$MATLABVERSION" == "" ]; then
    "$MATLABDIR"/bin/matlab.exe -wait -nodesktop -nosplash -r "fid=fopen('matlab_version.txt','wt'); fprintf(fid,'MATLAB_VERSION=%s\n', version); exit;" 
    MATLABVERSION="`grep MATLAB_VERSION matlab_version.txt | sed 's/.*(//' | sed 's/)//'`"
  fi
  cd "$NFFTBUILDDIR"
  "$NFFTDIR/configure" --enable-all $OMPFLAG --with-fftw3-libdir="$FFTWBUILDDIR"/.libs --with-fftw3-includedir="$FFTWDIR"/api --with-matlab="$MATLABDIR" "$MATLABARCHFLAG" --with-gcc-arch=$GCCARCH --disable-static --enable-shared --disable-applications --disable-examples --enable-exhaustive-unit-tests
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
mkdir "$MEXDIR"
rsync -rLt --exclude='Makefile*' --exclude='doxygen*' --exclude='*.c.in' --exclude='*.c' --exclude='*.h' "$NFFTDIR/matlab/" "$MEXDIR"
rsync -rLt --exclude='Makefile*' --exclude='.deps' --exclude='.libs' --exclude='*.la' --exclude='*.lo' --exclude='*.o' --exclude='*.c' "matlab/" "$MEXDIR"
for SUBDIR in nfft nfsft nfsoft nnfft fastsum nfct nfst infft1d fpt ; do
  cd "$MEXDIR/$SUBDIR"
  if [ -f simple_test.m ] ; then
  for TESTFILE in *test*.m
    do
    if [ "$SUBDIR" != "infft1d" ] ; then
      "$OCTAVEDIR"/bin/octave-cli.exe --no-window-system --eval="run('$TESTFILE')"
    fi
    if [ -f "$MATLABDIR"/bin/matlab.exe ] ; then
      PATH=/c/Windows/System32 "$MATLABDIR"/bin/matlab.exe -wait -nodesktop -nosplash -r "run('$TESTFILE'); exit"
    fi
  done
  fi
  cd "$NFFTBUILDDIR"
done

cd "$NFFTBUILDDIR"
cp "$NFFTDIR"/COPYING "$MEXDIR"/COPYING
echo 'This archive contains the Matlab and Octave interface of NFFT '$NFFTVERSION'
compiled for '$ARCH'-bit Windows using GCC' $GCCVERSION $ARCHNAME'-w64-mingw32
with -march='$GCCARCH$' and Matlab '$MATLABVERSION' and Octave '$OCTAVEVERSION'.
'"$BINARIES_ARCH_README""$READMECONTENT""$FFTWREADME" > "$MEXDIR"/readme-matlab.txt
unix2dos "$MEXDIR"/readme-matlab.txt
rm -f "$HOMEDIR/$MEXDIR".zip
7z a -r "$HOMEDIR/$MEXDIR".zip "$MEXDIR"


# Build applications
cd "$NFFTBUILDDIR"
CPPFLAGS="--static" LDFLAGS="--static" "$NFFTDIR/configure" --enable-all $OMPFLAG --with-fftw3-libdir="$FFTWBUILDDIR-static"/.libs --with-fftw3-includedir="$FFTWDIR"/api --with-gcc-arch=$GCCARCH --enable-static --disable-shared --enable-exhaustive-unit-tests
make all check

APPSDIR=nfft-"$NFFTVERSION"-applications-w$ARCH$OMPSUFFIX
mkdir "$APPSDIR"
cp "$NFFTDIR"/COPYING "$APPSDIR"/COPYING
rsync -rLt --exclude='Makefile*' --exclude='.deps' --exclude='.libs' --exclude='*.la' --exclude='*.lo' --exclude='*.o' 'applications/' "$APPSDIR"
rsync -rLt --exclude='Makefile*' --exclude='doxygen*' --exclude='*.c.in' "$NFFTDIR/applications/" "$APPSDIR"

echo 'This archive contains the NFFT' $NFFTVERSION 'applications.
The NFFT library was compiled with double precision support for' $ARCH'-bit Windows
using GCC' $GCCVERSION $ARCHNAME'-w64-mingw32 with march='$GCCARCH 'and FFTW' $FFTWVERSION'.
'"$BINARIES_ARCH_README""$READMECONTENT""$FFTWREADME" > "$APPSDIR"/readme.txt
unix2dos "$APPSDIR"/readme.txt
rm -f "$HOMEDIR/$APPSDIR".zip
7z a -r "$HOMEDIR/$APPSDIR".zip "$APPSDIR"

done
