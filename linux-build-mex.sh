#!/bin/bash

# This script builds Octave / Matlab interfaces for Linux.
# A Matlab installation must be specified in order to build the Matlab interface.
# The paths should not contain spaces!
# 
# The script is known to work on Ubuntu 18.10. At least the following packages
# are required:
# fftw3 libfftw3-dev libcunit1-dev make gcc octave liboctave-dev
#
# For running ./bootstrap.sh, the following packages are required:
# autoconf automake libtool
#
#
# Example call:
# ./linux-build-mex.sh --matlab=/path/to/matlab
# 
#
# For compiling with recent Octave version, one can use flatpak package org.octave.Octave.
# Then, the following RSYNC variable should be used.
RSYNC="flatpak-spawn --host rsync"
# Otherwise, the RSYNC variable should be set to
#RSYNC=rsync
#
# Possible flatpak call:
#   flatpak --command="bash" run org.octave.Octave
#   bash linux-build-mex.sh -o /app

# Any subsequent commands which fail will cause the shell script to exit immediately
set -ex

FFTWVERSION=3.3.8
GCCARCH=core2

# default values (to be overwritten if respective parameters are set)
OCTAVEDIR=/usr

# read the options
TEMP=`getopt -o o:m:f: --long octave:,matlab:,fftw: -n 'linux-build-mex.sh' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -o|--octave)
            case "$2" in
                "")  shift 2 ;;
                *) OCTAVEDIR=$2; shift 2 ;;
            esac ;;
        -m|--matlab)
            case "$2" in
                "")  shift 2 ;;
                *) MATLABDIR=$2; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

NFFTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
HOMEDIR="$NFFTDIR"/linux-build-mex
mkdir -p "$HOMEDIR"
cd "$HOMEDIR"
GCCVERSION=$(gcc -dumpversion)
OCTAVEVERSION=`"$OCTAVEDIR"/bin/octave-cli --eval "fprintf('OCTAVE_VERSION=%s\n', version); exit;" | grep OCTAVE_VERSION | sed 's/OCTAVE_VERSION=//'`


FFTWDIR=$HOMEDIR/fftw-$FFTWVERSION
# Build FFTW
if [ ! -f "$FFTWDIR/build-success" ]; then
  rm -rf "$FFTWDIR"
  curl "http://fftw.org/fftw-$FFTWVERSION.tar.gz" --output "fftw-$FFTWVERSION.tar.gz"
  tar -zxf "fftw-$FFTWVERSION.tar.gz"
  rm "fftw-$FFTWVERSION.tar.gz"
  cd "$FFTWDIR"

  mkdir build
  cd build
  ../configure --enable-static --enable-shared --enable-threads --with-pic --enable-sse2 --enable-avx --enable-avx2 --disable-fortran
  make -j4
  touch "$FFTWDIR/build-success"
  cd "$HOMEDIR"
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
make distclean || true

for OMPYN in 0 1
do
if [ $OMPYN = 1 ]; then
  NFFTBUILDDIR="$HOMEDIR/build-openmp"
  OMPFLAG="--enable-openmp"
  OMPLIBS="-fopenmp -static-libgcc"
  THREADSSUFFIX="_threads"
  OMPSUFFIX="-openmp"
  FFTWLIBSTATIC="$FFTWDIR/build/threads/.libs/libfftw3_threads.a $FFTWDIR/build/.libs/libfftw3.a"
else
  NFFTBUILDDIR="$HOMEDIR/build"
  OMPFLAG=""
  OMPLIBS=""
  THREADSSUFFIX=""
  OMPSUFFIX=""
  FFTWLIBSTATIC="$FFTWDIR/build/.libs/libfftw3.a"
fi

rm -f -r "$NFFTBUILDDIR"
mkdir "$NFFTBUILDDIR"
cd "$NFFTBUILDDIR"

LDFLAGS="-L$FFTWDIR/build/threads/.libs -L$FFTWDIR/build/.libs"
CPPFLAGS="-I$FFTWDIR/api"
"$NFFTDIR/configure" --enable-all $OMPFLAG --with-octave="$OCTAVEDIR" --with-gcc-arch="$GCCARCH" --disable-static --enable-shared
make
make check

NFFTVERSION=$( grep 'Version: ' nfft3.pc | cut -c10-)

# Create archive for Julia interface
cd julia
for LIB in nf*t
do
  cd "$LIB"
  gcc -shared  -fPIC -DPIC  .libs/lib"$LIB"julia.o  -Wl,--whole-archive ../../.libs/libnfft3_julia.a $FFTWLIBSTATIC -Wl,--no-whole-archive  $OMPLIBS -lm  -O3 -malign-double -march="$GCCARCH"   -Wl,-soname -Wl,lib"$LIB"julia.so -o .libs/lib"$LIB"julia.so
  cd ..
done
cd "$NFFTBUILDDIR"

ARCH=$(uname -m)
JULIADIR=nfft-"$NFFTVERSION"-julia-linux_$ARCH$OMPSUFFIX
mkdir "$JULIADIR"
$RSYNC -rLt --exclude='Makefile*' --exclude='doxygen*' --exclude='*.c.in' --exclude='*.c' --exclude='*.h' "$NFFTDIR/julia/" "$JULIADIR"
$RSYNC -rLt --exclude='Makefile*' --exclude='.deps' --exclude='.libs' --exclude='*.la' --exclude='*.lo' --exclude='*.o' --exclude='*.c' 'julia/' "$JULIADIR"
echo 'This archive contains the Julia interface of NFFT '$NFFTVERSION'
compiled for '$ARCH' Linux using GCC '$GCCVERSION' and FFTW '$FFTWVERSION'.
' "$READMECONTENT" "$FFTWREADME" > "$JULIADIR"/readme.txt
tar czf ../"$JULIADIR".tar.gz --owner=0 --group=0 "$JULIADIR"
# End of Julia interface


# Create Matlab/Octave release
DIR=nfft-$NFFTVERSION-mexa64$OMPSUFFIX
mkdir $DIR
$RSYNC -rLt --exclude='Makefile*' --exclude='doxygen*' --exclude='*.c.in' --exclude='*.c' --exclude='*.h' "$NFFTDIR/matlab/" "$DIR"
$RSYNC -rLt --exclude='Makefile*' --exclude='.deps' --exclude='.libs' --exclude='*.la' --exclude='*.lo' --exclude='*.o' --exclude='*.c' "matlab/" "$DIR"

# Compile with Matlab
if [ -n "$MATLABDIR" ]; then
  MATLABVERSION=`"$MATLABDIR"/bin/matlab -nodisplay -r "fprintf('MATLAB_VERSION=%s\n', version); exit;" | grep MATLAB_VERSION | sed 's/.*(//' | sed 's/)//'`
  cd "$NFFTBUILDDIR"
  make clean
  "$NFFTDIR/configure" --enable-all $OMPFLAG --with-matlab="$MATLABDIR" --with-gcc-arch="$GCCARCH" --disable-static --enable-shared
  make
  make check
fi

for SUBDIR in nfft nfsft nfsoft nnfft fastsum nfct nfst fpt
  do
  cp -f -L -r matlab/$SUBDIR/*.mex* "$DIR"/$SUBDIR/
done

cd "$NFFTBUILDDIR"
cp "$NFFTDIR"/COPYING "$DIR"/COPYING
if [ -n "$MATLABDIR" ]; then
echo 'This archive contains the Matlab and Octave interface of NFFT '$NFFTVERSION'
compiled for '$ARCH' Linux using GCC '$GCCVERSION' and Matlab '$MATLABVERSION' and Octave '$OCTAVEVERSION'.
' "$READMECONTENT" > "$DIR"/readme-matlab.txt
else
echo 'This archive contains the Octave interface of NFFT '$NFFTVERSION' compiled for
64-bit Linux using GCC '$GCCVERSION' and Octave '$OCTAVEVERSION'.
' "$READMECONTENT" > "$DIR"/readme-matlab.txt
fi
tar czf ../"$DIR".tar.gz --owner=0 --group=0 "$DIR"

done
