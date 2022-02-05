#!/bin/bash

# This script builds Octave / Matlab interfaces for Linux.
# A Matlab installation must be specified in order to build the Matlab interface.
# The paths should not contain spaces!
# 
# The script is known to work on Ubuntu 18.10. At least the following packages
# are required:
# libcunit1-dev make gcc octave liboctave-dev
#
# For running ./bootstrap.sh, the following packages are required:
# autoconf automake libtool
#
#
# For compiling with recent Octave version, one can use flatpak package org.octave.Octave.
# Example call:
#   flatpak --command="bash" run org.octave.Octave
#   bash linux-build-mex.sh -o /app --matlab=/path/to/matlab
# 
#
# When using flatpak, the following RSYNC variable should be used.
RSYNC="flatpak-spawn --host rsync"
# Otherwise, the RSYNC variable should be set to
RSYNC=rsync
#

# Any subsequent commands which fail will cause the shell script to exit immediately
set -ex

# Write log file
exec > >(tee linux-build-mex.log)
exec 2>&1

FFTWVERSION=3.3.10
GCCVERSION=11.2.0
GCCARCH=haswell
BINARIES_ARCH_README='
Please note that since the binaries were compiled with gcc flag -march=haswell,
they may not work on older CPUs (below Intel i3/i5/i7-4xxx or
AMD Excavator/4th gen Bulldozer) as well as on some Intel Atom/Pentium CPUs.
'
MPFRVERSION=4.0.1
MPCVERSION=1.1.0

# default values (to be overwritten if respective parameters are set)
OCTAVEDIR=/usr

JULIA_BIN=julia/julia-1.6.2/bin/julia
JULIA_ARCHIVE=julia-1.6.2-linux-x86_64.tar.gz
JULIA_URL=https://julialang-s3.julialang.org/bin/linux/x64/1.6/$JULIA_ARCHIVE

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
HOMEDIR=$(pwd)/linux-build-mex
mkdir -p "$HOMEDIR"
cd "$HOMEDIR"
#GCCVERSION=$(gcc -dumpfullversion)
OCTAVEVERSION=`"$OCTAVEDIR"/bin/octave-cli --eval "fprintf('OCTAVE_VERSION=%s\n', version); exit;" | grep OCTAVE_VERSION | sed 's/OCTAVE_VERSION=//'`


MPFRBUILDDIR=$HOMEDIR/mpfr-$MPFRVERSION
MPFRINSTALLDIR=$HOMEDIR/mpfr-$MPFRVERSION-install
# Build MPFR for GCC
if [ ! -f "$MPFRINSTALLDIR/build-success" ]; then
  rm -rf "$MPFRBUILDDIR"
  rm -rf "$MPFRINSTALLDIR"
  curl "https://ftp.gnu.org/gnu/mpfr/mpfr-$MPFRVERSION.tar.gz" --output "mpfr-$MPFRVERSION.tar.gz"
  tar -zxf "mpfr-$MPFRVERSION.tar.gz"
  rm "mpfr-$MPFRVERSION.tar.gz"
  cd $MPFRBUILDDIR
  ./configure --prefix="$MPFRINSTALLDIR"
  make -j4
  make install
  touch "$MPFRINSTALLDIR/build-success"
  cd $HOMEDIR
fi

MPCBUILDDIR=$HOMEDIR/mpc-$MPCVERSION
MPCINSTALLDIR=$HOMEDIR/mpc-$MPCVERSION-install
# Build MPC for GCC
if [ ! -f "$MPCINSTALLDIR/build-success" ]; then
  rm -rf "$MPCBUILDDIR"
  rm -rf "$MPCINSTALLDIR"
  curl "https://ftp.gnu.org/gnu/mpc/mpc-$MPCVERSION.tar.gz" --output "mpc-$MPCVERSION.tar.gz"
  tar -zxf "mpc-$MPCVERSION.tar.gz"
  rm "mpc-$MPCVERSION.tar.gz"
  cd $MPCBUILDDIR
  ./configure --prefix="$MPCINSTALLDIR" --with-mpfr="$MPFRINSTALLDIR"
  make -j4
  make install
  touch "$MPCINSTALLDIR/build-success"
  cd $HOMEDIR
fi

export LD_LIBRARY_PATH="$MPCINSTALLDIR/lib:$MPFRINSTALLDIR/lib:$LD_LIBRARY_PATH"

GCCBUILDDIR="$HOMEDIR/gcc-$GCCVERSION"
GCCINSTALLDIR="$HOMEDIR/gcc-$GCCVERSION-install"
# Build GCC
if [ ! -f "$GCCINSTALLDIR/build-success" ]; then
  rm -rf "$GCCBUILDDIR"
  rm -rf "$GCCINSTALLDIR"
  curl "https://ftp.gnu.org/gnu/gcc/gcc-$GCCVERSION/gcc-$GCCVERSION.tar.gz" --output "gcc-$GCCVERSION.tar.gz"
  tar -zxf "gcc-$GCCVERSION.tar.gz"
  rm "gcc-$GCCVERSION.tar.gz"
  cd $GCCBUILDDIR
  CFLAGS=-fPIC CXXFLAGS=-fPIC LDFLAGS=-fPIC ./configure -enable-threads=posix --enable-checking=release --with-system-zlib --enable-__cxa_atexit --enable-languages=c,lto --disable-multilib --disable-nls --enable-bootstrap --prefix="$GCCINSTALLDIR" --with-mpc="$MPCINSTALLDIR" --with-mpfr="$MPFRINSTALLDIR" --program-suffix="-$GCCVERSION"
  make -j4
  make install
  touch "$GCCINSTALLDIR/build-success"
  cd $HOMEDIR
fi


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
  CC="$GCCINSTALLDIR/bin/gcc-$GCCVERSION" ../configure --enable-static --enable-shared --enable-threads --with-pic --enable-sse2 --enable-avx --enable-avx2 --disable-fortran
  make -j4
  touch "$FFTWDIR/build-success"
  cd "$HOMEDIR"
fi


# Get Julia
if [ ! -f $JULIA_BIN ]; then
  rm -f -r julia
  mkdir julia
  cd julia
  curl "$JULIA_URL" --output "$JULIA_ARCHIVE"
  tar xzf $JULIA_ARCHIVE
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

for OMPYN in 1
do
if [ $OMPYN = 1 ]; then
  NFFTBUILDDIR="$HOMEDIR/build-openmp"
  OMPFLAG="--enable-openmp"
  OMPLIBS="-fopenmp -static-libgcc"
  THREADSSUFFIX="_threads"
  OMPSUFFIX="-openmp"
  FFTWLIBSTATIC="$FFTWDIR/build/threads/.libs/libfftw3_threads.a -pthread $FFTWDIR/build/.libs/libfftw3.a -lm"
  GOMPLIBSTATIC="$GCCINSTALLDIR/lib64/libgomp.a"
else
  NFFTBUILDDIR="$HOMEDIR/build"
  OMPFLAG=""
  OMPLIBS=""
  THREADSSUFFIX=""
  OMPSUFFIX=""
  FFTWLIBSTATIC="$FFTWDIR/build/.libs/libfftw3.a -lm"
  GOMPLIBSTATIC=""
fi

rm -f -r "$NFFTBUILDDIR"
mkdir "$NFFTBUILDDIR"
cd "$NFFTBUILDDIR"

LDFLAGS="-L$FFTWDIR/build/threads/.libs -L$FFTWDIR/build/.libs"
CPPFLAGS="-I$FFTWDIR/api"
CC="$GCCINSTALLDIR/bin/gcc-$GCCVERSION" "$NFFTDIR/configure" --enable-all $OMPFLAG --with-octave="$OCTAVEDIR" --with-gcc-arch="$GCCARCH" --disable-static --enable-shared
make
make check

NFFTVERSION=$( grep 'Version: ' nfft3.pc | cut -c10-)

# Create archive for Julia interface
cd julia
for LIB in nf*t
do
  cd "$LIB"
  "$GCCINSTALLDIR/bin/gcc-$GCCVERSION" -shared  -fPIC -DPIC  .libs/lib"$LIB"julia.o  -Wl,--whole-archive ../../.libs/libnfft3_julia.a $FFTWLIBSTATIC $GOMPLIBSTATIC -Wl,--no-whole-archive -O3 -malign-double -march="$GCCARCH" -Wl,-soname -Wl,lib"$LIB"julia.so -o .libs/lib"$LIB"julia.so
  cd ..
done
for LIB in fastsum
do
  cd "$LIB"
  "$GCCINSTALLDIR/bin/gcc-$GCCVERSION" -shared  -fPIC -DPIC  .libs/lib"$LIB"julia.o  -Wl,--whole-archive ../../applications/fastsum/.libs/libfastsum$THREADSSUFFIX.a ../../applications/fastsum/.libs/libkernels.a ../../.libs/libnfft3_julia.a $FFTWLIBSTATIC $GOMPLIBSTATIC -Wl,--no-whole-archive -O3 -malign-double -march="$GCCARCH" -Wl,-soname -Wl,lib"$LIB"julia.so -o .libs/lib"$LIB"julia.so
  cd ..
done
cd "$NFFTBUILDDIR"

ARCH=$(uname -m)
JULIADIR=nfft-"$NFFTVERSION"-julia-linux_$ARCH$OMPSUFFIX
mkdir "$JULIADIR"
$RSYNC -rLt --exclude='Makefile*' --exclude='doxygen*' --exclude='*.c.in' --exclude='*.c' --exclude='*.h' --exclude='*.so' "$NFFTDIR/julia/" "$JULIADIR"
$RSYNC -rLt --exclude='Makefile*' --exclude='.deps' --exclude='.libs' --exclude='*.la' --exclude='*.lo' --exclude='*.o' --exclude='*.c' 'julia/' "$JULIADIR"
for DIR in $JULIADIR/nf*t $JULIADIR/fastsum; do cd $DIR; for NAME in simple_test*.jl; do $HOMEDIR/$JULIA_BIN "$NAME"; done; cd "$NFFTBUILDDIR"; done;

echo 'This archive contains the NFFT' $NFFTVERSION 'Julia interface.
The NFFT library was compiled with double precision support for '$ARCH' Linux
using GCC '$GCCVERSION' with -march='$GCCARCH' and FFTW '$FFTWVERSION'.
'"$BINARIES_ARCH_README""$READMECONTENT""$FFTWREADME" > "$JULIADIR"/readme-julia.txt
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
  CC="$GCCINSTALLDIR/bin/gcc-$GCCVERSION" "$NFFTDIR/configure" --enable-all $OMPFLAG --with-matlab="$MATLABDIR" --with-gcc-arch="$GCCARCH" --disable-static --enable-shared --enable-exhaustive-unit-tests
  make
  if [ $OMPYN = 1 ]; then
    cd matlab
    for SUBDIR in nfft nfsft nfsoft nnfft fastsum nfct nfst fpt
    do
      cd "$SUBDIR"
      "$GCCINSTALLDIR/bin/gcc-$GCCVERSION" -shared  -fPIC -DPIC  .libs/lib"$SUBDIR"_la-"$SUBDIR"mex.o  -Wl,--whole-archive ../../.libs/libnfft3_matlab.a ../../matlab/.libs/libmatlab.a $GOMPLIBSTATIC -Wl,--no-whole-archive  -L$MATLABDIR/bin/glnxa64 -l:libmwfftw3.so.3 -lm -lmx -lmex -lmat -O3 -malign-double -march="$GCCARCH" -Wl,-soname -Wl,lib$SUBDIR.mexa64 -o .libs/lib$SUBDIR.mexa64
      cd ..
    done
    cd "$NFFTBUILDDIR"
  fi
  make check
fi

for SUBDIR in nfft nfsft nfsoft nnfft fastsum nfct nfst fpt
do
  cp -f -L -r matlab/$SUBDIR/*.mex* "$DIR"/$SUBDIR/
done

for SUBDIR in nfft nfsft nfsoft nnfft fastsum nfct nfst infft1d fpt ; do
  cd "$DIR/$SUBDIR"
  if [ -f simple_test.m ] ; then
  for TESTFILE in *test*.m
    do
    if [ "$SUBDIR" != "infft1d" ] ; then
      "$OCTAVEDIR"/bin/octave-cli --no-window-system --eval="run('$TESTFILE')"
    fi
     if [ -n "$MATLABDIR" ]; then
      "$MATLABDIR"/bin/matlab -nodisplay -r "run('$TESTFILE'); exit"
    fi
  done
  fi
  cd "$NFFTBUILDDIR"
done


cd "$NFFTBUILDDIR"
cp "$NFFTDIR"/COPYING "$DIR"/COPYING
if [ -n "$MATLABDIR" ]; then
echo 'This archive contains the Matlab and Octave interface of NFFT '$NFFTVERSION'
compiled for '$ARCH' Linux using GCC '$GCCVERSION' with -march='$GCCARCH'
and Matlab '$MATLABVERSION' and Octave '$OCTAVEVERSION'.
'"$BINARIES_ARCH_README""$READMECONTENT""$FFTWREADME" > "$DIR"/readme-matlab.txt
else
echo 'This archive contains the Octave interface of NFFT '$NFFTVERSION'
compiled for '$ARCH' Linux using GCC '$GCCVERSION' with -march='$GCCARCH'
and Octave '$OCTAVEVERSION'.
'"$BINARIES_ARCH_README""$READMECONTENT""$FFTWREADME" > "$DIR"/readme-matlab.txt
fi
tar czf ../"$DIR".tar.gz --owner=0 --group=0 "$DIR"

done
