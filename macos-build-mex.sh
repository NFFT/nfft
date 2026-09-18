#!/bin/bash

# This script builds Octave / Matlab interfaces for macOS.
# A Matlab installation must be specified in order to build the
# Matlab interface. The paths should not contain spaces!
#
# The script is known to work on macOS 11 Big Sur with Homebrew.
#
# At least the following packages are required:
# fftw automake autoconf gcc@12 gnu-sed cunit gnu-sed julia octave
#
# Example call:
# ./macos-build-mex.sh --matlab=/path/to/matlab
# 

# Any subsequent commands which fail will cause the shell script to exit immediately
set -ex

# Determine Homebrew prefix
HOMEBREW_PREFIX=$(brew --prefix)

# Pin the macOS deployment target to the running OS version.
# Homebrew GCC derives the deployment target from the Darwin kernel version.
# Since the macOS 15 -> 26 renumbering (Darwin 25), GCC 12 mis-derives this as
# "16.0", stamping objects with minos 16.0 while the SDK and the prebuilt
# Homebrew FFTW/Octave archives are stamped 26.0. The linker then rejects the
# mismatch ("object file was built for newer macOS version than being linked").
# Forcing the real major version keeps objects and libraries in sync, and is a
# no-op on macOS 15 and earlier where GCC already derives it correctly.
export MACOSX_DEPLOYMENT_TARGET="$(sw_vers -productVersion | cut -d. -f1).0"

# Architecture-dependent settings, auto-detected from the host CPU.
# The x86_64 path keeps the original Haswell baseline; the arm64 path mirrors it
# with apple-m1, the floor for all Apple Silicon Macs. Note that -march/-malign-double
# are x86-only: on arm64 GCC they are hard errors, so the arm64 branch uses -mcpu and
# drops -malign-double. MEXARCH/MATLABBIN follow MathWorks' per-arch naming.
HOSTARCH=$(uname -m)
if [ "$HOSTARCH" = "arm64" ]; then
  GCCARCH=apple-m1
  ARCHFLAGS="-mcpu=$GCCARCH"
  MEXARCH=mexmaca64
  MATLABBIN=maca64
  BINARIES_ARCH_NOTE='
Please note that the binaries were compiled with gcc flag -mcpu=apple-m1 and
therefore require an Apple Silicon (M1 or newer) Mac.
'
else
  GCCARCH=haswell
  ARCHFLAGS="-march=$GCCARCH -malign-double"
  MEXARCH=mexmaci64
  MATLABBIN=maci64
  BINARIES_ARCH_NOTE='
Please note that since the binaries were compiled with gcc flag -march=haswell,
they may not work on older CPUs (below Intel i3/i5/i7-4xxx or
AMD Excavator/4th gen Bulldozer) as well as on some Intel Atom/Pentium CPUs.
'
fi

FFTWDIR="$HOMEBREW_PREFIX"
GCC="gcc-12"

# default values (to be overwritten if respective parameters are set)
OCTAVEDIR="$HOMEBREW_PREFIX"
MATLABDIR=/Applications/MATLAB_R2021b.app
# read the options (pure POSIX)
while [ $# -gt 0 ]; do
    echo "opt: $1"
    case "$1" in
        -o|--octave)
            if [ -n "$2" ]; then
                OCTAVEDIR="$2"
                shift 2
            else
                shift
            fi
            ;;
        -o=*|--octave=*)
            OCTAVEDIR="${1#*=}"
            shift
            ;;
        -m|--matlab)
            if [ -n "$2" ]; then
                MATLABDIR="$2"
                shift 2
            else
                shift
            fi
            ;;
        -m=*|--matlab=*)
            MATLABDIR="${1#*=}"
            shift
            ;;
        -v|--matlab-version)
            if [ -n "$2" ]; then
                MATLABVERSION="$2"
                shift 2
            else
                shift
            fi
            ;;
        -v=*|--matlab-version=*)
            MATLABVERSION="${1#*=}"
            shift
            ;;
        -f|--fftw)
            if [ -n "$2" ]; then
                FFTWDIR="$2"
                shift 2
            else
                shift
            fi
            ;;
        -f=*|--fftw=*)
            FFTWDIR="${1#*=}"
            shift
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Internal error!"
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

NFFTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
HOMEDIR="$NFFTDIR"/macos-build-mex
mkdir -p "$HOMEDIR"
cd "$HOMEDIR"
GCCVERSION=`$GCC -dumpversion`
FFTWVERSION=`fftw-wisdom | grep fftw- | gsed 's/(fftw-//' | gsed 's/ fftw_wisdom.*//'`
OCTAVEVERSION=`"$OCTAVEDIR"/bin/octave-cli --eval "fprintf('OCTAVE_VERSION=%s\n', version); exit;" | grep OCTAVE_VERSION | gsed 's/OCTAVE_VERSION=//'`

# Build NFFT
READMECONTENT="
$(gsed -e '/^\[!/d' -e '/Directory structure/Q' $NFFTDIR/README.md) 
"
FFTWREADME='
FFTW
----
The compiled NFFT files contain parts of the FFTW library (http://www.fftw.org)
Copyright (c) 2003, 2007-14 Matteo Frigo
Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology'
BINARIES_ARCH_README="$BINARIES_ARCH_NOTE"

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
  FFTW_LINK_COMMAND="-Wl,-force_load,$FFTWDIR/lib/libfftw3_threads.a -Wl,-force_load,$FFTWDIR/lib/libfftw3.a"
else
  NFFTBUILDDIR="$HOMEDIR/build"
  OMPFLAG=""
  OMPLIBS="-static-libgcc"
  THREADSSUFFIX=""
  OMPSUFFIX=""
  FFTW_LINK_COMMAND="-Wl,-force_load,$FFTWDIR/lib/libfftw3.a"
fi

rm -f -r "$NFFTBUILDDIR"
mkdir "$NFFTBUILDDIR"
cd "$NFFTBUILDDIR"

CC=$GCC CPPFLAGS=-I"$FFTWDIR"/include LDFLAGS=-L"$FFTWDIR"/lib "$NFFTDIR/configure" --enable-all $OMPFLAG --with-octave="$OCTAVEDIR" --with-gcc-arch=$GCCARCH --disable-static --enable-shared --disable-examples --enable-applications
make
make check

NFFTVERSION=$( grep 'Version: ' nfft3.pc | cut -c10-)

# Create archive for Julia interface
cd julia
for LIB in nf*t
do
  cd "$LIB"
  $GCC -dynamiclib -o lib"$LIB"julia.dylib .libs/lib"$LIB"julia.o -Wl,-force_load,../../.libs/libnfft3_julia.a $FFTW_LINK_COMMAND -lm -O3 $ARCHFLAGS $OMPLIBS
  cd ..
done

cd fastsum
$GCC -dynamiclib -o libfastsumjulia.dylib .libs/libfastsumjulia.o -Wl,-force_load,../../.libs/libnfft3_julia.a $FFTW_LINK_COMMAND -Wl,-force_load,../../applications/fastsum/.libs/libfastsum$THREADSSUFFIX.a -Wl,-force_load,../../applications/fastsum/.libs/libkernels.a -lm -O3 $ARCHFLAGS $OMPLIBS
cd ..

cd "$NFFTBUILDDIR"

ARCH=$(uname -m)
JULIADIR=nfft-"$NFFTVERSION"-julia-macos_$ARCH$OMPSUFFIX
mkdir "$JULIADIR"
rsync -rLt --exclude='Makefile*' --exclude='doxygen*' --exclude='*.c.in' --exclude='*.c' --exclude='*.h' "$NFFTDIR/julia/" "$JULIADIR"
rsync -rLt --exclude='Makefile*' --exclude='.deps' --exclude='.libs' --exclude='*.la' --exclude='*.lo' --exclude='*.o' --exclude='*.c' 'julia/' "$JULIADIR"

for DIR in $JULIADIR/nf*t $JULIADIR/fastsum; do cd $DIR; for NAME in simple_test*.jl; do julia "$NAME"; done; cd "$NFFTBUILDDIR"; done;

echo 'This archive contains the NFFT' $NFFTVERSION 'Julia interface.
The NFFT library was compiled with double precision support for '$ARCH' macOS
using GCC '$GCCVERSION' with flags '$ARCHFLAGS' and FFTW '$FFTWVERSION'.
'"$BINARIES_ARCH_README""$READMECONTENT""$FFTWREADME" > "$JULIADIR"/readme-julia.txt
zip -9 -r ../"$JULIADIR".zip "$JULIADIR"
# End of Julia interface


# Create Matlab/Octave release
for LIB in nfft nfsft nfsoft nnfft fastsum nfct nfst fpt
do
  cd matlab/"$LIB"
  # The MEX entry points (mexFunction, mxGetPr, ...) live in Octave's liboctmex and
  # are resolved at load time from the host Octave process. The modern macOS linker no
  # longer tolerates undefined symbols in a -bundle by default, so make it explicit.
  $GCC -o .libs/lib"$LIB".mex -bundle -Wl,-undefined,dynamic_lookup  .libs/lib"$LIB"_la-"$LIB"mex.o -Wl,-force_load,../../.libs/libnfft3_matlab.a -Wl,-force_load,../../matlab/.libs/libmatlab.a -L"$OCTAVEDIR"/lib/octave/"$OCTAVEVERSION" $FFTW_LINK_COMMAND -lm -loctinterp -loctave -O3 $ARCHFLAGS $OMPLIBS
  cd ../..
done

DIR=nfft-$NFFTVERSION-$MEXARCH$OMPSUFFIX
mkdir $DIR
rsync -rLt --exclude='Makefile*' --exclude='doxygen*' --exclude='*.c.in' --exclude='*.c' --exclude='*.h' "$NFFTDIR/matlab/" "$DIR"
rsync -rLt --exclude='Makefile*' --exclude='.deps' --exclude='.libs' --exclude='*.la' --exclude='*.lo' --exclude='*.o' --exclude='*.c' "matlab/" "$DIR"

# Compile with Matlab
# Guard: the MEX must link against MATLAB libraries of the build architecture.
# MathWorks ships those under bin/<arch> (maca64 for Apple Silicon, maci64 for
# Intel). If the matching libdir is missing, skip the MATLAB build (and its 
# tests/readme) instead of failing, leaving the Octave interface intact.
if [ -n "$MATLABDIR" ] && [ ! -d "$MATLABDIR/bin/$MATLABBIN" ]; then
  echo "WARNING: $MATLABDIR has no bin/$MATLABBIN ($HOSTARCH) libraries;" \
       "skipping the MATLAB interface. Use an Apple Silicon MATLAB (R2023b+)" \
       "to build it, or pass --matlab to point at one." >&2
  MATLABDIR=""
fi
if [ -n "$MATLABDIR" ]; then
  if [ -z "$MATLABVERSION" ]; then
    MATLABVERSION=`"$MATLABDIR"/bin/matlab -wait -nodesktop -nosplash -r "fprintf('MATLAB_VERSION=%s\n', version); exit;" | grep MATLAB_VERSION | gsed 's/.*(//' | gsed 's/)//'`
  fi
  MATLABSTRING="and Matlab $MATLABVERSION "
  cd "$NFFTBUILDDIR"
  make clean
  CC=$GCC CPPFLAGS=-I"$FFTWDIR"/include LDFLAGS=-L"$FFTWDIR"/lib "$NFFTDIR/configure" --enable-all $OMPFLAG --with-matlab="$MATLABDIR" --with-gcc-arch=$GCCARCH --disable-static --enable-shared --disable-examples --disable-applications
  make
  make check
  for LIB in nfft nfsft nfsoft nnfft fastsum nfct nfst fpt
  do
    cd matlab/"$LIB"
    $GCC -o .libs/lib"$LIB".$MEXARCH -bundle  .libs/lib"$LIB"_la-"$LIB"mex.o   -Wl,-force_load,../../.libs/libnfft3_matlab.a -Wl,-force_load,../../matlab/.libs/libmatlab.a  -L"$MATLABDIR"/bin/$MATLABBIN -lm -lmwfftw3 -lmx -lmex -lmat -O3 $ARCHFLAGS $OMPLIBS
    cd ../..
  done
fi

for SUBDIR in nfft nfsft nfsoft nnfft fastsum nfct nfst fpt
  do
  cp -f -L matlab/$SUBDIR/*.mex* "$DIR"/$SUBDIR/
done

for SUBDIR in nfft nfsft nfsoft nnfft fastsum nfct nfst infft1d fpt ; do
  cd "$DIR/$SUBDIR"
  if [ -f simple_test.m ] ; then
  for TESTFILE in *test*.m
    do
    if [ "$SUBDIR" != "infft1d" ] ; then
      "$OCTAVEDIR"/bin/octave-cli --eval="run('$TESTFILE')"
    fi
     if [ -n "$MATLABDIR" ]; then
      "$MATLABDIR"/bin/matlab -wait -nodesktop -nosplash -r "run('$TESTFILE'); exit"
    fi
  done
  fi
  cd "$NFFTBUILDDIR"
done

cd "$NFFTBUILDDIR"
cp "$NFFTDIR"/COPYING "$DIR"/COPYING
if [ -n "$MATLABDIR" ]; then
echo 'This archive contains the Matlab and Octave interface of NFFT '$NFFTVERSION'
compiled for '$ARCH' macOS using GCC '$GCCVERSION' with flags '$ARCHFLAGS'
'$MATLABSTRING'and Octave '$OCTAVEVERSION'.
'"$BINARIES_ARCH_README""$READMECONTENT""$FFTWREADME" > "$DIR"/readme-matlab.txt
else
echo 'This archive contains the Octave interface of NFFT '$NFFTVERSION'
compiled for '$ARCH' macOS using GCC '$GCCVERSION' with flags '$ARCHFLAGS'
and Octave '$OCTAVEVERSION'.
'"$BINARIES_ARCH_README""$READMECONTENT""$FFTWREADME" > "$DIR"/readme-matlab.txt
fi

zip -9 -r ../"$DIR".zip "$DIR"

done
