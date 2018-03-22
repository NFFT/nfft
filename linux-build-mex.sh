#!/bin/bash

# This script builds Octave / Matlab interfaces for Linux.
# A Matlab installation must be specified in order to build the Matlab interface.
# The paths should not contain spaces!
# 
# Example call:
# ./linux-build-mex.sh --matlab=/path/to/matlab
# 

# Any subsequent commands which fail will cause the shell script to exit immediately
set -ex

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

# Build NFFT
READMECONTENT="
$(sed '/Directory structure/Q' $NFFTDIR/README) 
"

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
else
  NFFTBUILDDIR="$HOMEDIR/build"
  OMPFLAG=""
  OMPLIBS=""
  THREADSSUFFIX=""
  OMPSUFFIX=""
fi

rm -f -r "$NFFTBUILDDIR"
mkdir "$NFFTBUILDDIR"
cd "$NFFTBUILDDIR"

"$NFFTDIR/configure" --enable-all $OMPFLAG --with-octave="$OCTAVEDIR" --with-gcc-arch=core2 --disable-static --enable-shared
make
make check

NFFTVERSION=$( grep 'Version: ' nfft3.pc | cut -c10-)
DIR=nfft-$NFFTVERSION-mexa64$OMPSUFFIX

# Create Matlab/Octave release
for SUBDIR in nfft nfsft nfsoft nnfft fastsum nfct nfst infft1d nfsft/@f_hat
  do
  mkdir -p "$DIR"/$SUBDIR
  cp -f -L -r matlab/$SUBDIR/*.mex* "$DIR"/$SUBDIR/ || true
  cp -f -L -r "$NFFTDIR"/matlab/$SUBDIR/README "$DIR"/$SUBDIR/ || true
  cp -r "$NFFTDIR"/matlab/$SUBDIR/*.m "$DIR"/$SUBDIR/
done

# Compile with Matlab
if [ -n "$MATLABDIR" ]; then
  MATLABVERSION=`"$MATLABDIR"/bin/matlab -nodisplay -r "fprintf('MATLAB_VERSION=%s\n', version); exit;" | grep MATLAB_VERSION | sed 's/.*(//' | sed 's/)//'`
  cd "$NFFTBUILDDIR"
  make clean
  "$NFFTDIR/configure" --enable-all $OMPFLAG --with-matlab="$MATLABDIR" --with-gcc-arch=core2 --disable-static --enable-shared
  make
  make check
fi

for SUBDIR in nfft nfsft nfsoft nnfft fastsum nfct nfst
  do
  cp -f -L -r matlab/$SUBDIR/*.mex* "$DIR"/$SUBDIR/
done

cd "$NFFTBUILDDIR"
cp "$NFFTDIR"/COPYING "$DIR"/COPYING
if [ -n "$MATLABDIR" ]; then
echo 'This archive contains the Matlab and Octave interface of NFFT '$NFFTVERSION' compiled for
64-bit Linux using GCC '$GCCVERSION' and Matlab '$MATLABVERSION' and Octave '$OCTAVEVERSION'.
' "$READMECONTENT" > "$DIR"/readme-matlab.txt
else
echo 'This archive contains the Matlab and Octave interface of NFFT '$NFFTVERSION' compiled for
64-bit Linux using GCC '$GCCVERSION' and Octave '$OCTAVEVERSION'.
' "$READMECONTENT" > "$DIR"/readme-matlab.txt
fi
tar czf ../"$DIR".tar.gz "$DIR"

done
