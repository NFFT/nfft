#!/bin/bash

# This script builds Octave / Matlab interfaces for macOS.
# A Matlab installation must be specified in order to build the
# Matlab interface. The paths should not contain spaces!
#
# The script is known to work on macOS 10.5 Catalina with Homebrew.
#
# At least the following packages are required:
# octave gnu-sed cunit
#
# 
# Example call:
# ./macos-build-mex.sh --matlab=/path/to/matlab
# 

# Any subsequent commands which fail will cause the shell script to exit immediately
set -ex

GCCARCH=haswell
FFTWDIR=/usr/local
GCC="gcc-9 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk"

# default values (to be overwritten if respective parameters are set)
OCTAVEDIR=/usr/local

# read the options
TEMP=`getopt -o o:m:f:v: --long octave:,matlab:,matlab-version:,fftw: -n 'macos-build-mex.sh' -- "$@"`
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
        -v|--matlab-version)
            case "$2" in
                "")  shift 2 ;;
                *) MATLABVERSION=$2; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
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
$(gsed -e '/^\[!/d' -e '/Directory structure/Q' $NFFTDIR/README) 
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
  $GCC -o .libs/lib"$LIB"julia.so -bundle .libs/lib"$LIB"julia.o -Wl,-force_load,../../.libs/libnfft3_julia.a $FFTW_LINK_COMMAND -lm -O3 -malign-double -march=$GCCARCH $OMPLIBS
  cd ..
done

cd fastsum
$GCC -o .libs/libfastsumjulia.so -bundle .libs/libfastsumjulia.o -Wl,-force_load,../../.libs/libnfft3_julia.a $FFTW_LINK_COMMAND -Wl,-force_load,../../applications/fastsum/.libs/libfastsum$THREADSSUFFIX.a -Wl,-force_load,../../applications/fastsum/.libs/libkernels.a -lm -O3 -malign-double -march=$GCCARCH $OMPLIBS
cd ..

cd "$NFFTBUILDDIR"

ARCH=$(uname -m)
JULIADIR=nfft-"$NFFTVERSION"-julia-macos_$ARCH$OMPSUFFIX
mkdir "$JULIADIR"
rsync -rLt --exclude='Makefile*' --exclude='doxygen*' --exclude='*.c.in' --exclude='*.c' --exclude='*.h' "$NFFTDIR/julia/" "$JULIADIR"
rsync -rLt --exclude='Makefile*' --exclude='.deps' --exclude='.libs' --exclude='*.la' --exclude='*.lo' --exclude='*.o' --exclude='*.c' 'julia/' "$JULIADIR"
echo 'This archive contains the Julia interface of NFFT '$NFFTVERSION'
compiled for '$ARCH' macOS using GCC '$GCCVERSION' and FFTW '$FFTWVERSION'.
' "$READMECONTENT" "$FFTWREADME" > "$JULIADIR"/readme.txt
zip -9 -r ../"$JULIADIR".zip "$JULIADIR"
# End of Julia interface


# Create Matlab/Octave release
for LIB in nfft nfsft nfsoft nnfft fastsum nfct nfst fpt
do
  cd matlab/"$LIB"
  $GCC -o .libs/lib"$LIB".mex -bundle  .libs/lib"$LIB"_la-"$LIB"mex.o -Wl,-force_load,../../.libs/libnfft3_matlab.a -Wl,-force_load,../../matlab/.libs/libmatlab.a -L"$OCTAVEDIR"/lib/octave/"$OCTAVEVERSION" $FFTW_LINK_COMMAND -lm -loctinterp -loctave -O3 -malign-double -march=$GCCARCH $OMPLIBS
  cd ../..
done

DIR=nfft-$NFFTVERSION-mexmaci64$OMPSUFFIX
mkdir $DIR
rsync -rLt --exclude='Makefile*' --exclude='doxygen*' --exclude='*.c.in' --exclude='*.c' --exclude='*.h' "$NFFTDIR/matlab/" "$DIR"
rsync -rLt --exclude='Makefile*' --exclude='.deps' --exclude='.libs' --exclude='*.la' --exclude='*.lo' --exclude='*.o' --exclude='*.c' "matlab/" "$DIR"

# Compile with Matlab
if [ -n "$MATLABDIR" ]; then
  if [ -z "$MATLABVERSION" ]; then
    MATLABVERSION=`"$MATLABDIR"/bin/matlab -wait -nodesktop -nosplash -r "fprintf('MATLAB_VERSION=%s\n', version); exit;" | grep MATLAB_VERSION | gsed 's/.*(//' | gsed 's/)//'`
  fi
  cd "$NFFTBUILDDIR"
  make clean
  CC=$GCC CPPFLAGS=-I"$FFTWDIR"/include LDFLAGS=-L"$FFTWDIR"/lib "$NFFTDIR/configure" --enable-all $OMPFLAG --with-matlab="$MATLABDIR" --with-gcc-arch=$GCCARCH --disable-static --enable-shared --disable-examples --disable-applications --enable-exhaustive-unit-tests
  make
  make check
  for LIB in nfft nfsft nfsoft nnfft fastsum nfct nfst fpt
  do
    cd matlab/"$LIB"
    $GCC -o .libs/lib"$LIB".mexmaci64 -bundle  .libs/lib"$LIB"_la-"$LIB"mex.o   -Wl,-force_load,../../.libs/libnfft3_matlab.a -Wl,-force_load,../../matlab/.libs/libmatlab.a  -L"$MATLABDIR"/bin/maci64 -lm -lmwfftw3 -lmx -lmex -lmat -O3 -malign-double -march=$GCCARCH $OMPLIBS
    cd ../..
  done
fi

for SUBDIR in nfft nfsft nfsoft nnfft fastsum nfct nfst fpt
  do
  cp -f -L matlab/$SUBDIR/*.mex* "$DIR"/$SUBDIR/
done

cd "$NFFTBUILDDIR"
cp "$NFFTDIR"/COPYING "$DIR"/COPYING
if [ -n "$MATLABDIR" ]; then
echo 'This archive contains the Matlab and Octave interface of NFFT '$NFFTVERSION'
compiled for '$ARCH' macOS using GCC '$GCCVERSION' and Matlab '$MATLABVERSION'
and Octave '$OCTAVEVERSION' and FFTW '$FFTWVERSION'.
' "$READMECONTENT" "$FFTWREADME" > "$DIR"/readme-matlab.txt
else
echo 'This archive contains the Octave interface of NFFT '$NFFTVERSION' compiled
for '$ARCH' macOS using GCC '$GCCVERSION' and Octave '$OCTAVEVERSION' and FFTW '$FFTWVERSION'.
' "$READMECONTENT" "$FFTWREADME" > "$DIR"/readme-matlab.txt
fi

zip -9 -r ../"$DIR".zip "$DIR"

done
