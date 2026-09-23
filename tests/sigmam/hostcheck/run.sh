#!/bin/sh
# Build NFFT3 in double, then time nfft_sigma_m_1d's grid against the legacy
# grid on this host. See README.md.
#
# usage: tests/sigmam/hostcheck/run.sh [quick|full] [out-dir]
#
# Environment: CC, CFLAGS, LDFLAGS, and FFTW_DIR (prefix of an FFTW3 install
# that the compiler and CMake do not find by themselves).
set -e

MODE="${1:-quick}"
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$(cd "$HERE/../../.." && pwd)"
OUT="${2:-$SRC/build/hostcheck}"
BUILD="$OUT/nfft"
CASES="$HERE/cases-$MODE.txt"

[ -f "$CASES" ] || { echo "unknown mode '$MODE': quick or full" >&2; exit 1; }
for tool in cmake python3 "${CC:-cc}"; do
  command -v "$tool" >/dev/null 2>&1 || { echo "missing: $tool" >&2; exit 1; }
done
mkdir -p "$OUT"

FFTW_CPP=""
FFTW_LD=""
if [ -n "$FFTW_DIR" ]; then
  FFTW_CPP="-I$FFTW_DIR/include"
  FFTW_LD="-L$FFTW_DIR/lib -Wl,-rpath,$FFTW_DIR/lib"
  export CMAKE_PREFIX_PATH="$FFTW_DIR${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
elif command -v pkg-config >/dev/null 2>&1 && pkg-config --exists fftw3; then
  FFTW_CPP="$(pkg-config --cflags fftw3)"
  FFTW_LD="$(pkg-config --libs-only-L fftw3)"
fi

echo "building NFFT3 (double, serial) in $BUILD"
cmake -S "$SRC" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
  -DNFFT_ENABLE_OPENMP=OFF -DNFFT_ENABLE_TESTS=OFF \
  -DNFFT_ENABLE_EXAMPLES=OFF -DNFFT_ENABLE_APPLICATIONS=OFF \
  > "$OUT/cmake.log" 2>&1 || { echo "cmake failed, see $OUT/cmake.log" >&2; exit 1; }
cmake --build "$BUILD" -j > "$OUT/build.log" 2>&1 \
  || { echo "build failed, see $OUT/build.log" >&2; exit 1; }

LIB="$(dirname "$(find "$BUILD" -name 'libnfft3.*' \( -name '*.so' -o -name '*.dylib' \) | head -n 1)")"
${CC:-cc} -O2 -std=gnu99 $CFLAGS -I"$SRC/include" $FFTW_CPP \
  "$HERE/hostcheck.c" -o "$OUT/hostcheck" \
  -L"$LIB" -Wl,-rpath,"$LIB" -lnfft3 $FFTW_LD $LDFLAGS -lfftw3 -lm

RESULTS="$OUT/results-$MODE.csv"
{
  echo "# host $(uname -srm)"
  if command -v lscpu >/dev/null 2>&1; then
    echo "# cpu $(lscpu | sed -n 's/^Model name: *//p' | head -n 1)"
  elif command -v sysctl >/dev/null 2>&1; then
    echo "# cpu $(sysctl -n machdep.cpu.brand_string 2>/dev/null)"
  fi
  echo "# cc $(${CC:-cc} --version 2>/dev/null | head -n 1)"
  echo "# nfft $(git -C "$SRC" rev-parse --short HEAD 2>/dev/null || echo unknown)"
  echo "# mode $MODE, $(wc -l < "$CASES" | tr -d ' ') cases"
} > "$RESULTS"

echo "running $(wc -l < "$CASES" | tr -d ' ') cases ($MODE); results in $RESULTS"
"$OUT/hostcheck" < "$CASES" >> "$RESULTS"

python3 "$HERE/summarize.py" "$RESULTS" | tee "$OUT/summary-$MODE.txt"
echo
echo "summary: $OUT/summary-$MODE.txt"
