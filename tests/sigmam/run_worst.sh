#!/bin/sh
# Regenerate the worst-case sweeps kbfit.py fits kbmodel.json to.
#
# usage: run_worst.sh <src-root> <tree-d> <tree-f> <tree-l> <work-dir>
#
# <tree-x> is a build tree of that precision whose kernel/ (CMake) or .libs/
# (Autotools) holds libnfft3{,f,l}. One transform per k, f_hat = e_k, which is
# where the measure attains its supremum; every k is probed up to FULL_K_MAX.
# Nodes carry each type's own mantissa. Long double is quad on aarch64 and
# 80-bit on x86: rerun there and store the result as a separate type.
# Then: uv run --with scipy python kbfit.py [--write]
set -e
SRC="$1"; OUT="$5"
HERE="$(cd "$(dirname "$0")" && pwd)"
SIG=1.25,1.3,1.375,1.5,1.625,1.75,1.875,2.0,2.125,2.25,2.375,2.5
mkdir -p "$OUT"

gen() { # file bits NS FULL_K_MAX
  NS="$3" FULL_K_MAX="$4" uv run -q --with mpmath==1.3.0 \
    python "$HERE/floorgen.py" "$OUT/$1" "$2" 200
}
gen ref-d.txt 53 32,128,256,512,1024,4096 4096
gen ref-f.txt 24 32,128,256,512,1024,4096 4096
gen ref-l.txt 113 32,128,256,1024 256

for spec in "d 1 $2 24" "f 0 $3 14" "l 2 $4 28"; do
  set -- $spec
  suf=$1; [ $1 = d ] && suf=
  lib="$3/kernel"; [ -e "$lib/libnfft3$suf.so" ] || lib="$3/.libs"
  gcc -O2 -std=c99 -I"$SRC/include" -DSWEEP_PREC=$2 "$HERE/floor.c" \
      -o "$OUT/floor-$1" -L"$lib" -lnfft3$suf -lfftw3$suf -lm
  LD_LIBRARY_PATH="$lib" "$OUT/floor-$1" "$OUT/ref-$1.txt" $4 $SIG \
      | gzip > "$HERE/data/floor-worst-width-$1.csv.gz"
  echo "done $1"
done

uv run -q --with mpmath==1.3.0 --with numpy python "$HERE/exact_scan.py" \
  | gzip > "$HERE/data/trunc-exact-width.csv.gz"
echo "done exact scan"
