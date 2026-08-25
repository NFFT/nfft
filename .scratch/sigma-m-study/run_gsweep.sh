#!/bin/sh
# Build and run the geometry sweep in all three precisions.
# Usage: run_gsweep.sh <worktree-root> <out-dir> [m_max] [trials]
set -e
ROOT="$1"
OUT="$2"
MMAX="${3:-32}"
TRIALS="${4:-5}"
mkdir -p "$OUT"

for p in d f l; do
  case "$p" in
    d) prec=1; lib=nfft3;  tree=bx-d ;;
    f) prec=0; lib=nfft3f; tree=bx-f ;;
    l) prec=2; lib=nfft3l; tree=bx-l ;;
  esac
  gcc -O2 -std=c99 -I"$ROOT/include" -I"$ROOT/$tree" -DSWEEP_PREC=$prec \
      "$ROOT/.scratch/sigma-m-study/gsweep.c" -o "$OUT/gsweep-$p" \
      -L"$ROOT/$tree/kernel" -l$lib -lm
  LD_LIBRARY_PATH="$ROOT/$tree/kernel" "$OUT/gsweep-$p" "$MMAX" "$TRIALS" \
      > "$OUT/gsweep-$p.csv" 2> "$OUT/gsweep-$p.err"
  echo "done $p: $(wc -l < "$OUT/gsweep-$p.csv") rows"
done
