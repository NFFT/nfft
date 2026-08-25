#!/bin/sh
# Build and run the cost-weight timing sweep in all three precisions.
# Usage: run_costfit.sh <worktree-root> <out-dir> [seconds-per-point]
set -e
ROOT="$1"
OUT="$2"
BUDGET="${3:-0.02}"
mkdir -p "$OUT"

for p in d f l; do
  case "$p" in
    d) prec=1; lib=nfft3;  tree=bx-d ;;
    f) prec=0; lib=nfft3f; tree=bx-f ;;
    l) prec=2; lib=nfft3l; tree=bx-l ;;
  esac
  gcc -O2 -std=gnu99 -I"$ROOT/include" -I"$ROOT/$tree" -DSWEEP_PREC=$prec \
      "$ROOT/.scratch/sigma-m-study/costfit.c" -o "$OUT/costfit-$p" \
      -L"$ROOT/$tree/kernel" -l$lib -lfftw3 -lm
  LD_LIBRARY_PATH="$ROOT/$tree/kernel" "$OUT/costfit-$p" "$BUDGET" \
      > "$OUT/costfit-$p.csv" 2> "$OUT/costfit-$p.err"
  echo "done $p: $(wc -l < "$OUT/costfit-$p.csv") rows"
done
