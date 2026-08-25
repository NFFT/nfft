#!/bin/sh
# Build and run the tuned-vs-legacy comparison in all three precisions.
# Usage: run_compare.sh <worktree-root> <out-dir> [reps] [refine]
set -e
ROOT="$1"
OUT="$2"
REPS="${3:-50}"
REFINE="${4:-0}"
mkdir -p "$OUT"

for p in d f l; do
  case "$p" in
    d) prec=1; lib=nfft3;  tree=bx-d ;;
    f) prec=0; lib=nfft3f; tree=bx-f ;;
    l) prec=2; lib=nfft3l; tree=bx-l ;;
  esac
  gcc -O2 -std=gnu99 -I"$ROOT/include" -I"$ROOT/$tree" -DSWEEP_PREC=$prec \
      "$ROOT/.scratch/sigma-m-study/compare.c" -o "$OUT/compare-$p" \
      -L"$ROOT/$tree/kernel" -l$lib -lfftw3 -lm
  LD_LIBRARY_PATH="$ROOT/$tree/kernel" "$OUT/compare-$p" "$REPS" "$REFINE" \
      > "$OUT/compare-$p.csv" 2> "$OUT/compare-$p.err"
  echo "done $p: $(wc -l < "$OUT/compare-$p.csv") rows"
done
