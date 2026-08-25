#!/bin/sh
# Build and run the geometry sweep in all three precisions.
# Usage: run_gsweep.sh <worktree-root> <out-dir> [m_max] [trials] [Ns] [Mfacs] [sigmas]
# The three list arguments are comma-separated; an empty one keeps the default.
# <out-dir> also names the output files, so a second grid does not overwrite
# the first: pass TAG=<name> to change the gsweep-<p>.csv stem.
set -e
ROOT="$1"
OUT="$2"
MMAX="${3:-32}"
TRIALS="${4:-5}"
NS="${5:-}"
MFACS="${6:-}"
SIGMAS="${7:-}"
TAG="${TAG:-gsweep}"
mkdir -p "$OUT"

for p in d f l; do
  case "$p" in
    d) prec=1; lib=nfft3;  tree=bx-d ;;
    f) prec=0; lib=nfft3f; tree=bx-f ;;
    l) prec=2; lib=nfft3l; tree=bx-l ;;
  esac
  gcc -O2 -std=c99 -I"$ROOT/include" -I"$ROOT/$tree" -DSWEEP_PREC=$prec \
      "$ROOT/.scratch/sigma-m-study/gsweep.c" -o "$OUT/$TAG-$p" \
      -L"$ROOT/$tree/kernel" -l$lib -lm
  LD_LIBRARY_PATH="$ROOT/$tree/kernel" "$OUT/$TAG-$p" "$MMAX" "$TRIALS" "$NS" "$MFACS" "$SIGMAS" \
      > "$OUT/$TAG-$p.csv" 2> "$OUT/$TAG-$p.err"
  echo "done $p: $(wc -l < "$OUT/$TAG-$p.csv") rows"
done
