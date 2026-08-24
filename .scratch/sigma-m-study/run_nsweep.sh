#!/bin/sh
# Build and run the fine-grained n sweep (double precision).
# Usage: run_nsweep.sh <worktree-root> <out-dir> [m] [reps]
set -e
ROOT="$1"
OUT="$2"
M="${3:-6}"
REPS="${4:-100}"
mkdir -p "$OUT"
gcc -O2 -std=gnu99 -I"$ROOT/include" -I"$ROOT/bx-d" \
    "$ROOT/.scratch/sigma-m-study/nsweep.c" -o "$OUT/nsweep" \
    -L"$ROOT/bx-d/kernel" -lnfft3 -lfftw3 -lm
LD_LIBRARY_PATH="$ROOT/bx-d/kernel" "$OUT/nsweep" "$M" "$REPS" \
    > "$OUT/nsweep.csv" 2> "$OUT/nsweep.err"
echo "nsweep done: $(wc -l < "$OUT/nsweep.csv") rows"
