#!/usr/bin/env bash
# Build all three precision trees (float · double · long double) in one measurement mode.
#
# Deterministic Phase-A / Phase-E build step. Clears a stale Autotools config.h first
# (else the float/long-double trees mis-link), then configures + builds build-cmake{,-f,-l}
# with the CI-aligned flags. See ../details/phase-a-baseline.md and ../details/precision-matrix.md.
#
# Usage:  perf-build.sh [walltime|simulation] [extra cmake args...]
#   walltime  (default) — THE metric, local + CI (offline, no account); use this for the loop
#   simulation          — only for the optional callgrind tie-breaker (layout vs real work on an
#                         untouched control case); not a routine metric. See ../details/measurement-modes.md
# Optional env:  PERF_CODSPEED_SRC=/path/to/codspeed-cpp  reuses one checkout across trees.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$HERE/lib.sh"

MODE="${1:-walltime}"
case "$MODE" in
  walltime|simulation) ;;
  *) echo "usage: $(basename "$0") [walltime|simulation] [extra cmake args...]" >&2; exit 2 ;;
esac
shift || true

ROOT="$(perf_root)"; cd "$ROOT"

if [ -f include/config.h ]; then
  echo ">> clearing stale include/config.h (Autotools leftover) -> include/config.h.autotools.bak"
  mv include/config.h include/config.h.autotools.bak
fi

CFLAGS_STR="$(perf_cflags)"
# Array keeps the quoted CMAKE_C_FLAGS a SINGLE argument (unquoted it would split on spaces).
FLAGS=(-DNFFT_BENCHMARK_MODE="$MODE" -DNFFT_ENABLE_OPENMP=ON -DCMAKE_C_FLAGS="$CFLAGS_STR")

for p in "${PERF_PRECISIONS[@]}"; do
  t="$(perf_tree "$p")"
  extra=()
  case "$p" in
    f) extra+=(-DNFFT_ENABLE_FLOAT=ON) ;;
    l) extra+=(-DNFFT_ENABLE_LONG_DOUBLE=ON) ;;
  esac
  [ -n "${PERF_CODSPEED_SRC:-}" ] && extra+=(-DFETCHCONTENT_SOURCE_DIR_CODSPEED="$PERF_CODSPEED_SRC")
  echo ">> configuring $t  ($(perf_prec_name "$p"), mode=$MODE)"
  cmake -S . -B "$t" "${FLAGS[@]}" "${extra[@]}" "$@"
  echo ">> building $t"
  cmake --build "$t" -j
done

echo ">> built all precisions: trees build-cmake{,-f,-l}, mode=$MODE"
