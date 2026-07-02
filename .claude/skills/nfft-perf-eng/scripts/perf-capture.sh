#!/usr/bin/env bash
# Capture the FULL test + benchmark state for all precisions into a task dir's artifacts/.
#
# Deterministic Phase-A (baseline) and Phase-E (final) capture. For each precision tree it
# runs the full ctest suite (tee'd to a log) and the benchmark binary (all cases), then
# collates the per-process codspeed JSON into one flat array per precision. See
# ../details/phase-a-baseline.md, ../details/phase-f-exit-gate.md, ../details/precision-matrix.md.
#
# Usage:  perf-capture.sh <baseline|final|confirm> <task-dir> [--bench-only|--tests-only]
#                         [--prec d|f|l] [--filter REGEX]
#   e.g.  perf-capture.sh baseline .perfeng
#         perf-capture.sh final    .perfeng --bench-only   # re-measure cleanly
# Flags:
#   --bench-only   skip ctest (tests are deterministic — re-measure benchmarks without the slow
#                  long-double suite; the recommended way to re-run after a noisy capture)
#   --tests-only   skip the benchmark run
#   --prec P       restrict to one precision tree (d|f|l)
#   --filter RE    pass --benchmark_filter=RE to the benchmark binary (subset of cases)
# Writes:  <task-dir>/artifacts/<phase>-tests-{d,f,l}.log
#          <task-dir>/artifacts/<phase>-bench-{d,f,l}.json   (flat array of benchmark objects)
# Each precision logs an ISO timestamp + 1-min loadavg (quiescence signal) and WARNs if the
# machine is loaded — walltime medians are inflated by concurrent work, so capture on a quiet box.
# Exit code 0 = every (selected) precision green and captured; 1 = a precision was not fully
# green or could not be captured (details on stderr) — the caller decides what that means
# (Phase A: stop; Phase E: a real regression vs a missing precision).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$HERE/lib.sh"

PHASE="${1:-}"; TASKDIR="${2:-}"
case "$PHASE" in
  baseline|final|confirm) ;;
  *) echo "usage: $(basename "$0") <baseline|final|confirm> <task-dir> [--bench-only|--tests-only] [--prec d|f|l] [--filter RE]" >&2
     echo "  e.g. $(basename "$0") baseline .perfeng" >&2; exit 2 ;;
esac
[ -n "$TASKDIR" ] || { echo "error: task dir required" >&2; exit 2; }
shift 2 2>/dev/null || true

DO_TESTS=1; DO_BENCH=1; ONLY_PREC=""; FILTER=""
while [ $# -gt 0 ]; do
  case "$1" in
    --bench-only) DO_TESTS=0 ;;
    --tests-only) DO_BENCH=0 ;;
    --prec) shift; ONLY_PREC="${1:-}" ;;
    --filter) shift; FILTER="${1:-}" ;;
    *) echo "error: unknown option '$1'" >&2; exit 2 ;;
  esac
  shift
done
{ [ "$DO_TESTS" -eq 1 ] || [ "$DO_BENCH" -eq 1 ]; } || {
  echo "error: --bench-only and --tests-only are mutually exclusive" >&2; exit 2; }

ROOT="$(perf_root)"; cd "$ROOT"
[ -d "$TASKDIR" ] || { echo "error: task dir '$TASKDIR' not found (run perf-init.sh first)" >&2; exit 2; }
command -v jq >/dev/null || { echo "error: jq not found (needed to collate benchmark JSON)" >&2; exit 2; }

ART="$TASKDIR/artifacts"; mkdir -p "$ART"
SCRATCH="$(mktemp -d)"
trap 'rm -rf "$SCRATCH"' EXIT
green=1

for p in "${PERF_PRECISIONS[@]}"; do
  [ -n "$ONLY_PREC" ] && [ "$p" != "$ONLY_PREC" ] && continue
  t="$(perf_tree "$p")"
  if [ ! -d "$t" ]; then
    echo "WARN: [$p] tree $t missing — skipping (run perf-build.sh first)" >&2; green=0; continue
  fi

  # Quiescence signal — walltime medians are inflated by concurrent load. Record + warn.
  la="$(cut -d' ' -f1 /proc/loadavg 2>/dev/null || echo '?')"
  ncpu="$(nproc 2>/dev/null || echo '?')"
  echo ">> [$p] $(date -u +%FT%TZ)  loadavg(1m)=$la / ${ncpu} cpu" >&2
  case "$la$ncpu" in
    *'?'*) ;;
    *) if awk "BEGIN{exit !($la > 0.6*$ncpu)}" 2>/dev/null; then
         echo "WARN: [$p] load $la is high vs ${ncpu} cpu — quiesce the machine (no concurrent" \
              "builds/sweeps) before a capture; walltime medians may be inflated" >&2
       fi ;;
  esac

  if [ "$DO_TESTS" -eq 1 ]; then
    echo ">> [$p] ctest ($t)"
    if ! ctest --test-dir "$t" 2>&1 | tee "$ART/$PHASE-tests-$p.log"; then
      echo "WARN: [$p] ctest reported failures — see $ART/$PHASE-tests-$p.log" >&2; green=0
    fi
  fi

  if [ "$DO_BENCH" -eq 1 ]; then
    bin="$t/benchmarks/bench_nfft_direct"
    if [ ! -x "$bin" ]; then
      echo "WARN: [$p] benchmark binary not built ($bin)" >&2; green=0; continue
    fi
    echo ">> [$p] benchmark ($bin, ${FILTER:-all cases})"
    pf="$SCRATCH/$p"; mkdir -p "$pf"
    # provenance header into the (already-linked) bench log, so a noisy capture is identifiable
    echo "# $PHASE capture $(date -u +%FT%TZ) loadavg(1m)=$la/${ncpu}cpu filter=${FILTER:-<all>}" \
      > "$ART/$PHASE-bench-$p.log"
    fargs=(); [ -n "$FILTER" ] && fargs+=(--benchmark_filter="$FILTER")
    if ! CODSPEED_PROFILE_FOLDER="$pf" "$bin" "${fargs[@]}" 2>>"$ART/$PHASE-bench-$p.log"; then
      echo "WARN: [$p] benchmark binary failed — see $ART/$PHASE-bench-$p.log" >&2; green=0; continue
    fi
    if compgen -G "$pf/results/*.json" >/dev/null; then
      jq -s '[.[].benchmarks[]]' "$pf"/results/*.json > "$ART/$PHASE-bench-$p.json"
    else
      echo "WARN: [$p] no benchmark results JSON produced" >&2; green=0
    fi
  fi
done

echo ">> capture complete -> $ART/$PHASE-{tests,bench}-{d,f,l}.{log,json}"
if [ "$green" -eq 1 ]; then
  echo ">> all precisions green and captured"
else
  echo ">> NOTE: some precision not fully green/captured — see WARNs above" >&2
fi
exit $(( green == 1 ? 0 : 1 ))
