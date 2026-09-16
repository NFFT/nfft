#!/usr/bin/env bash
# Confirm (or dismiss) walltime regressions the noise rule flagged — the scripted form of
# "re-run before believing it" (../details/measurement-modes.md#the-noise-rule-this-is-the-metric).
#
# A single full capture can be inflated by transient background load, so perf-bench.py flags
# cases that may be pure noise. This re-measures — on the (presumably now quiet) machine — the
# benchmarks for every precision that had a flagged case, then re-applies the noise rule. A
# regression that shows up in BOTH the original final capture AND this independent re-run is
# real; one that evaporates was noise. (The deterministic clincher for an *untouched* control
# case is the instruction count — raw callgrind on a simulation build, the tie-breaker in
# ../details/measurement-modes.md — which is layout-insensitive; use it when the walltime story
# stays ambiguous.)
#
# Usage:  perf-confirm.sh <task-dir>
#   reads  <task-dir>/artifacts/{baseline,final}-bench-{d,f,l}.json
#   re-measures the affected precision(s) into a scratch dir (full bench, --bench-only style)
# Exit 0 = no flagged regression survived the re-run (all noise); 1 = at least one survived.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$HERE/lib.sh"

TASKDIR="${1:-}"
[ -n "$TASKDIR" ] || { echo "usage: $(basename "$0") <task-dir>" >&2; exit 2; }
ROOT="$(perf_root)"; cd "$ROOT"
ART="$TASKDIR/artifacts"
[ -d "$ART" ] || { echo "error: $ART not found (run a final capture first)" >&2; exit 2; }
command -v jq >/dev/null || { echo "error: jq not found" >&2; exit 2; }

BENCH="uv run python $HERE/perf-bench.py"
SCRATCH="$(mktemp -d)"; trap 'rm -rf "$SCRATCH"' EXIT

# Pass 1 — what the original final capture flagged.
$BENCH compare --taskdir "$TASKDIR" --emit-regressed > "$SCRATCH/pass1" || true
if [ ! -s "$SCRATCH/pass1" ]; then
  echo ">> no case flagged by the noise rule — nothing to confirm."
  exit 0
fi
n1=$(wc -l < "$SCRATCH/pass1")
echo ">> $n1 case(s) flagged in the final capture; re-measuring to see what survives:"
sed 's/^/   /' "$SCRATCH/pass1" >&2

# Re-measure each precision that had a flag (full bench, independent run), into scratch.
: > "$SCRATCH/pass2"
for p in $(cut -f1 "$SCRATCH/pass1" | sort -u); do
  t="$(perf_tree "$p")"; bin="$t/benchmarks/bench_nfft_direct"
  [ -x "$bin" ] || { echo "WARN: [$p] no benchmark binary ($bin) — cannot confirm" >&2; continue; }
  la="$(cut -d' ' -f1 /proc/loadavg 2>/dev/null || echo '?')"
  echo ">> [$p] re-measuring (loadavg(1m)=$la) ..." >&2
  pf="$SCRATCH/$p"; mkdir -p "$pf"
  CODSPEED_PROFILE_FOLDER="$pf" "$bin" >/dev/null 2>&1
  jq -s '[.[].benchmarks[]]' "$pf"/results/*.json > "$SCRATCH/confirm-$p.json"
  $BENCH compare --base "$ART/baseline-bench-$p.json" --final "$SCRATCH/confirm-$p.json" \
    --prec "$p" --emit-regressed >> "$SCRATCH/pass2" || true
done

# Survivors = flagged in BOTH passes.
sort -u "$SCRATCH/pass1" > "$SCRATCH/p1s"; sort -u "$SCRATCH/pass2" > "$SCRATCH/p2s"
comm -12 "$SCRATCH/p1s" "$SCRATCH/p2s" > "$SCRATCH/survivors"
comm -23 "$SCRATCH/p1s" "$SCRATCH/p2s" > "$SCRATCH/evaporated"

echo
if [ -s "$SCRATCH/evaporated" ]; then
  echo ">> evaporated on re-run (noise, not a real regression):"
  sed 's/^/   /' "$SCRATCH/evaporated"
fi
if [ -s "$SCRATCH/survivors" ]; then
  echo ">> CONFIRMED — regressed in both runs (real; investigate or attribute to code layout):"
  sed 's/^/   /' "$SCRATCH/survivors"
  echo ">> For an UNTOUCHED control case, confirm layout-vs-work with the deterministic"
  echo ">> instruction count (raw callgrind on a simulation build) — identical I refs ⇒ layout."
  exit 1
fi
echo ">> all flagged case(s) were noise — none survived the re-run."
exit 0
