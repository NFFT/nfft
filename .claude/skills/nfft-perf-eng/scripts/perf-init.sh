#!/usr/bin/env bash
# Step 0: bootstrap a task directory under docs/perfeng/ for one optimization.
#
# Picks the next NNNN sequence number, creates docs/perfeng/NNNN-<slug>/ with artifacts/,
# copies every deliverable template in (tracker -> README.md, the phase docs, and the
# error-analysis HTML), and stamps the baseline commit. It prints — but does NOT mutate —
# the index row to add to docs/perfeng/README.md (kept manual so the shared index is never
# silently rewritten). Deterministic scaffold. See ../details/deliverables.md.
#
# Usage:  perf-init.sh <target-slug>     e.g.  perf-init.sh trafo-direct
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$HERE/lib.sh"

SLUG="${1:-}"
[ -n "$SLUG" ] || { echo "usage: $(basename "$0") <target-slug>   e.g. trafo-direct" >&2; exit 2; }
case "$SLUG" in
  *[!a-z0-9-]*) echo "error: slug must be lower-case [a-z0-9-] only (got '$SLUG')" >&2; exit 2 ;;
esac

ROOT="$(perf_root)"; cd "$ROOT"
TPL="$HERE/../templates"
PERFDIR="docs/perfeng"
[ -d "$PERFDIR" ] || { echo "error: $PERFDIR missing" >&2; exit 2; }

# Next NNNN = highest existing zero-padded prefix + 1. Consider BOTH the task dirs AND any
# NNNN already referenced in the index README — a row may exist before its dir (e.g. a shipped
# worked-example row), and perf-init must never reuse such a number (that was the original
# collision bug). See ../details/deliverables.md.
n=0
for d in "$PERFDIR"/[0-9][0-9][0-9][0-9]-*/; do
  [ -d "$d" ] || continue
  b="$(basename "$d")"; num="${b%%-*}"; num=$((10#$num))
  [ "$num" -gt "$n" ] && n="$num"
done
IDX="$PERFDIR/README.md"
if [ -f "$IDX" ]; then
  maxidx="$(grep -oE '[0-9]{4}-[a-z0-9-]+' "$IDX" | grep -oE '^[0-9]{4}' | sort -rn | head -1 || true)"
  [ -n "$maxidx" ] && [ "$((10#$maxidx))" -gt "$n" ] && n="$((10#$maxidx))"
  if grep -qE "[0-9]{4}-$SLUG[]/]" "$IDX"; then
    echo "warning: index $IDX already references slug '$SLUG' — check for a stale/duplicate row" >&2
  fi
fi
NNNN="$(printf '%04d' $((n + 1)))"
TASKDIR="$PERFDIR/$NNNN-$SLUG"
[ -e "$TASKDIR" ] && { echo "error: $TASKDIR already exists" >&2; exit 2; }

mkdir -p "$TASKDIR/artifacts"
cp "$TPL/tracker.md" "$TASKDIR/README.md"
for f in phase-a-baseline phase-b-correctness-net phase-c-performance-metric \
         phase-d-error-analysis phase-e-inner-loop phase-f-exit-gate; do
  cp "$TPL/$f.md" "$TASKDIR/$f.md"
done
cp "$TPL/error-analysis.html" "$TASKDIR/error-analysis.html"

SHA="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
# Stamp the baseline commit where the templates leave a <sha> slot.
sed -i "s/<sha>/$SHA/g" "$TASKDIR/README.md" "$TASKDIR/phase-a-baseline.md" 2>/dev/null || true

echo ">> created $TASKDIR  (baseline commit $SHA)"
echo ">> scaffolded: README.md (tracker) + phase {a..f} deliverables + error-analysis.html + artifacts/"
echo ">> TODO (manual):"
echo "   1. set the Target line in $TASKDIR/README.md"
echo "   2. add this row to $PERFDIR/README.md index table:"
echo "      | [$NNNN-$SLUG]($NNNN-$SLUG/) | <target — file:line> | in-progress | — |"
