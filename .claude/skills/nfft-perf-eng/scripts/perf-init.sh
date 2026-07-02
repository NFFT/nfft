#!/usr/bin/env bash
# Step 0: bootstrap the single static, gitignored deliverables directory `.perfeng/` for
# one optimization run.
#
# Creates `.perfeng/` (with artifacts/), copies every deliverable template in (tracker ->
# README.md, the phase docs, and the error-analysis HTML), stamps the baseline commit, and
# records the squash base — the commit + branch HEAD points at right now — into `.perfeng/BASE`
# so the conclude step can collapse the run's intermediate commits deterministically.
#
# `.perfeng/` is a FIXED path, never committed (it is gitignored), and holds only the CURRENT
# run — its record ships as the run's squashed commit + PR + the zip attached to that PR, not
# in the repo tree. There is no NNNN sequence and no shared in-repo index. If `.perfeng/`
# already exists it is moved aside to `.perfeng.bak/` (one-level undo) and recreated fresh.
# Deterministic scaffold. See ../details/deliverables.md.
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
TASKDIR=".perfeng"

# Single static dir: preserve any prior run as a one-level undo, then start clean.
if [ -e "$TASKDIR" ]; then
  echo ">> $TASKDIR exists — moving it aside to $TASKDIR.bak (overwriting any previous .bak)" >&2
  rm -rf "$TASKDIR.bak"
  mv "$TASKDIR" "$TASKDIR.bak"
fi

mkdir -p "$TASKDIR/artifacts"
cp "$TPL/tracker.md" "$TASKDIR/README.md"
for f in phase-a-baseline phase-b-correctness-net phase-c-performance-metric \
         phase-d-error-analysis phase-e-inner-loop phase-f-exit-gate; do
  cp "$TPL/$f.md" "$TASKDIR/$f.md"
done
cp "$TPL/error-analysis.html" "$TASKDIR/error-analysis.html"

SHA="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
FULLSHA="$(git rev-parse HEAD 2>/dev/null || echo unknown)"
BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
# Stamp the baseline commit where the templates leave a <sha> slot.
sed -i "s/<sha>/$SHA/g" "$TASKDIR/README.md" "$TASKDIR/phase-a-baseline.md" 2>/dev/null || true
# Record the squash base for the conclude step (Phase G): all commits base..HEAD collapse to one.
{ echo "sha=$FULLSHA"; echo "branch=$BRANCH"; echo "slug=$SLUG"; } > "$TASKDIR/BASE"

echo ">> created $TASKDIR  (baseline / squash-base commit $SHA on $BRANCH)"
echo ">> scaffolded: README.md (tracker) + phase {a..f} deliverables + error-analysis.html + artifacts/ + BASE"
echo ">> $TASKDIR is gitignored — deliverables are never committed; commit only source changes as you go."
echo ">> TODO (manual): set the Target line in $TASKDIR/README.md"
