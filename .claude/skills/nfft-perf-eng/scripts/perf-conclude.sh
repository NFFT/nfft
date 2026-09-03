#!/usr/bin/env bash
# Phase G (conclude): the DETERMINISTIC wrap-up, in two steps that match the publish flow —
# because the PR number (which the deliverables archive is named for) only exists once the PR
# is created, i.e. after the squash. The interactive steps between them (push, open the PR,
# attach the zip) stay with the agent. See ../details/phase-g-conclude.md.
#
#   perf-conclude.sh squash -m "<message>"                 # after Phase F, BEFORE push/PR
#       Preflight — refuse unless: base recorded in .perfeng/BASE, working tree clean
#       (deliverables are gitignored, so a clean tree means all SOURCE work is committed), NOT on
#       the default branch, and the base is an ancestor of HEAD. Then collapse every commit
#       base..HEAD into ONE via `git reset --soft <base>` + commit (non-interactive; interactive
#       rebase is unavailable in this environment). The gitignored .perfeng/ never enters the commit.
#
#   perf-conclude.sh package <pr-number> [--zip-out <path>]  # AFTER the PR exists, number known
#       Rename .perfeng -> .perfeng-pr-<N> (no leading zeros) so the archive unpacks into a
#       PR-unique directory — unambiguous within this repo when several results zips are
#       downloaded — then zip it to perfeng-pr-<N>.zip OUTSIDE the tree (standard ZIP, for
#       maximum compatibility), ready to attach to the PR. Re-run it after any follow-up work
#       (post-conclude code/deliverable changes) to refresh the archive; when the dir is already
#       renamed it re-packages in place and tells you to OVERWRITE the zip in the existing PR
#       comment rather than attach a second one.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$HERE/lib.sh"
ROOT="$(perf_root)"; cd "$ROOT"
TASKDIR=".perfeng"

cmd_squash() {
  local MSG=""
  while [ $# -gt 0 ]; do
    case "$1" in
      -m|--message) MSG="${2:-}"; shift 2 ;;
      *) echo "usage: $(basename "$0") squash -m \"<commit message>\"" >&2; exit 2 ;;
    esac
  done
  [ -n "$MSG" ] || { echo "error: -m \"<commit message>\" is required" >&2; exit 2; }
  [ -d "$TASKDIR" ] || { echo "error: $TASKDIR missing — nothing to conclude (run perf-init first)" >&2; exit 2; }
  [ -f "$TASKDIR/BASE" ] || { echo "error: $TASKDIR/BASE missing — cannot determine the squash base" >&2; exit 2; }

  local BASE_SHA BASE_BRANCH BRANCH DEFAULT N_COMMITS
  BASE_SHA="$(sed -n 's/^sha=//p' "$TASKDIR/BASE")"
  BASE_BRANCH="$(sed -n 's/^branch=//p' "$TASKDIR/BASE")"
  [ -n "$BASE_SHA" ] || { echo "error: no sha= in $TASKDIR/BASE" >&2; exit 2; }
  BRANCH="$(git rev-parse --abbrev-ref HEAD)"
  DEFAULT="$(git rev-parse --abbrev-ref origin/HEAD 2>/dev/null | sed 's#^origin/##' || echo '')"

  # Preflight.
  if [ "$BRANCH" = "$DEFAULT" ] || [ "$BRANCH" = "HEAD" ]; then
    echo "error: on '$BRANCH' — squash/PR must run on a feature branch, not the default branch." >&2
    echo "       create a branch first (git switch -c perf/<slug>) and re-run." >&2
    exit 2
  fi
  if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "error: working tree not clean — commit (or revert) all SOURCE changes before concluding." >&2
    echo "       (deliverables under $TASKDIR/ are gitignored and never committed.)" >&2
    git status --short >&2; exit 2
  fi
  if ! git merge-base --is-ancestor "$BASE_SHA" HEAD; then
    echo "error: recorded base $BASE_SHA is not an ancestor of HEAD — refusing to squash." >&2; exit 2
  fi

  N_COMMITS="$(git rev-list --count "$BASE_SHA"..HEAD)"
  echo ">> branch '$BRANCH', base ${BASE_SHA:0:12} (recorded on '$BASE_BRANCH'), $N_COMMITS commit(s) since."
  if [ "$N_COMMITS" -eq 0 ]; then
    echo ">> no commits since base — nothing to squash. Ensure the optimization is committed." >&2; exit 2
  elif [ "$N_COMMITS" -eq 1 ]; then
    echo ">> already a single commit — amending its message only."
    git commit --amend -m "$MSG" >/dev/null
  else
    git reset --soft "$BASE_SHA"
    git commit -m "$MSG" >/dev/null
  fi
  echo ">> squashed to one commit: $(git rev-parse --short HEAD)  \"$MSG\""

  cat <<EOF

>> NEXT (interactive — confirm with the user first):
   push:    git push -u origin $BRANCH            # add --force-with-lease if the branch was pushed before the squash
   label:   gh label list | grep -q '^perf-eng' || gh label create perf-eng --description "nfft-perf-eng optimization run" --color 0e8a16
   PR:      gh pr create --base ${DEFAULT:-develop} --label perf-eng --title "<title>" --body "<summary>"
   package: $(basename "$0") package <pr-number>  # then rename .perfeng -> .perfeng-pr-<N> and zip it to attach
EOF
}

cmd_package() {
  local PR="${1:-}" ZIP_OUT=""
  shift 2>/dev/null || true
  while [ $# -gt 0 ]; do
    case "$1" in
      --zip-out) ZIP_OUT="${2:-}"; shift 2 ;;
      *) echo "usage: $(basename "$0") package <pr-number> [--zip-out <path>]" >&2; exit 2 ;;
    esac
  done
  case "$PR" in
    ''|*[!0-9]*) echo "error: package needs a numeric PR number (no leading zeros), e.g. '$(basename "$0") package 231'" >&2; exit 2 ;;
  esac
  PR=$((10#$PR))                       # normalise — strip any leading zeros
  local NAMED=".perfeng-pr-$PR" REPACKAGE=0

  if [ -d "$TASKDIR" ]; then
    [ -e "$NAMED" ] && { echo "error: $NAMED already exists — refusing to overwrite" >&2; exit 2; }
    mv "$TASKDIR" "$NAMED"
    echo ">> renamed $TASKDIR -> $NAMED"
  elif [ -d "$NAMED" ]; then
    REPACKAGE=1                         # already renamed once — this is a follow-up refresh
    echo ">> $NAMED already present — re-packaging (follow-up update)"
  else
    echo "error: neither $TASKDIR nor $NAMED exists — nothing to package" >&2; exit 2
  fi

  [ -n "$ZIP_OUT" ] || ZIP_OUT="$(dirname "$ROOT")/perfeng-pr-$PR.zip"
  rm -f "$ZIP_OUT"
  # Standard ZIP (default DEFLATE) so any unzip / Explorer / Finder opens it — no exotic codecs.
  ( cd "$ROOT" && zip -q -r "$ZIP_OUT" "$NAMED" )
  echo ">> deliverables archived: $ZIP_OUT   (standard zip; unzips into ./$NAMED/)"

  if [ "$REPACKAGE" -eq 1 ]; then
    cat <<EOF

>> UPDATE the deliverables on PR #$PR — the archive changed after follow-up work:
   edit the existing PR comment holding the old perfeng-pr-$PR.zip, remove that attachment, and
   drag-drop the new $ZIP_OUT into the SAME comment (overwrite it — same filename, one archive).
   (gh has no binary PR-attach, so this stays a manual web-UI step.)
EOF
  else
    cat <<EOF

>> attach $ZIP_OUT to PR #$PR — drag-drop it into a PR comment in the web UI
   (gh has no binary PR-attach; the deliverables are gitignored so they never ride in the commit).
EOF
  fi
}

SUB="${1:-}"; shift 2>/dev/null || true
case "$SUB" in
  squash)  cmd_squash "$@" ;;
  package) cmd_package "$@" ;;
  *) echo "usage: $(basename "$0") <squash -m \"MSG\" | package <pr-number>>" >&2; exit 2 ;;
esac
