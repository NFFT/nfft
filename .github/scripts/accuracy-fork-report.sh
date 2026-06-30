#!/usr/bin/env bash
# Fork-PR accuracy report (emoji comment + Check only), run from the trusted
# default branch via workflow_run. Best-effort: no failure here may matter.
set -uo pipefail
owner="${GITHUB_REPOSITORY%%/*}"; repo="${GITHUB_REPOSITORY##*/}"
baseline_raw="https://raw.githubusercontent.com/${owner}/${repo}/gh-pages/baseline"

# Resolve the PR number from the build's head SHA (workflow_run.pull_requests is
# unreliable for forks).
pr="$(gh api "repos/${GITHUB_REPOSITORY}/commits/${HEAD_SHA}/pulls" \
  --jq 'first(.[].number) // empty')"
[ -n "$pr" ] || { echo "no open PR for ${HEAD_SHA}; nothing to do"; exit 0; }

# Download the fork build's BMF artifacts from the triggering run (data only).
gh run download "${RUN_ID}" --pattern "accuracy-bmf-*" --dir bmf-artifacts \
  || { echo "no accuracy artifacts on run ${RUN_ID}; nothing to do"; exit 0; }

shopt -s nullglob
mkdir -p pr-bmf
for d in bmf-artifacts/accuracy-bmf-*/; do
  tb="$(basename "$d")"; tb="${tb#accuracy-bmf-}"
  [ -f "${d}accuracy.bmf.json" ] && cp "${d}accuracy.bmf.json" "pr-bmf/${tb}.bmf.json"
done
ls pr-bmf/*.bmf.json >/dev/null 2>&1 || { echo "no BMFs; nothing to do"; exit 0; }

# Fetch baseline (skip cleanly if develop has not published one yet).
mkdir -p base-bmf
for f in pr-bmf/*.bmf.json; do
  tb="$(basename "$f")"
  curl -fsSL --remove-on-error "${baseline_raw}/${tb}" -o "base-bmf/${tb}" \
    || echo "no baseline yet: ${tb}"
done

# Render WITHOUT png URLs -> the comment carries the inline emoji grid + lists.
if ls base-bmf/*.bmf.json >/dev/null 2>&1; then
  uv run --with matplotlib python tests/bench/pr_report.py pr-bmf base-bmf out
else
  uv run --with matplotlib python tests/bench/pr_report.py pr-bmf pr-bmf out --no-baseline
fi

# Post Check (PR head SHA) + upsert comment — non-fatal.
title="$(jq -r .title out/check.json)"; summary="$(jq -r .summary out/check.json)"
gh api -X POST "repos/${GITHUB_REPOSITORY}/check-runs" \
  -f name="accuracy" -f head_sha="${HEAD_SHA}" -f status="completed" \
  -f conclusion="neutral" -f output[title]="$title" -f output[summary]="$summary" \
  || echo "check post failed (non-fatal)"

body="$(cat out/comment.md)"
existing="$(gh api "repos/${GITHUB_REPOSITORY}/issues/${pr}/comments" \
  --jq 'first(.[] | select(.body | startswith("<!-- nfft-accuracy-report -->")) | .id) // empty')"
if [ -n "$existing" ]; then
  gh api -X PATCH "repos/${GITHUB_REPOSITORY}/issues/comments/${existing}" -f body="$body" \
    || echo "comment patch failed (non-fatal)"
else
  gh api -X POST "repos/${GITHUB_REPOSITORY}/issues/${pr}/comments" -f body="$body" \
    || echo "comment post failed (non-fatal)"
fi
