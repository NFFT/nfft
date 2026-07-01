#!/usr/bin/env bash
# Fork-PR accuracy report, run from the trusted default branch via workflow_run. 
# No gh-pages write (read-only token / no contents:write), so the report link 
# points at the develop baseline dashboard.
set -uo pipefail
owner="${GITHUB_REPOSITORY%%/*}"; repo="${GITHUB_REPOSITORY##*/}"
baseline_raw="https://raw.githubusercontent.com/${owner}/${repo}/gh-pages/baseline"
report_url="https://${owner}.github.io/${repo}/"   # develop dashboard (not per-PR)

# Resolve the PR number from the build's head SHA.
pr="$(gh api "repos/${GITHUB_REPOSITORY}/commits/${HEAD_SHA}/pulls" \
  --jq 'first(.[].number) // empty')"
[ -n "$pr" ] || { echo "no open PR for ${HEAD_SHA}; nothing to do"; exit 0; }

# Download the fork build's BMF artifacts from the triggering run.
gh run download "${RUN_ID}" --pattern "accuracy-bmf-*" --dir bmf-artifacts \
  || { echo "no accuracy artifacts on run ${RUN_ID}; nothing to do"; exit 0; }

shopt -s nullglob
mkdir -p pr-bmf
for d in bmf-artifacts/accuracy-bmf-*/; do
  tb="$(basename "$d")"; tb="${tb#accuracy-bmf-}"
  [ -f "${d}accuracy.bmf.json" ] && cp "${d}accuracy.bmf.json" "pr-bmf/${tb}.bmf.json"
done

# nullglob-safe: empty array when no BMFs.
pr_files=(pr-bmf/*.bmf.json)
[ "${#pr_files[@]}" -gt 0 ] || { echo "no BMFs; nothing to do"; exit 0; }

# Fetch baseline.
mkdir -p base-bmf
for f in "${pr_files[@]}"; do
  tb="$(basename "$f")"
  curl -fsSL --remove-on-error "${baseline_raw}/${tb}" -o "base-bmf/${tb}" \
    || echo "no baseline yet: ${tb}"
done

# nullglob-safe: empty array when no BMFs.
base_files=(base-bmf/*.bmf.json)
if [ "${#base_files[@]}" -gt 0 ]; then
  uv run python -m tests.accuracy.pr_report pr-bmf base-bmf out --report-url "$report_url"
else
  uv run python -m tests.accuracy.pr_report pr-bmf pr-bmf out --no-baseline --report-url "$report_url"
fi

# Post the Check on the PR HEAD commit.
title="$(jq -r .title out/check.json)"; summary="$(jq -r .summary out/check.json)"
gh api -X POST "repos/${GITHUB_REPOSITORY}/check-runs" \
  -f name="accuracy" -f head_sha="${HEAD_SHA}" -f status="completed" \
  -f conclusion="neutral" -f output[title]="$title" -f output[summary]="$summary" \
  || echo "check post failed (non-fatal)"

# Label the link so a fork reviewer doesn't mistake the develop dashboard for the
# PR's own report.
body="$(cat out/comment.md)"
body="${body/📊 \[Full accuracy report\]/📊 \[develop baseline dashboard (per-PR report unavailable for forks)\]}"

# --paginate so the marker is found even on PRs with >30 comments.
existing="$(gh api --paginate "repos/${GITHUB_REPOSITORY}/issues/${pr}/comments" \
  --jq 'first(.[] | select(.body | startswith("<!-- nfft-accuracy-report -->")) | .id) // empty')"
if [ -n "$existing" ]; then
  gh api -X PATCH "repos/${GITHUB_REPOSITORY}/issues/comments/${existing}" -f body="$body" \
    || echo "comment patch failed (non-fatal)"
else
  gh api -X POST "repos/${GITHUB_REPOSITORY}/issues/${pr}/comments" -f body="$body" \
    || echo "comment post failed (non-fatal)"
fi
