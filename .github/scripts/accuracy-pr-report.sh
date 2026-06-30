#!/usr/bin/env bash
# Post the accuracy Check + upserted comment for a same-repo PR, and archive its
# heatmaps to the gh-pages branch (permanent; linked via raw.githubusercontent
# URLs). Fork PRs never reach this script (the accuracy-report job is skipped for
# them — GitHub gives fork pull_request runs a read-only token). The job step runs
# with continue-on-error, so any failure here is informational and never reds CI.
set -uo pipefail
owner="${GITHUB_REPOSITORY%%/*}"; repo="${GITHUB_REPOSITORY##*/}"
baseline_raw="https://raw.githubusercontent.com/${owner}/${repo}/gh-pages/baseline"
images_raw="https://raw.githubusercontent.com/${owner}/${repo}/gh-pages/pr"
pr="$(jq -r .number "$GITHUB_EVENT_PATH")"
abs_url="${images_raw}/${pr}/absolute.png"
rel_url="${images_raw}/${pr}/relative.png"

run_report() { uv run --with matplotlib python tests/bench/pr_report.py "$@"; }
publish() {  # publish given PNGs to gh-pages under pr/<n>/
  mkdir -p "site/pr/${pr}"; cp "$@" "site/pr/${pr}/"
  bash .github/scripts/gh-pages-publish.sh "site/pr/${pr}" "pr ${pr} heatmaps" "pr/${pr}"
}

# Fetch baseline BMFs (skip cleanly if develop has not published a baseline yet).
shopt -s nullglob
mkdir -p base-bmf
for f in pr-bmf/*.bmf.json; do
  tb="$(basename "$f")"
  curl -fsSL --remove-on-error "${baseline_raw}/${tb}" -o "base-bmf/${tb}" \
    || echo "no baseline yet: ${tb}"
done

# nullglob-safe existence check: an unmatched glob is an EMPTY array (NOT a
# literal), whereas `ls base-bmf/*.bmf.json` under nullglob would run `ls` with no
# args, list the cwd, and falsely succeed -> always diff vs an empty baseline.
base_files=(base-bmf/*.bmf.json)
if [ "${#base_files[@]}" -gt 0 ]; then
  run_report pr-bmf base-bmf out --abs-url "$abs_url" --rel-url "$rel_url"
  publish out/absolute.png out/relative.png
else
  # No develop baseline yet (e.g. the first PR): absolute-only, no misleading diff.
  run_report pr-bmf pr-bmf out --no-baseline --abs-url "$abs_url"
  publish out/absolute.png
fi

# Post the Check on the PR HEAD commit (not the ephemeral merge SHA). Never fails CI.
title="$(jq -r .title out/check.json)"; summary="$(jq -r .summary out/check.json)"
gh api -X POST "repos/${GITHUB_REPOSITORY}/check-runs" \
  -f name="accuracy" -f head_sha="${HEAD_SHA}" -f status="completed" \
  -f conclusion="neutral" -f output[title]="$title" -f output[summary]="$summary" \
  || echo "check post failed (non-fatal)"

# Upsert the comment (find by marker, else create).
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
