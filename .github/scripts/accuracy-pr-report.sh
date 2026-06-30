#!/usr/bin/env bash
# Post the accuracy Check + upserted comment for a PR. Same-repo PRs also archive
# their heatmaps to the gh-pages branch (permanent); the comment links them via
# raw.githubusercontent URLs (stable, render inline, work without Pages enabled).
# Fork PRs get no gh-pages write (read-only token) -> emoji/text only.
set -euo pipefail
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
mkdir -p base-bmf
for f in pr-bmf/*.bmf.json; do
  tb="$(basename "$f")"
  curl -fsSL "${baseline_raw}/${tb}" -o "base-bmf/${tb}" || echo "no baseline yet: ${tb}"
done

if ls base-bmf/*.bmf.json >/dev/null 2>&1; then
  if [ "${IS_FORK:-false}" = "true" ]; then
    run_report pr-bmf base-bmf out                       # emoji only, no archive
  else
    run_report pr-bmf base-bmf out --abs-url "$abs_url" --rel-url "$rel_url"
    publish out/absolute.png out/relative.png
  fi
else
  # No develop baseline yet (e.g. the first PR): absolute-only, no misleading diff.
  if [ "${IS_FORK:-false}" = "true" ]; then
    run_report pr-bmf pr-bmf out --no-baseline           # text only, no archive
  else
    run_report pr-bmf pr-bmf out --no-baseline --abs-url "$abs_url"
    publish out/absolute.png
  fi
fi

# Post the Check (never fails CI).
title="$(jq -r .title out/check.json)"; summary="$(jq -r .summary out/check.json)"
gh api -X POST "repos/${GITHUB_REPOSITORY}/check-runs" \
  -f name="accuracy" -f head_sha="${GITHUB_SHA}" -f status="completed" \
  -f conclusion="neutral" -f output[title]="$title" -f output[summary]="$summary"

# Upsert the comment (find by marker, else create).
body="$(cat out/comment.md)"
existing="$(gh api "repos/${GITHUB_REPOSITORY}/issues/${pr}/comments" \
  --jq '.[] | select(.body | startswith("<!-- nfft-accuracy-report -->")) | .id' | head -1)"
if [ -n "$existing" ]; then
  gh api -X PATCH "repos/${GITHUB_REPOSITORY}/issues/comments/${existing}" -f body="$body"
else
  gh api -X POST "repos/${GITHUB_REPOSITORY}/issues/${pr}/comments" -f body="$body"
fi
