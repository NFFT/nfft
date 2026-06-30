#!/usr/bin/env bash
set -euo pipefail
owner="${GITHUB_REPOSITORY%%/*}"; repo="${GITHUB_REPOSITORY##*/}"
pages_raw="https://raw.githubusercontent.com/${owner}/${repo}/gh-pages/baseline"
pr="$(jq -r .number "$GITHUB_EVENT_PATH")"

# Fetch baseline BMFs (skip cleanly if the dashboard hasn't published yet).
mkdir -p base-bmf
for f in pr-bmf/*.bmf.json; do
  tb="$(basename "$f")"
  curl -fsSL "${pages_raw}/${tb}" -o "base-bmf/${tb}" || echo "no baseline yet: ${tb}"
done

if [ -d base-bmf ] && ls base-bmf/*.bmf.json >/dev/null 2>&1; then
  if [ "${IS_FORK:-false}" = "true" ]; then
    uv run --with matplotlib python tests/bench/pr_report.py pr-bmf base-bmf out
  else
    abs="https://${owner}.github.io/${repo}/pr/${pr}/absolute.png"
    rel="https://${owner}.github.io/${repo}/pr/${pr}/relative.png"
    uv run --with matplotlib python tests/bench/pr_report.py pr-bmf base-bmf out \
      --abs-url "$abs" --rel-url "$rel"
    mkdir -p site/pr/${pr}; cp out/absolute.png out/relative.png site/pr/${pr}/
    bash .github/scripts/gh-pages-publish.sh site/pr/${pr} "pr ${pr} heatmaps" "pr/${pr}"
  fi
else
  # No baseline published yet: still render PR-only heatmaps + a flat-ish comment.
  uv run --with matplotlib python tests/bench/pr_report.py pr-bmf pr-bmf out
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
