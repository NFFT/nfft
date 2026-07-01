#!/usr/bin/env bash
# Post the accuracy Check + upserted comment for a same-repo PR, and publish its
# HTML report to the gh-pages branch under pr/<n>/ (served by GitHub Pages).
set -uo pipefail
owner="${GITHUB_REPOSITORY%%/*}"; repo="${GITHUB_REPOSITORY##*/}"
baseline_raw="https://raw.githubusercontent.com/${owner}/${repo}/gh-pages/baseline"
pr="$(jq -r .number "$GITHUB_EVENT_PATH")"
report_url="https://${owner}.github.io/${repo}/pr/${pr}/"

run_report() { uv run python -m tests.accuracy.pr_report "$@"; }

# Fetch baseline BMFs.
shopt -s nullglob
mkdir -p base-bmf
for f in pr-bmf/*.bmf.json; do
  tb="$(basename "$f")"
  curl -fsSL --remove-on-error "${baseline_raw}/${tb}" -o "base-bmf/${tb}" \
    || echo "no baseline yet: ${tb}"
done

# nullglob-safe: empty array when no BMFs.
base_files=(base-bmf/*.bmf.json)
if [ "${#base_files[@]}" -gt 0 ]; then
  run_report pr-bmf base-bmf out --report-url "$report_url"
else
  # No develop baseline yet (e.g. the first PR): absolute-only, no misleading diff.
  run_report pr-bmf pr-bmf out --no-baseline --report-url "$report_url"
fi

# Publish the HTML report to gh-pages/pr/<n>/. Track whether it succeeded so the
# comment can degrade gracefully if Pages will lag / the publish failed.
published=1
mkdir -p "site/pr/${pr}"; cp out/index.html "site/pr/${pr}/"
bash .github/scripts/gh-pages-publish.sh "site/pr/${pr}" "pr ${pr} report" "pr/${pr}" \
  || { published=0; echo "gh-pages publish failed (non-fatal)"; }

# Post the Check on the PR HEAD commit.
title="$(jq -r .title out/check.json)"; summary="$(jq -r .summary out/check.json)"
gh api -X POST "repos/${GITHUB_REPOSITORY}/check-runs" \
  -f name="accuracy" -f head_sha="${HEAD_SHA}" -f status="completed" \
  -f conclusion="neutral" -f output[title]="$title" -f output[summary]="$summary" \
  || echo "check post failed (non-fatal)"

# Comment body; if the publish failed, note that the link may briefly lag.
body="$(cat out/comment.md)"
if [ "$published" -eq 0 ]; then
  body="${body}"$'\n\n_(report is still publishing — the link may 404 for a minute.)_'
fi

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
