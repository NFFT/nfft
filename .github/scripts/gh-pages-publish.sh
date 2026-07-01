#!/usr/bin/env bash
# Usage: gh-pages-publish.sh <src_dir> <commit_msg> [<dest_subdir>]
set -euo pipefail
src="$1"; msg="$2"; dest="${3:-}"
repo_url="https://x-access-token:${GH_TOKEN}@github.com/${GITHUB_REPOSITORY}.git"
tmp="$(mktemp -d)"
git clone --depth 1 --branch gh-pages "$repo_url" "$tmp" 2>/dev/null || {
  git clone --depth 1 "$repo_url" "$tmp"; git -C "$tmp" checkout --orphan gh-pages
  git -C "$tmp" rm -rf . >/dev/null 2>&1 || true
}
mkdir -p "$tmp/$dest"
cp -r "$src"/. "$tmp/$dest/"
git -C "$tmp" add -A
git -C "$tmp" -c user.name=ci -c user.email=ci@nfft commit -m "$msg" || { echo "no changes"; exit 0; }
git -C "$tmp" push origin gh-pages
