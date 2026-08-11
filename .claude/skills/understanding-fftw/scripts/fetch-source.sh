#!/usr/bin/env bash
# Shallow-clones (or updates) the upstream FFTW3 source tree into
# <project-root>/.fftw3-reference/fftw3-src so it can be browsed
# offline. Must be run from inside the target project's git repo.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_common.sh"

SRC_DIR="$REF_DIR/fftw3-src"
REPO_URL="https://github.com/FFTW/fftw3.git"

if [ -d "$SRC_DIR/.git" ]; then
  echo "Updating existing clone at $SRC_DIR"
  git -C "$SRC_DIR" fetch --depth 1 origin
  git -C "$SRC_DIR" reset --hard FETCH_HEAD
else
  echo "Cloning $REPO_URL (shallow, depth=1) into $SRC_DIR"
  git clone --depth 1 "$REPO_URL" "$SRC_DIR"
fi

commit="$(git -C "$SRC_DIR" rev-parse --short HEAD)"
date="$(git -C "$SRC_DIR" log -1 --format=%cd --date=short)"
echo "fftw3-src is now at commit $commit ($date)"
