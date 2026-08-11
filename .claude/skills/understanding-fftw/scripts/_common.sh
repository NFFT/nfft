# Shared setup for fetch-source.sh / fetch-docs.sh.
# Resolves the current project's root and reference-cache directory, and
# makes sure that directory is excluded from git locally (via
# .git/info/exclude) without ever touching a tracked .gitignore file.
set -euo pipefail

PROJECT_ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || {
  echo "Error: not inside a git repository. cd into the target project first." >&2
  exit 1
}

REF_DIR="$PROJECT_ROOT/.fftw3-reference"
EXCLUDE_FILE="$PROJECT_ROOT/.git/info/exclude"
EXCLUDE_PATTERN="/.fftw3-reference/"

mkdir -p "$REF_DIR"

if [ -f "$EXCLUDE_FILE" ] && grep -qxF "$EXCLUDE_PATTERN" "$EXCLUDE_FILE"; then
  : # already excluded
else
  echo "$EXCLUDE_PATTERN" >> "$EXCLUDE_FILE"
  echo "Added '$EXCLUDE_PATTERN' to $EXCLUDE_FILE (local-only, not tracked)"
fi
