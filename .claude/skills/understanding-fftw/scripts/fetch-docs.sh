#!/usr/bin/env bash
# Downloads the FFTW3 reference manual PDF and converts it to
# <project-root>/.fftw3-reference/fftw3-docs.md so it can be
# grepped/read offline. Must be run from inside the target project's git repo.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_common.sh"

PDF="$REF_DIR/fftw3.pdf"
OUT="$REF_DIR/fftw3-docs.md"
PDF_URL="https://fftw.org/fftw3.pdf"

echo "Downloading $PDF_URL"
curl -fsSL -o "$PDF" "$PDF_URL"

echo "Converting to Markdown ($OUT)"
if command -v pandoc >/dev/null 2>&1; then
  pandoc "$PDF" -f pdf -t gfm -o "$OUT"
elif command -v uvx >/dev/null 2>&1; then
  # pymupdf4llm produces clean, well-structured Markdown (headings, code
  # blocks) from the manual; this is the preferred fallback when pandoc
  # is unavailable, as it is on the devcontainer.
  uvx --from pymupdf4llm python3 - "$PDF" "$OUT" <<'PYEOF'
import sys, pathlib, pymupdf4llm
pdf_path, out_path = sys.argv[1], sys.argv[2]
md = pymupdf4llm.to_markdown(pdf_path)
pathlib.Path(out_path).write_text(md)
print(f"wrote {len(md)} chars")
PYEOF
elif command -v pdftotext >/dev/null 2>&1; then
  echo "Warning: pandoc/uvx unavailable; falling back to pdftotext (plain text, no headings/tables)" >&2
  pdftotext -layout "$PDF" "$OUT"
else
  echo "Error: no PDF-to-text/markdown converter available." >&2
  echo "Install one of: pandoc, uv (https://docs.astral.sh/uv/), pdftotext (poppler-utils)." >&2
  rm -f "$PDF"
  exit 1
fi

rm -f "$PDF"
echo "fftw3-docs.md: $(wc -l < "$OUT") lines"
