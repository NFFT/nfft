#!/usr/bin/env bash
# Refreshes both the local FFTW3 source mirror and the converted manual.
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"$DIR/fetch-source.sh"
"$DIR/fetch-docs.sh"
