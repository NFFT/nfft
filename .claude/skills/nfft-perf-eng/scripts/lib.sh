# Shared helpers for the nfft-perf-eng deterministic scripts.
# Source this (`. lib.sh`); do not execute it directly.
#
# The precision matrix (float · double · long double) and the CI-aligned build flags
# live here so every script agrees on them. See ../details/precision-matrix.md and
# ../details/phase-a-baseline.md.

# Precision suffixes, in canonical order. One CMake tree per precision.
PERF_PRECISIONS=(d f l)

# perf_tree <d|f|l> -> the build-tree directory name for that precision.
perf_tree() {
  case "$1" in
    d) echo build-cmake ;;
    f) echo build-cmake-f ;;
    l) echo build-cmake-l ;;
    *) echo "perf_tree: bad precision '$1' (want d|f|l)" >&2; return 2 ;;
  esac
}

# perf_prec_name <d|f|l> -> human-readable precision name.
perf_prec_name() {
  case "$1" in
    d) echo "double" ;;
    f) echo "float" ;;
    l) echo "long double" ;;
    *) echo "perf_prec_name: bad precision '$1'" >&2; return 2 ;;
  esac
}

# CI-aligned C flags. These MUST match .github/workflows/bench-linux.yml — the
# -falign-functions=64 -falign-loops=32 pair pins code layout so an untouched
# neighbour can't show a phantom alignment regression (see ../details/caveats.md).
perf_cflags() {
  echo "-O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math -falign-functions=64 -falign-loops=32"
}

# Repository root (scripts are always run from inside the repo).
perf_root() { git rev-parse --show-toplevel; }
