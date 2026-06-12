# Adopt a CMake build alongside Autotools

**Status:** accepted

We will add a CMake build system that **coexists** with the existing Autotools
build rather than replacing it, because Autotools cannot easily deliver four
things we need: stable CodSpeed benchmark CI, Windows/MSVC builds,
`find_package()` consumption, and good IDE integration. The precipitating pain is
**intermittent CodSpeed CI failures**: CodSpeed's `codspeed-cpp` (a google_benchmark
fork) is built by its own CMake, and our Autotools build currently reaches into
that build tree with hardcoded include/library paths and hand-set instrumentation
defines (`m4/nfft_lib_codspeed.m4`). That cross-build-system mismatch is the
likely cause of the flakiness, and CodSpeed only documents a CMake (`FetchContent`)
integration path.

## Key decisions

- **Coexistence, not replacement.** Autotools remains the supported build,
  including for legacy/HPC toolchains. This lets the CMake floor be modern
  (`cmake_minimum_required(VERSION 3.20)`) without locking anyone out.
- **One precision per build tree** (the FFTW model). Precision is fixed at
  configure time via `NFFT_SINGLE`/`NFFT_LDOUBLE`; an `ENABLE_FLOAT`/
  `ENABLE_LONG_DOUBLE` → `PREC_SUFFIX` scheme drives the library name. The public
  header `nfft3.h` is precision-agnostic and installed once; only the compiled
  library and its `NFFT3{,f,l}Config.cmake` are per-precision.
- **Benchmarks-first sequencing.** The first CMake milestone is the core library
  plus the CodSpeed benchmark (double precision only, with the precision
  scaffolding in place), obtaining `codspeed-cpp` via `FetchContent` (pinned tag,
  with a cache/override variable). This delivers the actual driver — stable
  benchmark CI — in the smallest slice. Tests, examples/applications, install/
  export, and the MATLAB/Julia interfaces follow.
- **Explicit source lists, hybrid layout.** Sources are listed explicitly (only
  ~33 stable `kernel/*.c` files; transform modules are single-file), not globbed.
  A root `CMakeLists.txt` plus `add_subdirectory` per top-level area co-locates
  each list with its `Makefile.am`.
- **Vendored `cmake/FindFFTW3.cmake`.** Handles the per-precision base/`_omp`/
  `_threads` FFTW triplet that off-the-shelf modules don't model.

## Consequences

- **Dual maintenance.** Build logic exists twice. The source of truth for config
  symbols is the `AC_DEFINE` set in `configure.ac` (`config.h.in` itself is
  `autoheader`-generated and not checked in); `cmake.config.h.in` must mirror it.
  A **CI parity guard** diffs the macro symbol sets between the two and fails on
  divergence to prevent silent drift.
- CMake CI starts as a single lane (double, all modules, OpenMP) and expands to
  the full window × precision × OpenMP matrix only in the final phase.
