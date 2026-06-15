# AGENTS.md

Guidance for coding agents working in this repository.

## 1. Project overview

**NFFT3** is a C library for computing **nonequispaced fast Fourier transforms**
and related transforms. It is a mature, research-grade numerical library.

Transforms implemented (each is its own module under `kernel/` with a public
API in `include/nfft3.h`):

- **NFFT** — nonequispaced FFT (forward + adjoint), the core module (`kernel/nfft`)
- **NNFFT** — nonequispaced in both time *and* frequency (`kernel/nnfft`)
- **NFCT / NFST** — real-valued (co)sine transforms (`kernel/nfct`, `kernel/nfst`)
- **NFSFT** — transform on the sphere S² (`kernel/nfsft`)
- **NFSOFT** — transform on the rotation group SO(3) (`kernel/nfsoft`)
- **NSFFT** — sparse/hyperbolic-cross FFT (`kernel/nsfft`)
- **FPT** — fast polynomial transform (`kernel/fpt`)
- **solver** — generalised inverse transforms via iterative methods (CGNR/CGNE)
- shared helpers in `kernel/util`

Application/example programs live in `applications/` (fastsum, fast Gauss
transform, MRI, polar FFT, Radon, S² quadrature, …) and `examples/`.

### Tech stack

- **Language:** C (C99). Benchmarks are C++17.
- **Build system:** GNU Autotools (autoconf/automake/libtool) + `make`.
  Experimental CMake support.
- **Hard dependency:** [FFTW3](https://fftw.org) (`libfftw3-dev`). Single,
  double, and long-double FFTW variants are all used depending on precision.
- **Tests:** [CUnit](http://cunit.sourceforge.net) (`libcunit1-dev`).
- **Benchmarks:** [CodSpeed](https://codspeed.io) C++ integration built on
  google_benchmark (CI-oriented; see §4).
- **Interfaces:** Julia (`julia/`), Matlab/Octave (`matlab/`).
- The dev container (`.devcontainer/`) is Ubuntu 24.04 with gcc-14, clang,
  FFTW (all precisions), CUnit, Octave, Julia, autotools, and clang-format
  preinstalled. Assume these are available on the local host.

### Code conventions (important)

The library is **precision-agnostic**: the same C sources compile in float,
double, or long-double precision. Read `CONVENTIONS.md` before editing C code.
Key rules:

- **Name mangling:** use `Y(foo)` for names exported library-wide (expands to
  `nfftf_`/`nfft_`/`nfftl_` by precision), `X(foo)` for module-local names
  (e.g. `nfct_foo` inside the NFCT module), and `FFTW(foo)` for FFTW names.
  Never hard-code a `nfft_`/`nfftf_`/`nfftl_` prefix.
- Use the type aliases `R` (real), `E` (real, possibly extended precision),
  `C` (complex); `A(...)` for assert, `CK(...)` for check.
- Indentation is 2 spaces, BSD brace style. A `.clang-format` is provided.
- New code must keep the float/double/long-double build matrix working.

## 2. Building (local host)

This is a source checkout, so the configure script must be generated first.

```bash
# One-time (and after editing any Makefile.am / configure.ac):
./bootstrap.sh                      # runs libtoolize + autoreconf

# Configure. Build ALL modules, tests, examples, and applications:
./configure --enable-all --enable-openmp

# Compile:
make -j
```

Useful `./configure` flags:

- `--enable-all` — build all transform modules (plus examples & applications).
- `--enable-openmp` — multithreaded (OpenMP) build; also produces the
  `*_omp` library and the `checkall_threads` test binary.
- `--enable-float` / `--enable-long-double` — change precision (default is
  double). These are mutually exclusive. Build each precision in a *separate*
  configured tree.
- `--with-window=ARG` — window function: `kaiserbessel` (default), `gaussian`,
  `bspline`, `sinc`, `dirac`.
- `--enable-exhaustive-unit-tests` — larger, slower test cases (CI uses this).
- `--enable-julia` — build the Julia interface (double precision + shared libs only).
- `--with-matlab=/path` / `--with-octave=/path` — build the MATLAB/Octave interface.
- `--enable-benchmarks --with-codspeed=<path>` — **legacy** Autotools benchmark
  build; benchmarks are now built with CMake (see §4).
- `./configure --help` lists everything.

A representative "build everything" configure line matching CI:

```bash
./configure --with-window=kaiserbessel --enable-all \
            --enable-exhaustive-unit-tests --enable-openmp --enable-julia
make -j
```

The build produces `libnfft3<suffix>.la` (and `..._omp.la` with OpenMP) in the
top directory, where `<suffix>` is empty for double, `f` for float, `l` for
long-double.

## 3. Running the tests

Tests use **CUnit** and only build when CUnit is detected at configure time.

```bash
make check          # builds + runs the test suite, prints PASS/FAIL summary
```

What this runs:

- `tests/checkall` — the full single-threaded suite.
- `tests/checkall_threads` — the same suite against the OpenMP library
  (only when configured with `--enable-openmp`).

Checking which tests pass/fail:

- The `make check` console output gives the per-test and overall pass/fail
  summary, and writes `tests/checkall.log` / `tests/checkall.trs` (automake
  test logs).
- CUnit also emits a detailed XML report at
  `tests/CUnitAutomated-Results.xml` (and `..._threads-Results.xml`). This can
  be converted to JUnit with the bundled stylesheet `tests/cunit2junit.xsl`.
- To re-run a single binary directly and capture full output:

  ```bash
  tests/checkall            # or: tests/checkall_threads
  ```

The C tests (`tests/nfft.c`, `nfct.c`, `nfst.c`, `check_nfsft.c`, …) validate
transforms against reference data in `tests/data/` and against the direct
(slow) transforms. Use `--enable-exhaustive-unit-tests` for the thorough set.

## 4. Running the benchmarks

Benchmarks (`benchmarks/bench_nfft_direct.cpp`) use the **CodSpeed** C++ integration
and run in CI for continuous regression tracking. Build them with the **CMake** build
(self-contained: `FetchContent` fetches/builds codspeed-cpp, submodules and all). A
single knob, `-DNFFT_BENCHMARK_MODE=`, both enables the benchmark build and picks the
measurement mode — **baked in at build time**:

- `off` (default) — benchmarks not built.
- `simulation` — deterministic instruction count; the metric CI gates on. The binary
  only measures under callgrind (`valgrind --tool=callgrind …`, read `I refs`); run
  directly it is inert (*"unknown environment"*).
- `walltime` — wall-clock timing; the binary writes a local stats JSON to
  `$CODSPEED_PROFILE_FOLDER/results/<pid>.json` (offline, no valgrind/runner).

```bash
cmake -S . -B build-cmake -DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON \
      -DCMAKE_C_FLAGS="-O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math"
cmake --build build-cmake -j        # binaries in build-cmake/benchmarks/
CODSPEED_PROFILE_FOLDER=/tmp/wt build-cmake/benchmarks/bench_nfft_direct
```

Switching mode = reconfigure the tree with a different `-DNFFT_BENCHMARK_MODE` (or use
a second build dir). Other options: `-DBENCHMARKS_PREFIX=`, `-DNFFT_AGNOSTIC_BENCHMARKS=`
(`"window:1,openmp:0,precision:1"`), `-DFETCHCONTENT_SOURCE_DIR_CODSPEED=<path>` (reuse
a checkout / offline). `valgrind` and the `codspeed` CLI are preinstalled in the dev
container.

Continuous tracking happens via CodSpeed in GitHub Actions (same CMake build,
`.github/workflows/build-linux.yml`). For the agent optimization workflow built on
these benchmarks, see
[`docs/agents/performance-optimization-loop.md`](docs/agents/performance-optimization-loop.md).

> **Legacy (do not use):** the Autotools benchmark path (`./configure
> --enable-benchmarks --with-codspeed=<path>` + `make bench`) predates CodSpeed; it
> needs a hand-built `codspeed-cpp` and prints no measurements. Use the CMake build.

## Quick reference

| Task | Command |
|------|---------|
| Regenerate build system | `./bootstrap.sh` |
| Configure (everything) | `./configure --enable-all --enable-openmp` |
| Build | `make -j` |
| Run tests | `make check` |
| Test results (detail) | `tests/checkall`, `tests/CUnitAutomated-Results.xml` |
| Build benchmarks | `cmake -S . -B build-cmake -DNFFT_BENCHMARK_MODE=walltime && cmake --build build-cmake -j` (§4) |
| Measure (walltime) | `CODSPEED_PROFILE_FOLDER=/tmp/wt build-cmake/benchmarks/bench_nfft_direct` → `/tmp/wt/results/*.json` |
| Clean | `make clean` / full reset: `make distclean` |
| Format C code | `clang-format -i <file>` (uses repo `.clang-format`) |

CI mirrors these steps in `.github/workflows/build-linux.yml` and
`build-macos.yml` across the full window × precision × OpenMP matrix — consult
those files when reproducing a CI configuration exactly.

## Agent skills

### Issue tracker

Issues and PRDs live as markdown files under `.scratch/<feature>/`. See `docs/agents/issue-tracker.md`.

### Triage labels

Five canonical triage roles using their default strings (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`). See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: one `CONTEXT.md` + `docs/adr/` at the repo root. See `docs/agents/domain.md`.

