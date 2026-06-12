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
- `--enable-benchmarks --with-codspeed=<path>` — see §4.
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

Benchmarks (`benchmarks/bench_nfft_direct.cpp`) use the **CodSpeed** C++
integration library and are primarily intended to run in CI. To build/run them
locally you must first build `codspeed-cpp` from source.

```bash
# 1. Build the CodSpeed integration library (one-time):
git clone --depth 1 https://github.com/CodSpeedHQ/codspeed-cpp
cmake -S codspeed-cpp/google_benchmark -B codspeed-cpp/google_benchmark \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBENCHMARK_DOWNLOAD_DEPENDENCIES=ON \
      -DCODSPEED_MODE=simulation -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_EXTENSIONS=OFF \
      -DCMAKE_CXX_FLAGS="-DCODSPEED_ROOT_DIR=\\\"$(pwd)\\\""
cmake --build codspeed-cpp/google_benchmark

# 2. Configure NFFT with benchmarks enabled:
./configure --enable-all --enable-benchmarks \
            --with-codspeed="$(pwd)/codspeed-cpp"
make -j

# 3. Build and run all benchmarks:
cd benchmarks
make bench          # builds each benchmark, then runs it, printing results
```

Other benchmark configure flags:

- `--with-benchmarks-prefix=PREFIX` — prefix applied to benchmark names.
- `--with-agnostic-benchmarks="window:1,openmp:0,precision:1"` — comma-separated
  `param:flag` pairs selecting which benchmark variants get built.

Benchmark binaries are placed in `benchmarks/`. The built binaries are plain
google_benchmark executables, so you can also run one directly to see timing
output (e.g. `./benchmarks/bench_nfft_direct`). Continuous performance tracking
(regression detection) happens via CodSpeed in GitHub Actions, not locally.

## Quick reference

| Task | Command |
|------|---------|
| Regenerate build system | `./bootstrap.sh` |
| Configure (everything) | `./configure --enable-all --enable-openmp` |
| Build | `make -j` |
| Run tests | `make check` |
| Test results (detail) | `tests/checkall`, `tests/CUnitAutomated-Results.xml` |
| Build + run benchmarks | `./configure --enable-benchmarks --with-codspeed=<path>` then `cd benchmarks && make bench` |
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

