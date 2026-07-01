# NFFT3 — Context

NFFT3 is a precision-agnostic C library for nonequispaced fast Fourier transforms.
This glossary fixes the vocabulary to use inside the project.

## Language

### Build systems & precision

**Coexistence**:
The policy that the CMake build is added *alongside* Autotools, never replacing it.
A change to one build system must not break the other; shared C sources,
`configure.ac`, and `Makefile.am` are not edited to serve CMake.
_Avoid_: migration, port (they imply replacement).

**Precision suffix**:
The single letter fixing a build tree's floating-point precision in the library
name and mangling: `` (double, `nfft_`), `f` (float, `nfftf_`), `l` (long-double,
`nfftl_`). One precision per build tree (the FFTW model).
_Avoid_: precision flavor, variant.

### Benchmarks & CodSpeed

**Benchmark name**:
The google_benchmark/CodSpeed identifier
`benchmarks/bench_<module>_direct.cpp::<prefix><function>/<args>`, one source file
per direct-transform module (`bench_nfft_direct.cpp`, `bench_nfct_direct.cpp`,
`bench_nfst_direct.cpp`), each mirroring the same 1D/2D/3D forward/adjoint
structure. The leaf (`<prefix><function>`) comes from the `BENCH` macro in
`benchmarks/util.h`; since codspeed-cpp v2.x the `benchmarks/bench_<module>_direct.cpp::`
prefix is supplied by codspeed from `__FILE__` relative to `CODSPEED_ROOT_DIR`
(pinned to the repo root in `benchmarks/CMakeLists.txt`). The unit CodSpeed tracks;
Preserving benchmark names for the same unit of work while making changes is important.
_Avoid_: benchmark id, test name, label.

**Benchmarks prefix** (`BENCHMARKS_PREFIX`):
The per-cell `<BUILD_CONFIG>/` string embedded in every benchmark name, where
`BUILD_CONFIG = <os>_<compiler>_<window>_<precision>_<openmp>`. The only
per-matrix-cell-variable part of a benchmark name; both build systems receive it
as an identical string.
_Avoid_: benchmark prefix tag, name prefix.

**CodSpeed profile**:
The full set of benchmark names CodSpeed has on record over time. The cutover must
leave it unchanged — no names added, none removed.
_Avoid_: benchmark suite, CodSpeed dashboard, baseline.

**Benchmark-producing cell**:
A CI matrix cell that builds and runs the benchmark, currently exactly the `kaiserbessel`
window cells. Gated by `BENCHMARK_AGNOSTIC_WINDOW`.
_Avoid_: benchmark job, active cell.

**Agnostic benchmark flags**:
The `window` / `openmp` / `precision` : `0/1` selectors choosing which matrix
cells build a given benchmark (Autotools `--with-agnostic-benchmarks`, CMake
`NFFT_AGNOSTIC_BENCHMARKS`; same `param:flag` CSV format). Only `window` is used
today.
_Avoid_: benchmark filters, build flags.

### The two parity guards (distinct — do not conflate)

**Config-symbol parity guard**:
`cmake/config-parity-check.sh` — diffs the macro symbol *set* in Autotools
`include/config.h` against CMake `build/config.h`, failing on any symbol the CMake
header omits. Guards against config drift between the two build systems. (As of
Phase 4 it exists but is not yet wired into CI.)
_Avoid_: the parity guard (ambiguous), config check.

### Tests

**ctest suite**:
The CMake counterpart of the Autotools `make check`: the CUnit binaries `checkall`
(serial) and `checkall_threads` (OpenMP) registered with `add_test` and run via
`ctest`. Covers exactly the **util / nfft / nfct / nfst** suites — the same set the
Autotools suite builds. In CI it runs in every matrix cell (full window × precision ×
openmp × compiler parity), like `make check`. Distinct from an **interface test**:
the term covers the CUnit core only, never the opt-in MATLAB/Octave runners.
_Avoid_: unit tests (ambiguous with the per-assertion `CU_ASSERT` cases), test harness.

**Interface test**:
An opt-in `ctest` case that runs a MATLAB/Octave unit-test runner
(`nfftUnitTestsRunAndExit` etc.) via the configured toolbox, mirroring the
`check_*_octave.sh`/`check_*_matlab.sh` tests Autotools runs under `make check`.
Registered only when a mex backend is configured (e.g. `NFFT_WITH_OCTAVE`); the
Octave set is nfft + nfsft (`HAVE_NFSFT`) + nfsoft (`HAVE_NFSOFT`). Separate from
the **ctest suite** (the CUnit core) — the two are never conflated.
_Avoid_: ctest suite (that is the CUnit core), unit test (the per-assertion cases).

**Orphaned test**:
`tests/check_nfsft.c` — a test source present in the tree but compiled by **neither**
build system. The CUnit suite has no tests for the six double-only modules (nnfft,
nsfft, mri, fpt, nfsft, nfsoft); the tested modules nfct/nfst build in every precision,
so "gate the six modules' tests to double precision" is a **no-op**.
_Avoid_: disabled test, skipped test (those imply a wired-up test deliberately not run).

**FindCUnit**:
The vendored `cmake/FindCUnit.cmake` — locates CUnit (the `<CUnit/CUnit.h>` header and
the `cunit` library) and provides the `CUnit::CUnit` imported target; the CMake
counterpart of `m4/nfft_lib_cunit.m4`. Honors `CUnit_ROOT` / `CUnit_INCLUDEDIR` /
`CUnit_LIBDIR`.
_Avoid_: the CUnit finder (ambiguous with the Autotools macro).

**Reference-data generator**:
The Python + mpmath tool `tests/refgen/` that is the single source of truth for the
file-based check data: it emits `tests/data/*.txt`, the generated C headers
`tests/data/generated/<module>_testcases.h`, and the data `EXTRA_DIST` list. Not a
build dependency. See ADR-0002.
_Avoid_: data scripts, Mathematica notebooks (removed).

**File-based check / Online check**:
The two CUnit test classes. A **file-based check** reads a high-precision reference
from `tests/data/*.txt` and validates both the direct and fast transform against it.
An **online check** generates random input, builds the reference with the direct
transform, and validates the fast transform against it. See
`docs/agents/test-methodology.md`.
_Avoid_: reference test / accuracy test (ambiguous — name the class).

### Accuracy tracking (HTML report)

**Accuracy metric**:
One row key in the aggregated accuracy JSON ("BMF") — the **max error over all
bound-absorbed parameter values** for a fixed combination of **error-shaping
parameters**. Rendered into the in-tree HTML report (never gates CI). Distinct from
a CodSpeed **Benchmark name**, which tracks instruction count, not error. See
`docs/agents/accuracy-tracking.md` and ADR-0004.
_Avoid_: accuracy benchmark, error benchmark (conflates with CodSpeed).

**Error-shaping parameter**:
A test parameter that changes a transform's *achievable* accuracy and so earns its
own **accuracy metric**: window, precision, runtime (serial vs OpenMP — the parallel
reduction order perturbs the low bits), dimension, transform kind (direct/fast,
forward/adjoint), precompute/init variant, and the file-vs-online oracle. Window and
precision are carried by the **accuracy testbed**; the rest are in the metric name.
_Avoid_: error parameter (ambiguous).

**Bound-absorbed parameter**:
A test parameter whose effect is already captured by the analytic error bound — the
bandwidth `N` and node count `M`. Collapsed via `max` within an **accuracy metric**,
never given its own series.
_Avoid_: size parameter.

**Accuracy digits**:
The primary accuracy measure, `-log10(max err)` — the worst-case number of accurate
digits (higher is better; a regression lowers it). Log-scaled so it reads cleanly
across the ~14 orders of magnitude the raw error spans. The raw `max err` is kept as
the secondary **max-error** measure for the exact figure, and `-log10(bound)` as
**bound-digits** so the heatmap can color by margin.
_Avoid_: tightness ratio (`err/bound`, the superseded measure), error ratio.

**Accuracy testbed**:
The per-cell `BUILD_CONFIG` (`<os>_<compiler>_<window>_<precision>`), one
`<testbed>.bmf.json` file. The only per-matrix-cell-variable part of an **accuracy
metric**; the metric name is identical across testbeds so `diff.py` joins PR-vs-
baseline figures by name (metric-name stability is what makes the diff align).
_Avoid_: cell, matrix testbed.

### Language interfaces

**Interface kernel**:
The `kernel/*.c` sources recompiled into a PIC static library (`nfft3_iface`, plus
`nfft3_iface_omp` under OpenMP) that links **no FFTW of its own** — the FFTW flavour
(system vs MATLAB-bundled) is chosen at the *final* shared-object link. One pair
serves every interface backend, collapsing Autotools' separate
`libnfft3${s}_julia.la` / `libnfft3${s}_matlab.la` convenience libs (libtool bakes
FFTW into each; CMake static libs defer linking).
_Avoid_: julia/matlab kernel variant (there is only one), kernel object lib.

**Mex backend**:
The toolbox a mex module is built against: **octave** (`NFFT_WITH_OCTAVE`, via the
vendored `FindOctave.cmake` with system FFTW or **matlab** (`find_package(Matlab)` 
+ bundled `libmwfftw3`. The two are mutually exclusive in one build tree (they produce 
the same `<mod>mex` files). The loadable produced is named exactly `<mod>mex.<ext>` 
(`.mex` for Octave), what the `.m` wrapper calls.
_Avoid_: mex flavor, octave/matlab variant (use "backend").

**Stub-MATLAB fixture**:
A fake MATLAB tree (`cmake/test-fixtures/fake-matlab/`: stub `extern/include` headers,
stub `bin/<arch>` libs, `libmwfftw3` = a copy of the real system FFTW, the
`mexFunction.map` version script) used to verify the MATLAB backend's **CMake wiring**
with no MATLAB present — that `find_package(Matlab)` + `matlab_add_mex` + the
`libmwfftw3` finder resolve and produce a `<mod>mex.<ext>` exporting only
`mexFunction`. It verifies **wiring only**, not that the C matches real MATLAB headers
(self-authored stub) and not runtime/numerical correctness.
_Avoid_: mock MATLAB (implies behavioural emulation — it emulates nothing).

**Provisional backend**:
A backend whose CMake wiring is verified (e.g. MATLAB via the **stub-MATLAB fixture**)
but which has not yet been built+run against the real toolbox. It stays provisional
until its runtime exit criterion passes (the MATLAB-host smoke test). Distinct from the
**Octave** backend, which is run-verified locally.
_Avoid_: working backend, supported backend (until the runtime test passes).

### Examples & applications

**Program target**:
An example (`examples/`) or application (`applications/`) executable. In Autotools
these are `noinst_PROGRAMS` — built by `make`, **never installed**, and **never run**
by `make check` (several block on `getc(stdin)`). CMake parity bar is therefore
that every such binary **builds**, not that it runs; no `ctest` cases are added for
them. The double-only ones are gated on the kernel `HAVE_*` variables, so they vanish
in float/long-double exactly as `m4/ax_nfft_module.m4` arranges.
_Avoid_: test program (implies it is run), installed binary.

**Orphaned application**:
`applications/iterS2` — present in the tree but built by **neither** build system
(commented out in `configure.ac`, no `Makefile.in`, `DIR_ITERS2` hard-wired empty).
Kept in the repo for possible later re-integration; the CMake build omits it to stay
at parity with Autotools (building it would break coexistence symmetry). The same
class as the orphaned test `tests/check_nfsft.c`. (`examples/mri/Makefile.am` is a
separate case — present but **empty**, so it builds nothing in either system.)
_Avoid_: disabled app, removed app (the source is retained, just unbuilt).
