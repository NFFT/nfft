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
header omits. Guards against config drift between the two build systems.
_Avoid_: the parity guard (ambiguous), config check.

### Tests

**ctest suite**:
The CMake counterpart of the Autotools `make check`: the CUnit binaries `checkall`
(serial) and `checkall_threads` (OpenMP) registered with `add_test` and run via
`ctest`. Covers exactly the **util / planner / nplan / nfft / nfct / nfst** suites
— the same set the Autotools suite builds. In CI it runs in every matrix cell (full window ×
precision × openmp × compiler parity), like `make check`. Distinct from an **interface
test**: the term covers the CUnit core only, never the opt-in MATLAB/Octave runners.
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

**wisdom**:
The planner's memo: for a problem size class (keyed by a 128-bit MD5 signature),
which solver won a search conducted under given impatience bounds. A cache, never a
source of truth — any mismatch on lookup or import degrades to a miss, never a wrong
plan. Lives in `kernel/planner/`; exported/imported through the public
`nfft_export_wisdom_*` / `nfft_import_wisdom_*` functions.

**timelimit** (planner):
A global wall-clock budget on the planner's *measured* candidate race, set by
`nfft_set_timelimit(seconds)` (`< 0` = unlimited, the default). Before measuring
each candidate the guru checks the elapsed coarse wall clock; on expiry it stops.
If no candidate finished before expiry (or there is no usable clock), the
selection degrades to **estimate grade** (unblessed) — so a timelimit-induced
loser is never blessed. ESTIMATE planning is never bounded. The timelimit is a
runtime bound on `struct planner_s`, *not* part of the wisdom key (the reserved
`timelimit_imp` impatience field stays 0).

**blessing**:
Marking a wisdom entry (`PLNR_BLESSING`) as worth persisting: it lives in the
planner's *blessed* hash table, is exported to wisdom files, and survives
`PLNR_FORGET_UNBLESSED`. An entry without this flag is *unblessed* — session-only
search memoisation, dropped on forget.

**subsumption**:
The impatience-lattice rule by which one wisdom entry answers another query's
bounds — `LEQ(a.u, q.u) && LEQ(q.l, a.l)` for a feasible entry `a` against query
`q` (with a mirrored rule for infeasible entries). Lookup returns the subsuming
entry with the smallest upper bound; insert kills every existing entry the
newcomer subsumes.

**solver registry**:
The planner's per-problem-kind table of registered solver recipes
(`kernel/planner/solvtab.c`, `solver.c`), keyed by registrar name + per-registrar
ordinal id. Wisdom import resolves a serialised `(name, id)` pair back to a solver
index through this table; an unresolvable pair is treated as malformed wisdom.

**provenance** (planner):
`kernel/planner/` follows the architecture of FFTW3's planner (problem / solver /
plan, wisdom keyed by an MD5 signature, printer/scanner, timer). The
implementation was written for this repository; FFTW3 is GPL-2.0-or-later like
NFFT3. `include/cycle.h` is FFTW's MIT-licensed file (see next entry).

**cycle.h provenance** (planner timing):
`include/cycle.h` is a vendored MIT-licensed copy of FFTW's (Frigo / MIT) — the
**one FFTW-derived file** the planner timing path depends on, kept under FFTW's
copyright header. The NFFT copy adds the aarch64 `CNTVCT_EL0`/`PMCCNTR_EL0`
tick-counter stanza and a guarded `<stdint.h>` include.

**problem / solver / plan** (planner trinity):
The planner's core (`kernel/planner/{problem,plan}.c`, search in `planner.c`):
a *problem* is a hashable description of work, a *solver* a registered recipe
whose `mkplan` returns an executable *plan* (with `pcost`) or NULL when
inapplicable; `planner_mkplan` consults wisdom, else keeps the cheapest
candidate and memoises the outcome unblessed. `pcost` is an analytic estimate
in the estimate-mode search; in a measured race it is in
**arbitrary tick units** (tick mode, when a cycle counter is present) or wall
seconds (slow-timer fallback), compared only within one race — measured and
estimate costs are never compared across modes.

**plan_ng** (NFFT bundled plan):
The next-generation NFFT lifecycle (`kernel/nfft/plan.c`): one bundle owns
one forward problem and one planner-selected plan that serves both directions
via `apply_adjoint`. Every plan (fast or direct) is coreless and
planner-native; there is no shared or per-plan legacy core. Naming rule:
`_ng` marks colliding nouns only; lifecycle verbs stay bare (`nfft_execute`,
`nfft_precompute`, ...). The guru + verbs + `nfft_fprint_plan` + planning
flags + wisdom I/O are declared in `nfft3.h`; the wisdom/print/timelimit
wrappers live in `kernel/nfft/api.c`. The guru takes only `fftw_flags` +
`planning` (no `nfft_flags`); `x`/`f_hat`/`f` are required in every mode,
`x` is copied into the plan, `f_hat`/`f` are borrowed; there are no
`plan_ng_x`/`plan_ng_f_hat`/`plan_ng_f` accessors.

**plan_ng (NFCT/NFST)** — *deferred*. The real-valued kinds were prototyped
and then removed from the new planner API before review; the new API is
NDFT/NFFT-only for now. The legacy NFCT/NFST API (`nfct_*`/`nfst_*`) is unaffected.

**fast-path geometry guard** (`guards_ok`, `kernel/nfft/nfft-nd.c`):
The applicability check the planner-native fast NFFT solver runs per axis
before it competes with the direct NDFT solvers. It accepts every geometry the
direct solvers accept — odd `N`, and per-axis type-I/type-II (mixed across
axes) — and additionally requires oversampling `sigma = n/N > 1` strictly on
every axis; at `sigma == 1` there is no zero-pad for DECONV to deconvolve into,
and for even `N` the window's `phi_hut` is exactly zero at the band edge. A
geometry the guard declines for a size reason still gets a correct answer from
the direct NDFT, just without the fast path.

`sigma <= 1` is different: `nfft_plan_ng_guru` **rejects** it and returns NULL
rather than quietly planning a direct transform, since the guru does not tell
the caller which solver won and a silent fallback would cost O(N*M) in place of
O(n log n). The rejection applies only while the fast solver is in play —
`NFFT_NO_FAST_NATIVE` takes it out of the running, so nothing is lost
unintentionally and `sigma <= 1` becomes legal. Unit axes are exempt, being
elided. The guard keeps its own check for callers who build a problem through
`mkproblem_nfft` below the guru. Wisdom consequence: for an odd or type-II shape
blessed before this guard existed, the wisdom key already hashed `variant[]`
and carried `N` parity, so the entry stays valid — it just still names a direct
solver and stays correct, only slower. Call `nfft_forget_wisdom()` (or delete
the wisdom file) to let such a shape re-plan onto the fast path.
_Avoid_: even-only fast path (no longer true).

**unit-axis invariant** (solver applicability):
Every NFFT solver of rank >= 1 declines a problem carrying an axis with
`N_t == 1` (`problem_nfft_has_unit_axis`). Unit axes are elided at problem
construction, so one reaching a solver means compression was bypassed — a bug
to surface, not a case to serve, even where the algorithm could handle it. The
rank-0 base case is unaffected: it has no axes, and it is what serves a problem
whose axes were all elided away. Only `mkproblem_nfft`'s borrowed path
(`copy_x = 0`) keeps full rank and can build such a problem.
_Avoid_: unit-axis support, degenerate-axis handling.

**configuration signature** (wisdom):
The MD5 over `sizeof(R)` followed by, in registration order, every registered
solver's `reg_id` and registrar name (`kernel/planner/planner.c`). It is
**global across all problem kinds** (iterated via `FORALL_SOLVERS`, not
per-kind), and is written on wisdom export and checked on import: a roster
mismatch (different precision, or a different solver set) degrades import to a
clean rejection, never a wrong plan. Consequence: registering a new kind's
solvers changes the signature, so wisdom from a process with a different
roster fails import — safe by design (a cache, never authoritative).
_Avoid_: per-kind signature, solver fingerprint (it is the whole-roster digest).

**tensor / mvdim** (planner):
The planner's geometry descriptor (`kernel/planner/tensor.c`): a `tensor`
is an array of rectangular `mvdim` factors `{n_in, is, n_out, os}`, each an
`n_out x n_in` strided operator, denoting their Kronecker product.
Generalises FFTW's square `iodim` (recovered as the state `n_in == n_out`);
the adjoint is an involution (swap sizes and strides per factor);
canonicalisation is behavior-matched to FFTW's `tensor_compress`(`_contiguous`).
Describes the structured stages of a transform — the scattered-node side of
NFFT is not a Kronecker factor and stays a flat count in the problem type.
_Avoid_: iodim (that is the square state, not a type here).

**size class**:
The structural identity a wisdom entry is keyed on: problem kind, geometry
(tensors, strides), window cutoff `m`, window type, power-of-two-bucketed node
count `M`, algorithm-relevant flags, precision, and thread count — but never
node coordinates. Wisdom means "for this size class on this machine, solver X
won."
_Avoid_: problem size (underspecified), problem hash (that is the key, not the class).

**estimate mode / measured mode**:
The two planning modes. *Estimate mode* (`NFFT_ESTIMATE` set) selects by the
fixed analytic cost model; nodes may be supplied later. *Measured mode* (the
default, flag absent — FFTW's convention) races the applicable candidates on
the user's actual nodes at plan time and blesses the winner. In the impatience
lattice, `ESTIMATE` rides in the upper bound only: "this search was allowed to
settle for estimate-grade evidence" — so measured wisdom answers estimate
queries, never the reverse.
_Avoid_: patient/exhaustive (not implemented), optimize (there is no separate refinement call).

**value-blind race**:
How measured mode times candidates: the real nodes drive the true memory-access
pattern (window support indices, node sorting), while the genuine value vectors
(`f_hat`, `f`) are irrelevant to the timing. Unlike FFTW (which zeroes its
operands), the NFFT race does **not** touch `f_hat`/`f` — it times `apply()` on
the caller's real arrays as-is. This is safe because the trip counts and access
pattern derive entirely from the nodes `x`, not from the values: forward reads
`f_hat` and writes `f`, adjoint the reverse, all node-driven. Timing is
pattern-faithful and value-independent; numerical output during the race is
meaningless by design. The ψ table is **precomputed, not zeroed**: `psi` and
`index_x` are functions of the nodes, not values, and zeroing them would destroy
the access pattern the race exists to measure. Precompute cost itself stays out of the
measurement (plans are executed many times), so what is timed is the
steady-state trafo. The NFFT analog of FFTW's zero-twiddle measurement —
which works here only because the access pattern derives from the nodes,
not the table values. The timer uses a **cycle-counter
fine clock** — raw ticks in arbitrary units via `include/cycle.h`'s
`getticks`/`elapsed`, timed the way FFTW times — with a **separate wall
clock** (`clock_gettime` `CLOCK_REALTIME`) bounding only the per-candidate
budget; where no tick counter exists it falls back to wall-clock seconds
(FFTW's slow-timer path), and where no clock exists it returns −1.0
(estimate-grade fallback).
_Avoid_: synthetic-node measurement (the nodes are real), bogus precompute
(nothing is precomputed for the race).

**precompute (awakening)**:
The deferrable lifecycle step (`nfft_precompute`) that materialises the real
node-dependent window values by awakening the bundle's direction plans; fast
plans build ψ exactly once per bundle awakening, direct plans awaken vacuously.
A measured-mode plan leaves its plans asleep — executing before precompute is
a checked error for every bundle (uniform rule).
_Avoid_: init-time precompute (that is `c_phi_inv`, node-independent, built at
construction).

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
