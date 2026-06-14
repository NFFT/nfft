# Performance optimization loop (test-pinned, benchmark-measured)

How a coding agent optimizes a scoped code region (a function, a hot loop) in this
repository **without human intervention**, keeping two independent feedback signals
in hand the whole time:

- **Correctness** — the CUnit accuracy suites (see
  [`test-methodology.md`](test-methodology.md)).
- **Performance** — the CodSpeed / google_benchmark binaries, built by the **CMake**
  tree (not Autotools — see [Phase 0](#phase-0--build-tree--full-baseline-the-exit-reference)).

The core idea: before changing anything for speed, *prove* which tests catch a
regression in the target and which benchmarks measure its speed. Do that by
deliberately breaking the target and watching what reacts. Then you optimize with a
known, trustworthy net (tests) and a known yardstick (benchmarks).

The shape of the whole task:

```
Phase 0  full baseline  ── capture EVERY test + benchmark result (the exit reference)
Phase A  pin correctness net ── which tests guard the target?   [HARD GATE]
Phase B  pin performance metric ── which benchmark measures it?  [HARD GATE]
Phase C  inner loop ── optimize against the scoped net + metric
Phase D  exit gate ── re-run the FULL Phase-0 baseline; no failure, no regression
```

Phases A and B carry **hard gates**: if breaking the target fails *no* test, or if
*no* concrete benchmark metric can be obtained for it, the task **cannot proceed** —
that is a coverage gap to fix or escalate, not something to work around. The narrow
net/metric from A and B drive the fast inner loop (Phase C); the full Phase-0
baseline is the slow, authoritative check at the end (Phase D) that nothing outside
that narrow scope was broken or slowed down.

The worked example below uses `intprod()` (`kernel/nfft/nfft.c:47`) as the
hypothetical target. Read the [caveats](#caveats) — `intprod` turns out to be a
*poor* target (it trips both hard gates), and the loop is what reveals that.

---

## Phase 0 — build tree + full baseline (the exit reference)

**Use one CMake tree for the whole loop.** It is self-contained: `FetchContent` pulls
and builds codspeed-cpp (submodules and all) and produces the library, the CUnit
tests, **and** the benchmark binaries together. Do *not* use the Autotools benchmark
path (`./configure --enable-benchmarks --with-codspeed=<path>`) — it is a pre-CodSpeed
leftover that needs a hand-built codspeed-cpp *and* a separately built library, and is
the source of the [legacy-path gaps](#tooling-status-agent-operability).

```bash
cmake -S . -B build-cmake \
  -DNFFT_ENABLE_BENCHMARKS=ON -DNFFT_ENABLE_OPENMP=ON -DCODSPEED_MODE=simulation \
  -DCMAKE_C_FLAGS="-O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math"
cmake --build build-cmake -j   # libnfft3, checkall(_threads), bench_nfft_direct(_omp)
```

This tree gives both signals:

- **correctness** — `ctest --test-dir build-cmake`, or run
  `build-cmake/tests/checkall` directly for the granular `-> OK/FAIL` stdout (same
  CUnit suite as `make check`).
- **performance** — `build-cmake/benchmarks/bench_nfft_direct[_omp]` under callgrind.

For network-free rebuilds, reuse a checkout with
`-DFETCHCONTENT_SOURCE_DIR_CODSPEED=<path>`. (Autotools `make check` is still a valid,
CI-canonical correctness path — see [`test-methodology.md`](test-methodology.md) — but
it does not build the benchmarks, so the loop below stays in the CMake tree.)

Now record the complete state of **both** signals on the unmodified tree — the
contract the finished work is judged against in
[Phase D](#phase-d--exit-gate-full-baseline-re-check). It
must be the *full* suite, not the scoped subset used later.

```bash
ctest --test-dir build-cmake 2>&1 | tee /tmp/baseline-tests.log   # full suite
valgrind --tool=callgrind --callgrind-out-file=/tmp/baseline-%p.cg \
    build-cmake/benchmarks/bench_nfft_direct 2>/tmp/baseline-bench.log   # ALL cases
```

If the baseline is not fully green, **stop** — optimization starts only from a clean
tree. Keep these artifacts for the whole task.

## Phase A — pin the correctness net

**A1. Identify the target.** A specific function / region, e.g. `intprod()`.

**A2. Inject a fault.** Make the *smallest* edit that changes the target's behaviour
— flip an operator, drop a term. The goal is to make dependent tests fail, not to be
realistic. Example: in `intprod`, `p *= vec[t] - a;` → `p += vec[t] - a;`.

**A3. Rebuild and see what flips to FAIL.** This set is your correctness net.

```bash
cmake --build build-cmake -j >/dev/null            # rebuild the library
build-cmake/tests/checkall >/tmp/buggy.log 2>&1; echo "exit=$?"
grep -E '\-> (FAIL|ERROR)' /tmp/buggy.log          # granular failing cases
```

Each line names the suite/case, the measured error and the bound, e.g.
`nfft_1d_50_50.txt … trafo_direct -> FAIL 5.7e+14 ( 1.07e-14)`. The affected
**suites** (`nfft` / `nfct` / `nfst`, never `util` for `intprod`) are the ones to
run in the inner loop. The detailed machine-readable report is
`tests/CUnitAutomated-Results.xml` (and `…_threads-Results.xml`).

> **HARD GATE — no failing test ⇒ no coverage ⇒ stop.** If the injected fault leaves
> *every* test green (exit 0, zero `-> FAIL`), the target is **not covered** by the
> suite. This is a blocking condition, not a green light: you cannot safely optimize
> code whose behaviour no test pins. Either add a test that fails on the fault first
> (then resume), or stop and report the coverage gap. Never optimize an uncovered
> region. (Try a *more destructive* fault before concluding there is no coverage — a
> too-subtle edit may stay within tolerance.)

**A4. Revert and re-confirm green.** Restore the exact original, rebuild, re-run;
`build-cmake/tests/checkall` must exit 0 with zero FAIL lines and `git diff` must be
empty.

You now know the precise tests that guard this region.

## Phase B — pin the performance metric

**B1. Measure a baseline for the target's case(s).** The benchmark binary built in
[Phase 0](#phase-0--build-tree--full-baseline-the-exit-reference) is
CodSpeed-**instrumented**: run directly it measures *nothing* ("running in an unknown
environment"). It only produces numbers under Valgrind/callgrind, which emits a
deterministic **instruction count** per benchmark.

```bash
valgrind --tool=callgrind --callgrind-out-file=/tmp/base.cg \
    build-cmake/benchmarks/bench_nfft_direct --benchmark_filter='nfft_forward_direct_1d.*'
```

Each benchmark case prints `Collected : <N>` / `I refs: <N>` (instruction
references) on stderr — that integer is the metric. Instruction count is
*deterministic* (no timing noise), which is exactly what an autonomous loop wants;
the trade-off is it measures instructions, not wall-clock. Record the per-case
counts as the baseline. (The `codspeed` CLI — `codspeed exec -- <binary>` — wraps
the same callgrind run and emits structured output if you prefer it.)

**B2. Inject a slowdown** in the target — e.g. wrap its body in a `for` loop that
repeats the work N times. Keep results correct (so you isolate cost, not
correctness). More instructions ⇒ higher count, so callgrind detects it reliably.

**B3. Re-run and diff against the baseline.** Whichever benchmark cases' instruction
counts rise are the ones that exercise the target — your performance metric. A case
whose count is unchanged does **not** touch the target.

> **HARD GATE — no obtainable metric ⇒ no progress.** The agent must be able to
> produce a *concrete, comparable* number for the target — in either **simulation**
> (callgrind instruction count) **or wall-clock** mode. Two ways this gate fails:
> (a) the measurement tooling cannot run at all (see [Tooling status](#tooling-status-agent-operability)),
> or (b) the tooling runs but the deliberate slowdown moves *no* benchmark — meaning
> no benchmark exercises the target, so a real improvement would be equally
> invisible. Either way, **stop**: add a benchmark that covers the target (see
> [caveats](#caveats)) or pick a different target. Optimization without a metric is
> unverifiable and must not proceed.

**B4. Revert** the slowdown; `git diff` must be empty.

## Phase C — inner loop (optimize against the scoped net + metric)

This is the fast loop, run against the *narrow* subset identified in A and B (not the
full suite — that is Phase D). After every change:

1. `cmake --build build-cmake -j && build-cmake/tests/checkall` (add
   `build-cmake/tests/checkall_threads` if the change touches OpenMP) — the Phase-A
   net must stay green.
2. Re-run the Phase-B benchmark case(s) under callgrind — the per-case instruction
   count should drop, and **must not** rise, versus the saved baseline.

Iterate until the metric is satisfactory. The scoped checks here are *necessary but
not sufficient*: they are fast feedback, but they only see the narrow slice. The
authoritative verdict is Phase D.

## Phase D — exit gate (full baseline re-check)

A scoped optimization can have effects outside its scope — a shared helper, a header
change, a compiler-visibility shift, an aliasing assumption. So the task is **not
done** until the *entire* Phase-0 baseline is re-run and compared:

```bash
cmake --build build-cmake -j
ctest --test-dir build-cmake                 # FULL suite — must be all-pass
valgrind --tool=callgrind --callgrind-out-file=/tmp/final-%p.cg \
    build-cmake/benchmarks/bench_nfft_direct 2>/tmp/final-bench.log   # ALL cases
```

**Exit condition (all must hold):**

1. The full test suite passes exactly as in Phase 0 — no new failures, in either the
   single-threaded or the OpenMP (`checkall_threads`) library.
2. **No** benchmark regresses versus the Phase-0 baseline — not just the target's,
   *every* benchmark. The target's metric should be improved (or at least equal).
3. `git diff` contains only the intended optimization.

If any check fails, the optimization is **not complete**. The agent may loop back to
Phase C with further changes and re-evaluate this gate, or — if it cannot satisfy all
three — **give up and revert**, reporting why. A faster target bought with a
regression or a broken test elsewhere is not a success.

---

## Operational reference

| Step | Command |
|------|---------|
| Configure build tree | `cmake -S . -B build-cmake -DNFFT_ENABLE_BENCHMARKS=ON -DNFFT_ENABLE_OPENMP=ON -DCODSPEED_MODE=simulation` |
| Build (lib + tests + benchmarks) | `cmake --build build-cmake -j` |
| Run tests | `ctest --test-dir build-cmake` |
| Granular test run | `build-cmake/tests/checkall` / `…/checkall_threads` (exit code + stdout `-> OK/FAIL`) |
| Failing cases only | `build-cmake/tests/checkall 2>&1 \| grep -E '\-> (FAIL\|ERROR)'` |
| Measure one benchmark | `valgrind --tool=callgrind build-cmake/benchmarks/bench_nfft_direct --benchmark_filter='…'` → read `I refs` / `Collected` |
| ~~Run benchmark directly~~ | inert — CodSpeed instrumentation makes no measurement outside Valgrind |

## Caveats

- **Crash-class faults lose granularity.** The Phase-A `intprod` fault corrupts the
  allocation sizes `N_total`/`n_total` (`kernel/nfft/nfft.c:5962`), so multi-dim
  cases under-allocate and the binary `SIGABRT`s (exit 134) partway through. The
  `-> FAIL` lines printed *before* the abort are still a reliable net, but the XML
  may be truncated. Prefer a fault that yields wrong-but-finite results over one that
  aborts when you want a complete failure list.
- **OpenMP-only changes show only in `checkall_threads`.** That binary links the
  `_omp` library; `checkall` does not. A change guarded by `#ifdef _OPENMP` (or
  affecting only parallel scheduling) will leave `checkall` green and must be judged
  by `checkall_threads`.
- **`intprod` is a deliberately bad example.** It runs only in `*_init*` (plan
  setup), and the current benchmarks call `init_*` **outside** the timed
  `for (auto _ : state)` loop (`benchmarks/bench_nfft_direct.cpp:65`) — they time
  only `trafo_direct` / `adjoint_direct`. So Phase B will show *no* benchmark moving
  for an `intprod` change: it is not a meaningful performance target under the
  current suite. Real targets live inside `trafo*` / `adjoint*`.
- **Benchmark coverage is narrow.** Only the **direct** (slow, O(N·M)) transforms
  are benchmarked, forward and adjoint, in 1d/2d/3d (`bench_nfft_direct.cpp`). The
  **fast** NFFT path (`trafo`, `adjoint`, the `precompute_one_psi` strategies) and
  every non-NFFT module have no benchmark yet. Optimizing those requires adding a
  benchmark first (Phase B with no coverage is meaningless).

## Tooling status (agent-operability)

Verified end to end in the dev container on branch `feature/coding-agents`, using the
CMake tree throughout:

- **Phase A is fully agent-operable.** `cmake --build` + `ctest` / direct
  `build-cmake/tests/checkall`, the `-> FAIL` stdout signal, exit codes, and the
  CUnit XML all work with no human step. Fault-inject → observe → revert was
  exercised end to end (67 `-> FAIL` lines on the seeded fault; 0 after revert).
- **Phase B is agent-operable via the CMake build.** `cmake -B build-cmake
  -DNFFT_ENABLE_BENCHMARKS=ON …` configured, built, and ran end to end:
  `FetchContent` cloned codspeed-cpp `v2.3.0` **with submodules recursed
  automatically** (`GIT_SUBMODULES_RECURSE TRUE` in `benchmarks/CMakeLists.txt`), and
  `build-cmake/benchmarks/bench_nfft_direct` under `valgrind --tool=callgrind`
  produced per-case instruction counts. `valgrind` and the `codspeed` CLI are
  installed by the dev container (`.devcontainer/Dockerfile`: `valgrind` via apt,
  `codspeed` via `https://codspeed.io/install.sh` → `~/.cargo/bin`). The bench binary
  is still *inert when run directly* (CodSpeed instrumentation: "unknown
  environment") — it only measures under callgrind; this is expected, not a gap.

- **The Autotools benchmark path is legacy and should not be used.**
  `./configure --enable-benchmarks --with-codspeed=<path>` predates the move to
  CodSpeed: it requires a **hand-built** codspeed-cpp and links a **separately built**
  library, and `make bench` just runs the binaries directly (so it prints no numbers).
  The hand-build is also where `AGENTS.md` §4 is wrong — its `git clone --depth 1`
  omits `--recurse-submodules`, so `core/instrument-hooks/` is empty and configure
  fails with `fatal error: core.h: No such file or directory`. **Recommended:** point
  `AGENTS.md` §4 at the CMake benchmark build (`-DNFFT_ENABLE_BENCHMARKS=ON`) and
  retire / clearly mark the Autotools `--with-codspeed` path as legacy.

- **No benchmark covers the fast path or the example target** — the one substantive
  gap left for *real* work. See [caveats](#caveats): only `*_direct` transforms are
  benchmarked, and `init`-only code like `intprod` sits outside every measured region.
  Closing perf work on `trafo`/`adjoint` requires *adding* benchmarks first.
