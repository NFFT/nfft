# Performance optimization loop (test-pinned, benchmark-measured)

How a coding agent optimizes a scoped code region (a function, a hot loop) in this
repository **without human intervention**, keeping two independent feedback signals
in hand the whole time:

- **Correctness** — the CUnit accuracy suites (see
  [`test-methodology.md`](test-methodology.md)).
- **Performance** — the CodSpeed / google_benchmark binaries built by the **CMake**
  tree, measured locally in **walltime** mode (see
  [Measurement modes](#measurement-modes-walltime-locally-simulation-in-ci)).

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

**One self-contained CMake tree drives the whole loop, in `walltime` mode** — the only
mode that emits benchmark results locally and offline (a per-case stats JSON; no
runner, token, or upload). `FetchContent` fetches/builds codspeed-cpp (submodules and
all) and the tree produces the library, the CUnit tests, **and** the benchmark
binaries together. (Don't use the Autotools `--with-codspeed` path — it is legacy.)

```bash
cmake -S . -B build-cmake -DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON \
  -DCMAKE_C_FLAGS="-O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math"
cmake --build build-cmake -j   # libnfft3, checkall(_threads), bench_nfft_direct(_omp)
```

Both signals come from this tree:

- **correctness** — `ctest --test-dir build-cmake`, or run `build-cmake/tests/checkall`
  directly for the granular `-> OK/FAIL` stdout (same CUnit suite as `make check`).
- **performance** — `build-cmake/benchmarks/bench_nfft_direct[_omp]`, run directly
  (walltime). See [Measurement modes](#measurement-modes-walltime-locally-simulation-in-ci)
  for when to *also* consult the CI **simulation** metric.

For network-free rebuilds, reuse a checkout with
`-DFETCHCONTENT_SOURCE_DIR_CODSPEED=<path>`. (Autotools `make check` is a valid,
CI-canonical correctness path — see [`test-methodology.md`](test-methodology.md) — but
it does not build the benchmarks, so the loop stays in the CMake tree.)

Now record the complete state of **both** signals on the unmodified tree — the
contract the finished work is judged against in
[Phase D](#phase-d--exit-gate-full-baseline-re-check). It must be the *full* suite, not
the scoped subset used later.

```bash
ctest --test-dir build-cmake 2>&1 | tee /tmp/baseline-tests.log              # full suite
CODSPEED_PROFILE_FOLDER=/tmp/bench-base build-cmake/benchmarks/bench_nfft_direct \
    2>/tmp/baseline-bench.log                          # ALL cases → /tmp/bench-base/results/*.json
```

If the baseline is not fully green, **stop** — optimization starts only from a clean
tree. Keep these artifacts for the whole task.

## Measurement modes (walltime locally, simulation in CI)

Mode is **baked in at build time** via `-DNFFT_BENCHMARK_MODE=` (`off`, the default,
builds no benchmarks). The two measuring modes play different roles:

| | **walltime** — local workhorse | **simulation** — CI gate |
|---|---|---|
| Build | `-DNFFT_BENCHMARK_MODE=walltime` | `-DNFFT_BENCHMARK_MODE=simulation` |
| Run | binary directly, offline | under `valgrind --tool=callgrind` (or in CI) |
| Metric | wall-clock ns: min/mean/median/stdev/IQR | deterministic instruction count |
| Local result | JSON at `$CODSPEED_PROFILE_FOLDER/results/<pid>.json` | callgrind `I refs` on stderr |
| Authoritative source | the local JSON | **CodSpeed** (uploaded from CI) |
| Noise | real timing noise — compare medians | none (deterministic) |

**Use walltime for the local loop.** It is the only mode that writes results locally
and offline, so Phases 0/B/C/D below read its JSON. Switching mode is just a
reconfigure (`-DNFFT_BENCHMARK_MODE=…`); you do **not** need two build trees unless you
want both modes at once. Because walltime is noisy in a container, compare **medians**
over the (auto-tuned) rounds, not single runs.

**Simulation is the deterministic enhancement** (and what CI gates on). It is
*optional* — see [Working without CodSpeed](#working-without-codspeed). Three ways to
get simulation numbers, in increasing fidelity, none of the first two needing a
CodSpeed account:

1. **Raw `valgrind --tool=callgrind` on a `simulation` build** — quick, but the
   `I refs` it prints is a **process total**, not per-benchmark. Usable only by
   filtering to a *single* case (and it still includes process startup). Good for a
   fast before/after on one case.
2. **`codspeed run --skip-upload`** (the CLI, in the dev container) — runs offline, no
   token, and post-processes callgrind into clean **per-benchmark** instruction
   counts. This is the way to get per-case simulation locally.
3. **CodSpeed MCP / cloud** — the base branch's CI numbers, for cross-commit parity.
   Needs a CodSpeed account and the repo onboarded; see
   [Tooling status](#tooling-status-agent-operability). Use it in
   [Phase D](#phase-d--exit-gate-full-baseline-re-check) to confirm against CI.

Instruction count and wall-clock can diverge a lot (e.g. a change measured here was
2.3× fewer instructions but only 1.3× faster wall-clock — the loop became
memory/arithmetic-bound once the transcendentals went). The wall-clock figure is the
user-facing speedup; instruction count is the low-noise proxy.

### Working without CodSpeed

A coding agent may have **no CodSpeed access at all** (no account, repo not onboarded,
no MCP). That is a supported — if non-ideal — configuration: **the entire loop below
runs on walltime alone.** What you lose is the deterministic metric and the CI-history
comparison; what you keep is a complete, offline measure → change → re-measure cycle.
The cost is **timing noise**, which you must manage explicitly:

- Compare **medians**, never single runs.
- Treat a case as regressed only when its median rises **beyond noise** — past, say,
  `max(3·stdev, 2% of the median)`. Identical, untouched code routinely swings several
  percent here (worst on the few-iteration 2d/3d cases), so a hard "no case may be
  slower at all" rule produces **false regressions**.
- **Re-run** any case that trips the threshold before believing it; noise rarely
  survives a second run, a real regression does.

If `codspeed run --skip-upload` is available (it is in the dev container), prefer it
for the affected cases — deterministic instruction counts sidestep the noise entirely,
no account required.

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

**B1. Measure a baseline for the target's case(s)** with the walltime binary from
Phase 0 (run directly; it writes a per-case JSON):

```bash
CODSPEED_PROFILE_FOLDER=/tmp/b build-cmake/benchmarks/bench_nfft_direct \
    --benchmark_filter='nfft_forward_direct_1d.*'
# /tmp/b/results/<pid>.json → per case: stats.median_ns (use the median, it is robust)
```

**B2. Inject a slowdown** in the target — e.g. wrap its body in a `for` loop that
repeats the work N times. Keep results correct (so you isolate cost, not correctness).

**B3. Re-run and compare medians against the baseline.** Whichever benchmark cases'
`median_ns` rise clearly (well beyond `stdev_ns`) are the ones that exercise the
target — your performance metric. A case whose median is unchanged does **not** touch
the target.

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
2. Re-run the Phase-B benchmark case(s) (walltime) — `median_ns` should drop, and
   **must not** rise beyond noise, versus the saved baseline.

Iterate until the metric is satisfactory. The scoped checks here are *necessary but
not sufficient*: they are fast feedback, but they only see the narrow slice. The
authoritative verdict is Phase D.

## Phase D — exit gate (full baseline re-check)

A scoped optimization can have effects outside its scope — a shared helper, a header
change, a compiler-visibility shift, an aliasing assumption. So the task is **not
done** until the *entire* Phase-0 baseline is re-run and compared:

```bash
cmake --build build-cmake -j
ctest --test-dir build-cmake                            # FULL suite — must be all-pass
CODSPEED_PROFILE_FOLDER=/tmp/bench-final build-cmake/benchmarks/bench_nfft_direct \
    2>/tmp/final-bench.log                              # ALL cases (walltime) vs Phase-0
```

The full walltime pass is **slow** (the 2d/3d cases run for seconds each, twice — here
and in Phase 0). That cost is the price of catching out-of-scope regressions; if you
must trim it, run the affected families at full fidelity and the rest as a coarse
sanity pass, and **say so** in the report rather than silently skipping cases.

**Exit condition (all must hold):**

1. The full test suite passes exactly as in Phase 0 — no new failures, in either the
   single-threaded or the OpenMP (`checkall_threads`) library.
2. **No benchmark regresses beyond noise** versus the Phase-0 baseline — *every* case,
   not just the target's. Apply the noise rule from
   [Working without CodSpeed](#working-without-codspeed): a case counts as regressed
   only if its `median_ns` rises past `max(3·stdev, 2% of the median)`, **and** the
   rise survives a re-run. Do **not** fail the gate on raw walltime jitter — untouched
   code routinely swings a few percent. The target's own metric should improve (or be
   equal). *(With CodSpeed/simulation available, judge this on the deterministic
   instruction count instead — no noise rule needed.)*
3. **CI-parity (simulation), if available:** the deterministic instruction count must
   not regress either. Best: read the base branch's numbers from CodSpeed via the MCP
   server. Locally without CodSpeed: there is no simulation baseline from Phase 0 (it
   captured walltime), so reconstruct one — `git stash` the change, build a
   `simulation` tree, measure the affected cases, `git stash pop`, rebuild, re-measure,
   compare. If no simulation is available at all, this step is **skipped** — note in
   the report that the verdict rests on walltime only and CI is the final authority.
4. `git diff` contains only the intended optimization.

If any check fails, the optimization is **not complete**. The agent may loop back to
Phase C with further changes and re-evaluate this gate, or — if it cannot satisfy all
of them — **give up and revert**, reporting why. A faster target bought with a
regression or a broken test elsewhere is not a success.

---

## Operational reference

| Step | Command |
|------|---------|
| Configure (local loop) | `cmake -S . -B build-cmake -DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON` |
| Build (lib + tests + benchmarks) | `cmake --build build-cmake -j` |
| Run tests | `ctest --test-dir build-cmake` |
| Granular test run | `build-cmake/tests/checkall` / `…/checkall_threads` (exit code + stdout `-> OK/FAIL`) |
| Failing cases only | `build-cmake/tests/checkall 2>&1 \| grep -E '\-> (FAIL\|ERROR)'` |
| Measure (walltime, local) | `CODSPEED_PROFILE_FOLDER=/tmp/b build-cmake/benchmarks/bench_nfft_direct --benchmark_filter='…'` → `/tmp/b/results/*.json` (`median_ns`) |
| Per-case simulation (local, no account) | `codspeed run --skip-upload -- build-cmake/benchmarks/bench_nfft_direct --benchmark_filter='…'` |
| Quick simulation (single case) | reconfigure `-DNFFT_BENCHMARK_MODE=simulation`; `valgrind --tool=callgrind … --benchmark_filter='<one case>'` → process-total `I refs` |
| CI-history simulation | CodSpeed MCP (`https://mcp.codspeed.io/mcp`) — needs account + onboarded repo |

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
- **`[ERROR] instrument-hooks: failed to write environment.json` is benign.** A
  walltime run prints it when the **current working directory** (where the binary
  tries to drop `environment.json`) is not writable; the per-case results JSON is
  still written to `$CODSPEED_PROFILE_FOLDER/results/` and the exit code is 0. Run
  from a writable dir to silence it, or ignore it.
- **Walltime regressions on untouched code are noise, not bugs.** A change confined to
  one function (e.g. `trafo_direct`'s 1d branch) will still show neighbouring,
  *byte-identical* cases (`adjoint_direct`, other dims) swinging ±several percent run
  to run — worst on the few-iteration 2d/3d cases. This is why Phase D uses a noise
  threshold and a re-run, not a bare "nothing got slower" check.

## Tooling status (agent-operability)

Verified end to end in the dev container, using the CMake tree throughout. The full
five-phase loop was run on a real target (`X(trafo_direct)`, optimised by phase
recurrence): the fault enumerated a 149-case net, the slowdown isolated the
`forward_direct_1d` metric, and the change landed at ~1.3× wall-clock / ~2.3×
instructions with the net green — the gaps that surfaced are folded into the notes
above.

- **Phase A is fully agent-operable.** `cmake --build` + `ctest` / direct
  `build-cmake/tests/checkall`, the `-> FAIL` stdout signal, exit codes, and the CUnit
  XML all work with no human step. Fault-inject → observe → revert was exercised end to
  end (67 `-> FAIL` lines on the seeded fault; 0 after revert).
- **Phase B (local, walltime) is agent-operable.** `cmake -B build-cmake
  -DNFFT_BENCHMARK_MODE=walltime` built and ran offline end to end: `FetchContent`
  fetched codspeed-cpp `v2.3.0` (submodules auto-recursed), and the binary wrote a
  per-case stats JSON (`median_ns`, …) with no runner, token, or upload. The
  `simulation` build also measures under `valgrind --tool=callgrind`. `valgrind` and
  the `codspeed` CLI ship in the dev container (`.devcontainer/Dockerfile`).
- **Simulation locally needs no CodSpeed account.** Clean per-case instruction counts
  come from `codspeed run --skip-upload` (CLI, offline); a raw `valgrind` run gives a
  single-case process total. The loop is fully operable on **walltime alone** with no
  CodSpeed dependency at all (see
  [Working without CodSpeed](#working-without-codspeed)) — that is the baseline,
  accepting timing noise.
- **Optional — CodSpeed MCP for CI-history parity.** To compare against the base
  branch's CI numbers in Phase D, register the **CodSpeed MCP server**
  (`npx add-mcp https://mcp.codspeed.io/mcp --name CodSpeed`, or the Claude Code plugin
  `CodSpeedHQ/codspeed`). It exposes tools to list/compare runs and read flamegraphs.
  **Prerequisites the user must provide:** a CodSpeed account, the repo onboarded to
  CodSpeed, and CI uploading results — so this is an *enhancement*, not a requirement,
  and is **not wired by default**. Without it the Phase-D simulation check is the local
  `codspeed run --skip-upload` cross-check (or is skipped, walltime-only).
- **The Autotools benchmark path is legacy** (`./configure --enable-benchmarks
  --with-codspeed=<path>` + `make bench`): it needs a hand-built codspeed-cpp and
  prints no measurements. AGENTS.md §4 and the loop above use the CMake build instead.
- **No benchmark covers the fast path or the example target** — the substantive gap for
  *real* work. See [caveats](#caveats): only `*_direct` transforms are benchmarked, and
  `init`-only code like `intprod` sits outside every measured region. Closing perf work
  on `trafo`/`adjoint` requires *adding* benchmarks first.
