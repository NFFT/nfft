# Performance optimization loop (test-pinned, benchmark-measured)

Overview and map. Each phase has its own detail doc under [`details/`](details/) —
read the phase doc when you reach that phase ([SKILL.md](SKILL.md) is the checklist).

How a coding agent optimizes a scoped code region (a function, a hot loop) in this
repository **without human intervention**, keeping two independent feedback signals
in hand the whole time:

- **Correctness** — the CUnit accuracy suites (see
  [`test-methodology.md`](../../../docs/agents/test-methodology.md)).
- **Performance** — the CodSpeed / google_benchmark binaries built by the **CMake**
  tree, measured locally in **walltime** mode (see
  [Measurement modes](details/measurement-modes.md)).

The core idea: before changing anything for speed, *prove* which tests catch a
regression in the target and which benchmarks measure its speed. Do that by
deliberately breaking the target and watching what reacts. Then — still before any speed
edit — *analyze* the target's rounding error under the standard model, so the optimization
is held to the right accuracy (the derived bound, not merely a passing test) and an
*avoidable* error source is rooted out rather than preserved. Then you optimize with a
known, trustworthy net (tests), a known yardstick (benchmarks), and a known error model.

The shape of the whole task:

```
Step 0    open the task directory ── docs/perfeng/NNNN-<slug>/ + tracker README
Phase A   full baseline  ── build walltime trees; capture EVERY test + benchmark result (the exit reference)
Phase B   pin correctness net ── which tests guard the target?   [HARD GATE]
Phase C   pin performance metric ── which benchmark measures it?  [HARD GATE]
Phase D   rounding-error analysis ── standard-model bound; set the accuracy objective
Phase E   inner loop ── optimize toward that objective against the scoped net + metric
Phase F   exit gate ── re-run the FULL Phase-A baseline; no failure, no regression
```

**Front-to-back tracking.** The whole run lives in one directory per optimization
under [`docs/perfeng/`](../../../docs/perfeng/) — a tracker `README.md` plus one
deliverable per phase. Each
phase produces concrete, canonical-format deliverables, and *deliverables = exit
gate*: a phase isn't done until its deliverables are written and the tracker row is flipped. The
layout, tracker template, and snapshot formats are in
[`details/deliverables.md`](details/deliverables.md) — read it before Phase A.

Phases B and C carry **hard gates**: if breaking the target fails *no* test, or if
*no* concrete benchmark metric can be obtained for it, the task **cannot proceed** —
that is a coverage gap to fix or escalate, not something to work around. The narrow
net/metric from B and C drive the fast inner loop (Phase E); the full Phase-A
baseline is the slow, authoritative check at the end (Phase F) that nothing outside
that narrow scope was broken or slowed down. Between the gates and the inner loop,
Phase D analyzes the target's rounding error and sets the *accuracy* objective the inner
loop optimizes toward (root out avoidable error first, or hold the derived bound).

Even that authoritative check is *bounded* — it covers specific sizes and grids. So every
run also carries a **risk assessment** ([`details/risk-assessment.md`](details/risk-assessment.md)):
a survey of the accuracy side effects the tests in hand can't see (e.g. a precision loss
that only shows at larger transforms), proven or not, surfaced in the human summary. When a
risk is material and cheaply testable, the net is **extended**
([`details/extending-tests.md`](details/extending-tests.md)).

The worked example throughout is `X(trafo_direct)()` — the direct, O(N·M) NDFT
(`kernel/nfft/nfft.c:145`): a real transformation, benchmarked forward and adjoint in
1d/2d/3d. Picking a *wrong* target is possible (code no test pins, or that no benchmark
measures), but you don't have to detect that up front — the Phase B and C hard gates
surface it for you: no failing test, or no benchmark moves, means **stop**. See
[caveats](details/caveats.md) for the common ways a target turns out unmeasurable.

## Map — read these in order

| Step | Detail doc |
|------|-----------|
| **Step 0 / deliverables** — task directory, tracker, canonical formats | [`details/deliverables.md`](details/deliverables.md) |
| **Phase A** — build tree + full baseline | [`details/phase-a-baseline.md`](details/phase-a-baseline.md) |
| **Phase B** — pin the correctness net *[HARD GATE]* | [`details/phase-b-correctness-net.md`](details/phase-b-correctness-net.md) |
| **Phase C** — pin the performance metric *[HARD GATE]* | [`details/phase-c-performance-metric.md`](details/phase-c-performance-metric.md) |
| **Phase D** — rounding-error analysis | [`details/phase-d-error-analysis.md`](details/phase-d-error-analysis.md) |
| **Phase E** — inner loop | [`details/phase-e-inner-loop.md`](details/phase-e-inner-loop.md) |
| **Phase F** — exit gate | [`details/phase-f-exit-gate.md`](details/phase-f-exit-gate.md) |

Cross-cutting references, consult as needed:

| Topic | Doc |
|-------|-----|
| task directory layout, tracker, canonical deliverable formats | [`details/deliverables.md`](details/deliverables.md) |
| **float · double · long double** matrix (three trees; A/D/E run all) | [`details/precision-matrix.md`](details/precision-matrix.md) |
| fill-in skeletons to copy per phase (tracker + one per deliverable) | [`templates/`](templates/) |
| **walltime** is the metric (local + CI); **the noise rule**; the optional callgrind tie-breaker | [`details/measurement-modes.md`](details/measurement-modes.md) |
| pitfalls (crash faults, OpenMP-only changes, narrow benchmark coverage, benign errors) | [`details/caveats.md`](details/caveats.md) |
| **risk assessment** — what a green run can still hide; the risk table | [`details/risk-assessment.md`](details/risk-assessment.md) |
| **extending the net** — online + refgen test cases to prove/disprove a risk | [`details/extending-tests.md`](details/extending-tests.md) |
| what is agent-operable (the whole walltime loop, no account) | [`details/tooling-status.md`](details/tooling-status.md) |

## Operational reference

| Step | Command |
|------|---------|
| Configure (local loop) | `cmake -S . -B build-cmake -DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON` |
| Build (lib + tests + benchmarks) | `cmake --build build-cmake -j` |
| Run tests | `ctest --test-dir build-cmake` |
| Granular test run | `build-cmake/tests/checkall` / `…/checkall_threads` (exit code + stdout `-> OK/FAIL`) |
| Failing cases only | `build-cmake/tests/checkall 2>&1 \| grep -E '\-> (FAIL\|ERROR)'` |
| Measure (walltime — the metric, local + CI) | `CODSPEED_PROFILE_FOLDER=/tmp/b build-cmake/benchmarks/bench_nfft_direct --benchmark_filter='…'` → `/tmp/b/results/*.json` (`median_ns`) |
| Tie-breaker only (layout vs real work, no account) | build `-DNFFT_BENCHMARK_MODE=simulation`; `valgrind --tool=callgrind … --benchmark_filter='<one case>'` → process-total `I refs` delta ([measurement-modes](details/measurement-modes.md)) |
