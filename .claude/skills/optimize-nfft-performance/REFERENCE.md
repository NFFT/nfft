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
deliberately breaking the target and watching what reacts. Then you optimize with a
known, trustworthy net (tests) and a known yardstick (benchmarks).

The shape of the whole task:

```
Step 0    open the task directory ── docs/perfeng/NNNN-<slug>/ + tracker README
Preflight pick the measurement track ── CodSpeed available? → simulation, else walltime
Phase A   full baseline  ── capture EVERY test + benchmark result (the exit reference)
Phase B   pin correctness net ── which tests guard the target?   [HARD GATE]
Phase C   pin performance metric ── which benchmark measures it?  [HARD GATE]
Phase D   inner loop ── optimize against the scoped net + metric
Phase E   exit gate ── re-run the FULL Phase-A baseline; no failure, no regression
```

**Front-to-back tracking.** The whole run lives in one directory per optimization
under [`docs/perfeng/`](../../../docs/perfeng/) — a tracker `README.md` plus one
deliverable per phase, the way decisions accumulate as ADRs under `docs/adr/`. Each
phase produces a concrete, canonical-format deliverable, and *deliverable = exit
gate*: a phase isn't done until its file is written and the tracker row flipped. The
layout, tracker template, and snapshot formats are in
[`details/deliverables.md`](details/deliverables.md) — read it before Phase A.

Phases B and C carry **hard gates**: if breaking the target fails *no* test, or if
*no* concrete benchmark metric can be obtained for it, the task **cannot proceed** —
that is a coverage gap to fix or escalate, not something to work around. The narrow
net/metric from B and C drive the fast inner loop (Phase D); the full Phase-A
baseline is the slow, authoritative check at the end (Phase E) that nothing outside
that narrow scope was broken or slowed down.

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
| **Preflight** — pick the measurement track | [`details/preflight.md`](details/preflight.md) |
| **Phase A** — build tree + full baseline | [`details/phase-a-baseline.md`](details/phase-a-baseline.md) |
| **Phase B** — pin the correctness net *[HARD GATE]* | [`details/phase-b-correctness-net.md`](details/phase-b-correctness-net.md) |
| **Phase C** — pin the performance metric *[HARD GATE]* | [`details/phase-c-performance-metric.md`](details/phase-c-performance-metric.md) |
| **Phase D** — inner loop | [`details/phase-d-inner-loop.md`](details/phase-d-inner-loop.md) |
| **Phase E** — exit gate | [`details/phase-e-exit-gate.md`](details/phase-e-exit-gate.md) |

Cross-cutting references, consult as needed:

| Topic | Doc |
|-------|-----|
| task directory layout, tracker, canonical deliverable formats | [`details/deliverables.md`](details/deliverables.md) |
| fill-in skeletons to copy per phase (tracker + one per deliverable) | [`templates/`](templates/) |
| walltime vs simulation; **working without CodSpeed** (the noise rule) | [`details/measurement-modes.md`](details/measurement-modes.md) |
| pitfalls (crash faults, OpenMP-only changes, narrow benchmark coverage, benign errors) | [`details/caveats.md`](details/caveats.md) |
| what is agent-operable, and what needs a CodSpeed account | [`details/tooling-status.md`](details/tooling-status.md) |

## Operational reference

| Step | Command |
|------|---------|
| Configure (local loop) | `cmake -S . -B build-cmake -DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON` |
| Build (lib + tests + benchmarks) | `cmake --build build-cmake -j` |
| Run tests | `ctest --test-dir build-cmake` |
| Granular test run | `build-cmake/tests/checkall` / `…/checkall_threads` (exit code + stdout `-> OK/FAIL`) |
| Failing cases only | `build-cmake/tests/checkall 2>&1 \| grep -E '\-> (FAIL\|ERROR)'` |
| Measure (walltime, local) | `CODSPEED_PROFILE_FOLDER=/tmp/b build-cmake/benchmarks/bench_nfft_direct --benchmark_filter='…'` → `/tmp/b/results/*.json` (`median_ns`) |
| Rough simulation (no account) | reconfigure `-DNFFT_BENCHMARK_MODE=simulation`; `valgrind --tool=callgrind … --benchmark_filter='<one case>'` → process-total `I refs` (incl. startup) |
| Clean per-case simulation | `codspeed run` — **needs a CodSpeed token** (`codspeed auth login`) |
| CI-history simulation | CodSpeed MCP (`https://mcp.codspeed.io/mcp`) — needs account + onboarded repo |
