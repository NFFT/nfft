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
Preflight pick the measurement track ── CodSpeed available? → simulation, else walltime
Phase 0   full baseline  ── capture EVERY test + benchmark result (the exit reference)
Phase A   pin correctness net ── which tests guard the target?   [HARD GATE]
Phase B   pin performance metric ── which benchmark measures it?  [HARD GATE]
Phase C   inner loop ── optimize against the scoped net + metric
Phase D   exit gate ── re-run the FULL Phase-0 baseline; no failure, no regression
```

Phases A and B carry **hard gates**: if breaking the target fails *no* test, or if
*no* concrete benchmark metric can be obtained for it, the task **cannot proceed** —
that is a coverage gap to fix or escalate, not something to work around. The narrow
net/metric from A and B drive the fast inner loop (Phase C); the full Phase-0
baseline is the slow, authoritative check at the end (Phase D) that nothing outside
that narrow scope was broken or slowed down.

The worked example throughout uses `intprod()` (`kernel/nfft/nfft.c:47`) as the
hypothetical target. Read the [caveats](details/caveats.md) — `intprod` turns out to
be a *poor* target (it trips both hard gates), and the loop is what reveals that.

## Map — read these in order

| Step | Detail doc |
|------|-----------|
| **Preflight** — pick the measurement track | [`details/preflight.md`](details/preflight.md) |
| **Phase 0** — build tree + full baseline | [`details/phase-0-baseline.md`](details/phase-0-baseline.md) |
| **Phase A** — pin the correctness net *[HARD GATE]* | [`details/phase-a-correctness-net.md`](details/phase-a-correctness-net.md) |
| **Phase B** — pin the performance metric *[HARD GATE]* | [`details/phase-b-performance-metric.md`](details/phase-b-performance-metric.md) |
| **Phase C** — inner loop | [`details/phase-c-inner-loop.md`](details/phase-c-inner-loop.md) |
| **Phase D** — exit gate | [`details/phase-d-exit-gate.md`](details/phase-d-exit-gate.md) |

Cross-cutting references, consult as needed:

| Topic | Doc |
|-------|-----|
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
