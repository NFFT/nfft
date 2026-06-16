---
name: optimize-nfft-performance
description: Safely optimize a scoped C region (a function, a hot loop) in the NFFT3 library without human intervention, keeping correctness (CUnit) and performance (CodSpeed/walltime benchmarks) as independent signals. Use when asked to optimize, speed up, or reduce the cost of a specific function/loop/kernel in this repo, or when working through the performance-optimization loop.
---

# Optimize NFFT performance (test-pinned, benchmark-measured)

Optimize a scoped region with two independent signals in hand the whole time:
**correctness** (CUnit accuracy suites) and **performance** (CodSpeed/google_benchmark
binaries). Before changing anything for speed, *prove* which tests catch a regression
and which benchmark measures the speed — by deliberately breaking the target and
watching what reacts. Then optimize with a trustworthy net and a known yardstick.

**Full methodology:** [REFERENCE.md](REFERENCE.md) is the overview + map; each phase
below links to its own detail doc under [`details/`](details/). This SKILL is the
enforced checklist — open a phase's doc when you reach it, plus the cross-cutting
[measurement-modes](details/measurement-modes.md), [caveats](details/caveats.md)
and [tooling-status](details/tooling-status.md) docs as needed.

## The loop (do these in order)

Create one TodoWrite item per phase. **Phases A and B are HARD GATES — if they fail,
stop and report, do not work around them.**

- [ ] **[Preflight — pick the measurement track](details/preflight.md).** CodSpeed fully available (MCP present
  + `codspeed status` logged in + repo onboarded with a base-branch run)? → **simulation**
  track (deterministic instruction count, what CI gates on). Any no → **walltime** track
  (local, offline, noisy; compare medians). State which track you're on in the report.

- [ ] **[Phase 0 — baseline](details/phase-0-baseline.md).** Build the self-contained CMake tree in the chosen mode and
  capture the FULL test + benchmark state (the exit reference for Phase D):
  ```bash
  cmake -S . -B build-cmake -DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON \
    -DCMAKE_C_FLAGS="-O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math"
  cmake --build build-cmake -j
  ctest --test-dir build-cmake 2>&1 | tee /tmp/baseline-tests.log
  CODSPEED_PROFILE_FOLDER=/tmp/bench-base build-cmake/benchmarks/bench_nfft_direct 2>/tmp/baseline-bench.log
  ```
  Baseline not fully green → **stop**; optimization starts only from a clean tree.

- [ ] **[Phase A — pin the correctness net](details/phase-a-correctness-net.md) [HARD GATE].** Inject the smallest fault into
  the target (flip an operator, drop a term), rebuild, run `build-cmake/tests/checkall`,
  and record which cases flip to `-> FAIL`. That set is your net. **No test fails ⇒ the
  region is uncovered ⇒ stop and exit** (try a more destructive fault first; then report 
  the gap — never optimize uncovered code). Revert, confirm green, `git diff` empty.

- [ ] **[Phase B — pin the performance metric](details/phase-b-performance-metric.md) [HARD GATE].** Measure the target's
  benchmark case(s), inject a slowdown (wrap the body in an N-times loop, results still
  correct), re-measure. Cases whose `median_ns` rise clearly are your metric. **No
  benchmark moves (or tooling can't measure) ⇒ stop and exit**. Revert, `git diff` empty.

- [ ] **[Phase C — inner loop](details/phase-c-inner-loop.md).** Optimize against the *narrow* A-net + B-metric only.
  After every change: rebuild, `build-cmake/tests/checkall` (add `checkall_threads` if
  OpenMP touched) stays green; the B benchmark `median_ns` drops and never rises beyond
  noise. Fast feedback, but not authoritative — Phase D is.

- [ ] **[Phase D — exit gate](details/phase-d-exit-gate.md).** Re-run the ENTIRE Phase-0 baseline (full `ctest` + ALL
  benchmark cases). Exit only when: full suite passes as in Phase 0; no benchmark
  regresses beyond `max(3·stdev, 2% of median)` *and* the rise survives a re-run; the
  target's metric improved or equal; `git diff` is only the intended change. Any failure
  ⇒ loop back to C, or revert and report why. A faster target bought with a regression or
  a broken test elsewhere is **not** a success.

## Key rules

- **Compare medians, never single runs** (walltime is often noisy). Untouched,
  byte-identical cases swing several percent run to run — that's noise, not a regression.
- **One CMake tree drives everything**; switching measurement mode is just a reconfigure
  with a different `-DNFFT_BENCHMARK_MODE=`. Don't use the legacy Autotools `--with-codspeed` path.
- **OpenMP-only changes show only in `checkall_threads`** (it links the `_omp` library).
- Respect the precision-agnostic conventions (`Y()`/`X()`/`FFTW()` mangling, `R`/`C`/`E`
  types) — see `CONVENTIONS.md`. Keep the float/double/long-double matrix building.
