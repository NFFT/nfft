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
[deliverables](details/deliverables.md), [measurement-modes](details/measurement-modes.md),
[caveats](details/caveats.md) and [tooling-status](details/tooling-status.md) docs as needed.

**Every phase produces a deliverable.** The loop is tracked front-to-back in one
directory per optimization under [`docs/perfeng/`](../../../docs/perfeng/) (the way
decisions accumulate as ADRs under `docs/adr/`). *Deliverable = exit gate:* a phase
is not done until its deliverable file exists in the task directory and the tracker
row is updated. The format, layout, and canonical snapshot shapes are in
[deliverables.md](details/deliverables.md) — read it before Phase A. Don't hand-write
deliverables: each phase has a fill-in skeleton under [`templates/`](templates/);
copy the named one and fill its `<…>` placeholders.

## The loop (do these in order)

Create one TodoWrite item per phase. **Phases B and C are HARD GATES — if they fail,
stop and report, do not work around them** (record the gate failure as a `blocked`
deliverable — see [deliverables.md](details/deliverables.md)).

- [ ] **Step 0 — open the task directory.** Create `docs/perfeng/NNNN-<target-slug>/`,
  copy [`templates/tracker.md`](templates/tracker.md) into it as `README.md` and set
  the target, add an index row to `docs/perfeng/README.md`. Every subsequent phase
  copies its template here and flips its tracker row. → *Deliverable:* the task
  directory + initialized tracker.

- [ ] **[Preflight — pick the measurement track](details/preflight.md).** CodSpeed fully available (MCP present
  + `codspeed status` logged in + repo onboarded with a base-branch run)? → **simulation**
  track (deterministic instruction count, what CI gates on). Any no → **walltime** track
  (local, offline, noisy; compare medians). State which track you're on in the report.
  → *Deliverable:* `preflight.md` (chosen track + the evidence for it); tracker **Track**
  field and Preflight row set.

- [ ] **[Phase A — baseline](details/phase-a-baseline.md).** Build the self-contained CMake tree in the chosen mode and
  capture the FULL test + benchmark state (the exit reference for Phase E):
  ```bash
  cmake -S . -B build-cmake -DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON \
    -DCMAKE_C_FLAGS="-O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math"
  cmake --build build-cmake -j
  ctest --test-dir build-cmake 2>&1 | tee /tmp/baseline-tests.log
  CODSPEED_PROFILE_FOLDER=/tmp/bench-base build-cmake/benchmarks/bench_nfft_direct 2>/tmp/baseline-bench.log
  ```
  Baseline not fully green → **stop**; optimization starts only from a clean tree.
  → *Deliverable:* `phase-a-baseline.md` (build config + commit + baseline snapshot
  table) and the raw `artifacts/baseline-{tests.log,bench.json}` — this is the exit
  reference Phase E is judged against.

- [ ] **[Phase B — pin the correctness net](details/phase-b-correctness-net.md) [HARD GATE].** Inject the smallest fault into
  the target (flip an operator, drop a term), rebuild, run `build-cmake/tests/checkall`,
  and record which cases flip to `-> FAIL`. That set is your net. **No test fails ⇒ the
  region is uncovered ⇒ stop and exit** (try a more destructive fault first; then report 
  the gap — never optimize uncovered code). Revert, confirm green, `git diff` empty.
  → *Deliverable:* `phase-b-correctness-net.md` (the fault, the net table + suite + size,
  revert confirmed) and `artifacts/fault.diff`; on gate failure this is the `blocked`
  coverage-gap report instead.

- [ ] **[Phase C — pin the performance metric](details/phase-c-performance-metric.md) [HARD GATE].** Measure the target's
  benchmark case(s), inject a slowdown (wrap the body in an N-times loop, results still
  correct), re-measure. Cases whose `median_ns` rise clearly are your metric. **No
  benchmark moves (or tooling can't measure) ⇒ stop and exit**. Revert, `git diff` empty.
  → *Deliverable:* `phase-c-performance-metric.md` (target baseline snapshot, the metric
  case(s) with before/after, revert confirmed) and `artifacts/slowdown.diff`; on gate
  failure this is the `blocked` no-metric report instead.

- [ ] **[Phase D — inner loop](details/phase-d-inner-loop.md).** Optimize against the *narrow* B-net + C-metric only.
  After every change: rebuild, `build-cmake/tests/checkall` (add `checkall_threads` if
  OpenMP touched) stays green; the C benchmark `median_ns` drops and never rises beyond
  noise. Fast feedback, but not authoritative — Phase E is. → *Deliverable:*
  `phase-d-inner-loop.md` (iteration journal: per-change net result + metric before→after)
  and the current `artifacts/change.diff`.

- [ ] **[Phase E — exit gate](details/phase-e-exit-gate.md).** Re-run the ENTIRE Phase-A baseline (full `ctest` + ALL
  benchmark cases). Exit only when: full suite passes as in Phase A; no benchmark
  regresses beyond `max(3·stdev, 2% of median)` *and* the rise survives a re-run; the
  target's metric improved or equal; `git diff` is only the intended change. Any failure
  ⇒ loop back to D, or revert and report why. → *Deliverable:* `phase-e-exit-gate.md`
  (full comparison table over ALL cases + the four-point checklist + verdict) with raw
  `artifacts/final-{tests.log,bench.json}`; then close out the tracker (**Status** =
  complete | reverted, **Outcome** one-liner), update the `docs/perfeng/` index, and
  write the human report `summary.html` (the reviewer-facing walkthrough — produced on
  *every* exit, including a Phase B/C hard-gate block). A faster target bought with a
  regression or a broken test elsewhere is **not** a success.

## Key rules

- **Compare medians, never single runs** (walltime is often noisy). Untouched,
  byte-identical cases swing several percent run to run — that's noise, not a regression.
- **One CMake tree drives everything**; switching measurement mode is just a reconfigure
  with a different `-DNFFT_BENCHMARK_MODE=`. Don't use the legacy Autotools `--with-codspeed` path.
- **OpenMP-only changes show only in `checkall_threads`** (it links the `_omp` library).
- Respect the precision-agnostic conventions (`Y()`/`X()`/`FFTW()` mangling, `R`/`C`/`E`
  types) — see `CONVENTIONS.md`. Keep the float/double/long-double matrix building.
