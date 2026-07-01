---
name: nfft-perf-eng
description: Safely optimize a scoped C region in the NFFT3 library for performance and/or accuracy without human intervention. Tests (CUnit) and benchmarks (CodSpeed walltime) ensure correctness and performance as independent signals. Deterministic steps (build, capture, compare, trend fit) are driven by helper scripts. Use when asked to optimize, speed up, or reduce the cost of a specific function/loop/kernel in this repo.
---

# Optimize NFFT performance (test-pinned, benchmark-measured)

Optimize a scoped region with two independent signals in hand the whole time:
**correctness** (CUnit accuracy suites) and **performance** (CodSpeed/google_benchmark
binaries). Before changing anything for speed or accuracy, *prove* which tests catch 
a regression and which benchmark measures the speed — by deliberately breaking the target and
watching what reacts. Then optimize with a trustworthy net and a known yardstick.

**Full methodology:** [REFERENCE.md](REFERENCE.md) is the overview + map; each phase
below links to its own detail doc under [`details/`](details/). This SKILL is the
enforced checklist — open a phase's doc when you reach it, plus the cross-cutting
[deliverables](details/deliverables.md), [precision-matrix](details/precision-matrix.md),
[measurement-modes](details/measurement-modes.md), [caveats](details/caveats.md),
[risk-assessment](details/risk-assessment.md), [extending-tests](details/extending-tests.md)
and [tooling-status](details/tooling-status.md) docs as needed.

**Every phase produces deliverables.** The loop is tracked front-to-back in one
directory per optimization under [`docs/perfeng/`](../../../docs/perfeng/). 
*Deliverables = exit gate:* a phase is not done until its deliverables exist in the task 
directory and the tracker row is updated. The format, layout, and canonical snapshot shapes are in
[deliverables.md](details/deliverables.md) — read it before Phase A. Don't hand-write
deliverables: each phase has a fill-in skeleton under [`templates/`](templates/);
copy the named one and fill its `<…>` placeholders.

## The loop (do these in order)

Create one TodoWrite item per phase. **Phases B and C are HARD GATES — if they fail,
stop and report, do not work around them** (record the gate failure as a `blocked`
deliverable — see [deliverables.md](details/deliverables.md)).

- [ ] **Step 0 — open the task directory.** Run
  [`scripts/perf-init.sh <target-slug>`](scripts/perf-init.sh): it creates
  `docs/perfeng/NNNN-<target-slug>/`, copies every template in (tracker → `README.md`, the
  phase docs, `error-analysis.html`), stamps the commit, and prints the index row to add to
  `docs/perfeng/README.md`. Set the **Target** line; add that printed index row. → *Deliverable:*
  the task directory + initialized tracker.

- [ ] **[Phase A — baseline](details/phase-a-baseline.md).** Build all three precision trees in
  **walltime** mode and capture the FULL test + benchmark state (the exit reference for Phase F) — **in
  all three precisions** (the sources are precision-agnostic; see [precision-matrix](details/precision-matrix.md)).
  Both steps are deterministic:
  ```bash
  SCR=.claude/skills/nfft-perf-eng/scripts
  $SCR/perf-build.sh walltime                                   # configure+build build-cmake{,-f,-l}
  $SCR/perf-capture.sh baseline docs/perfeng/0001-trafo-direct  # full ctest + all bench cases → artifacts/, collated per precision
  ```
  `perf-capture.sh` exits non-zero if any precision's baseline isn't fully green → **stop**;
  optimization starts only from a clean tree. Build the snapshot table with
  `uv run python $SCR/perf-bench.py snapshot artifacts/baseline-bench-<p>.json --prec <p>`.
  → *Deliverable:* `phase-a-baseline.md` (build config + commit + baseline snapshot table with
  a `prec` column) and raw per-precision `artifacts/baseline-{tests,bench}-{d,f,l}.*` — the
  exit reference Phase F is judged against.

- [ ] **[Phase B — pin the correctness net](details/phase-b-correctness-net.md) [HARD GATE].** Inject the smallest fault into
  the target (flip an operator, drop a term — agent's call), rebuild, run `build-cmake/tests/checkall`;
  `perf-net.py` turns its output into the net (the cases that flip to `-> FAIL`). **No test fails ⇒ the
  region is uncovered ⇒ stop and exit** (try a more destructive fault first; then report 
  the gap — never optimize uncovered code). Revert, confirm green, `git diff` empty.
  → *Deliverable:* `phase-b-correctness-net.md` (the fault, the net table + suite + size,
  revert confirmed) and `artifacts/fault.diff`; on gate failure this is the `blocked`
  coverage-gap report instead.

- [ ] **[Phase C — pin the performance metric](details/phase-c-performance-metric.md) [HARD GATE].** Measure the target's
  benchmark case(s), inject a slowdown (e.g. wrap the body in an N-times loop but beware of compiler optimizations, results still
  correct), re-measure. Cases whose `median_ns` rise clearly are your metric. **No
  benchmark moves ⇒ the region is uncovered/tooling can't measure ⇒ stop and exit**. Revert, `git diff` empty.
  → *Deliverable:* `phase-c-performance-metric.md` (target baseline snapshot, the metric
  case(s) with before/after, revert confirmed) and `artifacts/slowdown.diff`; on gate
  failure this is the `blocked` no-metric report instead.

- [ ] **[Phase D — rounding-error analysis](details/phase-d-error-analysis.md).** *Before any speed
  edit*, scrutinize the target's rounding error under the **standard model**
  (`fl(a∘b)=(a∘b)(1+δ), |δ|≤u`): derive the forward-error bound with any parameter
  (`N`, `M`, or other) dependence explicit, and reconcile it with the Phase-A baseline errors (and a
  [`scripts/perf-trend.py`](scripts/perf-trend.py) error-vs-`N` fit if growth is in question).
  Conclude with a verdict that **sets the inner loop's objective** — `clean` (error at the best
  achievable floor → optimize for speed without regressing the derived bound) or `improve-first`
  (an *avoidable* `N`/`M` dependence, cancellation, or unstable recurrence → **root out the error
  first, then tune speed**; often both improve). Seed the risk table with the size/precision
  dependence the analysis exposes. → *Deliverable:* `error-analysis.html` (self-contained,
  MathJax — the full derivation) **and** `phase-d-error-analysis.md` (verdict + objective).

- [ ] **[Phase E — inner loop](details/phase-e-inner-loop.md).** Optimize *toward the Phase-D
  objective* against the *narrow* B-net + C-metric only. After every change: rebuild **all three
  precision trees** and run the B-net (`checkall`, add `checkall_threads` if OpenMP touched) in
  each — it must stay green in float, double, *and* long double; re-measure the C-metric in each —
  `median_ns` drops and never rises beyond noise. The net is narrow, so 3× is cheap and catches a
  float-only break early. Fast feedback, but not authoritative — Phase F is. Keep the **risk
  table** current: for each change, note the accuracy/correctness side effects this narrow net
  can't see (size-dependent, precision-specific, input-range — see
  [risk-assessment](details/risk-assessment.md)); when a risk is material and cheaply testable,
  extend the net ([extending-tests](details/extending-tests.md)) to prove or disprove it.
  → *Deliverable:* `phase-e-inner-loop.md` (iteration journal: per-change net result + metric
  before→after, per precision) and the current `artifacts/change.diff`.

- [ ] **[Phase F — exit gate](details/phase-f-exit-gate.md).** Re-run the ENTIRE Phase-A baseline
  (full `ctest` + ALL benchmark cases) **in all three precisions** — `perf-capture.sh final …`
  then `uv run python …/perf-bench.py compare --taskdir …` (deterministic noise-rule verdict).
  Exit only when, *in float, double, and long double*: the full suite passes as in Phase A; no
  benchmark regresses beyond `max(3·stdev, 2% of median)` *and* the rise survives a re-run; the
  target's metric improved or equal; **the Phase-D accuracy objective is honoured** (avoidable
  source removed / derived bound not regressed); `git diff` is only the intended change; and the
  **risk assessment** is complete and honest — every material risk `proven`, `retired`,
  `accepted`, or `residual` ([risk-assessment](details/risk-assessment.md)). A regression in
  *any* precision fails the gate. A `proven` material accuracy drop is a **hard no** — fix or
  revert, never ship; an unsettled `residual` material drop may land but only as a `partial` run
  that surfaces it (with the settling test proposed) for the human reviewer. Any failure ⇒ loop
  back to E, or revert and report why. → *Deliverable:* `phase-f-exit-gate.md` (per-precision
  comparison table over ALL cases + the six-point checklist + verdict) with raw
  `artifacts/final-{tests,bench}-{d,f,l}.*`; then close out the tracker (**Status** = complete |
  reverted, **Outcome** one-liner), update the `docs/perfeng/` index, and write the human report
  `summary.html` (the reviewer-facing walkthrough — produced on *every* exit, including a Phase
  B/C hard-gate block; it presents every phase's result + numbers, **embeds the required charts**
  via `perf-summary.py charts`, and links every deliverable + artifact — `perf-summary.py check`
  verifies nothing is orphaned). A faster target bought with a regression or a broken test
  elsewhere is **not** a success.

## Key rules

- **Measure all three precisions** (float · double · long double) at the baseline, the inner
  loop, and the exit gate — one CMake tree per precision (`build-cmake{,-f,-l}`). A
  precision-agnostic target can regress in just one precision; double-only misses it. See
  [precision-matrix](details/precision-matrix.md).
- **A green suite is narrow, not a guarantee.** Passing tests prove correctness only for
  the sizes, grids, and distributions the cases cover. Every run produces a **risk
  assessment** of the accuracy side effects that fall outside that coverage — proven or
  unproven, all surfaced in `summary.html` — and extends the net to settle the material
  ones ([risk-assessment](details/risk-assessment.md), [extending-tests](details/extending-tests.md)).
- **Accuracy is analyzed before speed, not just tested.** Phase D derives the target's
  rounding-error bound under the standard model *before any speed edit* and sets the objective:
  an `improve-first` verdict makes *removing avoidable error* the inner loop's first job; a
  `clean` verdict makes the derived bound — not merely test pass/fail — the accuracy yardstick the
  exit gate enforces ([phase-d-error-analysis](details/phase-d-error-analysis.md)).
- **Deterministic steps are scripted, not hand-run.** [`scripts/`](scripts/) drives the build,
  capture, JSON collation, test-net parsing, noise-rule comparison, trend fit, chart generation,
  and the write-up completeness check so they are reproducible and token-light; the judgement
  steps (the fault, the slowdown, the analysis, the optimization) stay with you. See
  [scripts/README.md](scripts/README.md).
- **The write-up shows the whole run.** `summary.html` is the reviewer's artifact: every phase's
  result with numbers, the **required charts** (performance speedup; error-vs-`N` trend), and
  links to *every* deliverable and raw artifact — generated and verified by `perf-summary.py
  charts`/`check`. See [deliverables.md](details/deliverables.md#required-visualizations).
- **Walltime is the only metric — local and CI both use it** (CodSpeed Macro Runners run
  walltime), so a local result previews the CI gate. No simulation/instruction-count gate,
  no CodSpeed account, token, MCP, or cloud upload in the loop.
- **Compare medians, never single runs** (walltime is noisy). Untouched, byte-identical
  cases swing several percent run to run — that's noise, not a regression; re-run to confirm.
- **One CMake tree per precision drives that precision** (tests + benchmark), built with
  `-DNFFT_BENCHMARK_MODE=walltime`. The `simulation` build is used *only* for the optional
  callgrind tie-breaker on an ambiguous untouched control case ([measurement-modes](details/measurement-modes.md)).
  Don't use the legacy Autotools `--with-codspeed` path.
- **OpenMP-only changes show only in `checkall_threads`** (it links the `_omp` library).
- Respect the precision-agnostic conventions (`Y()`/`X()`/`FFTW()` mangling, `R`/`C`/`E`
  types) — see `CONVENTIONS.md`. Keep the float/double/long-double matrix building.
