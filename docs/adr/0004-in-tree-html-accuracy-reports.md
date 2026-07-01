# ADR-0004: In-tree HTML accuracy reports

## Status
Accepted (2026-06-30).

## Context
The `nfft` / `nfct` / `nfst` CUnit suites already compute, per case, a relative
error `err` and compare it to an analytic `bound`. That pass/fail gate stays in C
(`err < bound`). What it does not give us is *visibility into the value*: slow drift
in `err` is invisible until a case finally crosses its bound, and a PR has no way to
show that it improved or regressed accuracy. We want to track those figures over time
and diff them per PR — without standing up an external service, a project secret, or
gated CI environments for what an in-tree, GitHub Pages-served report can archive itself.

## Decision

- **Aggregate, don't series-per-case.** A run produces ~1456 cases per matrix cell;
  most sit near the round-off floor and would drown a per-case view in noise. Instead
  we aggregate into one **accuracy metric** per combination of *error-shaping
  parameters* — module, runtime (serial/OpenMP), oracle (file/online), speed
  (direct/fast), direction (forward/adjoint), dimension, init variant — collapsing the
  *bound-absorbed parameters* (`N`, `M`, whose effect the analytic bound already
  captures) via `max`. The metric name is
  `<module>/<runtime>/<oracle>/<speed>/<direction>/<dim>d/<init-slug>`.

- **Keep file and online oracles separate.** Merging them lets a small file case mask
  a regression in the larger-`N`, more realistic *online* cases; online uses a
  different oracle (the C direct transform), so it is an error-shaping difference, not
  merely a larger `N`. Serial and OpenMP are likewise distinct metrics: the parallel
  reduction order perturbs the low bits, so `runtime` is an error-shaping axis.

- **`accuracy-digits` is the primary measure.** Raw errors span ~14 orders of
  magnitude (1e-17 .. 1e-3), which no linear display renders legibly. The primary
  measure is `accuracy-digits` = `-log10(max err)` — the worst-case accurate digits,
  which reads cleanly (~3 .. 18) and where a regression *lowers* the value. The exact
  worst error is kept as the secondary `max-error`, and `-log10(bound)` as
  `bound-digits` so a cell can be colored by its **margin** (digits beyond the bound).

- **The C harness stays dumb.** It emits one raw, structured NDJSON record per case
  behind the `NFFT_BENCH_OUT` env var (a no-op when unset, so an ordinary `make check`
  is unaffected). All grouping and `max` aggregation live in a Python converter
  (`tests/accuracy/ndjson_to_bmf.py` → the aggregated accuracy JSON, an in-tree format
  we label "BMF"), so the policy is tunable without touching C or re-running builds.

- **Render self-contained HTML, serve via GitHub Pages.** `tests/accuracy/htmlreport.py`
  (Python stdlib only) emits one HTML document: CSS-only window tabs, one table per
  module, per-window testbeds as columns, margin-colored cells. No matplotlib, no
  numpy, no JavaScript. The develop dashboard is published to the `gh-pages` branch
  root; each PR's report under `pr/<n>/`. The PR comment carries a per-module
  improved/regressed rollup plus a link to the report.

## Consequences

- **Metric-name stability matters.** `diff.py` joins PR-vs-baseline BMFs by metric
  name, so changing the error-shaping/bound-absorbed split later breaks the diff
  (renamed metrics show as spurious added/removed) — which is why the split is recorded
  here.
- **Emission piggybacks on the existing test run.** Rather than a dedicated workflow
  that re-builds and re-runs the gcc matrix, emission is switched on inside
  `build-linux.yml`'s existing `make check` (via `NFFT_BENCH_OUT`); each cell publishes
  its BMF as an artifact with no secret. Tests run exactly once.
- **Never gates CI.** The report is informational; CUnit's `err < bound` remains the
  sole gate. The report steps are `continue-on-error` and the Check is `neutral`, so a
  publish or API failure can never red CI.
- **Pages dependency.** The full report depends on GitHub Pages serving the `gh-pages`
  branch; if Pages is disabled the comment link 404s (the per-module rollup in the
  comment is Pages-independent and still conveys whether anything changed).
- **Fork PRs.** A fork's build gets a read-only token, so the in-build report job is
  skipped; a trusted `workflow_run` companion (`accuracy-report-fork.yml`) instead
  posts the aggregate comment + Check over the fork's BMF artifacts (data only, no
  `gh-pages` write; the link points at the develop dashboard).
- **Track-only, phased.** Pure recording first (no thresholds — they misfire on a cold
  baseline); hard CI gating, if ever, only by explicit later decision.
- **Scope** is `nfft` / `nfct` / `nfst` (the three structurally identical harnesses);
  other modules are follow-up. Window and precision are carried by the **accuracy
  testbed** (`<os>_<compiler>_<window>_<precision>`), not the metric name.
