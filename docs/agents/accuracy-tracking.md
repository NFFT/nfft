# Accuracy tracking (in-tree HTML reports)

Every CUnit case in the `nfft` / `nfct` / `nfst` suites computes a relative error
`err` and compares it to an analytic `bound`.
The pass/fail gate stays in C (`err < bound`); the accuracy report additionally
visualizes the *value* so slow drift is visible before a case crosses its bound,
and a PR diff flags improvements/regressions. See
[ADR-0004](../adr/0004-in-tree-html-accuracy-reports.md) and the
*Accuracy tracking* vocabulary in [`CONTEXT.md`](../../CONTEXT.md).

## Pipeline

1. **Emit (C).** `accuracy_log_append` (`tests/accuracy_log.c`) appends one raw
   NDJSON record per case to `$NFFT_BENCH_OUT`. No-op when unset, so an ordinary
   `make check` is unaffected. The serial `checkall` writes `$NFFT_BENCH_OUT`; the
   OpenMP `checkall_threads` writes `$NFFT_BENCH_OUT.threads` (separate file, no
   interleaving) and tags its records `"openmp": 1`.
2. **Aggregate (Python).** `tests/accuracy/ndjson_to_bmf.py` groups by the
   *error-shaping parameters* and collapses the *bound-absorbed* `N`/`M` via `max`,
   emitting the accuracy JSON (an in-tree format historically called "BMF") with
   three measures per metric: `accuracy-digits` = `-log10(max(err))` (primary,
   **higher = better**), `max-error` = `max(err)` (secondary, the exact worst
   error), and `bound-digits` = `-log10(bound)` of the worst-error record (so the
   heatmap can color by the *margin* beyond each case's own bound). The `-log10`
   floor sits below any representable bound so a perfect (`err == 0`) case never
   colors as a false failure. Metric name:
   `<module>/<runtime>/<oracle>/<speed>/<direction>/<dim>d/<init-slug>`
   (`runtime` = `serial` | `omp`).
3. **Publish (CI).** The existing `make check` in
   [`.github/workflows/build-linux.yml`](../../.github/workflows/build-linux.yml)
   runs with `NFFT_BENCH_OUT` set, so each gcc window×precision cell emits NDJSON
   as a byproduct, converts it, and publishes an `accuracy-bmf-<BUILD_CONFIG>`
   artifact — **no secret involved**. The `accuracy-report` job consumes those
   artifacts and renders the HTML report (below). No external service is involved.

## Reporting

Human-facing reporting is the `accuracy-report` job in `build-linux.yml`, built
from pure stdlib modules in `tests/accuracy/`:

- **`diff.py`** — compares two sets of per-testbed BMFs on `accuracy-digits`;
  a case is *changed* when `|Δ digits| ≥ 0.5` (configurable).
- **`htmlreport.py`** — renders one **self-contained HTML document** (stdlib only,
  no matplotlib/numpy, no JavaScript). CSS-only **window tabs**; within each tab,
  one HTML `<table>` per module; the per-window testbeds (today the three
  precisions) are columns; the within-module submetrics are rows. Each cell is
  CSS-colored by its **margin** (`accuracy-digits − bound-digits`, correct digits
  beyond the bound — red = over bound, yellow = barely passing, green ramp =
  healthy headroom, grey = missing, `∞` = exact). In PR mode a **Changes** section
  on top lists improved/regressed (capped at 10/group) + an added/removed count
  line, and changed cells are marked (▲/▼, Δ in the tooltip).
- **`report.py`** — the Check summary and the upserted PR comment body (a
  per-module improved/regressed rollup + a link to the full HTML report).
- **`dashboard.py` / `pr_report.py`** — CLIs the workflow runs.

On **develop** push: the HTML dashboard (`index.html`) + `baseline/*.bmf.json` are
published to the `gh-pages` branch root, served by **GitHub Pages** at
`https://<owner>.github.io/<repo>/`. On a **PR**: the baseline is fetched from
`gh-pages`, diffed against the PR's BMFs, the PR report is published to
`gh-pages/pr/<n>/` (served at `…github.io/<repo>/pr/<n>/`), and a non-failing
**Check** plus an always-**upserted** comment (per-module rollup + report link) are
posted.

- **No baseline yet** (first PR, or any PR before `develop` has published): the
  comment shows a "no baseline yet" note + report link — never a misleading
  "unchanged".
- **Fork PRs** get a read-only token, so the in-build `accuracy-report` job is
  **skipped** for them; a trusted `workflow_run` companion
  (`accuracy-report-fork.yml`) instead posts the **aggregate comment + Check** (no
  `gh-pages` archive; the link points at the develop dashboard). It runs
  default-branch code over the fork's BMF artifacts (data only).
- The report depends on GitHub Pages serving the `gh-pages` branch; if Pages is
  disabled the comment link 404s (the per-module rollup is Pages-independent).
- The report steps are `continue-on-error` and the Check is `neutral`, so a publish
  or API failure is **never** able to red CI.

## Run it locally

```bash
NFFT_BENCH_OUT="$PWD/tests/accuracy.ndjson" tests/checkall > /dev/null
# Optional: also capture the OpenMP build (writes the .threads file):
NFFT_BENCH_OUT="$PWD/tests/accuracy.ndjson" tests/checkall_threads > /dev/null
cat tests/accuracy.ndjson tests/accuracy.ndjson.threads > tests/accuracy.all.ndjson 2>/dev/null \
  || cp tests/accuracy.ndjson tests/accuracy.all.ndjson
mkdir -p /tmp/bmf
uv run python -m tests.accuracy.ndjson_to_bmf tests/accuracy.all.ndjson \
  /tmp/bmf/local_gcc_kaiserbessel_double.bmf.json
uv run python -m tests.accuracy.dashboard /tmp/bmf /tmp/site
# open /tmp/site/index.html
```

## Conventions / gotchas

- **`checkall` is a `check_PROGRAM`** — `make all` does not build it; use
  `make -C tests checkall` (or `make check`). `checkall_threads` needs
  `--enable-openmp` configured.
- **Never gates CI.** The report is informational; the C `err < bound` assertion
  is the only gate. `accuracy-digits` is higher-is-better.
- **Serial vs OpenMP.** Both are tracked as distinct metrics (the `runtime` axis);
  the OpenMP binary routes to `<file>.threads` and tags `openmp: 1`. The serial
  results are the deterministic baseline (fixed `SEED`).
- **Metric-name stability.** Changing the grouping key breaks the PR diff — it
  joins PR-vs-baseline BMFs by metric name (`diff.py`), so renamed metrics show as
  spurious added/removed instead of a comparison.
- **Speed axis granularity.** The `fast` speed axis does not distinguish the
  dimension-specialized kernels (`trafo_1d/2d/3d`, `adjoint_1d/2d/3d`) from the
  generic guru transform (`trafo`/`adjoint`); within a group their errors are
  `max`-merged. The `init` variant usually separates them in practice (the
  dedicated 1D path pairs with `init_1d`, the guru path with `init_guru …`).
- **Scope.** Only nfft/nfct/nfst emit today. To add a module, give its harness the
  same `#include "accuracy_log.h"` + one `accuracy_log_append(...)` call.
