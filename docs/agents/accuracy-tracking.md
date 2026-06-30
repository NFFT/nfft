# Accuracy tracking with Bencher

Every CUnit case in the `nfft` / `nfct` / `nfst` suites computes a relative error
`err` and compares it to an analytic `bound` (`tests/nfft.c:393` and the clones).
The pass/fail gate stays in C (`err < bound`); Bencher additionally tracks the
*value* over time so slow drift is visible before a case crosses its bound. See
[ADR-0004](../adr/0004-accuracy-tracking-with-bencher.md) and the
*Accuracy tracking* vocabulary in [`CONTEXT.md`](../../CONTEXT.md).

## Pipeline

1. **Emit (C).** `bench_emit_accuracy` (`tests/bench_emit.c`) appends one raw
   NDJSON record per case to `$NFFT_BENCH_OUT`. No-op when unset, so an ordinary
   `make check` is unaffected. The serial `checkall` writes `$NFFT_BENCH_OUT`; the
   OpenMP `checkall_threads` writes `$NFFT_BENCH_OUT.threads` (separate file, no
   interleaving) and tags its records `"openmp": 1`.
2. **Aggregate (Python).** `tests/bench/ndjson_to_bmf.py` groups by the
   *error-shaping parameters* and collapses the *bound-absorbed* `N`/`M` via `max`,
   emitting BMF with two measures per metric: `accuracy-digits` = `-log10(max(err))`
   (primary, **higher = better**) and `max-error` = `max(err)` (secondary, the exact
   worst error). Metric name:
   `<module>/<runtime>/<oracle>/<speed>/<direction>/<dim>d/<init-slug>`
   (`runtime` = `serial` | `omp`).
3. **Upload (CI).** The existing `make check` in
   [`.github/workflows/build-linux.yml`](../../.github/workflows/build-linux.yml)
   runs with `NFFT_BENCH_OUT` set, so each gcc window×precision cell emits NDJSON
   as a byproduct, converts it, and publishes a `accuracy-bmf-<BUILD_CONFIG>`
   artifact — **no secret involved**. A separate, environment-gated
   `bencher-upload` job downloads those artifacts and calls
   `bencher run --project nfft --adapter json --file …` track-only (no thresholds,
   no `--err`), `develop` baseline, start-point recipe on PRs. The tests are **not
   re-run** just to produce Bencher data; only the upload is gated.

   `BENCHER_API_TOKEN` is an **environment secret** (not repo-wide; it holds the
   project-scoped `bencher_run_*` key) held in two environments, so only this job
   can read it, only after the environment's rules pass. **Only fork PRs are
   gated** (the one untrusted upload source); everything else runs unattended:
   - **fork `pull_request`** → `benchmarks` (required reviewer → manual approval);
   - **everything else** (push to default branches, same-repo PR,
     `workflow_dispatch`) → `bencher-baseline` (no reviewer → unattended).

   `bencher-baseline` must **allow all branches** in its deployment-branch policy,
   or same-repo PR feature branches are blocked from it. Without the key a run
   falls back to `--dry-run`.

## Reporting

Bencher is now only the long-term archive. Human-facing reporting is the
`accuracy-report` job in `build-linux.yml`, built from pure modules in
`tests/bench/`:

- **`diff.py`** — compares two sets of per-testbed BMFs on `accuracy-digits`;
  a case is *changed* when `|Δ digits| ≥ 0.5` (configurable).
- **`heatmap.py`** — absolute and relative heatmap PNGs + an inline emoji grid.
- **`report.py`** — the Check summary and the upserted PR comment body.
- **`dashboard.py` / `pr_report.py`** — CLIs the workflow runs.

On **develop** push: the absolute heatmap + `baseline/*.bmf.json` are published
to the `gh-pages` branch (the standing dashboard). On a **PR**: the baseline is
fetched from `gh-pages`, diffed against the PR's BMFs, and a non-failing
**Check** plus an always-upserted **comment** (emoji grid + itemized
improvements/regressions, capped at 10/group, links to the absolute + relative
heatmap PNGs) are posted. Fork PRs get the emoji grid only (no `gh-pages` write).

Scope: P1. The convergence-curve view (err vs N) is P2 — see
[`docs/superpowers/specs/2026-06-30-accuracy-reporting-layer-design.md`](../superpowers/specs/2026-06-30-accuracy-reporting-layer-design.md).

## Run it locally

```bash
NFFT_BENCH_OUT="$PWD/tests/accuracy.ndjson" tests/checkall > /dev/null
# Optional: also capture the OpenMP build (writes the .threads file):
NFFT_BENCH_OUT="$PWD/tests/accuracy.ndjson" tests/checkall_threads > /dev/null
cat tests/accuracy.ndjson tests/accuracy.ndjson.threads > tests/accuracy.all.ndjson 2>/dev/null \
  || cp tests/accuracy.ndjson tests/accuracy.all.ndjson
uv run python tests/bench/ndjson_to_bmf.py tests/accuracy.all.ndjson tests/accuracy.bmf.json
bencher run --dry-run --project nfft --branch local \
  --testbed local --adapter json --file tests/accuracy.bmf.json
```

## Conventions / gotchas

- **`checkall` is a `check_PROGRAM`** — `make all` does not build it; use
  `make -C tests checkall` (or `make check`). `checkall_threads` needs
  `--enable-openmp` configured.
- **Track-only, phased.** No thresholds yet; never `--err`. `accuracy-digits` is
  higher-is-better → add a *lower* boundary when alerts come later (alert if the
  worst-case accurate digits drop).
- **Serial vs OpenMP.** Both are tracked as distinct metrics (the `runtime` axis);
  the OpenMP binary routes to `<file>.threads` and tags `openmp: 1`. The serial
  results are the deterministic baseline (fixed `SEED`).
- **Metric-name stability.** Changing the grouping key breaks Bencher history.
- **Speed axis granularity.** The `fast` speed axis does not distinguish the
  dimension-specialized kernels (`trafo_1d/2d/3d`, `adjoint_1d/2d/3d`) from the
  generic guru transform (`trafo`/`adjoint`); within a group their errors are
  `max`-merged. The `init` variant usually separates them in practice (the
  dedicated 1D path pairs with `init_1d`, the guru path with `init_guru …`).
- **Scope.** Only nfft/nfct/nfst emit today. To add a module, give its harness the
  same `#include "bench_emit.h"` + one `bench_emit_accuracy(...)` call.
- **Quota dial.** Merge file+online in `group_key` (180 → 108 metrics/cell) first.
