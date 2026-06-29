# Track test accuracy in Bencher, aggregated by error-shaping parameters

We track the per-case accuracy figures (`err` vs `bound`) that the `nfft`/`nfct`/`nfst`
CUnit suites already compute, uploading them to the **Bencher** project `nfft` (track-only,
separate from the CodSpeed instruction-count benchmarks). Rather than one Bencher series
per case (~1456/cell × 12 cells), we aggregate: one **accuracy metric** per combination of
**error-shaping parameters** (window, precision, runtime serial/OpenMP, dimension,
direct/fast, forward/adjoint, init variant, file/online), collapsing the **bound-absorbed
parameters** (`N`, `M`) via `max`. The uploaded measures are **accuracy-digits**
`-log10(max(err))` (primary, higher = better) plus **max-error** `max(err)` raw (secondary).

## Status

accepted

## Considered Options

- **Per-case series.** Finest resolution, but ~20k metrics/run on Bencher Cloud (quota
  pressure) and a dashboard dominated by ~190 near-`ε` round-off-floor noise lines per
  cell. Rejected: cost and illegibility for near-zero signal.
- **Aggregate, merge file+online.** ~108 metrics/cell. Rejected: `max` over the merged
  group lets a file case mask a regression in the larger-`N` **online** cases — the
  highest-value, most-realistic sizes — because online uses a different oracle (the C
  direct transform) and is therefore an *error-shaping* difference, not merely a larger
  `N`. We split file vs online (~180/cell) to preserve that signal.
- **Which value to upload.** Initially `max(err/bound)` (a "tightness ratio"). Replaced:
  the bound is a fixed analytic function, so per metric the ratio carries the same
  regression signal as the raw error (Bencher's per-metric thresholds don't need the
  normalization), and it was as illegible as the raw error in the UI. Raw errors span
  ~14 orders of magnitude (1e-17 .. 1e-3), which no *linear* display renders well (all
  round to `0.00`). So the primary is **accuracy-digits** = `-log10(max(err))` — the
  worst-case accurate digits, which reads cleanly (~3 .. 18) and where a regression
  *lowers* the value (Bencher **lower-boundary** threshold). The exact worst error is
  kept as the secondary **max-error** for the precise figure.

## Consequences

- **Metric-name stability matters** (as with CodSpeed benchmark names): the grouping key
  fixes the Bencher series identity, so changing the error-shaping/​bound-absorbed split
  later breaks historical continuity. This is why the split is recorded here.
- The C harness stays *dumb* — it emits one raw, structured NDJSON record per case behind
  the `NFFT_BENCH_OUT` env var (a no-op when unset). **All grouping and `max` aggregation
  live in a Python converter**, so the policy is tunable without touching C or re-running
  builds.
- **Emission piggybacks on the existing test run; only the upload is gated.** Rather than a
  dedicated workflow that re-builds and re-runs the 12-cell gcc matrix, the emission is
  switched on inside `build-linux.yml`'s existing `make check` (via `NFFT_BENCH_OUT`); each
  cell publishes its BMF as an artifact with no secret. A separate `bencher-upload` job —
  gated on the protected `benchmarks` environment for PRs, automatic on default-branch
  pushes — is the *only* step that needs the API key. Tests run exactly once.
- **Serial and OpenMP are both tracked, not duplicated.** `make check` runs `checkall`
  (serial) and `checkall_threads` (OpenMP); the parallel reduction order perturbs the low
  bits, so `runtime` (serial/omp) is an error-shaping axis. The OpenMP binary routes its
  records to a separate `<file>.threads` path (no interleaving), and the converter gives
  serial and omp distinct metric names.
- **Track-only, phased.** Pure recording first (no thresholds — they misfire on a cold
  baseline); non-blocking upper-boundary alerts on the ratio once `develop` has baseline
  history; hard CI gating (`--err`) only later, by explicit decision. CUnit's `err < bound`
  remains the sole gate meanwhile.
- Scope is `nfft`/`nfct`/`nfst` (the three structurally identical harnesses); other
  modules are follow-up. Window/precision are the **accuracy testbed**, not the metric
  name.
