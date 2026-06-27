# Track test accuracy in Bencher, aggregated by error-shaping parameters

We track the per-case accuracy figures (`err` vs `bound`) that the `nfft`/`nfct`/`nfst`
CUnit suites already compute, uploading them to **Bencher** (track-only, separate from
the CodSpeed instruction-count benchmarks). Rather than one Bencher series per case
(~1680/cell × 12 cells ≈ 20k metrics/run), we aggregate: one **accuracy metric** per
combination of **error-shaping parameters** (window, precision, dimension, direct/fast,
forward/adjoint, init variant, file/online), collapsing the **bound-absorbed parameters**
(`N`, `M`) via `max`. The uploaded value is the **tightness ratio** `max(err/bound)`
(primary) plus `max(err)` raw (secondary). This yields ~180 metrics/cell (~4.3k/run).

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
- **Raw `err` only / ratio only.** We upload both: the **tightness ratio** is primary
  (normalizes out the bound-absorbed `N`-dependence and reads as distance-to-failure),
  raw `err` secondary (absolute context). Same group count, so cheap.

## Consequences

- **Metric-name stability matters** (as with CodSpeed benchmark names): the grouping key
  fixes the Bencher series identity, so changing the error-shaping/​bound-absorbed split
  later breaks historical continuity. This is why the split is recorded here.
- The C harness stays *dumb* — it emits one raw, structured NDJSON record per case behind
  the `NFFT_BENCH_OUT` env var (a no-op when unset). **All grouping and `max` aggregation
  live in a Python converter**, so the policy is tunable without touching C or re-running
  builds.
- **Track-only, phased.** Pure recording first (no thresholds — they misfire on a cold
  baseline); non-blocking upper-boundary alerts on the ratio once `develop` has baseline
  history; hard CI gating (`--err`) only later, by explicit decision. CUnit's `err < bound`
  remains the sole gate meanwhile.
- Scope is `nfft`/`nfct`/`nfst` (the three structurally identical harnesses); other
  modules are follow-up. Window/precision are the **accuracy testbed**, not the metric
  name.
