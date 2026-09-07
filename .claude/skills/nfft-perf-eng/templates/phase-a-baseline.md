# Phase A — baseline (the exit reference)

- **Build:** `<cmake flags, e.g. -DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON …>` (three trees: `build-cmake{,-f,-l}`)
- **Baseline commit:** <sha>
- **Tests (per precision):** double <N/N> ✅ · float <N/N> ✅ · long double <N/N> ✅ — `artifacts/baseline-tests-{d,f,l}.log`
- **Benchmark cases captured:** <N> × 3 precisions — `artifacts/baseline-bench-{d,f,l}.json`

## Benchmark snapshot

<!-- Benchmark snapshot format — keep these columns; one row per case × precision, ALL cases.
     prec = d (double) / f (float) / l (long double). -->
| prec | case                       | median_ns | stdev_ns | rounds |
|------|----------------------------|-----------|----------|--------|
| d    | <nfft_forward_direct_1d/…> | <123456>  | <789>    | <50>   |
| f    | <nfft_forward_direct_1d/…> | <98765>   | <654>    | <50>   |
| l    | <nfft_forward_direct_1d/…> | <210987>  | <912>    | <50>   |

## Notes

<!-- anything notable about the baseline state (warnings, skipped cases, a precision that is
     correctness-only because its benchmark could not build/run); '—' if nothing -->
—
