# Phase A — baseline (the exit reference)

- **Build:** `<cmake flags, e.g. -DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON …>`
- **Baseline commit:** <sha>
- **Track:** <walltime | simulation>
- **Tests:** <N passed / N total> — green ✅ — `artifacts/baseline-tests.log`
- **Benchmark cases captured:** <N> — `artifacts/baseline-bench.json`

## Benchmark snapshot

<!-- Benchmark snapshot format — keep these columns; one row per case, ALL cases. -->
| case                       | median_ns | stdev_ns | rounds |
|----------------------------|-----------|----------|--------|
| <nfft_forward_direct_1d/…> | <123456>  | <789>    | <50>   |

## Notes

<!-- anything notable about the baseline state (warnings, skipped cases); '—' if nothing -->
—
