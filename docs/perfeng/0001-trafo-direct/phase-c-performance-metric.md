# Phase C — performance metric

- **Outcome:** metric pinned ✅
- **Benchmark filter:** `--benchmark_filter='nfft_forward_direct_1d.*'` (target = 1d forward; `nfft_adjoint_direct_1d.*` as control)
- **Slowdown injected:** 1d `trafo_direct` body wrapped in a 10× repeat loop (`v` reset each rep, results identical) — `artifacts/slowdown.diff`

## Target baseline (double)

| prec | case | median_ns | stdev_ns | rounds |
|------|------|-----------|----------|--------|
| d | nfft_forward_direct_1d[32/100]   | 7694    | 24  | 5 |
| d | nfft_forward_direct_1d[128/400]  | 116439  | 265 | 5 |
| d | nfft_forward_direct_1d[512/1600] | 1862450 | —   | 5 |

## Metric — cases that moved under the slowdown (double)

| case | median_ns base | median_ns slowed | moves the target? |
|------|----------------|------------------|-------------------|
| nfft_forward_direct_1d[32/100]   | 7694    | 72919    | yes — 9.5× |
| nfft_forward_direct_1d[64/200]   | 29324   | 292944   | yes — 10.0× |
| nfft_forward_direct_1d[128/400]  | 116439  | 1162350  | yes — 10.0× |
| nfft_forward_direct_1d[256/800]  | 465332  | 4639490  | yes — 10.0× |
| nfft_forward_direct_1d[512/1600] | 1862450 | 18496900 | yes — 9.9× |
| nfft_adjoint_direct_1d[*] (control) | (5 cases) | unchanged | no — 1.0× |

**Metric for the inner loop:** the five `nfft_forward_direct_1d[*]` cases, measured **in all
three precisions** in Phase D. Pinned in double; the metric is the same case set per precision
(the long-double baseline shows these cases are ~250× more expensive — the recurrence should pay
off most there).

## Revert confirmation

- `git diff` empty: **yes**
