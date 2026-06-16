# Phase A — baseline (the exit reference)

- **Build:** three trees `build-cmake{,-f,-l}` — `-DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON -DCMAKE_C_FLAGS="-O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math"`, plus the per-tree precision flag. **Stale in-source `include/config.h` cleared first** (else the float/long-double trees mis-link — see precision-matrix.md).
- **Baseline commit:** a692216c
- **Track:** walltime
- **Tests (per precision):** double 1854/1854 ✅ · float 1843/1843 ✅ · long double 1854/1854 ✅ (`checkall` + `checkall_threads`) — `artifacts/baseline-tests-{d,f,l}.log`
- **Benchmark cases captured:** 22 × 3 precisions — `artifacts/baseline-bench-{d,f,l}.json`

## Benchmark snapshot

<!-- Benchmark snapshot (canonical format) — all 22 cases × 3 precisions (d/f/l), walltime medians over 5 rounds. -->
| prec | case | median_ns | stdev_ns | rounds |
|------|------|-----------|----------|--------|
| d | nfft_forward_direct_1d[32/100] | 7694 | 24 | 5 |
| d | nfft_forward_direct_1d[64/200] | 29324 | 150 | 5 |
| d | nfft_forward_direct_1d[128/400] | 116439 | 265 | 5 |
| d | nfft_forward_direct_1d[256/800] | 465332 | 1207 | 5 |
| d | nfft_forward_direct_1d[512/1600] | 1862450 | 6586 | 5 |
| d | nfft_adjoint_direct_1d[32/100] | 8464 | 107 | 5 |
| d | nfft_adjoint_direct_1d[64/200] | 31395 | 88 | 5 |
| d | nfft_adjoint_direct_1d[128/400] | 124373 | 423 | 5 |
| d | nfft_adjoint_direct_1d[256/800] | 502903 | 1537 | 5 |
| d | nfft_adjoint_direct_1d[512/1600] | 2006180 | 8601 | 5 |
| d | nfft_forward_direct_2d[16/16/500] | 1312200 | 15962 | 5 |
| d | nfft_forward_direct_2d[32/32/1000] | 10305100 | 127132 | 5 |
| d | nfft_forward_direct_2d[64/64/2000] | 66904600 | 224405 | 5 |
| d | nfft_adjoint_direct_2d[16/16/500] | 1260930 | 33903 | 5 |
| d | nfft_adjoint_direct_2d[32/32/1000] | 9881200 | 41632 | 5 |
| d | nfft_adjoint_direct_2d[64/64/2000] | 65500100 | 199041 | 5 |
| d | nfft_forward_direct_3d[4/4/4/250] | 96383 | 1644 | 5 |
| d | nfft_forward_direct_3d[8/8/8/500] | 3209680 | 9131 | 5 |
| d | nfft_forward_direct_3d[16/16/16/1000] | 42093800 | 108844 | 5 |
| d | nfft_adjoint_direct_3d[4/4/4/250] | 98346 | 1056 | 5 |
| d | nfft_adjoint_direct_3d[8/8/8/500] | 3081320 | 20590 | 5 |
| d | nfft_adjoint_direct_3d[16/16/16/1000] | 40254800 | 165137 | 5 |
| f | nfft_forward_direct_1d[32/100] | 7779 | 56 | 5 |
| f | nfft_forward_direct_1d[64/200] | 32829 | 151 | 5 |
| f | nfft_forward_direct_1d[128/400] | 133304 | 322 | 5 |
| f | nfft_forward_direct_1d[256/800] | 541102 | 1403 | 5 |
| f | nfft_forward_direct_1d[512/1600] | 2222930 | 5656 | 5 |
| f | nfft_adjoint_direct_1d[32/100] | 7849 | 33 | 5 |
| f | nfft_adjoint_direct_1d[64/200] | 33177 | 113 | 5 |
| f | nfft_adjoint_direct_1d[128/400] | 134820 | 246 | 5 |
| f | nfft_adjoint_direct_1d[256/800] | 547651 | 969 | 5 |
| f | nfft_adjoint_direct_1d[512/1600] | 2348450 | 7476 | 5 |
| f | nfft_forward_direct_2d[16/16/500] | 426308 | 1957 | 5 |
| f | nfft_forward_direct_2d[32/32/1000] | 3382850 | 10489 | 5 |
| f | nfft_forward_direct_2d[64/64/2000] | 25346000 | 412201 | 5 |
| f | nfft_adjoint_direct_2d[16/16/500] | 436553 | 1053 | 5 |
| f | nfft_adjoint_direct_2d[32/32/1000] | 3514870 | 25921 | 5 |
| f | nfft_adjoint_direct_2d[64/64/2000] | 26436400 | 191134 | 5 |
| f | nfft_forward_direct_3d[4/4/4/250] | 61852 | 270 | 5 |
| f | nfft_forward_direct_3d[8/8/8/500] | 936773 | 7195 | 5 |
| f | nfft_forward_direct_3d[16/16/16/1000] | 13284100 | 71192 | 5 |
| f | nfft_adjoint_direct_3d[4/4/4/250] | 61621 | 385 | 5 |
| f | nfft_adjoint_direct_3d[8/8/8/500] | 952768 | 5189 | 5 |
| f | nfft_adjoint_direct_3d[16/16/16/1000] | 13699100 | 24147 | 5 |
| l | nfft_forward_direct_1d[32/100] | 1895590 | 20273 | 5 |
| l | nfft_forward_direct_1d[64/200] | 7985810 | 20594 | 5 |
| l | nfft_forward_direct_1d[128/400] | 32886100 | 88815 | 5 |
| l | nfft_forward_direct_1d[256/800] | 131978000 | 944229 | 5 |
| l | nfft_forward_direct_1d[512/1600] | 528143000 | 1973650 | 5 |
| l | nfft_adjoint_direct_1d[32/100] | 1919880 | 6589 | 5 |
| l | nfft_adjoint_direct_1d[64/200] | 8150290 | 65001 | 5 |
| l | nfft_adjoint_direct_1d[128/400] | 32911800 | 74119 | 5 |
| l | nfft_adjoint_direct_1d[256/800] | 131517000 | 334137 | 5 |
| l | nfft_adjoint_direct_1d[512/1600] | 527667000 | 2848540 | 5 |
| l | nfft_forward_direct_2d[16/16/500] | 82915800 | 446961 | 5 |
| l | nfft_forward_direct_2d[32/32/1000] | 668311000 | 2264120 | 5 |
| l | nfft_forward_direct_2d[64/64/2000] | 5356080000 | 23438200 | 5 |
| l | nfft_adjoint_direct_2d[16/16/500] | 81763400 | 1651460 | 5 |
| l | nfft_adjoint_direct_2d[32/32/1000] | 664817000 | 3323220 | 5 |
| l | nfft_adjoint_direct_2d[64/64/2000] | 5325840000 | 9551960 | 5 |
| l | nfft_forward_direct_3d[4/4/4/250] | 9549980 | 53051 | 5 |
| l | nfft_forward_direct_3d[8/8/8/500] | 161259000 | 283571 | 5 |
| l | nfft_forward_direct_3d[16/16/16/1000] | 2623360000 | 5096900 | 5 |
| l | nfft_adjoint_direct_3d[4/4/4/250] | 9409850 | 25130 | 5 |
| l | nfft_adjoint_direct_3d[8/8/8/500] | 160025000 | 555912 | 5 |
| l | nfft_adjoint_direct_3d[16/16/16/1000] | 2625900000 | 3605150 | 5 |

## Notes

- **Long double is ~250–280× slower than double** on this target (e.g. forward_direct_1d[32/100]:
  7694 ns d vs 1,895,590 ns l) because `COSL`/`SINL` are extended-precision and far costlier than
  the double libm calls. Float is marginally *slower* than double here (no float-specific trig win,
  same call count). This precision spread is exactly what the matrix exists to surface, and it makes
  the cos/sin → complex-multiply recurrence especially valuable in long double.
- Benchmark coverage is the **direct** transforms only (forward/adjoint, 1d/2d/3d). The fault/metric
  pinning (Phases B/C) is done in double and applies to all precisions (same suite, same cases).
- float reports 1843 test cases vs 1854 for double/long double — a handful of suite cases are
  precision-gated; all run cases pass.
