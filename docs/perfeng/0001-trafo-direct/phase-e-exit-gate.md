# Phase E — exit gate

- **Verdict:** ⚠ **partial / landed with caveat** — target improves in all three precisions and the net is green everywhere, but **one confirmed float regression** on an untouched control case. A reviewer must weigh ship-vs-fix-vs-revert.
- **Target metric (forward_direct_1d):** double ~1.3× · float ~1.6× · **long double ~3.9×** faster
- **Landed change:** `artifacts/change.diff` (kernel/nfft/nfft.c, the 1d phase recurrence)

## Comparison (baseline → final, per precision)

<!-- Comparison table (canonical format, prec column). Noise rule: regressed only if median
     rises past max(3·stdev, 2% of base) AND survives a re-run. The 3 float flags were re-run:
     only adjoint_1d[32/100] survived. long double = 1d family full + (2d/3d trimmed, see notes). -->
| prec | case | base median_ns | final median_ns | Δ% | verdict |
|------|------|----------------|-----------------|----|---------|
| d | nfft_adjoint_direct_1d[128/400] | 124373 | 124247 | -0.1% | noise |
| d | nfft_adjoint_direct_1d[256/800] | 502903 | 500677 | -0.4% | noise |
| d | nfft_adjoint_direct_1d[32/100] | 8464 | 8491 | +0.3% | noise |
| d | nfft_adjoint_direct_1d[512/1600] | 2006180 | 2003650 | -0.1% | noise |
| d | nfft_adjoint_direct_1d[64/200] | 31395 | 31379 | -0.1% | noise |
| d | nfft_adjoint_direct_2d[16/16/500] | 1260930 | 1251880 | -0.7% | noise |
| d | nfft_adjoint_direct_2d[32/32/1000] | 9881200 | 9918570 | +0.4% | noise |
| d | nfft_adjoint_direct_2d[64/64/2000] | 65500100 | 65184700 | -0.5% | noise |
| d | nfft_adjoint_direct_3d[16/16/16/1000] | 40254800 | 40063600 | -0.5% | noise |
| d | nfft_adjoint_direct_3d[4/4/4/250] | 98346 | 97820 | -0.5% | noise |
| d | nfft_adjoint_direct_3d[8/8/8/500] | 3081320 | 3043050 | -1.2% | noise |
| d | nfft_forward_direct_1d[128/400] | 116439 | 91781 | -21.2% | faster |
| d | nfft_forward_direct_1d[256/800] | 465332 | 365396 | -21.5% | faster |
| d | nfft_forward_direct_1d[32/100] | 7694 | 5204 | -32.4% | faster |
| d | nfft_forward_direct_1d[512/1600] | 1862450 | 1468760 | -21.1% | faster |
| d | nfft_forward_direct_1d[64/200] | 29324 | 22703 | -22.6% | faster |
| d | nfft_forward_direct_2d[16/16/500] | 1312200 | 1274770 | -2.9% | noise |
| d | nfft_forward_direct_2d[32/32/1000] | 10305100 | 9971300 | -3.2% | noise |
| d | nfft_forward_direct_2d[64/64/2000] | 66904600 | 65103600 | -2.7% | faster |
| d | nfft_forward_direct_3d[16/16/16/1000] | 42093800 | 40236200 | -4.4% | faster |
| d | nfft_forward_direct_3d[4/4/4/250] | 96383 | 94932 | -1.5% | noise |
| d | nfft_forward_direct_3d[8/8/8/500] | 3209680 | 3078380 | -4.1% | faster |
| f | nfft_adjoint_direct_1d[128/400] | 134820 | 134708 | -0.1% | noise |
| f | nfft_adjoint_direct_1d[256/800] | 547651 | 559795 | +2.2% | noise (re-run +0.2%, did not survive) |
| f | nfft_adjoint_direct_1d[32/100] | 7849 | 8508 | +8.4% | ⚠ REGRESSED (confirmed: stable +8% vs reverted baseline; layout) |
| f | nfft_adjoint_direct_1d[512/1600] | 2348450 | 2295190 | -2.3% | faster |
| f | nfft_adjoint_direct_1d[64/200] | 33177 | 33214 | +0.1% | noise |
| f | nfft_adjoint_direct_2d[16/16/500] | 436553 | 445324 | +2.0% | noise (re-run +0.5%, did not survive) |
| f | nfft_adjoint_direct_2d[32/32/1000] | 3514870 | 3507950 | -0.2% | noise |
| f | nfft_adjoint_direct_2d[64/64/2000] | 26436400 | 26338700 | -0.4% | noise |
| f | nfft_adjoint_direct_3d[16/16/16/1000] | 13699100 | 13757100 | +0.4% | noise |
| f | nfft_adjoint_direct_3d[4/4/4/250] | 61621 | 62228 | +1.0% | noise |
| f | nfft_adjoint_direct_3d[8/8/8/500] | 952768 | 956795 | +0.4% | noise |
| f | nfft_forward_direct_1d[128/400] | 133304 | 88442 | -33.7% | faster |
| f | nfft_forward_direct_1d[256/800] | 541102 | 357644 | -33.9% | faster |
| f | nfft_forward_direct_1d[32/100] | 7779 | 4604 | -40.8% | faster |
| f | nfft_forward_direct_1d[512/1600] | 2222930 | 1439690 | -35.2% | faster |
| f | nfft_forward_direct_1d[64/200] | 32829 | 20858 | -36.5% | faster |
| f | nfft_forward_direct_2d[16/16/500] | 426308 | 422429 | -0.9% | noise |
| f | nfft_forward_direct_2d[32/32/1000] | 3382850 | 3377610 | -0.2% | noise |
| f | nfft_forward_direct_2d[64/64/2000] | 25346000 | 25269700 | -0.3% | noise |
| f | nfft_forward_direct_3d[16/16/16/1000] | 13284100 | 13277100 | -0.1% | noise |
| f | nfft_forward_direct_3d[4/4/4/250] | 61852 | 61705 | -0.2% | noise |
| f | nfft_forward_direct_3d[8/8/8/500] | 936773 | 935945 | -0.1% | noise |
| l | nfft_adjoint_direct_1d[128/400] | 32911800 | 31152400 | -5.3% | faster |
| l | nfft_adjoint_direct_1d[256/800] | 131517000 | 124886000 | -5.0% | faster |
| l | nfft_adjoint_direct_1d[32/100] | 1919880 | 1817470 | -5.3% | faster |
| l | nfft_adjoint_direct_1d[512/1600] | 527667000 | 501199000 | -5.0% | faster |
| l | nfft_adjoint_direct_1d[64/200] | 8150290 | 7602170 | -6.7% | faster |
| l | nfft_forward_direct_1d[128/400] | 32886100 | 8716830 | -73.5% | faster |
| l | nfft_forward_direct_1d[256/800] | 131978000 | 34083800 | -74.2% | faster |
| l | nfft_forward_direct_1d[32/100] | 1895590 | 527033 | -72.2% | faster |
| l | nfft_forward_direct_1d[512/1600] | 528143000 | 132684000 | -74.9% | faster |
| l | nfft_forward_direct_1d[64/200] | 7985810 | 2192720 | -72.5% | faster |

## Exit conditions (float · double · long double)

1. **Full test suite passes as in Phase A, every precision:** ✅ PASS — `ctest` 100% in d/f/l (`checkall` + `checkall_threads`) — `artifacts/final-tests-{d,f,l}.log`
2. **No benchmark regresses beyond noise, every case/precision:** ⚠ **QUALIFIED FAIL** — double 0/22 and long-double 0/10 clean; **float has 1 confirmed regression**: `adjoint_direct_1d[32/100]` +8% (untouched function; stable ~8480 ns vs stable reverted baseline ~7840 ns over 4+3 runs ⇒ a code-layout / i-cache effect of enlarging the sibling `trafo_direct`, not an algorithmic change). Two other float flags (`adjoint_1d[256/800]` +2.2%, `adjoint_2d[16/16/500]` +2.0%) were **noise** — they did not survive a re-run.
3. **Deterministic cross-check (simulation):** N/A — no CodSpeed this session; walltime-only. CI's CodSpeed run is the final authority.
4. **`git diff` is only the intended change:** ✅ PASS — kernel/nfft/nfft.c, the 1d recurrence only.

## Notes / trims

- **long double 2d/3d were not benchmarked** (each ~250× the double cost → tens of minutes). The
  change is strictly in the 1d branch; the multivariate (2d/3d) path is untouched and is covered
  at full fidelity in double and float (0 regressions there). long-double **tests** ran in full
  (1854/1854). Disclosed per the skill's trim rule.
- The float `adjoint_1d[32/100]` regression is **precision-specific**: in double the same case is
  noise (±0%). Double-only validation would have shipped this blind — this is precisely the value
  of the precision matrix.

## Verdict

The phase recurrence is a large, real win on the target in **every** precision — and a **~3.9×**
win in long double, where `COSL`/`SINL` dominate — with the 116-case correctness net green in
float, double, and long double. The exit gate does **not** pass clean: a reproducible ~8% float
slowdown on the smallest *adjoint* (untouched) micro-benchmark, traced to code layout, trips
condition 2. This is a judgment call for a human reviewer: ship (the regressed case is an 8µs
operation; the target gains 1.3–3.9×), attempt a layout-neutral formulation, or revert. Recorded
as **landed-with-caveat**, not a clean pass.
