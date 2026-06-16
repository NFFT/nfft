# Phase D — inner loop (iteration journal)

- **B-net:** `nfft` (116 cases; `checkall` + `checkall_threads`, **all three precisions**) · **C-metric:** `nfft_forward_direct_1d[*]`
- **Current change:** `artifacts/change.diff`

## The idea (discovered here)

The 1d inner loop evaluates `e^{-i·omega_k}` with `omega_k = 2π(k − N/2)·x[j]` from scratch each
iteration — two transcendentals per term. `omega_k` advances by a constant `Δ = 2π·x[j]`, so the
kernel is a geometric sequence `w_{k+1} = w_k · e^{-iΔ}`: evaluate the two trig calls once per
node, carry the phase by one complex multiply. Pattern: **phase recurrence / strength reduction**.
Risk: rounding drift over `N` steps — guarded by the B-net tolerance, **checked in every precision**.

## Iterations

<!-- Iteration journal (canonical format). net = green in all three precisions. -->
| iter | change | net | metric median (ns) before→after |
|------|--------|-----|----------------------------------|
| 1 | phase recurrence: hoist `COS/SIN` per node, carry `w *= dw` in the 1d loop | green d/f/l (d 1854/0, f 1843/0, l 1854/0) | see per-precision table below |

### Per-precision speedup at iter 1 (forward_direct_1d median, baseline → after)

| prec | [32/100] | [128/400] | [512/1600] | speedup range |
|------|----------|-----------|------------|---------------|
| d (double)      | 7694 → 5276       | 116439 → 92213      | 1862450 → 1463150     | **1.26–1.46×** |
| f (float)       | 7779 → 4483       | 133304 → 86394      | 2222930 → 1417260     | **1.54–1.74×** |
| l (long double) | 1895590 → 517441  | 32886100 → 8698350  | 528143000 → 132877000 | **3.66–3.97×** |

**The win is precision-dependent** — and largest exactly where the baseline is worst: long-double
`COSL`/`SINL` are ~250× the double cost, so removing them per-iteration is a ~4× win, vs ~1.3× in
double. A double-only loop would have reported only the ~1.3× and missed this entirely. The drift
risk did **not** materialise: the net stays green at every size in all three precisions.

## Exit state

- B-net green at latest kept state, **all precisions**: **yes** (d 1854/0, f 1843/0, l 1854/0; `checkall_threads` clean)
- metric median dropped beyond noise: **yes** (d ~1.3×, f ~1.6×, l ~3.9×; stdev <1.5%)
- `artifacts/change.diff` matches latest kept state: **yes**
