# 02 — Take `M` and choose `(n, m)` by cost

Status: done; see also issue 03, which followed from it

## The defect

`X(tune_plan)` minimises the oversampled size first and the cut-off second.
That is a memory-first policy. It is the wrong policy for run time whenever the
node count `M` is large, because the two costs pull in opposite directions:
the FFT is `O(n log n)`, the node convolution is `O(M m)`.

Measured against the legacy geometry `n = 2*next_power_of_2(N)` with an oracle
cut-off search (see `docs/tuning-vs-legacy.md`, 288 configurations):

| shape | M | median speedup | tuner faster |
|---|---|---|---|
| N-dominated | N/4 | 1.31x | 95/96 |
| balanced | N | 0.96x | 42/96 |
| M-dominated | 4N | 0.77x | 3/96 |

The tuner lands at about **0.62×** the legacy `n` and pays with a cut-off **2
to 3 larger**. Accuracy is unaffected — 288/288 goals met either way.

## Change

Add `M` and pick the pair by predicted cost:

```c
int X(tune_plan)(NFFT_INT N, NFFT_INT M, int adjoint, R goal,
                 NFFT_INT *n, int *m, R *attained);
```

This breaks the current signature. `tune_plan` is new and unreleased, so change
it in place rather than adding a variant. `X(tune)` and `X(tune_sigma)` keep
their signatures — they answer questions that genuinely do not involve `M`.

### Algorithm

1. Cap the goal as today: `cap = floor at the top of the band`,
   `a = max(goal, cap)`.
2. Enumerate the candidate sizes: every **even 5-smooth** `n` in
   `[tune_band_lo(N), 4*N]`. There are only tens of these
   (5-smooth numbers up to `x` number `~(ln x)^3 / (6 ln2 ln3 ln5)`), so
   enumerate them all rather than bisecting. `tune_next_smooth` in
   `kernel/nfft/tune.c` already generates them one at a time; generalise it to
   walk the set, or call it repeatedly.
3. For each candidate, the smallest feasible cut-off is `tune_scan(..., a, ...)`.
   Skip candidates where no `m <= tune_m_max(n)` reaches `a`.
4. Keep the candidate minimising the cost model below.
5. Return that pair, `*attained = E(m, n/N)`, and `1` when `*attained <= goal`.

`2*next_power_of_2(N)` is a power of two and lies in `[2N, 4N)`, so the legacy
size is **always already in the candidate set**. The search can therefore never
pick something the model rates worse than legacy's own geometry.

### The cost model

The repo already has one, in `kernel/nfft/nfft-nd.c:58` (`pcost` for the
native fast solver). For `d = 1` it reduces to

```
cost(n, m) = 2*N + 5*n*log2(n) + 2*M*(2m+2)
```

**Do not use those weights unchanged.** Fitting `t = a*n*log2(n) + b*M*(2m+2)`
to the 288 measured pairs gives `b/a ≈ 1.4` (double, long double) and `≈ 1.9`
(float), against the `2.0/5.0 = 0.4` the planner assumes. The planner
under-weights the node convolution by roughly 3.5×. Ranking accuracy on the
measured pairs:

| weights | ranks the measured pair correctly |
|---|---|
| planner's `b/a = 0.4` | 197/288 (68 %) — 99 % N-dominated, 44 % balanced, 62 % M-dominated |
| fitted `b/a ≈ 1.4` | 86–93 % per precision |

Correlation between predicted and measured speedup ratio is `r = 0.86` even
with the bad weights, so the model shape is right and only the ratio is off.
Only the ratio matters for ranking; the overall scale cancels.

Start from `b/a = 1.5` for all precisions (or 1.9 for float if you want the
per-precision value) and re-derive it properly with a dedicated calibration
sweep over `(n, m, M)` rather than reusing the 288 pairs, which were not
designed to separate the two terms. The fit script used for the estimate above
is worth rewriting from scratch; it was throwaway.

**The planner's own `pcost` is likely mis-weighted for the same reason.** That
affects which solver `NFFT_ESTIMATE` picks, which is a separate concern from
this issue — raise it as its own issue rather than changing it here.

## What is and is not achievable

**Achievable:** beating legacy comfortably when the FFT dominates, and closing
most of the gap when nodes dominate. With cost awareness the search can pay for
a larger `n` to buy a smaller `m`, which is exactly the trade the M-dominated
cases want.

**Not achievable analytically: always at least as good as the oracle.** The
baseline gives legacy the smallest `m` whose *measured* error meets the goal.
The tuner's model is an upper envelope over all measured geometries, so it
over-estimates the error of any particular one and hands out about **one extra
`m`** — median headroom 27x against the oracle's 3x, and `ln(27)/D ≈ 0.75` of
an `m` at `sigma = 2`. Even choosing legacy's exact `n`, the tuner would use
`m+1` and be `(2m+4)/(2m+2)` slower on the convolution — about 1.2× at `m = 4`.

Do not try to close that by shrinking the 1.25 safety margin: `ln(1.25)/D` is
0.05 of an `m`. The gap is worst-case-versus-actual, not margin.

Two honest ways to get a real guarantee, both out of scope here unless asked:

- An opt-in **measured refinement**: after picking `(n, m)`, run the transform
  against a direct NDFT and step `m` down while the goal still holds. That is
  what the oracle does; it costs `O(N*M)` and belongs behind a flag.
- A **per-geometry** rather than worst-case error model, which is a different
  and much harder fit.

Note also that against the legacy API *as it actually ships* — no cut-off
estimator at all — the tuner wins by default. The oracle baseline is a
deliberately unfair, and therefore informative, opponent.

## Acceptance

Re-run the head-to-head with the new tuner:

```sh
sh .scratch/sigma-m-study/run_compare.sh "$PWD" /tmp/cmp 50
cd .scratch/sigma-m-study && python3 compare_report.py \
    ../../docs/tuning-vs-legacy.md /tmp/cmp/compare-*.csv
```

`compare.c` calls `tune_plan`; update that call for the new signature. Targets:

- Accuracy unchanged: **288/288 goals met**, no regressions in the `tune` suite
  in any of the three precision trees.
- **No shape has a median speedup below 0.95×** (today: 0.77× M-dominated).
- **Overall median speedup ≥ 1.05×** (today: 0.97×).
- N-dominated stays at or above today's 1.31×.
- Worst single case no worse than 0.90× (today: 0.60×).

If the M-dominated median lands between 0.95× and 1.0×, that is the expected
residual from the extra `m` described above — record it in
`docs/tuning-vs-legacy.md` rather than chasing it.

Also add a test asserting the cost-ordering property directly, so the policy is
pinned independently of timing noise: for a fixed `(N, goal)`, the `n` returned
must be non-decreasing in `M`.

## Resolution

`X(tune_plan)` now takes `M` and walks every even 5-smooth `n` with sigma in
[5/4, 4], taking the smallest sufficient cut-off at each and keeping the pair
that minimises `n*log2(n) + w*M*(2m+2)`.

**The weight is `w = 4/5`, an operation count, not a fitted constant.** The
issue asked for a calibration sweep and a measured weight. The measurement was
done — `costfit.c`, `costfit.py`, 27600 timings over an `(N, n, m, M)` grid in
three precisions, data in `data/costfit-*.csv.gz` — but the constant it
produces is a property of the host: the fitted weight ranges 0.93 to 1.7 by
precision and direction, and the excess over the operation count is scattered
grid access, plus the adjoint scatter costing more than the forward gather.
That is cache behaviour. A library constant may use properties of the
floating-point type but not of the machine, so the sweep only checks the
ranking: 4/5 orders 92.6 % of candidate pairs as measured, the best weight
there 92.7 %, the planner's 2/5 90.4 %. The ranking is flat over [0.8, 1.5].

`check_tune_plan_cost` in `tests/tune_ng.c` pins `n` non-decreasing and `m`
non-increasing in `M`, free of timing noise.

### Against the acceptance targets

| target | result |
|---|---|
| accuracy 288/288, no test regressions | met — 288/288, `checkall_ng` clean in all three precisions |
| no shape median below 0.95x | **missed** — M-dominated 0.91x (was 0.77x) |
| overall median >= 1.05x | **missed** — 1.02x (was 0.97x) |
| N-dominated at or above 1.31x | met — 1.35x |
| worst single case no worse than 0.90x | **missed** — 0.68x (was 0.60x) |

The residual is the one this issue predicted, one `m` more than the oracle,
landing slightly below the 0.95x-to-1.0x band it expected. It is not the
weight: at `w = 1.5` on the double tree the M-dominated median moves by less
than 0.01x, the N-dominated median loses 0.10x, and the worst case gets worse.
The worst cases share one shape — the model demands `m+1` where the oracle
measured `m` to be enough, and with many nodes the convolution is the whole
bill. Closing it needs a per-geometry error model or an opt-in measured
refinement, both out of scope here.


## Follow-up: the geometry-aware model and the measured refinement

Recorded in [03](03-geometry-aware-error-model.md). The short version: the
cost policy above was not what held the M-dominated shape back. Two other
things were, and one of them was a defect.
