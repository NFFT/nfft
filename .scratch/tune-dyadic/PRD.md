# Dyadic NFFT parameter tuning — design

Written 2026-08-25. Branch `feature/tune-m`, on top of `ea6912c5`.

A second tuner for the 1-D Kaiser-Bessel NFFT, beside the existing
`X(tune_plan)` and not replacing it. It offers three oversampled grid sizes
instead of every even 5-smooth size in a band, and all three are powers of
two.

## 1. The idea

`next_power_of_2(2N) = 2 * next_power_of_2(N)` for every `N`, so the legacy
grid size is one member of a dyadic ladder. Write

    2^k = next_power_of_2(N)        t = 2^k / N,   t in (1, 2]

and the ladder is

| j | n | sigma | note |
|---|---|---|---|
| 0 | `2^k` | `t` | legal only when `t >= 5/4` |
| 1 | `2^(k+1)` | `2t` | the legacy choice |
| 2 | `2^(k+2)` | `4t` | |

One continuous parameter `t` and one discrete index `j` describe the whole
search. `j = 0` halves the FFT work and pays with a larger cut-off; `j = 2`
doubles it to buy a smaller one.

`t >= 5/4` is the tuner's existing hard oversampling floor (`4*n >= 5*N`)
written for `j = 0`. It is the only threshold: everything else is decided by
the accuracy floor and by cost.

## 2. Why this is worth building

The existing `X(tune_plan)` walks every even 5-smooth `n` with `sigma` in
`[5/4, 4]`. Against the legacy grid with an oracle cut-off it reaches an
overall median of 1.04x but only 0.95x when nodes dominate
(`docs/tuning-vs-legacy.md`). That shortfall splits as:

| cause | cases | median |
|---|---|---|
| cut-off already the oracle's, `n` a power of two | 15 | 0.99x |
| cut-off already the oracle's, `n` not a power of two | 18 | 0.94x |
| cut-off above the oracle's | 57 | 0.92x |

The dyadic ladder attacks the middle row structurally — every candidate is a
power of two, so the radix-2 penalty that makes `n = 480` lose to `n = 512`
cannot arise, and the `tune_twos` tie-break is not needed. It attacks the
bottom row indirectly but decisively: when nodes dominate, the ladder's
cheapest rung *is* the legacy grid, and there the median case already gets the
oracle's cut-off, so the median outcome is an exact tie rather than a loss.

### Measured gate

`.scratch/sigma-m-study/dgate.py` replays the `j=0` against `j=1` choice on
the existing `gsweep` measurements, using the sigma pairs in ratio 2 that the
sweep already contains (`t` = 1.25, 1.5, 1.75). Legacy is given the oracle
cut-off; the dyadic side is given the model's and then checked against the
measurement.

6464 cases, **0 accuracy misses**, median speedups:

| population | cases | median | p10 | N-dom | balanced | M-dom |
|---|---|---|---|---|---|---|
| all | 6464 | 1.253x | 1.000x | 1.889x | 1.409x | 1.001x |
| `goal/eps >= 1e3` | 6048 | 1.294x | 1.000x | 1.907x | 1.426x | 1.020x |
| `goal/eps < 1e3` | 416 | 1.000x | 0.638x | 1.000x | 1.000x | 1.000x |

Every shape's median clears 1.00x on both populations. The split matters for
the tail, not the median; section 8 says why.

The replay is a cost-model comparison, not wall time, and it reproduces the
sigma ratio rather than literal power-of-two sizes — so FFT size quality is
absent from it, and in the real method that can only help.

## 3. The model

Same two branches and the same exact rates as `X(tune_plan)`. With
`u = 1 - 1/sigma`:

    b = 2*pi*(1 - 1/(2*sigma))          Kaiser-Bessel shape parameter
    D = 2*pi*sqrt(u)                    truncation decays as exp(-D*m)
    A = b - D                           roundoff amplifies as exp(A*m)

    T = a * u^p * m^r * exp(-D*m)      * N^tn * M^tm
    U = c * eps * u^q * exp(alpha*A*m) * N^un * M^um
    E = T + U

The sweep behind the fit draws its input with real and imaginary parts on
`[0, 1)`, matching `Y(vrand_unit_complex)`. It drew them centred on zero
before this work, which measures a forward error up to 2.6x smaller and made
the forward envelope optimistic for the data the tests actually present --
issue 05 has the numbers and the consequences.

The change to the model itself is that the prefactors are fitted **once per
rung**, on the sigma band that rung occupies:

| j | fitted over sigma |
|---|---|
| 0 | `[5/4, 2]` |
| 1 | `[2, 4]` |
| 2 | `[4, 8]` |

The shipped envelope is raised until it dominates every measurement anywhere
in `[5/4, 4]`. A per-band envelope only has to dominate its own band, so it
can only be tighter. Measured on the existing sweep:

| band | direction | shipped: misses / extra m / exact | per-band |
|---|---|---|---|
| `[1.25, 2]` | forward | **32** / +0.32 / 68 % | 0 / +0.43 / 57 % |
| `[1.25, 2]` | adjoint | 0 / +0.27 / 73 % | 0 / +0.31 / 69 % |
| `[2, 4]` | forward | **64** / +0.17 / 82 % | 0 / +0.27 / 73 % |
| `[2, 4]` | adjoint | 0 / +0.20 / 80 % | 0 / +0.24 / 76 % |
| `[4, 8]` | forward | **97** / +0.08 / 90 % | 0 / +0.17 / 83 % |
| `[4, 8]` | adjoint | 0 / +0.16 / 84 % | 0 / +0.18 / 82 % |

The per-band envelopes look dearer in cut-off than the shipped one. That is
not a like-for-like comparison: the shipped forward constants **miss the goal
193 times** on this data, so they are not a valid tuner on it. The per-band
fit is the first model here that holds, and holding is what the extra cut-off
buys. The adjoint, which the input distribution barely moves, is within 0.04
of an m either way.

On the widest band `alpha` is pinned at 1 rather than fitted. It corrects the
exactly derived rate `A = b - D`, and `A` falls from 0.96 at `sigma = 5/4` to
0.056 at 4 and 0.012 at 8, so over `m <= 32` the branch is flat and the least
squares reads noise. Fitting it there produced cut-offs up to **29** above the
measured need; pinning it brings the worst case back to +1.

At run time the coefficient set is selected by `j`, not by testing `sigma`
against a boundary. Band 0 and band 1 overlap at the single point `sigma = 2`,
which arises when `t = 1`; selecting by `j` makes that unambiguous.

## 4. The algorithm

```
tune_plan_dyadic(N, M, adjoint, goal):
    2^k = next_power_of_2(N)
    cap goal at the floor of rung 2          # the widest grid on offer
    for j in 0, 1, 2:
        n = 2^j * 2^k
        skip if n <= N                       # j = 0 when N is a power of two
        skip if 4*n < 5*N                    # the oversampling floor
        skip if m_max(n) < 1
        m = smallest cut-off with E_j(m) <= capped goal, or skip
        cost = n*log2(n) + (4/5)*M*(2m+2)
    return the cheapest rung
```

`m_max(n) = min(32, n/2 - 1)` as today. No tie-break: every candidate is a
power of two, so the ordering the tie-break exists to correct cannot occur.

The feasibility test is the accuracy floor, `min_m E_j(m) <= goal`, reusing
`tune_floor`. Its risk is asymmetric and the design exploits that: a false
"reachable" misses the accuracy goal, a false "unreachable" only drops a rung
and falls back to a wider grid. It is therefore calibrated to be conservative,
gated on zero false-reachable calls across the sweep.

Rung 1 is always legal (`n = 2^(k+1) >= 2N > N`, and `4n >= 5N`), so the
search never comes back empty and the answer is never rated worse than the
legacy grid by the cost model.

## 5. The API

Declared in `include/nfft3.h` inside `NFFT_DEFINE_PLANNER_API`, implemented in
`kernel/nfft/tune.c` beside the existing entry points. 1-D and Kaiser-Bessel
only, analytic, microseconds.

```c
/* Cheapest (n, m) over the dyadic ladder n = 2^j * next_power_of_2(N),
 * j in {0, 1, 2}. Same contract as X(tune_plan).
 * 1 = goal met, 0 = goal was below the reachable floor and the capped goal
 * was met instead, -1 = bad args. */
int X(tune_plan_dyadic)(NFFT_INT N, NFFT_INT M, int adjoint, R goal,
                        NFFT_INT *n, int *m, R *attained);
```

`X(tune_plan)` keeps its current behaviour and constants. `X(tune)`,
`X(tune_sigma)` and `X(tune_refine)` are unchanged; `X(tune_refine)` accepts
any `n`, so it composes with the new entry point without modification.

## 6. Acceptance bar

Against the legacy grid `n = 2*next_power_of_2(N)` with a cut-off found by
searching upward for the smallest **measured** sufficient one:

1. zero accuracy misses,
2. no shape's median speedup below 1.00x,
3. overall median speedup above 1.00x.

## 7. Validation

Three arms on one N set, all through the same `plan_ng` path so only the
parameter pair differs: legacy with the oracle cut-off, the existing
`X(tune_plan)`, and `X(tune_plan_dyadic)`.

The existing 288-case set in `compare.c` **cannot be reused as it stands**.
Its N values are {243, 250, 251, 255, 256, 512}, every one with `t <= 1.06`,
so rung 0 is illegal in all 288 cases and the dyadic tuner would return the
legacy grid every time. The N set must span `t` across (1, 2]: N just above a
power of two (`t` near 2), near the 5/4 threshold, and at `t = 1` as a
control.

## 8. Known limits

- Rung 2 has no measured data yet. The sweep stops at `sigma = 4`, so band
  `(4, 8]` must be measured before its constants exist (issue 01).
- For `t < 5/4` the ladder reduces to rungs 1 and 2, and for most such N the
  answer is the legacy grid. The method's wins live in `t >= 5/4`, about 60 %
  of a dyadic octave.
- **Rung quantisation amplifies the envelope's caution near the accuracy
  floor.** The model errs towards calling a just-reachable goal unreachable.
  `X(tune_plan)` answers that by nudging `n` up a few percent; the ladder can
  only answer by climbing a whole rung, doubling or quadrupling the grid.
  Every one of the 58 replayed cases below 0.80x took rung 2 for this reason,
  and all 58 sit at `goal/eps < 100`. Above `goal/eps = 1e3` nothing falls
  below 0.809x. This is inherent to a ladder and is the price of its other
  properties.
- The remaining losses at ordinary goals are `+1 m` at rung 1: the envelope's
  cost, unchanged by this work.
- 1-D only, Kaiser-Bessel only, as with the existing tuner.
- The model is calibrated against random input of the shape the library's own
  tests present, not against the supremum over all inputs. That supremum
  exists and is finite, and no random draw reaches it. Issue 05 records what
  this costs and what deciding otherwise would take.

## 9. The work

| issue | title | status |
|---|---|---|
| [01](issues/01-dyadic-band-sweep.md) | Measure the `(4, 8]` oversampling band | done |
| [02](issues/02-per-band-constants.md) | Fit per-band constants and emit them | done |
| [03](issues/03-tune-plan-dyadic.md) | `X(tune_plan_dyadic)`, API, build, tests | done |
| [04](issues/04-three-way-validation.md) | Three-way validation on an extended N set | done |
| [05](issues/05-input-distribution-in-the-fit.md) | The fit's input distribution did not match the library's | fixed here, open for `X(tune_plan)` |

## 10. Measured result

`docs/tuning-dyadic.md`, 576 configurations, wall clock, three arms through
the same `plan_ng` path. Median speedup against the legacy grid with an oracle
cut-off:

| population | dyadic | `tune_plan` |
|---|---|---|
| all 576 | 1.01x | 1.16x |
| N-dominated | 1.11x | 1.50x |
| balanced | 1.04x | 1.16x |
| M-dominated | 1.00x | 1.00x |
| **`t >= 5/4`** (288) | **1.43x** | 1.36x |
| `t < 5/4` (288) | 0.98x | 1.04x |
| goal met | **576/576** | 575/576 |

The acceptance bar in section 6 is met: no accuracy miss, no shape median
below 1.00x, overall median above it.

Read past the bar, the result is sharper than the bar is. The ladder is
built for the case where `N` sits more than a quarter above a power of two,
and there it is the fastest of the three -- 1.43x, ahead of `tune_plan`'s
1.36x, on 288 cases. Where `N` sits just below a power of two the ladder has
no rung under the legacy grid, must take the legacy grid itself, and pays the
envelope's cut-off: 0.98x. `tune_plan` is free to pick a 5-smooth size a few
percent above `5N/4` there, and gets 1.04x.

So overall `tune_plan` is faster on this N set, 1.16x against 1.01x, and the
reason is the mix: half the set has no headroom for the ladder. **The two
tuners are complementary rather than ranked**, and which to call is decided by
`t` alone, which is free to compute. Whether the library should choose between
them for the caller is a design question this work does not settle.

The earlier gate replay predicted 1.25x where the wall clock says 1.01x. That
is not a miss: the replay's population was entirely `t >= 5/4`, and on that
population the wall clock says 1.43x, so the replay was conservative. The gap
is the mix, not the model.

`tune_plan`'s single accuracy miss is the input-distribution defect of issue
05, not a new fault.

### What the medians are worth

In 147 of the 576 cases the ladder returns the legacy pair itself, so both
arms run identical parameters and the measured ratio is pure timing noise.
It reads: median 1.000, mean 0.999, p05 0.938, p95 1.075, min 0.804, max
1.467.

The median is exact and the mean unbiased, so the medians above carry signal.
Individual rows do not: a case at 0.80x is inside what identical parameters
produce. Counts of rows below 1.00x, and the minimum columns in
`docs/tuning-dyadic.md`, are noise-dominated.
