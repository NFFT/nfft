# Kaiser-Bessel sigma/m study

How the cut-off `m` returned by `Y(tune)` is derived, and what it is worth.

## The question

For a one-dimensional Kaiser-Bessel NFFT with bandwidth `N` and oversampled
size `n`, pick the smallest window cut-off `m` that reaches a requested
accuracy, in the measure the test suite uses:

    forward   max_j |f_j - s_j| / ||f_hat||_1
    adjoint   max_k |fhat_k - s_k| / ||f||_1

## The model

Two effects oppose each other, both set by `sigma = n/N` alone. With

    b = 2*pi*(1 - 1/(2*sigma))          Kaiser-Bessel shape parameter
    D = 2*pi*sqrt(1 - 1/sigma)          truncation decay rate
    A = b - D                           roundoff amplification rate

the window is truncated to `2m+2` samples, and what is cut off decays like
`exp(-D*m)`. The deconvolution divides by the window's Fourier transform,
smallest at the band edge, which amplifies roundoff like `exp(A*m)`. Both
rates come from the identity

    b^2 - (pi/sigma)^2 = D^2.

The fitted formula keeps those two rates exactly and fits only prefactors:

    T = a * u^p * m^r * exp(-D*m)      * N^tn * M^tm    u = 1 - 1/sigma
    U = c * e * u^q * exp(alpha*A*m)   * N^un * M^um    e = machine epsilon
    E = T + U

The `N` and `M` powers are the measure's own normalisation, not curve
fitting. The forward error is a max over `M` nodes of a sum of `N` terms whose
phases cancel, divided by an `l1` norm over `N` terms, so it falls roughly
like `N^-1/2` and rises slowly with `M`; the adjoint swaps the two. Fitted:
forward `N^-0.58 M^+0.05`, adjoint `M^-0.51 N^+0.04` on the truncation branch.

Leaving them out is not merely loose. `sweep.c` ties `M` to `2N`, so a model
fitted on it cannot see the dependence at all, and away from that line it
**misses the goal**: 185 of 15226 swept geometries, nearly all adjoint with
`M < 2N`, where the true error exceeds the envelope.

`E` falls, reaches a minimum, and rises again. `tune` scans `m` and returns
the first `m` at or below the goal, or the argmin when the goal is out of
reach.

## The accuracy cap

`b^2 - (pi/sigma)^2 = D^2` forces `b > D`, so `A > 0` at every finite sigma
and the minimum always exists. Balancing the two terms puts it at

    E_min ~ eps^gamma,   gamma = D/b = sqrt(1-1/sigma) / (1 - 1/(2*sigma)).

`gamma < 1` for every finite sigma and reaches 1 only as sigma goes to
infinity. Oversampling buys back mantissa bits but never all of them:

| sigma | gamma | float | double | long double |
|-------|-------|-------|--------|-------------|
| 1.25 | 0.745 | 6.9e-06 | 2.2e-12 | 7.4e-26 |
| 1.50 | 0.866 | 1.0e-06 | 2.8e-14 | 6.3e-30 |
| 2.00 | 0.943 | 3.0e-07 | 1.7e-15 | 1.6e-32 |
| 3.00 | 0.980 | 1.7e-07 | 4.6e-16 | 9.2e-34 |
| 4.00 | 0.990 | 1.4e-07 | 3.2e-16 | 4.3e-34 |

Run `cap.py` to reproduce it. The exponent was checked against the sweep by
regressing each geometry's measured minimum against `eps` across the three
precisions; measured `gamma` tracks `D/b` closely (at sigma = 1.25, 0.752
measured against 0.745 predicted).

In exact arithmetic the cap disappears -- `T` alone decays without bound, so
any accuracy is reachable at any `sigma > 1`. The cap is purely a
floating-point effect.

`tune` floors sigma at `5/4` and refuses anything below, in exact integer
arithmetic (`4*n >= 5*N`).

## The experiment

`gsweep.c` measures the error of a Kaiser-Bessel `plan_ng` over a grid of
`sigma` in `[1.25, 4]`, `N` in `{32 … 1024}` crossed with `M/N` in
`{1/4, 1, 2, 8}`, `m` in `[1, 32]`, five random trials per cell, taking the
worst. The reference is a long-double direct NDFT; since it depends on the
nodes and data but not on `sigma` or `m`, it is computed once per
`(N, M, trial)` and reused across the whole `sigma`/`m` plane, which is what
makes the grid affordable. One binary per precision.

    sh run_gsweep.sh <worktree-root> <out-dir> 32 5

Optional third and fourth arguments override the `N` list and the `M` factors,
which is how the held-out validation at `N` of 2048 and 4096 was run.

`sweep.c` and `run_sweeps.sh` are the older, `M = 2N` version, kept because
the shipped constants before the geometry terms came from it.

## The fit

`gfit.py` splits each geometry's series at its measured minimum, regresses the
truncation branch for `a, p, r, tn, tm` and the roundoff branch for
`c, alpha, q, un, um`, then raises `a` and `c` to the smallest pair that
dominates *every* measurement, choosing among the dominating pairs the one
that keeps the envelope tightest. A 1.25x margin is applied on top.

    python3 gfit.py <out-dir>/gsweep-*.csv [--holdout <csv>...]

The constants it prints go into `kernel/nfft/tune.c`. It also validates the
shipped constants on the same data for a like-for-like baseline, and anything
after `--holdout` is excluded from the fit and validated separately -- the
geometry exponents are the part most at risk of not extrapolating, so the
holdout should sit outside the fitted box.

Held out at `N` of 2048 and 4096, four times the top of the fitted range:

| | misses | mean extra m | picks the measured m |
|---|---|---|---|
| the model this replaced | 1 | +0.74 forward, +0.47 adjoint | 26 %, 53 % |
| geometry-aware | 0 | +0.26 forward, +0.24 adjoint | 74 %, 76 % |

`fit.py` is the older fit, for the `M = 2N` sweep and the model without the
geometry terms.

## What the fit is worth

`fit.py` closes with the only check that matters for a tuner: over a ladder of
goals crossed with every geometry, it compares the `m` the model picks against
the smallest `m` the *measurement* shows to be sufficient. Reported per run:

- cases where the returned `m` misses the goal in measurement (must be 0),
- how much extra `m` the safety envelope costs,
- how pessimistic the reported attainable floor is.

Results are printed by each run; `constants.py` reports them alongside the
constants it emits.

## The three entry points

`Y(tune)` answers "what m at this geometry", `Y(tune_sigma)` answers "what
oversampling for this accuracy", and `Y(tune_plan)` answers "what pair for
this problem", which needs the node count M as well:

1. cap the goal at the best accuracy the window can reach -- the floor at the
   top of the band;
2. walk every even 5-smooth n with sigma in [5/4, 4] -- only tens of them, and
   the FFT runs at full speed on such a size (`nsweep.c` measures why: a size
   with a prime factor above 13 costs about 3x one without);
3. at each, take the smallest m reaching the capped goal, and keep the pair
   that costs least.

The cost is `n*log2(n) + (4/5)*M*(2m+2)`: the FFT is O(n log n) and the node
convolution O(M*(2m+2)), so more nodes argue for a larger grid and a smaller
cut-off. The weight 4/5 is an operation count -- about `5*n*log2(n)` real
operations for the FFT, four per window sample -- not a measurement, so the
library does not carry a constant fitted to one machine. `costfit.c` checks
the ranking it produces; see below.

The legacy size `2*next_power_of_2(N)` is a power of two in [2N, 4N), so it is
always one of the candidates and the answer is never rated worse than it.

Between two grids that need the same cut-off and cost within 10 % of each
other, the one with more factors of two wins. That is a tie-break, not a term
in the cost: FFT implementations optimise radix-2 hardest, so `n = 480` loses
to `n = 512` despite being smaller, and in float `n = 432` takes 1.7x the time
of `n = 512`. No operation count reproduces that ordering -- counted in
arithmetic, a radix-3 or radix-5 stage comes out cheaper per point. Widening
the window to chase the larger float spread does not pay: at 1.25 every
shape's median speed falls.

## The measured refinement

`Y(tune_refine)` is the opt-in half. It measures the transform against the
planner's own direct NDFT on the caller's nodes and steps the cut-off down
while the goal still holds, which is the only way to beat an envelope fitted
across geometries -- at a fixed geometry the error still varies by a factor of
a few from one set of nodes to another, and no formula can see that.

Two things bound what it can return:

- **Quantisation.** Dropping one cut-off multiplies the error by `exp(D)`,
  which is 30 to 90 across the band. A cut-off is only ever wasted when the
  headroom exceeds that, so the model is already exact in about four cases in
  five, and the refinement fires in the rest.
- **The input draw.** At fixed nodes the error still spans a median 1.55x and
  up to 6x across random inputs. Shaving on a single probe measurably misses
  the goal on the next vector -- 6 of 288 head-to-head cases. The refinement
  therefore takes the worst of eight probes and keeps a factor of two in hand;
  over 3840 fresh draws on 96 refined geometries, none exceeded the goal.

Step 2 uses `tune_next_smooth`, which enumerates every 3^b*5^c up to the bound
and rounds each up by powers of two -- O(log^3 n), not an upward scan. Scanning
would be near-linear: 5-smooth numbers thin out multiplicatively, so gaps near
k grow like k/log^3 k.

## The cost weight

`costfit.c` times one transform over a grid of (N, n, m, M) that moves the two
terms independently -- the head-to-head data cannot, since n and m move
together there. `costfit.py` reports both a least-squares fit of the weight and
the weight that ranks the most candidate pairs correctly.

    sh run_costfit.sh <worktree-root> <out-dir> 0.02
    cd .scratch/sigma-m-study && python3 costfit.py <out-dir>/costfit-*.csv

27600 timings on one aarch64 host, per precision and direction, put the fitted
weight between 0.93 and 1.7, and the best-ranking weight between 0.6 and 1.6.
The ranking is flat across that whole range:

| weight | pairs ordered as measured |
|---|---|
| 0.4 (the planner's `pcost`) | 90.4 % |
| **0.8 (the operation count, used)** | **92.6 %** |
| 1.0 | 92.7 % |
| 1.5 | 92.2 % |
| 3.0 | 89.8 % |

Measured weights above the operation count are the scattered grid access
costing more than its operation count says, and the adjoint scatter costing
more than the forward gather. Both are properties of one machine's memory
system, so neither is baked in.

## The n sweep

`nsweep.c` sweeps every integer n from 2N to 4N for six representative N,
measuring accuracy and mean transform time. `nanalyze.py` and `nsmooth.py`
summarise it; `nreport.py` emits the JSON behind the HTML report. It is what
justifies the 5-smooth rounding in `tune_plan`, and it shows the legacy
`2 * next_power_of_2(N)` rule overshooting when N sits just above a power of
two.

## Files

| File | What it is |
|---|---|
| `gsweep.c`, `run_gsweep.sh` | the sigma/m/N/M accuracy sweep, one binary per precision |
| `gfit.py` | fits the geometry-aware model, validates it, prints the constants |
| `sweep.c`, `run_sweeps.sh`, `fit.py` | the older M = 2N sweep and its fit |
| `selectform.py` | picks the model form by validation, not by in-sample fit |
| `constants.py` | prints the constants that go into `kernel/nfft/tune.c` |
| `cap.py`, `floors.py` | the attainable-accuracy cap, predicted and measured |
| `nsweep.c`, `run_nsweep.sh` | the fine-grained n sweep, accuracy and time |
| `nanalyze.py`, `nsmooth.py` | summarise it; compare the ways of choosing n |
| `nreport.py` | emits the JSON behind `nsweep-report.html` |
| `costfit.c`, `run_costfit.sh` | the (N, n, m, M) timing grid behind the cost weight |
| `gsweep-report` data | `data/gsweep-*.csv.gz`, `data/big-*.csv.gz` (the held-out N) |
| `costfit.py` | fits the weight and scores how well it ranks |
| `compare.c`, `run_compare.sh` | the head-to-head against the legacy parameter choice |
| `compare_report.py` | renders it into `docs/tuning-vs-legacy.md` |
| `plantable.c`, `table.c` | print what the tuning entry points return |
| `run_tests.sh` | build and run `checkall_ng` in one precision tree |
| `data/*.csv.gz` | the measurements behind the constants, the weight and the comparison |
| `nsweep-report.html` | the rendered n-sweep summary |

Regenerating everything, from the worktree root:

```sh
sh .scratch/sigma-m-study/run_sweeps.sh "$PWD" /tmp/sweeps 40
cd .scratch/sigma-m-study && python3 constants.py /tmp/sweeps/sweep-*.csv
sh .scratch/sigma-m-study/run_nsweep.sh "$PWD" /tmp/nsw 6 100
sh .scratch/sigma-m-study/run_costfit.sh "$PWD" /tmp/cost 0.02
cd .scratch/sigma-m-study && python3 costfit.py /tmp/cost/costfit-*.csv
sh .scratch/sigma-m-study/run_gsweep.sh "$PWD" /tmp/gs 32 5
cd .scratch/sigma-m-study && python3 gfit.py /tmp/gs/gsweep-*.csv
```
