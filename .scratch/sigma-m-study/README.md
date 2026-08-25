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

    T(m,s)   = a * u^p * m^r * exp(-D*m)                u = 1 - 1/sigma
    U(m,s,e) = c * e * u^q * m^w * exp(alpha*A*m)       e = machine epsilon
    E(m,s,e) = T + U

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

`sweep.c` measures the error of a Kaiser-Bessel `plan_ng` over a grid of
`sigma` in `[1.25, 4]`, `N` in `{64, 256, 1024}`, `m` in `[1, 40]`, with
`M = 2N` nodes and three random trials per cell, taking the worst trial. The
reference is a long-double direct NDFT computed in the driver, so the oracle
is well below the measured errors. One binary per precision.

    sh run_sweeps.sh <worktree-root> <out-dir> 40

## The fit

`fit.py` splits each geometry's series at its measured minimum, regresses the
truncation branch for `a, p, r` and the roundoff branch for `c, alpha, q, w`,
then raises `a` and `c` to the smallest pair that dominates *every*
measurement, choosing among the dominating pairs the one that keeps the
envelope tightest. A 1.25x margin is applied on top.

    python3 fit.py <out-dir>/sweep-d.csv <out-dir>/sweep-f.csv <out-dir>/sweep-l.csv

The constants it prints go into `kernel/nfft/tune.c`.

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
| `sweep.c`, `run_sweeps.sh` | the sigma/m accuracy sweep, one binary per precision |
| `fit.py` | fits the model and reports what the fit is worth |
| `selectform.py` | picks the model form by validation, not by in-sample fit |
| `constants.py` | prints the constants that go into `kernel/nfft/tune.c` |
| `cap.py`, `floors.py` | the attainable-accuracy cap, predicted and measured |
| `nsweep.c`, `run_nsweep.sh` | the fine-grained n sweep, accuracy and time |
| `nanalyze.py`, `nsmooth.py` | summarise it; compare the ways of choosing n |
| `nreport.py` | emits the JSON behind `nsweep-report.html` |
| `costfit.c`, `run_costfit.sh` | the (N, n, m, M) timing grid behind the cost weight |
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
```
