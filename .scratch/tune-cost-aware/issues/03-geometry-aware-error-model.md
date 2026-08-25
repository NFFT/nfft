# 03 — A geometry-aware error model, and an opt-in measured refinement

Status: done

Followed from issue 02, whose acceptance targets the cost policy alone could
not reach. The question asked was how close an error model could get the
cut-off to the one a measured search finds, and whether the M-dominated shape
could be brought to parity with the legacy geometry.

## What was wrong

**The shipped model had no geometry term at all**, because `sweep.c` ties
`M = 2N` and cannot separate `N` from `M`. That is not merely loose. Both
error measures are relative — the forward one divides by an `l1` norm over
`N` terms and takes a max over `M` nodes, the adjoint the other way round —
so the error genuinely moves with the geometry, and an envelope fitted along
one line misses off it. Measured on the crossed grid: **185 of 15226
geometries exceed the shipped envelope**, nearly all adjoint with `M < 2N`,
where the true error is larger than anything the old fit saw.

The over-prediction was also strongly systematic in `N` — 3.0x at `N = 64`,
5.9x at 256, 9.2x at 1024 — which is exactly the `N^-1/2` the measure
implies, and cost a median 0.44 extra cut-off.

## What was done

**`gsweep.c`** crosses `N` in {32 … 1024} with `M/N` in {1/4, 1, 2, 8}, five
draws, `m` to 32, all three precisions. The direct NDFT reference depends on
the nodes and data but not on `sigma` or `m`, so it is hoisted out of the
`sigma`/`m` plane; that is what makes the grid affordable.

**`gfit.py`** fits both branches with `N` and `M` powers, keeps the same
envelope discipline, and validates against the shipped constants on the same
data. Held out at `N` of 2048 and 4096 — four times the top of the fitted
range — the exponents extrapolate:

| | misses | mean extra m | picks the measured m |
|---|---|---|---|
| shipped | 1 | +0.74 fwd, +0.47 adj | 26 %, 53 % |
| geometry-aware | 0 | +0.26 fwd, +0.24 adj | 74 %, 76 % |

**`X(tune_refine)`** is the opt-in measured half: it runs the transform
against the planner's own direct NDFT on the caller's nodes and steps the
cut-off down while the goal holds. `X(tune)` and `X(tune_sigma)` now take `M`
too, since the model needs it.

**A tie-break** prefers the grid with more factors of two between two that
need the same cut-off and cost within 10 % of each other.

## What it is worth

Same 288 configurations, same host:

| | cost-blind | + tie-break | + geometry | + refinement |
|---|---|---|---|---|
| overall median | 0.97x | 1.04x | 1.02x | 1.04x |
| N-dominated | 1.31x | 1.33x | 1.28x | 1.28x |
| balanced | 0.96x | 1.03x | 1.03x | 1.05x |
| M-dominated | 0.77x | 0.92x | 0.94x | 0.95x |
| goal met | 288/288 | 288/288 | 288/288 | 288/288 |

## Two limits worth writing down

**The cut-off is quantised.** Dropping one multiplies the error by `exp(D)`,
30 to 90 across the band. A cut-off is only ever wasted when the headroom
exceeds that, so the model is already exact in about four cases in five and
no refinement can improve those. It also means a safety margin is cheap: it
raises the bar from `exp(D)` to `2*exp(D)`.

**A single measurement is not safe.** At fixed nodes the error still spans a
median 1.55x and up to 6x across input draws. A refinement shaving on one
probe missed the goal in 6 of 288 head-to-head cases. Taking the worst of
eight probes with a factor of two in hand does not: 0 of 3840 fresh draws on
96 refined geometries exceeded the goal.

## What is left, and where it belongs

The M-dominated residual is no longer the cut-off. Split by cause:

| M-dominated subset | cases | median |
|---|---|---|
| oracle's `m` already, `n` a power of two | 15 | 1.00x |
| oracle's `m` already, other `n` | 18 | 0.94x |
| more `m`, at a different `n` | 57 | 0.92x |

Where the tuner already has the oracle's cut-off and lands on a power of two
it is exactly on par. The rest is FFT size quality and the relative wall-clock
price of an FFT against a node convolution — `n = 432` takes 1.7x the time of
`n = 512` in float though it is the smaller grid, and no operation count
predicts that ordering, since counted in arithmetic a radix-3 or radix-5 stage
comes out cheaper per point.

That is machine knowledge, and it must not go into the analytic model: a
constant fitted to one host would travel to every other. A wider power-of-two
preference was tried and does not pay — at a 1.25 window every shape's median
falls. **Choosing among candidate pairs on measured time is measured
planning's job**, and handing the planner a shortlist of pairs to time is the
natural next step. Raise it as its own issue.
