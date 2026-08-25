# 05 — The error model is fitted on a different input distribution than it is used on

Status: ready-for-human

## What was found

`gsweep.c` drew both input vectors with real and imaginary parts uniform on
`[-1/2, 1/2)`. `Y(vrand_unit_complex)` — `kernel/util/rand.c:44`, used by the
library's own accuracy tests — draws them on
`[0, 1)`. The two are not interchangeable for the forward measure.

Measured over 200 fresh draws per cell, double precision:

| geometry | forward, `[0,1)` | forward, centred | ratio |
|---|---|---|---|
| N=128, sigma=2, m=1 | 9.81e-03 | 3.77e-03 | 2.60 |
| N=128, sigma=2, m=3 | 1.14e-06 | 6.19e-07 | 1.83 |
| N=100, sigma≈1.28, m=4 | 1.61e-06 | 1.48e-06 | 1.09 |
| N=256, sigma=2, m=5 | 9.74e-11 | 5.52e-11 | 1.77 |

The adjoint measure is within 10 % either way.

Why the forward measure moves and the adjoint does not: uncentred data has a
coherent component, so the forward transform peaks where the phases line up
and the error peaks with it, while the `l1` norm it is divided by does not
concentrate. The adjoint's input sits at the nodes and its norm runs over the
same index, so a constant offset scales both halves of the ratio alike.

## Why it matters

The fitted constants are raised until they dominate every measured point, and
that is the whole basis for the promise that the returned cut-off reaches the
goal. Fitted on the centred draw, the forward envelope is optimistic by up to
2.6x for the uncentred one — the very distribution the tests present.

It surfaced through the dyadic per-band fit (issue 02), which is tighter than
the shipped single-band envelope and so has less accidental slack to absorb
the gap. At N=128, M=256, forward, goal 1e-2 the ladder returned rung 1 with
m=1, predicted 9.97e-3, measured 1.07e-2.

A single envelope over sigma in [5/4, 4] carries the same exposure and hides
it: it happens to leave 2x or more of slack at the geometries the tests
exercise. That is luck, not design, and the per-band fit removed the luck.

## Done for this branch

`gsweep.c` now draws on `[0, 1)`, matching the library. The dyadic per-band
constants (issue 02) are fitted on that data.

## Not done — the real question

Neither draw is a bound. The measure

    max_j |f_j - s_j| / ||fhat||_1

has a finite supremum over all inputs — it is the largest absolute entry of
the error matrix `A - A~`, which does not depend on the input at all — and
random data approaches it from below whatever the distribution is. A tuner
that promises an accuracy for the caller's data, which is arbitrary, should
be calibrated against that supremum rather than against a draw.

It is computable: `N` transforms of unit vectors give the error matrix
exactly, which is affordable at the bandwidths the sweep uses. The cost is
that the envelope becomes markedly more conservative, and every head-to-head
against a legacy cut-off found by *measuring* one random draw would then
lose by construction. So it is a change of what the tuner promises, not a
bug fix, and it needs a decision before it is made.

Two follow-ups, neither taken here:

1. Decide whether `X(tune)` and friends promise "meets the goal for random
   input of this shape" or "meets the goal for any input", and say so in
   `nfft3.h`. Today the header says neither.
2. If the answer is "any input", refit against the error matrix supremum and
   re-run the head-to-head with a matching oracle.

## Applies to any future fit

Any model fitted against this measure must use the same draw as the code it
will be judged by, or state which draw it assumes. A fit on centred data and a
test on uncentred data disagree by up to 2.6x in the forward direction.
