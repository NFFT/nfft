# Tuned parameters against the legacy choice

Does `nfft_tune_plan` pick better `(n, m)` than the legacy default?

**Legacy** is `n = 2 * next_power_of_2(N)` with the cut-off found by
searching upward for the smallest `m` whose *measured* error meets the
goal. The legacy API has no cut-off estimator, so this hands it the best
`m` it could ever have used -- an oracle it does not actually possess.
**New** is whatever `nfft_tune_plan` returns, unaided.

Both run through the same `plan_ng` path, so the only thing that differs
is the parameter pair. Error is the standard measure against a
long-double direct NDFT; time is the mean forward or adjoint execution
over 50 repetitions, planning excluded.

## Headline

- 288 configurations: 6 bandwidths x 3 shapes x 2-3 goals x 2 directions x 3 precisions.
- Goal met: **new 288/288**, legacy-with-oracle-m 288/288.
- New is faster in **155/288** (54 %) of them.
- Median speedup **1.02x**.

## Verdict

**Accuracy: yes, consistently.** The tuner meets the goal in every one
of the 288 configurations, unaided. The legacy geometry also meets it in
every case, but only because it was handed the oracle cut-off; the legacy
API ships no way to find that `m`, so in practice the choice is between a
tuner that lands the accuracy and a guess that might not.

**Speed: it depends on how many nodes there are.**
With few nodes the tuner wins clearly (1.28x median, 85/96 cases); with many nodes it still trails (0.94x median, 19/96 wins).

`tune_plan` takes the node count and picks the pair its cost model
rates cheapest, `n*log2(n) + (4/5)*M*(2m+2)`, over every even 5-smooth
`n` with sigma in [5/4, 4]. The legacy size is a power of two in that
range, so it is always among the candidates.

The error model is geometry-aware: both measures are relative, so the
forward error falls like `N^-0.58` and rises slowly with `M`, and the
adjoint the other way round. Without those terms the envelope is not
just loose but wrong away from `M = 2N` -- the model this replaced
misses the goal outright in 185 of 15226 swept geometries.

### Where the remaining M-dominated gap sits

| M-dominated subset | cases | median |
|---|---|---|
| oracle's `m` already, `n` a power of two | 15 | 0.99x |
| oracle's `m` already, other `n` | 18 | 0.94x |
| more `m`, at a different `n` | 57 | 0.92x |

Both residuals are the same thing: an analytic model cannot know this
machine's relative speed of an FFT against a node convolution. Where
the tuner already has the oracle's cut-off and lands on a power of
two, it is exactly on par. Where it lands on another 5-smooth size it
pays for that size's codelets, which no operation count predicts --
`n = 432` measures 1.7x the time of `n = 512` in float though it is
the smaller grid. Where it takes an extra cut-off, the cost model
preferred a smaller grid on an operation count that undervalues the
convolution in wall time.

Neither is fixable in the model without fitting it to one machine.
A wider power-of-two preference was tried and costs more than it
returns: at a 1.25 window every shape's median falls. Choosing among
candidate pairs on measured time is what measured planning is for.

## How it got here

Four steps, all on the same 288 configurations and the same host. The
first column is the original policy: minimise `n`, then take whatever
cut-off that grid needs, with no `M` anywhere.

| | cost-blind | + tie-break | + geometry model | + refinement |
|---|---|---|---|---|
| overall median | 0.97x | 1.04x | 1.02x | 1.04x |
| N-dominated | 1.31x | 1.33x | 1.28x | 1.28x |
| balanced | 0.96x | 1.03x | 1.03x | 1.05x |
| M-dominated | 0.77x | 0.92x | 0.94x | 0.95x |
| goal met | 288/288 | 288/288 | 288/288 | 288/288 |

The tie-break prefers a power-of-two-richer grid when two need the
same cut-off and cost within 10 % of each other. The geometry model
adds the `N` and `M` powers. The refinement is opt-in and measured.

## With the measured refinement

`nfft_tune_refine` measures the transform against a direct NDFT on
the caller's own nodes and steps the cut-off down while the goal
still holds. It is opt-in and costs one O(N*M) probe, so the driver
calls it only where the cost model says the node convolution is at
least 30 % of the bill -- below that a needless cut-off is too
cheap to be worth removing.

| | model only | with refinement |
|---|---|---|
| overall median | 1.02x | 1.04x |
| N-dominated | 1.28x | 1.28x |
| balanced | 1.03x | 1.05x |
| M-dominated | 0.94x | 0.95x |
| goal met | 288/288 | 288/288 |

It is judged on eight probe vectors and keeps a factor of two in
hand, because at fixed nodes the error still spans a median 1.55x
and up to 6x across input draws -- shaving to a single probe does
miss the goal on the next vector. Over 3840 fresh draws on 96
refined geometries, none exceeded the goal.

## By problem shape

This is where the answer splits. `M` is the node count.

| shape | M | median speedup | new faster | median n_new / n_leg | median m_new - m_leg |
|---|---|---|---|---|---|
| N-dominated | N/4 | 1.28x | 85/96 | 0.62 | +2 |
| balanced | N | 1.03x | 51/96 | 0.75 | +1 |
| M-dominated | 4N | 0.94x | 19/96 | 0.78 | +1 |

## By precision

| precision | median speedup | new faster | goal met (new) | goal met (legacy) |
|---|---|---|---|---|
| float | 0.95x | 30/72 | 72/72 | 72/72 |
| double | 1.06x | 67/108 | 108/108 | 108/108 |
| long double | 1.01x | 58/108 | 108/108 | 108/108 |

## By bandwidth

| N | factors | median speedup | median n_new | n_leg | median m_new | median m_leg |
|---|---|---|---|---|---|---|
| 243 | 3·3·3·3·3 | 1.08x | 320 | 512 | 6 | 4 |
| 250 | 2·5·5·5 | 1.03x | 384 | 512 | 6 | 4 |
| 251 | 251 | 1.04x | 384 | 512 | 6 | 4 |
| 255 | 3·5·17 | 0.98x | 384 | 512 | 6 | 4 |
| 256 | 2·2·2·2·2·2·2·2 | 1.01x | 384 | 512 | 6 | 4 |
| 512 | 2·2·2·2·2·2·2·2·2 | 1.01x | 720 | 1024 | 6 | 4 |

## Accuracy

Both sides are asked for the same goal, so the question is not who is
more accurate but whether either misses. Cases where the goal was not
met:

- **New: none.**
- **Legacy (with oracle m): none.**

Median headroom below the goal: new 11.8x, legacy 3.1x.
Both overshoot the goal rather than miss it; the new side does so by
more, which is the cost of an upper-bound model against a measured search.

## Full results

| precision | N | M | shape | goal | dir | n_new | m_new | err_new | t_new (µs) | n_leg | m_leg | err_leg | t_leg (µs) | speedup |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| float | 243 | 60 | N-dominated | 1e-02 | adjoint | 320 | 2 | 2.13e-03 | 0.46 | 512 | 1 | 9.06e-03 | 0.62 | 1.36x |
| float | 243 | 60 | N-dominated | 1e-02 | forward | 320 | 2 | 2.87e-04 | 0.40 | 512 | 1 | 2.36e-03 | 0.57 | 1.42x |
| float | 243 | 60 | N-dominated | 1e-04 | adjoint | 400 | 3 | 8.74e-06 | 0.57 | 512 | 3 | 1.59e-06 | 0.66 | 1.17x |
| float | 243 | 60 | N-dominated | 1e-04 | forward | 320 | 3 | 1.15e-05 | 0.40 | 512 | 2 | 3.22e-05 | 0.66 | 1.63x |
| float | 243 | 243 | balanced | 1e-02 | adjoint | 320 | 2 | 1.38e-03 | 0.84 | 512 | 1 | 6.88e-03 | 0.91 | 1.08x |
| float | 243 | 243 | balanced | 1e-02 | forward | 320 | 2 | 4.07e-04 | 0.67 | 512 | 1 | 2.96e-03 | 0.75 | 1.12x |
| float | 243 | 243 | balanced | 1e-04 | adjoint | 384 | 3 | 6.56e-06 | 1.23 | 512 | 3 | 1.29e-06 | 1.15 | 0.93x |
| float | 243 | 243 | balanced | 1e-04 | forward | 320 | 3 | 1.37e-05 | 0.69 | 512 | 2 | 3.12e-05 | 0.84 | 1.22x |
| float | 243 | 972 | M-dominated | 1e-02 | adjoint | 320 | 2 | 3.98e-04 | 2.17 | 512 | 1 | 3.16e-03 | 2.05 | 0.94x |
| float | 243 | 972 | M-dominated | 1e-02 | forward | 320 | 2 | 3.39e-04 | 1.67 | 512 | 1 | 3.36e-03 | 1.64 | 0.98x |
| float | 243 | 972 | M-dominated | 1e-04 | adjoint | 384 | 3 | 2.55e-06 | 3.09 | 512 | 2 | 4.24e-05 | 2.36 | 0.76x |
| float | 243 | 972 | M-dominated | 1e-04 | forward | 320 | 3 | 1.69e-05 | 1.76 | 512 | 2 | 3.95e-05 | 1.81 | 1.03x |
| float | 250 | 62 | N-dominated | 1e-02 | adjoint | 320 | 2 | 4.46e-03 | 0.45 | 512 | 2 | 2.18e-04 | 0.64 | 1.42x |
| float | 250 | 62 | N-dominated | 1e-02 | forward | 320 | 2 | 3.80e-04 | 0.40 | 512 | 1 | 2.64e-03 | 0.58 | 1.45x |
| float | 250 | 62 | N-dominated | 1e-04 | adjoint | 400 | 4 | 1.78e-06 | 0.77 | 512 | 3 | 2.09e-06 | 0.67 | 0.87x |
| float | 250 | 62 | N-dominated | 1e-04 | forward | 320 | 3 | 1.68e-05 | 0.40 | 512 | 2 | 3.09e-05 | 0.60 | 1.51x |
| float | 250 | 250 | balanced | 1e-02 | adjoint | 320 | 2 | 1.72e-03 | 0.86 | 512 | 1 | 4.72e-03 | 0.93 | 1.08x |
| float | 250 | 250 | balanced | 1e-02 | forward | 320 | 2 | 4.15e-04 | 0.67 | 512 | 1 | 2.90e-03 | 0.78 | 1.15x |
| float | 250 | 250 | balanced | 1e-04 | adjoint | 384 | 3 | 7.37e-06 | 1.26 | 512 | 2 | 7.76e-05 | 1.00 | 0.79x |
| float | 250 | 250 | balanced | 1e-04 | forward | 324 | 3 | 1.67e-05 | 1.07 | 512 | 2 | 3.85e-05 | 0.85 | 0.80x |
| float | 250 | 1000 | M-dominated | 1e-02 | adjoint | 320 | 2 | 9.33e-04 | 2.21 | 512 | 1 | 3.13e-03 | 2.05 | 0.93x |
| float | 250 | 1000 | M-dominated | 1e-02 | forward | 320 | 2 | 4.88e-04 | 1.69 | 512 | 1 | 2.93e-03 | 1.82 | 1.08x |
| float | 250 | 1000 | M-dominated | 1e-04 | adjoint | 384 | 3 | 4.48e-06 | 3.01 | 512 | 2 | 5.26e-05 | 2.38 | 0.79x |
| float | 250 | 1000 | M-dominated | 1e-04 | forward | 384 | 3 | 3.66e-06 | 2.07 | 512 | 2 | 3.83e-05 | 1.93 | 0.93x |
| float | 251 | 62 | N-dominated | 1e-02 | adjoint | 320 | 2 | 2.07e-03 | 0.42 | 512 | 2 | 1.39e-04 | 0.79 | 1.86x |
| float | 251 | 62 | N-dominated | 1e-02 | forward | 320 | 2 | 3.43e-04 | 0.42 | 512 | 1 | 2.46e-03 | 0.61 | 1.44x |
| float | 251 | 62 | N-dominated | 1e-04 | adjoint | 400 | 4 | 1.48e-06 | 0.59 | 512 | 3 | 2.04e-06 | 0.64 | 1.09x |
| float | 251 | 62 | N-dominated | 1e-04 | forward | 324 | 3 | 1.42e-05 | 0.67 | 512 | 2 | 3.10e-05 | 0.62 | 0.94x |
| float | 251 | 251 | balanced | 1e-02 | adjoint | 320 | 2 | 1.56e-03 | 0.86 | 512 | 1 | 4.13e-03 | 0.94 | 1.10x |
| float | 251 | 251 | balanced | 1e-02 | forward | 320 | 2 | 4.69e-04 | 0.70 | 512 | 1 | 2.41e-03 | 0.78 | 1.13x |
| float | 251 | 251 | balanced | 1e-04 | adjoint | 384 | 3 | 7.32e-06 | 1.41 | 512 | 2 | 7.45e-05 | 1.00 | 0.71x |
| float | 251 | 251 | balanced | 1e-04 | forward | 324 | 3 | 1.57e-05 | 0.98 | 512 | 2 | 3.59e-05 | 0.85 | 0.87x |
| float | 251 | 1004 | M-dominated | 1e-02 | adjoint | 320 | 2 | 8.55e-04 | 2.37 | 512 | 1 | 2.78e-03 | 2.16 | 0.91x |
| float | 251 | 1004 | M-dominated | 1e-02 | forward | 320 | 2 | 4.08e-04 | 1.68 | 512 | 1 | 2.51e-03 | 1.60 | 0.95x |
| float | 251 | 1004 | M-dominated | 1e-04 | adjoint | 384 | 3 | 3.76e-06 | 3.33 | 512 | 2 | 5.01e-05 | 2.40 | 0.72x |
| float | 251 | 1004 | M-dominated | 1e-04 | forward | 384 | 3 | 3.62e-06 | 2.16 | 512 | 2 | 3.88e-05 | 1.89 | 0.87x |
| float | 255 | 63 | N-dominated | 1e-02 | adjoint | 324 | 2 | 5.56e-03 | 0.73 | 512 | 2 | 2.17e-04 | 0.65 | 0.89x |
| float | 255 | 63 | N-dominated | 1e-02 | forward | 320 | 2 | 4.80e-04 | 0.43 | 512 | 1 | 2.32e-03 | 0.65 | 1.50x |
| float | 255 | 63 | N-dominated | 1e-04 | adjoint | 432 | 3 | 5.95e-06 | 1.28 | 512 | 3 | 1.56e-06 | 0.69 | 0.53x |
| float | 255 | 63 | N-dominated | 1e-04 | forward | 320 | 4 | 1.31e-06 | 0.64 | 512 | 2 | 3.45e-05 | 0.67 | 1.04x |
| float | 255 | 255 | balanced | 1e-02 | adjoint | 320 | 2 | 2.41e-03 | 0.88 | 512 | 1 | 4.37e-03 | 0.92 | 1.05x |
| float | 255 | 255 | balanced | 1e-02 | forward | 320 | 2 | 4.46e-04 | 1.18 | 512 | 1 | 2.90e-03 | 0.78 | 0.66x |
| float | 255 | 255 | balanced | 1e-04 | adjoint | 384 | 3 | 5.93e-06 | 1.27 | 512 | 2 | 7.58e-05 | 1.12 | 0.89x |
| float | 255 | 255 | balanced | 1e-04 | forward | 384 | 3 | 3.06e-06 | 1.00 | 512 | 2 | 3.62e-05 | 0.85 | 0.85x |
| float | 255 | 1020 | M-dominated | 1e-02 | adjoint | 512 | 1 | 3.23e-03 | 2.10 | 512 | 1 | 3.23e-03 | 2.04 | 0.97x |
| float | 255 | 1020 | M-dominated | 1e-02 | forward | 320 | 2 | 6.12e-04 | 1.78 | 512 | 1 | 3.23e-03 | 1.60 | 0.90x |
| float | 255 | 1020 | M-dominated | 1e-04 | adjoint | 384 | 3 | 4.13e-06 | 3.22 | 512 | 2 | 5.99e-05 | 2.38 | 0.74x |
| float | 255 | 1020 | M-dominated | 1e-04 | forward | 384 | 3 | 3.84e-06 | 2.03 | 512 | 2 | 4.09e-05 | 1.88 | 0.93x |
| float | 256 | 64 | N-dominated | 1e-02 | adjoint | 324 | 2 | 2.13e-03 | 0.72 | 512 | 2 | 2.03e-04 | 0.68 | 0.94x |
| float | 256 | 64 | N-dominated | 1e-02 | forward | 320 | 2 | 5.60e-04 | 0.42 | 512 | 1 | 2.29e-03 | 0.58 | 1.38x |
| float | 256 | 64 | N-dominated | 1e-04 | adjoint | 432 | 3 | 4.06e-06 | 1.15 | 512 | 3 | 2.22e-06 | 0.73 | 0.64x |
| float | 256 | 64 | N-dominated | 1e-04 | forward | 320 | 4 | 1.53e-06 | 0.72 | 512 | 2 | 3.23e-05 | 0.62 | 0.86x |
| float | 256 | 256 | balanced | 1e-02 | adjoint | 320 | 2 | 1.70e-03 | 0.88 | 512 | 1 | 4.25e-03 | 0.90 | 1.02x |
| float | 256 | 256 | balanced | 1e-02 | forward | 320 | 2 | 4.65e-04 | 0.66 | 512 | 1 | 2.41e-03 | 0.77 | 1.16x |
| float | 256 | 256 | balanced | 1e-04 | adjoint | 384 | 3 | 1.55e-05 | 1.27 | 512 | 2 | 9.44e-05 | 0.97 | 0.76x |
| float | 256 | 256 | balanced | 1e-04 | forward | 384 | 3 | 2.84e-06 | 0.99 | 512 | 2 | 3.52e-05 | 0.83 | 0.84x |
| float | 256 | 1024 | M-dominated | 1e-02 | adjoint | 512 | 1 | 2.31e-03 | 2.11 | 512 | 1 | 2.31e-03 | 2.10 | 1.00x |
| float | 256 | 1024 | M-dominated | 1e-02 | forward | 320 | 2 | 5.44e-04 | 1.72 | 512 | 1 | 3.09e-03 | 1.59 | 0.93x |
| float | 256 | 1024 | M-dominated | 1e-04 | adjoint | 384 | 3 | 4.47e-06 | 3.11 | 512 | 2 | 3.92e-05 | 2.42 | 0.78x |
| float | 256 | 1024 | M-dominated | 1e-04 | forward | 384 | 3 | 4.25e-06 | 2.06 | 512 | 2 | 3.87e-05 | 1.88 | 0.91x |
| float | 512 | 128 | N-dominated | 1e-02 | adjoint | 640 | 2 | 3.40e-03 | 0.96 | 1024 | 1 | 8.45e-03 | 1.36 | 1.42x |
| float | 512 | 128 | N-dominated | 1e-02 | forward | 640 | 2 | 3.31e-04 | 0.85 | 1024 | 1 | 1.68e-03 | 1.24 | 1.46x |
| float | 512 | 128 | N-dominated | 1e-04 | adjoint | 864 | 3 | 4.68e-06 | 2.73 | 1024 | 3 | 2.27e-06 | 1.94 | 0.71x |
| float | 512 | 128 | N-dominated | 1e-04 | forward | 640 | 3 | 2.14e-05 | 0.84 | 1024 | 2 | 2.42e-05 | 1.30 | 1.54x |
| float | 512 | 512 | balanced | 1e-02 | adjoint | 640 | 2 | 1.66e-03 | 1.65 | 1024 | 1 | 4.58e-03 | 1.90 | 1.15x |
| float | 512 | 512 | balanced | 1e-02 | forward | 640 | 2 | 4.02e-04 | 1.36 | 1024 | 1 | 2.13e-03 | 1.61 | 1.18x |
| float | 512 | 512 | balanced | 1e-04 | adjoint | 768 | 3 | 1.06e-05 | 2.21 | 1024 | 2 | 6.27e-05 | 2.10 | 0.95x |
| float | 512 | 512 | balanced | 1e-04 | forward | 648 | 3 | 1.64e-05 | 1.97 | 1024 | 2 | 2.87e-05 | 1.80 | 0.92x |
| float | 512 | 2048 | M-dominated | 1e-02 | adjoint | 864 | 1 | 2.73e-03 | 5.00 | 1024 | 1 | 2.20e-03 | 4.28 | 0.86x |
| float | 512 | 2048 | M-dominated | 1e-02 | forward | 640 | 2 | 4.32e-04 | 3.44 | 1024 | 1 | 2.13e-03 | 3.34 | 0.97x |
| float | 512 | 2048 | M-dominated | 1e-04 | adjoint | 768 | 3 | 2.69e-06 | 5.85 | 1024 | 2 | 4.28e-05 | 4.90 | 0.84x |
| float | 512 | 2048 | M-dominated | 1e-04 | forward | 648 | 3 | 2.12e-05 | 4.22 | 1024 | 2 | 3.24e-05 | 3.80 | 0.90x |
| double | 243 | 60 | N-dominated | 1e-04 | adjoint | 320 | 4 | 4.30e-06 | 0.96 | 512 | 3 | 1.15e-06 | 1.28 | 1.32x |
| double | 243 | 60 | N-dominated | 1e-04 | forward | 320 | 3 | 1.13e-05 | 0.74 | 512 | 2 | 3.22e-05 | 1.26 | 1.71x |
| double | 243 | 60 | N-dominated | 1e-08 | adjoint | 320 | 7 | 2.29e-10 | 1.13 | 512 | 5 | 1.24e-10 | 1.51 | 1.33x |
| double | 243 | 60 | N-dominated | 1e-08 | forward | 320 | 6 | 8.63e-10 | 0.82 | 512 | 4 | 2.84e-09 | 1.17 | 1.42x |
| double | 243 | 60 | N-dominated | 1e-11 | adjoint | 384 | 8 | 5.42e-14 | 1.26 | 512 | 6 | 1.30e-12 | 1.42 | 1.13x |
| double | 243 | 60 | N-dominated | 1e-11 | forward | 324 | 8 | 7.45e-13 | 0.95 | 512 | 6 | 4.33e-13 | 1.23 | 1.30x |
| double | 243 | 243 | balanced | 1e-04 | adjoint | 320 | 4 | 2.86e-06 | 1.79 | 512 | 3 | 1.05e-06 | 2.07 | 1.16x |
| double | 243 | 243 | balanced | 1e-04 | forward | 320 | 3 | 1.38e-05 | 1.06 | 512 | 2 | 3.12e-05 | 1.52 | 1.43x |
| double | 243 | 243 | balanced | 1e-08 | adjoint | 320 | 7 | 1.33e-10 | 2.39 | 512 | 5 | 1.07e-10 | 2.53 | 1.06x |
| double | 243 | 243 | balanced | 1e-08 | forward | 320 | 6 | 9.60e-10 | 1.37 | 512 | 4 | 3.65e-09 | 1.60 | 1.17x |
| double | 243 | 243 | balanced | 1e-11 | adjoint | 384 | 7 | 1.12e-12 | 2.38 | 512 | 6 | 1.05e-12 | 2.53 | 1.06x |
| double | 243 | 243 | balanced | 1e-11 | forward | 324 | 8 | 1.03e-12 | 1.98 | 512 | 6 | 3.70e-13 | 1.84 | 0.93x |
| double | 243 | 972 | M-dominated | 1e-04 | adjoint | 320 | 3 | 1.32e-05 | 4.60 | 512 | 2 | 4.25e-05 | 4.16 | 0.90x |
| double | 243 | 972 | M-dominated | 1e-04 | forward | 320 | 3 | 1.71e-05 | 2.63 | 512 | 2 | 3.95e-05 | 2.51 | 0.96x |
| double | 243 | 972 | M-dominated | 1e-08 | adjoint | 384 | 5 | 1.30e-09 | 6.26 | 512 | 4 | 3.85e-09 | 6.00 | 0.96x |
| double | 243 | 972 | M-dominated | 1e-08 | forward | 384 | 5 | 9.95e-10 | 3.55 | 512 | 4 | 3.39e-09 | 3.30 | 0.93x |
| double | 243 | 972 | M-dominated | 1e-11 | adjoint | 512 | 6 | 5.33e-13 | 7.48 | 512 | 6 | 5.33e-13 | 7.24 | 0.97x |
| double | 243 | 972 | M-dominated | 1e-11 | forward | 512 | 6 | 3.76e-13 | 4.60 | 512 | 6 | 3.76e-13 | 4.38 | 0.95x |
| double | 250 | 62 | N-dominated | 1e-04 | adjoint | 320 | 4 | 1.61e-05 | 0.91 | 512 | 3 | 2.63e-06 | 1.28 | 1.40x |
| double | 250 | 62 | N-dominated | 1e-04 | forward | 320 | 3 | 1.66e-05 | 0.74 | 512 | 2 | 3.08e-05 | 1.27 | 1.71x |
| double | 250 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 3.33e-09 | 0.96 | 512 | 5 | 4.40e-10 | 1.38 | 1.44x |
| double | 250 | 62 | N-dominated | 1e-08 | forward | 324 | 6 | 1.77e-09 | 0.84 | 512 | 4 | 3.12e-09 | 1.15 | 1.36x |
| double | 250 | 62 | N-dominated | 1e-11 | adjoint | 384 | 8 | 1.29e-13 | 1.24 | 512 | 6 | 5.76e-12 | 1.28 | 1.03x |
| double | 250 | 62 | N-dominated | 1e-11 | forward | 320 | 9 | 4.38e-13 | 0.91 | 512 | 6 | 4.76e-13 | 1.20 | 1.32x |
| double | 250 | 250 | balanced | 1e-04 | adjoint | 384 | 3 | 7.52e-06 | 1.79 | 512 | 2 | 7.73e-05 | 1.84 | 1.03x |
| double | 250 | 250 | balanced | 1e-04 | forward | 320 | 3 | 1.96e-05 | 1.11 | 512 | 2 | 3.88e-05 | 1.51 | 1.37x |
| double | 250 | 250 | balanced | 1e-08 | adjoint | 384 | 6 | 9.38e-11 | 2.39 | 512 | 5 | 1.27e-10 | 2.46 | 1.03x |
| double | 250 | 250 | balanced | 1e-08 | forward | 324 | 6 | 1.62e-09 | 1.53 | 512 | 4 | 3.84e-09 | 1.71 | 1.12x |
| double | 250 | 250 | balanced | 1e-11 | adjoint | 400 | 7 | 1.90e-12 | 2.76 | 512 | 6 | 1.70e-12 | 2.62 | 0.95x |
| double | 250 | 250 | balanced | 1e-11 | forward | 384 | 7 | 9.34e-13 | 1.61 | 512 | 6 | 4.63e-13 | 1.98 | 1.23x |
| double | 250 | 1000 | M-dominated | 1e-04 | adjoint | 384 | 3 | 4.63e-06 | 4.66 | 512 | 2 | 5.26e-05 | 4.23 | 0.91x |
| double | 250 | 1000 | M-dominated | 1e-04 | forward | 320 | 3 | 2.49e-05 | 2.45 | 512 | 2 | 3.83e-05 | 2.42 | 0.99x |
| double | 250 | 1000 | M-dominated | 1e-08 | adjoint | 400 | 5 | 2.34e-09 | 6.38 | 512 | 4 | 5.02e-09 | 5.80 | 0.91x |
| double | 250 | 1000 | M-dominated | 1e-08 | forward | 384 | 5 | 1.75e-09 | 3.66 | 512 | 4 | 5.57e-09 | 3.27 | 0.89x |
| double | 250 | 1000 | M-dominated | 1e-11 | adjoint | 512 | 6 | 7.71e-13 | 7.38 | 512 | 6 | 7.71e-13 | 7.37 | 1.00x |
| double | 250 | 1000 | M-dominated | 1e-11 | forward | 512 | 6 | 5.19e-13 | 4.33 | 512 | 6 | 5.19e-13 | 4.44 | 1.03x |
| double | 251 | 62 | N-dominated | 1e-04 | adjoint | 320 | 4 | 5.95e-06 | 0.93 | 512 | 3 | 1.45e-06 | 1.30 | 1.40x |
| double | 251 | 62 | N-dominated | 1e-04 | forward | 320 | 3 | 1.81e-05 | 0.75 | 512 | 2 | 3.10e-05 | 1.11 | 1.48x |
| double | 251 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 6.15e-10 | 1.08 | 512 | 5 | 1.58e-10 | 1.38 | 1.28x |
| double | 251 | 62 | N-dominated | 1e-08 | forward | 324 | 6 | 1.55e-09 | 0.90 | 512 | 4 | 3.42e-09 | 1.17 | 1.29x |
| double | 251 | 62 | N-dominated | 1e-11 | adjoint | 384 | 8 | 9.40e-14 | 1.28 | 512 | 6 | 1.41e-12 | 1.43 | 1.12x |
| double | 251 | 62 | N-dominated | 1e-11 | forward | 320 | 9 | 2.24e-13 | 0.95 | 512 | 6 | 3.99e-13 | 1.29 | 1.35x |
| double | 251 | 251 | balanced | 1e-04 | adjoint | 384 | 3 | 7.49e-06 | 1.68 | 512 | 2 | 7.46e-05 | 1.79 | 1.06x |
| double | 251 | 251 | balanced | 1e-04 | forward | 320 | 3 | 1.89e-05 | 1.20 | 512 | 2 | 3.60e-05 | 1.56 | 1.31x |
| double | 251 | 251 | balanced | 1e-08 | adjoint | 384 | 6 | 1.11e-10 | 2.38 | 512 | 4 | 8.46e-09 | 2.09 | 0.88x |
| double | 251 | 251 | balanced | 1e-08 | forward | 384 | 5 | 1.33e-09 | 1.33 | 512 | 4 | 5.49e-09 | 1.58 | 1.19x |
| double | 251 | 251 | balanced | 1e-11 | adjoint | 400 | 7 | 8.81e-13 | 2.66 | 512 | 6 | 9.73e-13 | 2.67 | 1.00x |
| double | 251 | 251 | balanced | 1e-11 | forward | 384 | 7 | 6.60e-13 | 1.65 | 512 | 6 | 4.63e-13 | 1.92 | 1.17x |
| double | 251 | 1004 | M-dominated | 1e-04 | adjoint | 384 | 3 | 3.64e-06 | 4.67 | 512 | 2 | 4.99e-05 | 4.01 | 0.86x |
| double | 251 | 1004 | M-dominated | 1e-04 | forward | 324 | 3 | 1.95e-05 | 2.40 | 512 | 2 | 3.88e-05 | 2.48 | 1.03x |
| double | 251 | 1004 | M-dominated | 1e-08 | adjoint | 400 | 5 | 1.90e-09 | 6.21 | 512 | 4 | 4.42e-09 | 5.68 | 0.92x |
| double | 251 | 1004 | M-dominated | 1e-08 | forward | 384 | 5 | 1.32e-09 | 3.81 | 512 | 4 | 3.98e-09 | 3.19 | 0.84x |
| double | 251 | 1004 | M-dominated | 1e-11 | adjoint | 512 | 6 | 4.47e-13 | 7.14 | 512 | 6 | 4.47e-13 | 7.16 | 1.00x |
| double | 251 | 1004 | M-dominated | 1e-11 | forward | 512 | 6 | 4.36e-13 | 3.93 | 512 | 6 | 4.36e-13 | 4.17 | 1.06x |
| double | 255 | 63 | N-dominated | 1e-04 | adjoint | 320 | 4 | 7.87e-06 | 0.96 | 512 | 3 | 1.46e-06 | 1.27 | 1.33x |
| double | 255 | 63 | N-dominated | 1e-04 | forward | 324 | 3 | 2.22e-05 | 0.79 | 512 | 2 | 3.45e-05 | 1.09 | 1.39x |
| double | 255 | 63 | N-dominated | 1e-08 | adjoint | 324 | 7 | 1.53e-09 | 1.25 | 512 | 5 | 2.36e-10 | 1.50 | 1.20x |
| double | 255 | 63 | N-dominated | 1e-08 | forward | 320 | 7 | 1.88e-10 | 0.84 | 512 | 4 | 3.71e-09 | 1.22 | 1.45x |
| double | 255 | 63 | N-dominated | 1e-11 | adjoint | 360 | 9 | 3.50e-14 | 1.44 | 512 | 6 | 2.72e-12 | 1.44 | 1.01x |
| double | 255 | 63 | N-dominated | 1e-11 | forward | 324 | 9 | 3.46e-13 | 1.07 | 512 | 6 | 4.83e-13 | 1.27 | 1.19x |
| double | 255 | 255 | balanced | 1e-04 | adjoint | 384 | 3 | 6.00e-06 | 1.77 | 512 | 2 | 7.59e-05 | 1.90 | 1.07x |
| double | 255 | 255 | balanced | 1e-04 | forward | 324 | 3 | 1.81e-05 | 1.15 | 512 | 2 | 3.63e-05 | 1.35 | 1.17x |
| double | 255 | 255 | balanced | 1e-08 | adjoint | 384 | 6 | 8.43e-11 | 2.37 | 512 | 4 | 7.42e-09 | 2.21 | 0.93x |
| double | 255 | 255 | balanced | 1e-08 | forward | 400 | 5 | 1.03e-09 | 1.52 | 512 | 4 | 4.22e-09 | 1.57 | 1.03x |
| double | 255 | 255 | balanced | 1e-11 | adjoint | 384 | 8 | 7.76e-14 | 2.75 | 512 | 6 | 8.22e-13 | 2.64 | 0.96x |
| double | 255 | 255 | balanced | 1e-11 | forward | 384 | 7 | 1.11e-12 | 1.60 | 512 | 6 | 4.95e-13 | 1.86 | 1.16x |
| double | 255 | 1020 | M-dominated | 1e-04 | adjoint | 384 | 3 | 4.34e-06 | 4.94 | 512 | 2 | 5.99e-05 | 4.11 | 0.83x |
| double | 255 | 1020 | M-dominated | 1e-04 | forward | 384 | 3 | 3.84e-06 | 2.60 | 512 | 2 | 4.10e-05 | 2.50 | 0.96x |
| double | 255 | 1020 | M-dominated | 1e-08 | adjoint | 400 | 5 | 1.14e-09 | 6.44 | 512 | 4 | 5.71e-09 | 5.76 | 0.89x |
| double | 255 | 1020 | M-dominated | 1e-08 | forward | 400 | 5 | 1.11e-09 | 3.71 | 512 | 4 | 5.29e-09 | 3.46 | 0.93x |
| double | 255 | 1020 | M-dominated | 1e-11 | adjoint | 512 | 6 | 5.97e-13 | 7.47 | 512 | 6 | 5.97e-13 | 7.39 | 0.99x |
| double | 255 | 1020 | M-dominated | 1e-11 | forward | 512 | 6 | 7.14e-13 | 4.52 | 512 | 6 | 7.14e-13 | 4.43 | 0.98x |
| double | 256 | 64 | N-dominated | 1e-04 | adjoint | 320 | 4 | 1.41e-05 | 0.92 | 512 | 3 | 2.10e-06 | 1.26 | 1.37x |
| double | 256 | 64 | N-dominated | 1e-04 | forward | 324 | 3 | 2.63e-05 | 1.23 | 512 | 2 | 3.23e-05 | 1.10 | 0.90x |
| double | 256 | 64 | N-dominated | 1e-08 | adjoint | 324 | 7 | 1.02e-09 | 1.07 | 512 | 5 | 1.95e-10 | 1.32 | 1.23x |
| double | 256 | 64 | N-dominated | 1e-08 | forward | 320 | 7 | 2.87e-10 | 0.84 | 512 | 4 | 3.99e-09 | 1.08 | 1.28x |
| double | 256 | 64 | N-dominated | 1e-11 | adjoint | 360 | 9 | 1.31e-13 | 1.38 | 512 | 6 | 2.62e-12 | 1.39 | 1.01x |
| double | 256 | 64 | N-dominated | 1e-11 | forward | 324 | 9 | 5.94e-13 | 0.92 | 512 | 6 | 5.27e-13 | 1.26 | 1.37x |
| double | 256 | 256 | balanced | 1e-04 | adjoint | 384 | 3 | 1.51e-05 | 1.75 | 512 | 2 | 9.43e-05 | 1.82 | 1.04x |
| double | 256 | 256 | balanced | 1e-04 | forward | 384 | 3 | 2.87e-06 | 1.27 | 512 | 2 | 3.49e-05 | 1.38 | 1.08x |
| double | 256 | 256 | balanced | 1e-08 | adjoint | 384 | 6 | 2.96e-10 | 2.33 | 512 | 5 | 2.10e-10 | 2.43 | 1.04x |
| double | 256 | 256 | balanced | 1e-08 | forward | 400 | 5 | 1.05e-09 | 1.47 | 512 | 4 | 4.30e-09 | 1.46 | 0.99x |
| double | 256 | 256 | balanced | 1e-11 | adjoint | 384 | 8 | 2.46e-13 | 2.73 | 512 | 6 | 2.96e-12 | 2.61 | 0.95x |
| double | 256 | 256 | balanced | 1e-11 | forward | 384 | 7 | 1.32e-12 | 1.61 | 512 | 6 | 5.83e-13 | 1.82 | 1.13x |
| double | 256 | 1024 | M-dominated | 1e-04 | adjoint | 384 | 3 | 4.36e-06 | 4.75 | 512 | 2 | 3.93e-05 | 4.41 | 0.93x |
| double | 256 | 1024 | M-dominated | 1e-04 | forward | 384 | 3 | 3.91e-06 | 2.56 | 512 | 2 | 3.86e-05 | 2.49 | 0.97x |
| double | 256 | 1024 | M-dominated | 1e-08 | adjoint | 400 | 5 | 3.44e-09 | 6.45 | 512 | 4 | 4.28e-09 | 5.91 | 0.92x |
| double | 256 | 1024 | M-dominated | 1e-08 | forward | 400 | 5 | 1.35e-09 | 3.69 | 512 | 4 | 4.85e-09 | 3.37 | 0.91x |
| double | 256 | 1024 | M-dominated | 1e-11 | adjoint | 512 | 6 | 5.31e-13 | 7.66 | 512 | 6 | 5.31e-13 | 8.27 | 1.08x |
| double | 256 | 1024 | M-dominated | 1e-11 | forward | 512 | 6 | 6.43e-13 | 4.34 | 512 | 6 | 6.43e-13 | 4.76 | 1.10x |
| double | 512 | 128 | N-dominated | 1e-04 | adjoint | 640 | 4 | 1.63e-05 | 2.06 | 1024 | 3 | 2.01e-06 | 2.83 | 1.38x |
| double | 512 | 128 | N-dominated | 1e-04 | forward | 640 | 3 | 2.14e-05 | 2.06 | 1024 | 2 | 2.42e-05 | 2.69 | 1.30x |
| double | 512 | 128 | N-dominated | 1e-08 | adjoint | 648 | 7 | 2.05e-09 | 2.53 | 1024 | 5 | 3.69e-10 | 3.12 | 1.23x |
| double | 512 | 128 | N-dominated | 1e-08 | forward | 640 | 7 | 2.33e-10 | 1.96 | 1024 | 4 | 3.17e-09 | 2.57 | 1.31x |
| double | 512 | 128 | N-dominated | 1e-11 | adjoint | 720 | 9 | 1.15e-13 | 2.90 | 1024 | 6 | 5.02e-12 | 3.14 | 1.09x |
| double | 512 | 128 | N-dominated | 1e-11 | forward | 648 | 9 | 3.51e-13 | 2.98 | 1024 | 6 | 3.78e-13 | 2.73 | 0.91x |
| double | 512 | 512 | balanced | 1e-04 | adjoint | 640 | 4 | 6.36e-06 | 3.79 | 1024 | 2 | 6.25e-05 | 3.96 | 1.04x |
| double | 512 | 512 | balanced | 1e-04 | forward | 640 | 3 | 2.02e-05 | 2.28 | 1024 | 2 | 2.87e-05 | 2.97 | 1.30x |
| double | 512 | 512 | balanced | 1e-08 | adjoint | 640 | 7 | 1.91e-09 | 5.01 | 1024 | 4 | 6.93e-09 | 4.68 | 0.93x |
| double | 512 | 512 | balanced | 1e-08 | forward | 768 | 5 | 1.73e-09 | 3.17 | 1024 | 4 | 3.73e-09 | 3.36 | 1.06x |
| double | 512 | 512 | balanced | 1e-11 | adjoint | 800 | 7 | 1.16e-12 | 5.62 | 1024 | 6 | 1.01e-12 | 5.48 | 0.97x |
| double | 512 | 512 | balanced | 1e-11 | forward | 768 | 7 | 9.88e-13 | 3.67 | 1024 | 6 | 5.31e-13 | 4.08 | 1.11x |
| double | 512 | 2048 | M-dominated | 1e-04 | adjoint | 1024 | 2 | 4.31e-05 | 8.70 | 1024 | 2 | 4.31e-05 | 8.63 | 0.99x |
| double | 512 | 2048 | M-dominated | 1e-04 | forward | 640 | 3 | 2.51e-05 | 5.30 | 1024 | 2 | 3.19e-05 | 5.24 | 0.99x |
| double | 512 | 2048 | M-dominated | 1e-08 | adjoint | 800 | 5 | 1.28e-09 | 13.46 | 1024 | 4 | 6.71e-09 | 11.76 | 0.87x |
| double | 512 | 2048 | M-dominated | 1e-08 | forward | 768 | 5 | 1.73e-09 | 7.82 | 1024 | 4 | 4.14e-09 | 7.08 | 0.91x |
| double | 512 | 2048 | M-dominated | 1e-11 | adjoint | 1024 | 6 | 7.75e-13 | 15.17 | 1024 | 6 | 7.75e-13 | 15.03 | 0.99x |
| double | 512 | 2048 | M-dominated | 1e-11 | forward | 1024 | 6 | 5.64e-13 | 8.93 | 1024 | 6 | 5.64e-13 | 8.80 | 0.99x |
| long double | 243 | 60 | N-dominated | 1e-08 | adjoint | 320 | 7 | 2.29e-10 | 93.77 | 512 | 5 | 1.24e-10 | 124.79 | 1.33x |
| long double | 243 | 60 | N-dominated | 1e-08 | forward | 320 | 6 | 8.63e-10 | 86.23 | 512 | 4 | 2.84e-09 | 124.77 | 1.45x |
| long double | 243 | 60 | N-dominated | 1e-14 | adjoint | 324 | 11 | 9.55e-16 | 115.19 | 512 | 8 | 1.48e-16 | 143.37 | 1.24x |
| long double | 243 | 60 | N-dominated | 1e-14 | forward | 320 | 11 | 1.41e-16 | 104.43 | 512 | 7 | 3.21e-15 | 129.96 | 1.24x |
| long double | 243 | 60 | N-dominated | 1e-20 | adjoint | 320 | 16 | 3.59e-22 | 121.68 | 512 | 11 | 1.90e-22 | 152.40 | 1.25x |
| long double | 243 | 60 | N-dominated | 1e-20 | forward | 320 | 15 | 8.35e-22 | 123.37 | 512 | 10 | 3.07e-21 | 149.37 | 1.21x |
| long double | 243 | 243 | balanced | 1e-08 | adjoint | 320 | 7 | 1.33e-10 | 183.41 | 512 | 5 | 1.07e-10 | 223.50 | 1.22x |
| long double | 243 | 243 | balanced | 1e-08 | forward | 320 | 6 | 9.60e-10 | 179.88 | 512 | 4 | 3.65e-09 | 197.09 | 1.10x |
| long double | 243 | 243 | balanced | 1e-14 | adjoint | 384 | 9 | 5.89e-16 | 258.84 | 512 | 7 | 7.20e-15 | 222.37 | 0.86x |
| long double | 243 | 243 | balanced | 1e-14 | forward | 384 | 9 | 1.37e-16 | 273.26 | 512 | 7 | 3.28e-15 | 255.41 | 0.93x |
| long double | 243 | 243 | balanced | 1e-20 | adjoint | 384 | 13 | 1.94e-22 | 301.26 | 512 | 11 | 1.12e-22 | 294.44 | 0.98x |
| long double | 243 | 243 | balanced | 1e-20 | forward | 400 | 12 | 3.97e-22 | 345.03 | 512 | 10 | 3.16e-21 | 286.54 | 0.83x |
| long double | 243 | 972 | M-dominated | 1e-08 | adjoint | 384 | 5 | 1.30e-09 | 507.24 | 512 | 4 | 3.85e-09 | 472.44 | 0.93x |
| long double | 243 | 972 | M-dominated | 1e-08 | forward | 384 | 5 | 9.95e-10 | 523.45 | 512 | 4 | 3.39e-09 | 477.60 | 0.91x |
| long double | 243 | 972 | M-dominated | 1e-14 | adjoint | 576 | 7 | 9.69e-16 | 734.92 | 512 | 7 | 5.32e-15 | 738.55 | 1.00x |
| long double | 243 | 972 | M-dominated | 1e-14 | forward | 576 | 7 | 8.15e-16 | 766.98 | 512 | 7 | 3.47e-15 | 772.71 | 1.01x |
| long double | 243 | 972 | M-dominated | 1e-20 | adjoint | 576 | 10 | 8.15e-22 | 963.90 | 512 | 10 | 4.07e-21 | 901.15 | 0.93x |
| long double | 243 | 972 | M-dominated | 1e-20 | forward | 576 | 10 | 3.67e-22 | 1028.73 | 512 | 10 | 3.77e-21 | 1018.46 | 0.99x |
| long double | 250 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 3.33e-09 | 90.97 | 512 | 5 | 4.40e-10 | 122.79 | 1.35x |
| long double | 250 | 62 | N-dominated | 1e-08 | forward | 324 | 6 | 1.77e-09 | 101.53 | 512 | 4 | 3.12e-09 | 134.42 | 1.32x |
| long double | 250 | 62 | N-dominated | 1e-14 | adjoint | 320 | 12 | 1.73e-15 | 109.71 | 512 | 8 | 8.25e-16 | 146.60 | 1.34x |
| long double | 250 | 62 | N-dominated | 1e-14 | forward | 320 | 11 | 9.79e-16 | 107.51 | 512 | 7 | 4.64e-15 | 138.08 | 1.28x |
| long double | 250 | 62 | N-dominated | 1e-20 | adjoint | 320 | 17 | 8.50e-22 | 127.78 | 512 | 11 | 1.34e-21 | 145.26 | 1.14x |
| long double | 250 | 62 | N-dominated | 1e-20 | forward | 320 | 16 | 4.50e-22 | 131.32 | 512 | 10 | 7.38e-21 | 143.65 | 1.09x |
| long double | 250 | 250 | balanced | 1e-08 | adjoint | 384 | 6 | 9.38e-11 | 209.65 | 512 | 5 | 1.27e-10 | 215.40 | 1.03x |
| long double | 250 | 250 | balanced | 1e-08 | forward | 324 | 6 | 1.62e-09 | 205.68 | 512 | 4 | 3.84e-09 | 178.40 | 0.87x |
| long double | 250 | 250 | balanced | 1e-14 | adjoint | 400 | 9 | 9.77e-16 | 284.71 | 512 | 8 | 2.43e-16 | 279.45 | 0.98x |
| long double | 250 | 250 | balanced | 1e-14 | forward | 384 | 9 | 4.69e-16 | 275.65 | 512 | 7 | 4.82e-15 | 229.23 | 0.83x |
| long double | 250 | 250 | balanced | 1e-20 | adjoint | 384 | 13 | 6.87e-22 | 333.20 | 512 | 11 | 3.93e-22 | 289.08 | 0.87x |
| long double | 250 | 250 | balanced | 1e-20 | forward | 400 | 12 | 1.69e-21 | 361.62 | 512 | 10 | 6.12e-21 | 287.59 | 0.80x |
| long double | 250 | 1000 | M-dominated | 1e-08 | adjoint | 400 | 5 | 2.34e-09 | 533.40 | 512 | 4 | 5.02e-09 | 506.44 | 0.95x |
| long double | 250 | 1000 | M-dominated | 1e-08 | forward | 384 | 5 | 1.75e-09 | 555.61 | 512 | 4 | 5.57e-09 | 474.05 | 0.85x |
| long double | 250 | 1000 | M-dominated | 1e-14 | adjoint | 576 | 7 | 3.84e-15 | 753.11 | 512 | 7 | 9.32e-15 | 712.63 | 0.95x |
| long double | 250 | 1000 | M-dominated | 1e-14 | forward | 576 | 7 | 1.41e-15 | 824.08 | 512 | 7 | 5.79e-15 | 753.12 | 0.91x |
| long double | 250 | 1000 | M-dominated | 1e-20 | adjoint | 576 | 10 | 3.02e-21 | 969.34 | 512 | 11 | 1.80e-22 | 1013.47 | 1.05x |
| long double | 250 | 1000 | M-dominated | 1e-20 | forward | 576 | 10 | 9.46e-22 | 1056.55 | 512 | 10 | 9.36e-21 | 1049.85 | 0.99x |
| long double | 251 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 6.15e-10 | 90.91 | 512 | 5 | 1.58e-10 | 123.09 | 1.35x |
| long double | 251 | 62 | N-dominated | 1e-08 | forward | 324 | 6 | 1.55e-09 | 100.48 | 512 | 4 | 3.42e-09 | 134.94 | 1.34x |
| long double | 251 | 62 | N-dominated | 1e-14 | adjoint | 320 | 12 | 3.67e-16 | 110.56 | 512 | 8 | 2.32e-16 | 145.80 | 1.32x |
| long double | 251 | 62 | N-dominated | 1e-14 | forward | 320 | 11 | 6.10e-16 | 107.33 | 512 | 7 | 4.68e-15 | 133.07 | 1.24x |
| long double | 251 | 62 | N-dominated | 1e-20 | adjoint | 320 | 17 | 1.61e-22 | 129.15 | 512 | 11 | 2.85e-22 | 157.23 | 1.22x |
| long double | 251 | 62 | N-dominated | 1e-20 | forward | 320 | 16 | 3.49e-22 | 130.45 | 512 | 10 | 5.17e-21 | 149.33 | 1.14x |
| long double | 251 | 251 | balanced | 1e-08 | adjoint | 384 | 6 | 1.11e-10 | 204.39 | 512 | 4 | 8.46e-09 | 180.88 | 0.88x |
| long double | 251 | 251 | balanced | 1e-08 | forward | 384 | 5 | 1.33e-09 | 191.70 | 512 | 4 | 5.49e-09 | 184.30 | 0.96x |
| long double | 251 | 251 | balanced | 1e-14 | adjoint | 400 | 9 | 4.66e-16 | 251.36 | 512 | 7 | 7.51e-15 | 259.76 | 1.03x |
| long double | 251 | 251 | balanced | 1e-14 | forward | 384 | 9 | 3.44e-16 | 283.04 | 512 | 7 | 4.78e-15 | 270.97 | 0.96x |
| long double | 251 | 251 | balanced | 1e-20 | adjoint | 384 | 13 | 6.71e-22 | 304.27 | 512 | 11 | 1.48e-22 | 327.04 | 1.07x |
| long double | 251 | 251 | balanced | 1e-20 | forward | 384 | 13 | 1.34e-22 | 337.41 | 512 | 10 | 5.38e-21 | 300.01 | 0.89x |
| long double | 251 | 1004 | M-dominated | 1e-08 | adjoint | 400 | 5 | 1.90e-09 | 548.24 | 512 | 4 | 4.42e-09 | 553.82 | 1.01x |
| long double | 251 | 1004 | M-dominated | 1e-08 | forward | 384 | 5 | 1.32e-09 | 558.65 | 512 | 4 | 3.98e-09 | 562.42 | 1.01x |
| long double | 251 | 1004 | M-dominated | 1e-14 | adjoint | 576 | 7 | 7.66e-16 | 761.10 | 512 | 7 | 4.88e-15 | 711.73 | 0.94x |
| long double | 251 | 1004 | M-dominated | 1e-14 | forward | 576 | 7 | 1.15e-15 | 801.08 | 512 | 7 | 5.38e-15 | 751.28 | 0.94x |
| long double | 251 | 1004 | M-dominated | 1e-20 | adjoint | 576 | 10 | 8.47e-22 | 998.43 | 512 | 11 | 9.87e-23 | 1038.92 | 1.04x |
| long double | 251 | 1004 | M-dominated | 1e-20 | forward | 576 | 10 | 7.27e-22 | 1099.00 | 512 | 10 | 7.31e-21 | 1014.47 | 0.92x |
| long double | 255 | 63 | N-dominated | 1e-08 | adjoint | 324 | 7 | 1.53e-09 | 105.14 | 512 | 5 | 2.36e-10 | 122.95 | 1.17x |
| long double | 255 | 63 | N-dominated | 1e-08 | forward | 320 | 7 | 1.88e-10 | 92.89 | 512 | 4 | 3.71e-09 | 134.14 | 1.44x |
| long double | 255 | 63 | N-dominated | 1e-14 | adjoint | 324 | 12 | 1.13e-15 | 125.56 | 512 | 8 | 4.52e-16 | 147.89 | 1.18x |
| long double | 255 | 63 | N-dominated | 1e-14 | forward | 324 | 11 | 7.83e-16 | 125.81 | 512 | 7 | 4.88e-15 | 136.10 | 1.08x |
| long double | 255 | 63 | N-dominated | 1e-20 | adjoint | 320 | 17 | 6.20e-22 | 131.92 | 512 | 11 | 4.62e-22 | 153.77 | 1.17x |
| long double | 255 | 63 | N-dominated | 1e-20 | forward | 324 | 16 | 5.01e-22 | 149.12 | 512 | 10 | 7.43e-21 | 147.51 | 0.99x |
| long double | 255 | 255 | balanced | 1e-08 | adjoint | 384 | 6 | 8.43e-11 | 203.73 | 512 | 4 | 7.42e-09 | 179.23 | 0.88x |
| long double | 255 | 255 | balanced | 1e-08 | forward | 400 | 5 | 1.03e-09 | 202.10 | 512 | 4 | 4.22e-09 | 178.60 | 0.88x |
| long double | 255 | 255 | balanced | 1e-14 | adjoint | 400 | 9 | 8.14e-16 | 275.07 | 512 | 8 | 1.20e-16 | 280.96 | 1.02x |
| long double | 255 | 255 | balanced | 1e-14 | forward | 384 | 9 | 5.75e-16 | 285.52 | 512 | 7 | 6.26e-15 | 269.30 | 0.94x |
| long double | 255 | 255 | balanced | 1e-20 | adjoint | 432 | 12 | 3.53e-22 | 338.93 | 512 | 11 | 1.80e-22 | 310.94 | 0.92x |
| long double | 255 | 255 | balanced | 1e-20 | forward | 384 | 13 | 2.47e-22 | 337.95 | 512 | 10 | 7.00e-21 | 334.71 | 0.99x |
| long double | 255 | 1020 | M-dominated | 1e-08 | adjoint | 400 | 5 | 1.14e-09 | 536.32 | 512 | 4 | 5.71e-09 | 486.97 | 0.91x |
| long double | 255 | 1020 | M-dominated | 1e-08 | forward | 400 | 5 | 1.11e-09 | 561.75 | 512 | 4 | 5.29e-09 | 490.05 | 0.87x |
| long double | 255 | 1020 | M-dominated | 1e-14 | adjoint | 576 | 7 | 1.23e-15 | 772.34 | 512 | 7 | 7.22e-15 | 718.03 | 0.93x |
| long double | 255 | 1020 | M-dominated | 1e-14 | forward | 576 | 7 | 1.55e-15 | 835.46 | 512 | 7 | 6.52e-15 | 772.57 | 0.92x |
| long double | 255 | 1020 | M-dominated | 1e-20 | adjoint | 576 | 10 | 7.33e-22 | 990.84 | 512 | 11 | 1.52e-22 | 1051.35 | 1.06x |
| long double | 255 | 1020 | M-dominated | 1e-20 | forward | 576 | 10 | 9.49e-22 | 1082.77 | 512 | 10 | 9.17e-21 | 1024.10 | 0.95x |
| long double | 256 | 64 | N-dominated | 1e-08 | adjoint | 324 | 7 | 1.02e-09 | 103.70 | 512 | 5 | 1.95e-10 | 126.59 | 1.22x |
| long double | 256 | 64 | N-dominated | 1e-08 | forward | 320 | 7 | 2.87e-10 | 110.38 | 512 | 4 | 3.99e-09 | 129.90 | 1.18x |
| long double | 256 | 64 | N-dominated | 1e-14 | adjoint | 324 | 12 | 6.49e-16 | 127.80 | 512 | 8 | 4.08e-16 | 134.20 | 1.05x |
| long double | 256 | 64 | N-dominated | 1e-14 | forward | 324 | 11 | 2.48e-15 | 133.59 | 512 | 7 | 6.97e-15 | 153.84 | 1.15x |
| long double | 256 | 64 | N-dominated | 1e-20 | adjoint | 324 | 17 | 5.04e-22 | 146.08 | 512 | 11 | 7.69e-22 | 151.34 | 1.04x |
| long double | 256 | 64 | N-dominated | 1e-20 | forward | 324 | 16 | 1.31e-21 | 148.97 | 512 | 11 | 1.51e-22 | 152.33 | 1.02x |
| long double | 256 | 256 | balanced | 1e-08 | adjoint | 384 | 6 | 2.96e-10 | 204.21 | 512 | 5 | 2.10e-10 | 217.23 | 1.06x |
| long double | 256 | 256 | balanced | 1e-08 | forward | 400 | 5 | 1.05e-09 | 212.30 | 512 | 4 | 4.30e-09 | 181.35 | 0.85x |
| long double | 256 | 256 | balanced | 1e-14 | adjoint | 400 | 9 | 9.32e-16 | 272.02 | 512 | 8 | 4.68e-16 | 280.59 | 1.03x |
| long double | 256 | 256 | balanced | 1e-14 | forward | 384 | 9 | 9.01e-16 | 287.30 | 512 | 7 | 7.25e-15 | 241.17 | 0.84x |
| long double | 256 | 256 | balanced | 1e-20 | adjoint | 432 | 12 | 1.02e-21 | 373.91 | 512 | 11 | 8.82e-22 | 307.93 | 0.82x |
| long double | 256 | 256 | balanced | 1e-20 | forward | 384 | 13 | 4.51e-22 | 361.14 | 512 | 11 | 1.31e-22 | 330.49 | 0.92x |
| long double | 256 | 1024 | M-dominated | 1e-08 | adjoint | 400 | 5 | 3.44e-09 | 559.75 | 512 | 4 | 4.28e-09 | 544.67 | 0.97x |
| long double | 256 | 1024 | M-dominated | 1e-08 | forward | 400 | 5 | 1.35e-09 | 568.43 | 512 | 4 | 4.85e-09 | 502.59 | 0.88x |
| long double | 256 | 1024 | M-dominated | 1e-14 | adjoint | 640 | 7 | 3.73e-16 | 786.44 | 512 | 7 | 7.31e-15 | 715.50 | 0.91x |
| long double | 256 | 1024 | M-dominated | 1e-14 | forward | 576 | 7 | 1.89e-15 | 814.53 | 512 | 7 | 8.76e-15 | 825.42 | 1.01x |
| long double | 256 | 1024 | M-dominated | 1e-20 | adjoint | 576 | 10 | 1.82e-21 | 990.32 | 512 | 11 | 1.61e-22 | 1046.56 | 1.06x |
| long double | 256 | 1024 | M-dominated | 1e-20 | forward | 576 | 10 | 1.46e-21 | 1088.10 | 512 | 11 | 1.55e-22 | 1149.74 | 1.06x |
| long double | 512 | 128 | N-dominated | 1e-08 | adjoint | 648 | 7 | 2.05e-09 | 336.51 | 1024 | 5 | 3.69e-10 | 385.74 | 1.15x |
| long double | 512 | 128 | N-dominated | 1e-08 | forward | 640 | 7 | 2.33e-10 | 256.91 | 1024 | 4 | 3.17e-09 | 422.07 | 1.64x |
| long double | 512 | 128 | N-dominated | 1e-14 | adjoint | 648 | 12 | 1.48e-15 | 365.73 | 1024 | 8 | 7.98e-16 | 447.57 | 1.22x |
| long double | 512 | 128 | N-dominated | 1e-14 | forward | 648 | 11 | 1.33e-15 | 371.55 | 1024 | 7 | 5.62e-15 | 373.85 | 1.01x |
| long double | 512 | 128 | N-dominated | 1e-20 | adjoint | 640 | 17 | 4.90e-21 | 292.47 | 1024 | 11 | 1.55e-21 | 384.30 | 1.31x |
| long double | 512 | 128 | N-dominated | 1e-20 | forward | 648 | 16 | 8.30e-22 | 431.51 | 1024 | 10 | 7.38e-21 | 491.50 | 1.14x |
| long double | 512 | 512 | balanced | 1e-08 | adjoint | 640 | 7 | 1.91e-09 | 445.96 | 1024 | 4 | 6.93e-09 | 614.16 | 1.38x |
| long double | 512 | 512 | balanced | 1e-08 | forward | 768 | 5 | 1.73e-09 | 479.22 | 1024 | 4 | 3.73e-09 | 614.06 | 1.28x |
| long double | 512 | 512 | balanced | 1e-14 | adjoint | 800 | 9 | 7.72e-16 | 659.03 | 1024 | 7 | 8.35e-15 | 615.39 | 0.93x |
| long double | 512 | 512 | balanced | 1e-14 | forward | 768 | 9 | 8.05e-16 | 626.65 | 1024 | 7 | 4.99e-15 | 642.49 | 1.03x |
| long double | 512 | 512 | balanced | 1e-20 | adjoint | 768 | 13 | 2.89e-21 | 733.64 | 1024 | 11 | 1.82e-22 | 789.35 | 1.08x |
| long double | 512 | 512 | balanced | 1e-20 | forward | 864 | 12 | 1.55e-22 | 860.56 | 1024 | 10 | 9.06e-21 | 783.12 | 0.91x |
| long double | 512 | 2048 | M-dominated | 1e-08 | adjoint | 800 | 5 | 1.28e-09 | 1328.81 | 1024 | 4 | 6.71e-09 | 1267.51 | 0.95x |
| long double | 512 | 2048 | M-dominated | 1e-08 | forward | 768 | 5 | 1.73e-09 | 1359.48 | 1024 | 4 | 4.14e-09 | 1238.56 | 0.91x |
| long double | 512 | 2048 | M-dominated | 1e-14 | adjoint | 1152 | 7 | 1.41e-15 | 1776.10 | 1024 | 7 | 7.48e-15 | 1670.59 | 0.94x |
| long double | 512 | 2048 | M-dominated | 1e-14 | forward | 864 | 8 | 1.89e-15 | 1961.03 | 1024 | 7 | 6.18e-15 | 1820.94 | 0.93x |
| long double | 512 | 2048 | M-dominated | 1e-20 | adjoint | 1152 | 10 | 1.13e-21 | 2215.95 | 1024 | 11 | 1.33e-22 | 2225.68 | 1.00x |
| long double | 512 | 2048 | M-dominated | 1e-20 | forward | 1152 | 10 | 9.91e-22 | 2438.20 | 1024 | 10 | 8.39e-21 | 2276.31 | 0.93x |

`!` marks a goal that was not met.

