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
With few nodes the tuner wins clearly (1.35x median, 96/96 cases); with many nodes it still trails (0.91x median, 10/96 wins).

`tune_plan` takes the node count and picks the pair its cost model rates
cheapest, `n*log2(n) + (4/5)*M*(2m+2)`, over every even 5-smooth `n` with
sigma in [5/4, 4]. The legacy size is a power of two in that range, so it
is always among the candidates. More nodes therefore buy a larger grid
and a smaller cut-off, which is what the M-dominated shape wants.

What remains is not the choice of pair but the error model. It is an
upper envelope over every measured geometry, while the oracle uses the
actual error of the one geometry in front of it, so the tuner hands out
about one more `m` than needed -- see the headroom figures below. At the
same `n` that alone costs `(2m+4)/(2m+2)`, and that is what is left in
the M-dominated shape, where the convolution is the whole bill.
Raising the cost weight does not recover it. On the double tree, a weight
of 1.5 instead of 4/5 moves the M-dominated median by less than 0.01x,
costs 0.10x in the N-dominated shape, and makes the worst case worse.

## Against the cost-blind policy

The first version of `tune_plan` did not take `M`. It minimised `n` and
then took the smallest `m` that worked there, which is a memory-first
policy. Same 288 configurations, same host:

| | cost-blind | cost-aware |
|---|---|---|
| overall median | 0.97x | 1.02x |
| N-dominated | 1.31x | 1.35x |
| balanced | 0.96x | 1.01x |
| M-dominated | 0.77x | 0.91x |
| worst single case | 0.60x | 0.68x |
| goal met | 288/288 | 288/288 |

## By problem shape

This is where the answer splits. `M` is the node count.

| shape | M | median speedup | new faster | median n_new / n_leg | median m_new - m_leg |
|---|---|---|---|---|---|
| N-dominated | N/4 | 1.35x | 96/96 | 0.62 | +2 |
| balanced | N | 1.01x | 49/96 | 0.70 | +1 |
| M-dominated | 4N | 0.91x | 10/96 | 0.78 | +1 |

## By precision

| precision | median speedup | new faster | goal met (new) | goal met (legacy) |
|---|---|---|---|---|
| float | 1.10x | 40/72 | 72/72 | 72/72 |
| double | 1.00x | 53/108 | 108/108 | 108/108 |
| long double | 1.02x | 62/108 | 108/108 | 108/108 |

## By bandwidth

| N | factors | median speedup | median n_new | n_leg | median m_new | median m_leg |
|---|---|---|---|---|---|---|
| 243 | 3·3·3·3·3 | 1.09x | 320 | 512 | 6 | 4 |
| 250 | 2·5·5·5 | 1.01x | 360 | 512 | 6 | 4 |
| 251 | 251 | 0.99x | 360 | 512 | 6 | 4 |
| 255 | 3·5·17 | 0.99x | 360 | 512 | 6 | 4 |
| 256 | 2·2·2·2·2·2·2·2 | 1.01x | 360 | 512 | 6 | 4 |
| 512 | 2·2·2·2·2·2·2·2·2 | 1.06x | 640 | 1024 | 7 | 4 |

## Accuracy

Both sides are asked for the same goal, so the question is not who is
more accurate but whether either misses. Cases where the goal was not
met:

- **New: none.**
- **Legacy (with oracle m): none.**

Median headroom below the goal: new 13.7x, legacy 3.1x.
Both overshoot the goal rather than miss it; the new side does so by
more, which is the cost of an upper-bound model against a measured search.

## Full results

| precision | N | M | shape | goal | dir | n_new | m_new | err_new | t_new (µs) | n_leg | m_leg | err_leg | t_leg (µs) | speedup |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| float | 243 | 60 | N-dominated | 1e-02 | adjoint | 320 | 2 | 2.13e-03 | 0.46 | 512 | 1 | 9.06e-03 | 0.56 | 1.22x |
| float | 243 | 60 | N-dominated | 1e-02 | forward | 320 | 2 | 2.87e-04 | 0.40 | 512 | 1 | 2.36e-03 | 0.57 | 1.42x |
| float | 243 | 60 | N-dominated | 1e-04 | adjoint | 320 | 4 | 4.78e-06 | 0.49 | 512 | 3 | 1.59e-06 | 0.60 | 1.22x |
| float | 243 | 60 | N-dominated | 1e-04 | forward | 320 | 3 | 1.15e-05 | 0.36 | 512 | 2 | 3.22e-05 | 0.57 | 1.58x |
| float | 243 | 243 | balanced | 1e-02 | adjoint | 320 | 2 | 1.38e-03 | 0.80 | 512 | 1 | 6.88e-03 | 0.91 | 1.14x |
| float | 243 | 243 | balanced | 1e-02 | forward | 320 | 2 | 4.07e-04 | 0.66 | 512 | 1 | 2.96e-03 | 0.76 | 1.16x |
| float | 243 | 243 | balanced | 1e-04 | adjoint | 320 | 4 | 2.68e-06 | 0.98 | 512 | 3 | 1.29e-06 | 1.10 | 1.12x |
| float | 243 | 243 | balanced | 1e-04 | forward | 320 | 3 | 1.37e-05 | 0.67 | 512 | 2 | 3.12e-05 | 0.84 | 1.25x |
| float | 243 | 972 | M-dominated | 1e-02 | adjoint | 320 | 2 | 3.98e-04 | 1.99 | 512 | 1 | 3.16e-03 | 1.95 | 0.98x |
| float | 243 | 972 | M-dominated | 1e-02 | forward | 320 | 2 | 3.39e-04 | 1.61 | 512 | 1 | 3.36e-03 | 1.57 | 0.98x |
| float | 243 | 972 | M-dominated | 1e-04 | adjoint | 360 | 3 | 4.51e-06 | 3.35 | 512 | 2 | 4.24e-05 | 2.40 | 0.72x |
| float | 243 | 972 | M-dominated | 1e-04 | forward | 320 | 3 | 1.69e-05 | 1.66 | 512 | 2 | 3.95e-05 | 1.82 | 1.09x |
| float | 250 | 62 | N-dominated | 1e-02 | adjoint | 320 | 2 | 4.46e-03 | 0.45 | 512 | 2 | 2.18e-04 | 0.70 | 1.55x |
| float | 250 | 62 | N-dominated | 1e-02 | forward | 320 | 2 | 3.80e-04 | 0.40 | 512 | 1 | 2.64e-03 | 0.57 | 1.42x |
| float | 250 | 62 | N-dominated | 1e-04 | adjoint | 320 | 4 | 1.63e-05 | 0.45 | 512 | 3 | 2.09e-06 | 0.60 | 1.34x |
| float | 250 | 62 | N-dominated | 1e-04 | forward | 320 | 4 | 9.95e-07 | 0.41 | 512 | 2 | 3.09e-05 | 0.59 | 1.42x |
| float | 250 | 250 | balanced | 1e-02 | adjoint | 320 | 2 | 1.72e-03 | 0.89 | 512 | 1 | 4.72e-03 | 0.91 | 1.01x |
| float | 250 | 250 | balanced | 1e-02 | forward | 320 | 2 | 4.15e-04 | 0.69 | 512 | 1 | 2.90e-03 | 0.77 | 1.12x |
| float | 250 | 250 | balanced | 1e-04 | adjoint | 360 | 3 | 1.19e-05 | 1.18 | 512 | 2 | 7.76e-05 | 0.90 | 0.77x |
| float | 250 | 250 | balanced | 1e-04 | forward | 360 | 3 | 4.74e-06 | 0.99 | 512 | 2 | 3.85e-05 | 0.87 | 0.88x |
| float | 250 | 1000 | M-dominated | 1e-02 | adjoint | 320 | 2 | 9.33e-04 | 2.31 | 512 | 1 | 3.13e-03 | 2.04 | 0.88x |
| float | 250 | 1000 | M-dominated | 1e-02 | forward | 320 | 2 | 4.88e-04 | 1.68 | 512 | 1 | 2.93e-03 | 1.56 | 0.93x |
| float | 250 | 1000 | M-dominated | 1e-04 | adjoint | 360 | 3 | 5.38e-06 | 3.12 | 512 | 2 | 5.26e-05 | 2.36 | 0.76x |
| float | 250 | 1000 | M-dominated | 1e-04 | forward | 360 | 3 | 6.65e-06 | 2.14 | 512 | 2 | 3.83e-05 | 1.88 | 0.88x |
| float | 251 | 62 | N-dominated | 1e-02 | adjoint | 320 | 2 | 2.07e-03 | 0.48 | 512 | 2 | 1.39e-04 | 0.62 | 1.29x |
| float | 251 | 62 | N-dominated | 1e-02 | forward | 320 | 2 | 3.43e-04 | 0.42 | 512 | 1 | 2.46e-03 | 0.58 | 1.37x |
| float | 251 | 62 | N-dominated | 1e-04 | adjoint | 320 | 4 | 6.15e-06 | 0.55 | 512 | 3 | 2.04e-06 | 0.60 | 1.10x |
| float | 251 | 62 | N-dominated | 1e-04 | forward | 320 | 4 | 1.06e-06 | 0.45 | 512 | 2 | 3.10e-05 | 0.61 | 1.36x |
| float | 251 | 251 | balanced | 1e-02 | adjoint | 320 | 2 | 1.56e-03 | 0.76 | 512 | 1 | 4.13e-03 | 0.87 | 1.14x |
| float | 251 | 251 | balanced | 1e-02 | forward | 320 | 2 | 4.69e-04 | 0.68 | 512 | 1 | 2.41e-03 | 0.81 | 1.20x |
| float | 251 | 251 | balanced | 1e-04 | adjoint | 360 | 3 | 1.74e-05 | 1.30 | 512 | 2 | 7.45e-05 | 1.04 | 0.80x |
| float | 251 | 251 | balanced | 1e-04 | forward | 360 | 3 | 5.24e-06 | 1.04 | 512 | 2 | 3.59e-05 | 0.89 | 0.86x |
| float | 251 | 1004 | M-dominated | 1e-02 | adjoint | 320 | 2 | 8.55e-04 | 2.31 | 512 | 1 | 2.78e-03 | 2.03 | 0.88x |
| float | 251 | 1004 | M-dominated | 1e-02 | forward | 320 | 2 | 4.08e-04 | 1.73 | 512 | 1 | 2.51e-03 | 1.55 | 0.90x |
| float | 251 | 1004 | M-dominated | 1e-04 | adjoint | 360 | 3 | 8.52e-06 | 3.08 | 512 | 2 | 5.01e-05 | 2.34 | 0.76x |
| float | 251 | 1004 | M-dominated | 1e-04 | forward | 360 | 3 | 5.62e-06 | 2.12 | 512 | 2 | 3.88e-05 | 1.82 | 0.86x |
| float | 255 | 63 | N-dominated | 1e-02 | adjoint | 320 | 2 | 2.88e-03 | 0.49 | 512 | 2 | 2.17e-04 | 0.67 | 1.39x |
| float | 255 | 63 | N-dominated | 1e-02 | forward | 320 | 2 | 4.80e-04 | 0.41 | 512 | 1 | 2.32e-03 | 0.60 | 1.46x |
| float | 255 | 63 | N-dominated | 1e-04 | adjoint | 320 | 4 | 7.73e-06 | 0.56 | 512 | 3 | 1.56e-06 | 0.91 | 1.63x |
| float | 255 | 63 | N-dominated | 1e-04 | forward | 320 | 4 | 1.31e-06 | 0.43 | 512 | 2 | 3.45e-05 | 0.58 | 1.37x |
| float | 255 | 255 | balanced | 1e-02 | adjoint | 320 | 2 | 2.41e-03 | 0.91 | 512 | 1 | 4.37e-03 | 0.95 | 1.04x |
| float | 255 | 255 | balanced | 1e-02 | forward | 320 | 2 | 4.46e-04 | 0.71 | 512 | 1 | 2.90e-03 | 0.80 | 1.13x |
| float | 255 | 255 | balanced | 1e-04 | adjoint | 360 | 3 | 1.69e-05 | 1.33 | 512 | 2 | 7.58e-05 | 0.93 | 0.70x |
| float | 255 | 255 | balanced | 1e-04 | forward | 360 | 3 | 6.01e-06 | 1.05 | 512 | 2 | 3.62e-05 | 0.88 | 0.84x |
| float | 255 | 1020 | M-dominated | 1e-02 | adjoint | 320 | 2 | 1.21e-03 | 2.38 | 512 | 1 | 3.23e-03 | 2.12 | 0.89x |
| float | 255 | 1020 | M-dominated | 1e-02 | forward | 320 | 2 | 6.12e-04 | 1.81 | 512 | 1 | 3.23e-03 | 1.47 | 0.81x |
| float | 255 | 1020 | M-dominated | 1e-04 | adjoint | 360 | 3 | 1.05e-05 | 3.19 | 512 | 2 | 5.99e-05 | 2.42 | 0.76x |
| float | 255 | 1020 | M-dominated | 1e-04 | forward | 360 | 3 | 7.71e-06 | 1.98 | 512 | 2 | 4.09e-05 | 1.86 | 0.94x |
| float | 256 | 64 | N-dominated | 1e-02 | adjoint | 320 | 2 | 5.68e-03 | 0.45 | 512 | 2 | 2.03e-04 | 0.61 | 1.35x |
| float | 256 | 64 | N-dominated | 1e-02 | forward | 320 | 2 | 5.60e-04 | 0.37 | 512 | 1 | 2.29e-03 | 0.52 | 1.41x |
| float | 256 | 64 | N-dominated | 1e-04 | adjoint | 320 | 4 | 1.48e-05 | 0.50 | 512 | 3 | 2.22e-06 | 0.72 | 1.44x |
| float | 256 | 64 | N-dominated | 1e-04 | forward | 320 | 4 | 1.53e-06 | 0.42 | 512 | 2 | 3.23e-05 | 0.62 | 1.48x |
| float | 256 | 256 | balanced | 1e-02 | adjoint | 320 | 2 | 1.70e-03 | 0.97 | 512 | 1 | 4.25e-03 | 0.92 | 0.95x |
| float | 256 | 256 | balanced | 1e-02 | forward | 320 | 2 | 4.65e-04 | 0.67 | 512 | 1 | 2.41e-03 | 0.78 | 1.16x |
| float | 256 | 256 | balanced | 1e-04 | adjoint | 320 | 4 | 7.61e-06 | 1.10 | 512 | 2 | 9.44e-05 | 1.00 | 0.91x |
| float | 256 | 256 | balanced | 1e-04 | forward | 360 | 3 | 6.11e-06 | 1.02 | 512 | 2 | 3.52e-05 | 0.83 | 0.81x |
| float | 256 | 1024 | M-dominated | 1e-02 | adjoint | 320 | 2 | 8.84e-04 | 2.40 | 512 | 1 | 2.31e-03 | 2.12 | 0.89x |
| float | 256 | 1024 | M-dominated | 1e-02 | forward | 320 | 2 | 5.44e-04 | 1.97 | 512 | 1 | 3.09e-03 | 1.58 | 0.80x |
| float | 256 | 1024 | M-dominated | 1e-04 | adjoint | 384 | 3 | 4.47e-06 | 3.15 | 512 | 2 | 3.92e-05 | 2.30 | 0.73x |
| float | 256 | 1024 | M-dominated | 1e-04 | forward | 360 | 3 | 6.99e-06 | 1.96 | 512 | 2 | 3.87e-05 | 1.76 | 0.90x |
| float | 512 | 128 | N-dominated | 1e-02 | adjoint | 640 | 2 | 3.40e-03 | 0.99 | 1024 | 1 | 8.45e-03 | 1.37 | 1.39x |
| float | 512 | 128 | N-dominated | 1e-02 | forward | 640 | 2 | 3.31e-04 | 0.86 | 1024 | 1 | 1.68e-03 | 1.30 | 1.51x |
| float | 512 | 128 | N-dominated | 1e-04 | adjoint | 640 | 4 | 1.62e-05 | 1.14 | 1024 | 3 | 2.27e-06 | 1.46 | 1.28x |
| float | 512 | 128 | N-dominated | 1e-04 | forward | 640 | 4 | 1.86e-06 | 0.86 | 1024 | 2 | 2.42e-05 | 1.34 | 1.56x |
| float | 512 | 512 | balanced | 1e-02 | adjoint | 640 | 2 | 1.66e-03 | 1.67 | 1024 | 1 | 4.58e-03 | 1.84 | 1.10x |
| float | 512 | 512 | balanced | 1e-02 | forward | 640 | 2 | 4.02e-04 | 1.37 | 1024 | 1 | 2.13e-03 | 1.65 | 1.20x |
| float | 512 | 512 | balanced | 1e-04 | adjoint | 640 | 4 | 6.65e-06 | 2.06 | 1024 | 2 | 6.27e-05 | 2.24 | 1.09x |
| float | 512 | 512 | balanced | 1e-04 | forward | 640 | 4 | 1.81e-06 | 1.48 | 1024 | 2 | 2.87e-05 | 1.64 | 1.11x |
| float | 512 | 2048 | M-dominated | 1e-02 | adjoint | 640 | 2 | 9.04e-04 | 4.55 | 1024 | 1 | 2.20e-03 | 4.27 | 0.94x |
| float | 512 | 2048 | M-dominated | 1e-02 | forward | 640 | 2 | 4.32e-04 | 3.49 | 1024 | 1 | 2.13e-03 | 3.25 | 0.93x |
| float | 512 | 2048 | M-dominated | 1e-04 | adjoint | 750 | 3 | 4.02e-06 | 6.21 | 1024 | 2 | 4.28e-05 | 4.84 | 0.78x |
| float | 512 | 2048 | M-dominated | 1e-04 | forward | 720 | 3 | 6.41e-06 | 4.34 | 1024 | 2 | 3.24e-05 | 3.59 | 0.83x |
| double | 243 | 60 | N-dominated | 1e-04 | adjoint | 320 | 4 | 4.30e-06 | 0.93 | 512 | 3 | 1.15e-06 | 1.29 | 1.38x |
| double | 243 | 60 | N-dominated | 1e-04 | forward | 320 | 3 | 1.13e-05 | 0.75 | 512 | 2 | 3.22e-05 | 1.28 | 1.70x |
| double | 243 | 60 | N-dominated | 1e-08 | adjoint | 320 | 7 | 2.29e-10 | 1.18 | 512 | 5 | 1.24e-10 | 1.36 | 1.15x |
| double | 243 | 60 | N-dominated | 1e-08 | forward | 320 | 6 | 8.63e-10 | 0.82 | 512 | 4 | 2.84e-09 | 1.26 | 1.53x |
| double | 243 | 60 | N-dominated | 1e-11 | adjoint | 320 | 9 | 5.09e-13 | 1.07 | 512 | 6 | 1.30e-12 | 1.36 | 1.27x |
| double | 243 | 60 | N-dominated | 1e-11 | forward | 320 | 9 | 6.28e-14 | 0.90 | 512 | 6 | 4.33e-13 | 1.30 | 1.45x |
| double | 243 | 243 | balanced | 1e-04 | adjoint | 320 | 4 | 2.86e-06 | 1.78 | 512 | 3 | 1.05e-06 | 1.96 | 1.10x |
| double | 243 | 243 | balanced | 1e-04 | forward | 320 | 3 | 1.38e-05 | 1.04 | 512 | 2 | 3.12e-05 | 1.26 | 1.21x |
| double | 243 | 243 | balanced | 1e-08 | adjoint | 320 | 7 | 1.33e-10 | 2.36 | 512 | 5 | 1.07e-10 | 2.23 | 0.94x |
| double | 243 | 243 | balanced | 1e-08 | forward | 320 | 6 | 9.60e-10 | 1.37 | 512 | 4 | 3.65e-09 | 1.52 | 1.11x |
| double | 243 | 243 | balanced | 1e-11 | adjoint | 384 | 7 | 1.12e-12 | 2.45 | 512 | 6 | 1.05e-12 | 2.44 | 1.00x |
| double | 243 | 243 | balanced | 1e-11 | forward | 384 | 7 | 3.63e-13 | 1.54 | 512 | 6 | 3.70e-13 | 1.69 | 1.09x |
| double | 243 | 972 | M-dominated | 1e-04 | adjoint | 360 | 3 | 4.51e-06 | 4.71 | 512 | 2 | 4.25e-05 | 4.11 | 0.87x |
| double | 243 | 972 | M-dominated | 1e-04 | forward | 320 | 3 | 1.71e-05 | 2.41 | 512 | 2 | 3.95e-05 | 2.46 | 1.02x |
| double | 243 | 972 | M-dominated | 1e-08 | adjoint | 400 | 5 | 8.18e-10 | 6.18 | 512 | 4 | 3.85e-09 | 5.74 | 0.93x |
| double | 243 | 972 | M-dominated | 1e-08 | forward | 384 | 5 | 9.95e-10 | 3.67 | 512 | 4 | 3.39e-09 | 3.18 | 0.87x |
| double | 243 | 972 | M-dominated | 1e-11 | adjoint | 486 | 6 | 9.74e-13 | 7.18 | 512 | 6 | 5.33e-13 | 7.51 | 1.05x |
| double | 243 | 972 | M-dominated | 1e-11 | forward | 480 | 6 | 7.30e-13 | 4.40 | 512 | 6 | 3.76e-13 | 4.31 | 0.98x |
| double | 250 | 62 | N-dominated | 1e-04 | adjoint | 320 | 4 | 1.61e-05 | 0.93 | 512 | 3 | 2.63e-06 | 1.25 | 1.35x |
| double | 250 | 62 | N-dominated | 1e-04 | forward | 324 | 3 | 1.51e-05 | 0.80 | 512 | 2 | 3.08e-05 | 1.21 | 1.51x |
| double | 250 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 3.33e-09 | 0.97 | 512 | 5 | 4.40e-10 | 1.35 | 1.40x |
| double | 250 | 62 | N-dominated | 1e-08 | forward | 320 | 7 | 1.24e-10 | 0.78 | 512 | 4 | 3.12e-09 | 1.25 | 1.61x |
| double | 250 | 62 | N-dominated | 1e-11 | adjoint | 324 | 9 | 4.86e-12 | 1.23 | 512 | 6 | 5.76e-12 | 1.46 | 1.19x |
| double | 250 | 62 | N-dominated | 1e-11 | forward | 320 | 9 | 4.38e-13 | 0.89 | 512 | 6 | 4.76e-13 | 1.31 | 1.47x |
| double | 250 | 250 | balanced | 1e-04 | adjoint | 360 | 3 | 1.18e-05 | 1.97 | 512 | 2 | 7.73e-05 | 1.71 | 0.87x |
| double | 250 | 250 | balanced | 1e-04 | forward | 324 | 3 | 1.67e-05 | 1.17 | 512 | 2 | 3.88e-05 | 1.37 | 1.17x |
| double | 250 | 250 | balanced | 1e-08 | adjoint | 360 | 6 | 2.29e-10 | 2.45 | 512 | 5 | 1.27e-10 | 2.37 | 0.97x |
| double | 250 | 250 | balanced | 1e-08 | forward | 384 | 5 | 1.54e-09 | 1.49 | 512 | 4 | 3.84e-09 | 1.58 | 1.06x |
| double | 250 | 250 | balanced | 1e-11 | adjoint | 360 | 8 | 2.25e-13 | 2.95 | 512 | 6 | 1.70e-12 | 2.61 | 0.88x |
| double | 250 | 250 | balanced | 1e-11 | forward | 384 | 7 | 9.34e-13 | 1.67 | 512 | 6 | 4.63e-13 | 1.79 | 1.07x |
| double | 250 | 1000 | M-dominated | 1e-04 | adjoint | 360 | 3 | 5.20e-06 | 5.00 | 512 | 2 | 5.26e-05 | 3.97 | 0.79x |
| double | 250 | 1000 | M-dominated | 1e-04 | forward | 324 | 3 | 2.13e-05 | 2.60 | 512 | 2 | 3.83e-05 | 2.32 | 0.89x |
| double | 250 | 1000 | M-dominated | 1e-08 | adjoint | 432 | 5 | 3.86e-10 | 6.32 | 512 | 4 | 5.02e-09 | 5.87 | 0.93x |
| double | 250 | 1000 | M-dominated | 1e-08 | forward | 384 | 5 | 1.75e-09 | 3.45 | 512 | 4 | 5.57e-09 | 3.44 | 1.00x |
| double | 250 | 1000 | M-dominated | 1e-11 | adjoint | 500 | 6 | 1.23e-12 | 7.38 | 512 | 6 | 7.71e-13 | 7.20 | 0.98x |
| double | 250 | 1000 | M-dominated | 1e-11 | forward | 480 | 6 | 1.26e-12 | 4.54 | 512 | 6 | 5.19e-13 | 4.29 | 0.95x |
| double | 251 | 62 | N-dominated | 1e-04 | adjoint | 320 | 4 | 5.95e-06 | 0.86 | 512 | 3 | 1.45e-06 | 1.39 | 1.60x |
| double | 251 | 62 | N-dominated | 1e-04 | forward | 320 | 4 | 8.89e-07 | 0.79 | 512 | 2 | 3.10e-05 | 1.15 | 1.46x |
| double | 251 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 6.15e-10 | 1.10 | 512 | 5 | 1.58e-10 | 1.36 | 1.24x |
| double | 251 | 62 | N-dominated | 1e-08 | forward | 320 | 7 | 8.76e-11 | 0.88 | 512 | 4 | 3.42e-09 | 1.19 | 1.35x |
| double | 251 | 62 | N-dominated | 1e-11 | adjoint | 324 | 10 | 1.32e-13 | 1.28 | 512 | 6 | 1.41e-12 | 1.39 | 1.08x |
| double | 251 | 62 | N-dominated | 1e-11 | forward | 320 | 9 | 2.24e-13 | 0.94 | 512 | 6 | 3.99e-13 | 1.30 | 1.37x |
| double | 251 | 251 | balanced | 1e-04 | adjoint | 360 | 3 | 1.73e-05 | 1.91 | 512 | 2 | 7.46e-05 | 1.81 | 0.94x |
| double | 251 | 251 | balanced | 1e-04 | forward | 360 | 3 | 5.07e-06 | 1.43 | 512 | 2 | 3.60e-05 | 1.40 | 0.97x |
| double | 251 | 251 | balanced | 1e-08 | adjoint | 360 | 6 | 3.10e-10 | 2.46 | 512 | 4 | 8.46e-09 | 1.99 | 0.81x |
| double | 251 | 251 | balanced | 1e-08 | forward | 400 | 5 | 7.09e-10 | 1.46 | 512 | 4 | 5.49e-09 | 1.46 | 1.00x |
| double | 251 | 251 | balanced | 1e-11 | adjoint | 400 | 7 | 8.81e-13 | 2.67 | 512 | 6 | 9.73e-13 | 2.56 | 0.96x |
| double | 251 | 251 | balanced | 1e-11 | forward | 384 | 7 | 6.60e-13 | 1.53 | 512 | 6 | 4.63e-13 | 1.66 | 1.09x |
| double | 251 | 1004 | M-dominated | 1e-04 | adjoint | 360 | 3 | 8.73e-06 | 4.80 | 512 | 2 | 4.99e-05 | 4.05 | 0.84x |
| double | 251 | 1004 | M-dominated | 1e-04 | forward | 360 | 3 | 5.66e-06 | 2.61 | 512 | 2 | 3.88e-05 | 2.31 | 0.88x |
| double | 251 | 1004 | M-dominated | 1e-08 | adjoint | 432 | 5 | 3.35e-10 | 6.52 | 512 | 4 | 4.42e-09 | 6.38 | 0.98x |
| double | 251 | 1004 | M-dominated | 1e-08 | forward | 400 | 5 | 7.69e-10 | 3.73 | 512 | 4 | 3.98e-09 | 3.23 | 0.87x |
| double | 251 | 1004 | M-dominated | 1e-11 | adjoint | 500 | 6 | 9.11e-13 | 7.55 | 512 | 6 | 4.47e-13 | 7.03 | 0.93x |
| double | 251 | 1004 | M-dominated | 1e-11 | forward | 480 | 6 | 8.89e-13 | 4.42 | 512 | 6 | 4.36e-13 | 4.13 | 0.94x |
| double | 255 | 63 | N-dominated | 1e-04 | adjoint | 320 | 4 | 7.87e-06 | 1.03 | 512 | 3 | 1.46e-06 | 1.17 | 1.14x |
| double | 255 | 63 | N-dominated | 1e-04 | forward | 320 | 4 | 1.17e-06 | 0.83 | 512 | 2 | 3.45e-05 | 1.15 | 1.38x |
| double | 255 | 63 | N-dominated | 1e-08 | adjoint | 320 | 7 | 8.64e-10 | 1.11 | 512 | 5 | 2.36e-10 | 1.36 | 1.23x |
| double | 255 | 63 | N-dominated | 1e-08 | forward | 320 | 7 | 1.88e-10 | 0.88 | 512 | 4 | 3.71e-09 | 1.26 | 1.43x |
| double | 255 | 63 | N-dominated | 1e-11 | adjoint | 360 | 8 | 1.18e-12 | 1.27 | 512 | 6 | 2.72e-12 | 1.31 | 1.04x |
| double | 255 | 63 | N-dominated | 1e-11 | forward | 324 | 9 | 3.46e-13 | 1.02 | 512 | 6 | 4.83e-13 | 1.32 | 1.30x |
| double | 255 | 255 | balanced | 1e-04 | adjoint | 360 | 3 | 1.68e-05 | 1.92 | 512 | 2 | 7.59e-05 | 2.48 | 1.29x |
| double | 255 | 255 | balanced | 1e-04 | forward | 360 | 3 | 6.00e-06 | 1.40 | 512 | 2 | 3.63e-05 | 1.38 | 0.99x |
| double | 255 | 255 | balanced | 1e-08 | adjoint | 360 | 6 | 4.95e-10 | 2.27 | 512 | 4 | 7.42e-09 | 1.99 | 0.88x |
| double | 255 | 255 | balanced | 1e-08 | forward | 400 | 5 | 1.03e-09 | 1.60 | 512 | 4 | 4.22e-09 | 1.38 | 0.86x |
| double | 255 | 255 | balanced | 1e-11 | adjoint | 360 | 8 | 8.52e-13 | 2.59 | 512 | 6 | 8.22e-13 | 2.35 | 0.91x |
| double | 255 | 255 | balanced | 1e-11 | forward | 384 | 7 | 1.11e-12 | 2.03 | 512 | 6 | 4.95e-13 | 1.59 | 0.79x |
| double | 255 | 1020 | M-dominated | 1e-04 | adjoint | 360 | 3 | 1.05e-05 | 4.86 | 512 | 2 | 5.99e-05 | 4.08 | 0.84x |
| double | 255 | 1020 | M-dominated | 1e-04 | forward | 360 | 3 | 7.74e-06 | 2.91 | 512 | 2 | 4.10e-05 | 2.42 | 0.83x |
| double | 255 | 1020 | M-dominated | 1e-08 | adjoint | 432 | 5 | 4.01e-10 | 6.59 | 512 | 4 | 5.71e-09 | 6.17 | 0.94x |
| double | 255 | 1020 | M-dominated | 1e-08 | forward | 400 | 5 | 1.11e-09 | 3.98 | 512 | 4 | 5.29e-09 | 3.17 | 0.80x |
| double | 255 | 1020 | M-dominated | 1e-11 | adjoint | 512 | 6 | 5.97e-13 | 7.52 | 512 | 6 | 5.97e-13 | 7.12 | 0.95x |
| double | 255 | 1020 | M-dominated | 1e-11 | forward | 480 | 6 | 1.49e-12 | 4.39 | 512 | 6 | 7.14e-13 | 4.34 | 0.99x |
| double | 256 | 64 | N-dominated | 1e-04 | adjoint | 320 | 4 | 1.41e-05 | 0.89 | 512 | 3 | 2.10e-06 | 1.25 | 1.41x |
| double | 256 | 64 | N-dominated | 1e-04 | forward | 320 | 4 | 1.37e-06 | 0.79 | 512 | 2 | 3.23e-05 | 0.99 | 1.27x |
| double | 256 | 64 | N-dominated | 1e-08 | adjoint | 320 | 7 | 4.00e-09 | 1.03 | 512 | 5 | 1.95e-10 | 1.34 | 1.30x |
| double | 256 | 64 | N-dominated | 1e-08 | forward | 320 | 7 | 2.87e-10 | 0.79 | 512 | 4 | 3.99e-09 | 1.05 | 1.34x |
| double | 256 | 64 | N-dominated | 1e-11 | adjoint | 360 | 8 | 3.69e-12 | 1.30 | 512 | 6 | 2.62e-12 | 1.40 | 1.07x |
| double | 256 | 64 | N-dominated | 1e-11 | forward | 324 | 9 | 5.94e-13 | 1.05 | 512 | 6 | 5.27e-13 | 1.66 | 1.57x |
| double | 256 | 256 | balanced | 1e-04 | adjoint | 360 | 3 | 3.14e-05 | 1.86 | 512 | 2 | 9.43e-05 | 1.63 | 0.88x |
| double | 256 | 256 | balanced | 1e-04 | forward | 360 | 3 | 6.06e-06 | 1.32 | 512 | 2 | 3.49e-05 | 1.24 | 0.94x |
| double | 256 | 256 | balanced | 1e-08 | adjoint | 360 | 6 | 9.53e-10 | 2.21 | 512 | 5 | 2.10e-10 | 2.19 | 0.99x |
| double | 256 | 256 | balanced | 1e-08 | forward | 400 | 5 | 1.05e-09 | 1.39 | 512 | 4 | 4.30e-09 | 1.41 | 1.01x |
| double | 256 | 256 | balanced | 1e-11 | adjoint | 360 | 8 | 8.31e-13 | 2.58 | 512 | 6 | 2.96e-12 | 2.36 | 0.91x |
| double | 256 | 256 | balanced | 1e-11 | forward | 400 | 7 | 4.99e-13 | 1.50 | 512 | 6 | 5.83e-13 | 1.60 | 1.07x |
| double | 256 | 1024 | M-dominated | 1e-04 | adjoint | 360 | 3 | 1.47e-05 | 5.07 | 512 | 2 | 3.93e-05 | 4.17 | 0.82x |
| double | 256 | 1024 | M-dominated | 1e-04 | forward | 360 | 3 | 7.27e-06 | 2.70 | 512 | 2 | 3.86e-05 | 2.40 | 0.89x |
| double | 256 | 1024 | M-dominated | 1e-08 | adjoint | 432 | 5 | 6.32e-10 | 6.53 | 512 | 4 | 4.28e-09 | 5.83 | 0.89x |
| double | 256 | 1024 | M-dominated | 1e-08 | forward | 400 | 5 | 1.35e-09 | 3.75 | 512 | 4 | 4.85e-09 | 3.24 | 0.86x |
| double | 256 | 1024 | M-dominated | 1e-11 | adjoint | 512 | 6 | 5.31e-13 | 6.97 | 512 | 6 | 5.31e-13 | 6.79 | 0.97x |
| double | 256 | 1024 | M-dominated | 1e-11 | forward | 480 | 6 | 1.67e-12 | 4.54 | 512 | 6 | 6.43e-13 | 4.36 | 0.96x |
| double | 512 | 128 | N-dominated | 1e-04 | adjoint | 640 | 4 | 1.63e-05 | 1.89 | 1024 | 3 | 2.01e-06 | 2.88 | 1.53x |
| double | 512 | 128 | N-dominated | 1e-04 | forward | 640 | 4 | 1.13e-06 | 1.80 | 1024 | 2 | 2.42e-05 | 2.55 | 1.42x |
| double | 512 | 128 | N-dominated | 1e-08 | adjoint | 640 | 7 | 4.90e-09 | 2.16 | 1024 | 5 | 3.69e-10 | 2.96 | 1.37x |
| double | 512 | 128 | N-dominated | 1e-08 | forward | 640 | 7 | 2.33e-10 | 2.06 | 1024 | 4 | 3.17e-09 | 2.96 | 1.44x |
| double | 512 | 128 | N-dominated | 1e-11 | adjoint | 720 | 8 | 2.84e-12 | 2.67 | 1024 | 6 | 5.02e-12 | 3.23 | 1.21x |
| double | 512 | 128 | N-dominated | 1e-11 | forward | 648 | 9 | 3.51e-13 | 2.28 | 1024 | 6 | 3.78e-13 | 2.79 | 1.23x |
| double | 512 | 512 | balanced | 1e-04 | adjoint | 640 | 4 | 6.36e-06 | 3.74 | 1024 | 2 | 6.25e-05 | 3.85 | 1.03x |
| double | 512 | 512 | balanced | 1e-04 | forward | 640 | 4 | 1.19e-06 | 2.44 | 1024 | 2 | 2.87e-05 | 3.01 | 1.23x |
| double | 512 | 512 | balanced | 1e-08 | adjoint | 640 | 7 | 1.91e-09 | 4.97 | 1024 | 4 | 6.93e-09 | 4.67 | 0.94x |
| double | 512 | 512 | balanced | 1e-08 | forward | 640 | 7 | 2.05e-10 | 3.26 | 1024 | 4 | 3.73e-09 | 3.45 | 1.06x |
| double | 512 | 512 | balanced | 1e-11 | adjoint | 720 | 8 | 8.28e-13 | 5.70 | 1024 | 6 | 1.01e-12 | 5.55 | 0.97x |
| double | 512 | 512 | balanced | 1e-11 | forward | 720 | 8 | 1.65e-13 | 3.85 | 1024 | 6 | 5.31e-13 | 3.99 | 1.04x |
| double | 512 | 2048 | M-dominated | 1e-04 | adjoint | 720 | 3 | 7.36e-06 | 9.58 | 1024 | 2 | 4.31e-05 | 8.95 | 0.93x |
| double | 512 | 2048 | M-dominated | 1e-04 | forward | 720 | 3 | 5.57e-06 | 5.70 | 1024 | 2 | 3.19e-05 | 5.35 | 0.94x |
| double | 512 | 2048 | M-dominated | 1e-08 | adjoint | 864 | 5 | 7.38e-10 | 13.39 | 1024 | 4 | 6.71e-09 | 11.32 | 0.85x |
| double | 512 | 2048 | M-dominated | 1e-08 | forward | 800 | 5 | 8.06e-10 | 8.01 | 1024 | 4 | 4.14e-09 | 7.20 | 0.90x |
| double | 512 | 2048 | M-dominated | 1e-11 | adjoint | 1024 | 6 | 7.75e-13 | 14.95 | 1024 | 6 | 7.75e-13 | 14.89 | 1.00x |
| double | 512 | 2048 | M-dominated | 1e-11 | forward | 960 | 6 | 1.32e-12 | 9.13 | 1024 | 6 | 5.64e-13 | 9.05 | 0.99x |
| long double | 243 | 60 | N-dominated | 1e-08 | adjoint | 320 | 7 | 2.29e-10 | 89.53 | 512 | 5 | 1.24e-10 | 124.18 | 1.39x |
| long double | 243 | 60 | N-dominated | 1e-08 | forward | 320 | 6 | 8.63e-10 | 86.37 | 512 | 4 | 2.84e-09 | 122.08 | 1.41x |
| long double | 243 | 60 | N-dominated | 1e-14 | adjoint | 320 | 11 | 1.09e-15 | 103.21 | 512 | 8 | 1.48e-16 | 131.55 | 1.27x |
| long double | 243 | 60 | N-dominated | 1e-14 | forward | 320 | 11 | 1.41e-16 | 103.81 | 512 | 7 | 3.21e-15 | 130.37 | 1.26x |
| long double | 243 | 60 | N-dominated | 1e-20 | adjoint | 320 | 16 | 3.59e-22 | 121.46 | 512 | 11 | 1.90e-22 | 144.30 | 1.19x |
| long double | 243 | 60 | N-dominated | 1e-20 | forward | 320 | 15 | 8.35e-22 | 143.15 | 512 | 10 | 3.07e-21 | 149.39 | 1.04x |
| long double | 243 | 243 | balanced | 1e-08 | adjoint | 320 | 7 | 1.33e-10 | 190.56 | 512 | 5 | 1.07e-10 | 213.80 | 1.12x |
| long double | 243 | 243 | balanced | 1e-08 | forward | 320 | 6 | 9.60e-10 | 173.12 | 512 | 4 | 3.65e-09 | 175.13 | 1.01x |
| long double | 243 | 243 | balanced | 1e-14 | adjoint | 384 | 9 | 5.89e-16 | 249.27 | 512 | 7 | 7.20e-15 | 251.32 | 1.01x |
| long double | 243 | 243 | balanced | 1e-14 | forward | 384 | 9 | 1.37e-16 | 272.05 | 512 | 7 | 3.28e-15 | 226.39 | 0.83x |
| long double | 243 | 243 | balanced | 1e-20 | adjoint | 384 | 13 | 1.94e-22 | 321.81 | 512 | 11 | 1.12e-22 | 328.13 | 1.02x |
| long double | 243 | 243 | balanced | 1e-20 | forward | 360 | 13 | 6.10e-22 | 330.09 | 512 | 10 | 3.16e-21 | 287.06 | 0.87x |
| long double | 243 | 972 | M-dominated | 1e-08 | adjoint | 400 | 5 | 8.18e-10 | 524.61 | 512 | 4 | 3.85e-09 | 462.72 | 0.88x |
| long double | 243 | 972 | M-dominated | 1e-08 | forward | 384 | 5 | 9.95e-10 | 523.15 | 512 | 4 | 3.39e-09 | 476.06 | 0.91x |
| long double | 243 | 972 | M-dominated | 1e-14 | adjoint | 450 | 8 | 4.24e-16 | 838.99 | 512 | 7 | 5.32e-15 | 700.89 | 0.84x |
| long double | 243 | 972 | M-dominated | 1e-14 | forward | 576 | 7 | 8.15e-16 | 774.16 | 512 | 7 | 3.47e-15 | 714.15 | 0.92x |
| long double | 243 | 972 | M-dominated | 1e-20 | adjoint | 576 | 10 | 8.15e-22 | 939.85 | 512 | 10 | 4.07e-21 | 889.36 | 0.95x |
| long double | 243 | 972 | M-dominated | 1e-20 | forward | 540 | 10 | 1.21e-21 | 1039.78 | 512 | 10 | 3.77e-21 | 961.70 | 0.92x |
| long double | 250 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 3.33e-09 | 93.39 | 512 | 5 | 4.40e-10 | 128.69 | 1.38x |
| long double | 250 | 62 | N-dominated | 1e-08 | forward | 320 | 7 | 1.24e-10 | 96.49 | 512 | 4 | 3.12e-09 | 128.50 | 1.33x |
| long double | 250 | 62 | N-dominated | 1e-14 | adjoint | 320 | 12 | 1.73e-15 | 112.08 | 512 | 8 | 8.25e-16 | 145.91 | 1.30x |
| long double | 250 | 62 | N-dominated | 1e-14 | forward | 320 | 11 | 9.79e-16 | 111.34 | 512 | 7 | 4.64e-15 | 134.30 | 1.21x |
| long double | 250 | 62 | N-dominated | 1e-20 | adjoint | 324 | 16 | 4.81e-21 | 136.56 | 512 | 11 | 1.34e-21 | 160.33 | 1.17x |
| long double | 250 | 62 | N-dominated | 1e-20 | forward | 320 | 16 | 4.50e-22 | 133.75 | 512 | 10 | 7.38e-21 | 150.43 | 1.12x |
| long double | 250 | 250 | balanced | 1e-08 | adjoint | 360 | 6 | 2.29e-10 | 210.13 | 512 | 5 | 1.27e-10 | 213.48 | 1.02x |
| long double | 250 | 250 | balanced | 1e-08 | forward | 384 | 5 | 1.54e-09 | 181.45 | 512 | 4 | 3.84e-09 | 184.76 | 1.02x |
| long double | 250 | 250 | balanced | 1e-14 | adjoint | 360 | 10 | 2.20e-16 | 286.82 | 512 | 8 | 2.43e-16 | 244.41 | 0.85x |
| long double | 250 | 250 | balanced | 1e-14 | forward | 384 | 9 | 4.69e-16 | 269.51 | 512 | 7 | 4.82e-15 | 229.22 | 0.85x |
| long double | 250 | 250 | balanced | 1e-20 | adjoint | 384 | 13 | 6.87e-22 | 326.01 | 512 | 11 | 3.93e-22 | 349.99 | 1.07x |
| long double | 250 | 250 | balanced | 1e-20 | forward | 384 | 13 | 1.59e-22 | 409.14 | 512 | 10 | 6.12e-21 | 279.98 | 0.68x |
| long double | 250 | 1000 | M-dominated | 1e-08 | adjoint | 432 | 5 | 3.86e-10 | 565.05 | 512 | 4 | 5.02e-09 | 473.85 | 0.84x |
| long double | 250 | 1000 | M-dominated | 1e-08 | forward | 384 | 5 | 1.75e-09 | 535.46 | 512 | 4 | 5.57e-09 | 478.55 | 0.89x |
| long double | 250 | 1000 | M-dominated | 1e-14 | adjoint | 480 | 8 | 2.86e-16 | 790.66 | 512 | 7 | 9.32e-15 | 742.02 | 0.94x |
| long double | 250 | 1000 | M-dominated | 1e-14 | forward | 576 | 7 | 1.41e-15 | 825.65 | 512 | 7 | 5.79e-15 | 712.98 | 0.86x |
| long double | 250 | 1000 | M-dominated | 1e-20 | adjoint | 600 | 10 | 3.42e-22 | 1004.42 | 512 | 11 | 1.80e-22 | 1005.40 | 1.00x |
| long double | 250 | 1000 | M-dominated | 1e-20 | forward | 576 | 10 | 9.46e-22 | 1049.73 | 512 | 10 | 9.36e-21 | 1009.20 | 0.96x |
| long double | 251 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 6.15e-10 | 92.46 | 512 | 5 | 1.58e-10 | 125.77 | 1.36x |
| long double | 251 | 62 | N-dominated | 1e-08 | forward | 320 | 7 | 8.76e-11 | 93.60 | 512 | 4 | 3.42e-09 | 128.54 | 1.37x |
| long double | 251 | 62 | N-dominated | 1e-14 | adjoint | 320 | 12 | 3.67e-16 | 111.08 | 512 | 8 | 2.32e-16 | 137.11 | 1.23x |
| long double | 251 | 62 | N-dominated | 1e-14 | forward | 320 | 11 | 6.10e-16 | 110.53 | 512 | 7 | 4.68e-15 | 139.47 | 1.26x |
| long double | 251 | 62 | N-dominated | 1e-20 | adjoint | 324 | 16 | 1.44e-21 | 141.40 | 512 | 11 | 2.85e-22 | 149.56 | 1.06x |
| long double | 251 | 62 | N-dominated | 1e-20 | forward | 320 | 16 | 3.49e-22 | 131.80 | 512 | 10 | 5.17e-21 | 146.49 | 1.11x |
| long double | 251 | 251 | balanced | 1e-08 | adjoint | 360 | 6 | 3.10e-10 | 209.67 | 512 | 4 | 8.46e-09 | 207.77 | 0.99x |
| long double | 251 | 251 | balanced | 1e-08 | forward | 400 | 5 | 7.09e-10 | 197.79 | 512 | 4 | 5.49e-09 | 183.80 | 0.93x |
| long double | 251 | 251 | balanced | 1e-14 | adjoint | 400 | 9 | 4.66e-16 | 257.74 | 512 | 7 | 7.51e-15 | 277.82 | 1.08x |
| long double | 251 | 251 | balanced | 1e-14 | forward | 384 | 9 | 3.44e-16 | 246.10 | 512 | 7 | 4.78e-15 | 238.85 | 0.97x |
| long double | 251 | 251 | balanced | 1e-20 | adjoint | 384 | 13 | 6.71e-22 | 307.75 | 512 | 11 | 1.48e-22 | 294.04 | 0.96x |
| long double | 251 | 251 | balanced | 1e-20 | forward | 384 | 13 | 1.34e-22 | 368.16 | 512 | 10 | 5.38e-21 | 336.20 | 0.91x |
| long double | 251 | 1004 | M-dominated | 1e-08 | adjoint | 432 | 5 | 3.35e-10 | 572.61 | 512 | 4 | 4.42e-09 | 467.76 | 0.82x |
| long double | 251 | 1004 | M-dominated | 1e-08 | forward | 400 | 5 | 7.69e-10 | 558.33 | 512 | 4 | 3.98e-09 | 498.98 | 0.89x |
| long double | 251 | 1004 | M-dominated | 1e-14 | adjoint | 480 | 8 | 1.95e-16 | 782.68 | 512 | 7 | 4.88e-15 | 672.68 | 0.86x |
| long double | 251 | 1004 | M-dominated | 1e-14 | forward | 576 | 7 | 1.15e-15 | 797.80 | 512 | 7 | 5.38e-15 | 789.93 | 0.99x |
| long double | 251 | 1004 | M-dominated | 1e-20 | adjoint | 600 | 10 | 6.51e-22 | 1008.49 | 512 | 11 | 9.87e-23 | 1043.54 | 1.03x |
| long double | 251 | 1004 | M-dominated | 1e-20 | forward | 576 | 10 | 7.27e-22 | 1041.62 | 512 | 10 | 7.31e-21 | 1044.64 | 1.00x |
| long double | 255 | 63 | N-dominated | 1e-08 | adjoint | 320 | 7 | 8.64e-10 | 92.08 | 512 | 5 | 2.36e-10 | 125.22 | 1.36x |
| long double | 255 | 63 | N-dominated | 1e-08 | forward | 320 | 7 | 1.88e-10 | 93.06 | 512 | 4 | 3.71e-09 | 134.68 | 1.45x |
| long double | 255 | 63 | N-dominated | 1e-14 | adjoint | 320 | 12 | 8.43e-16 | 109.49 | 512 | 8 | 4.52e-16 | 133.75 | 1.22x |
| long double | 255 | 63 | N-dominated | 1e-14 | forward | 320 | 12 | 9.96e-17 | 114.01 | 512 | 7 | 4.88e-15 | 134.04 | 1.18x |
| long double | 255 | 63 | N-dominated | 1e-20 | adjoint | 320 | 17 | 6.20e-22 | 128.98 | 512 | 11 | 4.62e-22 | 147.32 | 1.14x |
| long double | 255 | 63 | N-dominated | 1e-20 | forward | 324 | 16 | 5.01e-22 | 141.73 | 512 | 10 | 7.43e-21 | 144.61 | 1.02x |
| long double | 255 | 255 | balanced | 1e-08 | adjoint | 360 | 6 | 4.95e-10 | 212.82 | 512 | 4 | 7.42e-09 | 183.06 | 0.86x |
| long double | 255 | 255 | balanced | 1e-08 | forward | 400 | 5 | 1.03e-09 | 199.94 | 512 | 4 | 4.22e-09 | 201.28 | 1.01x |
| long double | 255 | 255 | balanced | 1e-14 | adjoint | 400 | 9 | 8.14e-16 | 254.68 | 512 | 8 | 1.20e-16 | 284.02 | 1.12x |
| long double | 255 | 255 | balanced | 1e-14 | forward | 384 | 9 | 5.75e-16 | 281.13 | 512 | 7 | 6.26e-15 | 248.07 | 0.88x |
| long double | 255 | 255 | balanced | 1e-20 | adjoint | 432 | 12 | 3.53e-22 | 338.48 | 512 | 11 | 1.80e-22 | 307.31 | 0.91x |
| long double | 255 | 255 | balanced | 1e-20 | forward | 384 | 13 | 2.47e-22 | 335.81 | 512 | 10 | 7.00e-21 | 356.42 | 1.06x |
| long double | 255 | 1020 | M-dominated | 1e-08 | adjoint | 432 | 5 | 4.01e-10 | 615.72 | 512 | 4 | 5.71e-09 | 478.21 | 0.78x |
| long double | 255 | 1020 | M-dominated | 1e-08 | forward | 400 | 5 | 1.11e-09 | 561.11 | 512 | 4 | 5.29e-09 | 550.79 | 0.98x |
| long double | 255 | 1020 | M-dominated | 1e-14 | adjoint | 480 | 8 | 2.26e-16 | 776.41 | 512 | 7 | 7.22e-15 | 700.19 | 0.90x |
| long double | 255 | 1020 | M-dominated | 1e-14 | forward | 576 | 7 | 1.55e-15 | 805.96 | 512 | 7 | 6.52e-15 | 744.21 | 0.92x |
| long double | 255 | 1020 | M-dominated | 1e-20 | adjoint | 600 | 10 | 4.94e-22 | 1032.03 | 512 | 11 | 1.52e-22 | 1048.86 | 1.02x |
| long double | 255 | 1020 | M-dominated | 1e-20 | forward | 576 | 10 | 9.49e-22 | 1088.25 | 512 | 10 | 9.17e-21 | 1048.76 | 0.96x |
| long double | 256 | 64 | N-dominated | 1e-08 | adjoint | 320 | 7 | 4.00e-09 | 93.86 | 512 | 5 | 1.95e-10 | 127.52 | 1.36x |
| long double | 256 | 64 | N-dominated | 1e-08 | forward | 320 | 7 | 2.87e-10 | 90.62 | 512 | 4 | 3.99e-09 | 122.33 | 1.35x |
| long double | 256 | 64 | N-dominated | 1e-14 | adjoint | 320 | 12 | 3.73e-15 | 112.46 | 512 | 8 | 4.08e-16 | 141.09 | 1.25x |
| long double | 256 | 64 | N-dominated | 1e-14 | forward | 320 | 12 | 2.95e-16 | 113.26 | 512 | 7 | 6.97e-15 | 136.49 | 1.21x |
| long double | 256 | 64 | N-dominated | 1e-20 | adjoint | 320 | 17 | 3.65e-21 | 131.59 | 512 | 11 | 7.69e-22 | 156.21 | 1.19x |
| long double | 256 | 64 | N-dominated | 1e-20 | forward | 324 | 16 | 1.31e-21 | 149.49 | 512 | 11 | 1.51e-22 | 163.60 | 1.09x |
| long double | 256 | 256 | balanced | 1e-08 | adjoint | 360 | 6 | 9.53e-10 | 216.64 | 512 | 5 | 2.10e-10 | 218.01 | 1.01x |
| long double | 256 | 256 | balanced | 1e-08 | forward | 400 | 5 | 1.05e-09 | 184.94 | 512 | 4 | 4.30e-09 | 221.10 | 1.20x |
| long double | 256 | 256 | balanced | 1e-14 | adjoint | 400 | 9 | 9.32e-16 | 253.72 | 512 | 8 | 4.68e-16 | 281.40 | 1.11x |
| long double | 256 | 256 | balanced | 1e-14 | forward | 384 | 9 | 9.01e-16 | 290.69 | 512 | 7 | 7.25e-15 | 249.09 | 0.86x |
| long double | 256 | 256 | balanced | 1e-20 | adjoint | 432 | 12 | 1.02e-21 | 324.42 | 512 | 11 | 8.82e-22 | 298.85 | 0.92x |
| long double | 256 | 256 | balanced | 1e-20 | forward | 384 | 13 | 4.51e-22 | 316.65 | 512 | 11 | 1.31e-22 | 374.37 | 1.18x |
| long double | 256 | 1024 | M-dominated | 1e-08 | adjoint | 432 | 5 | 6.32e-10 | 632.44 | 512 | 4 | 4.28e-09 | 485.62 | 0.77x |
| long double | 256 | 1024 | M-dominated | 1e-08 | forward | 400 | 5 | 1.35e-09 | 555.56 | 512 | 4 | 4.85e-09 | 565.15 | 1.02x |
| long double | 256 | 1024 | M-dominated | 1e-14 | adjoint | 480 | 8 | 9.57e-16 | 783.07 | 512 | 7 | 7.31e-15 | 696.56 | 0.89x |
| long double | 256 | 1024 | M-dominated | 1e-14 | forward | 600 | 7 | 1.30e-15 | 882.66 | 512 | 7 | 8.76e-15 | 733.68 | 0.83x |
| long double | 256 | 1024 | M-dominated | 1e-20 | adjoint | 600 | 10 | 1.48e-21 | 1048.09 | 512 | 11 | 1.61e-22 | 1023.88 | 0.98x |
| long double | 256 | 1024 | M-dominated | 1e-20 | forward | 576 | 10 | 1.46e-21 | 1067.87 | 512 | 11 | 1.55e-22 | 1152.67 | 1.08x |
| long double | 512 | 128 | N-dominated | 1e-08 | adjoint | 640 | 7 | 4.90e-09 | 253.86 | 1024 | 5 | 3.69e-10 | 376.02 | 1.48x |
| long double | 512 | 128 | N-dominated | 1e-08 | forward | 640 | 7 | 2.33e-10 | 268.08 | 1024 | 4 | 3.17e-09 | 425.69 | 1.59x |
| long double | 512 | 128 | N-dominated | 1e-14 | adjoint | 640 | 12 | 4.92e-15 | 296.41 | 1024 | 8 | 7.98e-16 | 391.71 | 1.32x |
| long double | 512 | 128 | N-dominated | 1e-14 | forward | 640 | 12 | 1.57e-16 | 308.69 | 1024 | 7 | 5.62e-15 | 358.66 | 1.16x |
| long double | 512 | 128 | N-dominated | 1e-20 | adjoint | 640 | 17 | 4.90e-21 | 299.32 | 1024 | 11 | 1.55e-21 | 382.24 | 1.28x |
| long double | 512 | 128 | N-dominated | 1e-20 | forward | 648 | 16 | 8.30e-22 | 429.78 | 1024 | 10 | 7.38e-21 | 509.84 | 1.19x |
| long double | 512 | 512 | balanced | 1e-08 | adjoint | 640 | 7 | 1.91e-09 | 460.57 | 1024 | 4 | 6.93e-09 | 516.91 | 1.12x |
| long double | 512 | 512 | balanced | 1e-08 | forward | 640 | 7 | 2.05e-10 | 489.83 | 1024 | 4 | 3.73e-09 | 614.01 | 1.25x |
| long double | 512 | 512 | balanced | 1e-14 | adjoint | 720 | 10 | 1.07e-15 | 616.36 | 1024 | 7 | 8.35e-15 | 621.11 | 1.01x |
| long double | 512 | 512 | balanced | 1e-14 | forward | 768 | 9 | 8.05e-16 | 615.19 | 1024 | 7 | 4.99e-15 | 643.29 | 1.05x |
| long double | 512 | 512 | balanced | 1e-20 | adjoint | 864 | 12 | 3.52e-22 | 795.85 | 1024 | 11 | 1.82e-22 | 764.41 | 0.96x |
| long double | 512 | 512 | balanced | 1e-20 | forward | 768 | 13 | 3.31e-22 | 797.56 | 1024 | 10 | 9.06e-21 | 776.25 | 0.97x |
| long double | 512 | 2048 | M-dominated | 1e-08 | adjoint | 864 | 5 | 7.38e-10 | 1330.82 | 1024 | 4 | 6.71e-09 | 1205.39 | 0.91x |
| long double | 512 | 2048 | M-dominated | 1e-08 | forward | 800 | 5 | 8.06e-10 | 1391.34 | 1024 | 4 | 4.14e-09 | 1332.86 | 0.96x |
| long double | 512 | 2048 | M-dominated | 1e-14 | adjoint | 960 | 8 | 3.06e-16 | 1773.25 | 1024 | 7 | 7.48e-15 | 1818.85 | 1.03x |
| long double | 512 | 2048 | M-dominated | 1e-14 | forward | 900 | 8 | 8.02e-16 | 2035.88 | 1024 | 7 | 6.18e-15 | 1770.34 | 0.87x |
| long double | 512 | 2048 | M-dominated | 1e-20 | adjoint | 1200 | 10 | 1.05e-21 | 2380.32 | 1024 | 11 | 1.33e-22 | 2333.36 | 0.98x |
| long double | 512 | 2048 | M-dominated | 1e-20 | forward | 1152 | 10 | 9.91e-22 | 2441.12 | 1024 | 10 | 8.39e-21 | 2288.35 | 0.94x |

`!` marks a goal that was not met.

