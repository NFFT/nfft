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
- New is faster in **140/288** (49 %) of them.
- Median speedup **0.97x**.

## Verdict

**Accuracy: yes, consistently.** The tuner meets the goal in every one
of the 288 configurations, unaided. The legacy geometry also meets it in
every case, but only because it was handed the oracle cut-off; the legacy
API ships no way to find that `m`, so in practice the choice is between a
tuner that lands the accuracy and a guess that might not.

**Speed: no, not consistently -- it depends on how many nodes there are.**
With few nodes the tuner wins clearly (1.31x median, 95/96 cases); with many nodes it loses (0.77x median, only 3/96 wins).

The reason is visible in the table below: the tuner drives `n` down to
about 0.62 of the legacy size, and pays for it with a cut-off 2 to 3
larger. The FFT costs `O(n log n)` and the node convolution `O(M m)`, so
shrinking `n` while growing `m` is a good trade exactly while the FFT
dominates. Once `M` outgrows `N` it is the wrong trade.

`tune_plan` minimises the oversampling first and the cut-off second,
which is a memory-first policy, not a time-first one. Making it
time-optimal means giving it `M` and minimising `c1 n log n + c2 M m`
over the pair instead of taking the smallest `n` that works.

## By problem shape

This is where the answer splits. `M` is the node count.

| shape | M | median speedup | new faster | median n_new / n_leg | median m_new - m_leg |
|---|---|---|---|---|---|
| N-dominated | N/4 | 1.31x | 95/96 | 0.62 | +2 |
| balanced | N | 0.96x | 42/96 | 0.62 | +3 |
| M-dominated | 4N | 0.77x | 3/96 | 0.62 | +3 |

## By precision

| precision | median speedup | new faster | goal met (new) | goal met (legacy) |
|---|---|---|---|---|
| float | 1.12x | 44/72 | 72/72 | 72/72 |
| double | 0.98x | 52/108 | 108/108 | 108/108 |
| long double | 0.90x | 44/108 | 108/108 | 108/108 |

## By bandwidth

| N | factors | median speedup | median n_new | n_leg | median m_new | median m_leg |
|---|---|---|---|---|---|---|
| 243 | 3·3·3·3·3 | 1.09x | 320 | 512 | 7 | 4 |
| 250 | 2·5·5·5 | 0.93x | 320 | 512 | 7 | 4 |
| 251 | 251 | 0.94x | 320 | 512 | 7 | 4 |
| 255 | 3·5·17 | 1.03x | 320 | 512 | 7 | 4 |
| 256 | 2·2·2·2·2·2·2·2 | 0.97x | 320 | 512 | 7 | 4 |
| 512 | 2·2·2·2·2·2·2·2·2 | 1.05x | 640 | 1024 | 7 | 4 |

## Accuracy

Both sides are asked for the same goal, so the question is not who is
more accurate but whether either misses. Cases where the goal was not
met:

- **New: none.**
- **Legacy (with oracle m): none.**

Median headroom below the goal: new 27.4x, legacy 3.1x.
Both overshoot the goal rather than miss it; the new side does so by
more, which is the cost of an upper-bound model against a measured search.

## Full results

| precision | N | M | shape | goal | dir | n_new | m_new | err_new | t_new (µs) | n_leg | m_leg | err_leg | t_leg (µs) | speedup |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| float | 243 | 60 | N-dominated | 1e-02 | adjoint | 320 | 2 | 2.13e-03 | 0.45 | 512 | 1 | 9.06e-03 | 0.56 | 1.25x |
| float | 243 | 60 | N-dominated | 1e-02 | forward | 320 | 2 | 2.87e-04 | 0.40 | 512 | 1 | 2.36e-03 | 0.58 | 1.44x |
| float | 243 | 60 | N-dominated | 1e-04 | adjoint | 320 | 4 | 4.78e-06 | 0.49 | 512 | 3 | 1.59e-06 | 0.61 | 1.23x |
| float | 243 | 60 | N-dominated | 1e-04 | forward | 320 | 3 | 1.15e-05 | 0.48 | 512 | 2 | 3.22e-05 | 0.64 | 1.33x |
| float | 243 | 243 | balanced | 1e-02 | adjoint | 320 | 2 | 1.38e-03 | 0.80 | 512 | 1 | 6.88e-03 | 1.01 | 1.26x |
| float | 243 | 243 | balanced | 1e-02 | forward | 320 | 2 | 4.07e-04 | 0.65 | 512 | 1 | 2.96e-03 | 0.77 | 1.19x |
| float | 243 | 243 | balanced | 1e-04 | adjoint | 320 | 4 | 2.68e-06 | 1.10 | 512 | 3 | 1.29e-06 | 0.98 | 0.90x |
| float | 243 | 243 | balanced | 1e-04 | forward | 320 | 3 | 1.37e-05 | 0.65 | 512 | 2 | 3.12e-05 | 0.84 | 1.29x |
| float | 243 | 972 | M-dominated | 1e-02 | adjoint | 320 | 2 | 3.98e-04 | 2.38 | 512 | 1 | 3.16e-03 | 2.06 | 0.86x |
| float | 243 | 972 | M-dominated | 1e-02 | forward | 320 | 2 | 3.39e-04 | 1.52 | 512 | 1 | 3.36e-03 | 1.74 | 1.14x |
| float | 243 | 972 | M-dominated | 1e-04 | adjoint | 320 | 4 | 9.30e-07 | 3.49 | 512 | 2 | 4.24e-05 | 2.59 | 0.74x |
| float | 243 | 972 | M-dominated | 1e-04 | forward | 320 | 3 | 1.69e-05 | 1.68 | 512 | 2 | 3.95e-05 | 1.89 | 1.13x |
| float | 250 | 62 | N-dominated | 1e-02 | adjoint | 320 | 2 | 4.46e-03 | 0.55 | 512 | 2 | 2.18e-04 | 0.64 | 1.15x |
| float | 250 | 62 | N-dominated | 1e-02 | forward | 320 | 2 | 3.80e-04 | 0.40 | 512 | 1 | 2.64e-03 | 0.57 | 1.42x |
| float | 250 | 62 | N-dominated | 1e-04 | adjoint | 320 | 4 | 1.63e-05 | 0.51 | 512 | 3 | 2.09e-06 | 0.66 | 1.30x |
| float | 250 | 62 | N-dominated | 1e-04 | forward | 320 | 4 | 9.95e-07 | 0.39 | 512 | 2 | 3.09e-05 | 0.59 | 1.53x |
| float | 250 | 250 | balanced | 1e-02 | adjoint | 320 | 2 | 1.72e-03 | 0.92 | 512 | 1 | 4.72e-03 | 0.83 | 0.91x |
| float | 250 | 250 | balanced | 1e-02 | forward | 320 | 2 | 4.15e-04 | 0.68 | 512 | 1 | 2.90e-03 | 0.78 | 1.16x |
| float | 250 | 250 | balanced | 1e-04 | adjoint | 320 | 4 | 3.10e-06 | 1.09 | 512 | 2 | 7.76e-05 | 0.93 | 0.85x |
| float | 250 | 250 | balanced | 1e-04 | forward | 320 | 4 | 1.01e-06 | 0.78 | 512 | 2 | 3.85e-05 | 0.88 | 1.13x |
| float | 250 | 1000 | M-dominated | 1e-02 | adjoint | 320 | 2 | 9.33e-04 | 2.30 | 512 | 1 | 3.13e-03 | 2.06 | 0.90x |
| float | 250 | 1000 | M-dominated | 1e-02 | forward | 320 | 2 | 4.88e-04 | 1.69 | 512 | 1 | 2.93e-03 | 1.55 | 0.92x |
| float | 250 | 1000 | M-dominated | 1e-04 | adjoint | 320 | 4 | 3.78e-06 | 3.29 | 512 | 2 | 5.26e-05 | 2.23 | 0.68x |
| float | 250 | 1000 | M-dominated | 1e-04 | forward | 320 | 4 | 1.79e-06 | 2.05 | 512 | 2 | 3.83e-05 | 1.80 | 0.88x |
| float | 251 | 62 | N-dominated | 1e-02 | adjoint | 320 | 2 | 2.07e-03 | 0.42 | 512 | 2 | 1.39e-04 | 0.66 | 1.59x |
| float | 251 | 62 | N-dominated | 1e-02 | forward | 320 | 2 | 3.43e-04 | 0.36 | 512 | 1 | 2.46e-03 | 0.57 | 1.59x |
| float | 251 | 62 | N-dominated | 1e-04 | adjoint | 320 | 4 | 6.15e-06 | 0.49 | 512 | 3 | 2.04e-06 | 0.60 | 1.22x |
| float | 251 | 62 | N-dominated | 1e-04 | forward | 320 | 4 | 1.06e-06 | 0.42 | 512 | 2 | 3.10e-05 | 0.58 | 1.39x |
| float | 251 | 251 | balanced | 1e-02 | adjoint | 320 | 2 | 1.56e-03 | 0.82 | 512 | 1 | 4.13e-03 | 0.83 | 1.01x |
| float | 251 | 251 | balanced | 1e-02 | forward | 320 | 2 | 4.69e-04 | 0.61 | 512 | 1 | 2.41e-03 | 0.78 | 1.27x |
| float | 251 | 251 | balanced | 1e-04 | adjoint | 320 | 4 | 2.78e-06 | 0.96 | 512 | 2 | 7.45e-05 | 0.91 | 0.95x |
| float | 251 | 251 | balanced | 1e-04 | forward | 320 | 4 | 1.59e-06 | 0.71 | 512 | 2 | 3.59e-05 | 0.87 | 1.23x |
| float | 251 | 1004 | M-dominated | 1e-02 | adjoint | 320 | 2 | 8.55e-04 | 2.26 | 512 | 1 | 2.78e-03 | 2.06 | 0.91x |
| float | 251 | 1004 | M-dominated | 1e-02 | forward | 320 | 2 | 4.08e-04 | 1.68 | 512 | 1 | 2.51e-03 | 1.58 | 0.94x |
| float | 251 | 1004 | M-dominated | 1e-04 | adjoint | 320 | 4 | 2.43e-06 | 3.24 | 512 | 2 | 5.01e-05 | 2.28 | 0.70x |
| float | 251 | 1004 | M-dominated | 1e-04 | forward | 320 | 4 | 1.65e-06 | 2.04 | 512 | 2 | 3.88e-05 | 1.83 | 0.90x |
| float | 255 | 63 | N-dominated | 1e-02 | adjoint | 320 | 2 | 2.88e-03 | 0.55 | 512 | 2 | 2.17e-04 | 0.60 | 1.08x |
| float | 255 | 63 | N-dominated | 1e-02 | forward | 320 | 2 | 4.80e-04 | 0.40 | 512 | 1 | 2.32e-03 | 0.58 | 1.44x |
| float | 255 | 63 | N-dominated | 1e-04 | adjoint | 320 | 4 | 7.73e-06 | 0.49 | 512 | 3 | 1.56e-06 | 0.66 | 1.36x |
| float | 255 | 63 | N-dominated | 1e-04 | forward | 320 | 4 | 1.31e-06 | 0.39 | 512 | 2 | 3.45e-05 | 0.60 | 1.55x |
| float | 255 | 255 | balanced | 1e-02 | adjoint | 320 | 2 | 2.41e-03 | 0.81 | 512 | 1 | 4.37e-03 | 0.94 | 1.16x |
| float | 255 | 255 | balanced | 1e-02 | forward | 320 | 2 | 4.46e-04 | 0.69 | 512 | 1 | 2.90e-03 | 0.76 | 1.10x |
| float | 255 | 255 | balanced | 1e-04 | adjoint | 320 | 4 | 7.85e-06 | 1.05 | 512 | 2 | 7.58e-05 | 0.94 | 0.89x |
| float | 255 | 255 | balanced | 1e-04 | forward | 320 | 4 | 1.83e-06 | 0.79 | 512 | 2 | 3.62e-05 | 0.86 | 1.09x |
| float | 255 | 1020 | M-dominated | 1e-02 | adjoint | 320 | 2 | 1.21e-03 | 2.31 | 512 | 1 | 3.23e-03 | 2.06 | 0.89x |
| float | 255 | 1020 | M-dominated | 1e-02 | forward | 320 | 2 | 6.12e-04 | 1.73 | 512 | 1 | 3.23e-03 | 1.68 | 0.97x |
| float | 255 | 1020 | M-dominated | 1e-04 | adjoint | 320 | 4 | 2.87e-06 | 3.38 | 512 | 2 | 5.99e-05 | 2.43 | 0.72x |
| float | 255 | 1020 | M-dominated | 1e-04 | forward | 320 | 4 | 2.01e-06 | 2.07 | 512 | 2 | 4.09e-05 | 1.71 | 0.82x |
| float | 256 | 64 | N-dominated | 1e-02 | adjoint | 320 | 2 | 5.68e-03 | 0.42 | 512 | 2 | 2.03e-04 | 0.60 | 1.43x |
| float | 256 | 64 | N-dominated | 1e-02 | forward | 320 | 2 | 5.60e-04 | 0.37 | 512 | 1 | 2.29e-03 | 0.52 | 1.40x |
| float | 256 | 64 | N-dominated | 1e-04 | adjoint | 320 | 4 | 1.48e-05 | 0.46 | 512 | 3 | 2.22e-06 | 0.68 | 1.47x |
| float | 256 | 64 | N-dominated | 1e-04 | forward | 320 | 4 | 1.53e-06 | 0.39 | 512 | 2 | 3.23e-05 | 0.54 | 1.41x |
| float | 256 | 256 | balanced | 1e-02 | adjoint | 320 | 2 | 1.70e-03 | 0.81 | 512 | 1 | 4.25e-03 | 0.85 | 1.05x |
| float | 256 | 256 | balanced | 1e-02 | forward | 320 | 2 | 4.65e-04 | 0.65 | 512 | 1 | 2.41e-03 | 0.77 | 1.19x |
| float | 256 | 256 | balanced | 1e-04 | adjoint | 320 | 4 | 7.61e-06 | 1.17 | 512 | 2 | 9.44e-05 | 1.06 | 0.90x |
| float | 256 | 256 | balanced | 1e-04 | forward | 320 | 4 | 1.71e-06 | 0.76 | 512 | 2 | 3.52e-05 | 0.85 | 1.13x |
| float | 256 | 1024 | M-dominated | 1e-02 | adjoint | 320 | 2 | 8.84e-04 | 2.33 | 512 | 1 | 2.31e-03 | 2.26 | 0.97x |
| float | 256 | 1024 | M-dominated | 1e-02 | forward | 320 | 2 | 5.44e-04 | 1.59 | 512 | 1 | 3.09e-03 | 1.51 | 0.95x |
| float | 256 | 1024 | M-dominated | 1e-04 | adjoint | 320 | 4 | 3.65e-06 | 3.44 | 512 | 2 | 3.92e-05 | 2.41 | 0.70x |
| float | 256 | 1024 | M-dominated | 1e-04 | forward | 320 | 4 | 1.79e-06 | 2.11 | 512 | 2 | 3.87e-05 | 1.72 | 0.81x |
| float | 512 | 128 | N-dominated | 1e-02 | adjoint | 640 | 2 | 3.40e-03 | 0.88 | 1024 | 1 | 8.45e-03 | 1.26 | 1.42x |
| float | 512 | 128 | N-dominated | 1e-02 | forward | 640 | 2 | 3.31e-04 | 0.99 | 1024 | 1 | 1.68e-03 | 1.11 | 1.12x |
| float | 512 | 128 | N-dominated | 1e-04 | adjoint | 640 | 4 | 1.62e-05 | 0.97 | 1024 | 3 | 2.27e-06 | 1.48 | 1.52x |
| float | 512 | 128 | N-dominated | 1e-04 | forward | 640 | 4 | 1.86e-06 | 0.88 | 1024 | 2 | 2.42e-05 | 1.14 | 1.30x |
| float | 512 | 512 | balanced | 1e-02 | adjoint | 640 | 2 | 1.66e-03 | 1.73 | 1024 | 1 | 4.58e-03 | 1.89 | 1.10x |
| float | 512 | 512 | balanced | 1e-02 | forward | 640 | 2 | 4.02e-04 | 1.40 | 1024 | 1 | 2.13e-03 | 1.62 | 1.15x |
| float | 512 | 512 | balanced | 1e-04 | adjoint | 640 | 4 | 6.65e-06 | 2.17 | 1024 | 2 | 6.27e-05 | 1.90 | 0.87x |
| float | 512 | 512 | balanced | 1e-04 | forward | 640 | 4 | 1.81e-06 | 1.65 | 1024 | 2 | 2.87e-05 | 1.77 | 1.07x |
| float | 512 | 2048 | M-dominated | 1e-02 | adjoint | 640 | 2 | 9.04e-04 | 4.34 | 1024 | 1 | 2.20e-03 | 3.97 | 0.92x |
| float | 512 | 2048 | M-dominated | 1e-02 | forward | 640 | 2 | 4.32e-04 | 3.40 | 1024 | 1 | 2.13e-03 | 3.42 | 1.01x |
| float | 512 | 2048 | M-dominated | 1e-04 | adjoint | 640 | 4 | 3.14e-06 | 6.23 | 1024 | 2 | 4.28e-05 | 4.83 | 0.78x |
| float | 512 | 2048 | M-dominated | 1e-04 | forward | 640 | 4 | 2.32e-06 | 4.18 | 1024 | 2 | 3.24e-05 | 3.94 | 0.94x |
| double | 243 | 60 | N-dominated | 1e-04 | adjoint | 320 | 4 | 4.30e-06 | 0.85 | 512 | 3 | 1.15e-06 | 1.16 | 1.36x |
| double | 243 | 60 | N-dominated | 1e-04 | forward | 320 | 3 | 1.13e-05 | 0.75 | 512 | 2 | 3.22e-05 | 1.20 | 1.60x |
| double | 243 | 60 | N-dominated | 1e-08 | adjoint | 320 | 7 | 2.29e-10 | 0.97 | 512 | 5 | 1.24e-10 | 1.28 | 1.31x |
| double | 243 | 60 | N-dominated | 1e-08 | forward | 320 | 6 | 8.63e-10 | 0.82 | 512 | 4 | 2.84e-09 | 1.04 | 1.26x |
| double | 243 | 60 | N-dominated | 1e-11 | adjoint | 320 | 9 | 5.09e-13 | 1.04 | 512 | 6 | 1.30e-12 | 1.32 | 1.26x |
| double | 243 | 60 | N-dominated | 1e-11 | forward | 320 | 9 | 6.28e-14 | 0.83 | 512 | 6 | 4.33e-13 | 1.08 | 1.30x |
| double | 243 | 243 | balanced | 1e-04 | adjoint | 320 | 4 | 2.86e-06 | 1.63 | 512 | 3 | 1.05e-06 | 1.88 | 1.15x |
| double | 243 | 243 | balanced | 1e-04 | forward | 320 | 3 | 1.38e-05 | 1.07 | 512 | 2 | 3.12e-05 | 1.25 | 1.17x |
| double | 243 | 243 | balanced | 1e-08 | adjoint | 320 | 7 | 1.33e-10 | 2.33 | 512 | 5 | 1.07e-10 | 2.13 | 0.92x |
| double | 243 | 243 | balanced | 1e-08 | forward | 320 | 6 | 9.60e-10 | 1.38 | 512 | 4 | 3.65e-09 | 1.44 | 1.04x |
| double | 243 | 243 | balanced | 1e-11 | adjoint | 320 | 9 | 2.83e-13 | 2.77 | 512 | 6 | 1.05e-12 | 2.52 | 0.91x |
| double | 243 | 243 | balanced | 1e-11 | forward | 320 | 9 | 7.09e-14 | 1.82 | 512 | 6 | 3.70e-13 | 1.70 | 0.93x |
| double | 243 | 972 | M-dominated | 1e-04 | adjoint | 320 | 4 | 5.49e-07 | 5.20 | 512 | 2 | 4.25e-05 | 3.85 | 0.74x |
| double | 243 | 972 | M-dominated | 1e-04 | forward | 320 | 3 | 1.71e-05 | 2.45 | 512 | 2 | 3.95e-05 | 2.34 | 0.96x |
| double | 243 | 972 | M-dominated | 1e-08 | adjoint | 320 | 7 | 5.78e-11 | 7.28 | 512 | 4 | 3.85e-09 | 5.55 | 0.76x |
| double | 243 | 972 | M-dominated | 1e-08 | forward | 320 | 6 | 1.03e-09 | 3.80 | 512 | 4 | 3.39e-09 | 3.10 | 0.82x |
| double | 243 | 972 | M-dominated | 1e-11 | adjoint | 320 | 9 | 1.07e-13 | 8.44 | 512 | 6 | 5.33e-13 | 6.89 | 0.82x |
| double | 243 | 972 | M-dominated | 1e-11 | forward | 320 | 9 | 5.67e-14 | 5.08 | 512 | 6 | 3.76e-13 | 3.97 | 0.78x |
| double | 250 | 62 | N-dominated | 1e-04 | adjoint | 320 | 4 | 1.61e-05 | 0.92 | 512 | 3 | 2.63e-06 | 1.14 | 1.25x |
| double | 250 | 62 | N-dominated | 1e-04 | forward | 320 | 4 | 8.85e-07 | 0.77 | 512 | 2 | 3.08e-05 | 1.09 | 1.41x |
| double | 250 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 3.33e-09 | 0.96 | 512 | 5 | 4.40e-10 | 1.36 | 1.42x |
| double | 250 | 62 | N-dominated | 1e-08 | forward | 320 | 7 | 1.24e-10 | 0.85 | 512 | 4 | 3.12e-09 | 1.24 | 1.46x |
| double | 250 | 62 | N-dominated | 1e-11 | adjoint | 320 | 10 | 6.36e-13 | 1.08 | 512 | 6 | 5.76e-12 | 1.28 | 1.18x |
| double | 250 | 62 | N-dominated | 1e-11 | forward | 320 | 9 | 4.38e-13 | 0.86 | 512 | 6 | 4.76e-13 | 1.29 | 1.51x |
| double | 250 | 250 | balanced | 1e-04 | adjoint | 320 | 4 | 3.40e-06 | 1.84 | 512 | 2 | 7.73e-05 | 1.65 | 0.90x |
| double | 250 | 250 | balanced | 1e-04 | forward | 320 | 4 | 9.13e-07 | 1.16 | 512 | 2 | 3.88e-05 | 1.25 | 1.07x |
| double | 250 | 250 | balanced | 1e-08 | adjoint | 320 | 7 | 4.95e-10 | 2.52 | 512 | 5 | 1.27e-10 | 2.47 | 0.98x |
| double | 250 | 250 | balanced | 1e-08 | forward | 320 | 7 | 1.24e-10 | 1.50 | 512 | 4 | 3.84e-09 | 1.41 | 0.93x |
| double | 250 | 250 | balanced | 1e-11 | adjoint | 320 | 10 | 9.97e-14 | 2.85 | 512 | 6 | 1.70e-12 | 2.40 | 0.84x |
| double | 250 | 250 | balanced | 1e-11 | forward | 320 | 9 | 3.71e-13 | 1.54 | 512 | 6 | 4.63e-13 | 1.66 | 1.08x |
| double | 250 | 1000 | M-dominated | 1e-04 | adjoint | 320 | 4 | 3.44e-06 | 5.17 | 512 | 2 | 5.26e-05 | 4.09 | 0.79x |
| double | 250 | 1000 | M-dominated | 1e-04 | forward | 320 | 4 | 1.06e-06 | 2.71 | 512 | 2 | 3.83e-05 | 2.23 | 0.82x |
| double | 250 | 1000 | M-dominated | 1e-08 | adjoint | 320 | 7 | 7.07e-10 | 7.07 | 512 | 4 | 5.02e-09 | 5.42 | 0.77x |
| double | 250 | 1000 | M-dominated | 1e-08 | forward | 320 | 7 | 1.52e-10 | 4.38 | 512 | 4 | 5.57e-09 | 2.93 | 0.67x |
| double | 250 | 1000 | M-dominated | 1e-11 | adjoint | 320 | 10 | 7.90e-14 | 9.73 | 512 | 6 | 7.71e-13 | 7.11 | 0.73x |
| double | 250 | 1000 | M-dominated | 1e-11 | forward | 320 | 9 | 4.77e-13 | 5.22 | 512 | 6 | 5.19e-13 | 4.14 | 0.79x |
| double | 251 | 62 | N-dominated | 1e-04 | adjoint | 320 | 4 | 5.95e-06 | 0.85 | 512 | 3 | 1.45e-06 | 1.19 | 1.39x |
| double | 251 | 62 | N-dominated | 1e-04 | forward | 320 | 4 | 8.89e-07 | 0.71 | 512 | 2 | 3.10e-05 | 1.05 | 1.48x |
| double | 251 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 6.15e-10 | 1.13 | 512 | 5 | 1.58e-10 | 1.24 | 1.09x |
| double | 251 | 62 | N-dominated | 1e-08 | forward | 320 | 7 | 8.76e-11 | 0.88 | 512 | 4 | 3.42e-09 | 1.15 | 1.30x |
| double | 251 | 62 | N-dominated | 1e-11 | adjoint | 324 | 10 | 1.32e-13 | 1.20 | 512 | 6 | 1.41e-12 | 1.28 | 1.07x |
| double | 251 | 62 | N-dominated | 1e-11 | forward | 320 | 9 | 2.24e-13 | 1.08 | 512 | 6 | 3.99e-13 | 1.20 | 1.12x |
| double | 251 | 251 | balanced | 1e-04 | adjoint | 320 | 4 | 2.66e-06 | 1.78 | 512 | 2 | 7.46e-05 | 1.66 | 0.94x |
| double | 251 | 251 | balanced | 1e-04 | forward | 320 | 4 | 1.07e-06 | 1.27 | 512 | 2 | 3.60e-05 | 1.34 | 1.06x |
| double | 251 | 251 | balanced | 1e-08 | adjoint | 320 | 7 | 4.12e-10 | 2.46 | 512 | 4 | 8.46e-09 | 2.12 | 0.86x |
| double | 251 | 251 | balanced | 1e-08 | forward | 320 | 7 | 1.03e-10 | 1.46 | 512 | 4 | 5.49e-09 | 1.47 | 1.00x |
| double | 251 | 251 | balanced | 1e-11 | adjoint | 324 | 10 | 6.30e-14 | 3.26 | 512 | 6 | 9.73e-13 | 2.43 | 0.75x |
| double | 251 | 251 | balanced | 1e-11 | forward | 320 | 9 | 2.77e-13 | 1.90 | 512 | 6 | 4.63e-13 | 1.67 | 0.88x |
| double | 251 | 1004 | M-dominated | 1e-04 | adjoint | 320 | 4 | 2.32e-06 | 5.40 | 512 | 2 | 4.99e-05 | 3.92 | 0.73x |
| double | 251 | 1004 | M-dominated | 1e-04 | forward | 320 | 4 | 1.04e-06 | 2.77 | 512 | 2 | 3.88e-05 | 2.17 | 0.78x |
| double | 251 | 1004 | M-dominated | 1e-08 | adjoint | 320 | 7 | 2.26e-10 | 7.67 | 512 | 4 | 4.42e-09 | 5.36 | 0.70x |
| double | 251 | 1004 | M-dominated | 1e-08 | forward | 320 | 7 | 1.06e-10 | 4.24 | 512 | 4 | 3.98e-09 | 2.99 | 0.71x |
| double | 251 | 1004 | M-dominated | 1e-11 | adjoint | 324 | 10 | 3.25e-14 | 10.08 | 512 | 6 | 4.47e-13 | 6.99 | 0.69x |
| double | 251 | 1004 | M-dominated | 1e-11 | forward | 320 | 9 | 3.48e-13 | 5.16 | 512 | 6 | 4.36e-13 | 4.24 | 0.82x |
| double | 255 | 63 | N-dominated | 1e-04 | adjoint | 320 | 4 | 7.87e-06 | 0.86 | 512 | 3 | 1.46e-06 | 1.23 | 1.42x |
| double | 255 | 63 | N-dominated | 1e-04 | forward | 320 | 4 | 1.17e-06 | 0.69 | 512 | 2 | 3.45e-05 | 1.09 | 1.58x |
| double | 255 | 63 | N-dominated | 1e-08 | adjoint | 320 | 7 | 8.64e-10 | 0.98 | 512 | 5 | 2.36e-10 | 1.26 | 1.29x |
| double | 255 | 63 | N-dominated | 1e-08 | forward | 320 | 7 | 1.88e-10 | 0.84 | 512 | 4 | 3.71e-09 | 1.06 | 1.26x |
| double | 255 | 63 | N-dominated | 1e-11 | adjoint | 360 | 8 | 1.18e-12 | 1.22 | 512 | 6 | 2.72e-12 | 1.42 | 1.16x |
| double | 255 | 63 | N-dominated | 1e-11 | forward | 324 | 9 | 3.46e-13 | 0.99 | 512 | 6 | 4.83e-13 | 1.17 | 1.18x |
| double | 255 | 255 | balanced | 1e-04 | adjoint | 320 | 4 | 7.95e-06 | 1.87 | 512 | 2 | 7.59e-05 | 1.91 | 1.03x |
| double | 255 | 255 | balanced | 1e-04 | forward | 320 | 4 | 1.11e-06 | 1.25 | 512 | 2 | 3.63e-05 | 1.48 | 1.19x |
| double | 255 | 255 | balanced | 1e-08 | adjoint | 320 | 7 | 9.51e-10 | 2.36 | 512 | 4 | 7.42e-09 | 2.27 | 0.96x |
| double | 255 | 255 | balanced | 1e-08 | forward | 320 | 7 | 2.09e-10 | 1.54 | 512 | 4 | 4.22e-09 | 1.68 | 1.09x |
| double | 255 | 255 | balanced | 1e-11 | adjoint | 360 | 8 | 8.52e-13 | 2.70 | 512 | 6 | 8.22e-13 | 2.56 | 0.95x |
| double | 255 | 255 | balanced | 1e-11 | forward | 324 | 9 | 2.91e-13 | 1.70 | 512 | 6 | 4.95e-13 | 1.91 | 1.12x |
| double | 255 | 1020 | M-dominated | 1e-04 | adjoint | 320 | 4 | 2.96e-06 | 5.36 | 512 | 2 | 5.99e-05 | 4.13 | 0.77x |
| double | 255 | 1020 | M-dominated | 1e-04 | forward | 320 | 4 | 1.92e-06 | 2.84 | 512 | 2 | 4.10e-05 | 2.23 | 0.79x |
| double | 255 | 1020 | M-dominated | 1e-08 | adjoint | 320 | 7 | 5.40e-10 | 7.63 | 512 | 4 | 5.71e-09 | 5.31 | 0.70x |
| double | 255 | 1020 | M-dominated | 1e-08 | forward | 320 | 7 | 2.13e-10 | 4.29 | 512 | 4 | 5.29e-09 | 3.39 | 0.79x |
| double | 255 | 1020 | M-dominated | 1e-11 | adjoint | 360 | 8 | 2.61e-13 | 8.81 | 512 | 6 | 5.97e-13 | 6.99 | 0.79x |
| double | 255 | 1020 | M-dominated | 1e-11 | forward | 324 | 9 | 3.16e-13 | 5.37 | 512 | 6 | 7.14e-13 | 4.26 | 0.79x |
| double | 256 | 64 | N-dominated | 1e-04 | adjoint | 320 | 4 | 1.41e-05 | 0.88 | 512 | 3 | 2.10e-06 | 1.19 | 1.34x |
| double | 256 | 64 | N-dominated | 1e-04 | forward | 320 | 4 | 1.37e-06 | 0.72 | 512 | 2 | 3.23e-05 | 1.05 | 1.46x |
| double | 256 | 64 | N-dominated | 1e-08 | adjoint | 320 | 7 | 4.00e-09 | 1.04 | 512 | 5 | 1.95e-10 | 1.28 | 1.23x |
| double | 256 | 64 | N-dominated | 1e-08 | forward | 320 | 7 | 2.87e-10 | 0.79 | 512 | 4 | 3.99e-09 | 1.11 | 1.41x |
| double | 256 | 64 | N-dominated | 1e-11 | adjoint | 360 | 8 | 3.69e-12 | 1.29 | 512 | 6 | 2.62e-12 | 1.33 | 1.03x |
| double | 256 | 64 | N-dominated | 1e-11 | forward | 324 | 9 | 5.94e-13 | 0.89 | 512 | 6 | 5.27e-13 | 1.17 | 1.31x |
| double | 256 | 256 | balanced | 1e-04 | adjoint | 320 | 4 | 7.64e-06 | 1.83 | 512 | 2 | 9.43e-05 | 1.73 | 0.94x |
| double | 256 | 256 | balanced | 1e-04 | forward | 320 | 4 | 1.33e-06 | 1.10 | 512 | 2 | 3.49e-05 | 1.24 | 1.13x |
| double | 256 | 256 | balanced | 1e-08 | adjoint | 320 | 7 | 2.38e-09 | 2.42 | 512 | 5 | 2.10e-10 | 2.27 | 0.94x |
| double | 256 | 256 | balanced | 1e-08 | forward | 320 | 7 | 2.80e-10 | 1.47 | 512 | 4 | 4.30e-09 | 1.74 | 1.18x |
| double | 256 | 256 | balanced | 1e-11 | adjoint | 360 | 8 | 8.31e-13 | 2.62 | 512 | 6 | 2.96e-12 | 2.44 | 0.93x |
| double | 256 | 256 | balanced | 1e-11 | forward | 324 | 9 | 6.91e-13 | 1.80 | 512 | 6 | 5.83e-13 | 1.61 | 0.89x |
| double | 256 | 1024 | M-dominated | 1e-04 | adjoint | 320 | 4 | 3.72e-06 | 5.25 | 512 | 2 | 3.93e-05 | 3.91 | 0.75x |
| double | 256 | 1024 | M-dominated | 1e-04 | forward | 320 | 4 | 1.73e-06 | 3.02 | 512 | 2 | 3.86e-05 | 2.52 | 0.83x |
| double | 256 | 1024 | M-dominated | 1e-08 | adjoint | 320 | 7 | 1.12e-09 | 7.78 | 512 | 4 | 4.28e-09 | 5.42 | 0.70x |
| double | 256 | 1024 | M-dominated | 1e-08 | forward | 320 | 7 | 3.40e-10 | 4.11 | 512 | 4 | 4.85e-09 | 3.27 | 0.80x |
| double | 256 | 1024 | M-dominated | 1e-11 | adjoint | 360 | 8 | 1.16e-12 | 8.80 | 512 | 6 | 5.31e-13 | 7.05 | 0.80x |
| double | 256 | 1024 | M-dominated | 1e-11 | forward | 324 | 9 | 7.86e-13 | 5.50 | 512 | 6 | 6.43e-13 | 4.30 | 0.78x |
| double | 512 | 128 | N-dominated | 1e-04 | adjoint | 640 | 4 | 1.63e-05 | 1.89 | 1024 | 3 | 2.01e-06 | 2.63 | 1.39x |
| double | 512 | 128 | N-dominated | 1e-04 | forward | 640 | 4 | 1.13e-06 | 1.50 | 1024 | 2 | 2.42e-05 | 2.30 | 1.53x |
| double | 512 | 128 | N-dominated | 1e-08 | adjoint | 640 | 7 | 4.90e-09 | 2.27 | 1024 | 5 | 3.69e-10 | 2.92 | 1.29x |
| double | 512 | 128 | N-dominated | 1e-08 | forward | 640 | 7 | 2.33e-10 | 1.85 | 1024 | 4 | 3.17e-09 | 2.42 | 1.31x |
| double | 512 | 128 | N-dominated | 1e-11 | adjoint | 720 | 8 | 2.84e-12 | 2.71 | 1024 | 6 | 5.02e-12 | 2.95 | 1.09x |
| double | 512 | 128 | N-dominated | 1e-11 | forward | 648 | 9 | 3.51e-13 | 2.07 | 1024 | 6 | 3.78e-13 | 2.48 | 1.19x |
| double | 512 | 512 | balanced | 1e-04 | adjoint | 640 | 4 | 6.36e-06 | 3.57 | 1024 | 2 | 6.25e-05 | 3.75 | 1.05x |
| double | 512 | 512 | balanced | 1e-04 | forward | 640 | 4 | 1.19e-06 | 2.39 | 1024 | 2 | 2.87e-05 | 2.76 | 1.16x |
| double | 512 | 512 | balanced | 1e-08 | adjoint | 640 | 7 | 1.91e-09 | 4.86 | 1024 | 4 | 6.93e-09 | 4.50 | 0.92x |
| double | 512 | 512 | balanced | 1e-08 | forward | 640 | 7 | 2.05e-10 | 3.09 | 1024 | 4 | 3.73e-09 | 3.49 | 1.13x |
| double | 512 | 512 | balanced | 1e-11 | adjoint | 720 | 8 | 8.28e-13 | 5.44 | 1024 | 6 | 1.01e-12 | 5.17 | 0.95x |
| double | 512 | 512 | balanced | 1e-11 | forward | 648 | 9 | 4.54e-13 | 3.70 | 1024 | 6 | 5.31e-13 | 3.61 | 0.98x |
| double | 512 | 2048 | M-dominated | 1e-04 | adjoint | 640 | 4 | 3.22e-06 | 10.71 | 1024 | 2 | 4.31e-05 | 9.01 | 0.84x |
| double | 512 | 2048 | M-dominated | 1e-04 | forward | 640 | 4 | 1.31e-06 | 6.45 | 1024 | 2 | 3.19e-05 | 5.21 | 0.81x |
| double | 512 | 2048 | M-dominated | 1e-08 | adjoint | 640 | 7 | 8.06e-10 | 15.53 | 1024 | 4 | 6.71e-09 | 11.63 | 0.75x |
| double | 512 | 2048 | M-dominated | 1e-08 | forward | 640 | 7 | 2.71e-10 | 9.18 | 1024 | 4 | 4.14e-09 | 6.99 | 0.76x |
| double | 512 | 2048 | M-dominated | 1e-11 | adjoint | 720 | 8 | 3.30e-13 | 16.78 | 1024 | 6 | 7.75e-13 | 14.96 | 0.89x |
| double | 512 | 2048 | M-dominated | 1e-11 | forward | 648 | 9 | 4.94e-13 | 11.28 | 1024 | 6 | 5.64e-13 | 8.86 | 0.78x |
| long double | 243 | 60 | N-dominated | 1e-08 | adjoint | 320 | 7 | 2.29e-10 | 87.99 | 512 | 5 | 1.24e-10 | 133.22 | 1.51x |
| long double | 243 | 60 | N-dominated | 1e-08 | forward | 320 | 6 | 8.63e-10 | 84.52 | 512 | 4 | 2.84e-09 | 122.94 | 1.45x |
| long double | 243 | 60 | N-dominated | 1e-14 | adjoint | 320 | 11 | 1.09e-15 | 102.89 | 512 | 8 | 1.48e-16 | 132.86 | 1.29x |
| long double | 243 | 60 | N-dominated | 1e-14 | forward | 320 | 11 | 1.41e-16 | 104.72 | 512 | 7 | 3.21e-15 | 126.42 | 1.21x |
| long double | 243 | 60 | N-dominated | 1e-20 | adjoint | 320 | 16 | 3.59e-22 | 120.06 | 512 | 11 | 1.90e-22 | 144.07 | 1.20x |
| long double | 243 | 60 | N-dominated | 1e-20 | forward | 320 | 15 | 8.35e-22 | 114.76 | 512 | 10 | 3.07e-21 | 143.47 | 1.25x |
| long double | 243 | 243 | balanced | 1e-08 | adjoint | 320 | 7 | 1.33e-10 | 195.06 | 512 | 5 | 1.07e-10 | 212.30 | 1.09x |
| long double | 243 | 243 | balanced | 1e-08 | forward | 320 | 6 | 9.60e-10 | 178.53 | 512 | 4 | 3.65e-09 | 173.81 | 0.97x |
| long double | 243 | 243 | balanced | 1e-14 | adjoint | 320 | 11 | 7.90e-16 | 251.29 | 512 | 7 | 7.20e-15 | 224.03 | 0.89x |
| long double | 243 | 243 | balanced | 1e-14 | forward | 320 | 11 | 1.49e-16 | 290.00 | 512 | 7 | 3.28e-15 | 228.80 | 0.79x |
| long double | 243 | 243 | balanced | 1e-20 | adjoint | 320 | 16 | 2.52e-22 | 322.03 | 512 | 11 | 1.12e-22 | 289.76 | 0.90x |
| long double | 243 | 243 | balanced | 1e-20 | forward | 320 | 15 | 7.11e-22 | 417.48 | 512 | 10 | 3.16e-21 | 279.18 | 0.67x |
| long double | 243 | 972 | M-dominated | 1e-08 | adjoint | 320 | 7 | 5.78e-11 | 623.09 | 512 | 4 | 3.85e-09 | 542.19 | 0.87x |
| long double | 243 | 972 | M-dominated | 1e-08 | forward | 320 | 6 | 1.03e-09 | 578.78 | 512 | 4 | 3.39e-09 | 534.63 | 0.92x |
| long double | 243 | 972 | M-dominated | 1e-14 | adjoint | 320 | 11 | 2.55e-16 | 929.20 | 512 | 7 | 5.32e-15 | 685.65 | 0.74x |
| long double | 243 | 972 | M-dominated | 1e-14 | forward | 320 | 11 | 1.53e-16 | 988.44 | 512 | 7 | 3.47e-15 | 735.13 | 0.74x |
| long double | 243 | 972 | M-dominated | 1e-20 | adjoint | 320 | 16 | 4.98e-23 | 1265.73 | 512 | 10 | 4.07e-21 | 904.96 | 0.71x |
| long double | 243 | 972 | M-dominated | 1e-20 | forward | 320 | 15 | 7.43e-22 | 1300.07 | 512 | 10 | 3.77e-21 | 1002.78 | 0.77x |
| long double | 250 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 3.33e-09 | 109.28 | 512 | 5 | 4.40e-10 | 132.24 | 1.21x |
| long double | 250 | 62 | N-dominated | 1e-08 | forward | 320 | 7 | 1.24e-10 | 93.23 | 512 | 4 | 3.12e-09 | 130.22 | 1.40x |
| long double | 250 | 62 | N-dominated | 1e-14 | adjoint | 320 | 12 | 1.73e-15 | 107.90 | 512 | 8 | 8.25e-16 | 131.09 | 1.21x |
| long double | 250 | 62 | N-dominated | 1e-14 | forward | 320 | 11 | 9.79e-16 | 104.03 | 512 | 7 | 4.64e-15 | 127.53 | 1.23x |
| long double | 250 | 62 | N-dominated | 1e-20 | adjoint | 320 | 17 | 8.50e-22 | 128.59 | 512 | 11 | 1.34e-21 | 150.10 | 1.17x |
| long double | 250 | 62 | N-dominated | 1e-20 | forward | 320 | 16 | 4.50e-22 | 122.24 | 512 | 10 | 7.38e-21 | 147.81 | 1.21x |
| long double | 250 | 250 | balanced | 1e-08 | adjoint | 320 | 7 | 4.95e-10 | 186.77 | 512 | 5 | 1.27e-10 | 219.84 | 1.18x |
| long double | 250 | 250 | balanced | 1e-08 | forward | 320 | 7 | 1.24e-10 | 202.05 | 512 | 4 | 3.84e-09 | 172.72 | 0.85x |
| long double | 250 | 250 | balanced | 1e-14 | adjoint | 320 | 12 | 1.91e-16 | 294.42 | 512 | 8 | 2.43e-16 | 273.64 | 0.93x |
| long double | 250 | 250 | balanced | 1e-14 | forward | 320 | 11 | 1.38e-15 | 304.83 | 512 | 7 | 4.82e-15 | 256.77 | 0.84x |
| long double | 250 | 250 | balanced | 1e-20 | adjoint | 320 | 17 | 8.91e-23 | 344.23 | 512 | 11 | 3.93e-22 | 297.97 | 0.87x |
| long double | 250 | 250 | balanced | 1e-20 | forward | 320 | 16 | 4.73e-22 | 398.58 | 512 | 10 | 6.12e-21 | 344.08 | 0.86x |
| long double | 250 | 1000 | M-dominated | 1e-08 | adjoint | 320 | 7 | 7.07e-10 | 650.13 | 512 | 4 | 5.02e-09 | 489.14 | 0.75x |
| long double | 250 | 1000 | M-dominated | 1e-08 | forward | 320 | 7 | 1.52e-10 | 795.44 | 512 | 4 | 5.57e-09 | 480.53 | 0.60x |
| long double | 250 | 1000 | M-dominated | 1e-14 | adjoint | 320 | 12 | 3.72e-16 | 1001.70 | 512 | 7 | 9.32e-15 | 695.24 | 0.69x |
| long double | 250 | 1000 | M-dominated | 1e-14 | forward | 320 | 11 | 1.37e-15 | 1007.95 | 512 | 7 | 5.79e-15 | 736.91 | 0.73x |
| long double | 250 | 1000 | M-dominated | 1e-20 | adjoint | 320 | 17 | 1.90e-22 | 1383.65 | 512 | 11 | 1.80e-22 | 995.00 | 0.72x |
| long double | 250 | 1000 | M-dominated | 1e-20 | forward | 320 | 16 | 5.57e-22 | 1413.57 | 512 | 10 | 9.36e-21 | 1021.16 | 0.72x |
| long double | 251 | 62 | N-dominated | 1e-08 | adjoint | 320 | 7 | 6.15e-10 | 89.41 | 512 | 5 | 1.58e-10 | 135.76 | 1.52x |
| long double | 251 | 62 | N-dominated | 1e-08 | forward | 320 | 7 | 8.76e-11 | 91.01 | 512 | 4 | 3.42e-09 | 133.61 | 1.47x |
| long double | 251 | 62 | N-dominated | 1e-14 | adjoint | 320 | 12 | 3.67e-16 | 146.74 | 512 | 8 | 2.32e-16 | 134.53 | 0.92x |
| long double | 251 | 62 | N-dominated | 1e-14 | forward | 320 | 11 | 6.10e-16 | 105.38 | 512 | 7 | 4.68e-15 | 130.69 | 1.24x |
| long double | 251 | 62 | N-dominated | 1e-20 | adjoint | 320 | 17 | 1.61e-22 | 148.88 | 512 | 11 | 2.85e-22 | 151.67 | 1.02x |
| long double | 251 | 62 | N-dominated | 1e-20 | forward | 320 | 16 | 3.49e-22 | 121.48 | 512 | 10 | 5.17e-21 | 151.42 | 1.25x |
| long double | 251 | 251 | balanced | 1e-08 | adjoint | 320 | 7 | 4.12e-10 | 206.37 | 512 | 4 | 8.46e-09 | 176.43 | 0.85x |
| long double | 251 | 251 | balanced | 1e-08 | forward | 320 | 7 | 1.03e-10 | 194.90 | 512 | 4 | 5.49e-09 | 220.14 | 1.13x |
| long double | 251 | 251 | balanced | 1e-14 | adjoint | 320 | 12 | 1.62e-16 | 295.37 | 512 | 7 | 7.51e-15 | 226.63 | 0.77x |
| long double | 251 | 251 | balanced | 1e-14 | forward | 320 | 11 | 7.02e-16 | 295.33 | 512 | 7 | 4.78e-15 | 233.32 | 0.79x |
| long double | 251 | 251 | balanced | 1e-20 | adjoint | 320 | 17 | 7.79e-23 | 351.18 | 512 | 11 | 1.48e-22 | 300.79 | 0.86x |
| long double | 251 | 251 | balanced | 1e-20 | forward | 320 | 16 | 4.54e-22 | 403.88 | 512 | 10 | 5.38e-21 | 295.98 | 0.73x |
| long double | 251 | 1004 | M-dominated | 1e-08 | adjoint | 320 | 7 | 2.26e-10 | 644.67 | 512 | 4 | 4.42e-09 | 484.25 | 0.75x |
| long double | 251 | 1004 | M-dominated | 1e-08 | forward | 320 | 7 | 1.06e-10 | 717.28 | 512 | 4 | 3.98e-09 | 480.43 | 0.67x |
| long double | 251 | 1004 | M-dominated | 1e-14 | adjoint | 320 | 12 | 1.44e-16 | 1008.00 | 512 | 7 | 4.88e-15 | 712.89 | 0.71x |
| long double | 251 | 1004 | M-dominated | 1e-14 | forward | 320 | 11 | 6.25e-16 | 1106.41 | 512 | 7 | 5.38e-15 | 910.81 | 0.82x |
| long double | 251 | 1004 | M-dominated | 1e-20 | adjoint | 320 | 17 | 6.66e-23 | 1393.53 | 512 | 11 | 9.87e-23 | 1009.17 | 0.72x |
| long double | 251 | 1004 | M-dominated | 1e-20 | forward | 320 | 16 | 2.98e-22 | 1383.03 | 512 | 10 | 7.31e-21 | 1003.26 | 0.73x |
| long double | 255 | 63 | N-dominated | 1e-08 | adjoint | 320 | 7 | 8.64e-10 | 90.25 | 512 | 5 | 2.36e-10 | 129.81 | 1.44x |
| long double | 255 | 63 | N-dominated | 1e-08 | forward | 320 | 7 | 1.88e-10 | 89.59 | 512 | 4 | 3.71e-09 | 128.30 | 1.43x |
| long double | 255 | 63 | N-dominated | 1e-14 | adjoint | 320 | 12 | 8.43e-16 | 111.73 | 512 | 8 | 4.52e-16 | 132.95 | 1.19x |
| long double | 255 | 63 | N-dominated | 1e-14 | forward | 320 | 12 | 9.96e-17 | 110.14 | 512 | 7 | 4.88e-15 | 135.13 | 1.23x |
| long double | 255 | 63 | N-dominated | 1e-20 | adjoint | 320 | 17 | 6.20e-22 | 128.82 | 512 | 11 | 4.62e-22 | 151.59 | 1.18x |
| long double | 255 | 63 | N-dominated | 1e-20 | forward | 320 | 17 | 6.24e-23 | 131.87 | 512 | 10 | 7.43e-21 | 145.08 | 1.10x |
| long double | 255 | 255 | balanced | 1e-08 | adjoint | 320 | 7 | 9.51e-10 | 200.96 | 512 | 4 | 7.42e-09 | 173.66 | 0.86x |
| long double | 255 | 255 | balanced | 1e-08 | forward | 320 | 7 | 2.09e-10 | 187.95 | 512 | 4 | 4.22e-09 | 220.52 | 1.17x |
| long double | 255 | 255 | balanced | 1e-14 | adjoint | 320 | 12 | 1.11e-15 | 272.93 | 512 | 8 | 1.20e-16 | 282.28 | 1.03x |
| long double | 255 | 255 | balanced | 1e-14 | forward | 320 | 12 | 1.04e-16 | 321.61 | 512 | 7 | 6.26e-15 | 234.80 | 0.73x |
| long double | 255 | 255 | balanced | 1e-20 | adjoint | 320 | 17 | 6.27e-22 | 348.88 | 512 | 11 | 1.80e-22 | 295.51 | 0.85x |
| long double | 255 | 255 | balanced | 1e-20 | forward | 320 | 17 | 8.25e-23 | 430.58 | 512 | 10 | 7.00e-21 | 343.56 | 0.80x |
| long double | 255 | 1020 | M-dominated | 1e-08 | adjoint | 320 | 7 | 5.40e-10 | 711.17 | 512 | 4 | 5.71e-09 | 489.34 | 0.69x |
| long double | 255 | 1020 | M-dominated | 1e-08 | forward | 320 | 7 | 2.13e-10 | 707.34 | 512 | 4 | 5.29e-09 | 528.35 | 0.75x |
| long double | 255 | 1020 | M-dominated | 1e-14 | adjoint | 320 | 12 | 3.65e-16 | 1013.74 | 512 | 7 | 7.22e-15 | 723.53 | 0.71x |
| long double | 255 | 1020 | M-dominated | 1e-14 | forward | 320 | 12 | 1.65e-16 | 1088.21 | 512 | 7 | 6.52e-15 | 768.89 | 0.71x |
| long double | 255 | 1020 | M-dominated | 1e-20 | adjoint | 320 | 17 | 2.44e-22 | 1404.24 | 512 | 11 | 1.52e-22 | 1020.90 | 0.73x |
| long double | 255 | 1020 | M-dominated | 1e-20 | forward | 320 | 17 | 1.00e-22 | 1542.64 | 512 | 10 | 9.17e-21 | 1035.75 | 0.67x |
| long double | 256 | 64 | N-dominated | 1e-08 | adjoint | 320 | 7 | 4.00e-09 | 93.66 | 512 | 5 | 1.95e-10 | 141.66 | 1.51x |
| long double | 256 | 64 | N-dominated | 1e-08 | forward | 320 | 7 | 2.87e-10 | 93.95 | 512 | 4 | 3.99e-09 | 131.44 | 1.40x |
| long double | 256 | 64 | N-dominated | 1e-14 | adjoint | 320 | 12 | 3.73e-15 | 105.07 | 512 | 8 | 4.08e-16 | 144.05 | 1.37x |
| long double | 256 | 64 | N-dominated | 1e-14 | forward | 320 | 12 | 2.95e-16 | 118.15 | 512 | 7 | 6.97e-15 | 276.05 | 2.34x |
| long double | 256 | 64 | N-dominated | 1e-20 | adjoint | 320 | 17 | 3.65e-21 | 131.36 | 512 | 11 | 7.69e-22 | 152.84 | 1.16x |
| long double | 256 | 64 | N-dominated | 1e-20 | forward | 320 | 17 | 2.15e-22 | 138.67 | 512 | 11 | 1.51e-22 | 171.96 | 1.24x |
| long double | 256 | 256 | balanced | 1e-08 | adjoint | 320 | 7 | 2.38e-09 | 202.33 | 512 | 5 | 2.10e-10 | 216.68 | 1.07x |
| long double | 256 | 256 | balanced | 1e-08 | forward | 320 | 7 | 2.80e-10 | 184.97 | 512 | 4 | 4.30e-09 | 222.16 | 1.20x |
| long double | 256 | 256 | balanced | 1e-14 | adjoint | 320 | 12 | 2.39e-15 | 297.33 | 512 | 8 | 4.68e-16 | 282.03 | 0.95x |
| long double | 256 | 256 | balanced | 1e-14 | forward | 320 | 12 | 2.52e-16 | 323.14 | 512 | 7 | 7.25e-15 | 236.31 | 0.73x |
| long double | 256 | 256 | balanced | 1e-20 | adjoint | 320 | 17 | 2.13e-21 | 356.05 | 512 | 11 | 8.82e-22 | 342.93 | 0.96x |
| long double | 256 | 256 | balanced | 1e-20 | forward | 320 | 17 | 2.24e-22 | 398.29 | 512 | 11 | 1.31e-22 | 323.31 | 0.81x |
| long double | 256 | 1024 | M-dominated | 1e-08 | adjoint | 320 | 7 | 1.12e-09 | 664.05 | 512 | 4 | 4.28e-09 | 549.68 | 0.83x |
| long double | 256 | 1024 | M-dominated | 1e-08 | forward | 320 | 7 | 3.40e-10 | 733.75 | 512 | 4 | 4.85e-09 | 482.88 | 0.66x |
| long double | 256 | 1024 | M-dominated | 1e-14 | adjoint | 320 | 12 | 1.14e-15 | 1011.63 | 512 | 7 | 7.31e-15 | 718.91 | 0.71x |
| long double | 256 | 1024 | M-dominated | 1e-14 | forward | 320 | 12 | 2.86e-16 | 1110.45 | 512 | 7 | 8.76e-15 | 794.79 | 0.72x |
| long double | 256 | 1024 | M-dominated | 1e-20 | adjoint | 320 | 17 | 1.09e-21 | 1391.44 | 512 | 11 | 1.61e-22 | 1026.50 | 0.74x |
| long double | 256 | 1024 | M-dominated | 1e-20 | forward | 320 | 17 | 2.44e-22 | 1527.44 | 512 | 11 | 1.55e-22 | 1093.71 | 0.72x |
| long double | 512 | 128 | N-dominated | 1e-08 | adjoint | 640 | 7 | 4.90e-09 | 248.13 | 1024 | 5 | 3.69e-10 | 380.77 | 1.53x |
| long double | 512 | 128 | N-dominated | 1e-08 | forward | 640 | 7 | 2.33e-10 | 252.43 | 1024 | 4 | 3.17e-09 | 417.69 | 1.65x |
| long double | 512 | 128 | N-dominated | 1e-14 | adjoint | 640 | 12 | 4.92e-15 | 291.71 | 1024 | 8 | 7.98e-16 | 380.46 | 1.30x |
| long double | 512 | 128 | N-dominated | 1e-14 | forward | 640 | 12 | 1.57e-16 | 307.50 | 1024 | 7 | 5.62e-15 | 370.50 | 1.20x |
| long double | 512 | 128 | N-dominated | 1e-20 | adjoint | 640 | 17 | 4.90e-21 | 300.31 | 1024 | 11 | 1.55e-21 | 477.13 | 1.59x |
| long double | 512 | 128 | N-dominated | 1e-20 | forward | 640 | 17 | 1.33e-22 | 357.93 | 1024 | 10 | 7.38e-21 | 421.69 | 1.18x |
| long double | 512 | 512 | balanced | 1e-08 | adjoint | 640 | 7 | 1.91e-09 | 451.62 | 1024 | 4 | 6.93e-09 | 513.07 | 1.14x |
| long double | 512 | 512 | balanced | 1e-08 | forward | 640 | 7 | 2.05e-10 | 470.51 | 1024 | 4 | 3.73e-09 | 611.23 | 1.30x |
| long double | 512 | 512 | balanced | 1e-14 | adjoint | 640 | 12 | 1.90e-15 | 636.20 | 1024 | 7 | 8.35e-15 | 628.26 | 0.99x |
| long double | 512 | 512 | balanced | 1e-14 | forward | 640 | 12 | 1.81e-16 | 687.38 | 1024 | 7 | 4.99e-15 | 651.63 | 0.95x |
| long double | 512 | 512 | balanced | 1e-20 | adjoint | 640 | 17 | 1.96e-21 | 809.97 | 1024 | 11 | 1.82e-22 | 768.71 | 0.95x |
| long double | 512 | 512 | balanced | 1e-20 | forward | 640 | 17 | 1.45e-22 | 884.35 | 1024 | 10 | 9.06e-21 | 772.40 | 0.87x |
| long double | 512 | 2048 | M-dominated | 1e-08 | adjoint | 640 | 7 | 8.06e-10 | 1467.68 | 1024 | 4 | 6.71e-09 | 1200.60 | 0.82x |
| long double | 512 | 2048 | M-dominated | 1e-08 | forward | 640 | 7 | 2.71e-10 | 1562.39 | 1024 | 4 | 4.14e-09 | 1240.06 | 0.79x |
| long double | 512 | 2048 | M-dominated | 1e-14 | adjoint | 640 | 12 | 8.48e-16 | 2200.66 | 1024 | 7 | 7.48e-15 | 1644.28 | 0.75x |
| long double | 512 | 2048 | M-dominated | 1e-14 | forward | 640 | 12 | 2.27e-16 | 2352.63 | 1024 | 7 | 6.18e-15 | 1769.69 | 0.75x |
| long double | 512 | 2048 | M-dominated | 1e-20 | adjoint | 640 | 17 | 7.73e-22 | 2932.27 | 1024 | 11 | 1.33e-22 | 2224.73 | 0.76x |
| long double | 512 | 2048 | M-dominated | 1e-20 | forward | 640 | 17 | 1.61e-22 | 3222.18 | 1024 | 10 | 8.39e-21 | 2300.18 | 0.71x |

`!` marks a goal that was not met.

