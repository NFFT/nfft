# 02 — Fit per-band constants and emit them

Status: done

## Problem

`kernel/nfft/tune.c` carries one fitted envelope per direction, raised until
it dominates every measurement anywhere in `sigma` in [5/4, 4]. The dyadic
tuner uses three disjoint sigma bands, so each rung can have its own envelope,
which only has to dominate its own band.

## What to do

`dfit.py` already fits and validates per band. Extend it to:

- take the band `(4, 8]` seriously once issue 01 lands (today it collapses,
  because only `sigma = 4` falls in it and the least-squares design matrix is
  then singular — the `u` column is constant). Make the band bounds
  `[1.25, 2]`, `[2, 4]`, `[4, 8]` and report a clear message rather than a
  traceback when a band has fewer than two distinct sigma values.
- emit the six coefficient sets as C, in the form
  `kernel/nfft/tune.c` can paste, named `tune_dyadic_fwd[3]` and
  `tune_dyadic_adj[3]` indexed by rung.
- validate each band against the shipped constants on the same
  band-restricted population, as it does now.

## Acceptance

Per band and direction, against the measured smallest sufficient cut-off over
a ladder of goals:

- **misses must be 0**, as for the shipped constants;
- mean extra cut-off no worse than the shipped constants on the same
  population.

Recorded on the existing data for bands 0 and 1:

| band | direction | shipped | per-band |
|---|---|---|---|
| `[1.25, 2]` | forward | +0.44, 56 % exact | +0.38, 62 % |
| `[1.25, 2]` | adjoint | +0.27, 73 % | +0.27, 73 % |
| `[2, 4]` | forward | +0.30, 70 % | +0.16, 84 % |
| `[2, 4]` | adjoint | +0.20, 80 % | +0.19, 81 % |

## Note

`selectform.py` exists because a form that fits better in sample can pick a
worse cut-off. Do not change the functional form here — only the fitting
population changes.
