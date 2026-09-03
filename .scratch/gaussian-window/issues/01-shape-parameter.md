# 01 — Set the Gaussian `b` for the truncation width the stencil actually has

Status: ready-for-human

Implemented on branch `feature/gaussian-shape-parameter` (`d5b233b`), branched
from `develop` so it can be filed on its own — it is a semantic change and wants
reviewing apart from the `FG_PSI` correctness fix in
[`02-centred-fast-gaussian.md`](02-centred-fast-gaussian.md). `make check` is
green there in float, double and long double with no test-data changes.

Note for whoever reviews the two together: the wider `b` lowers the `FG_PSI`
running factor, so in float the m=9 failure drops from every node to 187 of
2000 while m >= 11 stays wholly non-finite. Landing this first masks the bug a
little; it does not fix it.

## Problem

`WINDOW_HELP_INIT` (`include/infft.h`, `GAUSSIAN` branch) sets

```c
ths->b[i] = (K(2.0)*ths->sigma[i]) / (K(2.0)*ths->sigma[i] - K(1.0)) * (((R)ths->m) / KPI);
```

That is the classical balance between truncation and aliasing error for a
Gaussian truncated at half-width `m`. The library's stencil is not truncated at
`m`: `uo()` places the run at `u = floor(n x) - m` over `2m + 2` points, so the
largest `|n x - l|` in a run lies between `m` and `m + 1`. `b` is therefore too
small and the window is pinched harder than the truncation requires.

## Change

```c
ths->b[i] = (K(2.0)*ths->sigma[i]) / (K(2.0)*ths->sigma[i] - K(1.0))
          * (((R)ths->m + K(1.0)) / KPI);
```

## Justification

The minimiser of the exact single-frequency error, computed in mpmath at 50
digits over the same `2m+2` stencil, tracks `(m + c)` with `c` in 0.9–1.2 for
every `sigma` in [1.25, 3] and `m` in [2, 16]. A scan of `c` in the library
(0.0, 0.8, 1.0, 1.2, 1.4) puts the optimum at 1.0–1.2 in every configuration.
`c = 1` is within a few percent of the true minimum everywhere and is the value
with a structural reason behind it.

Measured `E_inf` (N = 256, `PRE_PSI`, deterministic sub-cell sweep):

| case             | today    | `m+1`    | gain |
|------------------|----------|----------|------|
| sigma=2, m=4     | 7.83e-06 | 1.14e-06 | 6.9x |
| sigma=2, m=7     | 1.10e-08 | 1.43e-09 | 7.7x |
| sigma=2, m=10    | 1.70e-11 | 2.79e-12 | 6.1x |
| sigma=2, m=13    | 2.69e-14 | 4.53e-15 | 5.9x |
| sigma=1.25, m=8  | 3.64e-06 | 1.57e-06 | 2.3x |
| sigma=1.5, m=8   | 6.57e-08 | 1.38e-08 | 4.8x |
| sigma=3, m=8     | 5.73e-11 | 6.76e-12 | 8.5x |

No configuration regresses. Cost is zero: `b` is computed once per axis per plan.

## Scope and risk

This is a semantic change — every Gaussian result moves. It should be its own
commit.

- `make check` passes unchanged with the shift applied (2355 assertions).
- `ths->b` is read by the `FG_PSI` / `PRE_FG_PSI` paths in `nfft.c`, `nfct.c`
  and `nfst.c`; they consume whatever value `WINDOW_HELP_INIT` sets, so nothing
  there needs touching.
- `kernel/nsfft` requires the Gaussian window but goes through ordinary NFFT
  plans and never reads `ths->b` itself.
- Worth re-running the accuracy report (`NFFT_BENCH_OUT`) so the improvement is
  recorded rather than read as drift.

## Tests

A `tests/window.c` check that `phi` and `phi_hut` agree with an mpmath reference
generated under `tests/windowref/` for a fixed `(n, N, m)`, per #263's pattern.
The reference has to be generated from the new `b`, so it lands with this change.
