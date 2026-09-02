# 02 — Recentre the fast-Gaussian recurrence (fixes NaN output in single precision)

Status: needs-triage

## Problem

`FG_PSI` and `PRE_FG_PSI` factorise the run as

```
phi(nx0 - l) = phi(nx0) * q^l * e(l),   q = exp(2 nx0 / b),  e(l) = exp(-l^2/b)
```

anchored at `l = 0`, the far left end of the run, where `phi` is smallest. The
running factor has to climb the window's whole dynamic range: at sigma=2, m=13
it reaches `q^(2m+1) = 3.2e59` while `e(l)` falls to `4.2e-58`.

Single precision tops out at 3.4e38. From `m = 8` the factor overflows to `+inf`
and the table underflows to zero; the product is NaN.

Reproduced against `libnfft3f` (float, Gaussian, sigma=2, N=256, M=2000):

```
  FG_PSI       m= 7   non-finite outputs: 0/2000
  FG_PSI       m= 9   non-finite outputs: 2000/2000
  PRE_FG_PSI   m= 9   non-finite outputs: 2000/2000
  FG_PSI       m=13   non-finite outputs: 2000/2000
```

CI does not catch it because the tests only instantiate
`WINDOW_HELP_ESTIMATE_m`, which is 5 for float.

Inside its range the anchoring still costs accuracy — window values against a
long double reference, worst over the run, relative to the run peak:

| precision | m  | as implemented | centred | direct  |
|-----------|----|----------------|---------|---------|
| float     | 7  | 1.2e-06        | 5.6e-08 | 5.6e-08 |
| double    | 8  | 3.3e-15        | 1.6e-16 | 1.6e-16 |
| double    | 13 | 1.9e-14        | 8.3e-17 | 1.2e-16 |

At m=13 that is the same order as the whole transform error. `make check` shows
it today: the `nfst_online` 2-D case reports `PRE FG PSI` at 5.04e-14 against
`PRE PSI` at 1.45e-14 on identical inputs.

## Change

Anchor at the grid point nearest the node. With `c = LRINT(nx0)` and
`r = nx0 - c`, `|r| <= 1/2`:

```
phi(nx0 - l) = phi(r) * q^(l-c) * e(|l-c|)
```

The running factor then stays inside `exp(+-(m+1)/b)` — 12 to 19 for the `b`
values this library uses — bounded in every precision at every `m`, and it is
stepped at most `m+1` times from the peak of the run instead of `2m+1` times
from its smallest value.

Equivalently and more cheaply, walk the neighbour ratio outward from `c`:
`w *= g; g *= exp(-2/b)` in each direction, with `exp(-2/b)` a per-axis constant
and `g` seeded at `exp((2r-1)/b)` upward and `exp(-2/b)/exp((2r-1)/b)` downward.
Two `EXP`s and one division for the whole run, two multiplies per point, both
loops branch-free. Measured at ~1 ulp; see `04-phi-run.md`.

## Sites

- `kernel/nfft/nfft.c`: `nfft_1d_init_fg_exp_l` and every `FG_PSI` /
  `PRE_FG_PSI` block (1-D, 2-D, 3-D, d-D, serial and OpenMP), and
  `X(precompute_fg_psi)`.
- `kernel/nfct/nfct.c`, `kernel/nfst/nfst.c`: the same blocks.

`PRE_FG_PSI` can keep its two-reals-per-node layout: store `phi(r)` and the
upward seed, from which `c` is recoverable (`c = m` when the seed is >= the
downward one, `m + 1` otherwise) — or store `r` directly if that reads better.

## Tests

A `tests/window.c` check that the run agrees with per-point direct evaluation at
every offset of the node inside a cell and at cut-offs well past the default,
and that no entry is non-finite — the last part is what would have caught the
float bug. #263's `check_bspline_run` is the model.
