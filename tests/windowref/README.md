# Window reference values

Prints the reference tables that `tests/window.c` bounds the B-spline and
sinc-power window evaluations against:

```bash
uv run --with mpmath==1.3.0 python -m tests.windowref.generate
```

The output goes to stdout for pasting into `tests/window.c`. Nothing in the
tree is written, and there is no generated header, so the values stay visible
next to the test that uses them, as the Kaiser-Bessel tables do.

## What each table covers

All four run at `n = 512`, `N = 256`, so `sigma = 2`.

| table | evaluation | how |
|---|---|---|
| `bspline_phi_hut_ref_512_256_11` | `(sin t/t)^(2m)/n`, `t = k pi/n` | mpmath, 60 digits |
| `bspline_phi_ref_512_11` | `B_2m(t-l)/n` over one full run | exact rational |
| `sincpow_phi_ref_512_8` | `w (sin a/a)^(2m)` over one full run | mpmath, 60 digits |
| `sincpow_phi_hut_ref_512_256_8` | `B_2m(k/(w n) + m)` | exact rational |

## Two choices worth knowing

**The B-spline values are exact, not `mpmath`.** The cardinal B-spline comes
from the truncated power basis,

    B_k(x) = 1/(k-1)! sum_j (-1)^j C(k,j) (x-j)_+^(k-1)

evaluated over `Fraction`. That formula cancels catastrophically in floating
point, which is why `kernel/util/bspline.c` uses de Boor instead, but over the
rationals it is exact and it shares nothing with the scheme under test.

**The run offset is `5/16` of a grid cell, not `0.3`.** A dyadic fraction is
exact in every precision, so no part of the measured error is the rounding of
the offset.

The sinc-power tables run at `m = 8` rather than the default cutoff, because
there `w = (2 sigma - 1)/(2 m sigma)` is `3/32` and exactly representable. At
`m = 11` it is `3/44`, and the rounding of the constant would be charged to the
window. The B-spline tables have no such constant and run at the default 11.

The remaining rounding the tests do carry is `KPI` and, for
`sincpow_phi_hut`, the stored `1/w`. Both are noted at the tests that carry
them.
