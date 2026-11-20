# Polynomial sinc-power degree table

Produces `kernel/util/sincpowpoly_data.h`, the table `Y(sincpow_poly_degree)`
(`kernel/util/kbpoly.c`) reads to choose the degree of the polynomial
sinc-power window, for every floating-point format the library builds for:

```bash
uv run --with mpmath==1.3.0 python -m gen.sincpowpolygen.generate
```

Tests:

```bash
uv run --with mpmath==1.3.0 --with pytest python -m pytest gen/sincpowpolygen/tests -q
```

## The rule

The window is `phi(nu) = w sinc(pi w nu)^(2m)`, `w = (2 sigma - 1)/(2 m sigma)`,
with peak `w`. The plan fits one polynomial per unit cell of `|nu|`, the
`m + 2` cells `[m + 1 - c, m + 2 - c)`. The table holds the lowest `d` with

    2 sum_{j>d} max_cells |c_j| <= eps * w

`c_j` the exact Chebyshev coefficients of the window on a cell. There is no
cut at the transform error `E`: deconvolution amplifies the window error at
the band edge, most at low sigma, so a cut at `E` loses accuracy there.

The plan raises an odd `d` by one (`CELL_POLY_ROWS` in `include/infft.h`).

## Layout

Per format and `m`: the degree at sigma 1.1 and the sigma thresholds, in
units of 1e-4, from each of which the degree is one higher. The generator
stops if the degree ever falls or jumps in sigma. Thresholds are rounded
down, so a lookup never gives less than the exact degree.

Rows stop at `m = 32`. The degree falls with `m`; the generator checks at
`m = 64` that no point needs more than row 32, and the lookup reads row
`min(m, 32)`.

Double goes last, as the `#else`, so a format the table does not know about
builds against the double table.
