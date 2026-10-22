# Polynomial Kaiser-Bessel degree table

Produces `kernel/util/kbpoly_data.h`, the table `Y(kb_poly_degree)`
(`kernel/util/kbpoly.c`) reads to choose the degree of the polynomial
Kaiser-Bessel window, for every floating-point format the library builds for:

```bash
uv run --with mpmath==1.3.0 python -m tests.kbpolygen.generate
```

About 50 s on 10 cores. Tests:

```bash
uv run --with mpmath==1.3.0 --with pytest python -m pytest tests/kbpolygen/tests -q
```

## The rule

The plan fits one polynomial per unit cell of the window. Truncating its
Chebyshev series after degree `d` leaves an error that is smooth in the offset,
so it adds up over the `2m + 2` taps, unlike rounding noise. The table holds
the lowest `d` with

    2 sum_{j>d} max_cells |c_j| <= max(0.05 pi e^(-b R), eps) * peak

`R = m + 1` the reach, `b = pi (2 - 1/sigma)`.

- `0.05 pi e^(-b R)` keeps the transform within 5% of its own worst-case error
  `E = c / I0(2 pi R sqrt(1 - 1/sigma))`. Certified and sampled values put `c`
  in [1.10, 2.81] over sigma 1.1 .. 4; the rule takes `c = 1`.
- `eps` keeps the polynomial below the transform's rounding noise where `E` is
  smaller still.

The plan raises an odd `d` by one: at even degree the leading interpolation
error is odd in the cell and has no mean, which would otherwise add up over
the taps (`KB_POLY_ROWS` in `include/infft.h`).

The coefficients are computed in mpmath, so the table is free of the rounding
noise a fit in working precision sees; in float that noise moved the degree by
one at random.

## Layout

Per format and `m`: the degree at sigma 1.1 and the sigma thresholds, in
units of 1e-4, from each of which the degree is one higher. The generator
stops if the degree ever falls or jumps in sigma. Thresholds are rounded
down, so a lookup never gives less than the exact degree.

## Rows stop at m_opt

Past `m_opt`, larger `m` no longer lowers the transform error: truncation
`c / I0(a R)` meets the round-off floor `u I0(b R) / I0(a R)` where
`I0(b R) = c / u`, which gives one constant per format,

    KB_POLY_BR_OPT = L + ln(2 pi L) / 2,  L = ln(2.81 / u),
    m_opt = ceil(KB_POLY_BR_OPT / b) - 1.

Measured transform errors reach their floor at `m_opt` or one below it. Past
`m_opt` the degree does not rise, so the lookup reads row `min(m, m_opt)`; a
plan with `m > m_opt` gets the degree of `m_opt`, up to 4 more than it needs.

Double goes last, as the `#else`, so a format the table does not know about
builds against the double table.
