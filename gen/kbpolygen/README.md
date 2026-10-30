# Polynomial Kaiser-Bessel degree table

Produces `kernel/util/kbpoly_data.h`, the table `Y(kb_poly_degree)`
(`kernel/util/kbpoly.c`) reads to choose the degree of the polynomial
Kaiser-Bessel window, for every floating-point format the library builds for:

```bash
uv run --with mpmath==1.3.0 --with numpy python -m gen.kbpolygen.generate
```

About 50 s on 10 cores. Tests:

```bash
uv run --with mpmath==1.3.0 --with numpy --with pytest python -m pytest gen/kbpolygen/tests -q
```

## The rule

The plan fits one polynomial per unit cell of |nx| and the transform sees
it only through the 2m + 2 taps at a node's offset and the deconvolution
`C(xi)`. The relative error the interpolation adds to the transform is

    poly_part = max_{xi, t} |D_e(xi, t)| / C(xi),

`D_e` the tap sum of the interpolation error (`functional.py`). The table
holds, per format, the lowest half degree `h` (degree `2h`, even: at even
degree the leading interpolation error is odd in the cell and has no mean)
with

    poly_part <= max(0.05 * 1.10 / I0(a R), 2 eps),   a = 2 pi sqrt(1 - 1/sigma), R = m + 1.

- `1.10 / I0(a R)` is the transform's own worst-case error `c / I0(a R)` at
  the certified lower end `c = 1.10` of `c` over sigma 1.1 .. 4, so the
  polynomial adds at most 5% of it.
- `2 eps` keeps the polynomial below the transform's rounding noise where
  the error is smaller still.

The interpolation error coefficients are exact (mpmath); the functional's
sum, which cancels nothing, runs in double from them and is accurate for
every format. Compared with the earlier rule, which held the Chebyshev tail
summed over the cells to 5% of `pi e^(-bR)`, this measures the quantity the
transform sees and uses the actual error as the reference; it is one to
two `h` lower in most of the table.

## The correction guard

`KB_POLY_BR_LS_<mant>` is the `b R` up to which the optional plan-time
correction of `Y(kb_poly_correct)` pays in a type of that mantissa width:
the truncation error `1.10 / I0(a R)` is still 30 times the round-off floor
`u I0(b R) / I0(a R)`, the `m_opt` rule below with a margin of 30. Past it
the correction would fit rounding noise; `Y(kb_poly_init)` skips it there.

## Layout

Per format and `m`: `h` at sigma 1.1 and the sigma thresholds, in units of
1e-4, from each of which `h` is one higher (a step of two repeats the
threshold). The generator stops if `h` ever falls in sigma. Thresholds are
rounded down, so a lookup never gives less than the exact `h`.

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
