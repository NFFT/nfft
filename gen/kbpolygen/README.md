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
optimization of `Y(kb_poly_correct)` can pay in a type of that mantissa
width: the truncation error `1.10 / I0(a R)` is still 100 times the
round-off floor `u I0(b R) / I0(a R)`, the `m_opt` rule below with a margin
of 100. Past it the optimization would fit rounding noise;
`Y(kb_poly_init)` skips it there. The margin is a coarse first filter only:
the optimized window lowers the error up to 100 times, so near the floor
it can still fit noise and come out worse, which `Y(kb_poly_correct)`
catches itself by validating the result on a hold-out sample set and
keeping the plain table where the optimized one does not beat it.

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

# Optimized window table

`optima.py` produces `kernel/util/kbpoly_opt_data.h`, the polynomial window
of `Y(kb_poly_init)` with the `OPT_WINDOW` correction of
`Y(kb_poly_correct)` already applied, for sigma in {1.1, 1.25, 1.5, 1.75,
2, 2.5, 3, 4}, `m` from 1 to the `KB_POLY_BR_LS_113` guard and `h` of the
113-bit degree table:

```bash
uv run --with python-flint==0.9.0 --with numpy --with mpmath==1.3.0 python -m gen.kbpolygen.optima --mmax 11 -j 5
```

It runs the algorithm of `kbpoly.c` (design,
first least-squares pass, Dinkelbach steps and bisection on cone
subproblems, exchange on the 401 x 65 fine grid) in 200-bit ball arithmetic
of `python-flint==0.9.0`, which is deterministic for a pinned version, and
solves to 0.1% instead of 1%. The initial design is that of `kbpoly.c` for
`h <= 8` and `m <= 9`; above, it has `2 h` offsets and, for `m > 9`,
`3 (m + 1)` frequencies, since the 17 x 24 design of `kbpoly.c` is too
coarse for the 113-bit cells. Per cell the header holds the monomial table,
`e`, the fine-grid maximum of `|E|` with local refinement, and up to 400
design points as `(round(65535 xi / ximax), round(65535 t))`. The script
re-evaluates `e` from the emitted 36 digits and stops if it moved by more
than 0.1%. `--cell SIGMA M H` solves one cell and writes nothing, `--mmax M`
writes a partial header with every `KB_POLY_OPT_MMAX` capped at `M`.
The committed header is such a partial header, `--mmax 11` (88 cells, 56
minutes with `-j 5`, about 0.7 GB per worker): it covers float, double and
80-bit long double, but not the 113-bit cells `m = 12 .. 21`, which take
hours per cell. The exchange stops after 10 rounds; the log gives per cell
the lower bound `lo` of the last design problem next to `e`. The log of the committed run is
`.scratch/kbpoly-optimize/results/optima_log.txt`.
