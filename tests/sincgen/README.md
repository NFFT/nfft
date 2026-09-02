# log|sinc| coefficient generator

Produces `kernel/util/sinc_data.h`, the coefficient table behind `Y(log_sinc)`
(`kernel/util/sinc.c`), for every floating-point format the library builds for.
Built on the minimax fitter in `tests/besselgen/remez.py`.

## Why the logarithm

Two windows raise sinc to the `2m`-th power: `PHI_HUT` for the B-spline window
and `PHI` for the sinc-power window. Over the range that carries the weight,
sinc is a number near 1, and rounding it to an `R` throws away exactly the part
that survives the exponentiation. `POW` then multiplies that loss by the
condition number `2m`, so the error of the window value grows linearly in `m`
whatever the rest of the code does.

Returning `log|sinc(x)|` keeps the small quantity small. The caller forms
`EXP(2m * log_sinc(x))`, whose relative error is `2m` times the *absolute*
error of the logarithm; branch 1 below holds that at `|L| * eps` instead of
`eps`. For the B-spline `PHI_HUT` at the default oversampling this is a factor
of about 4 in double precision, and the same factor in single precision, where
the window's own truncation error is far below the evaluation error.

The gain is tied to `|L| < 1`, that is `sinc > 1/e`. Past that the direct
`POW` is the better form, but there the window value is already negligible
(`1e-48` at `m = 11`).

The absolute value is deliberate: the callers' exponent `2m` is even, so
`POW(sinc, 2m)` is positive where sinc is negative, and `PHI` for the
sinc-power window is sampled past the first zero of sinc when sigma is large.

## Dependencies (dev only; not a build dependency)

Run via [uv](https://docs.astral.sh/uv/), the same way as `tests/besselgen`.

## Regenerate

```bash
# from the repository root
uv run --with mpmath==1.3.0 python -m tests.sincgen.generate
```

Deterministic: same `SPEC` in, same bytes out. Nothing here reads a clock or a
random seed.

## Run the generator's own tests

```bash
uv run --with mpmath==1.3.0 --with pytest python -m pytest tests/sincgen/tests -q
```

## Options

- `--out P` output header (default `kernel/util/sinc_data.h`)
- `--precisions L` comma-separated `MANT_DIG` values (default `24,53,64,113`).
  A subset writes an incomplete `#if` chain, so use it only for experiments.
- `--dps N` mpmath working digits (default 60, or 90 for quadruple)
- `--search` re-derive the degrees from the accuracy target instead of using
  the pinned ones, and report them. Use this after changing the split, then
  move the reported degrees into `SPEC` deliberately rather than letting the
  emitted header drift.
- `--split S` / `--n N` override the split or the degree for a single
  precision, for experiments.

## What is fitted

One range per format; the other needs no coefficients.

| range | form | evaluated by |
|---|---|---|
| `\|x\| <= split` | `log\|sinc(x)\| = y * Q(y)`, `y = x*x` | Horner in `y` |
| `\|x\| >  split` | `log\|sinc(x)\| = LOG(FABS(SIN(x)/x))` | directly |

`Q` is a **weighted** minimax fit, the weight chosen so what equioscillates is
the *relative* error of `y*Q(y)`. `fq` is bounded away from zero on the whole
branch, running from `-1/6` at `y = 0` to about `-0.197` at `y = 4`, so unlike
the I0 branch-1 weight this one vanishes nowhere and needs no stepped-in search
interval.

## The split

`NFFT_LOG_SINC_SPLIT` is 2 in every format, and it is a conditioning choice
rather than an accuracy one. The fitted coefficients start at the Taylor series
of `log sinc` (`-1/6, -1/180, -1/2835, -1/37800`), all negative, and at a split
of 2 the Horner growth factor `max sum|c_j| y^j / |sum c_j y^j|` is 1.000 — the
chain cannot cancel. At 2.5 the same degree still fits but the growth factor is
already 53. `scheme.horner_growth` measures it and `--search` prints it.

The split also puts the whole B-spline `PHI_HUT` range, `|t| <= pi/2`, inside
the polynomial branch.

## Pinned degrees

| MANT_DIG | degree | design error | measured runtime error |
|---|---|---|---|
| 24 | 9 | 1.35e-10 | 2.1 ulp |
| 53 | 19 | 9.12e-20 | 1.2 ulp |
| 64 | 22 | 1.72e-22 | not measured (no 80-bit host to hand) |
| 113 | 39 | 7.8e-38 | 1.8 ulp |

## Notes for anyone changing this

- **Report the error of the form you ship.** The fit runs in the Chebyshev
  basis, but the runtime evaluates a monomial Horner chain.
  `scheme.verify_branch` re-measures the shipped form in exact arithmetic, and
  that is what the generator prints.
- Every error reported here is a dense-scan measurement, never the
  equioscillation level the linear system returns.
- Widening the split past 2 buys nothing: the direct branch is already at its
  own floor there, and the polynomial branch starts to cancel.
