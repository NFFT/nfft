# Bessel I0 coefficient generator

Produces `kernel/util/bessel_i0_data.h`, the coefficient tables for the modified
Bessel function I0 in every floating-point format the library builds for. This
replaces the Mathematica notebook (`kernel/util/bessel_i0.nb` / `.m`) that used
to be the only way to regenerate them.

## Dependencies (dev only; not a build dependency)

Run via [uv](https://docs.astral.sh/uv/), the same way as `tests/refgen`. The
only dependency, `mpmath`, is declared inline on the command; there is no
`pyproject.toml` or lockfile to maintain.

## Regenerate

```bash
# from the repository root
uv run --with mpmath==1.3.0 python -m tests.besselgen.generate
```

Deterministic: same `SPEC` in, same bytes out. Nothing here reads a clock or a
random seed. Takes a few minutes, mostly in the quadruple-precision fits.

## Run the generator's own tests

```bash
uv run --with mpmath==1.3.0 --with pytest python -m pytest tests/besselgen/tests -q
```

## Options

- `--out P` output header (default `kernel/util/bessel_i0_data.h`)
- `--precisions L` comma-separated `MANT_DIG` values (default `24,53,64,113`).
  A subset writes an incomplete `#if` chain, so use it only for experiments.
- `--dps N` mpmath working digits (default 60, or 90 for quadruple)
- `--search` re-derive the polynomial degrees from the accuracy target instead
  of using the pinned ones, and report them. Use this after changing a split,
  then move the reported degrees into `SPEC` deliberately rather than letting
  the emitted header drift.

## What is fitted

Two ranges per format, split at `NFFT_I0_ASYMP_SPLIT`:

| range | form | evaluated by |
|---|---|---|
| `x <= split` | `I0(x) = 1 + y*P1(y)`, `y = (x/2)^2` | Horner in `y` |
| `x > split` | `I0(x) = exp(x)/sqrt(x) * P2(t)`, `t = 2*split/x - 1` | Clenshaw in `t` |

Both are **weighted** minimax fits (`remez.py`), the weight chosen so the
quantity that equioscillates is the *relative* error of the final expression.
That is what a floating-point kernel needs; an unweighted fit would spend its
accuracy where `I0` is largest.

The two forms are chosen for conditioning rather than for the smallest degree:

- `1 + y*P1(y)` keeps every `P1` coefficient positive, so the Horner chain has
  no cancellation, and the leading `1` makes the form exact at `x = 0`. `P1`
  starts at the Taylor coefficients `1, 1/4, 1/36, ...` with a small minimax
  correction on top.
- `P2` is fitted against `sqrt(x)*exp(-x)*I0(x)`, which is smooth and bounded
  away from zero on `(0, 1/split]`, and the argument is mapped into `[-1, 1)`
  so Clenshaw stays in its stable domain.

`SPEC` pins the split and the two degrees per format. `docs/bessel-i0-accuracy.md`
records the sweep those values came from and the resulting ULP measurements.

## Notes for anyone changing this

- **Report the error of the form you ship.** The fit runs in the Chebyshev
  basis, but branch 1 ships monomial coefficients evaluated by Horner.
  `scheme.verify_branch1` / `verify_branch2` re-measure the shipped form in
  exact arithmetic, and that is what the generator prints. Reporting the fit's
  own residual instead once hid a basis-conversion error floor that was
  harmless at double and fatal at quadruple precision.
- The branch-1 weight vanishes at `y = 0`, which would pin the equioscillation
  level to zero if a reference point landed there. The search interval steps
  just inside; the Chebyshev basis interval must **not** move with it.
- The equioscillation level that falls out of the linear system is only
  meaningful once the reference has converged. Every error reported by this
  package is a dense-scan measurement instead.
