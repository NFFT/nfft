# ADR-0003: 128-bit quadruple precision readiness

## Status
Accepted (data-pipeline readiness); the `__float128` build is scoped, not implemented.

## Context
We want the test reference-data pipeline ready for 128-bit quadruple precision and,
later, even wider types. The reference-data generator (`tests/refgen/`, see
[ADR-0002](0002-python-reference-data-generator.md)) now emits ≥64-digit values and
takes `--precision`, so the **data** is quad-ready. Wiring an actual `__float128`
*build* is a separate, larger effort with these touch points:

- `configure.ac`: a new precision (`NFFT_QUAD`, `PRECISION=q`, `PREC_SUFFIX=q`,
  `NFFT_PRECISION_MACRO=NFFT_PRECISION_QUAD`), conflict checks vs single/long-double,
  and `-lquadmath` detection. `MANT_DIG==113` is already accepted at `configure.ac`.
- FFTW: `m4/nfft_lib_fftw3.m4` is fully parameterized by `$PREC_SUFFIX`, so it finds
  `fftw3q`/`fftwq_*` automatically once `PREC_SUFFIX=q` — plus libquadmath in the link.
- `include/infft.h`: an `#elif defined(NFFT_QUAD)` arm for `R`/`C`, the `K(x)`
  `Q`-suffix branch, the `EPSILON`/`MANT_DIG` arm (`FLT128_*`), and the full quad
  math-function macro set (`fabsq`/`sqrtq`/`expq`/`cosq`/`csqrtq`/…). This macro set
  is the largest mechanical change.
- I/O: `printf`/`scanf` cannot handle `__float128`. The `__FI__`-based `fscanf`
  (`tests/*.c`) and `__FE__`-based prints must switch to `strtoflt128` /
  `quadmath_snprintf` (`%Qe`). The current `MANT_DIG==113 → "...LE"` format arms are
  written for 113-bit *long double*, not `__float128`, and would misbehave for true
  quad.
- Public API: `include/nfft3.h` needs `NFFT_MANGLE_QUAD` + a 4th `*_DEFINE_API`
  instantiation per module; `include/nfft3mp.h` needs an `NFFT_PRECISION_QUAD` arm.
- Test bounds: the `MANT_DIG==113`/`64` fast-bound `a`/`b` tables in
  `tests/{nfft,nfct,nfst}.c` carry `// TODO`s and need tuning once a quad build runs.

## Decision
Make the reference-data pipeline quad-ready now (done: generator `--precision`,
≥64-digit files, precision-portable format). Defer the `__float128` build to a
dedicated follow-on, tracked against the touch-point list above.

## Verification evidence
The regenerated data was validated across the full **4 windows
(kaiserbessel/gaussian/bspline/sinc) × 3 precisions (float/double/long-double)**
matrix — the three `MANT_DIG` bound branches (24/53/64) every supported window can
hit. (`dirac` has no CUnit suite: `#error Unsupported window function`.) Any fast-bound
`a`/`b` re-tunes made for the new draw are listed here:

```
PASS  kaiserbessel  float
PASS  kaiserbessel  double
PASS  kaiserbessel  long-double
PASS  gaussian  float
PASS  gaussian  double
PASS  gaussian  long-double
PASS  bspline  float
PASS  bspline  double
PASS  bspline  long-double
PASS  sinc  float
PASS  sinc  double
PASS  sinc  long-double
```

All 12 cells PASSED with the regenerated reference data. **No fast-bound `a`/`b`
constants were re-tuned** — the existing tolerances in `tests/{nfft,nfct,nfst}.c`
covered the new random draw unchanged, confirming the bounds are draw-independent as
designed.

The `MANT_DIG==113` (quad) branch cannot be exercised until the `__float128` build
above exists; its `a`/`b` constants remain `// TODO` and must be tuned then.

## Consequences
- Existing double/long-double trees are unaffected.
- A future quad build reuses the same `tests/data/*.txt` files (the C reader rounds
  on input); only the build plumbing and I/O path above must be added.
- The `__FI__`/`__FE__` `__float128` I/O work is the genuine blocker for *running*
  the suite at quad and must be done before the bound tables can be tuned.
