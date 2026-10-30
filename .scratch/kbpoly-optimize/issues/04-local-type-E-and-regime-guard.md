# 04: The local type `E` and the regime guard table

Status: done (feature/kbpoly-optimize)

Part of: PRD R3, G5, D5, D6, Q2. Depends on: issue 01 (functional) for the
thresholds. Background: `analysis.md` section 10.

## Task

### `E`

`CONVENTIONS.md` lists `E` as "real type for local variables (possibly extra
precision)" and no definition exists under `include/` on this branch. Decide
with the maintainer whether it already exists elsewhere (grep the whole tree,
including `config.h` and the Autotools and CMake configuration); if not,
define it in `include/infft.h`: `long double` when the build is double and
`HAVE_LONG_DOUBLE`-equivalent support is configured, otherwise `R`. Provide
the handful of math wrappers the fit needs in `E` (sqrt, cos, sin, exp, fabs)
in the style of the `K()`/`SQRT`/`COS` macros, and `MANT_DIG_E`. The float
build gets `E = R` unless the maintainer wants `double` there; record the
decision in the issue.

### Guard

In `tests/kbpolygen/generate.py`, for each format `MANT_DIG` in `{24, 53,
64, 113}` compute the largest `bR` (equivalently the largest `m` per `sigma`)
at which the exact Kaiser-Bessel error `c / I0(aR)` is still at least
`K eps`, `K = 30` for the least-squares step (and `K = 100` for the Lawson
step of issue 05), using the functional of issue 01 for `c`. Emit it to
`kernel/util/kbpoly_data.h` as `KB_POLY_BR_LS` and `KB_POLY_BR_MM` per
format, in the style of `KB_POLY_BR_OPT`, together with a short derivation
comment. In `kernel/util/kbpoly.c` add `Y(kb_poly_step_ok)(b, m, which)` that
selects the row by `MANT_DIG_E` (falling back to `R`'s row when `E == R`) and
returns whether `b (m + 1) <= KB_POLY_BR_*`.

## Acceptance

- The table rows reproduce the analysis's float32 verdicts: at `sigma = 2`,
  `m = 2` admitted, `m = 3` rejected for least squares in float; at
  `sigma = 1.5, m = 4` rejected; in double `sigma = 2, m = 5` admitted.
- Every build (float, double, long double; gcc and clang; Autotools and CMake)
  compiles with `E` defined, and a build with `E` forced to `R` (a `-D`
  override documented in the issue) also compiles and passes `make check`.
- A generator test pins the four verdicts above.
