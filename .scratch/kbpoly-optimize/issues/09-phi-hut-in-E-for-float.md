# 09: Evaluate the deconvolution factors to one ulp in `R`

Status: ready-for-agent

Part of: PRD addendum. Background: `analysis.md` sections 14 and 15.1.

## Task

The deconvolution factors `PHI_HUT` (`Y(kb_phi_hut)`: an exponential, the
scaled Bessel function `Y(bessel_i0_exp_scaled)`, a reciprocal) carry an
error of a few ulp that lands on the result unattenuated. Whether that shows
depends on the format's rounding floor, not on the format: in double and
long double the floor is 8 to 30 ulp from the FFT and the gather and the
factors do not register; in float the floor is about 14 ulp at the band
edge and exact factors lower the error 1.17x there and 1.3x to 1.5x in RMS
at `sigma = 2, m = 3, 4` (`results/floor_float_s2.txt`, "ownc" rows).

Forming the exponential in a wider type was tried on this branch and
reverted: it recovered almost nothing, the remaining error is in the float
Bessel evaluation, and a wider type is a special case that only the float
build can have. The precision-agnostic statement of the goal is an accuracy
target per format: evaluate `PHI_HUT` to about one ulp in `R`.

- Measure the ulp error of `Y(bessel_i0_exp_scaled)` and of `Y(kb_phi_hut)`
  per format against mpmath over the arguments the band uses (`m ra` from
  `m sqrt(b^2 - pi^2/sigma^2)` to `m b`); `tests/bessel.c` has the harness.
- Tighten the float fit in `gen/besselgen` (more terms, or a split that
  keeps the argument reduction exact) until `kb_phi_hut` is within one ulp
  in float; check that double and long double already are, and tighten them
  if not.
- Keep the evaluation in `R` and the code the same in every format.

## Acceptance

- Float build, `kbgain2.c`-style measurement at `sigma = 2, m = 3, 4`: the
  band-edge error and RMS within 5 percent of the "ownc" rows of
  `results/floor_float_s2.txt`.
- Double and long double: `make check` output unchanged in the digits the
  tests print.
- No change to the run-time instruction count of the transform; the factors
  are `O(N)` per plan or per call.
