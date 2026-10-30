# 09: Evaluate the deconvolution factors in `E` in the float build

Status: done (feature/kbpoly-optimize)

Part of: PRD addendum. Background: `analysis.md` section 14.

## Task

In the float build the deconvolution factors `PHI_HUT` (`Y(kb_phi_hut)`:
`exp`, `bessel_i0_exp_scaled`, a reciprocal) are evaluated in float and are off
by a few eps, and that error lands on the result unattenuated. Measured on
develop at `sigma = 2, N = 2048`: replacing them by accurate values lowers the
band-edge error 1.17x at `m = 3` and `4` and the RMS error 1.3x to 1.5x
(`results/floor_float_s2.txt`, rows "ownc" against "stock"). In double the
same replacement changes nothing.

Evaluate `Y(kb_phi_hut)` in `E` (issue 04 defines it; `double` in the float
build on every platform the library supports) and round once to `R`, both on
the `PRE_PHI_HUT` table path (`precompute_phi_hut`) and the per-call path.
`bessel_i0_exp_scaled` then needs an `E` entry point, or the table of
`kernel/util/bessel_i0_data.h` evaluated in `E`; check which the generator in
`gen/besselgen` already provides for the double build and reuse it.

## Acceptance

- Float build, `tests/kbgain2.c`-style measurement (port the driver's
  single-frequency check into `tests/kbpoly.c` or a new test): band-edge
  error at `sigma = 2, m = 3, 4` within 5 percent of the "ownc" rows of
  `results/floor_float_s2.txt`.
- Double and long double builds bit-identical in `make check` output.
- No change to CodSpeed benchmark names; the transform's instruction count is
  unchanged (the change is in precompute or in the per-`k` factor evaluation,
  which is `O(N)`).
