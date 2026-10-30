# PRD: Fit the polynomial Kaiser-Bessel window to the transform error

Status: implemented on `feature/kbpoly-optimize` (issues 01 to 04, 06 and 09 in one pass; 05 folded into 03 as the bounded reweighting; 07 and 08 open). Measurements in analysis section 15.1. Background and measurements: `analysis.md` in this
directory (sections referenced below as A.n).

## Problem

The branch `feature/kb-poly-window` evaluates the Kaiser-Bessel window from one
polynomial per unit cell, fitted at plan time by Chebyshev interpolation of the
window function, with the degree read from a generated table. Two things are
left on the table:

1. The degree table holds the polynomial to 5 percent of `pi e^(-bR)`, a lower
   bound of the transform error that sits 10x to 100x below the actual
   Kaiser-Bessel error. The polynomial is therefore held to 0.05 to 0.5 percent
   of the error it is meant to be negligible against, which costs 0 to 3 degrees
   of Horner work per tap (A.9).
2. The polynomial reproduces the Kaiser-Bessel function, which is itself only
   near-optimal for the 2m+2-tap stencil. Fitted to the transform's own error
   functional instead, the same degree gives 1.3x to 3.1x lower RMS error
   (least squares) or 3.4x to 8x lower worst-case error (minimax), at zero
   run-time cost (A.6).

## Goals

- G1. Regenerate `kernel/util/kbpoly_data.h` from the transform's error
  functional, so the degree is the smallest at which the polynomial's
  contribution to the transform error is within 5 percent of the actual
  Kaiser-Bessel error.
- G2. Add an optional plan-time correction step after the Chebyshev
  interpolation in `Y(kb_poly_fit)` that fits the difference to the interpolant
  by weighted least squares against the error functional.
- G3. Add an optional Lawson (iteratively reweighted least squares) step for the
  worst-case objective on top of G2.
- G4. Keep every run-time path unchanged: the coefficient table has the same
  layout, `kb_poly_run` and `kb_poly_phi` do not change, `PHI_HUT` does not
  change.
- G5. Keep the float, double and long-double builds working, with the fit in
  the build's own precision; use the local type `E` where it is wider, without
  depending on it.

## Non-goals

- Re-optimizing the deconvolution factors `PHI_HUT` (measured: no gain, A.5.4).
- Changing the reach `R = m + 1`, the stencil, or the `e^(-bR)` rate.
- A coefficient table in the header (table-size and `sigma`-grid cost, A.5).
- Node-adaptive weighting of the fit (A.5.7); a later feature on the same
  functional.
- Windows other than Kaiser-Bessel.

## Requirements

### Functional

- R1. The generator (`tests/kbpolygen`) gains an exact evaluation of the error
  functional of a grid-aligned piecewise polynomial window (A.3): `D(xi, t)`,
  `C(xi)`, `J_inf`, `J_2`, and the polynomial part `max |D_poly - D_KB| / |C|`,
  for a given `(sigma, m, d)` and format. It must reproduce the Kaiser-Bessel
  constant `c` in `[1.10, 2.81]` for the truncated window (A.3) and the
  prototype's figures in `results/` to three digits.
- R2. The degree rule becomes: smallest `d` with polynomial part `<= 0.05 x
  J_inf(truncated KB)` (the "strict" criterion, A.9), or below the format's
  rounding floor as today. The measured criterion (transform `J_inf` within 5
  percent of Kaiser-Bessel's) is computed and reported by the generator but
  does not set the table, because it certifies the NFFT functional only (A.9).
- R3. The generator emits, per format, a regime threshold for the correction
  step: the step is taken only when the actual Kaiser-Bessel error is at least
  `K eps` of the type the fit runs in, `K = 30` for least squares and `K = 100`
  for Lawson (A.10). The threshold is stored as a bound on `bR` per format and
  `sigma` range, in the style of `KB_POLY_BR_OPT`, and the code selects the row
  by the mantissa of `E`.
- R4. `Y(kb_poly_fit)` gets an optional correction step: unknowns are the
  Chebyshev coefficients of the correction per cell for the stencil cells
  `-R .. R-1` (the two outer columns of `KB_POLY_COLS(m)` keep their Chebyshev
  interpolant), evenness halves the unknowns, the residual `D_cheb / C - 1`
  and the solve are formed in `E`, the solve is a Householder QR of the
  stacked real and imaginary rows on a `40 xi x 32 t` sample grid, with
  Tikhonov regularization toward zero correction, and the correction is added
  to the interpolant and rounded to `R` before the existing monomial
  conversion and refinement step (A.5.1, A.10).
- R5. A plan flag selects the step. Bit 14 of the plan flags is free; the name
  is to be chosen (proposal: `OPTIMIZED_WINDOW`). The step is opt-in at first;
  flipping the default is a separate decision after the accuracy report shows
  the change (see Open questions). The flag is rejected by `X(check)` with
  `ANALYTIC_WINDOW` and on builds whose window is not Kaiser-Bessel, like
  `PRE_POLY_PSI` was.
- R6. The Lawson step (G3) is a second level of the same flag or a second flag
  (to be decided in its issue), runs on a `120 xi x 96 t` grid for at most 20
  iterations keeping the best iterate, reuses the QR of R4, and is guarded by
  its own threshold (`K = 100`).
- R7. Every module that calls `WINDOW_HELP_POLY_INIT` (NFFT, NFCT, NFST, NNFFT)
  gets the same behavior, since all anchor the run below the node and so share
  the stencil `-R .. R-1`. Any module that reads the outer columns
  (`PRE_LIN_PSI`) reads the uncorrected interpolant there.

### Precision

- P1. All plan-time arithmetic is in `R` or `E`; no `double` literals or types
  in the fit. With `E == R` the behavior is the float32 picture of A.10 and
  must pass the tests.
- P2. The correction must never make the transform error worse than the
  uncorrected polynomial's by more than 5 percent on any `(sigma, m)` the
  table admits; the regime threshold of R3 is what guarantees it, and the
  tests check it.

### Performance

- F1. No change to any CodSpeed benchmark name or to run-time instruction
  counts of the transform (the benchmarks measure the transform, not the
  plan).
- F2. Plan-time cost of the least-squares step at `m = 11, d = 12` in double
  under 2 ms per axis, of the Lawson step under 20 ms per axis, measured and
  recorded in the issue.

### Compatibility

- C1. Autotools and CMake builds both build and test the change (Coexistence
  policy, `CONTEXT.md`).
- C2. The Julia and MATLAB/Octave interfaces do not mirror `ANALYTIC_WINDOW`
  today; the new flag is likewise not exposed until asked for.
- C3. The generated header carries the `do not edit` banner and the generator's
  tests cover the new rule.

## Design decisions (taken)

- D1. Fit the correction, not the window, so that only the residual needs
  precision (A.10).
- D2. Chebyshev basis for the solve, monomials only for the run time, through
  the existing conversion and refinement (A.5.1).
- D3. Least squares before minimax: deterministic, closed form, and the
  worst-case gain of minimax needs a dense grid to be real (A.6.2).
- D4. Strict degree criterion for the table (R2), because it carries over to
  every consumer of the table (A.9).
- D5. The guard lives in the generated table, not in a plan-time estimate,
  because a quantity computed at the noise floor cannot judge its own noise
  (A.10).
- D6. `E` for the residual and the solve, `R` for the stored coefficients.

## Open questions

- Q1. Default on or off for the correction step once it exists. Proposal: off
  until one release of accuracy-report data; the analysis measured no case
  where the guarded step makes the error worse, but the measurement is six
  `(sigma, m)` pairs.
- Q2. Whether `E` exists in `include/infft.h` on this branch. `CONVENTIONS.md`
  lists it; no definition was found under `include/`. Issue 04 settles it.
- Q3. Whether the NFCT and NFST transform error with the corrected window
  should also be measured by the generator (their functional combines
  `E(xi, t)` with the reflected offset). Not required for correctness, since
  the strict criterion bounds the polynomial part for every consumer, but it
  would show how much of the gain they get.

## Milestones

1. Generator: error functional (issue 01), new degree rule and table (02).
2. Library: `E`, flag, least-squares correction and guard (03, 04).
3. Lawson step (05).
4. Tests, measurements, docs (06, 07, 08).

## Success metrics

- The regenerated table lowers the double-precision degree by at least one at
  `sigma = 2, m = 2` and `sigma = 2, m = 3` (A.9: 7 -> 6, 9 -> 8) with the
  transform tests passing at the same bounds.
- With the correction step on, the accuracy report's `nfft/.../fast/...`
  metrics improve by at least 1.3x in `max err` for the Kaiser-Bessel cells in
  double, and no cell regresses by more than 5 percent.
- With `E == R` (float build) the same holds for the cases the guard admits,
  and the cases it rejects are unchanged.

## Addendum (2026-10-30): the rounding-aware objective

Analysis section 14 tried goal 4 (aliasing plus rounding amplification in one
objective). It is numerically out of reach for double-precision solvers and,
where a well-conditioned substitute could be measured, buys about 1.1x at the
band edge at the float floor. It is dropped from the goals. Two things from
that work stay in scope: evaluating `PHI_HUT` in `E` in the float build
(issue 09, measured 1.17x band edge and 1.3x to 1.5x RMS at the float floor),
and the note that `kbgain2.c` is the harness to use for floor measurements,
since a `cexp` of a large argument in the reference is itself a floor.
