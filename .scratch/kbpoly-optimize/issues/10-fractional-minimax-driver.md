# 10: Fractional minimax driver with a dense cone solver

Status: done (feature/kbpoly-optimize, c37bfec8 to 272d85c8)

Part of: PRD addendum (2026-11-06). Replaces issue 05. Background:
`analysis.md` section 19.

## Task

Replace the Lawson loop of `Y(kb_poly_correct)` by a solver of the discrete
ratio problem `min_c max_i |g_i c| / (b_i c)`:

- `kernel/util/socp3.c`, `Y(socp3_solve)`: dense primal-dual interior-point
  method for `M` cones `Q^3`, Nesterov-Todd scaling, Mehrotra
  predictor-corrector, normal equations by Cholesky with a ridge of
  `100 eps max diag`. Early stop on a sign decision. NaN check by `ilogb`
  (`-ffast-math` folds `isfinite`).
- Pivoted Householder QR in `E` as the preconditioner, truncated at
  `500 eps |R_00|`.
- Design points as a growable list; scale fixed by `sum_i b_i c` over the
  initial design, unknown with the largest `|Bsum|` eliminated.
- Safeguarded Dinkelbach steps (accepted on the true ratio, kept while a step
  gains 1.3x) and bisection at `sqrt(lo hi)`; stop at `hi <= 1.01 lo` or 30
  solves.
- All constants dimensionless or multiples of `E_EPSILON`; the guard
  `Y(kb_poly_correct_ok)` keyed to `R`.

## Result

Appendix A of the plan, 24 rows in double: 0.989x to 1.017x of the reference
optimum after the odd-degree fix (09373f17). Tests: `socp3`, `kbpoly_qr`,
`kbpoly_driver` in all three trees.
