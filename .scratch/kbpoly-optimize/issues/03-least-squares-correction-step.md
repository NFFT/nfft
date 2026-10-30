# 03: Plan-time least-squares correction of the polynomial window

Status: done (feature/kbpoly-optimize)

Part of: PRD G2, R4, R5, R7, P1, D1, D2, D6. Depends on: issue 04 (for `E`
and the guard) can be developed in parallel but must land first. Background:
`analysis.md` sections 4, 5.1, 10; prototype `kbopt4.py` function
`fit_residual`.

## Task

In `kernel/util/kbpoly.c`, after `cheb_of` and before `cheb_to_monomial` in
`Y(kb_poly_fit)` (or as a new function `Y(kb_poly_correct)` called from
`Y(kb_poly_init)`), add the step:

1. Sample grid: `n_xi = 40` points uniform in `[0, 1/(2 sigma)]` (evenness of
   the window gives the negative half), `n_t = 32` midpoints in `[0, 1)`.
   `sigma` follows from `b` as in `Y(kb_poly_degree)`.
2. Residual: for each sample, `r = D_cheb(xi, t) / C(xi) - 1` with `D_cheb`
   evaluated from the Chebyshev coefficients of the stencil cells
   `-R .. R-1` (table columns `1 .. 2m + 2`) and `C` from `Y(kb_phi_hut)`
   with the plan's constants. Form `D_cheb`, `C` and `r` in `E`.
3. Unknowns: Chebyshev coefficients `delta[a][j]` of the correction, `a` over
   the stencil cells, `j = 0 .. d`, with evenness `delta[-1-a][j] =
   (-1)^j delta[a][j]` so that `(m + 1)(d + 1)` unknowns remain.
4. System: rows are the real and imaginary parts of
   `T_j(2t - 1) e^(-2 pi i xi (a + t)) / C(xi)` (with the mirrored cell folded
   in), right-hand side `-r`; append `sqrt(lambda)` times the identity for the
   Tikhonov term with `lambda` chosen so that it is negligible at the table
   degree and dominant one degree above it (start from `lambda = 1e-2` times
   the mean squared row norm and tune against the prototype). Solve by
   Householder QR in `E`; no normal equations.
5. Add `delta` to the Chebyshev coefficients, round to `R`, and continue with
   `cheb_to_monomial` and the refinement step as today.
6. Take the step only when the plan flag of PRD R5 is set and the guard of
   issue 04 admits `(b, m)`; otherwise the function is unchanged.

Memory: one `(2 n_xi n_t + n_u) x n_u` matrix in `E`, `n_u <= 13 * 17` in
double; allocate with `Y(malloc)` once per fit and free it. Cost target: under
2 ms per axis at `m = 11, d = 12` in double on the dev container; record the
measurement in this issue.

Flag: take bit 14 (`1U<<14`) in `include/nfft3.h`, name per PRD R5, document
next to `ANALYTIC_WINDOW`. Reject the combination with `ANALYTIC_WINDOW` and
on non-Kaiser-Bessel builds in each module's `X(check)`, the way
`ANALYTIC_WINDOW` and `PRE_POLY_PSI` were handled on this branch (see the
commit that added `PRE_POLY_PSI`).

## Acceptance

- With the flag off, `kernel/util/kbpoly.c` produces bit-identical tables to
  today's (test: compare against a table built with the flag off before the
  change, or keep the old path and diff).
- With the flag on, in double at `sigma = 2, m = 2, d = 6`, the transform
  error against the direct transform on the `tests/nfft.c` online checks drops
  by at least 1.3x in RMS over the existing `poly_trafo` cases, and no case's
  `err` exceeds its `bound`.
- `tests/kbpoly.c` gains a check (see issue 06) that the corrected window's
  measured `J_inf`, computed in the test from the coefficient table by the
  formula of section 3, is at most that of the uncorrected window for every
  `(sigma, m)` the guard admits, in every precision.
- Float build with `E == R`: the same check passes and the guard rejects every
  `(sigma, m)` where the analysis's float32 table shows a loss (section 10).
- No run-time file changes except the flag definition and `X(check)`.
