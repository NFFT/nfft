# ADR-0005: Fractional minimax window correction with shipped optima

## Status
Accepted (2026-11-06).

## Context
`OPT_WINDOW` corrects the polynomial Kaiser-Bessel window at plan time so
that the transform's relative error `E(xi, t) = (D - C) / C` is small over the
band and the node offsets, with `C` the window's own `psi_hat`. On a design of
points this is `min_c max_i |g_i c| / (b_i c)`, a ratio of linear forms. Its
level sets are intersections of second-order cones, so the problem is
quasiconvex.

The Lawson loop (iteratively reweighted least squares) froze the denominator
at the current window. The frozen problem is a different problem: at sigma
1.5, library m 5, deg 8 its optimum is 2.66e-9, the window that attains it has
a true ratio of 6.7e-8, and the true optimum is 1.57e-8. Gap of the loop to the
optimum: median 4 percent, worst 95 percent.

The optimum depends on `(sigma, m, h)` only, not on `N`, `d`, the nodes or the
format. From scratch it costs, in double, a median of 74.7 ms and at most
277 ms on the `kbfrac` cell set (sigma 1.25 to 2.1, m 3 to 6) and 430 to
473 ms on the Appendix A rows at sigma 1.125. In software binary128 the
trafo test did not finish its m 16 plan from scratch in 24 min, and one
m 21 generator cell needs hours.

## Decision

- **Dense cone solver in `E` replaces IRLS.** `Y(socp3_solve)` (primal-dual
  interior point, Nesterov-Todd scaling, Mehrotra predictor-corrector) solves
  the level-set subproblem. A safeguarded Dinkelbach driver with bisection
  brackets the optimum, and an exchange step adds the fine-grid maxima to the
  design.
- **Constants are multiples of `E_EPSILON` or dimensionless.** No literal below
  `1e-4` in the solver path of the correction or in `socp3.c`; the one exception
  is the `1e-9` relative tolerance of the exact-hit match `same_sigma`, which
  compares sigma values, not errors. `E` is never assumed wider than `R`: the
  double and long double builds run with `E = R` and are the test matrix; the
  float build has `E = double`.
- **The guard is keyed to `R`.** `Y(kb_poly_correct_ok)` reads
  `KB_POLY_BR_LS_<MANT_DIG>`, so the float build optimizes only cells whose
  plain error is 100x above the float floor.
- **Hold-out acceptance at 0.9 of plain** (`KB_POLY_CORR_ACCEPT`). The
  hold-out evaluates the `R`-rounded window with the hut table as the transform
  does. At 0.8 it rejected the plan's own reference optimum (0.836x of plain at
  sigma 1.5, library m 5, deg 8).
- **One shared table from a deterministic Python generator.**
  `gen/kbpolygen/optima.py` runs the same algorithm in 200-bit `python-flint`
  ball arithmetic and writes `kernel/util/kbpoly_opt_data.h`: per
  `(sigma, m, h)` the 36-digit table, its error `e` and up to 400 exchange
  points. One cell per library degree of each format at each `(sigma, m)`.
  Each build compiles the rows up to its own m cap.
- **Exchange gate 0.5 percent in the generator.** 80 points per round, up to
  40 rounds, until `e <= 1.005 lo`; a cell that misses it fails the run. The
  library's exchange runs 40 points times 10 rounds to 2 percent.
- **The table is used on exact `(sigma, degree)` hits only.** The plan copies
  the cell, computes the hut series and runs the hold-out; no solve. Off the
  table it solves from scratch. The nearest-below cell start stays in
  `Y(kb_poly_correct)` but `Y(kb_poly_init)` does not use it: it was 0.87x to
  2.6x of the from-scratch time at equal quality.
- **No test-only library symbols.** The only observable of the correction is
  the return value of `Y(kb_poly_correct)`; costs are measured by wall time in
  scratch harnesses.

## Alternatives rejected

- Pure mpmath generator: 1.4 us per multiply-add, 60 h for the largest cell.
- A C harness as generator: the table would depend on the host's long double
  and libm.
- Gauss-Newton on the ratio: linear convergence, denominators change by order
  one along weak directions.
- Dinkelbach alone from a far start: the steps stall far above the optimum,
  2x to 190x off at sigma 1.0625.
- Weighted L2 ratio by a generalized eigenproblem: windows with sign-changing
  spectra.
- Newton on `tau` from the dual slope: the same degeneracy as Dinkelbach.
- Preconditioner truncated at `1e-8` of the largest singular value: loses the
  optimum at large `m`; `500 eps` is right.
- Interpolation or extrapolation of optimal coefficients in sigma: 1.2x to 8x
  worse than either endpoint; a window for a larger sigma is not designed for
  the wider band.
- A rounding-aware objective: the cone program cannot resolve order-one
  differences to the floor while it weights a rounding term over the same
  quantities; 1.1x at best.

## Consequences

- Appendix A of the plan, 24 rows in double: 0.989x to 1.017x of the published
  optimum.
- An exact hit costs 0.34 to 0.40 ms per axis in double. From scratch: see
  Context; the median is 24 percent over the plan's 60 ms estimate
  (issue 16).
- The header is about 0.9 MB of source. It is regenerated, never hand-edited;
  the `--mmax 11 -j 3` run takes about 3 h on 3 cores.
- The 113-bit cells `m >= 12` are missing (issue 14); 113-bit plans there solve
  from scratch.
- In double at small sigma and large m the hut table's rounding at the band
  edge, not the window, limits the gain (issue 15).
- Near the rounding floor of an `E = R` format the driver cannot resolve
  levels (issue 13); there only the table helps.
