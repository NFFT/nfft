# 12: Shipped optima for a fixed sigma set

Status: done (feature/kbpoly-optimize, d4f25f4c to 5d348d64)

Part of: PRD addendum (2026-11-06). Depends on: issues 10, 11. Background:
`analysis.md` section 19, ADR-0005.

## Task

- `gen/kbpolygen/optima.py`: the driver and exchange of issues 10 and 11 in
  200-bit `python-flint` arithmetic. One cell per `(sigma, m, h)` with
  `sigma in {1.1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4}` and `h` the library degree
  of each format whose m cap admits `m`. Exchange 80 points per round, up to
  40 rounds, until `e <= 1.005 lo`. Header `kernel/util/kbpoly_opt_data.h`:
  36-digit monomial tables, `e`, up to 400 exchange points as 16-bit pairs.
- `Y(kb_poly_table)`: nearest cell below, nearest above, exact hit.
- `Y(kb_poly_init)`: on an exact `(sigma, h)` hit copy the cell, compute the
  hut series, run the hold-out; otherwise solve from scratch.
- Harness `kbfrac.c`, `kbinit.c`: Appendix A rows and plan-time cost.

## Result

- 104 cells for `m <= 11` (`--mmax 11 -j 3`, 10527 s). The 113-bit cells
  `m >= 12` are issue 14.
- Exact hit 0.34 to 0.40 ms per axis in double. From scratch: median
  74.7 ms, max 277 ms (issue 16).
- `kbpoly_table` asserts an exact hit at each format's own degree for every
  `(sigma, m <= KB_POLY_OPT_MMAX)`: double 72, float 24, 113-bit 88 cells.
