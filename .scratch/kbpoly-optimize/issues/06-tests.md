# 06: Tests for the corrected window

Status: ready-for-agent

Part of: PRD P2, C3. Depends on: issues 02, 03, 04 (05 for the Lawson cases).

## Task

`tests/kbpoly.c`:

- `check_kaiser_bessel_poly` (fit error against the Kaiser-Bessel function at
  the table degree) keeps testing the uncorrected interpolant only; make that
  explicit by building the table with the correction flag off, and relax its
  target to the new rule of issue 02 (the function error is no longer held to
  `0.05 pi e^(-bR)`; compute the expectation from the same functional the
  generator uses, or pin per-case values the generator prints).
- New `check_kaiser_bessel_poly_corrected`: for every `(sigma, m)` the guard
  admits and every precision, build both tables, evaluate `J_inf` of each by
  the section 3 formula in the test (a direct double sum over a modest grid is
  fine; it is a test), and assert `J_inf(corrected) <= J_inf(uncorrected)`.
  For `(sigma, m)` the guard rejects, assert the tables are identical.
- `check_kaiser_bessel_poly_trafo`: add the correction flag (and the Lawson
  flag once issue 05 lands) to the delegate matrix, so the transform tests run
  the corrected window against the direct transform under the existing
  bounds, in the NFFT, NFCT and NFST suites.

`tests/kbpolygen/tests/`: tests from issues 01, 02, 04 and 05 as specified
there.

Accuracy report: the Kaiser-Bessel cells will show a step down; no change to
the report tooling is needed, but the PR must link the before and after
report for the record (PRD success metrics).

## Acceptance

- `make check` and `ctest` pass in float, double, long double, with and
  without OpenMP, with the new cases registered in `tests/check.c` and in the
  CMake test list.
- The `_threads` binary passes too (the correction is per plan, not per
  thread, so no new parallel code; the check is that nothing regressed).
