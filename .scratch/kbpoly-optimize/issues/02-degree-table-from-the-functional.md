# 02: Derive the degree table from the error functional

Status: ready-for-agent

Part of: PRD R2, D4. Depends on: issue 01. Background: `analysis.md` sections
6.1 and 9.

## Task

Change the rule in `tests/kbpolygen/generate.py` from the Chebyshev-tail bound
to the measured polynomial part: the degree of `(sigma, m, format)` is the
smallest `d` such that

    poly_part(d) <= 0.05 * J_inf(exact truncated KB)   or   <= eps_format * peak-relative floor,

where `poly_part` and `J_inf` come from issue 01 and the rounding-floor branch
keeps its present form (`2^(1 - mant)`). Evaluate every candidate degree in
order; do not assume the polynomial part is monotone in `d` (section 6.2 shows
it need not be). Keep the `sigma` bisection of the thresholds, the
rounding-down of thresholds, the monotonicity-in-`sigma` check and the
`m_opt` row limit as they are.

Also compute and print, per `(sigma, m, format)`, the degree under the
"measured" criterion (transform `J_inf` within 5 percent of Kaiser-Bessel's),
as information only; it must not set the table (see PRD D4 for why).

Regenerate `kernel/util/kbpoly_data.h`. Update `tests/kbpolygen/README.md`
("The rule" section) and the comment block above `Y(kb_poly_degree)` in
`include/infft.h`, and `AGENTS.md` section 5 if the command changes.

## Acceptance

- Double-precision degrees at `sigma = 2`: `m = 2 -> 6`, `m = 3 -> 8`; at
  `sigma = 4, m = 2 -> 6`; at `sigma = 1.5, m = 4 -> 9` (section 9 of the
  analysis, strict reading).
- No degree in any format rises relative to the current table.
- `make check` (Autotools) and `ctest` (CMake) pass in float, double and
  long double, Kaiser-Bessel window, with and without OpenMP, with the
  transform tests' bounds unchanged.
- The generator's existing test `test_degrees_at_known_points` is updated to
  the new values and a test pins the four cases above.
- The generator's run time stays under five minutes on the dev container.
