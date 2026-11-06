# 15: Rounding of the hut table at the band edge

Status: needs-triage

Part of: PRD addendum (2026-11-06). Related: issue 09.

## Problem

The hold-out and the transform read `psi_hat` from the 64-term Chebyshev hut
table. Its rounding is about `eps C(0)` absolute, so at the band edge the
relative error is `eps C(0) / C(edge)`. Measured in double:

| case | hold-out, hut table | hold-out, exact `psi_hat` | cell `e` |
|---|---|---|---|
| sigma 1.1, m 8 | 7.58e-9 | 1.26e-9 | 1.18e-9 |
| sigma 1.25, m 8 | 2.58e-11 | 1.12e-12 | 8.17e-13 |
| sigma 2, m 6 (kbfrac cell set) | 4.97x `e` | 1.39x `e` | |

The `kbpoly_table` cell check carries a slack `4 (2m + 2) eps_R / cmin` for
this. Consequence: the optimized window's gain at small sigma and large m in
double is capped by the hut table, not by the window.

## Task

Evaluate `psi_hat` at the band edge to a few ulp relative: for example a
series in a variable that resolves the edge, a relative rather than an
absolute truncation rule, or the series of `psi_hat / C_KB` against the
Kaiser-Bessel formula.

## Acceptance

Double, sigma 1.1 and 1.25 at m 8: hold-out with the table within 1.2x of the
exact `psi_hat`; the slack in `kbpoly_table` drops; no change of the
transform's instruction count.
