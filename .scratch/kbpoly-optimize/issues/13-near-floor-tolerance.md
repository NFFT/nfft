# 13: Tolerance of the driver near the rounding floor

Status: needs-triage

Part of: PRD addendum (2026-11-06). Background: `analysis.md` section 19.

## Problem

The residual `g_i c` is a cancelling sum of O(1) terms with an absolute
noise of about `10 eps` (measured). The driver decides levels at 1 percent,
which it can resolve only above about `1000 eps`. The transform itself
realizes errors down to `4 to 8 eps`. Between the two lies a band of two
decades where a better window still pays but the driver cannot tell windows
apart:

| format | band of the error | m at sigma 1.5 and up |
|---|---|---|
| float (`E = double`) | 6e-5 to 6e-7 | 2 to 3 |
| double | 2e-13 to 2e-15 | 8 to 9 |
| 64-bit | 2e-16 to 2e-18 | 10 to 11 |
| 113-bit | 2e-32 to 2e-34 | 20 to 22 |

In float `E = double` covers the band. In the `E = R` builds the band is
given up off the table and served on the table by the generator's optimum.

## Idea

Widen the relative tolerance with the level: `tol = max(1 percent,
30 eps / level)`. Open question: the feasibility margins of `Y(socp3_solve)`
(`DELTA / 10` on the objectives) and the floor guard (`500 eps`) must move
with it, or a widened tolerance stops the bisection on a noise decision.

## Acceptance

Off-table cells in the band of each `E = R` format reach at least the
from-scratch value of today, and none gets worse than plain.
