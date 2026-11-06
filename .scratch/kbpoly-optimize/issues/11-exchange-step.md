# 11: Exchange step on a fine grid

Status: done (feature/kbpoly-optimize, 150956e9)

Part of: PRD addendum (2026-11-06). Depends on: issue 10. Background:
`analysis.md` section 19.

## Task

After the driver, add the 40 largest local maxima of `|E|` on a `401 x 65`
fine grid (each refined on an `11 x 11` patch) to the design and run the
driver again, `lo` carried over, until the fine maximum is within 2 percent
of `hi`, after 10 rounds, or after two rounds below `500 eps_E`. Keep the
round with the smallest fine maximum.

## Result

`kbpoly_driver_reference` 2.1821e-8 against 2.183e-8 (1.000x, was 1.086x
without the exchange). `kbpoly_exchange`: sigma 1.25, m 4, fine `J_inf`
1.014x of the returned maximum. Double minimax `opt` rows within 1.02x of
the Lawson rows or better, one floor row at 1.033x (150 eps).
