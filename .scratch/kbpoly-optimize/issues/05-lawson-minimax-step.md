# 05: Lawson minimax step on top of the least-squares correction

Status: done in part (feature/kbpoly-optimize): three floored reweightings inside Y(kb_poly_correct), accepted only when the RMS stays within a tenth of least squares; a full minimax level remains open

Part of: PRD G3, R6. Depends on: issues 03 and 04. Background: `analysis.md`
sections 5.2, 6.2; prototypes `kbopt3.py`, `kbopt4.py`.

## Task

Extend the correction step of issue 03 with Lawson's iteration for the
worst-case objective:

- Grid `120 xi x 96 t` (the `40 x 32` grid lets the minimax solution exploit
  the gaps, section 6.2).
- Weights `w_i` start uniform; each iteration solves the weighted system
  (rows scaled by `sqrt(w_i)`, same QR as issue 03), computes `|r_i|` of the
  new iterate, keeps the iterate with the smallest maximum, and updates
  `w_i <- w_i |r_i| / sum w |r|`. At most 20 iterations; stop early when the
  maximum has not improved for 5 iterations.
- Guarded by `KB_POLY_BR_MM` (`K = 100`, issue 04).
- Selected by a second flag bit or a value of the flag of issue 03; decide
  with the maintainer and document.

Cost target: under 20 ms per axis at `m = 11, d = 12` in double; measure and
record. If the dense grid makes the QR too slow, fold the `t` dimension: the
rows for one `xi` share the Chebyshev factor, so the Gram contribution per
`xi` is a Kronecker product; this is the optimization of section 5.1 and may
be deferred to a follow-up if the target is met without it.

The generator certifies the step: for the four cases of section 6.2 it
computes the LP optimum on the design grid (`scipy.optimize.linprog`, HiGHS,
32-gon modulus) as a lower bound and checks that the C implementation's
fine-grid `J_inf` is within 25 percent of it.

## Acceptance

- Double, `sigma = 2, m = 2, d = 6`: fine-grid `J_inf` at most `8e-6`
  (section 6.2: `7.5e-6`); `sigma = 4, m = 2, d = 6`: at most `3e-7`.
- The RMS error of the Lawson result is within 2x of the least-squares result
  (section 6.2 saw 1.3x to 1.9x); if it is worse, the best-iterate rule is
  not working.
- Float with `E == R`: the guard rejects every case of the float32 table in
  section 10 where Lawson lost its gain, and the admitted cases match double
  to 10 percent.
