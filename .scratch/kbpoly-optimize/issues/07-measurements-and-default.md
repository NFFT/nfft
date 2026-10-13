# 07: Measure plan-time cost and accuracy, decide the default

Status: ready-for-human

Part of: PRD F2, Q1, success metrics. Depends on: issues 03, 05, 06.

## Task

- Measure plan-time cost of `Y(kb_poly_init)` with the flag off, with the
  least-squares step, and with the Lawson step, for `m = 2, 5, 8, 11` at
  `sigma = 2` in double and float, on x86-64 and aarch64 (the dev container
  and the CI runners), and record the table here. The targets are PRD F2.
- Compare the accuracy report (`NFFT_BENCH_OUT`, see
  `docs/agents/accuracy-tracking.md`) for the Kaiser-Bessel cells before and
  after, with the flag on, in every precision and runtime.
- Confirm that no CodSpeed benchmark name changed and that the transform
  instruction counts are unchanged (they must be; the run time does not see
  the fit).
- Decide PRD Q1: default on or off. Record the decision and the evidence in
  `PRD.md`, and if the default flips, implement it as a follow-up issue (the
  flag then becomes an opt-out like `ANALYTIC_WINDOW`).

## Acceptance

- Tables of cost and accuracy in this file, the Q1 decision recorded in
  `PRD.md`.
