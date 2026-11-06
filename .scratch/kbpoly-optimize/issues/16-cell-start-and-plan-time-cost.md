# 16: Cell start and plan-time cost of the correction

Status: needs-triage

Part of: PRD addendum (2026-11-06). Depends on: issue 12.

## Problem

`Y(kb_poly_init)` uses the table on exact hits only. The nearest-below cell
start in `Y(kb_poly_correct)` stays for a later quick mode: it was measured
slower than from scratch at equal quality (0.87x to 2.6x of the from-scratch
time, against the plan's target of 0.5x). Causes:

- The 400 cell points join the 216-point base design; 616 points triple the
  IPM cost per iteration.
- The lower bound `e_above / 1.05` is weak: the error falls steeply with
  sigma, so `lo` sits 4x to 5x below `hi` and round 0 runs a full bisection.
- The start saves only 1 to 2 exchange rounds.

Quality off the table: sigma 1.25, m 8 reached only 1.63x of the cell's `e`
with the exact `psi_hat` (8.17e-13 against 1.33e-12).

From scratch in double (host 10 cores, median of 7): `kbinit` OPT_WINDOW
minus plain init, median 73.8 ms (29 to 188 ms); `kbfrac` cell set median
74.7 ms, max 277 ms, against the plan's 60 ms median and 500 ms max. The IPM
takes 75 to 80 percent of the time; round 0 spends 9 to 11 solves (2
Dinkelbach steps plus a bisection from `hi / 4`).

## Options

- Round 0 by Dinkelbach steps only.
- Drop the base grid under a cell start.
- A subset of the cell points.

## Acceptance

From scratch median at most 60 ms on the `kbinit` set; a cell start at most
0.5x of from scratch at equal quality, then used by `Y(kb_poly_init)`.
