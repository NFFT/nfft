# 04 — Validation on an extended N set

Status: done

## Problem

`compare.c` ran over N in {243, 250, 251, 255, 256, 512}. Every one of those
has `t = next_power_of_2(N)/N <= 1.06`, so the dyadic ladder's rung 0 is
illegal in all 288 cases and the tuner would return the legacy grid every
time. That data cannot decide anything about this work.

## What to do

Point `compare.c` at `X(tune_plan_dyadic)` and extend the N set to span `t`
across (1, 2]:

| N | `t` | why |
|---|---|---|
| 256, 512 | 1.00 | control: rung 0 illegal, must not regress |
| 250, 251 | 1.02 | control, near the current set |
| 200 | 1.28 | just above the 5/4 threshold |
| 160, 320 | 1.60 | mid band |
| 140 | 1.83 | near the top of the band |
| 1100 | 1.86 | the worked example in the request |
| 600 | 1.71 | second octave |

Keep the rest of the protocol: three shapes (`M = N/4`, `N`, `4N`), the same
goal ladder, both directions, all three precisions, everything through the
same `plan_ng` path.

## Report

`dcompare_report.py` renders `docs/tuning-dyadic.md`, split by shape and by
`t` band (`t < 5/4` against `t >= 5/4`). The split is the whole point: below
the threshold the ladder can only return the legacy grid.

## Acceptance

For `X(tune_plan_dyadic)` against legacy with the oracle cut-off:

1. zero accuracy misses,
2. no shape's median speedup below 1.00x,
3. overall median speedup above 1.00x.

The control rows (`t <= 1.06`) must sit at 1.00x within noise. A control row
below 1.00x means the cut-off at rung 1 is worse than the shipped model's,
which would be a regression in issue 02's constants rather than in the ladder.

## Note

Long double runs an `O(N*M)` direct NDFT per case and is software binary128
here. `N = 1100` with `M = 4N` in long double is the most expensive cell in
the matrix — time one before launching the full run.
