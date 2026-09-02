# 05 — Build `fg_exp_l` directly instead of by nested recurrence

Status: ready-for-agent

## Problem

`nfft_1d_init_fg_exp_l` (`kernel/nfft/nfft.c`) and the inline copies in the d-D
`MACRO_B` blocks and in `kernel/nfct/nfct.c` / `kernel/nfst/nfst.c` build
`exp(-l^2/b)` as

```c
  fg_exp_b0 = EXP(K(-1.0)/b);
  fg_exp_b0_sq = fg_exp_b0*fg_exp_b0;
  ...
  fg_exp_b2 = fg_exp_b1*fg_exp_b0;
  fg_exp_b1 *= fg_exp_b0_sq;
  fg_exp_l[l] = fg_exp_l[l-1]*fg_exp_b2;
```

Two chained recurrences: the relative error accumulates roughly `2l` ulp, about
50 ulp at the tail of an `m = 13` table.

## Change

`EXP(-(R)l * (R)l / b)`, or `EXP(-(R)(l*l) * b_inv)` once `b_inv` exists
(`03-hoist-constants.md`). The table is `2m+2` entries built once per transform,
not per node, so `2m+2` `EXP`s are not measurable against the `M`-sized loop that
follows, and every entry is then correct to 1 ulp.

If `02-centred-fast-gaussian.md` lands as the neighbour-ratio form, this table
disappears from the run path entirely and the change is moot for those call
sites — but `PRE_FG_PSI` still needs a per-axis table, so do it there.

## Tests

Covered by the run check in `02` / `04`.
