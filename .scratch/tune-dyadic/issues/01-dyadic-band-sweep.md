# 01 — Measure the `(4, 8]` oversampling band

Status: done

## Problem

The dyadic ladder's third rung has `sigma = 4t` with `t` in (1, 2], so it
occupies `sigma` in (4, 8]. `gsweep.c` stops at `sigma = 4`, so that band has
no measurements and rung 2 has no constants.

## What to do

`gsweep.c` hardcodes its sigma list at line 79. Make it an optional fifth
argument, comma-separated, defaulting to the current list, and pass it through
`run_gsweep.sh`. Then measure

    sigma in {5, 6, 7, 8}

which completes the ratio-2 triples `(t, 2t, 4t)` for `t` in
{1.25, 1.5, 1.75, 2} against the sigma values the sweep already holds. Same N
list, M factors, cut-off cap and trial count as the shipped `gsweep-*.csv.gz`,
so the two merge.

`sigma = 4.25` and similar intermediate values are not needed: the fit only
has to cover the band, and four values across (4, 8] with the full N x M x m
grid behind each is the same density the other bands have.

## Done when

- `run_gsweep.sh <root> <out> 32 5 "" "" "5,6,7,8"` produces one CSV per
  precision that merges with the existing data.
- `data/gsweep2-{d,f,l}.csv.gz` committed.
- `dfit.py` reports a non-empty band `(4, 8]` with a non-singular fit.

## Note

Long double is software binary128 here and the existing 13-sigma sweep took
about 40 minutes. Four sigmas should be proportionally shorter, but `n = 8N`
at `N = 1024` is the largest grid the study has yet transformed — check the
run time before assuming it fits.
