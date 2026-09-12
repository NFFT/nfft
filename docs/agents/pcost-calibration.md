# pcost calibration

How the planner's prune gate is calibrated against measurement, and how to
redo it. The sweep program is deliberately *not* in the tree: it is a one-off
calibration harness, not a test. Recreate it from this document when the
operation-count estimator changes, when the fast solver's inner loops change,
or when the ratio needs re-deriving on a different host.

Last run: 2026-09-12, dev container (Ubuntu 24.04, gcc-14), double precision,
`kaiserbessel`, serial, 976 geometries.

## What is being calibrated

Two solvers compete for the same NFFT problem: direct NDFT and planner-native
fast NFFT. Both now report an operation count (`opcnt {add, mul, fma, other;}`,
`include/iplanner.h`) filled in at `mkplan` time -- `other` is loads, stores
and index arithmetic, which no cost model counted before. `Y(iestimate_cost)`
(`kernel/planner/planner.c`) collapses that record into one number with
FFTW's own weights: `add + mul + (HAVE_FMA ? 1 : 2)*fma + other`. The fast
solver's estimate is the sum of its DECONV child's count, its CONV child's
count, and `FFTW(flops)` on the FFTW plan it actually built -- not a fitted
formula.

`pcost` holds that estimate in estimate mode, and a real measurement (or the
estimate as a fallback) otherwise. It gates the measured race through
`PLNR_PRUNE_RATIO` (`include/iplanner.h`): a candidate whose estimate exceeds
that multiple of the cheapest candidate's is never timed. That gate is only
as good as the two solvers' estimates being on a *common* scale.

Unlike FFTW, we cannot rule the direct solver out with a fixed size constant
(FFTW's `GENERIC_MAX_SLOW` and friends, `kernel/ifftw.h`). FFTW's slow path is
reachable only at prime `n` and its fast path is exact. Ours is applicable at
every geometry and the fast path's cost carries `m` and `sigma`, so the
crossover moves: the direct wins at `d=1 N=8192 M=16` and at `d=3 N=32 M=64`.
A model gate is required; a constant cannot express that.

Both solvers used to carry their own hand-fitted cost scalar (`Y(nfft_ndft_pcost)`
with a fitted `NDFT_COST_UNIT`, and a `2*Ntot + 5*ntot*log2(ntot) +
2*M*(2m+2)^d` formula for the fast solver) instead of counting operations.
Both, and their constants, are gone; this sweep is what replaced them and
what re-derives the gate against the replacement.

## The sweep

One C program, built against the installed public headers, ~180 lines.

For each geometry it builds two estimate-mode plans over the same random
`x`/`f_hat`:

```c
p_dir  = NFFT(plan_ng_guru)(d, Nv, NULL, nv, M, m, window, x, f_hat, f_d, 0u,
                            NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
p_fast = NFFT(plan_ng_guru)(d, Nv, NULL, nv, M, m, window, x, f_hat, f_f, 0u,
                            NFFT_ESTIMATE | NFFT_NO_DIRECT);
```

A NULL `p_fast` means the fast guards refused the geometry (`N > m`,
`n > 2m+2`, `n > N`); skip the row.

Then it `NFFT(precompute)`s both, reads each plan's estimate out of `pcost`,
and times one `NFFT(execute)` of each.

**Reading pcost.** There is no accessor. `NFFT(fprint_plan)` prints
`pcost=` and the parent prints its own figure before descending into
children, so the first occurrence is the top-level one -- an estimate-mode
plan's own `Y(iestimate_cost)` result. Print into an `open_memstream` buffer
and `strtod` past `strstr(buf, "pcost=") + 6`. Verified by hand before the
sweep ran: at `N=64, n=128, M=16, m=6` the direct plan prints
`pcost=16480` and the fast plan's top line prints `pcost=4576` (its
DECONV/CONV children print `pcost=0` -- only the composite parent's estimate
is written, which is exactly the figure this sweep reads).

**Timing.** Re-run `NFFT(execute)` until a wall budget elapses
(`NFFT(clock_gettime_seconds)()`, 0.2 s per step is enough), report the mean.

**Grid.**

| rank | N | M | m | sigma |
|---|---|---|---|---|
| 1 | 16 … 8192, powers of 2 | 16, 64, 256, 1024, 4096 | 2, 4, 8, 12 | 1.25, 1.5, 2.0 |
| 2 | 8 … 128 | same | same | same |
| 3 | 4 … 32 | same | same | same |

`n` is the smallest even integer strictly above `sigma*N`. Rows with
`Ntot*M > 5e7` are skipped to keep one direct apply affordable; rows the fast
guards refuse are also skipped. 1140 candidate geometries, 2 dropped by the
size cap, 140 by the fast guards, 976 rows survive.

Emit TSV: `d N M m sigma n pcost_dir pcost_fast r_model t_dir t_fast r_meas`.

## The analysis

The single number that matters is **seconds per model cost unit** for each
solver separately, `t/pcost`. If the estimator's shape is right, that number
is constant across the whole grid.

Results:

| solver | median s/unit | p05-p95 spread |
|---|---|---|
| direct | 1.05e-10 | 1.32x |
| fast | 7.34e-11 | 3.38x |

**The direct estimate is good.** 1.32x over three ranks and four decades of
size -- essentially unchanged from the hand-written model's 1.28x, which makes
sense: a direct NDFT's operation mix does not change shape when it is counted
instead of formula-fit.

**The fast estimate's spread shrank a lot: 3.38x, against 5.08x before.**
Splitting the grid by `t_fast` into its smallest and largest deciles gives a
median s/unit of 6.79e-11 at the small end and 7.80e-11 at the large end -- a
1.15x drift with problem size, against a 3x drift measured under the old
formula. `other` (loads, stores, index arithmetic) evidently accounts for
most of what the old fitted formula missed. What is left is presumably still
cache and memory bandwidth, which no static operation count can see; there is
no formula here to refit, since the estimate is an exact sum of FFTW's own
`opcnt`, not a fitted shape.

**The defect the original sweep found** -- the direct and fast models sitting
on different unit scales, which `NDFT_COST_UNIT` existed to correct -- is
much smaller now but has not vanished: the two medians differ by **1.43x**
(1.05e-10 / 7.34e-11), against 2.31x under the hand-written models. Sharing
FFTW's own weighting between both solvers, instead of fitting each
separately, removed most of the scale mismatch but not all of it -- the two
solvers' operation mixes (a direct NDFT's trig evaluations and complex
multiplies vs. the fast solver's FFTW butterfly, DECONV window evaluation and
CONV gather/scatter) still cost different amounts of wall-clock per counted
operation on this host.

## Judging the gate

Do not score the gate on "was the model ratio close to the measured ratio".
The two errors that cost anything are:

- **harmful prune** — the direct was actually faster (`r_meas < 1`) but its
  estimate exceeded the ratio, so it was never timed. The planner ships the
  wrong solver.
- **waste** — wall seconds spent timing a direct that then lost.

| ratio | harmful prunes | total waste | worst single |
|---|---|---|---|
| 1.5 | 3 | 87 ms | 32 ms |
| 1.6 (observed floor) | 0 | 96 ms | 32 ms |
| 2.0 | 0 | 228 ms | 32 ms |
| **4.0 (kept)** | **0** | **316 ms** | **32 ms** |
| 8.0 | 0 | 694 ms | 59 ms |

The largest `r_model` among rows where the direct solver actually won
(`r_meas < 1`) is 1.565; below that, some ratios prune a winner (1.5 already
does, with 3 harmful prunes). 4.0 keeps a 2.5x margin over that observed
floor -- more margin than the old sweep's own 4-over-a-2-floor choice -- while
the wasted timing it costs is trivial in absolute terms (316 ms summed across
the whole 976-geometry grid, 32 ms worst case). Nothing in this sweep
justifies moving off 4.0.

## What was applied

`PLNR_PRUNE_RATIO` stays at 4.0 in `include/iplanner.h`. This is a
confirmation, not a no-op: the estimator underneath it changed completely
(operation counts instead of two separately-fitted formulas), the unit-scale
gap it has to cover shrank from 2.31x to 1.43x, and the ratio was re-measured
against the new numbers rather than carried over.

Three `tests/nplan.c` cases need a geometry the gate admits, because they
assert on a two-candidate race actually happening: `destructive_default`,
`timelimit_tight_adopts_untimed`, and `awake_zero_internal`. Their fixture
(`N = 64, n = 128, M = 16, m = 6`) was not retuned: under the new estimator it
prints `pcost=16480` direct vs `pcost=4576` fast, `r_model = 3.60`, still
below the 4.0 ratio, so the direct candidate still survives the gate and the
race still happens. Any future change to the ratio or to either solver's
operation count has to re-check this; the probe is two `NFFT(plan_ng_guru)`
calls with `NFFT_NO_FAST_NATIVE` and `NFFT_NO_DIRECT` plus the `fprint_plan`
pcost read described above.

This check ran in double precision only, and the answer is
precision-dependent: the fast estimate carries `FFTW(flops)` on the plan the
precision's own FFTW library built, and on this host `fftw3f` plans `n = 128`
cheaper than `fftw3`/`fftw3l` (2504 against 4576), which puts the same
geometry at `r_model = 6.58` in float — above the gate, so the direct
candidate is pruned and those tests do not race in a float build. See
ADR-0006.

The ratio is a host calibration, not a law: it is seconds per counted
operation on this machine, and it moves with the machine. Re-run the sweep if
either solver's operation count changes, or if the ratio is ever suspected of
being wrong on a different host.
