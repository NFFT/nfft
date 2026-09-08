# pcost calibration

How the planner's two analytic cost models were calibrated against measurement,
and how to redo it. The sweep program is deliberately *not* in the tree: it is a
one-off calibration harness, not a test. Recreate it from this document when a
cost model changes, when the fast solver's inner loops change, or when the
constants below need re-deriving on a different host.

Last run: 2026-09-08, dev container (Ubuntu 24.04, gcc-14), double precision,
`kaiserbessel`, serial, 961 geometries.

## What is being calibrated

Two solvers compete for the same NFFT problem:

- direct NDFT, `Y(nfft_ndft_pcost)` in `kernel/nfft/ndft-1d.c` — `Ntot * M`
  times a constant.
- planner-native fast NFFT, `pcost()` in `kernel/nfft/nfft-nd.c` —
  `2*Ntot + 5*ntot*log2(ntot) + 2*M*(2m+2)^d`.

`pcost` has no unit. It only has to rank candidates, and it gates the measured
race through `PLNR_PRUNE_RATIO` (`include/iplanner.h`): a candidate whose pcost
exceeds that multiple of the cheapest candidate's is never timed. That gate is
only as good as the two models being on a *common* scale.

Unlike FFTW, we cannot rule the direct solver out with a fixed size constant
(FFTW's `GENERIC_MAX_SLOW` and friends, `kernel/ifftw.h`). FFTW's slow path is
reachable only at prime `n` and its fast path is exact. Ours is applicable at
every geometry and the fast path's cost carries `m` and `sigma`, so the
crossover moves: the direct wins at `d=1 N=8192 M=16` and at `d=3 N=32 M=64`.
A model gate is required; a constant cannot express that.

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

Then it `NFFT(precompute)`s both, reads each plan's analytic pcost, and times
one `NFFT(execute)` of each.

**Reading pcost.** There is no accessor. `NFFT(fprint_plan)` prints
`pcost=` and the parent prints its own figure before descending into children,
so the first occurrence is the top-level one. Print into an `open_memstream`
buffer and `strtod` past `strstr(buf, "pcost=") + 6`.

**Timing.** Re-run `NFFT(execute)` until a wall budget elapses
(`NFFT(clock_gettime_seconds)()`, 0.2 s per step is enough), report the mean.

**Grid.**

| rank | N | M | m | sigma |
|---|---|---|---|---|
| 1 | 16 … 8192, powers of 2 | 16, 64, 256, 1024, 4096 | 2, 4, 8, 12 | 1.25, 1.5, 2.0 |
| 2 | 8 … 128 | same | same | same |
| 3 | 4 … 32 | same | same | same |

`n` is the smallest even integer strictly above `sigma*N`. Rows with
`Ntot*M > 5e7` are skipped to keep one direct apply affordable.

Emit TSV: `d N M m sigma n pcost_dir pcost_fast r_model t_dir t_fast r_meas`.

## The analysis

The single number that matters is **seconds per model cost unit** for each
solver separately, `t/pcost`. If a model's shape is right, that number is
constant across the whole grid.

Results:

| solver | median s/unit | p05-p95 spread |
|---|---|---|
| direct | 2.22e-10 | 1.28x |
| fast | 9.63e-11 | 5.08x |

**The direct model is good.** 1.28x over three ranks and four decades of size.

**The fast model has a 5x spread that is not a term-weighting error.** A
least-squares refit of `a*Ntot + b*ntot*log2(ntot) + c*M*(2m+2)^d` moves the
spread only to 4.81x and wants a negative `a`. The residual is a monotone drift
with problem size instead — 7.4e-11 s/unit at `t_fast ~ 4e-7 s`, rising to
2.2e-10 at `t_fast ~ 2e-3 s`. That is cache and memory bandwidth, which an op
count cannot see. Do not try to fix it analytically; finding it is what the
measured race is for.

**The defect the sweep found** is that the two medians differ by 2.31x, i.e.
the models were on different unit scales, which made the gate systematically
2.3x too permissive toward the direct solver.

## Judging the gate

Do not score the gate on "was the model ratio close to the measured ratio". The
two errors that cost anything are:

- **harmful prune** — the direct was actually faster (`r_meas < 1`) but its
  pcost exceeded the ratio, so it was never timed. The planner ships the wrong
  solver.
- **waste** — wall seconds spent timing a direct that then lost.

| | harmful prunes | total waste | worst single |
|---|---|---|---|
| before: direct unscaled, ratio 8 | 0 | 765 ms | 59 ms |
| direct x2.31, ratio 8 | 0 | 334 ms | 32 ms |
| **after: direct x2.31, ratio 4** | **0** | **81 ms** | **32 ms** |
| direct x2.31, ratio 2 | 4 | 17 ms | 8 ms |

The old gate was already safe; it was wasteful. Ratio 2 is where the fast
model's residual spread starts pruning winners, so 4 is the floor with margin.

## What was applied

- `NDFT_COST_UNIT 2.31` in `kernel/nfft/ndft-1d.c`, putting the direct model on
  the fast model's unit.
- `PLNR_PRUNE_RATIO` 8.0 -> 4.0 in `include/iplanner.h`.

Three `tests/nplan.c` cases need a geometry the gate admits, because they assert
on a two-candidate race actually happening: `destructive_default`,
`timelimit_tight_degrades_to_estimate`, `awake_zero_internal`. Their fixture
moved from `M = 64` to `M = 16` at `N = 64, n = 128, m = 6` (pcost ratio 11.2 ->
3.5). Any future change to the ratio or to either model has to retune them; the
probe is two `NFFT(plan_ng_guru)` calls with `NFFT_NO_FAST_NATIVE` and
`NFFT_NO_DIRECT` plus the `fprint_plan` pcost read described above.

Both numbers are host calibrations, not laws. 2.31 is the ratio of two
op-count-to-seconds conversions and moves with the machine; the sweep is how to
re-derive it. Re-run it if either cost model is edited.
