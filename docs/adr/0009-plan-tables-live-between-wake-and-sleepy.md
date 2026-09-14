# ADR-0009: Plan tables live between wake and SLEEPY

## Status
Accepted (2026-09-14).

## Context
Every DECONV and CONV solver allocated its tables in `mkplan`. The race loops
build every candidate before timing (`Y(planner_candidates)`), and a candidate
returned to `SLEEPY` freed nothing, so a measured race held all candidates'
tables at once. CONV's `psi` is `M*d*(2m+2)` reals. Serial against threaded at
PATIENT already doubled it; a sorted CONV family (ADR-0008) would have made it
4x, 2.7 GiB at 3D, `m = 6`, M = 2^21. The guru's own race also kept its
running best at `AWAKE_ZERO`, so its peak was the best plus the current
candidate.

## Decision
Follow FFTW unless NFFT genuinely differs. FFTW, from the local mirror:
- tables are allocated when a plan leaves `SLEEPY` and freed on return
  (`dft/bluestein.c:130-146`, `kernel/twiddle.c:245`);
- `AWAKE_ZERO` still allocates, with cheap placeholder content
  (`kernel/trig.c:212`);
- every transition has `SLEEPY` on exactly one side (`kernel/plan.c:64`);
- a timed candidate goes back to `SLEEPY`, the best included
  (`kernel/timer.c:148-176`).

So:
- DECONV and CONV, serial and threaded, allocate `phi_hut_inv`, `psi`, `u`,
  the threaded band arrays and `perm` in `awake` on leaving `SLEEPY`, and free
  them on return. `AWAKE_ZERO` fills zeros where filling costs a window
  evaluation, and true values for indices.
- `Y(plan_awake)` ignores a request for the current state and asserts `SLEEPY`
  on exactly one side of every real change. The solvers' `level` fields are
  gone.
- The guru race puts every timed candidate back to `SLEEPY`, as
  `Y(planner_mkplan)` already did. `Y(precompute)` wakes the winner.
- `tests/alloc_count.h` counts live bytes through the malloc and free hooks.
  `check_nfast_conv_race_memory` and `check_ng_omp_race_memory` require a
  measured CONV race to peak below 1.5 times one candidate's tables.

## Deliberate divergences from FFTW
1. **The grids `g1`/`g2` stay allocated in `fast_native`'s `mkplan`.** FFTW's
   child plans are bound to them, and allocating hundreds of MiB on every
   execute would fault them in each time. There is one `fast_native` per race
   level, so they are never duplicated.
2. **Both race loops build every candidate before timing, and destroy
   nothing early.** The estimate prune (`PLNR_PRUNE_RATIO`) needs every
   estimate first, as do the lone-survivor adoption, the timelimit split and
   the guru's fallback. With tables in `awake` a sleeping candidate holds only
   small arrays. An FFTW-shaped search is recorded as its own work item
   (`.scratch/planner-fftw-search/PRD.md`).

## Consequences
- Peak table memory in a measured race is one candidate's tables.
- A race winner rebuilds its index tables at `precompute`: `u`, the bands, and
  for a sorted CONV plan `perm`. `psi` is filled at `AWAKE` either way.
- A direct `AWAKE_ZERO` to `AWAKE` move is a checked error under
  `--enable-debug`.
