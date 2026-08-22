# Planning modes, the race & flags

How the planner chooses an algorithm. Ground truth: `kernel/nfft/plan_ng.c`
(the NFFT guru), `kernel/planner/planner.c` (search + lattice),
`kernel/planner/timer.c` (measurement), `include/iplanner.h`.

## Two modes

**Measured (default, `NFFT_MEASURE == 0`).** The guru races the applicable
candidate plans **on your actual nodes** at plan time and blesses the winner
into wisdom. The realistic node distribution drives the true memory-access
pattern, so the measurement reflects real performance. This is FFTW's
default-measure convention.

**Estimate (`NFFT_ESTIMATE`, `1<<0`).** No race. Each solver's `mkplan`
returns a plan with an analytic `pcost` (fixed default constants); the planner
keeps the cheapest applicable one and memoises it **unblessed**. Instant, and
**never** time-bounded. Use when you cannot afford planning time or want
deterministic selection.

There is **no separate `nfft_optimize()` verb** — measurement is always at plan
time. Nodes arrive at the guru; you do not "optimize later".

## The measured race (what actually happens)

Per direction (forward only today), under measured bounds `{l=u=F}` where `F` is
the `PLNR_*` image of your gate flags:

1. **Wisdom lookup** (`planner_hlookup`). Infeasible hit → this direction is
   absent. Feasible hit → re-run that solver's `mkplan` and adopt the plan
   without racing (stale hit → fall through to the race, which heals it).
2. **Enumerate candidates** (`planner_candidates`, cap 8) — every applicable
   solver of the kind under the current bounds. Zero → direction absent. One →
   adopt untimed and **bless** immediately.
3. **Two or more → race.** For each candidate: `plan_awake(AWAKE_ZERO)` (builds
   ψ, but may use bogus values), `plan_measure_cost(...)`, `plan_awake(SLEEPY)` 
   (invalidate ψ, restore for the next candidate). **Strictly smallest cost wins**; ties keep the
   earlier-encountered candidate. Losers destroyed; winner **blessed**, left
   `SLEEPY` (precompute awakens it).
4. **No usable clock** (`plan_measure_cost` returns exactly `-1.0`) → tear down
   the raced candidates and **restart in estimate mode** (unblessed). A
   timelimit-induced or clockless loser is never blessed.

The race is **value-blind**: it does *not* zero or otherwise touch `f_hat`/`f` —
it times `apply()` on whatever they hold. This is safe because NFFT trip counts
and access patterns derive entirely from the nodes `x`, not from the values. ψ
is *precomputed, not zeroed* (zeroing it would destroy the access pattern the
race exists to measure). Numerical output during the race is meaningless by
design.

### The timer (`kernel/planner/timer.c`)

`Y(plan_measure_cost)(plan, problem)` returns a cost in **arbitrary units**
(comparable only within one process) or exactly `-1.0` (no clock). Two clocks:

- **fine clock** — a high-resolution CPU tick counter (`getticks`/`elapsed` via
  `include/cycle.h`), used for the actual measurement. `cycle.h` is the *one*
  FFTW-derived (MIT-licensed) file the timing path depends on; `timer.c` itself
  is clean-room.
- **budget clock** — wall seconds (`clock_gettime`), bounds only the total
  measurement time.

Loop constants (`iplanner.h`): `PLNR_TIME_REPEAT` = 8 best-of batches per
doubling level; `PLNR_TIME_MIN_TICKS` = 100 (tick-mode acceptance floor);
`PLNR_TIME_MIN_SLOW_SECONDS` = 1e-3 (slow-timer floor);
`PLNR_TIME_LIMIT_SECONDS` = 2.0 (per-candidate wall budget). Where only
`clock_gettime` exists it falls back to timing in wall seconds (FFTW's
slow-timer path); where no clock exists it returns `-1.0` → estimate-grade.
The timer never touches awake state, never allocates, and is not thread-safe.

## The impatience / subsumption lattice

Wisdom entries and queries carry a pair of bitsets `(l, u)` over the internal
`PLNR_*` flags (`iplanner.h`). More bits set = **more impatient** = fewer
algorithms considered. The lattice order is subset inclusion:

```c
#define LEQ(x, y) (((x) & (y)) == (x))   /* x <= y  iff  x is a subset of y */
```

A solution `(l, u)` asserts: every component of the winning plan is at least as
impatient as `l`, and the search tried every solver at least as impatient as `u`
(smaller `u` = searched harder). Solutions maintain `LEQ(l, u)`.

**Subsumption law** — a stored solution `a` answers a query `q` on the same key:
- `a` feasible: `LEQ(a.u, q.u) && LEQ(q.l, a.l)`.
- `a` infeasible: `LEQ(a.l, q.l) && a.timelimit_imp <= q.timelimit_imp`.

Lookup returns the subsuming entry with the smallest `u`; insert kills every
live same-key entry the newcomer subsumes.

### Evidence grades (why measured wisdom answers estimate queries, not vice versa)

`PLNR_ESTIMATE` rides in the **upper bound only**:

- estimate query/memo: `{l = F, u = PLNR_ESTIMATE | F}`
- measured query/bless: `{l = u = F}`

So a measured (blessed) entry subsumes a same-`F` estimate query (you get the
better decision for free), while a *feasible* estimate memo never answers a
measured query (no poisoning of measured planning with estimate-grade
evidence). `PLNR_ESTIMATE` must **never** appear in `l`. Estimate-grade results
are **never blessed** (`planner_bless` asserts no `PLNR_ESTIMATE` in `u`).

## The public planning flags → internal images

Public `NFFT_*` (in `nfft3.h`) map to internal `PLNR_*` (in `iplanner.h`) via
`map_planning_flags` in `plan_ng.c`:

| Public flag | Value | Internal | Effect |
|-------------|-------|----------|--------|
| `NFFT_MEASURE` | `0` | — | Default. Race + bless. |
| `NFFT_ESTIMATE` | `1<<0` | `PLNR_ESTIMATE` (in `u`) | Analytic pick, unblessed, no race. |
| `NFFT_NO_DIRECT` | `1<<1` | `PLNR_NO_DIRECT` | Forbid O(N·M) direct/NDFT solvers. |
| `NFFT_NO_FAST_NATIVE` | `1<<4` | `PLNR_NO_FAST_NATIVE` | Forbid the native fast NFFT. |

A solver reads the current bounds (`PLNR_L(pl)`) inside its `mkplan` and returns
`NULL` when a gate forbids it — impatience gating and plain inapplicability are
the same mechanism.

`fftw_flags` is a **separate** parameter controlling the internal FFTW child
plans (`0` is fine; the child defaults to `FFTW_ESTIMATE`). Its planning-rigor
bits (`FFTW_ESTIMATE`/`MEASURE`/`PATIENT`/...) participate in the wisdom key;
the input-preservation bits (`FFTW_DESTROY_INPUT`/`FFTW_PRESERVE_INPUT`) are
stripped from the key (`keyable_fftw_flags`) so a caller's spelling never
fragments the cache — every native candidate is out-of-place anyway.

> There is **no patience ladder**. `NFFT_PATIENT` / `NFFT_EXHAUSTIVE` were
> dropped; do not expect them. See [history-and-drift.md](history-and-drift.md).

## Timelimit

`nfft_set_timelimit(double seconds)` sets a **process-global wall-clock budget**
on the *measured* race (`< 0` = unlimited, the default). It is declared under
`nfft_*` only but the global store is shared library-wide.

- Before measuring each candidate the guru checks the elapsed coarse wall clock;
  on expiry it stops measuring further candidates.
- If **no** candidate finished before expiry (or there is no usable clock),
  selection degrades to **estimate grade (unblessed)** — a timelimit-induced
  loser is never blessed. A candidate that *did* finish before expiry is a
  legitimate blessed winner.
- **Estimate mode is never time-bounded.**
- The timelimit is a runtime `double` on `struct planner_s`; it is **not part of
  the wisdom key** (the reserved `timelimit_imp` field stays 0).
