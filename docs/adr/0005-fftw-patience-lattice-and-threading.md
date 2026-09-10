# ADR-0005: FFTW's patience lattice and threading model in `plan_ng`

## Status
Accepted (2026-09-10).

## Context
The next-generation planner (`plan_ng`) already carried an FFTW-shaped
impatience lattice (`flags_t {l, u}`, `LEQ`, subsumption) and a thread count
(`planner.nthr`, hashed into the wisdom key). What it did not carry was FFTW's
*policy*: a patience level above `MEASURE`, the rule that decides whether a
serial solver may compete against a threaded one, and the packaging that makes
that rule safe to ship. Three further defects sharpened the case for doing this
now, before any threaded solver lands: the flags handed to the internal FFTW
child plans ignored the NFFT patience level the caller asked for; FFTW's own
thread count could change which child plan gets built without changing our
wisdom key; and only the top-level solution was ever blessed, so an exported
wisdom file described a fraction of the plan tree — the DECONV/CONV children
that made the winner fast were missing from it.

## Decision
Mirror FFTW rather than invent a new scheme.

- **Patience is the absence of restriction.** `NFFT_PATIENT` widens the search
  by *lifting* restrictions (`PLNR_NO_NONTHREADED`, `PLNR_BELIEVE_PCOST`) rather
  than adding a "try harder" bit of its own. A mapping stage
  (`kernel/nfft/mapflags.c`) turns the public planning word into the internal
  `PLNR_*` image and into the derived child FFTW flags; nothing downstream
  reads the public word directly.
- **Serial solvers decline below `PATIENT`.** `NO_NONTHREADEDP(pl)` —
  `(PLNR_L(pl) & PLNR_NO_NONTHREADED) && pl->nthr > 1` — is checked by every
  solver that does not parallelise internally. Below `NFFT_PATIENT` a serial
  solver simply does not offer a plan when more than one thread was requested,
  so the threaded plan is chosen without ever being timed against the serial
  one. At `NFFT_PATIENT` and above the restriction lifts and the race decides.
- **The thread count is an explicit parameter, never searched.** `pl->nthr` is
  set by `X(plan_with_nthreads)`, saved and restored around every `mkplan` call
  the search makes (`invoke_solver`), and folded into the wisdom key. Nothing
  in the planner tries several thread counts and picks the fastest.
- **The thread-count API and the threaded roster ship in an add-on library**,
  linked *in addition* to `libnfft3<suffix>`, so asking for threads is what
  installs them — a program that never calls `X(plan_with_nthreads)` cannot
  raise the count above 1, and a program that never links the add-on cannot
  call it at all.
- **Blessing reaches the winning tree by a second pass**, not by a
  child-enumeration hook. After a measured race picks a winner (or, in
  estimate mode, after the single applicable candidate is chosen),
  `plan.c`'s `bless_whole_tree` re-runs `Y(planner_mkplan)` on the same
  top-level problem with `PLNR_BLESSING` set. Every node the first pass
  memoised is hit again by this second pass, so the whole tree the winner is
  built from gets blessed with no need to know what that tree contains.
- **Wisdom-only is a planner state outside the flags word**
  (`wisdom_state_t`, `struct planner_s::wisdom_state`), not a `PLNR_*` bit —
  it must never participate in `LEQ` subsumption or key a wisdom entry. It is
  set at the top of the guru from `NFFT_WISDOM_ONLY`, reset at every return
  path, so a refusal never survives past the call that produced it.

### Why the packaging matters
This is the whole safety argument. `NO_NONTHREADEDP` checks two conditions,
not FFTW's on-paper three, and carries no marker on `solver_adt` for "this
solver only kicks in with threads." That is safe only because a program that
has not linked the add-on library cannot raise `pl->nthr` above 1 in the first
place — the third condition FFTW would need is supplied by the link line
instead of by a runtime check. `X(init_threads)` returning 0 on an empty
threaded roster closes the one remaining window: until a threaded solver is
registered, `X(plan_with_nthreads)` cannot leave the count above 1 either, so
there is never a state with `nthr > 1` and nothing to serve it.

### The five deliberate divergences from FFTW
1. `NO_NONTHREADEDP` has two conditions, not three, because packaging supplies
   the invariant FFTW gets at runtime from `fftw_plan_with_nthreads`.
2. No `relax_tab`. FFTW's `search()` widens `PLNR_NO_UGLY`/`PLNR_NO_SLOW`
   across a stalled search; no NFFT solver reads either bit, so porting the
   table would be dead code. Wisdom still partitions across patience levels,
   exactly as in FFTW.
3. No patience-escalation loop under a timelimit
   (FFTW's `api/apiplan.c:108-125`). This is the main reason our patience
   levels do less work than FFTW's: a stalled measured race degrades straight
   to estimate grade rather than retrying at a lower patience first.
4. FFTW's own thread count is observed only when the add-on library is
   linked — through a null-by-default hook the add-on installs, because
   `fftw_planner_nthreads` lives in FFTW's threads library and `libnfft3`
   links `@fftw3_LIBS@` only.
5. `NFFT_EXHAUSTIVE` is reserved (`1U << 6`) rather than exposed, so its
   mapping branch in `levels()` ships dead and uncovered until the flag is
   made public.

### Naming
`libnfft3<suffix>_ng_omp` is the FFTW-shaped add-on introduced here: threading
support for `plan_ng`, opt-in by linking, empty roster until a threaded solver
exists. `libnfft3<suffix>_omp` remains what it always was: a whole-kernel
OpenMP rebuild of the legacy API, linked *instead of* the serial library, and
untouched by this work. Two libraries, two models, one letter apart — a
deliberate but unavoidable naming collision.

### The FFTW version floor
`fftw_planner_nthreads` — the function the add-on's hook calls — arrived in
FFTW 3.3.9. `configure.ac` probes for it with `AC_CHECK_DECL`; the add-on
library (`ENABLE_NG_OMP`) is built only when the declaration is found (and
`--enable-openmp` and FFTW's own threads library are both present). An older
FFTW leaves `plan_ng` fully functional, just without the add-on: no
`X(plan_with_nthreads)`, no threaded roster, always one thread.

## Divergences from the plan's text found during implementation
Four points where the implementation deviates from the original plan document,
recorded here because they change externally observable behaviour.

- **`pl->nthr` no longer follows `Y(get_num_threads)()`.** The line
  `pl->nthr = (int)Y(get_num_threads)();` was deleted from `Y(plan_ng_guru)`.
  The thread count defaults to 1 and is written only by the add-on's
  `X(plan_with_nthreads)`, so on an OpenMP build the wisdom key no longer
  tracks the OpenMP thread count — a caller who wants more than one thread
  must ask for it explicitly.
- **`FFTW_WISDOM_ONLY` is kept out of the key by `hash()`, not stripped at
  construction.** `keyable_fftw_flags` (`kernel/nfft/plan.c`) strips the
  preservation bits only; the wisdom-only bit stays in the word stored on the
  problem so the *child* FFTW plans built in `nfft-nd.c` are planned
  wisdom-only too. `hash()` in `kernel/nfft/problem.c` strips it from the key
  instead, so it shapes what gets planned but never what gets remembered.
  Stripping it at construction, as the plan's text called for, would have kept
  it out of the key but left the child FFTW plans un-gated.
- **The estimate path blesses more than the measured path's invariant
  promises.** `plan.c` sets `PLNR_BLESSING` around a full estimate search, not
  just around the chosen winner, so in estimate mode a *losing* candidate's
  children would also be memoised blessed if a losing candidate had any. "A
  losing candidate's children are not blessed" therefore holds in measured
  mode only — vacuous today because `nfft_solver_fast_native` is the only
  NFFT-kind solver with children (DECONV/CONV), but not actually enforced. See
  Known gaps.
- **The `fftw_planner_nthreads` configure probe is CPPFLAGS-safe.** The
  `AC_CHECK_DECL` for it is wrapped to save `CPPFLAGS`, prepend
  `$nfft_fftw3_CPPFLAGS`, and restore afterward — the same treatment
  `NFFT_LIB_FFTW3` gives its own probes. Without it, a build configured with
  `--with-fftw3-includedir` would probe against the wrong (or no) `fftw3.h`,
  silently conclude the declaration is absent, and drop the whole add-on
  library rather than fail loudly.

## Consequences
- **Every stored wisdom file is invalidated, for every user.** The
  configuration signature carries a vocabulary tag (`"plnr-flags-v2"` in
  `config_signature`, `kernel/planner/planner.c`), bumped because the
  `PLNR_*` bit vocabulary and the problem hash both changed. A mismatch is a
  clean import rejection, never a wrong plan.
- **Wisdom files grow roughly threefold.** Blessing now reaches DECONV and
  CONV children alongside the top-level NFFT solution, so a single measured
  win writes about three entries where it used to write one.
- **Estimate-only users begin writing wisdom.** Estimate solutions are now
  blessed and exported (FFTW parity); a program that never measures still
  accumulates a wisdom file, whose entries can only ever answer other estimate
  queries.
- **On OpenMP builds the key changes again, independently of everything else
  here**, because `nthr` no longer follows the OpenMP thread count — a program
  that only ever ran single-threaded before this branch sees the same key it
  always did; one that relied on the implicit OpenMP thread count no longer
  gets it for free.

## Known gaps
- `NFFT_EXHAUSTIVE`'s branch in `levels()` (`kernel/nfft/mapflags.c`) ships
  dead and uncovered: the reserved bit that would trigger it is not exposed.
- `PLNR_NO_UGLY`, `PLNR_NO_SLOW`, `PLNR_ALLOW_PRUNING` and
  `PLNR_BELIEVE_PCOST` are set correctly by the mapping stage and read by
  nothing — there is no relaxation search and no cost-model consumer for them
  yet.
- The child FFTW plans' out-of-place and destroy-input invariants at
  `kernel/nfft/nfft-nd.c:209-217` have no test; they are not independently
  observable from outside the plan.
- The guru's `NULL` return conflates "no wisdom" (a wisdom-only miss) with
  "bad arguments" (a release-safe guard). There is no error channel that
  distinguishes them.
- `X(plan_with_nthreads)` must be called before any other NFFT routine — it
  destroys the planner and its in-memory wisdom — and nothing enforces that
  ordering beyond documentation, mirroring FFTW's own contract for
  `fftw_plan_with_nthreads`.
- The estimate path blesses a losing candidate's children whenever a losing
  candidate has any: "losing candidates are not blessed" holds for the
  measured race only, and is unenforced (currently vacuous, since no losing
  NFFT-kind candidate has children).

## Not decided here
- Racing several thread counts against each other rather than taking the
  count as given.
- A pthreads (non-OpenMP) threading variant.
- FFTW's patience-escalation loop under a timelimit
  (`api/apiplan.c:108-125`) — a stalled measured race degrades straight to
  estimate grade instead of retrying at a lower patience first.
- An error channel that separates "no wisdom" from "bad arguments" in the
  guru's `NULL` return.
