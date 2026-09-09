# Spec: adopt FFTW's threading model and impatience lattice in `plan_ng`

**Status:** approved for planning, 2026-09-09
**Applies to:** the next-generation planner (`plan_ng`) only. Legacy
`nfft_plan` is untouched.

## Why

The planner already carries an FFTW-shaped impatience lattice
(`flags_t {l, u}`, `LEQ`, subsumption) and a thread count (`planner.nthr`,
hashed into the wisdom key). What it does not carry is FFTW's *policy*: the
public patience levels above `MEASURE`, and the rule that decides whether a
serial plan may compete against a threaded one.

Threaded solvers are coming (the convolution first). Transplanting the policy
now means the first threaded solver lands into machinery that already knows
how to choose against it, rather than forcing a redesign at that moment.

## The FFTW behaviour to mirror

Verified against the FFTW 3.3.11 source and manual.

1. **Four patience levels.** `ESTIMATE`, `MEASURE`, `PATIENT`, `EXHAUSTIVE`.
   `EXHAUSTIVE` implies `PATIENT`. `ESTIMATE` implies not `PATIENT`.

2. **Patience is expressed by the absence of `NO_*` bits.** There is no
   internal `PATIENT` bit. `api/mapflags.c:108-115` sets a block of `NO_*`
   flags whenever the request is *below* `PATIENT`:

   ```c
   IMPLIES(NO(FFTW_PATIENT),
           YES(FFTW_NO_VRECURSE | FFTW_NO_RANK_SPLITS | FFTW_NO_VRANK_SPLITS
               | FFTW_NO_NONTHREADED | FFTW_NO_DFT_R2HC
               | FFTW_NO_FIXED_RADIX_LARGE_N | FFTW_BELIEVE_PCOST))
   ```

   and `NO_UGLY` whenever below `EXHAUSTIVE`.

3. **The serial/threaded decision.** `kernel/ifftw.h:692`:

   ```c
   #define NO_NONTHREADEDP(plnr) \
     ((PLNR_L(plnr) & NO_NONTHREADED) && (plnr)->nthr > 1)
   ```

   When true the **serial** solvers decline the problem — `dft/ct.c:132` and
   `dft/vrank-geq1.c:139` return a null plan with the comment
   *"prefer threaded version"*. The serial plan is never built, so it is never
   timed. Only at `PATIENT` is the bit clear and both compete.

4. **The thread count is a parameter, never searched.**
   `fftw_plan_with_nthreads(n)` sets `plnr->nthr = imax(1, n)`
   (`threads/api.c:70-81`); the default is 1. It is a maximum: a threaded
   solver may use fewer and passes the remaining budget to its children.
   `fftw_planner_nthreads()` reads it back. FFTW never races thread counts.

5. **`nthr` is dynamically scoped across `mkplan`.**
   `kernel/planner.c:476-489` saves `ego->nthr` and `ego->flags` before
   `s->adt->mkplan(...)` and restores them after, so a solver may lower the
   budget for its children without leaking the change.

6. **`nthr` is part of the wisdom key.** `kernel/planner.c:170-177` feeds
   `plnr->nthr` into the md5 before the problem's own hash, so each thread
   count keeps its own entries.

## Requirements

**R1. Public patience flags.** Add `NFFT_PATIENT` and `NFFT_EXHAUSTIVE`
alongside `NFFT_MEASURE` and `NFFT_ESTIMATE`. Existing flag values must not
move: `NFFT_MEASURE = 0`, `NFFT_ESTIMATE = 1<<0`, `NFFT_NO_DIRECT = 1<<1`,
`NFFT_NO_FAST_NATIVE = 1<<4` stay as they are. New flags take fresh bits above
the highest in use; retired bits are not reused.

**R2. Public threading override.** Add `NFFT_NO_NONTHREADED` as a public
beyond-guru flag, mirroring FFTW's, so a caller can force the threaded
decision independently of patience.

**R3. Flag mapping stage.** The public planning word is normalised through an
`IMPLIES`-style table before it becomes `PLNR_*` bits, in its own translation
unit so it can be unit-tested without constructing a plan. Mirrors FFTW's
`api/mapflags.c`.

**R4. Thread count API.** Add `X(plan_with_nthreads)(int)` and
`X(planner_nthreads)(void)`. The planner default stays 1. The guru must stop
deriving `nthr` from `Y(get_num_threads)()`; the count becomes explicit, as in
FFTW.

**R5. Serial solvers decline.** Every currently registered solver that
performs no internal parallelism declines when `NO_NONTHREADEDP` holds. This
covers the NFFT kind (`fast_native`, `ndft_1d`, `ndft_nd`) and the DECONV and
CONV kinds. The rank-0 base case is exempt: there is nothing to parallelise
and it must remain the terminal solver for a fully elided problem.

**R6. Dynamic scoping.** `nthr` and the flag bounds are saved and restored
around every `mkplan` call the search makes, so a future threaded solver can
lower its children's budget.

**R7. Wisdom safety.** Entries written before this change carry `l`/`u` words
whose bit meanings have changed. Imports of such entries must be rejected
rather than misread. The existing configuration signature is the vehicle.

**R8. Observable behaviour today.** With no threaded solvers registered:
- `nfft_plan_with_nthreads(1)` — everything behaves exactly as before.
- `nfft_plan_with_nthreads(n > 1)` with `MEASURE` or `ESTIMATE` — every solver
  declines and `nfft_plan_ng_guru` returns `NULL`.
- `nfft_plan_with_nthreads(n > 1)` with `NFFT_PATIENT` or
  `NFFT_EXHAUSTIVE` — the serial plan is found and used.

R8 is the proof the machinery is wired. It is a deliberate, documented, tested
state, not an accident.

## Out of scope

- Any threaded solver, OpenMP pragma or thread pool.
- Racing thread counts. FFTW does not; we may later, and the wisdom key
  already admits it.
- `NO_UGLY` / `NO_SLOW` / `ALLOW_PRUNING` / `BELIEVE_PCOST` consumers. These
  bits exist in `iplanner.h` and no solver reads them. The mapping sets them
  per the lattice; giving them consumers is separate work.
- Legacy `nfft_plan`, `nfft_set_num_threads`, `nfsft`, `fpt`.
- Passing a patience level down to the child FFTW plans. Today
  `kernel/nfft/nfft-nd.c:210` forwards the caller's `fftw_flags` verbatim and
  that stays true.

## Acceptance

`make check` passes. New CUnit cases in the `planner` suite of
`tests/checkall_ng` cover the mapping table, the thread-count API, the decline
rule, dynamic scoping, wisdom keying by `nthr`, and every bullet of R8.
