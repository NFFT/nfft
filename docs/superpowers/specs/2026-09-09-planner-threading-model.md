# Spec: adopt FFTW's threading model and impatience lattice in `plan_ng`

**Status:** approved for planning, revised 2026-09-09 after review
**Applies to:** the next-generation planner (`plan_ng`) only. Legacy
`nfft_plan` is untouched.

## Why

The planner already carries an FFTW-shaped impatience lattice
(`flags_t {l, u}`, `LEQ`, subsumption) and a thread count (`planner.nthr`,
hashed into the wisdom key). What it does not carry is FFTW's *policy*: the
public patience levels above `MEASURE`, the rule that decides whether a serial
plan may compete against a threaded one, and the packaging that makes that rule
safe.

Threaded solvers are coming, the convolution first. Transplanting the policy
now means the first threaded solver lands into machinery that already knows how
to choose against it.

A second, independent problem surfaced during review: the flags we hand to the
child FFTW plans are inconsistent with the patience the caller asked for, and
FFTW's own thread count can change our child plans without changing our wisdom
key. Both are fixed here, because both become worse once threading exists.

## The FFTW behaviour to mirror

Verified against the FFTW 3.3.11 source and manual.

1. **Four patience levels.** `ESTIMATE`, `MEASURE`, `PATIENT`, `EXHAUSTIVE`.
   `EXHAUSTIVE` implies `PATIENT`; `ESTIMATE` implies not `PATIENT`.

2. **Patience is the absence of `NO_*` bits.** There is no internal `PATIENT`
   bit. `api/mapflags.c:108-115` sets a block of restrictions whenever the
   request is *below* `PATIENT`:

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

   When true the serial solvers decline — `dft/ct.c:132`,
   `dft/vrank-geq1.c:139`, comment *"prefer threaded version"*. The serial plan
   is never built, so it is never timed. `NO_NONTHREADED` is **not** in
   `search()`'s `relax_tab` (`kernel/planner.c:578-585`), so no fallback
   relaxation rescues it.

4. **Threading ships as a separate library that is linked in addition.**
   `--enable-openmp` builds `libfftw3_omp`, linked as `-lfftw3_omp -lfftw3`.
   `--enable-threads` builds `libfftw3_threads` (pthreads). A program links at
   most one of them. `--with-combined-threads` folds threading into the main
   library instead.

5. **Threaded code enters through two doors, both from the threads library.**
   `X(threads_conf_standard)` registers extra solvers (`threads/conf.c`), and
   `threads_register_hooks()` sets `X(mksolver_ct_hook)` and
   `X(mksolver_hc2hc_hook)` (`threads/api.c:26-30`) so the main library's
   existing radix registration sites produce threaded solvers.

6. **Asking for threads is what installs them.** `fftw_plan_with_nthreads`
   calls `X(cleanup)()` and then `X(init_threads)()` when threads are not yet
   initialised (`threads/api.c:70-81`), so FFTW can never reach a state where
   `nthr > 1` and no threaded solver is registered. `init_threads` returns zero
   on failure. The `cleanup()` is deliberate: the roster is about to change, so
   the planner and its wisdom are wiped.

7. **The thread count is a parameter, never searched.** It is a maximum; a
   threaded solver may use fewer and passes the remainder to its children.
   `fftw_planner_nthreads()` reads it back (FFTW 3.3.9 and later).

8. **`nthr` is dynamically scoped and is part of the key.**
   `kernel/planner.c:476-489` saves and restores `ego->nthr` and `ego->flags`
   around `mkplan`; `kernel/planner.c:170-177` feeds `plnr->nthr` into the md5
   before the problem's own hash.

## Requirements

### Patience levels

**R1. Public patience flags.** Add `NFFT_PATIENT` and `NFFT_EXHAUSTIVE`.
Existing values do not move: `NFFT_MEASURE = 0`, `NFFT_ESTIMATE = 1<<0`,
`NFFT_NO_DIRECT = 1<<1`, `NFFT_NO_FAST_NATIVE = 1<<4`. New flags take fresh
bits above the highest in use; retired bits `1<<2` and `1<<3` stay unused.

**R2. Public threading override.** Add `NFFT_NO_NONTHREADED`, mirroring FFTW's
beyond-guru flag, so a caller can force the decision independently of patience.

**R3. Flag mapping stage.** The public planning word is normalised through an
`IMPLIES`-style table before it becomes `PLNR_*` bits, in its own translation
unit so it is unit-testable without constructing a plan. Mirrors
`api/mapflags.c`.

**R4. The decline rule.** Add `PLNR_NO_NONTHREADED` and the macro
`NO_NONTHREADEDP(pl)` in FFTW's exact two-condition form. Every registered
solver that performs no internal parallelism declines when it holds: the NFFT
kind (`fast_native`, `ndft_1d`, `ndft_nd`), DECONV and CONV. The rank-0 base
case is exempt — nothing to parallelise, and it must stay the terminal solver
for a fully elided problem.

**R5. Dynamic scoping.** `nthr` and the flag bounds are saved and restored
around every `mkplan` the search performs, so a future threaded solver can
lower its children's budget without leaking the change to the next candidate.

### Packaging

**R6. A separate add-on library.** Threaded planner code lives in
`libnfft3<suffix>_ng_omp`, built only under `--enable-openmp`, linked **in
addition** to `libnfft3<suffix>`. It contains the thread-count API and the
threaded solver roster and nothing else.

The existing `libnfft3<suffix>_omp` is a whole-kernel OpenMP rebuild that
callers link *instead of* the serial library. It keeps that meaning and is not
touched. The two models coexist; the documentation must say which is which.

**R7. Asking for threads installs them.** `X(plan_with_nthreads)`,
`X(planner_nthreads)`, `X(init_threads)` and `X(cleanup_threads)` are defined
**only** in the add-on library, declared unconditionally in `nfft3.h`.
`X(plan_with_nthreads)` calls `X(init_threads)` when threads are not yet
initialised, and `X(init_threads)` forgets wisdom and destroys the planner
before registering the threaded roster, because the roster change invalidates
both.

This reproduces FFTW's invariant: a program linked against the main library
alone cannot raise `nthr` above 1, so `NO_NONTHREADEDP` can never fire without
threaded solvers present. No `threaded` marker on `solver_adt` and no per-kind
counter are needed.

**R8. An empty roster is a reported failure, not a trap.** Until the first
threaded solver exists the roster registers nothing. `X(init_threads)` must
then return 0, as FFTW's does on failure, and `X(plan_with_nthreads)` must
leave the count at 1. The day a threaded solver is added, the roster is
non-empty and the machinery switches on with no further change.

### Child FFTW plans

**R9. Derive the child flags from patience.** `problem_nfft.fftw_flags`
already documents `0 = derive` (`iplanner.h:411`) and nothing derives:
`nfft-nd.c:210` passes the word through, so `0` reaches FFTW as
`FFTW_MEASURE`, while the planner skill documents `FFTW_ESTIMATE`. Three
sources, three answers. Implement the documented contract:

| NFFT level | derived child flags |
|---|---|
| `NFFT_ESTIMATE` | `FFTW_ESTIMATE` |
| `NFFT_MEASURE` | `FFTW_MEASURE` |
| `NFFT_PATIENT` | `FFTW_PATIENT` |
| `NFFT_EXHAUSTIVE` | `FFTW_EXHAUSTIVE` |

always with `| FFTW_DESTROY_INPUT` and `& ~FFTW_PRESERVE_INPUT`, as now. A
non-zero caller word is honoured verbatim after the same normalisation, so the
escape hatch survives. This is backward compatible for the common case: `0`
under `NFFT_MEASURE` still yields `FFTW_MEASURE`.

Without this an `NFFT_PATIENT` plan hands `FFTW_MEASURE` to the FFT child,
which is the exact trap the patience lattice exists to avoid.

**R10. Observe FFTW's thread count in the key.** The problem hash records
`pl->nthr` and `ego->fftw_flags` but not FFTW's own global count, so a caller's
`fftw_plan_with_nthreads(8)` changes our child FFT plans without changing our
key. Feed FFTW's count into the NFFT problem hash. We never write FFTW's
global state.

`fftw_planner_nthreads` is declared in `fftw3.h` but **defined only in FFTW's
threads library**, and `libnfft3` links `@fftw3_LIBS@` alone. The main library
therefore cannot call it directly without breaking the serial link. Use a
function-pointer hook, the same device FFTW uses for `mksolver_ct_hook`
(`threads/api.c:26-30`): the main library owns a null-by-default hook and
hashes 1 when it is null; the add-on library sets it to
`FFTW(planner_nthreads)` in `X(init_threads)` and clears it in
`X(cleanup_threads)`.

Consequence to document: FFTW's count is observed only when
`libnfft3_ng_omp` is linked. Driving FFTW's thread count without it is
unsupported and the key will read 1.

`fftw_planner_nthreads` is FFTW 3.3.9 and later, so it needs a configure probe.
When it is absent the add-on library is not built and the hook stays null.

## Out of scope

- Any threaded solver, OpenMP pragma or thread pool. The roster ships empty.
- Racing thread counts. FFTW does not; the wisdom key already admits it.
- Consumers for `PLNR_NO_UGLY`, `PLNR_NO_SLOW`, `PLNR_ALLOW_PRUNING` and
  `PLNR_BELIEVE_PCOST`. They exist in `iplanner.h` and no solver reads them.
  The mapping sets them per the lattice; giving them consumers is separate work.
- Intercepting our own solver construction the way FFTW's `mksolver_ct_hook`
  does. We have no factory-registered radix machinery for it to intercept; the
  hook in R10 is a different device with the same shape.
- A pthreads variant. OpenMP only, as agreed.
- Legacy `nfft_plan`, `Y(set_num_threads)`, `libnfft3_omp`, `nfsft`, `fpt`.

## Acceptance

`make check` passes across the double, float, long-double and OpenMP configure
matrix. New CUnit cases cover the mapping table, the decline rule, dynamic
scoping, the derived child flags, the extended problem key, and the add-on
library's API. A new `checkall_ng_ngomp` binary, built only under
`--enable-openmp`, links `libnfft3 + libnfft3_ng_omp` and covers R7 and R8.

## Observable behaviour on delivery

- A program linked against `libnfft3` alone: unchanged in every respect except
  the child FFTW flags now following the patience level (R9) and the wisdom key
  now including FFTW's thread count (R10). Both invalidate old wisdom, which
  the configuration signature turns into a clean miss.
- A program also linked against `libnfft3_ng_omp`: `nfft_init_threads()`
  returns 0 and `nfft_plan_with_nthreads(n)` leaves the count at 1, because the
  threaded roster is empty. Every plan is serial, as before.
