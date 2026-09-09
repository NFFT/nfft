# Spec: adopt FFTW's threading model, patience lattice and wisdom discipline in `plan_ng`

**Status:** settled after adversarial review, 2026-09-09
**Applies to:** the next-generation planner (`plan_ng`) only. Legacy
`nfft_plan` is untouched.

## Why

The planner already carries an FFTW-shaped impatience lattice
(`flags_t {l, u}`, `LEQ`, subsumption) and a thread count (`planner.nthr`,
hashed into the wisdom key). What it does not carry is FFTW's *policy*: the
patience level above `MEASURE`, the rule that decides whether a serial plan may
compete against a threaded one, and the packaging that makes that rule safe.

Threaded solvers are coming, the convolution first. Transplanting the policy now
means the first threaded solver lands into machinery that already knows how to
choose against it.

Review surfaced three further defects that become worse once threading exists,
so they are fixed here: the flags handed to the child FFTW plans ignore the
patience the caller asked for, FFTW's own thread count can change our child
plans without changing our wisdom key, and only the top-level solution is ever
blessed, so an exported wisdom file describes a fraction of the plan tree.

## The FFTW behaviour to mirror

Verified against the FFTW 3.3.11 source, not from memory.

1. **Patience is the absence of restriction.** No internal `PATIENT` bit.
   `api/mapflags.c:108-115` sets a block of `NO_*` flags whenever the request is
   below `PATIENT`, and `NO_UGLY` below `EXHAUSTIVE`.

2. **The serial/threaded decision.** `kernel/ifftw.h:692`:
   `NO_NONTHREADEDP(plnr)` is `(PLNR_L(plnr) & NO_NONTHREADED) && plnr->nthr > 1`.
   Serial solvers then decline (`dft/ct.c:132`, `dft/vrank-geq1.c:139`), so the
   serial plan is never built and never timed. `NO_NONTHREADED` is absent from
   `search()`'s `relax_tab` (`kernel/planner.c:578-585`), so nothing rescues it.

3. **Threading is a separate library linked in addition.** `--enable-openmp`
   builds `libfftw3_omp`, linked `-lfftw3_omp -lfftw3`. Its solvers register
   through `X(threads_conf_standard)` (`threads/conf.c`).

4. **Asking for threads installs them.** `fftw_plan_with_nthreads` calls
   `X(cleanup)()` then `X(init_threads)()` when threads are not yet initialised
   (`threads/api.c:70-81`), so FFTW cannot reach a state with `nthr > 1` and no
   threaded solver. `init_threads` returns zero on failure. The `cleanup()` is
   deliberate: the roster is about to change.

5. **The thread count is a parameter, never searched.** A maximum, dynamically
   scoped across `mkplan` (`kernel/planner.c:476-489`) and part of the wisdom
   key (`kernel/planner.c:170-177`).

6. **Wisdom does not span impatience levels.** `mkplan` performs exactly one
   `hlookup` with the caller's flags; `search()`'s relaxation drives the
   *solver enumeration*, never the lookup (`kernel/planner.c:518-615`). There is
   no `bless` that widens `l` — blessing is only the info bit choosing a table
   (`kernel/planner.c:392`). So FFTW partitions levels exactly as we do.

7. **Blessing reaches the whole winning tree by a second pass.**
   `api/apiplan.c:147` re-creates the finished plan with `BLESSING` set; every
   lookup on that pass hits the first pass's memo, and the bit is dynamically
   scoped down through `invoke_solver`. Losing candidates' children are not
   blessed. No child-enumeration machinery exists.

8. **Blessing is independent of patience.** `mkplan0` sets
   `plnr->flags.hash_info = BLESSING` whatever the level, and `exprt`
   (`kernel/planner.c`) writes every live blessed entry with no filter. Estimate
   plans are blessed and exported.

9. **Wisdom-only is a planner state, not an impatience bit.**
   `wisdom_state_t` (`kernel/ifftw.h:707-718`) sits outside the flags word. It
   is set per top-level planning call by `mkplan0`, so it does not persist
   across calls; within a call it is sticky and unwinds through children. A miss
   at `do_search:` jumps to `wisdom_is_bogus` and returns 0. Wisdom-only takes
   its own branch (`api/apiplan.c:102-107`) that bypasses `mkplan`'s
   forget-everything recovery, and plans with `hash_info = 0`, so it writes
   nothing back.

## Requirements

### Patience

**R1.** Add `NFFT_PATIENT` = `1U << 5`. Reserve `1U << 6` for a future
`NFFT_EXHAUSTIVE`, which is **not** defined publicly. Existing values do not
move: `NFFT_MEASURE = 0`, `NFFT_ESTIMATE = 1<<0`, `NFFT_NO_DIRECT = 1<<1`,
`NFFT_NO_FAST_NATIVE = 1<<4`. Bits `1<<2` and `1<<3` belonged to retired flags
and stay unused.

**R2.** Add `NFFT_NO_NONTHREADED` = `1U << 7`, mirroring FFTW's beyond-guru
flag.

**R3.** The public planning word is normalised through an `IMPLIES`-style table
in its own translation unit, unit-testable without constructing a plan. The
exhaustive branch is written and reads the reserved bit; it ships unreachable
and uncovered until the flag is exposed. That is a known, accepted cost.

**R4.** Add `PLNR_NO_NONTHREADED` and `NO_NONTHREADEDP(pl)` in FFTW's exact
two-condition form. Every solver that does not parallelise internally declines
when it holds: `fast_native`, `ndft_1d`, `ndft_nd`, and every DECONV and CONV
solver. `rnk0` is exempt.

**R5.** `nthr` and the flag bounds are saved and restored around every `mkplan`
the search performs.

**R6. No relaxation table.** No NFFT solver reads `PLNR_NO_UGLY` or
`PLNR_NO_SLOW`, so a port of `relax_tab` would be dead code. Wisdom stays
partitioned across patience levels, as in FFTW. Record, do not fix.

### Packaging

**R7.** Threaded planner code lives in `libnfft3<suffix>_ng_omp`, built only
under `--enable-openmp`, linked **in addition** to `libnfft3<suffix>`. The
existing `libnfft3<suffix>_omp` is a whole-kernel OpenMP rebuild linked
*instead of* the serial library; it keeps that meaning and is not touched.

**R8.** `X(init_threads)`, `X(cleanup_threads)`, `X(plan_with_nthreads)` and
`X(planner_nthreads)` are declared in `nfft3.h` and defined **only** in the
add-on. `X(plan_with_nthreads)` calls `X(init_threads)` first;
`X(init_threads)` destroys the planner before registering the threaded roster.
A program linked against the serial library alone cannot raise the count, which
is why `NO_NONTHREADEDP` needs no third condition.

Destroying the planner discards in-memory wisdom. It is verified not to dangle
live plans: `Y(planner_destroy)` frees `slvdescs` and the tables, and a live
plan holds only file-scope `static const` adt pointers and no `solver*`.
Document `X(plan_with_nthreads)` as a call to make before any other NFFT
routine, as FFTW does.

**R9.** Until a threaded solver exists the roster registers nothing;
`X(init_threads)` must then return 0 and `X(plan_with_nthreads)` must leave the
count at 1.

### Child FFTW plans

**R10.** `problem_nfft.fftw_flags` documents `0 = derive` (`iplanner.h:411`) and
nothing derives. Implement it, three rows:

| NFFT level | derived child flags |
|---|---|
| `NFFT_ESTIMATE` | `FFTW_ESTIMATE` |
| `NFFT_MEASURE` | `FFTW_MEASURE` |
| `NFFT_PATIENT` | `FFTW_PATIENT` |

A non-zero caller word is used as given. `FFTW_DESTROY_INPUT` is forced and
`FFTW_PRESERVE_INPUT` stripped where the child plan is built
(`nfft-nd.c:209-217`), unchanged, and the child plans stay out of place because
`g1 != g2`. Those two properties are not observable from a test and are left to
code review; `nfft-nd.c:210` ships untested.

`FFTW_PATIENT` is what makes FFTW compare its own threaded and serial
candidates, so this is the one place a patience level does real work today.

**R11.** Feed FFTW's thread count into the NFFT problem hash through a
null-by-default hook that the add-on installs, because
`fftw_planner_nthreads` is defined in FFTW's threads library and `libnfft3`
links `@fftw3_LIBS@` only. Probe it with `AC_CHECK_DECL` using the
precision-mangled spelling; without it the add-on is not built and the hook
stays null.

### Wisdom and blessing

**R12.** Blessing reaches the whole winning tree by FFTW's mechanism, not by a
child-enumeration hook. `search_flags()` carries `pl->flags.info & PLNR_BLESSING`
into its inserts. In measured mode a second `Y(planner_mkplan)` pass runs after
the top-level bless with the bit set; every node hits the first pass's memo. In
estimate mode there is no race, so the single pass runs with the bit set.

**R13.** DECONV and CONV children are blessed and exported. The wisdom grammar
is kind-agnostic — entries are keyed by md5 and name a solver by registrar name
and ordinal (`planner.c:307-340`) — so **no format change is required**.

**R14.** Estimate solutions are blessed and exported, reversing
`plan.c:125`. FFTW parity. Estimate-only users start writing wisdom files whose
entries can only ever answer other estimate queries.

**R15.** The configuration signature carries a vocabulary tag, bumped whenever a
`PLNR_*` bit or a problem-hash field changes.

### Wisdom-only

**R16.** Add `NFFT_WISDOM_ONLY` = `1U << 8`, and a `wisdom_state` field on
`planner_s` — **not** a `PLNR_*` bit, because it must not participate in `LEQ`
subsumption or key an entry.

**R17.** Semantics, mirroring FFTW:
- set from the planning word at the top of the guru, reset at the guru
  boundary, so it never persists across calls;
- sticky within one call: a miss sets the bogus state, returns null, and unwinds
  through parents;
- gates every problem kind, children included, which works because R12 and R13
  put children in wisdom;
- takes a branch that bypasses any forget-everything recovery, so a miss can
  never wipe the store;
- plans unblessed and writes nothing back;
- returns `NULL` from the guru on a miss, in every mode.

**R18.** `Y(nfft_derive_fftw_flags)` sets `FFTW_WISDOM_ONLY` in the child word
when `NFFT_WISDOM_ONLY` is given and clears it otherwise, on both the derived
and the caller-supplied path. `FFTW_WISDOM_ONLY` inside `fftw_flags` stops
working; `tests/nplan.c:1495` moves to the planning word.

**R19.** `keyable_fftw_flags` (`plan.c:39-42`) strips `FFTW_WISDOM_ONLY` as well
as the preservation bits. It says how hard to look for a plan, not which plan is
wanted; leaving it in the key would make a wisdom-only attempt look under a
different key from the entry it needs and never find it.

## Deliberate divergences from FFTW

Each must appear in the ADR.

1. `NO_NONTHREADEDP` has two conditions, not three, because packaging supplies
   the invariant FFTW gets from `fftw_plan_with_nthreads`.
2. No `relax_tab`, because no solver reads `NO_UGLY` or `NO_SLOW`.
3. No patience-escalation loop under a timelimit (`api/apiplan.c:108-125`).
   This is the main reason our patience levels do less than FFTW's.
4. FFTW's thread count is observed only when the add-on is linked.
5. `NFFT_EXHAUSTIVE` is reserved rather than exposed, so its mapping branch
   ships dead.

## Out of scope

- Any threaded solver, OpenMP pragma or thread pool. The roster ships empty.
- Racing thread counts; a pthreads variant; the patience-escalation loop.
- Consumers for `PLNR_NO_UGLY`, `PLNR_NO_SLOW`, `PLNR_ALLOW_PRUNING`,
  `PLNR_BELIEVE_PCOST`.
- An error channel distinguishing "no wisdom" from "bad arguments"; the guru's
  `NULL` stays overloaded.
- Legacy `nfft_plan`, `Y(set_num_threads)`, `libnfft3<suffix>_omp`, `nfsft`,
  `fpt`.

## Acceptance

`make check` passes across double, float, long-double and OpenMP configurations.
New CUnit cases cover the mapping table, the derived child flags, the decline
rule, dynamic scoping, the extended problem key, tree blessing, and every clause
of R17. A `checkall_ngomp` binary, built only under `--enable-openmp`, links
`libnfft3 + libnfft3_ng_omp` and covers R8 and R9.
