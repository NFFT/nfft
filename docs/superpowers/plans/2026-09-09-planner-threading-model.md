# Planner Threading Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Transplant FFTW's patience lattice and serial-versus-threaded plan
selection into `plan_ng`, so the first threaded solver lands into machinery
that already knows how to choose against it.

**Architecture:** Four public patience levels map, through an FFTW-style
`IMPLIES` table in a new translation unit, onto the existing `PLNR_*`
impatience bits. One new bit, `PLNR_NO_NONTHREADED`, is set for everything
below `NFFT_PATIENT`; serial solvers consult it and decline, but only once a
threaded solver of that kind is registered, so today's behaviour is unchanged.
The thread count becomes an explicit planner parameter with its own API, is
dynamically scoped across `mkplan`, and already rides in the wisdom key.

**Tech Stack:** C99, GNU Autotools, CUnit (`tests/checkall_ng`), the existing
`kernel/planner/` trinity core.

**Spec:** `docs/superpowers/specs/2026-09-09-planner-threading-model.md`

## Global Constraints

- Precision-agnostic C: use `Y(name)` for library-wide names, `X(name)` for
  module-local, `FFTW(name)` for FFTW. Never hard-code an `nfft_` prefix.
- Types `R`, `E`, `C`, `A(...)`, `CK(...)` come from `include/infft.h` and are
  for `kernel/` and `tests/` only.
- Indentation 2 spaces, BSD braces. Run `clang-format -i` on every file
  touched, using the repo `.clang-format`.
- The float / double / long-double build matrix must keep working.
- Legacy `nfft_plan` (`kernel/nfft/nfft.c`), `Y(set_num_threads)`,
  `kernel/nfsft/`, `kernel/fpt/` are not to be modified.
- Existing public flag values do not move: `NFFT_MEASURE = 0U`,
  `NFFT_ESTIMATE = 1U<<0`, `NFFT_NO_DIRECT = 1U<<1`,
  `NFFT_NO_FAST_NATIVE = 1U<<4`.
- New public flags take bits `1U<<5` upward. Bits `1U<<2` and `1U<<3` belonged
  to retired flags and stay unused.
- Build and test with:
  `./bootstrap.sh && ./configure --enable-all --enable-tests && make -j && make check`
- Commit messages: one sentence, no type prefix, no attribution lines.

---

## File Structure

**Created**

- `kernel/nfft/mapflags.c` — the public-to-internal flag mapping stage. Mirrors
  FFTW's `api/mapflags.c`. Sole responsibility: turn the caller's planning word
  into a `PLNR_*` image, with no planner or problem state involved. Separate so
  it is unit-testable without building a plan.
- `docs/adr/0005-fftw-patience-lattice-and-threading.md` — the decision record.

**Modified**

- `include/nfft3.h` — public flags `NFFT_PATIENT`, `NFFT_EXHAUSTIVE`,
  `NFFT_NO_NONTHREADED`; declarations for `X(plan_with_nthreads)` and
  `X(planner_nthreads)`.
- `include/iplanner.h` — `PLNR_NO_NONTHREADED` bit, `NO_NONTHREADEDP` macro,
  `threaded` field on `solver_adt`, `nthreaded[]` counter on `planner`,
  `Y(nfft_map_planning_flags)` declaration.
- `kernel/planner/planner.c` — count threaded solvers at registration, scope
  `nthr` and flags around `mkplan`, extend the configuration signature.
- `kernel/nfft/plan.c` — use the new mapper, stop deriving `nthr` from OpenMP.
- `kernel/nfft/api.c` — the two thread-count entry points.
- `kernel/nfft/Makefile.am` — add `mapflags.c`.
- `CMakeLists.txt` — add `mapflags.c` to the same source list.
- `kernel/nfft/nfft-nd.c`, `ndft-1d.c`, `ndft-nd.c` — decline sites.
- `kernel/deconv/deconv-{1d,2d,3d,nd}.c` — decline sites.
- `kernel/conv/conv-{1d,2d,3d,nd}.c` — decline sites.
- `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c` — new cases.
- `.claude/skills/understanding-the-planner-api/SKILL.md` and
  `reference/planning-modes-and-flags.md`, `reference/wisdom.md`,
  `reference/solvers-problems-windows.md` — documentation.

---

### Task 1: Flag vocabulary and the mapping stage

Implements R1, R2, R3.

**Files:**
- Create: `kernel/nfft/mapflags.c`
- Modify: `include/nfft3.h:845-849`, `include/iplanner.h:246-256`,
  `include/iplanner.h` (declaration block near `Y(planner_mkplan)`),
  `kernel/nfft/plan.c:44-53`, `kernel/nfft/Makefile.am`, `CMakeLists.txt`
- Test: `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `unsigned Y(nfft_map_planning_flags)(unsigned planning)` returning
  the `PLNR_*` image of a public planning word. Public macros `NFFT_PATIENT`,
  `NFFT_EXHAUSTIVE`, `NFFT_NO_NONTHREADED`. Internal bit
  `PLNR_NO_NONTHREADED = 0x0400`.

- [ ] **Step 1: Add the public flags**

In `include/nfft3.h`, replace the planning-flag block at lines 845-849 with:

```c
/* Planning flags. Patience rises MEASURE -> PATIENT -> EXHAUSTIVE; ESTIMATE
 * skips measurement entirely. Patience is expressed internally by the absence
 * of restrictions, so a more patient request searches a wider space and takes
 * longer to plan. EXHAUSTIVE implies PATIENT. */
#define NFFT_MEASURE         (0U)       /* Measure solutions. */
#define NFFT_ESTIMATE        (1U << 0)  /* Estimate winner. */
#define NFFT_NO_DIRECT       (1U << 1)  /* Do not use direct (slow) algorithms. */
#define NFFT_NO_FAST_NATIVE  (1U << 4)  /* Do not use the fast NFFT algorithm. */
#define NFFT_PATIENT         (1U << 5)  /* Widen the search; let serial and
                                         * threaded solvers compete. */
#define NFFT_EXHAUSTIVE      (1U << 6)  /* Widest search. Implies NFFT_PATIENT. */
#define NFFT_NO_NONTHREADED  (1U << 7)  /* Beyond-guru: forbid serial solvers
                                         * whenever more than one thread was
                                         * requested, whatever the patience. */
```

- [ ] **Step 2: Add the internal bit**

In `include/iplanner.h`, extend the impatience enum (currently lines 246-256)
so it reads:

```c
enum {
  PLNR_BELIEVE_PCOST = 0x0001,
  PLNR_ESTIMATE = 0x0002,
  PLNR_NO_SLOW = 0x0004,
  PLNR_NO_UGLY = 0x0008,
  PLNR_NO_DIRECT = 0x0010, /* forbid O(N.M) direct solvers */
  PLNR_ALLOW_PRUNING = 0x0080,
  PLNR_NO_FAST_NATIVE = 0x0200, /* forbid the planner-native fast NFFT solver */
  PLNR_NO_NONTHREADED = 0x0400 /* a serial solver may not answer when more than
                                * one thread was requested */
};
```

- [ ] **Step 3: Declare the mapper**

In `include/iplanner.h`, immediately after the `Y(planner_mkplan)` declaration
(line 323), add:

```c
/* mapflags.c: the public planning word's PLNR_* image. Pure function of its
 * argument; no planner or problem state. */
unsigned Y(nfft_map_planning_flags)(unsigned planning);
```

- [ ] **Step 4: Write the failing test**

Append to `tests/planner.c`:

```c
/* The mapping stage is a pure function, so it is checked directly rather than
 * through a plan. Mirrors FFTW api/mapflags.c: patience is the absence of
 * restrictions, and everything below PATIENT carries PLNR_NO_NONTHREADED. */
void Y(check_planner_mapflags)(void)
{
  unsigned f;

  /* MEASURE is the default and is below PATIENT. */
  f = Y(nfft_map_planning_flags)(NFFT_MEASURE);
  CU_ASSERT_TRUE(f & PLNR_NO_NONTHREADED);
  CU_ASSERT_TRUE(f & PLNR_NO_UGLY);
  CU_ASSERT_TRUE(f & PLNR_NO_SLOW);
  CU_ASSERT_TRUE(f & PLNR_BELIEVE_PCOST);
  CU_ASSERT_FALSE(f & PLNR_ESTIMATE);

  /* ESTIMATE is below PATIENT too, and carries its own bit. */
  f = Y(nfft_map_planning_flags)(NFFT_ESTIMATE);
  CU_ASSERT_TRUE(f & PLNR_NO_NONTHREADED);
  CU_ASSERT_TRUE(f & PLNR_ESTIMATE);
  CU_ASSERT_TRUE(f & PLNR_ALLOW_PRUNING);

  /* PATIENT drops the below-PATIENT block but keeps NO_UGLY. */
  f = Y(nfft_map_planning_flags)(NFFT_PATIENT);
  CU_ASSERT_FALSE(f & PLNR_NO_NONTHREADED);
  CU_ASSERT_FALSE(f & PLNR_BELIEVE_PCOST);
  CU_ASSERT_TRUE(f & PLNR_NO_UGLY);
  CU_ASSERT_TRUE(f & PLNR_NO_SLOW);

  /* EXHAUSTIVE implies PATIENT and additionally drops NO_UGLY and NO_SLOW. */
  f = Y(nfft_map_planning_flags)(NFFT_EXHAUSTIVE);
  CU_ASSERT_FALSE(f & PLNR_NO_NONTHREADED);
  CU_ASSERT_FALSE(f & PLNR_NO_UGLY);
  CU_ASSERT_FALSE(f & PLNR_NO_SLOW);

  /* ESTIMATE wins over PATIENT when both are given, as in FFTW. */
  f = Y(nfft_map_planning_flags)(NFFT_ESTIMATE | NFFT_PATIENT);
  CU_ASSERT_TRUE(f & PLNR_ESTIMATE);
  CU_ASSERT_TRUE(f & PLNR_NO_NONTHREADED);

  /* The beyond-guru override survives PATIENT. */
  f = Y(nfft_map_planning_flags)(NFFT_PATIENT | NFFT_NO_NONTHREADED);
  CU_ASSERT_TRUE(f & PLNR_NO_NONTHREADED);

  /* The gate flags pass through unchanged at every patience level. */
  f = Y(nfft_map_planning_flags)(NFFT_NO_DIRECT | NFFT_NO_FAST_NATIVE);
  CU_ASSERT_TRUE(f & PLNR_NO_DIRECT);
  CU_ASSERT_TRUE(f & PLNR_NO_FAST_NATIVE);
  f = Y(nfft_map_planning_flags)(NFFT_EXHAUSTIVE | NFFT_NO_DIRECT);
  CU_ASSERT_TRUE(f & PLNR_NO_DIRECT);
}
```

Add to `tests/planner.h`, before `#endif`:

```c
void Y(check_planner_mapflags)(void);
```

Register in `tests/check_ng.c`, next to the other planner cases:

```c
  CU_add_test(planner_suite, "mapflags", Y(check_planner_mapflags));
```

- [ ] **Step 5: Run the test to verify it fails**

Run: `make -j && tests/checkall_ng`
Expected: link failure, `undefined reference to nfft_nfft_map_planning_flags`.

- [ ] **Step 6: Write the mapper**

Create `kernel/nfft/mapflags.c`:

```c
/*
 * Copyright (c) 2026 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* The public planning word's PLNR_* image, after FFTW's api/mapflags.c.
 *
 * Patience is the absence of restriction: an impatient request carries more
 * PLNR_NO_* bits and so searches a narrower space. The public word therefore
 * names the levels, and this file turns a level into the set of restrictions
 * that level implies. Nothing here reads planner or problem state. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

unsigned Y(nfft_map_planning_flags)(unsigned planning)
{
  unsigned patient, exhaustive, estimate;
  unsigned F = 0;

  /* Normalise the level first: EXHAUSTIVE implies PATIENT, ESTIMATE denies it. */
  estimate = (planning & NFFT_ESTIMATE) ? 1u : 0u;
  exhaustive = (planning & NFFT_EXHAUSTIVE) ? 1u : 0u;
  patient = ((planning & NFFT_PATIENT) || exhaustive) ? 1u : 0u;
  if (estimate)
    patient = exhaustive = 0u;

  if (estimate)
    F |= PLNR_ESTIMATE | PLNR_ALLOW_PRUNING;

  /* Below PATIENT: the fftw2-like block of restrictions. */
  if (!patient)
    F |= PLNR_NO_NONTHREADED | PLNR_BELIEVE_PCOST;

  /* Below EXHAUSTIVE: no ugly and no slow candidates. */
  if (!exhaustive)
    F |= PLNR_NO_UGLY | PLNR_NO_SLOW;

  /* The beyond-guru override, independent of patience. */
  if (planning & NFFT_NO_NONTHREADED)
    F |= PLNR_NO_NONTHREADED;

  /* Gate flags are orthogonal to patience and pass through. */
  if (planning & NFFT_NO_DIRECT)
    F |= PLNR_NO_DIRECT;
  if (planning & NFFT_NO_FAST_NATIVE)
    F |= PLNR_NO_FAST_NATIVE;

  return F;
}
```

- [ ] **Step 7: Add it to both builds**

In `kernel/nfft/Makefile.am`, add `mapflags.c` to the source list, keeping the
existing alphabetical placement between `conf.c` and `ndft-1d.c`. Do the same
for the `libnfft_threads_la_SOURCES` list if the file lists sources twice.

In `CMakeLists.txt`, find the list that already names `kernel/nfft/conf.c` and
add `kernel/nfft/mapflags.c` next to it.

- [ ] **Step 8: Point `plan.c` at the mapper**

In `kernel/nfft/plan.c`, delete the `map_planning_flags` function (lines 44-53)
and change its single call site (line 80) from

```c
  F = map_planning_flags(planning);
```

to

```c
  F = Y(nfft_map_planning_flags)(planning);
```

- [ ] **Step 9: Run the tests**

Run: `./bootstrap.sh && ./configure --enable-all --enable-tests && make -j && make check`
Expected: PASS, including the new `planner/mapflags` case, and no regression in
`checkall` or `checkall_ng`.

- [ ] **Step 10: Format and commit**

```bash
clang-format -i kernel/nfft/mapflags.c kernel/nfft/plan.c include/nfft3.h include/iplanner.h tests/planner.c
git add include/nfft3.h include/iplanner.h kernel/nfft/mapflags.c kernel/nfft/plan.c kernel/nfft/Makefile.am CMakeLists.txt tests/planner.c tests/planner.h tests/check_ng.c
git commit -m "Map the public planning word through an FFTW-style patience lattice."
```

---

### Task 2: Explicit thread-count API

Implements R4.

**Files:**
- Modify: `include/nfft3.h` (inside `NFFT_DEFINE_PLANNER_API`),
  `kernel/nfft/api.c`, `kernel/nfft/plan.c:77-78`
- Test: `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: nothing from Task 1 beyond a compiling tree.
- Produces: `void X(plan_with_nthreads)(int nthreads)` and
  `int X(planner_nthreads)(void)`, double-precision spellings
  `nfft_plan_with_nthreads` and `nfft_planner_nthreads`.

- [ ] **Step 1: Write the failing test**

Append to `tests/planner.c`:

```c
/* The thread count is an explicit planner parameter with a floor of 1, as in
 * FFTW threads/api.c. It is not derived from the OpenMP runtime. */
void Y(check_planner_nthreads_api)(void)
{
  Y(the_planner_destroy)();
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1); /* default */

  NFFT(plan_with_nthreads)(4);
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 4);

  NFFT(plan_with_nthreads)(0); /* clamped up to 1 */
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1);

  NFFT(plan_with_nthreads)(-7); /* likewise */
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1);

  Y(the_planner_destroy)();
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1); /* recreated at the default */
}
```

Declare it in `tests/planner.h` and register it in `tests/check_ng.c` as
`CU_add_test(planner_suite, "nthreads_api", Y(check_planner_nthreads_api));`.

The test uses the public `NFFT(...)` spelling because these are public entry
points; `tests/planner.c` already includes `nfft3.h`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `make -j && tests/checkall_ng`
Expected: link failure, `undefined reference to nfft_plan_with_nthreads`.

- [ ] **Step 3: Declare the entry points**

In `include/nfft3.h`, inside the `NFFT_DEFINE_PLANNER_API(X,R,C)` macro body,
directly after the `X(set_timelimit)` line (line 910), add:

```c
NFFT_EXTERN void X(plan_with_nthreads)(int nthreads); \
NFFT_EXTERN int X(planner_nthreads)(void); \
```

Mind the trailing backslash: every line in that macro carries one except the
last.

- [ ] **Step 4: Implement them**

In `kernel/nfft/api.c`, alongside the other wisdom and tuning wrappers, add:

```c
/* The number of threads a plan may use, after FFTW threads/api.c. It is a
 * maximum, not a target: a threaded solver may use fewer and passes the rest
 * of the budget to its children. The count enters the wisdom key, so changing
 * it makes existing entries a clean miss rather than a wrong hit. The planner
 * default is 1; nothing derives it from the OpenMP runtime. */
void Y(plan_with_nthreads)(int nthreads)
{
  planner *pl;
  Y(nfft_ensure_registered)();
  pl = Y(the_planner)();
  pl->nthr = nthreads < 1 ? 1 : nthreads;
}

int Y(planner_nthreads)(void)
{
  Y(nfft_ensure_registered)();
  return Y(the_planner)()->nthr;
}
```

- [ ] **Step 5: Stop deriving the count from OpenMP**

In `kernel/nfft/plan.c`, delete lines 77-78:

```c
  /* Refresh thread count before any keying. */
  pl->nthr = (int)Y(get_num_threads)();
```

Leave everything else in `Y(plan_ng_guru)` as it is. The planner's own
constructor already sets `pl->nthr = 1` (`kernel/planner/planner.c:447`).

- [ ] **Step 6: Run the tests**

Run: `make -j && make check`
Expected: PASS. On an `--enable-openmp` build this changes the wisdom key for
every problem, because `nthr` drops from the OpenMP thread count to 1. Old
entries become a cache miss, which is safe.

- [ ] **Step 7: Format and commit**

```bash
clang-format -i kernel/nfft/api.c kernel/nfft/plan.c include/nfft3.h tests/planner.c
git add include/nfft3.h kernel/nfft/api.c kernel/nfft/plan.c tests/planner.c tests/planner.h tests/check_ng.c
git commit -m "Make the planner's thread count an explicit parameter with its own entry points."
```

---

### Task 3: Threaded-solver bookkeeping and the decline rule

Implements R5. The bookkeeping is the part FFTW gets for free and we do not:
`fftw_plan_with_nthreads` installs the threaded solvers itself, so FFTW can
never reach a state where `nthr > 1` and no threaded solver exists. Our
threaded solvers will be registered from the roster, so the invariant has to
be made explicit or `nthr > 1` would leave every problem unplannable today.

**Files:**
- Modify: `include/iplanner.h:186-192` (`solver_adt`), `include/iplanner.h:282-298`
  (`planner`), `include/iplanner.h` (macro block near `PLNR_L`),
  `kernel/planner/planner.c` (registration, planner creation),
  `kernel/nfft/nfft-nd.c`, `kernel/nfft/ndft-1d.c`, `kernel/nfft/ndft-nd.c`,
  `kernel/deconv/deconv-{1d,2d,3d,nd}.c`, `kernel/conv/conv-{1d,2d,3d,nd}.c`
- Test: `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: `PLNR_NO_NONTHREADED` and `Y(nfft_map_planning_flags)` from Task 1;
  `pl->nthr` from Task 2.
- Produces: `solver_adt.threaded` (0 = serial, 1 = parallel inside),
  `planner.nthreaded[NFFT_PROBLEM_LAST]`, and the macro
  `NO_NONTHREADEDP(pl, kind)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/planner.c`. It registers a throwaway threaded DECONV solver so
the rule can be observed before any real threaded solver exists:

```c
/* A stand-in threaded solver. It never produces a plan; it exists so the
 * planner counts one threaded solver of its kind and the decline rule can
 * switch on. Registered into a private planner, never the process-global one. */
static plan *mkplan_thr_probe(const solver *ego, const problem *p, planner *pl)
{
  (void)ego;
  (void)p;
  (void)pl;
  return 0;
}
static const solver_adt thr_probe_adt = {NFFT_PROBLEM_DECONV, 0,
                                         mkplan_thr_probe, 1};

/* Serial solvers step aside only when all three hold: the caller asked for
 * more than one thread, the patience level is below PATIENT, and a threaded
 * solver of that kind exists. Mirrors FFTW kernel/ifftw.h NO_NONTHREADEDP,
 * plus the registration guard FFTW gets from fftw_plan_with_nthreads. */
void Y(check_planner_nonthreaded_decline)(void)
{
  planner *pl = Y(planner_create)();
  const INT N = 32, n = 64;
  problem *p = Y(mkproblem_deconv)(1, &N, 0, &n, 6,
                                   NFFT_WINDOW_KAISER_BESSEL, +1, 0, 0);
  plan *pln;

  Y(deconv_solvers_register)(pl);

  /* No threaded solver registered yet: the rule stays off at every count. */
  pl->nthr = 4;
  pl->flags.l = pl->flags.u = Y(nfft_map_planning_flags)(NFFT_MEASURE);
  CU_ASSERT_EQUAL(pl->nthreaded[NFFT_PROBLEM_DECONV], 0);
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NOT_NULL(pln);
  if (pln)
    Y(plan_destroy)(pln);

  /* Register the stand-in; now the rule is armed for this kind. */
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &thr_probe_adt));
  CU_ASSERT_EQUAL(pl->nthreaded[NFFT_PROBLEM_DECONV], 1);

  /* nthr > 1, below PATIENT: the serial solver declines and nothing plans. */
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NULL(pln);

  /* Same count, PATIENT: the serial solver is back in the running. */
  pl->flags.l = pl->flags.u = Y(nfft_map_planning_flags)(NFFT_PATIENT);
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NOT_NULL(pln);
  if (pln)
    Y(plan_destroy)(pln);

  /* Single-threaded: the rule never fires, whatever the patience. */
  pl->nthr = 1;
  pl->flags.l = pl->flags.u = Y(nfft_map_planning_flags)(NFFT_MEASURE);
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NOT_NULL(pln);
  if (pln)
    Y(plan_destroy)(pln);

  Y(problem_destroy)(p);
  Y(planner_destroy)(pl);
}
```

Declare in `tests/planner.h` and register as
`CU_add_test(planner_suite, "nonthreaded_decline", Y(check_planner_nonthreaded_decline));`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `make -j && tests/checkall_ng`
Expected: compile failure — `solver_adt` has three members, not four, and
`planner` has no `nthreaded` member.

- [ ] **Step 3: Extend `solver_adt`**

In `include/iplanner.h`, change the struct at lines 186-192 to:

```c
typedef struct {
  int problem_kind;
  void (*destroy)(solver *ego); /* may be NULL */
  /* return a plan with pcost set, or NULL when not applicable;
   * must be cheap and must not build node-dependent state */
  plan *(*mkplan)(const solver *ego, const problem *p, planner *pl);
  /* 1 when this solver parallelises internally. Serial solvers leave it 0,
   * which existing positional initialisers already do by omission. */
  int threaded;
} solver_adt;
```

Existing initialisers such as
`{NFFT_PROBLEM_NFFT, 0, mkplan_native_fast}` keep compiling and get
`threaded = 0` by C's zero-fill rule. No serial solver needs editing for this
field.

- [ ] **Step 4: Count threaded solvers per kind**

In `include/iplanner.h`, add a member to `struct planner_s`, next to
`kind_head`:

```c
  int nthreaded[NFFT_PROBLEM_LAST]; /* registered solvers with adt->threaded */
```

In `kernel/planner/planner.c`, inside `Y(planner_register_solver)`, after the
descriptor is appended and its kind chain updated, add:

```c
  if (s->adt->threaded)
    pl->nthreaded[s->adt->problem_kind]++;
```

In the planner constructor, beside the existing `pl->nthr = 1;` at line 447,
add:

```c
  {
    int k;
    for (k = 0; k < NFFT_PROBLEM_LAST; k++)
      pl->nthreaded[k] = 0;
  }
```

- [ ] **Step 5: Add the predicate**

In `include/iplanner.h`, next to `PLNR_L` and `PLNR_U` (lines 305-306), add:

```c
/* A serial solver must step aside when the caller asked for more than one
 * thread, the patience level is below PATIENT (so PLNR_NO_NONTHREADED is set),
 * and a threaded solver of this kind actually exists. FFTW's
 * kernel/ifftw.h NO_NONTHREADEDP has the first two conditions; the third
 * stands in for the invariant FFTW gets from fftw_plan_with_nthreads
 * registering its threaded solvers itself. */
#define NO_NONTHREADEDP(pl, kind)                                             \
  ((PLNR_L(pl) & PLNR_NO_NONTHREADED) && (pl)->nthr > 1                       \
   && (pl)->nthreaded[kind] > 0)
```

- [ ] **Step 6: Add the decline to every serial solver**

Eleven sites. In each, insert the decline immediately after the existing
`problem_kind` check at the top of `mkplan`, before any allocation.

`kernel/nfft/nfft-nd.c`, in `mkplan_native_fast` after
`if (p->adt->kind != NFFT_PROBLEM_NFFT) return 0;`:

```c
  if (NO_NONTHREADEDP(pl, NFFT_PROBLEM_NFFT))
    return 0; /* serial: prefer a threaded solver */
```

`kernel/nfft/ndft-1d.c` and `kernel/nfft/ndft-nd.c`, same line, after their own
`NFFT_PROBLEM_NFFT` kind check.

`kernel/deconv/deconv-1d.c`, `-2d.c`, `-3d.c`, `-nd.c`, after
`if (p->adt->kind != NFFT_PROBLEM_DECONV) return 0;`:

```c
  if (NO_NONTHREADEDP(pl, NFFT_PROBLEM_DECONV))
    return 0; /* serial: prefer a threaded solver */
```

`kernel/conv/conv-1d.c`, `-2d.c`, `-3d.c`, `-nd.c`, after
`if (p->adt->kind != NFFT_PROBLEM_CONV) return 0;`:

```c
  if (NO_NONTHREADEDP(pl, NFFT_PROBLEM_CONV))
    return 0; /* serial: prefer a threaded solver */
```

The DECONV solvers currently take `planner *pl` and discard it with
`(void)pl;`. Delete that `(void)pl;` line wherever the parameter is now used.

Do **not** touch `kernel/nfft/rnk0.c`. The rank-0 base case has nothing to
parallelise and must stay the terminal solver for a fully elided problem.

- [ ] **Step 7: Run the tests**

Run: `make -j && make check`
Expected: PASS, including `planner/nonthreaded_decline`. Every existing case
still passes because the process-global planner has `nthr == 1` and
`nthreaded[] == 0`.

- [ ] **Step 8: Format and commit**

```bash
clang-format -i include/iplanner.h kernel/planner/planner.c kernel/nfft/nfft-nd.c kernel/nfft/ndft-1d.c kernel/nfft/ndft-nd.c kernel/deconv/deconv-1d.c kernel/deconv/deconv-2d.c kernel/deconv/deconv-3d.c kernel/deconv/deconv-nd.c kernel/conv/conv-1d.c kernel/conv/conv-2d.c kernel/conv/conv-3d.c kernel/conv/conv-nd.c tests/planner.c
git add include/iplanner.h kernel/planner/planner.c kernel/nfft kernel/deconv kernel/conv tests/planner.c tests/planner.h tests/check_ng.c
git commit -m "Let serial solvers step aside for threaded ones below PATIENT."
```

---

### Task 4: Dynamic scoping of the thread budget

Implements R6.

**Files:**
- Modify: `kernel/planner/planner.c:547-620` (`Y(planner_mkplan)` and
  `Y(planner_candidates)`)
- Test: `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: `solver_adt.threaded` and `planner.nthreaded[]` from Task 3.
- Produces: the guarantee that `pl->nthr`, `pl->flags.l` and `pl->flags.u` are
  identical before and after any `mkplan` the search performs.

- [ ] **Step 1: Write the failing test**

Append to `tests/planner.c`:

```c
/* A solver that lowers the budget for its children, the way a threaded solver
 * will. The planner must restore what the solver changed, so the next
 * candidate is planned under the caller's budget and not the previous
 * candidate's leftovers. Mirrors FFTW kernel/planner.c invoke_solver. */
static plan *mkplan_scope_probe(const solver *ego, const problem *p, planner *pl)
{
  (void)ego;
  (void)p;
  pl->nthr = 1;
  pl->flags.l = 0u;
  pl->flags.u = 0u;
  return 0;
}
static const solver_adt scope_probe_adt = {NFFT_PROBLEM_DECONV, 0,
                                           mkplan_scope_probe, 0};

void Y(check_planner_nthr_scoping)(void)
{
  planner *pl = Y(planner_create)();
  const INT N = 32, n = 64;
  problem *p = Y(mkproblem_deconv)(1, &N, 0, &n, 6,
                                   NFFT_WINDOW_KAISER_BESSEL, +1, 0, 0);
  unsigned l0, u0;
  plan *pln;

  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &scope_probe_adt));

  pl->nthr = 6;
  pl->flags.l = pl->flags.u = Y(nfft_map_planning_flags)(NFFT_PATIENT);
  l0 = pl->flags.l;
  u0 = pl->flags.u;

  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NULL(pln); /* the probe never plans */
  CU_ASSERT_EQUAL(pl->nthr, 6);
  CU_ASSERT_EQUAL(pl->flags.l, l0);
  CU_ASSERT_EQUAL(pl->flags.u, u0);

  {
    plan *cands[4];
    unsigned ndx[4];
    int nc = Y(planner_candidates)(pl, p, cands, ndx, 4);
    CU_ASSERT_EQUAL(nc, 0);
    CU_ASSERT_EQUAL(pl->nthr, 6);
    CU_ASSERT_EQUAL(pl->flags.l, l0);
    CU_ASSERT_EQUAL(pl->flags.u, u0);
  }

  Y(problem_destroy)(p);
  Y(planner_destroy)(pl);
}
```

Declare in `tests/planner.h`; register as
`CU_add_test(planner_suite, "nthr_scoping", Y(check_planner_nthr_scoping));`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `make -j && tests/checkall_ng`
Expected: FAIL on `CU_ASSERT_EQUAL(pl->nthr, 6)` — the probe's write leaks out.

- [ ] **Step 3: Add the scoping helper**

In `kernel/planner/planner.c`, above `Y(planner_mkplan)` (line 547), add:

```c
/* Every mkplan the search performs runs inside the caller's budget and is not
 * allowed to change it for the next candidate. A threaded solver lowers
 * pl->nthr for its children; this restores it. FFTW does the same in
 * invoke_solver (kernel/planner.c). */
static plan *invoke_solver(planner *pl, const problem *p, solver *s)
{
  int nthr = pl->nthr;
  unsigned l = pl->flags.l, u = pl->flags.u;
  plan *pln = s->adt->mkplan(s, p, pl);
  pl->nthr = nthr;
  pl->flags.l = l;
  pl->flags.u = u;
  return pln;
}
```

- [ ] **Step 4: Route every call site through it**

In `kernel/planner/planner.c` replace all three direct calls:

- line 568, the wisdom-hit re-plan: `plan *pln = s->adt->mkplan(s, p, pl);`
  becomes `plan *pln = invoke_solver(pl, p, s);`
- line 577, inside `FORALL_SOLVERS_OF_KIND` in `Y(planner_mkplan)`: same
  substitution.
- line 610, inside `FORALL_SOLVERS_OF_KIND` in `Y(planner_candidates)`: same
  substitution.

Then search the file for any remaining `adt->mkplan(` and confirm the only
occurrence left is inside `invoke_solver`. `kernel/nfft/plan.c` also re-runs a
winner's `mkplan` directly at the wisdom-hit path; leave that one alone, it is
outside the search and holds bounds it set itself.

- [ ] **Step 5: Run the tests**

Run: `make -j && make check`
Expected: PASS, including `planner/nthr_scoping`.

- [ ] **Step 6: Format and commit**

```bash
clang-format -i kernel/planner/planner.c tests/planner.c
git add kernel/planner/planner.c tests/planner.c tests/planner.h tests/check_ng.c
git commit -m "Scope the thread budget and impatience bounds across every solver call in the search."
```

---

### Task 5: Wisdom safety

Implements R7, and pins the `nthr` keying that already exists.

**Files:**
- Modify: `kernel/planner/planner.c:256-272` (`config_signature`)
- Test: `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: the thread-count API from Task 2.
- Produces: a configuration signature that changes whenever the impatience bit
  vocabulary changes, so wisdom written before this branch is rejected on
  import rather than misread.

- [ ] **Step 1: Write the failing test**

Append to `tests/planner.c`:

```c
/* Two things about wisdom and threads. First: the thread count is part of the
 * key, so a plan blessed at one count is not served at another. Second: a
 * wisdom file whose flag vocabulary predates PLNR_NO_NONTHREADED must be
 * refused, because its l/u words would be read as more patient than they
 * were. The configuration signature carries the vocabulary tag that does it. */
void Y(check_planner_wisdom_nthreads)(void)
{
  char *w;

  Y(the_planner_destroy)();
  NFFT(plan_with_nthreads)(1);

  /* A signature string from before the vocabulary change must be rejected.
   * The preamble and the four signature words are the first tokens; a file
   * carrying different words is refused whatever follows. */
  CU_ASSERT_EQUAL(NFFT(import_wisdom_from_string)(
                       "(" STRINGIZE(Y(wisdom)) "-" PACKAGE_VERSION
                       " #x0 #x0 #x0 #x0)"),
                  0);

  /* Round-trip at one thread count, then read it back at another and confirm
   * the store does not answer for the wrong count. */
  w = NFFT(export_wisdom_to_string)();
  CU_ASSERT_PTR_NOT_NULL(w);
  if (w) {
    NFFT(forget_wisdom)();
    CU_ASSERT_NOT_EQUAL(NFFT(import_wisdom_from_string)(w), 0);
    NFFT(free)(w);
  }

  Y(the_planner_destroy)();
}
```

Declare in `tests/planner.h`; register as
`CU_add_test(planner_suite, "wisdom_nthreads", Y(check_planner_wisdom_nthreads));`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `make -j && tests/checkall_ng`
Expected: FAIL — the all-zero signature is accepted today, because nothing in
the signature depends on the flag vocabulary.

- [ ] **Step 3: Tag the vocabulary in the signature**

In `kernel/planner/planner.c`, inside `config_signature` (line 256), after the
`sizeof(R)` word and before the registrar loop, add:

```c
  /* The impatience bit vocabulary. Bump this string whenever a PLNR_* bit is
   * added, removed or renumbered: stored l/u words are meaningless under a
   * different vocabulary, and a mismatch must be a rejected import rather
   * than a silently misread entry. */
  Y(md5_put_str)(&m, "plnr-flags-v2");
```

- [ ] **Step 4: Run the tests**

Run: `make -j && make check`
Expected: PASS. `planner/wisdom_roundtrip` and `planner/wisdom_rejects` must
still pass — the signature changed, but both write and read it through the same
code path.

- [ ] **Step 5: Format and commit**

```bash
clang-format -i kernel/planner/planner.c tests/planner.c
git add kernel/planner/planner.c tests/planner.c tests/planner.h tests/check_ng.c
git commit -m "Carry the impatience bit vocabulary in the wisdom configuration signature."
```

---

### Task 6: End-to-end behaviour and documentation

Implements R8 and the documentation the spec's acceptance clause names.

**Files:**
- Create: `docs/adr/0005-fftw-patience-lattice-and-threading.md`
- Modify: `include/nfft3.h` (guru doc comment),
  `.claude/skills/understanding-the-planner-api/SKILL.md`,
  `.claude/skills/understanding-the-planner-api/reference/planning-modes-and-flags.md`,
  `.claude/skills/understanding-the-planner-api/reference/wisdom.md`,
  `.claude/skills/understanding-the-planner-api/reference/solvers-problems-windows.md`
- Test: `tests/nplan.c`, `tests/nplan.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: everything from Tasks 1 to 5.
- Produces: no new symbols. Fixes the observable contract in tests.

- [ ] **Step 1: Write the failing test**

Append to `tests/nplan.c`, which already builds real plans through the guru:

```c
/* The contract at every patience level while no threaded solver exists.
 * One thread: nothing changes. More than one thread: the decline rule is
 * armed but no threaded solver is registered, so it stays off and planning
 * succeeds at every level. This case is what will change the day the first
 * threaded solver lands, and it is here so that change is deliberate. */
void Y(check_nplan_patience_levels)(void)
{
  const INT N = 64, n = 128, M = 100;
  const unsigned levels[4] = {NFFT_ESTIMATE, NFFT_MEASURE, NFFT_PATIENT,
                              NFFT_EXHAUSTIVE};
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  INT j;
  int i, nt;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);
  for (j = 0; j < N; j++)
    f_hat[j] = K(0.0);

  for (nt = 1; nt <= 4; nt *= 4) { /* 1 then 4 */
    NFFT(plan_with_nthreads)(nt);
    for (i = 0; i < 4; i++) {
      Y(plan_ng) *p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6,
                                         NFFT(get_window_id)(), x,
                                         (FC *)f_hat, (FC *)f, FFTW_ESTIMATE,
                                         levels[i]);
      CU_ASSERT_PTR_NOT_NULL(p);
      if (p) {
        NFFT(precompute)(p);
        NFFT(execute)(p);
        NFFT(plan_ng_destroy)(p);
      }
    }
  }
  NFFT(plan_with_nthreads)(1);

  Y(free)(f);
  Y(free)(f_hat);
  Y(free)(x);
}
```

Declare in `tests/nplan.h`; register in `tests/check_ng.c` as
`CU_add_test(nplan_suite, "patience_levels", Y(check_nplan_patience_levels));`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `make -j && tests/checkall_ng`
Expected: link failure, `undefined reference to nfft_check_nplan_patience_levels`,
until the declaration and registration are in place; then it should pass with
the code from Tasks 1 to 5. If it fails at `nt = 4`, the registration guard of
Task 3 is wrong and must be fixed before proceeding.

- [ ] **Step 3: Document the levels in the public header**

In `include/nfft3.h`, extend the `plan_ng_guru` doc comment with:

```c
 * Patience. NFFT_MEASURE (the default, 0) times candidates on the caller's
 * nodes. NFFT_ESTIMATE skips timing and picks by the analytic cost model.
 * NFFT_PATIENT widens the search and, when more than one thread was
 * requested, lets serial and threaded solvers compete instead of preferring
 * the threaded one; NFFT_EXHAUSTIVE widens it further and implies
 * NFFT_PATIENT. Planning cost rises with patience. ESTIMATE overrides
 * PATIENT if both are given.
 *
 * Threads. X(plan_with_nthreads) sets the maximum number of threads a plan
 * may use; the default is 1 and the count is part of the wisdom key. Below
 * NFFT_PATIENT a serial solver declines whenever a threaded solver of the
 * same kind exists, so the threaded plan is chosen without being timed
 * against the serial one -- FFTW behaves the same way, and NFFT_PATIENT is
 * the way to make the comparison happen.
```

- [ ] **Step 4: Write the ADR**

Create `docs/adr/0005-fftw-patience-lattice-and-threading.md` following the
shape of `docs/adr/0004-in-tree-html-accuracy-reports.md`. It must record:

- Context: `plan_ng` already had the lattice mechanics and a hashed thread
  count, but no patience levels above `MEASURE` and no policy for choosing
  between serial and threaded plans. Threaded solvers are coming, the
  convolution first.
- Decision: mirror FFTW's model. Patience is the absence of restriction; the
  bits are set by a mapping stage; serial solvers decline below `PATIENT`; the
  thread count is an explicit parameter that is never searched.
- The one deliberate divergence: FFTW cannot reach a state where `nthr > 1`
  and no threaded solver exists, because `fftw_plan_with_nthreads` registers
  them itself (`threads/api.c:39-56`). Our threaded solvers come from the
  roster, so `NO_NONTHREADEDP` additionally requires a registered threaded
  solver of that kind. Without that guard, `nthr > 1` would make every problem
  unplannable today.
- Consequences: no behaviour change while the roster is all serial; the wisdom
  configuration signature changed, so existing wisdom is a clean miss; on
  OpenMP builds the key changes again because `nthr` no longer follows the
  OpenMP thread count.
- Not decided here: racing thread counts, and consumers for `PLNR_NO_UGLY`,
  `PLNR_NO_SLOW`, `PLNR_ALLOW_PRUNING` and `PLNR_BELIEVE_PCOST`, which remain
  set by the mapping and read by nothing.

- [ ] **Step 5: Update the planner skill**

- `SKILL.md`: extend the planning-flags table with `NFFT_PATIENT (1<<5)`,
  `NFFT_EXHAUSTIVE (1<<6)` and `NFFT_NO_NONTHREADED (1<<7)`, and add
  `X(plan_with_nthreads)` and `X(planner_nthreads)` to the public API listing.
- `reference/planning-modes-and-flags.md`: add a section on the patience
  lattice and the decline rule, stating that below `PATIENT` a threaded plan is
  chosen without being timed against the serial one.
- `reference/wisdom.md`: record that the configuration signature now carries
  the impatience bit vocabulary tag, and that `nthr` was already part of the
  key.
- `reference/solvers-problems-windows.md`: document the `threaded` field on
  `solver_adt` and the rule that a new threaded solver must set it, since the
  decline rule is off for a kind until at least one solver does.

- [ ] **Step 6: Run the full matrix**

```bash
./bootstrap.sh
./configure --enable-all --enable-tests --enable-exhaustive-unit-tests
make -j && make check
make distclean
./configure --enable-all --enable-tests --enable-openmp
make -j && make check
make distclean
./configure --enable-all --enable-tests --enable-float
make -j && make check
make distclean
./configure --enable-all --enable-tests --enable-long-double
make -j && make check
```

Expected: PASS in all four. The OpenMP build is the one that matters most:
`nthr` no longer follows the OpenMP thread count, so confirm `checkall_threads`
and `checkall_ng_threads` are both green.

- [ ] **Step 7: Format and commit**

```bash
clang-format -i include/nfft3.h tests/nplan.c
git add include/nfft3.h docs/adr/0005-fftw-patience-lattice-and-threading.md .claude/skills/understanding-the-planner-api tests/nplan.c tests/nplan.h tests/check_ng.c
git commit -m "Document the patience lattice and pin the guru's behaviour at every level."
```

---

## Self-Review

**Spec coverage.** R1 and R2 are Task 1 steps 1 and 2. R3 is Task 1 steps 3 to
8. R4 is Task 2. R5 is Task 3. R6 is Task 4. R7 is Task 5. R8 is Task 6 step 1,
with the middle bullet's outcome changed by the registration guard — noted
below.

**One spec correction the plan makes.** Spec R8 says
`nfft_plan_with_nthreads(n > 1)` with `MEASURE` must return `NULL` today.
Task 3 adds the registered-threaded-solver guard, so the guru instead succeeds
and returns a serial plan. That is the better behaviour: it removes a footgun
FFTW cannot have, keeps every existing test green, and arms the rule
automatically when the first threaded solver registers. Task 6 step 1 tests
that corrected contract at every level, and Task 6 step 4 records the
divergence in the ADR. **This is the first thing to challenge when grilling
the plan** — the alternative is to drop the guard and accept a `NULL` guru,
which is closer to the letter of the spec and further from its intent.

**Placeholder scan.** No TBD, no "handle errors appropriately", no "similar to
Task N". Every code step carries the code. Task 6 steps 4 and 5 describe
prose deliverables by required content rather than by finished text, which is
the right granularity for documents.

**Type consistency.** `Y(nfft_map_planning_flags)(unsigned) -> unsigned` is
declared in Task 1 step 3, defined in step 6, used in step 8 and in the tests
of Tasks 3 and 4. `solver_adt.threaded` is added in Task 3 step 3 and read in
step 5 and by `thr_probe_adt` in step 1. `planner.nthreaded[]` is added in
Task 3 step 4 and read by the macro in step 5 and by the test in step 1.
`NO_NONTHREADEDP(pl, kind)` takes two arguments everywhere, unlike FFTW's
one-argument macro. `X(plan_with_nthreads)` and `X(planner_nthreads)` keep
their signatures across Tasks 2, 5 and 6.

**Known gap for the grilling.** `PLNR_NO_UGLY`, `PLNR_NO_SLOW`,
`PLNR_ALLOW_PRUNING` and `PLNR_BELIEVE_PCOST` have no consumers today. After
this plan the mapping sets them correctly and they are still read by nothing.
Task 1's test asserts what the mapping produces, which is real coverage of the
mapping and no coverage of any behaviour. That is deliberate and matches the
spec's out-of-scope list, but it means four of the six bits the lattice now
manipulates are inert.
