# Planner Threading Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Transplant FFTW's patience lattice, serial-versus-threaded plan
selection, add-on threading library, tree-wide blessing and wisdom-only state
into `plan_ng`, so the first threaded solver lands into finished machinery.

**Architecture:** `NFFT_PATIENT` maps, through an `IMPLIES` table in a new
translation unit, onto `PLNR_*` bits including a new `PLNR_NO_NONTHREADED`;
serial solvers consult it and decline. Safety comes from packaging: the
thread-count API lives only in `libnfft3_ng_omp`, so a program linked against
the serial library cannot raise the count — FFTW's own invariant. The same unit
derives the child FFTW flags from patience. Blessing reaches the whole winning
tree by FFTW's second-pass trick rather than any new hook, which is what makes a
recursive `wisdom_state` gate workable.

**Tech Stack:** C99, GNU Autotools, CUnit (`tests/checkall_ng`), the
`kernel/planner/` trinity core, FFTW3 (3.3.9+ for the add-on).

**Spec:** `docs/superpowers/specs/2026-09-09-planner-threading-model.md`

## Global Constraints

- Precision-agnostic C: `Y(name)` library-wide, `X(name)` module-local,
  `FFTW(name)` for FFTW. Never hard-code an `nfft_` prefix.
- `R`, `E`, `C`, `A(...)`, `CK(...)` from `include/infft.h`, `kernel/` and
  `tests/` only.
- 2-space indent, BSD braces, `clang-format -i` every file touched.
- The float / double / long-double matrix must keep working.
- Do not modify `kernel/nfft/nfft.c`, `Y(set_num_threads)`,
  `libnfft3<suffix>_omp`, `kernel/nfsft/`, `kernel/fpt/`.
- Public flag values do not move: `NFFT_MEASURE = 0U`,
  `NFFT_ESTIMATE = 1U<<0`, `NFFT_NO_DIRECT = 1U<<1`,
  `NFFT_NO_FAST_NATIVE = 1U<<4`. `1U<<2` and `1U<<3` stay unused.
- `libnfft3<suffix>` links `@fftw3_LIBS@` only and must never reference a
  symbol from FFTW's threads library.
- Baseline: `./bootstrap.sh && ./configure --enable-all --enable-tests && make -j && make check`
- Commit messages: one sentence, no type prefix, no attribution lines.

## File Structure

**Created**

- `kernel/nfft/mapflags.c` — translation of the caller's request: the `PLNR_*`
  image of the planning word, and the child FFTW planner flags. Mirrors
  `api/mapflags.c`. Separate so both are unit-testable without a plan.
- `kernel/threads/api.c`, `kernel/threads/conf.c`, `kernel/threads/Makefile.am`
  — the add-on library.
- `tests/ngomp.c`, `tests/ngomp.h`, `tests/check_ngomp.c` — its test binary.
- `docs/adr/0005-fftw-patience-lattice-and-threading.md`

**Modified**

`include/nfft3.h`, `include/iplanner.h`, `kernel/planner/planner.c`,
`kernel/nfft/{plan.c,problem.c,nfft-nd.c,ndft-1d.c,ndft-nd.c,Makefile.am}`,
`kernel/deconv/deconv-{1d,2d,3d,nd}.c`, `kernel/conv/conv-{1d,2d,3d,nd}.c`,
`configure.ac`, `Makefile.am`, `kernel/Makefile.am`, `tests/Makefile.am`,
`CMakeLists.txt`, `tests/{planner.c,planner.h,nplan.c,nplan.h,check_ng.c}`,
and the planner skill under `.claude/skills/understanding-the-planner-api/`.

---

### Task 1: Flag vocabulary, the mapping stage, and the derived child flags

Implements R1, R2, R3, R10, R18, R19.

**Files:**
- Create: `kernel/nfft/mapflags.c`
- Modify: `include/nfft3.h:845-849`, `include/iplanner.h:246-256`,
  `include/iplanner.h:323`, `kernel/nfft/plan.c:34-53`, `:80`, `:105-107`,
  `kernel/nfft/Makefile.am`, `CMakeLists.txt`, `tests/nplan.c:1492-1510`
- Test: `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: nothing.
- Produces: `unsigned Y(nfft_map_planning_flags)(unsigned planning)`,
  `unsigned Y(nfft_derive_fftw_flags)(unsigned planning, unsigned fftw_flags)`,
  public `NFFT_PATIENT (1U<<5)`, `NFFT_NO_NONTHREADED (1U<<7)`,
  `NFFT_WISDOM_ONLY (1U<<8)`, internal `PLNR_NO_NONTHREADED = 0x0400`.

- [ ] **Step 1: Public flags**

Replace `include/nfft3.h:845-849` with:

```c
/* Planning flags. Patience rises MEASURE -> PATIENT; ESTIMATE skips
 * measurement entirely. Patience is expressed internally by the absence of
 * restrictions, so a more patient request searches a wider space and takes
 * longer to plan. ESTIMATE overrides PATIENT if both are given. */
#define NFFT_MEASURE         (0U)       /* Measure solutions. */
#define NFFT_ESTIMATE        (1U << 0)  /* Estimate winner. */
#define NFFT_NO_DIRECT       (1U << 1)  /* Do not use direct (slow) algorithms. */
#define NFFT_NO_FAST_NATIVE  (1U << 4)  /* Do not use the fast NFFT algorithm. */
#define NFFT_PATIENT         (1U << 5)  /* Widen the search; let serial and
                                         * threaded solvers compete. */
/* (1U << 6) is reserved for NFFT_EXHAUSTIVE. Do not reuse it. */
#define NFFT_NO_NONTHREADED  (1U << 7)  /* Beyond-guru: forbid serial solvers
                                         * whenever more than one thread was
                                         * requested, whatever the patience. */
#define NFFT_WISDOM_ONLY     (1U << 8)  /* Plan only from existing wisdom; the
                                         * guru returns NULL on a miss. */
```

- [ ] **Step 2: Internal bit**

`include/iplanner.h`, extend the impatience enum after
`PLNR_NO_FAST_NATIVE = 0x0200,`:

```c
  PLNR_NO_NONTHREADED = 0x0400 /* a serial solver may not answer when more than
                                * one thread was requested */
```

- [ ] **Step 3: Declare the mapping functions**

`include/iplanner.h`, after `Y(planner_mkplan)` (line 323):

```c
/* mapflags.c: translation of the caller's request. Both are pure functions of
 * their arguments; no planner or problem state is read. */
unsigned Y(nfft_map_planning_flags)(unsigned planning);
unsigned Y(nfft_derive_fftw_flags)(unsigned planning, unsigned fftw_flags);
```

- [ ] **Step 4: Write the failing tests**

Append to `tests/planner.c`:

```c
/* Both mapping functions are pure, so they are checked directly rather than
 * through a plan. Mirrors FFTW api/mapflags.c: patience is the absence of
 * restrictions, and everything below PATIENT carries PLNR_NO_NONTHREADED. */
void Y(check_planner_mapflags)(void)
{
  unsigned f;

  f = Y(nfft_map_planning_flags)(NFFT_MEASURE); /* the default, below PATIENT */
  CU_ASSERT_TRUE(f & PLNR_NO_NONTHREADED);
  CU_ASSERT_TRUE(f & PLNR_BELIEVE_PCOST);
  CU_ASSERT_TRUE(f & PLNR_NO_UGLY);
  CU_ASSERT_TRUE(f & PLNR_NO_SLOW);
  CU_ASSERT_FALSE(f & PLNR_ESTIMATE);

  f = Y(nfft_map_planning_flags)(NFFT_ESTIMATE);
  CU_ASSERT_TRUE(f & PLNR_NO_NONTHREADED);
  CU_ASSERT_TRUE(f & PLNR_ESTIMATE);
  CU_ASSERT_TRUE(f & PLNR_ALLOW_PRUNING);

  f = Y(nfft_map_planning_flags)(NFFT_PATIENT); /* drops the below-PATIENT block */
  CU_ASSERT_FALSE(f & PLNR_NO_NONTHREADED);
  CU_ASSERT_FALSE(f & PLNR_BELIEVE_PCOST);
  /* NO_UGLY and NO_SLOW stay set: only EXHAUSTIVE clears them, and its flag is
   * reserved rather than exposed. */
  CU_ASSERT_TRUE(f & PLNR_NO_UGLY);
  CU_ASSERT_TRUE(f & PLNR_NO_SLOW);

  f = Y(nfft_map_planning_flags)(NFFT_ESTIMATE | NFFT_PATIENT); /* ESTIMATE wins */
  CU_ASSERT_TRUE(f & PLNR_ESTIMATE);
  CU_ASSERT_TRUE(f & PLNR_NO_NONTHREADED);

  f = Y(nfft_map_planning_flags)(NFFT_PATIENT | NFFT_NO_NONTHREADED);
  CU_ASSERT_TRUE(f & PLNR_NO_NONTHREADED); /* the override survives PATIENT */

  f = Y(nfft_map_planning_flags)(NFFT_PATIENT | NFFT_NO_DIRECT
                                 | NFFT_NO_FAST_NATIVE);
  CU_ASSERT_TRUE(f & PLNR_NO_DIRECT); /* gates are orthogonal to patience */
  CU_ASSERT_TRUE(f & PLNR_NO_FAST_NATIVE);
}

/* iplanner.h documents problem_nfft.fftw_flags as "0 = derive". Zero takes the
 * patience level from the NFFT request; a non-zero word is the caller's own.
 * Wisdom-only is settled from the planning word on both paths. */
void Y(check_planner_derive_fftw_flags)(void)
{
  unsigned f;

  f = Y(nfft_derive_fftw_flags)(NFFT_ESTIMATE, 0u);
  CU_ASSERT_TRUE(f & FFTW_ESTIMATE);
  CU_ASSERT_FALSE(f & FFTW_PATIENT);

  f = Y(nfft_derive_fftw_flags)(NFFT_MEASURE, 0u); /* FFTW_MEASURE is zero */
  CU_ASSERT_FALSE(f & FFTW_ESTIMATE);
  CU_ASSERT_FALSE(f & FFTW_PATIENT);

  f = Y(nfft_derive_fftw_flags)(NFFT_PATIENT, 0u);
  CU_ASSERT_TRUE(f & FFTW_PATIENT);

  f = Y(nfft_derive_fftw_flags)(NFFT_ESTIMATE, FFTW_PATIENT); /* caller wins */
  CU_ASSERT_TRUE(f & FFTW_PATIENT);
  CU_ASSERT_FALSE(f & FFTW_ESTIMATE);

  CU_ASSERT_NOT_EQUAL(Y(nfft_derive_fftw_flags)(NFFT_MEASURE, 0u),
                      Y(nfft_derive_fftw_flags)(NFFT_PATIENT, 0u));
  CU_ASSERT_NOT_EQUAL(Y(nfft_derive_fftw_flags)(NFFT_ESTIMATE, 0u),
                      Y(nfft_derive_fftw_flags)(NFFT_MEASURE, 0u));

  f = Y(nfft_derive_fftw_flags)(NFFT_MEASURE | NFFT_WISDOM_ONLY, 0u);
  CU_ASSERT_TRUE(f & FFTW_WISDOM_ONLY);
  f = Y(nfft_derive_fftw_flags)(NFFT_PATIENT | NFFT_WISDOM_ONLY, FFTW_ESTIMATE);
  CU_ASSERT_TRUE(f & FFTW_WISDOM_ONLY);
  CU_ASSERT_TRUE(f & FFTW_ESTIMATE);
  f = Y(nfft_derive_fftw_flags)(NFFT_MEASURE, FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
  CU_ASSERT_FALSE(f & FFTW_WISDOM_ONLY); /* cleared even when the caller set it */
  CU_ASSERT_TRUE(f & FFTW_ESTIMATE);
}

/* Wisdom-only is a planning directive, not a property of the problem, so it
 * must not enter the key. If it did, a wisdom-only attempt would look under a
 * different key from the plan that wrote the entry and never find it. */
void Y(check_planner_wisdom_only_not_keyed)(void)
{
  planner *pl = Y(planner_create)();
  const INT N = 32, n = 64, M = 10;
  R x[10];
  md5sig a, b;
  problem *p, *q;
  INT j;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);

  p = Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, +1,
                        Y(nfft_derive_fftw_flags)(NFFT_MEASURE, 0u), x, 1, 0, 0);
  q = Y(mkproblem_nfft)(
       1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, +1,
       Y(nfft_derive_fftw_flags)(NFFT_MEASURE | NFFT_WISDOM_ONLY, 0u), x, 1, 0,
       0);
  Y(problem_md5)(pl, p, a);
  Y(problem_md5)(pl, q, b);
  CU_ASSERT_TRUE(a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3]);

  Y(problem_destroy)(q);
  Y(problem_destroy)(p);
  Y(planner_destroy)(pl);
}
```

Declare all three in `tests/planner.h`; register in `tests/check_ng.c`:

```c
  CU_add_test(planner_suite, "mapflags", Y(check_planner_mapflags));
  CU_add_test(planner_suite, "derive_fftw_flags",
              Y(check_planner_derive_fftw_flags));
  CU_add_test(planner_suite, "wisdom_only_not_keyed",
              Y(check_planner_wisdom_only_not_keyed));
```

- [ ] **Step 5: Run to verify failure**

Run: `make -j && tests/checkall_ng`
Expected: link failure, `undefined reference to nfft_nfft_map_planning_flags`.

- [ ] **Step 6: Write the mapping unit**

Create `kernel/nfft/mapflags.c` with the repo's GPL header, then:

```c
/* Translation of the caller's request, after FFTW's api/mapflags.c.
 *
 * Patience is the absence of restriction: an impatient request carries more
 * PLNR_NO_* bits and searches a narrower space. The public word names the
 * level; this file turns a level into the restrictions it implies, and into the
 * planner flags handed to the child FFTW plans. Nothing here reads planner or
 * problem state. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

/* Reserved public bit for a future NFFT_EXHAUSTIVE. The branch below is written
 * against it so that exposing the flag is a one-line change; until then it is
 * unreachable and uncovered, which is a known, accepted cost. */
#define NFFT_EXHAUSTIVE_RESERVED (1U << 6)

/* EXHAUSTIVE implies PATIENT; ESTIMATE denies both. */
static void levels(unsigned planning, unsigned *estimate, unsigned *patient,
                   unsigned *exhaustive)
{
  *estimate = (planning & NFFT_ESTIMATE) ? 1u : 0u;
  *exhaustive = (planning & NFFT_EXHAUSTIVE_RESERVED) ? 1u : 0u;
  *patient = ((planning & NFFT_PATIENT) || *exhaustive) ? 1u : 0u;
  if (*estimate)
    *patient = *exhaustive = 0u;
}

unsigned Y(nfft_map_planning_flags)(unsigned planning)
{
  unsigned estimate, patient, exhaustive;
  unsigned F = 0;

  levels(planning, &estimate, &patient, &exhaustive);

  if (estimate)
    F |= PLNR_ESTIMATE | PLNR_ALLOW_PRUNING;

  if (!patient) /* the fftw2-like block of restrictions */
    F |= PLNR_NO_NONTHREADED | PLNR_BELIEVE_PCOST;

  if (!exhaustive)
    F |= PLNR_NO_UGLY | PLNR_NO_SLOW;

  if (planning & NFFT_NO_NONTHREADED) /* beyond-guru, patience-independent */
    F |= PLNR_NO_NONTHREADED;

  if (planning & NFFT_NO_DIRECT) /* gates are orthogonal to patience */
    F |= PLNR_NO_DIRECT;
  if (planning & NFFT_NO_FAST_NATIVE)
    F |= PLNR_NO_FAST_NATIVE;

  return F;
}

/* The child FFTW plans' flags. Zero means derive from the NFFT patience level,
 * which is what problem_nfft.fftw_flags has always documented. FFTW_PATIENT is
 * the one that does real work: it makes FFTW compare its own threaded and
 * serial candidates instead of preferring the threaded one.
 *
 * FFTW_DESTROY_INPUT is added and FFTW_PRESERVE_INPUT stripped where the child
 * plan is built (kernel/nfft/nfft-nd.c), not here, so those bits stay out of
 * the wisdom key. */
unsigned Y(nfft_derive_fftw_flags)(unsigned planning, unsigned fftw_flags)
{
  unsigned estimate, patient, exhaustive;
  unsigned ff;

  if (fftw_flags != 0u)
    ff = fftw_flags;
  else {
    levels(planning, &estimate, &patient, &exhaustive);
    if (estimate)
      ff = (unsigned)FFTW_ESTIMATE;
    else if (exhaustive)
      ff = (unsigned)FFTW_EXHAUSTIVE;
    else if (patient)
      ff = (unsigned)FFTW_PATIENT;
    else
      ff = (unsigned)FFTW_MEASURE; /* zero */
  }

  /* One source of truth for wisdom-only, on both paths. */
  if (planning & NFFT_WISDOM_ONLY)
    ff |= (unsigned)FFTW_WISDOM_ONLY;
  else
    ff &= ~(unsigned)FFTW_WISDOM_ONLY;

  return ff;
}
```

- [ ] **Step 7: Add to both builds**

`kernel/nfft/Makefile.am`: add `mapflags.c` alphabetically between `conf.c` and
`ndft-1d.c`, in every source list the file defines.
`CMakeLists.txt`: add `kernel/nfft/mapflags.c` beside `kernel/nfft/conf.c`.

- [ ] **Step 8: Wire into `plan.c` and keep wisdom-only out of the key**

Replace `keyable_fftw_flags` (`plan.c:34-42`) with:

```c
/* Strip the bits that are planning directives rather than properties of the
 * problem, before fftw_flags reaches the wisdom key. The preservation bits go
 * because no planner-native candidate mutates its input in place, so the two
 * spellings must not key distinct entries. FFTW_WISDOM_ONLY goes because it
 * says how hard to look for a plan, not which plan is wanted: a wisdom-only
 * attempt must find the entry an ordinary plan wrote. */
static unsigned keyable_fftw_flags(unsigned fftw_flags)
{
  return fftw_flags
         & ~(unsigned)(FFTW_DESTROY_INPUT | FFTW_PRESERVE_INPUT
                       | FFTW_WISDOM_ONLY);
}
```

Delete `map_planning_flags` (lines 44-53) and change line 80 to
`F = Y(nfft_map_planning_flags)(planning);`.

Change the problem construction at lines 105-107 to pass the derived word:

```c
  /* The derived word is what the child FFTW plans use and what the wisdom key
   * records, so two patience levels never share a key while planning different
   * child FFTs. */
  p->prob[FWD] = Y(mkproblem_nfft)(
       d, N, variant, n, M, m, window, +1,
       keyable_fftw_flags(Y(nfft_derive_fftw_flags)(planning, fftw_flags)), x,
       /*copy_x=*/1, (C *)f_hat, (C *)f);
```

Grep for every other `Y(mkproblem_nfft)` in the file and apply the same change.

- [ ] **Step 9: Migrate the existing wisdom-only test**

`tests/nplan.c:1492-1510` asks for wisdom-only through `fftw_flags`, which
Step 6 now clears. Replace the guru call in
`Y(check_nplan_fftw_wisdom_only_declines)` with:

```c
  /* Wisdom-only now comes from the planning word; FFTW_WISDOM_ONLY inside
   * fftw_flags is cleared by Y(nfft_derive_fftw_flags). */
  CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(
       1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u,
       NFFT_ESTIMATE | NFFT_NO_DIRECT | NFFT_WISDOM_ONLY));

  /* The old spelling is inert: the child plans are built normally, so the guru
   * succeeds. */
  {
    Y(plan_ng) *p = Y(plan_ng_guru)(
         1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f,
         FFTW_WISDOM_ONLY | FFTW_ESTIMATE, NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p)
      Y(plan_ng_destroy)(p);
  }
```

Update the comment above the function; keep its registered name so the CUnit
history stays continuous. **The first assertion only passes once Task 6 lands.**
Until then this case is deliberately red; that is the plan's only cross-task red
state and Task 6 step 6 closes it.

- [ ] **Step 10: Run**

Run: `./bootstrap.sh && ./configure --enable-all --enable-tests && make -j && make check`
Expected: `planner/mapflags`, `planner/derive_fftw_flags` and
`planner/wisdom_only_not_keyed` pass. `nplan/fftw_wisdom_only_declines` fails on
its first assertion, as above.

- [ ] **Step 11: Format and commit**

```bash
clang-format -i kernel/nfft/mapflags.c kernel/nfft/plan.c include/nfft3.h include/iplanner.h tests/planner.c tests/nplan.c
git add include/nfft3.h include/iplanner.h kernel/nfft/mapflags.c kernel/nfft/plan.c kernel/nfft/Makefile.am CMakeLists.txt tests/planner.c tests/planner.h tests/nplan.c tests/check_ng.c
git commit -m "Map the public planning word through an FFTW-style patience lattice and make it the only source of the child FFTW flags."
```

**Not tested, by decision:** the `| FFTW_DESTROY_INPUT`, the
`& ~FFTW_PRESERVE_INPUT` and the out-of-place `g1 != g2` at
`kernel/nfft/nfft-nd.c:209-217`. None is observable from a test — FFTW's plan
description carries no placement information — so they are left to code review.
Do not change that block.

---

### Task 2: The decline rule

Implements R4.

**Files:**
- Modify: `include/iplanner.h:305-306`, `kernel/nfft/nfft-nd.c`,
  `kernel/nfft/ndft-1d.c`, `kernel/nfft/ndft-nd.c`,
  `kernel/deconv/deconv-{1d,2d,3d,nd}.c`, `kernel/conv/conv-{1d,2d,3d,nd}.c`
- Test: `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: `PLNR_NO_NONTHREADED`, `Y(nfft_map_planning_flags)` from Task 1.
- Produces: `NO_NONTHREADEDP(pl)`.

- [ ] **Step 1: Write the failing test**

`checkall_ng` links `libnfft3` alone, so `X(plan_with_nthreads)` is unavailable
to it. The test sets the count on a private planner directly. Append to
`tests/planner.c`:

```c
/* A serial solver steps aside when the caller asked for more than one thread
 * and the patience level is below PATIENT. Mirrors FFTW kernel/ifftw.h
 * NO_NONTHREADEDP and the declines at dft/ct.c:132, dft/vrank-geq1.c:139. */
void Y(check_planner_nonthreaded_decline)(void)
{
  planner *pl = Y(planner_create)();
  const INT N = 32, n = 64;
  problem *p = Y(mkproblem_deconv)(1, &N, 0, &n, 6,
                                   NFFT_WINDOW_KAISER_BESSEL, +1, 0, 0);
  plan *pln;

  Y(deconv_solvers_register)(pl);

  pl->nthr = 1; /* one thread: the rule never fires */
  pl->flags.l = pl->flags.u = Y(nfft_map_planning_flags)(NFFT_MEASURE);
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NOT_NULL(pln);
  if (pln)
    Y(plan_destroy)(pln);

  pl->nthr = 4; /* more than one, below PATIENT: the serial solver declines */
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NULL(pln);

  pl->flags.l = pl->flags.u = Y(nfft_map_planning_flags)(NFFT_PATIENT);
  pln = Y(planner_mkplan)(pl, p); /* back in the running */
  CU_ASSERT_PTR_NOT_NULL(pln);
  if (pln)
    Y(plan_destroy)(pln);

  pl->flags.l = pl->flags.u =
       Y(nfft_map_planning_flags)(NFFT_PATIENT | NFFT_NO_NONTHREADED);
  pln = Y(planner_mkplan)(pl, p); /* the override reinstates the decline */
  CU_ASSERT_PTR_NULL(pln);

  Y(problem_destroy)(p);
  Y(planner_destroy)(pl);
}
```

Declare in `tests/planner.h`; register as
`CU_add_test(planner_suite, "nonthreaded_decline", Y(check_planner_nonthreaded_decline));`.

- [ ] **Step 2: Run to verify failure**

Run: `make -j && tests/checkall_ng`
Expected: FAIL at the first `CU_ASSERT_PTR_NULL` — nothing declines yet.

- [ ] **Step 3: Add the predicate**

`include/iplanner.h`, next to `PLNR_L` and `PLNR_U` (lines 305-306):

```c
/* A serial solver must step aside when the caller asked for more than one
 * thread and the patience level is below PATIENT, so PLNR_NO_NONTHREADED is
 * set. FFTW's kernel/ifftw.h has the same two conditions. It is safe without a
 * third because raising the count requires the add-on library, and the add-on
 * library is where threaded solvers are registered. */
#define NO_NONTHREADEDP(pl)                                                   \
  ((PLNR_L(pl) & PLNR_NO_NONTHREADED) && (pl)->nthr > 1)
```

- [ ] **Step 4: Eleven decline sites**

In each `mkplan`, immediately after the existing `problem_kind` check and before
any allocation:

```c
  if (NO_NONTHREADEDP(pl))
    return 0; /* serial: prefer a threaded solver */
```

Sites: `kernel/nfft/nfft-nd.c` (`mkplan_native_fast`), `kernel/nfft/ndft-1d.c`,
`kernel/nfft/ndft-nd.c`, `kernel/deconv/deconv-{1d,2d,3d,nd}.c`,
`kernel/conv/conv-{1d,2d,3d,nd}.c`. Delete the `(void)pl;` line wherever the
parameter becomes used.

Do **not** touch `kernel/nfft/rnk0.c`: the rank-0 base case has nothing to
parallelise and must stay the terminal solver for a fully elided problem.

- [ ] **Step 5: Run**

Run: `make -j && make check`
Expected: PASS. Every existing case still passes because nothing linked against
`libnfft3` alone can raise `nthr` above 1.

- [ ] **Step 6: Format and commit**

```bash
clang-format -i include/iplanner.h kernel/nfft/nfft-nd.c kernel/nfft/ndft-1d.c kernel/nfft/ndft-nd.c kernel/deconv/*.c kernel/conv/*.c tests/planner.c
git add include/iplanner.h kernel/nfft kernel/deconv kernel/conv tests/planner.c tests/planner.h tests/check_ng.c
git commit -m "Let serial solvers step aside for threaded ones below PATIENT."
```

---

### Task 3: Dynamic scoping of the thread budget

Implements R5.

**Files:** Modify `kernel/planner/planner.c:547-620`. Test `tests/planner.c`,
`tests/planner.h`, `tests/check_ng.c`.

**Interfaces:**
- Consumes: nothing.
- Produces: the guarantee that `pl->nthr`, `pl->flags.l` and `pl->flags.u` are
  unchanged across any `mkplan` the search performs. Task 5 extends the same
  helper to a fourth field, so keep the save-and-restore in one place.

- [ ] **Step 1: Write the failing test**

Append to `tests/planner.c`:

```c
/* A solver that lowers the budget for its children, the way a threaded solver
 * will. The planner must restore what the solver changed, so the next candidate
 * is planned under the caller's budget. Mirrors FFTW's invoke_solver. */
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
                                           mkplan_scope_probe};

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
  CU_ASSERT_PTR_NULL(pln);
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

- [ ] **Step 2: Run to verify failure**

Run: `make -j && tests/checkall_ng`
Expected: FAIL at `CU_ASSERT_EQUAL(pl->nthr, 6)`.

- [ ] **Step 3: Add the helper**

`kernel/planner/planner.c`, above `Y(planner_mkplan)` (line 547):

```c
/* Every mkplan the search performs runs inside the caller's budget and may not
 * change it for the next candidate. A threaded solver lowers pl->nthr for its
 * children; this restores it. FFTW does the same in invoke_solver. */
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

- [ ] **Step 4: Route the search through it**

Replace the three direct calls with `invoke_solver(pl, p, s)`: line 568 (the
wisdom-hit re-plan), line 577 (`Y(planner_mkplan)`'s enumeration), line 610
(`Y(planner_candidates)`'s enumeration). Grep for `adt->mkplan(` and confirm the
only remaining occurrence is inside `invoke_solver`. Leave `plan.c`'s own
wisdom-hit re-plan alone; it is outside the search.

- [ ] **Step 5: Run and commit**

Run: `make -j && make check` — PASS including `planner/nthr_scoping`.

```bash
clang-format -i kernel/planner/planner.c tests/planner.c
git add kernel/planner/planner.c tests/planner.c tests/planner.h tests/check_ng.c
git commit -m "Scope the thread budget and impatience bounds across every solver call in the search."
```

---

### Task 4: FFTW's thread count in the problem key

Implements R11.

**Files:** Modify `configure.ac`, `include/iplanner.h:312`,
`kernel/planner/planner.c`, `kernel/nfft/problem.c:71`. Test `tests/planner.c`,
`tests/planner.h`, `tests/check_ng.c`.

**Interfaces:**
- Consumes: nothing.
- Produces: `extern int (*Y(fftw_nthreads_hook))(void);`, null by default.
  Task 7's `X(init_threads)` installs it.

- [ ] **Step 1: Write the failing test**

Append to `tests/planner.c`:

```c
static int fake_fftw_nthreads_value = 1;
static int fake_fftw_nthreads(void) { return fake_fftw_nthreads_value; }

/* FFTW's own thread count selects a different child FFT plan for the same
 * flags, so it must change our key. libnfft3 links @fftw3_LIBS@ only and
 * fftw_planner_nthreads lives in FFTW's threads library, so the value arrives
 * through a hook the add-on installs. A null hook means one thread. */
void Y(check_planner_fftw_nthreads_key)(void)
{
  planner *pl = Y(planner_create)();
  const INT N = 32, n = 64, M = 10;
  R x[10];
  md5sig a, b, c;
  problem *p;
  INT j;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);

  CU_ASSERT_PTR_NULL(Y(fftw_nthreads_hook));

  p = Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u,
                        x, 1, 0, 0);
  Y(problem_md5)(pl, p, a);

  Y(fftw_nthreads_hook) = fake_fftw_nthreads;
  fake_fftw_nthreads_value = 1; /* a hook reporting 1 must not move the key */
  Y(problem_md5)(pl, p, b);
  CU_ASSERT_TRUE(a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3]);

  fake_fftw_nthreads_value = 8; /* a different count must */
  Y(problem_md5)(pl, p, c);
  CU_ASSERT_FALSE(a[0] == c[0] && a[1] == c[1] && a[2] == c[2] && a[3] == c[3]);

  Y(fftw_nthreads_hook) = 0;
  Y(problem_destroy)(p);
  Y(planner_destroy)(pl);
}
```

Declare in `tests/planner.h`; register as
`CU_add_test(planner_suite, "fftw_nthreads_key", Y(check_planner_fftw_nthreads_key));`.

- [ ] **Step 2: Run to verify failure**

Run: `make -j && tests/checkall_ng`
Expected: link failure, `undefined reference to nfft_fftw_nthreads_hook`.

- [ ] **Step 3: Probe FFTW**

`configure.ac`, after the FFTW3 checks that set `fftw3_LIBS_omp` (around line
480). Use the precision-mangled spelling; `PREC_SUFFIX` is set at lines 138-142
and `m4/nfft_lib_fftw3.m4:144` already probes `fftw3${PREC_SUFFIX}_omp` the same
way.

```m4
# fftw_planner_nthreads arrived in FFTW 3.3.9. The add-on threading library
# reports FFTW's thread count into our wisdom key through it; without it that
# library is not built.
AC_CHECK_DECL([fftw${PREC_SUFFIX}_planner_nthreads],
  [AC_DEFINE([HAVE_FFTW_PLANNER_NTHREADS],[1],
     [Define to 1 if FFTW declares fftw_planner_nthreads.])
   nfft_have_fftw_planner_nthreads=yes],
  [nfft_have_fftw_planner_nthreads=no],
  [[#include <fftw3.h>]])
```

- [ ] **Step 4: Declare and define the hook**

`include/iplanner.h`, near line 312:

```c
/* FFTW's own planner thread count, reported by the add-on threading library.
 * libnfft3 links @fftw3_LIBS@ only and fftw_planner_nthreads is defined in
 * FFTW's threads library, so the main library cannot call it and takes the
 * value through this hook. Null means one thread. Same device as FFTW's
 * threads/api.c mksolver_ct_hook. */
extern int (*Y(fftw_nthreads_hook))(void);
```

`kernel/planner/planner.c`, at file scope:

```c
int (*Y(fftw_nthreads_hook))(void) = 0;
```

- [ ] **Step 5: Hash it**

`kernel/nfft/problem.c`, in `hash()`, after
`Y(md5_put_unsigned)(ctx, ego->fftw_flags);` (line 71):

```c
  /* FFTW's thread count selects a different child FFT plan for the same flags,
   * so it belongs in the key. Observed, never set. */
  Y(md5_put_int)(ctx, Y(fftw_nthreads_hook) ? Y(fftw_nthreads_hook)() : 1);
```

- [ ] **Step 6: Run and commit**

Run: `./bootstrap.sh && ./configure --enable-all --enable-tests && make -j && make check`
Expected: PASS. Wisdom written before this change is now a clean miss, which
Task 9's signature tag makes explicit.

```bash
clang-format -i include/iplanner.h kernel/planner/planner.c kernel/nfft/problem.c tests/planner.c
git add configure.ac include/iplanner.h kernel/planner/planner.c kernel/nfft/problem.c tests/planner.c tests/planner.h tests/check_ng.c
git commit -m "Take FFTW's planner thread count into the problem key through a hook."
```

---

### Task 5: Bless the whole winning tree

Implements R12, R13, R14. Prerequisite for Task 6: a recursive wisdom-only gate
only works if children are in wisdom.

**Files:** Modify `kernel/planner/planner.c` (`search_flags`, `invoke_solver`),
`kernel/nfft/plan.c:120-130`, `:342-355`, `:375-380`. Test `tests/nplan.c`,
`tests/nplan.h`, `tests/check_ng.c`.

**Interfaces:**
- Consumes: `invoke_solver` from Task 3.
- Produces: the invariant that after a successful guru call every node of the
  winning plan tree has a blessed wisdom entry, and no node of a losing
  candidate does. Also `int count_wisdom_entries(void)` in `tests/nplan.c`,
  declared in `tests/nplan.h` because Task 6 reuses it.

- [ ] **Step 1: Write the failing test**

Append to `tests/nplan.c`:

```c
/* Entry lines in an exported wisdom file are indented; the preamble and the
 * closing paren are not. Declared in nplan.h because the wisdom-only case
 * reuses it. */
int Y(count_wisdom_entries)(void)
{
  char *w = NFFT(export_wisdom_to_string)();
  int n = 0;
  if (w) {
    const char *s = w;
    while ((s = strchr(s, '\n')) != 0) {
      s++;
      if (*s == ' ')
        n++;
    }
    NFFT(free)(w);
  }
  return n;
}

/* After planning, an export must describe the whole winning tree: the NFFT
 * solution and its DECONV and CONV children. FFTW gets this by re-creating the
 * finished plan with BLESSING set so every lookup hits the first pass's memo
 * (api/apiplan.c:147); we do the same. Losing candidates' children stay out. */
void Y(check_nplan_blesses_whole_tree)(void)
{
  const INT N = 64, n = 128, M = 100;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  INT j;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);
  for (j = 0; j < N; j++)
    f_hat[j] = K(0.0);

  Y(the_planner_destroy)();
  NFFT(forget_wisdom)();
  {
    Y(plan_ng) *p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6,
                                       NFFT(get_window_id)(), x, (FC *)f_hat,
                                       (FC *)f, 0u, NFFT_MEASURE);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p)
      NFFT(plan_ng_destroy)(p);
  }
  CU_ASSERT_TRUE(Y(count_wisdom_entries)() >= 3); /* NFFT + DECONV + CONV */

  /* Estimate solutions are blessed too, FFTW parity: mkplan0 sets BLESSING
   * whatever the patience and exprt filters nothing. */
  Y(the_planner_destroy)();
  NFFT(forget_wisdom)();
  {
    Y(plan_ng) *p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6,
                                       NFFT(get_window_id)(), x, (FC *)f_hat,
                                       (FC *)f, 0u, NFFT_ESTIMATE);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p)
      NFFT(plan_ng_destroy)(p);
  }
  CU_ASSERT_TRUE(Y(count_wisdom_entries)() >= 3);

  Y(the_planner_destroy)();
  NFFT(forget_wisdom)();
  Y(free)(f);
  Y(free)(f_hat);
  Y(free)(x);
}
```

Declare `Y(count_wisdom_entries)` and `Y(check_nplan_blesses_whole_tree)` in
`tests/nplan.h`; register as
`CU_add_test(nplan_suite, "blesses_whole_tree", Y(check_nplan_blesses_whole_tree));`.

- [ ] **Step 2: Run to verify failure**

Run: `make -j && tests/checkall_ng`
Expected: FAIL — both counts are 1; only the top-level entry is blessed.

- [ ] **Step 3: Carry the blessing bit into the search's inserts**

`kernel/planner/planner.c`, in `search_flags`, replace `q.info = 0;` with:

```c
  /* Blessing is dynamically scoped: whatever the planner is currently blessing,
   * every solution memoised beneath it inherits. That is how FFTW gets the
   * whole winning tree blessed on its second pass without enumerating children
   * (api/apiplan.c:147, kernel/planner.c BLISS). */
  q.info = pl->flags.info & PLNR_BLESSING;
```

Extend `invoke_solver` from Task 3 to save and restore `pl->flags.info`, so a
solver cannot leak a blessing change sideways:

```c
static plan *invoke_solver(planner *pl, const problem *p, solver *s)
{
  int nthr = pl->nthr;
  unsigned l = pl->flags.l, u = pl->flags.u, info = pl->flags.info;
  plan *pln = s->adt->mkplan(s, p, pl);
  pl->nthr = nthr;
  pl->flags.l = l;
  pl->flags.u = u;
  pl->flags.info = info;
  return pln;
}
```

- [ ] **Step 4: Bless the tree in the measured path**

`kernel/nfft/plan.c`, the successful-race branch — currently the `else` of
`if (timed_out)` at line 354, holding
`Y(planner_bless)(pl, p->prob[dirn], cslvndx[dirn][winner_idx]);`. Replace that
one statement with:

```c
      else {
        Y(planner_bless)(pl, p->prob[dirn], cslvndx[dirn][winner_idx]);
        /* Second pass, FFTW's trick: re-plan the winner with the blessing bit
         * set. Every node hits the memo the race left behind, so the whole
         * winning tree is blessed and no losing candidate's children are. The
         * plan it returns is discarded; the raced winner is already in
         * p->dir[dirn] and replacing it would disturb the race's result. */
        pl->flags.info |= PLNR_BLESSING;
        {
          plan *bless_pass = Y(planner_mkplan)(pl, p->prob[dirn]);
          if (bless_pass)
            Y(plan_destroy)(bless_pass);
        }
        pl->flags.info &= ~(unsigned)PLNR_BLESSING;
      }
```

The timed-out branch is unchanged: a partial race stays unblessed, so its
children must stay unblessed too.

- [ ] **Step 5: Bless in the estimate path**

`kernel/nfft/plan.c:120-130`, the estimate branch. There is no race, so one pass
suffices. Set the bit around `select_estimate` and clear it on every exit,
keeping every existing cleanup statement:

```c
  if (is_estimate) {
    /* Bounds {l = F, u = PLNR_ESTIMATE | F}. Estimate solutions are blessed and
     * exported, matching FFTW: mkplan0 sets BLESSING whatever the patience and
     * exprt filters nothing. They can only ever answer another estimate query,
     * because PLNR_ESTIMATE rides in u and LEQ(a->u, q->u) then fails for a
     * measured query -- which is FFTW's behaviour too. */
    pl->flags.l = F;
    pl->flags.u = PLNR_ESTIMATE | F;
    pl->flags.info |= PLNR_BLESSING;

    if (!select_estimate(p, pl)) {
      pl->flags.info &= ~(unsigned)PLNR_BLESSING;
      /* ... existing cleanup and return 0, unchanged ... */
    }

    pl->flags.info &= ~(unsigned)PLNR_BLESSING;
    /* ... existing bound restore and return p, unchanged ... */
  }
```

Apply the same set-and-clear around the estimate-restart path near line 377.

- [ ] **Step 6: Run**

Run: `make -j && make check`
Expected: PASS including `nplan/blesses_whole_tree`. Watch
`planner/wisdom_roundtrip`, `planner/subsumption` and `planner/forget`: exported
files now carry more entries, so any test asserting an exact count needs its
expectation updated, not its assertion removed.

- [ ] **Step 7: Format and commit**

```bash
clang-format -i kernel/planner/planner.c kernel/nfft/plan.c tests/nplan.c
git add kernel/planner/planner.c kernel/nfft/plan.c tests/nplan.c tests/nplan.h tests/check_ng.c
git commit -m "Bless the whole winning plan tree with a second pass, as FFTW does."
```

---

### Task 6: Wisdom-only as a planner state

Implements R16, R17.

**Files:** Modify `include/iplanner.h` (enum + `planner_s`),
`kernel/planner/planner.c` (`Y(planner_mkplan)`, constructor),
`kernel/nfft/plan.c` (guru entry, measured path, estimate path). Test
`tests/nplan.c`, `tests/nplan.h`, `tests/check_ng.c`.

**Interfaces:**
- Consumes: tree blessing from Task 5, `NFFT_WISDOM_ONLY` from Task 1,
  `Y(count_wisdom_entries)` from Task 5.
- Produces: `wisdom_state_t` and `planner_s.wisdom_state`.

- [ ] **Step 1: Write the failing test**

Append to `tests/nplan.c`:

```c
/* Wisdom-only plans from the store or not at all, at every patience level and
 * for every kind, children included. It writes nothing back, and a miss must
 * never disturb the store -- FFTW gives it its own branch that bypasses the
 * forget-everything recovery for exactly that reason (api/apiplan.c:102-107). */
void Y(check_nplan_wisdom_only)(void)
{
  const INT N = 64, n = 128, M = 100;
  const INT N2 = 32, n2 = 64;
  const unsigned levels[3] = {NFFT_ESTIMATE, NFFT_MEASURE, NFFT_PATIENT};
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  INT j;
  int i;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);
  for (j = 0; j < N; j++)
    f_hat[j] = K(0.0);

  for (i = 0; i < 3; i++) {
    Y(plan_ng) *p;

    Y(the_planner_destroy)();
    NFFT(forget_wisdom)();

    /* Empty store: refused at every level, and the refusal writes nothing. */
    p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT(get_window_id)(), x,
                           (FC *)f_hat, (FC *)f, 0u,
                           levels[i] | NFFT_WISDOM_ONLY);
    CU_ASSERT_PTR_NULL(p);
    CU_ASSERT_EQUAL(Y(count_wisdom_entries)(), 0);

    /* Plan normally to fill the store, then wisdom-only must succeed. */
    p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT(get_window_id)(), x,
                           (FC *)f_hat, (FC *)f, 0u, levels[i]);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p)
      NFFT(plan_ng_destroy)(p);

    p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT(get_window_id)(), x,
                           (FC *)f_hat, (FC *)f, 0u,
                           levels[i] | NFFT_WISDOM_ONLY);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p) {
      NFFT(precompute)(p);
      NFFT(execute)(p);
      NFFT(plan_ng_destroy)(p);
    }

    /* A failed wisdom-only call must not poison the next ordinary call: the
     * state is reset at the guru boundary, as FFTW's mkplan0 does. */
    p = NFFT(plan_ng_guru)(1, &N2, 0, &n2, M, 6, NFFT(get_window_id)(), x,
                           (FC *)f_hat, (FC *)f, 0u,
                           levels[i] | NFFT_WISDOM_ONLY);
    CU_ASSERT_PTR_NULL(p); /* a different size class, not in the store */
    p = NFFT(plan_ng_guru)(1, &N2, 0, &n2, M, 6, NFFT(get_window_id)(), x,
                           (FC *)f_hat, (FC *)f, 0u, levels[i]);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p)
      NFFT(plan_ng_destroy)(p);
  }

  Y(the_planner_destroy)();
  NFFT(forget_wisdom)();
  Y(free)(f);
  Y(free)(f_hat);
  Y(free)(x);
}
```

Declare in `tests/nplan.h`; register as
`CU_add_test(nplan_suite, "wisdom_only", Y(check_nplan_wisdom_only));`.

- [ ] **Step 2: Run to verify failure**

Run: `make -j && tests/checkall_ng`
Expected: FAIL on the first `CU_ASSERT_PTR_NULL` — nothing gates the search yet.
`nplan/fftw_wisdom_only_declines` from Task 1 step 9 is still red here too.

- [ ] **Step 3: Add the state**

`include/iplanner.h`, beside the awake-state enum:

```c
/* Whether the planner may search. Deliberately outside flags_t: wisdom-only is
 * a directive about how hard to look, not a property of a stored solution, so
 * it must not take part in LEQ subsumption or key an entry. FFTW keeps it
 * outside its flags word for the same reason (kernel/ifftw.h wisdom_state_t). */
typedef enum {
  PLNR_WISDOM_NORMAL = 0, /* search freely */
  PLNR_WISDOM_ONLY = 1, /* answer from the store or fail */
  PLNR_WISDOM_IS_BOGUS = 2 /* a wisdom-only lookup missed; unwind */
} wisdom_state_t;
```

Add `wisdom_state_t wisdom_state;` to `struct planner_s` beside `flags`, and
initialise it to `PLNR_WISDOM_NORMAL` in the constructor beside
`pl->nthr = 1;` (`planner.c:447`).

- [ ] **Step 4: Gate the search**

`kernel/planner/planner.c`, in `Y(planner_mkplan)`, immediately after the
`hlookup` block and before the `FORALL_SOLVERS_OF_KIND` enumeration:

```c
  /* Wisdom-only: the store did not answer, so there is nothing to do. Setting
   * the bogus state makes the failure sticky for the rest of this planning
   * call, so a parent solver cannot quietly substitute a different child and
   * hand back a plan that was never in wisdom. FFTW unwinds the same way
   * (kernel/planner.c do_search). */
  if (pl->wisdom_state == PLNR_WISDOM_ONLY) {
    pl->wisdom_state = PLNR_WISDOM_IS_BOGUS;
    return 0;
  }
  if (pl->wisdom_state == PLNR_WISDOM_IS_BOGUS)
    return 0;
```

- [ ] **Step 5: Drive it from the guru**

`kernel/nfft/plan.c`, in `Y(plan_ng_guru)`, after `F` is computed (line 80):

```c
  /* Set per call and cleared on every exit, so a refusal never affects the next
   * call. FFTW's mkplan0 assigns the state on each top-level planning call for
   * the same reason. A wisdom-only plan is never blessed and writes nothing
   * back, so this path must not set PLNR_BLESSING. */
  pl->wisdom_state = (planning & NFFT_WISDOM_ONLY) ? PLNR_WISDOM_ONLY
                                                   : PLNR_WISDOM_NORMAL;
```

Then, on **every** return path of the guru — argument-validation failures, the
fast-guard refusal, the estimate branch, the measured branch, the restart path —
add `pl->wisdom_state = PLNR_WISDOM_NORMAL;` beside the existing
`pl->flags.l = saved_l; pl->flags.u = saved_u;` restores. Grep the function for
`saved_l` and add the line at each site; the early returns that precede the
`saved_l` capture need it too.

In the measured branch, skip the race and the blessing entirely when
`pl->wisdom_state != PLNR_WISDOM_NORMAL`: the lookup either answered or the
direction is absent. In the estimate branch, do not set `PLNR_BLESSING` when
wisdom-only is active.

- [ ] **Step 6: Run**

Run: `make -j && make check`
Expected: PASS including `nplan/wisdom_only`, and
`nplan/fftw_wisdom_only_declines` now green, closing Task 1's red state.

- [ ] **Step 7: Format and commit**

```bash
clang-format -i include/iplanner.h kernel/planner/planner.c kernel/nfft/plan.c tests/nplan.c
git add include/iplanner.h kernel/planner/planner.c kernel/nfft/plan.c tests/nplan.c tests/nplan.h tests/check_ng.c
git commit -m "Gate the search on a wisdom-only planner state that unwinds through children."
```

---

### Task 7: The add-on threading library

Implements R7, R8, R9.

**Files:** Create `kernel/threads/{api.c,conf.c,Makefile.am}`. Modify
`include/nfft3.h`, `include/iplanner.h`, `configure.ac`, `Makefile.am`,
`kernel/Makefile.am`. Tested by Task 8.

**Interfaces:**
- Consumes: `Y(fftw_nthreads_hook)` from Task 4.
- Produces: `int X(init_threads)(void)`, `void X(cleanup_threads)(void)`,
  `void X(plan_with_nthreads)(int)`, `int X(planner_nthreads)(void)`,
  `void Y(nfft_threads_conf_standard)(planner *)`,
  `int Y(nfft_threads_roster_size)(void)`.

- [ ] **Step 1: Declare the entry points**

`include/nfft3.h`, inside `NFFT_DEFINE_PLANNER_API(X,R,C)` after the
`X(set_timelimit)` line (910) — every line carries a trailing backslash except
the macro's last:

```c
NFFT_EXTERN int X(init_threads)(void); \
NFFT_EXTERN void X(cleanup_threads)(void); \
NFFT_EXTERN void X(plan_with_nthreads)(int nthreads); \
NFFT_EXTERN int X(planner_nthreads)(void); \
```

Above the macro, beside the thread-safety note at lines 838-843:

```c
/* Threading. These four are declared here but DEFINED ONLY in the add-on
 * library libnfft3<suffix>_ng_omp, linked in addition to libnfft3<suffix>:
 * `-lnfft3_ng_omp -lnfft3`. FFTW packages its threading the same way. A program
 * that links only libnfft3<suffix> cannot raise the thread count above 1, which
 * is exactly why a serial plan is always available there.
 *
 * X(plan_with_nthreads) initialises threading, which destroys the planner and
 * its in-memory wisdom because the solver roster is about to change. Call it
 * before any other NFFT routine, as FFTW documents for fftw_init_threads.
 *
 * Note libnfft3<suffix>_omp is a different thing: a whole-library OpenMP
 * rebuild for the legacy API, linked INSTEAD of libnfft3<suffix>. */
```

- [ ] **Step 2: Declare the roster**

`include/iplanner.h`, next to `Y(nfft_ensure_registered)` (line 449):

```c
/* kernel/threads/conf.c, in the add-on library only: the threaded solver
 * roster, and whether it registered anything. Empty until the first threaded
 * solver exists. */
void Y(nfft_threads_conf_standard)(planner *pl);
int Y(nfft_threads_roster_size)(void);
```

- [ ] **Step 3: Write the roster**

Create `kernel/threads/conf.c` with the GPL header, then:

```c
/* The threaded solver roster. This is the one place a threaded solver is
 * registered, and it exists only in the add-on library, so a program that does
 * not link that library has no threaded solvers and cannot raise the planner's
 * thread count. FFTW's threads/conf.c has the same job.
 *
 * The table is empty: no threaded solver has been written yet.
 * Y(init_threads) reports that by returning 0, so nothing can reach a state
 * with more than one thread and nothing to run on them. Adding the first entry
 * here is all that is needed to switch the machinery on. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

static const solvtab s = {SOLVTAB_END};

void Y(nfft_threads_conf_standard)(planner *pl) { Y(solvtab_exec)(s, pl); }

int Y(nfft_threads_roster_size)(void)
{
  int n = 0;
  while (s[n].reg != 0)
    n++;
  return n;
}
```

- [ ] **Step 4: Write the API**

Create `kernel/threads/api.c` with the GPL header, then:

```c
/* The add-on threading library's entry points, after FFTW's threads/api.c.
 *
 * Asking for threads is what installs them: X(plan_with_nthreads) initialises
 * threading first, and initialising registers the threaded roster. FFTW relies
 * on the same order, which is why its NO_NONTHREADEDP needs no third condition.
 * Registering solvers changes the configuration signature and so invalidates
 * every stored decision, so initialisation destroys the planner first, as
 * FFTW's X(cleanup)() does. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

#include <fftw3.h>

static int threads_inited = 0;

int Y(init_threads)(void)
{
  if (threads_inited)
    return 1;

  /* Nothing to register yet: report failure rather than leave a caller with a
   * raised thread count and no threaded solver to serve it. */
  if (Y(nfft_threads_roster_size)() == 0)
    return 0;

  Y(the_planner_destroy)(); /* the roster is about to change */
  Y(nfft_ensure_registered)();
  Y(nfft_threads_conf_standard)(Y(the_planner)());

  /* Let the main library observe FFTW's thread count for the wisdom key. */
  Y(fftw_nthreads_hook) = FFTW(planner_nthreads);

  threads_inited = 1;
  return 1;
}

void Y(cleanup_threads)(void)
{
  if (!threads_inited)
    return;
  Y(fftw_nthreads_hook) = 0;
  Y(the_planner_destroy)();
  threads_inited = 0;
}

/* The maximum number of threads a plan may use, not a target: a threaded solver
 * may use fewer and passes the rest of the budget to its children. The count is
 * part of the wisdom key, so changing it makes existing entries a clean miss.
 * Leaves the count at 1 while no threaded solver exists. */
void Y(plan_with_nthreads)(int nthreads)
{
  planner *pl;
  if (!Y(init_threads)())
    return;
  pl = Y(the_planner)();
  pl->nthr = nthreads < 1 ? 1 : nthreads;
}

int Y(planner_nthreads)(void)
{
  Y(nfft_ensure_registered)();
  return Y(the_planner)()->nthr;
}
```

- [ ] **Step 5: Build it**

Create `kernel/threads/Makefile.am`:

```make
AM_CPPFLAGS = -I$(top_srcdir)/include @fftw3_CPPFLAGS@

noinst_LTLIBRARIES = libthreads_ng.la

libthreads_ng_la_SOURCES = api.c conf.c
libthreads_ng_la_CFLAGS = $(OPENMP_CFLAGS)
```

`kernel/Makefile.am`: add `threads` to `SUBDIRS` under an `ENABLE_NG_OMP`
conditional, following the `nfsft` pattern, and do **not** add
`threads/libthreads_ng.la` to `libkernel.la` or `libkernel_threads.la`:

```make
if ENABLE_NG_OMP
  DIR_THREADS_NG = threads
else
  DIR_THREADS_NG =
endif
```

Top-level `Makefile.am`, beside the `ENABLE_OPENMP` blocks at lines 43-49 and
61-71:

```make
if ENABLE_NG_OMP
LIBNFFT3_NG_OMP_LA = libnfft3@PREC_SUFFIX@_ng_omp.la
else
LIBNFFT3_NG_OMP_LA =
endif
```

add `$(LIBNFFT3_NG_OMP_LA)` to `lib_LTLIBRARIES`, and:

```make
if ENABLE_NG_OMP
libnfft3@PREC_SUFFIX@_ng_omp_la_SOURCES =
libnfft3@PREC_SUFFIX@_ng_omp_la_LIBADD = kernel/threads/libthreads_ng.la libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @fftw3_LIBS_omp@
libnfft3@PREC_SUFFIX@_ng_omp_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@ $(OPENMP_CFLAGS) @fftw3_LDFLAGS@
libnfft3@PREC_SUFFIX@_ng_omp_la_CFLAGS = $(OPENMP_CFLAGS)
endif
```

The `LIBADD` on `libnfft3@PREC_SUFFIX@.la` is what makes this an add-on: the
serial library stays free of OpenMP and of FFTW's threads library.

`configure.ac`, after the Task 4 probe:

```m4
AM_CONDITIONAL(ENABLE_NG_OMP,
  test "x$enable_threads" = "xyes" \
    -a "x$nfft_have_fftw_planner_nthreads" = "xyes" \
    -a "x$nfft_fftw3_have_lib_omp" = "xyes")
```

and add `kernel/threads/Makefile` to `AC_CONFIG_FILES` beside the other
`kernel/*/Makefile` entries.

- [ ] **Step 6: Verify the linkage both ways**

```bash
./bootstrap.sh
./configure --enable-all --enable-tests --enable-openmp
make -j
nm -D .libs/libnfft3.so | grep -c "plan_with_nthreads"                # expect 0
nm -D .libs/libnfft3_ng_omp.so | grep -c "T nfft_plan_with_nthreads"  # expect 1
nm -D .libs/libnfft3.so | grep -c "fftw_planner_nthreads"             # expect 0
```

A non-zero last count means the Task 4 hook was bypassed somewhere. Then confirm
the serial build is unaffected:

```bash
make distclean && ./configure --enable-all --enable-tests && make -j && make check
```

- [ ] **Step 7: Commit**

```bash
clang-format -i kernel/threads/api.c kernel/threads/conf.c include/nfft3.h include/iplanner.h
git add configure.ac Makefile.am kernel/Makefile.am kernel/threads include/nfft3.h include/iplanner.h
git commit -m "Add the add-on threading library that owns the thread count and the threaded solver roster."
```

---

### Task 8: The add-on library's test binary

**Files:** Create `tests/ngomp.c`, `tests/ngomp.h`, `tests/check_ngomp.c`.
Modify `tests/Makefile.am`.

**Interfaces:**
- Consumes: the four entry points from Task 7.
- Produces: `checkall_ngomp`, built only under `ENABLE_NG_OMP`.

- [ ] **Step 1: Write the tests**

Create `tests/ngomp.c` with the GPL header, includes of `nfft3.h`, `infft.h`,
`iplanner.h`, `CUnit/CUnit.h` and `ngomp.h`, then:

```c
/* The add-on library's contract while its roster is empty. init_threads reports
 * failure, so plan_with_nthreads leaves the count at 1 and no caller can reach
 * a state with more than one thread and nothing to run on them. When the first
 * threaded solver is added to kernel/threads/conf.c these expectations invert,
 * deliberately and visibly. */
void Y(check_ngomp_empty_roster)(void)
{
  CU_ASSERT_EQUAL(Y(nfft_threads_roster_size)(), 0);
  CU_ASSERT_EQUAL(NFFT(init_threads)(), 0);

  NFFT(plan_with_nthreads)(4);
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1);

  NFFT(cleanup_threads)(); /* safe when init never succeeded */
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1);
}

/* Planning still works at every patience level through the add-on library. */
void Y(check_ngomp_plans_serially)(void)
{
  const INT N = 64, n = 128, M = 100;
  const unsigned levels[3] = {NFFT_ESTIMATE, NFFT_MEASURE, NFFT_PATIENT};
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  INT j;
  int i;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);
  for (j = 0; j < N; j++)
    f_hat[j] = K(0.0);

  NFFT(plan_with_nthreads)(4); /* refused; the count stays 1 */
  for (i = 0; i < 3; i++) {
    Y(plan_ng) *p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6,
                                       NFFT(get_window_id)(), x, (FC *)f_hat,
                                       (FC *)f, 0u, levels[i]);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p) {
      NFFT(precompute)(p);
      NFFT(execute)(p);
      NFFT(plan_ng_destroy)(p);
    }
  }

  Y(free)(f);
  Y(free)(f_hat);
  Y(free)(x);
}
```

Create `tests/ngomp.h` with the GPL header, guard `NGOMP_TEST_H`,
`#include "infft.h"`, and both declarations. Create `tests/check_ngomp.c`
following `tests/check_ng.c`: basic registry, one suite named `ngomp`, both
cases added, the same exit code and XML-report handling.

- [ ] **Step 2: Add the binary**

`tests/Makefile.am`, beside the `checkall_ng` rules:

```make
if ENABLE_NG_OMP
CHECK_NG_OMP = checkall_ngomp
else
CHECK_NG_OMP =
endif

checkall_ngomp_SOURCES = check_ngomp.c ngomp.c ngomp.h
checkall_ngomp_CFLAGS = $(OPENMP_CFLAGS)
checkall_ngomp_LDFLAGS = $(OPENMP_CFLAGS) @fftw3_LDFLAGS@ @cunit_LDFLAGS@
checkall_ngomp_LDADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_ng_omp.la $(top_builddir)/libnfft3@PREC_SUFFIX@.la @fftw3_LIBS@ @fftw3_LIBS_omp@ @cunit_LIBS@
```

Add `$(CHECK_NG_OMP)` to the `CHECK` list and to `check_PROGRAMS`, following how
`checkall_ng_threads` is wired.

- [ ] **Step 3: Run and commit**

```bash
./bootstrap.sh && ./configure --enable-all --enable-tests --enable-openmp
make -j && make check && tests/checkall_ngomp
```

A failure of `check_ngomp_empty_roster` on `init_threads` returning non-zero
means `kernel/threads/conf.c` acquired an entry it should not have.

```bash
clang-format -i tests/ngomp.c tests/check_ngomp.c
git add tests/ngomp.c tests/ngomp.h tests/check_ngomp.c tests/Makefile.am
git commit -m "Cover the add-on threading library's contract with its own test binary."
```

---

### Task 9: Wisdom configuration signature

Implements R15. Tasks 1, 4, 5 and 6 each change what a stored word or key means.

**Files:** Modify `kernel/planner/planner.c:256-272`. Test `tests/planner.c`,
`tests/planner.h`, `tests/check_ng.c`.

- [ ] **Step 1: Write the failing test**

Append to `tests/planner.c`:

```c
/* A wisdom file from before this branch must be refused: its l/u words predate
 * PLNR_NO_NONTHREADED, its keys were computed without FFTW's thread count, and
 * it carries only top-level entries. The configuration signature's vocabulary
 * tag is what refuses it. */
void Y(check_planner_wisdom_vocabulary)(void)
{
  NFFT(forget_wisdom)();

  CU_ASSERT_EQUAL(NFFT(import_wisdom_from_string)(
                       "(" STRINGIZE(Y(wisdom)) "-" PACKAGE_VERSION
                       " #x0 #x0 #x0 #x0)"),
                  0);

  {
    char *w = NFFT(export_wisdom_to_string)(); /* our own still round-trips */
    CU_ASSERT_PTR_NOT_NULL(w);
    if (w) {
      NFFT(forget_wisdom)();
      CU_ASSERT_NOT_EQUAL(NFFT(import_wisdom_from_string)(w), 0);
      NFFT(free)(w);
    }
  }
}
```

Declare in `tests/planner.h`; register as
`CU_add_test(planner_suite, "wisdom_vocabulary", Y(check_planner_wisdom_vocabulary));`.

- [ ] **Step 2: Run to verify failure**

Run: `make -j && tests/checkall_ng`
Expected: FAIL — the all-zero signature is accepted today.

- [ ] **Step 3: Tag the vocabulary**

`kernel/planner/planner.c`, in `config_signature` (line 256), after the
`sizeof(R)` word and before the registrar loop:

```c
  /* The impatience bit vocabulary and the key's composition. Bump this string
   * whenever a PLNR_* bit is added, removed or renumbered, or whenever the
   * problem hash gains or loses a field: stored words are meaningless under a
   * different vocabulary, and a mismatch must be a rejected import rather than
   * a silently misread entry. */
  Y(md5_put_str)(&m, "plnr-flags-v2");
```

- [ ] **Step 4: Run and commit**

Run: `make -j && make check` — `planner/wisdom_roundtrip` and
`planner/wisdom_rejects` must still pass.

```bash
clang-format -i kernel/planner/planner.c tests/planner.c
git add kernel/planner/planner.c tests/planner.c tests/planner.h tests/check_ng.c
git commit -m "Carry the impatience bit vocabulary in the wisdom configuration signature."
```

---

### Task 10: End-to-end behaviour and documentation

**Files:** Create `docs/adr/0005-fftw-patience-lattice-and-threading.md`. Modify
`include/nfft3.h`, `tests/nplan.c`, `tests/nplan.h`, `tests/check_ng.c`, and the
planner skill.

- [ ] **Step 1: Write the test**

Append to `tests/nplan.c`:

```c
/* Every patience level plans and executes, and each derives its own child FFTW
 * flags, which is what stops two levels sharing a wisdom key while planning
 * different child FFTs. */
void Y(check_nplan_patience_levels)(void)
{
  const INT N = 64, n = 128, M = 100;
  const unsigned levels[3] = {NFFT_ESTIMATE, NFFT_MEASURE, NFFT_PATIENT};
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  md5sig sig[3];
  planner *pl = Y(the_planner)();
  INT j;
  int i, k;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);
  for (j = 0; j < N; j++)
    f_hat[j] = K(0.0);

  for (i = 0; i < 3; i++) {
    Y(plan_ng) *p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6,
                                       NFFT(get_window_id)(), x, (FC *)f_hat,
                                       (FC *)f, 0u, levels[i]);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p) {
      NFFT(precompute)(p);
      NFFT(execute)(p);
      NFFT(plan_ng_destroy)(p);
    }
    {
      problem *q = Y(mkproblem_nfft)(
           1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, +1,
           Y(nfft_derive_fftw_flags)(levels[i], 0u), x, 1, 0, 0);
      Y(problem_md5)(pl, q, sig[i]);
      Y(problem_destroy)(q);
    }
  }

  for (i = 0; i < 3; i++)
    for (k = i + 1; k < 3; k++)
      CU_ASSERT_FALSE(sig[i][0] == sig[k][0] && sig[i][1] == sig[k][1]
                      && sig[i][2] == sig[k][2] && sig[i][3] == sig[k][3]);

  Y(free)(f);
  Y(free)(f_hat);
  Y(free)(x);
}
```

Declare in `tests/nplan.h`; register as
`CU_add_test(nplan_suite, "patience_levels", Y(check_nplan_patience_levels));`.

Run: `make -j && tests/checkall_ng` — PASS after Tasks 1 to 9. Colliding
signatures mean Task 1's derivation returns the same word for two levels.

- [ ] **Step 2: Document in the public header**

Extend the `plan_ng_guru` doc comment in `include/nfft3.h`:

```c
 * Patience. NFFT_MEASURE (the default, 0) times candidates on the caller's
 * nodes. NFFT_ESTIMATE skips timing and picks by the analytic cost model.
 * NFFT_PATIENT widens the search and, when more than one thread was requested,
 * lets serial and threaded solvers compete instead of preferring the threaded
 * one. Planning cost rises with patience. ESTIMATE overrides PATIENT if both
 * are given. Wisdom is not shared between levels.
 *
 * fftw_flags. Zero means derive the child FFTW plans' patience from the NFFT
 * level: ESTIMATE -> FFTW_ESTIMATE, MEASURE -> FFTW_MEASURE, PATIENT ->
 * FFTW_PATIENT. A non-zero word is used as given. Either way input preservation
 * is stripped and destruction forced, because the scratch grids belong to the
 * plan, and FFTW_WISDOM_ONLY is set or cleared from NFFT_WISDOM_ONLY, so
 * passing it in fftw_flags has no effect. The derived word is part of the
 * wisdom key, except for those bits, which are planning directives rather than
 * properties of the problem.
 *
 * NFFT_WISDOM_ONLY plans from the store or not at all, at every level and for
 * every internal stage, and returns NULL on a miss without writing or
 * discarding anything.
 *
 * Threads. X(plan_with_nthreads) sets the maximum number of threads a plan may
 * use; it lives in libnfft3<suffix>_ng_omp, so a program linked against
 * libnfft3<suffix> alone always plans at one thread. Both that count and FFTW's
 * own are part of the wisdom key. Below NFFT_PATIENT a serial solver declines
 * whenever more than one thread was requested, so the threaded plan is chosen
 * without being timed against the serial one -- FFTW behaves the same way, and
 * NFFT_PATIENT is how the comparison is made to happen.
```

- [ ] **Step 3: Write the ADR**

Create `docs/adr/0005-fftw-patience-lattice-and-threading.md`, shaped like
`docs/adr/0004-in-tree-html-accuracy-reports.md`, recording:

- Context: the lattice mechanics and hashed thread count existed; the policy,
  the packaging, consistent child FFTW flags, a key aware of FFTW's count, and
  tree-wide blessing did not.
- Decision: mirror FFTW. Patience is the absence of restriction; a mapping stage
  produces the bits; serial solvers decline below `PATIENT`; the thread count is
  an explicit parameter, never searched; the thread-count API and threaded
  roster ship in an add-on linked in addition, so asking for threads installs
  them; blessing reaches the winning tree by a second pass; wisdom-only is a
  planner state outside the flags word.
- Why the packaging matters: it is the whole safety argument.
  `NO_NONTHREADEDP` has FFTW's two conditions and no marker on `solver_adt`,
  because a program that has not linked the add-on cannot raise the count.
  `X(init_threads)` returning 0 on an empty roster closes the remaining window.
- The five deliberate divergences listed in the spec, each with its reason.
- Naming: `libnfft3<suffix>_ng_omp` is the FFTW-shaped add-on;
  `libnfft3<suffix>_omp` remains the legacy whole-kernel rebuild linked instead
  of the serial library. Two libraries, two models, one letter apart.
- The FFTW version floor: `fftw_planner_nthreads` is 3.3.9 and later.
- Consequences: stored wisdom is invalidated for every user, by the signature
  tag, and files now carry roughly three entries per problem instead of one.
  Estimate-only users begin writing wisdom. On OpenMP builds the key changes
  again because `nthr` no longer follows the OpenMP thread count.
- Known gaps: `NFFT_EXHAUSTIVE`'s mapping branch ships dead and uncovered;
  `PLNR_NO_UGLY`, `PLNR_NO_SLOW`, `PLNR_ALLOW_PRUNING` and `PLNR_BELIEVE_PCOST`
  are set and read by nothing; the child plans' out-of-place and destroy-input
  invariants at `nfft-nd.c:209-217` have no test; the guru's `NULL` conflates
  "no wisdom" with "bad arguments"; `X(plan_with_nthreads)` must precede
  planning and nothing enforces it.
- Not decided here: racing thread counts, a pthreads variant, the
  patience-escalation loop (`api/apiplan.c:108-125`), an error channel.

- [ ] **Step 4: Update the planner skill**

- `SKILL.md`: add `NFFT_PATIENT (1<<5)`, `NFFT_NO_NONTHREADED (1<<7)`,
  `NFFT_WISDOM_ONLY (1<<8)` to the flags table and note `1<<6` reserved; add the
  four add-on entry points with the `-lnfft3_ng_omp` requirement; correct the
  `fftw_flags` paragraph, which says `0` means `FFTW_MEASURE`, to the derive
  rule; state that `FFTW_WISDOM_ONLY` inside `fftw_flags` is ignored.
- `reference/planning-modes-and-flags.md`: the patience lattice, the decline
  rule, that below `PATIENT` a threaded plan is chosen without being timed
  against the serial one, and that wisdom is not shared between levels.
- `reference/wisdom.md`: the vocabulary tag; FFTW's thread count in the key via
  the hook; blessing now covering the whole winning tree and estimate solutions.
- `reference/solvers-problems-windows.md`: a threaded solver is registered in
  `kernel/threads/conf.c` and nowhere else, and adding the first entry is what
  makes `X(init_threads)` succeed and the decline rule take effect.
- `reference/building-testing-examples.md`: `libnfft3<suffix>_ng_omp`, how it
  differs from `libnfft3<suffix>_omp`, and `checkall_ngomp`.

- [ ] **Step 5: Run the full matrix**

```bash
./bootstrap.sh
./configure --enable-all --enable-tests --enable-exhaustive-unit-tests && make -j && make check && make distclean
./configure --enable-all --enable-tests --enable-openmp && make -j && make check && make distclean
./configure --enable-all --enable-tests --enable-float && make -j && make check && make distclean
./configure --enable-all --enable-tests --enable-long-double && make -j && make check
```

All four PASS. In the OpenMP configuration confirm `checkall`,
`checkall_threads`, `checkall_ng`, `checkall_ng_threads` and `checkall_ngomp`
are green.

- [ ] **Step 6: Format and commit**

```bash
clang-format -i include/nfft3.h tests/nplan.c
git add include/nfft3.h docs/adr/0005-fftw-patience-lattice-and-threading.md .claude/skills/understanding-the-planner-api tests/nplan.c tests/nplan.h tests/check_ng.c
git commit -m "Document the patience lattice and pin the guru's behaviour at every level."
```

---

## Self-Review

**Spec coverage.** R1-R3, R10, R18, R19 → Task 1. R4 → Task 2. R5 → Task 3.
R11 → Task 4. R12-R14 → Task 5. R16-R17 → Task 6. R7-R9 → Tasks 7 and 8.
R15 → Task 9. R6 and the divergences → Task 10's ADR.

**Ordering.** Task 5 precedes Task 6 because a recursive wisdom-only gate needs
children in wisdom. Task 4 precedes Task 7 because the add-on installs the hook.
Task 1 step 9 leaves `nplan/fftw_wisdom_only_declines` deliberately red until
Task 6 step 6; that is the only cross-task red state and it is called out at
both ends.

**Placeholder scan.** No TBD, no "handle errors appropriately", no "similar to
Task N". Every code step carries its code. Task 10 steps 3 and 4 specify prose
by required content, the right granularity for documents.

**Type consistency.** `Y(nfft_map_planning_flags)(unsigned) -> unsigned` and
`Y(nfft_derive_fftw_flags)(unsigned, unsigned) -> unsigned` are declared in
Task 1 step 3, defined in step 6, used in Tasks 1, 2, 3 and 10.
`int (*Y(fftw_nthreads_hook))(void)` is declared and defined in Task 4 and set
in Task 7. `Y(nfft_threads_roster_size)(void) -> int` is declared in Task 7
step 2, defined in step 3, read in step 4 and in Task 8.
`Y(count_wisdom_entries)(void) -> int` is defined in Task 5 step 1, declared in
`nplan.h` there, reused in Task 6 step 1. `NO_NONTHREADEDP(pl)` takes one
argument everywhere. `wisdom_state_t` and `PLNR_WISDOM_{NORMAL,ONLY,IS_BOGUS}`
are introduced in Task 6 step 3 and used in steps 4 and 5. `invoke_solver` is
introduced in Task 3 with three saved fields and extended in Task 5 to four;
both tasks show the whole function. `solver_adt` is unchanged, so
`scope_probe_adt` in Task 3 has three members.

**Known gaps, accepted during review.**
1. `NFFT_EXHAUSTIVE`'s branch in `levels()` ships dead and uncovered.
2. `PLNR_NO_UGLY`, `PLNR_NO_SLOW`, `PLNR_ALLOW_PRUNING` and
   `PLNR_BELIEVE_PCOST` are set correctly and read by nothing.
3. `nfft-nd.c:209-217` has no test; its invariants are unobservable.
4. Wisdom files grow roughly threefold, and estimate-only users start writing
   them.
5. Every existing wisdom file is invalidated, including for users who never
   touch threads.
6. `X(plan_with_nthreads)` destroys in-memory wisdom, so it must precede
   planning; nothing enforces that beyond documentation.
