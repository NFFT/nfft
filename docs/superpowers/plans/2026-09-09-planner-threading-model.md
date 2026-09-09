# Planner Threading Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Transplant FFTW's patience lattice, its serial-versus-threaded plan
selection, and its add-on threading library shape into `plan_ng`, and make the
child FFTW plans consistent with the patience the caller asked for.

**Architecture:** Four public patience levels map, through an FFTW-style
`IMPLIES` table in a new translation unit, onto the existing `PLNR_*` bits. One
new bit, `PLNR_NO_NONTHREADED`, is set below `NFFT_PATIENT`; serial solvers
consult it and decline. Safety comes from packaging rather than bookkeeping: the
thread-count API lives only in a new add-on library `libnfft3_ng_omp`, so a
program linked against `libnfft3` alone can never raise the count above 1, which
is FFTW's own invariant. The same translation unit derives the child FFTW
planner flags from the patience level, and a null-by-default hook lets the
add-on report FFTW's thread count into the wisdom key.

**Tech Stack:** C99, GNU Autotools, CUnit (`tests/checkall_ng`), the existing
`kernel/planner/` trinity core, FFTW3 (3.3.9+ for the add-on).

**Spec:** `docs/superpowers/specs/2026-09-09-planner-threading-model.md`

## Global Constraints

- Precision-agnostic C: `Y(name)` library-wide, `X(name)` module-local,
  `FFTW(name)` for FFTW. Never hard-code an `nfft_` prefix.
- `R`, `E`, `C`, `A(...)`, `CK(...)` come from `include/infft.h`, for `kernel/`
  and `tests/` only.
- Indentation 2 spaces, BSD braces. Run `clang-format -i` on every file touched.
- The float / double / long-double build matrix must keep working.
- Legacy `nfft_plan` (`kernel/nfft/nfft.c`), `Y(set_num_threads)`,
  `libnfft3<suffix>_omp`, `kernel/nfsft/`, `kernel/fpt/` are not to be modified.
- Existing public flag values do not move: `NFFT_MEASURE = 0U`,
  `NFFT_ESTIMATE = 1U<<0`, `NFFT_NO_DIRECT = 1U<<1`,
  `NFFT_NO_FAST_NATIVE = 1U<<4`. New flags take `1U<<5` upward; bits `1U<<2`
  and `1U<<3` belonged to retired flags and stay unused.
- `libnfft3<suffix>` links `@fftw3_LIBS@` only. It must never reference a
  symbol that lives in FFTW's threads library.
- Baseline build and test:
  `./bootstrap.sh && ./configure --enable-all --enable-tests && make -j && make check`
- Commit messages: one sentence, no type prefix, no attribution lines.

---

## File Structure

**Created**

- `kernel/nfft/mapflags.c` — translation of the caller's request. Two pure
  functions: the `PLNR_*` image of the planning word, and the child FFTW
  planner flags. Mirrors FFTW's `api/mapflags.c`. Separate so both are
  unit-testable without building a plan.
- `kernel/threads/api.c` — the add-on library's entry points:
  `X(init_threads)`, `X(cleanup_threads)`, `X(plan_with_nthreads)`,
  `X(planner_nthreads)`.
- `kernel/threads/conf.c` — the threaded solver roster.
  `Y(nfft_threads_conf_standard)`. Empty today; the single place a threaded
  solver gets registered.
- `kernel/threads/Makefile.am` — builds `libthreads_ng.la`.
- `tests/ngomp.c`, `tests/ngomp.h`, `tests/check_ngomp.c` — the add-on
  library's own CUnit binary.
- `docs/adr/0005-fftw-patience-lattice-and-threading.md`

**Modified**

- `include/nfft3.h` — `NFFT_PATIENT`, `NFFT_EXHAUSTIVE`, `NFFT_NO_NONTHREADED`;
  declarations for the four add-on entry points.
- `include/iplanner.h` — `PLNR_NO_NONTHREADED`, `NO_NONTHREADEDP`, the two
  `mapflags.c` declarations, the `Y(fftw_nthreads_hook)` declaration.
- `kernel/planner/planner.c` — `invoke_solver` scoping, the hook's definition,
  the configuration signature tag.
- `kernel/nfft/problem.c` — hash FFTW's thread count.
- `kernel/nfft/plan.c` — use both mapflags functions, drop the OpenMP-derived
  count.
- `kernel/nfft/nfft-nd.c`, `ndft-1d.c`, `ndft-nd.c` — decline sites.
- `kernel/deconv/deconv-{1d,2d,3d,nd}.c`, `kernel/conv/conv-{1d,2d,3d,nd}.c` —
  decline sites.
- `configure.ac` — probe `fftw_planner_nthreads`, add
  `kernel/threads/Makefile`, `AM_CONDITIONAL(ENABLE_NG_OMP, ...)`.
- `Makefile.am` — the `libnfft3<suffix>_ng_omp.la` target.
- `kernel/Makefile.am` — add the `threads` subdirectory without pulling it into
  `libkernel.la`.
- `kernel/nfft/Makefile.am`, `CMakeLists.txt` — add `mapflags.c`.
- `tests/Makefile.am` — the `checkall_ngomp` binary.
- `tests/planner.c`, `tests/planner.h`, `tests/nfast.c`, `tests/nfast.h`,
  `tests/nplan.c`, `tests/nplan.h`, `tests/check_ng.c` — new cases.
- `.claude/skills/understanding-the-planner-api/SKILL.md` and its
  `reference/planning-modes-and-flags.md`, `reference/wisdom.md`,
  `reference/solvers-problems-windows.md`, `reference/building-testing-examples.md`.

---

### Task 1: Flag vocabulary, the mapping stage, and the derived child flags

Implements R1, R2, R3, R9.

**Files:**
- Create: `kernel/nfft/mapflags.c`
- Modify: `include/nfft3.h:845-849`, `include/iplanner.h:246-256`,
  `include/iplanner.h` (declarations near `Y(planner_mkplan)`, line 323),
  `kernel/nfft/plan.c:44-53`, `kernel/nfft/plan.c:80`,
  `kernel/nfft/plan.c:105-107`, `kernel/nfft/Makefile.am`, `CMakeLists.txt`
- Test: `tests/planner.c`, `tests/planner.h`, `tests/nfast.c`, `tests/nfast.h`,
  `tests/check_ng.c`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `unsigned Y(nfft_map_planning_flags)(unsigned planning)`
  - `unsigned Y(nfft_derive_fftw_flags)(unsigned planning, unsigned fftw_flags)`
  - public macros `NFFT_PATIENT (1U<<5)`, `NFFT_EXHAUSTIVE (1U<<6)`,
    `NFFT_NO_NONTHREADED (1U<<7)`
  - internal bit `PLNR_NO_NONTHREADED = 0x0400`

- [ ] **Step 1: Add the public flags**

In `include/nfft3.h`, replace the planning-flag block at lines 845-849 with:

```c
/* Planning flags. Patience rises MEASURE -> PATIENT -> EXHAUSTIVE; ESTIMATE
 * skips measurement entirely. Patience is expressed internally by the absence
 * of restrictions, so a more patient request searches a wider space and takes
 * longer to plan. EXHAUSTIVE implies PATIENT; ESTIMATE overrides both. */
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

In `include/iplanner.h`, the impatience enum (lines 246-256) becomes:

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

- [ ] **Step 3: Declare both mapping functions**

In `include/iplanner.h`, immediately after the `Y(planner_mkplan)` declaration
(line 323), add:

```c
/* mapflags.c: translation of the caller's request. Both are pure functions of
 * their arguments; no planner or problem state is read. */
unsigned Y(nfft_map_planning_flags)(unsigned planning);
unsigned Y(nfft_derive_fftw_flags)(unsigned planning, unsigned fftw_flags);
```

- [ ] **Step 4: Write the failing test**

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
  CU_ASSERT_TRUE(f & PLNR_NO_UGLY);
  CU_ASSERT_TRUE(f & PLNR_NO_SLOW);
  CU_ASSERT_TRUE(f & PLNR_BELIEVE_PCOST);
  CU_ASSERT_FALSE(f & PLNR_ESTIMATE);

  f = Y(nfft_map_planning_flags)(NFFT_ESTIMATE);
  CU_ASSERT_TRUE(f & PLNR_NO_NONTHREADED);
  CU_ASSERT_TRUE(f & PLNR_ESTIMATE);
  CU_ASSERT_TRUE(f & PLNR_ALLOW_PRUNING);

  f = Y(nfft_map_planning_flags)(NFFT_PATIENT); /* drops the below-PATIENT block */
  CU_ASSERT_FALSE(f & PLNR_NO_NONTHREADED);
  CU_ASSERT_FALSE(f & PLNR_BELIEVE_PCOST);
  CU_ASSERT_TRUE(f & PLNR_NO_UGLY);
  CU_ASSERT_TRUE(f & PLNR_NO_SLOW);

  f = Y(nfft_map_planning_flags)(NFFT_EXHAUSTIVE); /* implies PATIENT */
  CU_ASSERT_FALSE(f & PLNR_NO_NONTHREADED);
  CU_ASSERT_FALSE(f & PLNR_NO_UGLY);
  CU_ASSERT_FALSE(f & PLNR_NO_SLOW);

  f = Y(nfft_map_planning_flags)(NFFT_ESTIMATE | NFFT_PATIENT); /* ESTIMATE wins */
  CU_ASSERT_TRUE(f & PLNR_ESTIMATE);
  CU_ASSERT_TRUE(f & PLNR_NO_NONTHREADED);

  f = Y(nfft_map_planning_flags)(NFFT_PATIENT | NFFT_NO_NONTHREADED);
  CU_ASSERT_TRUE(f & PLNR_NO_NONTHREADED); /* the override survives PATIENT */

  f = Y(nfft_map_planning_flags)(NFFT_EXHAUSTIVE | NFFT_NO_DIRECT
                                 | NFFT_NO_FAST_NATIVE);
  CU_ASSERT_TRUE(f & PLNR_NO_DIRECT); /* gates are orthogonal to patience */
  CU_ASSERT_TRUE(f & PLNR_NO_FAST_NATIVE);
}

/* iplanner.h documents problem_nfft.fftw_flags as "0 = derive". Zero means
 * take the patience level from the NFFT request; a non-zero word is the
 * caller's own and is honoured after the input-preservation normalisation. */
void Y(check_planner_derive_fftw_flags)(void)
{
  unsigned f;

  f = Y(nfft_derive_fftw_flags)(NFFT_ESTIMATE, 0u);
  CU_ASSERT_TRUE(f & FFTW_ESTIMATE);
  CU_ASSERT_FALSE(f & FFTW_PATIENT);

  f = Y(nfft_derive_fftw_flags)(NFFT_MEASURE, 0u);
  CU_ASSERT_FALSE(f & FFTW_ESTIMATE); /* FFTW_MEASURE is zero */
  CU_ASSERT_FALSE(f & FFTW_PATIENT);
  CU_ASSERT_FALSE(f & FFTW_EXHAUSTIVE);

  f = Y(nfft_derive_fftw_flags)(NFFT_PATIENT, 0u);
  CU_ASSERT_TRUE(f & FFTW_PATIENT);
  CU_ASSERT_FALSE(f & FFTW_EXHAUSTIVE);

  f = Y(nfft_derive_fftw_flags)(NFFT_EXHAUSTIVE, 0u);
  CU_ASSERT_TRUE(f & FFTW_EXHAUSTIVE);

  /* A caller word passes through, patience notwithstanding. */
  f = Y(nfft_derive_fftw_flags)(NFFT_ESTIMATE, FFTW_EXHAUSTIVE);
  CU_ASSERT_TRUE(f & FFTW_EXHAUSTIVE);
  CU_ASSERT_FALSE(f & FFTW_ESTIMATE);

  /* Deriving decides patience and nothing else. Input preservation is
   * normalised once, where the child plan is built. */
  f = Y(nfft_derive_fftw_flags)(NFFT_MEASURE, FFTW_PRESERVE_INPUT | FFTW_PATIENT);
  CU_ASSERT_TRUE(f & FFTW_PATIENT);

  /* Distinct patience levels must produce distinct words, or two plans with
   * different child FFTs would share a wisdom key. */
  CU_ASSERT_NOT_EQUAL(Y(nfft_derive_fftw_flags)(NFFT_MEASURE, 0u),
                      Y(nfft_derive_fftw_flags)(NFFT_PATIENT, 0u));
  CU_ASSERT_NOT_EQUAL(Y(nfft_derive_fftw_flags)(NFFT_ESTIMATE, 0u),
                      Y(nfft_derive_fftw_flags)(NFFT_MEASURE, 0u));
}
```

Declare both in `tests/planner.h` before `#endif`:

```c
void Y(check_planner_mapflags)(void);
void Y(check_planner_derive_fftw_flags)(void);
```

Register in `tests/check_ng.c` beside the other planner cases:

```c
  CU_add_test(planner_suite, "mapflags", Y(check_planner_mapflags));
  CU_add_test(planner_suite, "derive_fftw_flags",
              Y(check_planner_derive_fftw_flags));
```

- [ ] **Step 5: Run the test to verify it fails**

Run: `make -j && tests/checkall_ng`
Expected: link failure, `undefined reference to nfft_nfft_map_planning_flags`.

- [ ] **Step 6: Write the mapping unit**

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

/* Translation of the caller's request, after FFTW's api/mapflags.c.
 *
 * Patience is the absence of restriction: an impatient request carries more
 * PLNR_NO_* bits and searches a narrower space. The public word names the
 * level; this file turns a level into the restrictions it implies, and into
 * the planner flags handed to the child FFTW plans. Nothing here reads planner
 * or problem state. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

/* EXHAUSTIVE implies PATIENT; ESTIMATE denies both. */
static void levels(unsigned planning, unsigned *estimate, unsigned *patient,
                   unsigned *exhaustive)
{
  *estimate = (planning & NFFT_ESTIMATE) ? 1u : 0u;
  *exhaustive = (planning & NFFT_EXHAUSTIVE) ? 1u : 0u;
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

/* The child FFTW plans' flags. A zero fftw_flags means derive from the NFFT
 * patience level, which is what problem_nfft.fftw_flags has always documented;
 * a non-zero word is the caller's own and passes through. Input preservation
 * is stripped on both paths because the scratch grids belong to the plan, and
 * FFTW_DESTROY_INPUT is added where the plan is built (kernel/nfft/nfft-nd.c),
 * not here, so it stays out of the wisdom key. */
unsigned Y(nfft_derive_fftw_flags)(unsigned planning, unsigned fftw_flags)
{
  unsigned estimate, patient, exhaustive;

  if (fftw_flags != 0u)
    return fftw_flags;

  levels(planning, &estimate, &patient, &exhaustive);

  if (estimate)
    return (unsigned)FFTW_ESTIMATE;
  if (exhaustive)
    return (unsigned)FFTW_EXHAUSTIVE;
  if (patient)
    return (unsigned)FFTW_PATIENT;
  return (unsigned)FFTW_MEASURE; /* zero */
}
```

- [ ] **Step 7: Add it to both builds**

In `kernel/nfft/Makefile.am`, add `mapflags.c` to the source list in its
alphabetical place, between `conf.c` and `ndft-1d.c`. If the file lists sources
twice (a serial and a `_threads` variant), add it to both.

In `CMakeLists.txt`, find the list naming `kernel/nfft/conf.c` and add
`kernel/nfft/mapflags.c` next to it.

- [ ] **Step 8: Wire both functions into `plan.c`**

In `kernel/nfft/plan.c`:

Delete the `map_planning_flags` function at lines 44-53. Change its call site
at line 80 to:

```c
  F = Y(nfft_map_planning_flags)(planning);
```

Change the problem construction at lines 105-107 from

```c
  p->prob[FWD] = Y(mkproblem_nfft)(d, N, variant, n, M, m, window, +1,
                                keyable_fftw_flags(fftw_flags), x,
                                /*copy_x=*/1, (C *)f_hat, (C *)f);
```

to

```c
  /* The derived word is what the child FFTW plans use and what the wisdom key
   * records, so two patience levels never share a key while planning
   * different child FFTs. keyable_fftw_flags still drops the input-preservation
   * spellings, which do not change any measured cost. */
  p->prob[FWD] = Y(mkproblem_nfft)(
       d, N, variant, n, M, m, window, +1,
       keyable_fftw_flags(Y(nfft_derive_fftw_flags)(planning, fftw_flags)), x,
       /*copy_x=*/1, (C *)f_hat, (C *)f);
```

Apply the identical change to any other `Y(mkproblem_nfft)` call in the file;
grep for `mkproblem_nfft` and confirm each passes the derived word.

- [ ] **Step 9: Pin the child plan's flag invariants against regression**

Four properties hold today at `kernel/nfft/nfft-nd.c:209-217` and must still
hold afterwards. FFTW has no in-place or out-of-place flag: in-place is
expressed by passing the same pointer twice, so the first two are pointer
properties, not flag properties.

| invariant | enforced by |
|---|---|
| the child FFTs are out of place | `g1 != g2` in the two `FFTW(plan_dft)` calls |
| forward is `g1 -> g2`, backward is `g2 -> g1` | the same two calls |
| `FFTW_DESTROY_INPUT` always set | `\| FFTW_DESTROY_INPUT` at line 210 |
| `FFTW_PRESERVE_INPUT` never set | `& ~FFTW_PRESERVE_INPUT` at line 210 |

Append to `tests/nfast.c`, which already reaches inside the fast solver:

```c
/* The child FFTW plans must stay out of place with input destruction forced,
 * whatever the caller passed and whatever the patience level. FFTW's own
 * description of the plan names the direction and whether it is in place, so
 * it is the cheapest way to assert this from outside. */
void Y(check_nfast_child_fftw_flags)(void)
{
  const INT N = 64, n = 128, M = 64;
  const unsigned levels[4] = {NFFT_ESTIMATE, NFFT_MEASURE, NFFT_PATIENT,
                              NFFT_EXHAUSTIVE};
  /* A caller word that would break both invariants if it were honoured. */
  const unsigned hostile[2] = {0u, (unsigned)FFTW_PRESERVE_INPUT};
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  INT j;
  int i, h;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);
  for (j = 0; j < N; j++)
    f_hat[j] = K(0.0);

  for (h = 0; h < 2; h++)
    for (i = 0; i < 4; i++) {
      Y(plan_ng) *p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6,
                                         NFFT(get_window_id)(), x,
                                         (FC *)f_hat, (FC *)f, hostile[h],
                                         levels[i]);
      CU_ASSERT_PTR_NOT_NULL(p);
      if (p) {
        /* The plan tree splices in FFTW's own description of the child. An
         * in-place child would say so; an out-of-place one does not. */
        char buf[4096];
        FILE *s = tmpfile();
        size_t got;
        CU_ASSERT_PTR_NOT_NULL(s);
        NFFT(precompute)(p);
        NFFT(fprint_plan)(p, s);
        rewind(s);
        got = fread(buf, 1, sizeof buf - 1, s);
        buf[got] = '\0';
        fclose(s);
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "fftw"));
        CU_ASSERT_PTR_NULL(strstr(buf, "in-place"));
        NFFT(plan_ng_destroy)(p);
      }
    }

  Y(free)(f);
  Y(free)(f_hat);
  Y(free)(x);
}
```

Declare it in `tests/nfast.h`; register as
`CU_add_test(nfast_suite, "child_fftw_flags", Y(check_nfast_child_fftw_flags));`.

Run: `make -j && tests/checkall_ng`
Expected: PASS. If `in-place` appears, `nfft-nd.c` was changed and the two
`FFTW(plan_dft)` calls no longer use distinct buffers. If the guru returns
`NULL` for the `FFTW_PRESERVE_INPUT` word, the strip at line 210 was lost.

- [ ] **Step 10: Run the tests**

Run: `./bootstrap.sh && ./configure --enable-all --enable-tests && make -j && make check`
Expected: PASS, including `planner/mapflags`, `planner/derive_fftw_flags` and
`nfast/child_fftw_flags`.

- [ ] **Step 11: Format and commit**

```bash
clang-format -i kernel/nfft/mapflags.c kernel/nfft/plan.c include/nfft3.h include/iplanner.h tests/planner.c tests/nfast.c
git add include/nfft3.h include/iplanner.h kernel/nfft/mapflags.c kernel/nfft/plan.c kernel/nfft/Makefile.am CMakeLists.txt tests/planner.c tests/planner.h tests/nfast.c tests/nfast.h tests/check_ng.c
git commit -m "Map the public planning word through an FFTW-style patience lattice and derive the child FFTW flags from it."
```

---

### Task 2: The decline rule

Implements R4.

**Files:**
- Modify: `include/iplanner.h` (macro block near `PLNR_L`, lines 305-306),
  `kernel/nfft/nfft-nd.c`, `kernel/nfft/ndft-1d.c`, `kernel/nfft/ndft-nd.c`,
  `kernel/deconv/deconv-{1d,2d,3d,nd}.c`, `kernel/conv/conv-{1d,2d,3d,nd}.c`
- Test: `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: `PLNR_NO_NONTHREADED` and `Y(nfft_map_planning_flags)` from Task 1.
- Produces: `NO_NONTHREADEDP(pl)`, FFTW's two-condition macro.

- [ ] **Step 1: Write the failing test**

`checkall_ng` links `libnfft3` alone, so `X(plan_with_nthreads)` is not
available to it. The test sets the count on a private planner directly, which
is what the add-on library will do through the API. Append to
`tests/planner.c`:

```c
/* A serial solver steps aside when the caller asked for more than one thread
 * and the patience level is below PATIENT. Mirrors FFTW kernel/ifftw.h
 * NO_NONTHREADEDP and the declines at dft/ct.c:132, dft/vrank-geq1.c:139.
 * Reached here by setting the private planner's count directly; through the
 * public API only the add-on library can raise it. */
void Y(check_planner_nonthreaded_decline)(void)
{
  planner *pl = Y(planner_create)();
  const INT N = 32, n = 64;
  problem *p = Y(mkproblem_deconv)(1, &N, 0, &n, 6,
                                   NFFT_WINDOW_KAISER_BESSEL, +1, 0, 0);
  plan *pln;

  Y(deconv_solvers_register)(pl);

  /* One thread: the rule never fires, whatever the patience. */
  pl->nthr = 1;
  pl->flags.l = pl->flags.u = Y(nfft_map_planning_flags)(NFFT_MEASURE);
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NOT_NULL(pln);
  if (pln)
    Y(plan_destroy)(pln);

  /* More than one thread, below PATIENT: the serial solver declines. */
  pl->nthr = 4;
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NULL(pln);

  /* Same count at PATIENT: the serial solver is back in the running. */
  pl->flags.l = pl->flags.u = Y(nfft_map_planning_flags)(NFFT_PATIENT);
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NOT_NULL(pln);
  if (pln)
    Y(plan_destroy)(pln);

  /* And at EXHAUSTIVE, which implies PATIENT. */
  pl->flags.l = pl->flags.u = Y(nfft_map_planning_flags)(NFFT_EXHAUSTIVE);
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NOT_NULL(pln);
  if (pln)
    Y(plan_destroy)(pln);

  /* The beyond-guru override reinstates the decline at PATIENT. */
  pl->flags.l = pl->flags.u =
       Y(nfft_map_planning_flags)(NFFT_PATIENT | NFFT_NO_NONTHREADED);
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NULL(pln);

  Y(problem_destroy)(p);
  Y(planner_destroy)(pl);
}
```

Declare in `tests/planner.h`; register as
`CU_add_test(planner_suite, "nonthreaded_decline", Y(check_planner_nonthreaded_decline));`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `make -j && tests/checkall_ng`
Expected: compile failure, `NO_NONTHREADEDP` undeclared, once the solvers use
it; before that, FAIL at the first `CU_ASSERT_PTR_NULL` because nothing
declines.

- [ ] **Step 3: Add the predicate**

In `include/iplanner.h`, next to `PLNR_L` and `PLNR_U` (lines 305-306):

```c
/* A serial solver must step aside when the caller asked for more than one
 * thread and the patience level is below PATIENT, so PLNR_NO_NONTHREADED is
 * set. FFTW's kernel/ifftw.h has the same two conditions. It is safe without a
 * third because raising the count requires the add-on library, and the add-on
 * library is where threaded solvers are registered. */
#define NO_NONTHREADEDP(pl)                                                   \
  ((PLNR_L(pl) & PLNR_NO_NONTHREADED) && (pl)->nthr > 1)
```

- [ ] **Step 4: Add the decline to every serial solver**

Eleven sites. In each, insert immediately after the existing `problem_kind`
check at the top of `mkplan`, before any allocation.

`kernel/nfft/nfft-nd.c`, `mkplan_native_fast`, after
`if (p->adt->kind != NFFT_PROBLEM_NFFT) return 0;`:

```c
  if (NO_NONTHREADEDP(pl))
    return 0; /* serial: prefer a threaded solver */
```

`kernel/nfft/ndft-1d.c` and `kernel/nfft/ndft-nd.c`: the same two lines, after
their own `NFFT_PROBLEM_NFFT` kind check.

`kernel/deconv/deconv-1d.c`, `-2d.c`, `-3d.c`, `-nd.c`: the same two lines,
after `if (p->adt->kind != NFFT_PROBLEM_DECONV) return 0;`.

`kernel/conv/conv-1d.c`, `-2d.c`, `-3d.c`, `-nd.c`: the same two lines, after
`if (p->adt->kind != NFFT_PROBLEM_CONV) return 0;`.

The DECONV solvers currently discard the planner with `(void)pl;`. Delete that
line wherever the parameter is now used.

Do **not** touch `kernel/nfft/rnk0.c`. The rank-0 base case has nothing to
parallelise and must stay the terminal solver for a fully elided problem.

- [ ] **Step 5: Run the tests**

Run: `make -j && make check`
Expected: PASS. Every existing case still passes because nothing linked against
`libnfft3` alone can raise `nthr` above 1.

- [ ] **Step 6: Format and commit**

```bash
clang-format -i include/iplanner.h kernel/nfft/nfft-nd.c kernel/nfft/ndft-1d.c kernel/nfft/ndft-nd.c kernel/deconv/deconv-1d.c kernel/deconv/deconv-2d.c kernel/deconv/deconv-3d.c kernel/deconv/deconv-nd.c kernel/conv/conv-1d.c kernel/conv/conv-2d.c kernel/conv/conv-3d.c kernel/conv/conv-nd.c tests/planner.c
git add include/iplanner.h kernel/nfft kernel/deconv kernel/conv tests/planner.c tests/planner.h tests/check_ng.c
git commit -m "Let serial solvers step aside for threaded ones below PATIENT."
```

---

### Task 3: Dynamic scoping of the thread budget

Implements R5.

**Files:**
- Modify: `kernel/planner/planner.c:547-620`
- Test: `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: nothing beyond a compiling tree.
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
Expected: FAIL at `CU_ASSERT_EQUAL(pl->nthr, 6)` — the probe's write leaks out.

- [ ] **Step 3: Add the scoping helper**

In `kernel/planner/planner.c`, above `Y(planner_mkplan)` (line 547):

```c
/* Every mkplan the search performs runs inside the caller's budget and may not
 * change it for the next candidate. A threaded solver lowers pl->nthr for its
 * children; this restores it. FFTW does the same in invoke_solver
 * (kernel/planner.c). */
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

- [ ] **Step 4: Route every search call site through it**

In `kernel/planner/planner.c` replace all three direct calls with
`invoke_solver(pl, p, s)`:

- line 568, the wisdom-hit re-plan inside `Y(planner_mkplan)`
- line 577, inside `FORALL_SOLVERS_OF_KIND` in `Y(planner_mkplan)`
- line 610, inside `FORALL_SOLVERS_OF_KIND` in `Y(planner_candidates)`

Then grep the file for `adt->mkplan(` and confirm the only remaining
occurrence is inside `invoke_solver`. Leave `kernel/nfft/plan.c`'s own
wisdom-hit re-plan alone: it is outside the search and holds bounds it set
itself.

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

### Task 4: FFTW's thread count in the problem key

Implements R10.

**Files:**
- Modify: `configure.ac` (probe), `include/iplanner.h` (hook declaration),
  `kernel/planner/planner.c` (hook definition), `kernel/nfft/problem.c:53-76`
- Test: `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `extern int (*Y(fftw_nthreads_hook))(void);`, null by default. Task
  6's `X(init_threads)` sets it.

- [ ] **Step 1: Write the failing test**

Append to `tests/planner.c`:

```c
static int fake_fftw_nthreads_value = 1;
static int fake_fftw_nthreads(void) { return fake_fftw_nthreads_value; }

/* FFTW's own thread count changes the child FFT plans, so it must change our
 * key. libnfft3 links @fftw3_LIBS@ only and fftw_planner_nthreads lives in
 * FFTW's threads library, so the count arrives through a hook that the add-on
 * library installs. Null hook means one thread. */
void Y(check_planner_fftw_nthreads_key)(void)
{
  planner *pl = Y(planner_create)();
  const INT N = 32, n = 64;
  const INT M = 10;
  R x[10];
  md5sig a, b, c;
  problem *p;
  INT j;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);

  CU_ASSERT_PTR_NULL(Y(fftw_nthreads_hook)); /* default */

  p = Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u,
                        x, 1, 0, 0);
  Y(problem_md5)(pl, p, a);

  /* An installed hook reporting one thread must not move the key. */
  Y(fftw_nthreads_hook) = fake_fftw_nthreads;
  fake_fftw_nthreads_value = 1;
  Y(problem_md5)(pl, p, b);
  CU_ASSERT_TRUE(a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3]);

  /* A different count must move it. */
  fake_fftw_nthreads_value = 8;
  Y(problem_md5)(pl, p, c);
  CU_ASSERT_FALSE(a[0] == c[0] && a[1] == c[1] && a[2] == c[2] && a[3] == c[3]);

  Y(fftw_nthreads_hook) = 0;
  Y(problem_destroy)(p);
  Y(planner_destroy)(pl);
}
```

Declare in `tests/planner.h`; register as
`CU_add_test(planner_suite, "fftw_nthreads_key", Y(check_planner_fftw_nthreads_key));`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `make -j && tests/checkall_ng`
Expected: link failure, `undefined reference to nfft_fftw_nthreads_hook`.

- [ ] **Step 3: Probe FFTW for the accessor**

In `configure.ac`, after the FFTW3 library checks that set `fftw3_LIBS_omp`
(around line 480), add:

```m4
# fftw_planner_nthreads arrived in FFTW 3.3.9. The add-on threading library
# reports FFTW's thread count into our wisdom key through it; without it that
# library is not built.
AC_CHECK_DECL([fftw_planner_nthreads],
  [AC_DEFINE([HAVE_FFTW_PLANNER_NTHREADS],[1],
     [Define to 1 if FFTW declares fftw_planner_nthreads.])
   nfft_have_fftw_planner_nthreads=yes],
  [nfft_have_fftw_planner_nthreads=no],
  [[#include <fftw3.h>]])
```

- [ ] **Step 4: Declare and define the hook**

In `include/iplanner.h`, beside the other planner declarations (near line 312):

```c
/* FFTW's own planner thread count, reported by the add-on threading library.
 * libnfft3 links @fftw3_LIBS@ only and fftw_planner_nthreads is defined in
 * FFTW's threads library, so the main library cannot call it and takes the
 * value through this hook instead. Null means one thread. Same device as
 * FFTW's threads/api.c mksolver_ct_hook. */
extern int (*Y(fftw_nthreads_hook))(void);
```

In `kernel/planner/planner.c`, near the top with the other file-scope state:

```c
int (*Y(fftw_nthreads_hook))(void) = 0;
```

- [ ] **Step 5: Hash it**

In `kernel/nfft/problem.c`, inside `hash()`, immediately after
`Y(md5_put_unsigned)(ctx, ego->fftw_flags);` (line 71), add:

```c
  /* FFTW's thread count selects a different child FFT plan for the same
   * flags, so it belongs in the key. Observed, never set. */
  Y(md5_put_int)(ctx, Y(fftw_nthreads_hook) ? Y(fftw_nthreads_hook)() : 1);
```

- [ ] **Step 6: Run the tests**

Run: `./bootstrap.sh && ./configure --enable-all --enable-tests && make -j && make check`
Expected: PASS, including `planner/fftw_nthreads_key`. Wisdom written before
this change is now a clean miss, which Task 7's signature tag makes explicit.

- [ ] **Step 7: Format and commit**

```bash
clang-format -i include/iplanner.h kernel/planner/planner.c kernel/nfft/problem.c tests/planner.c
git add configure.ac include/iplanner.h kernel/planner/planner.c kernel/nfft/problem.c tests/planner.c tests/planner.h tests/check_ng.c
git commit -m "Take FFTW's planner thread count into the problem key through a hook."
```

---

### Task 5: The add-on threading library

Implements R6, R7, R8.

**Files:**
- Create: `kernel/threads/api.c`, `kernel/threads/conf.c`,
  `kernel/threads/Makefile.am`
- Modify: `include/nfft3.h` (inside `NFFT_DEFINE_PLANNER_API`),
  `include/iplanner.h` (roster declaration), `configure.ac`, `Makefile.am`,
  `kernel/Makefile.am`
- Test: created in Task 6

**Interfaces:**
- Consumes: `Y(fftw_nthreads_hook)` from Task 4.
- Produces: `int X(init_threads)(void)`, `void X(cleanup_threads)(void)`,
  `void X(plan_with_nthreads)(int nthreads)`, `int X(planner_nthreads)(void)`,
  all defined only in `libnfft3<suffix>_ng_omp`. Internal:
  `void Y(nfft_threads_conf_standard)(planner *pl)`.

- [ ] **Step 1: Declare the entry points**

In `include/nfft3.h`, inside the `NFFT_DEFINE_PLANNER_API(X,R,C)` macro body,
after the `X(set_timelimit)` line (line 910), add — every line carries a
trailing backslash except the macro's last:

```c
NFFT_EXTERN int X(init_threads)(void); \
NFFT_EXTERN void X(cleanup_threads)(void); \
NFFT_EXTERN void X(plan_with_nthreads)(int nthreads); \
NFFT_EXTERN int X(planner_nthreads)(void); \
```

Above the macro, next to the thread-safety note at lines 838-843, add:

```c
/* Threading. These four are declared here but DEFINED ONLY in the add-on
 * library libnfft3<suffix>_ng_omp, which is linked in addition to
 * libnfft3<suffix>: `-lnfft3_ng_omp -lnfft3`. FFTW packages its threading the
 * same way. A program that links only libnfft3<suffix> cannot raise the thread
 * count above 1, which is exactly why a serial plan is always available there.
 *
 * Note libnfft3<suffix>_omp is a different thing: a whole-library OpenMP
 * rebuild for the legacy API, linked INSTEAD of libnfft3<suffix>. */
```

- [ ] **Step 2: Declare the roster**

In `include/iplanner.h`, next to `Y(nfft_ensure_registered)` (line 449):

```c
/* kernel/threads/conf.c, in the add-on library only: the threaded solver
 * roster. Empty until the first threaded solver exists. */
void Y(nfft_threads_conf_standard)(planner *pl);
```

- [ ] **Step 3: Write the roster**

Create `kernel/threads/conf.c` with the standard GPL header, then:

```c
/* The threaded solver roster. This is the one place a threaded solver is
 * registered, and it exists only in the add-on library, so a program that does
 * not link that library has no threaded solvers and cannot raise the planner's
 * thread count. FFTW's threads/conf.c has the same job.
 *
 * The table is empty: no threaded solver has been written yet. Y(init_threads)
 * reports that by returning 0, so nothing can reach a state with more than one
 * thread and nothing to run on them. Adding the first entry here is all that is
 * needed to switch the machinery on. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

static const solvtab s = {SOLVTAB_END};

void Y(nfft_threads_conf_standard)(planner *pl) { Y(solvtab_exec)(s, pl); }

/* Whether the roster registered anything. Y(init_threads) fails when it did
 * not. Counting entries rather than testing a compile-time constant keeps this
 * honest when the table grows. */
int Y(nfft_threads_roster_size)(void)
{
  int n = 0;
  while (s[n].reg != 0)
    n++;
  return n;
}
```

Add to `include/iplanner.h` beside the previous declaration:

```c
int Y(nfft_threads_roster_size)(void);
```

- [ ] **Step 4: Write the API**

Create `kernel/threads/api.c` with the standard GPL header, then:

```c
/* The add-on threading library's entry points, after FFTW's threads/api.c.
 *
 * Asking for threads is what installs them: X(plan_with_nthreads) initialises
 * threading first, and initialising registers the threaded roster. FFTW relies
 * on the same order, which is why its NO_NONTHREADEDP needs no third
 * condition. Registering solvers changes the configuration signature and so
 * invalidates every stored decision, so initialisation destroys the planner
 * first, as FFTW's X(cleanup)() does. */

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

  /* The roster is about to change; the planner and its wisdom are stale. */
  Y(the_planner_destroy)();
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

/* The maximum number of threads a plan may use, not a target: a threaded
 * solver may use fewer and passes the rest of the budget to its children. The
 * count is part of the wisdom key, so changing it makes existing entries a
 * clean miss. Leaves the count at 1 while no threaded solver exists. */
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

- [ ] **Step 5: Build the library**

Create `kernel/threads/Makefile.am`:

```make
AM_CPPFLAGS = -I$(top_srcdir)/include @fftw3_CPPFLAGS@

noinst_LTLIBRARIES = libthreads_ng.la

libthreads_ng_la_SOURCES = api.c conf.c
libthreads_ng_la_CFLAGS = $(OPENMP_CFLAGS)
```

In `kernel/Makefile.am`, add `threads` to `SUBDIRS` **only** under the OpenMP
conditional, and do **not** add `threads/libthreads_ng.la` to `libkernel.la` or
`libkernel_threads.la`. Follow the existing conditional-subdirectory pattern
used for `nfsft`:

```make
if ENABLE_NG_OMP
  DIR_THREADS_NG = threads
else
  DIR_THREADS_NG =
endif
```

and append `$(DIR_THREADS_NG)` to `SUBDIRS`.

In the top-level `Makefile.am`, beside the existing `ENABLE_OPENMP` block at
lines 43-49 and 61-71:

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

In `configure.ac`, after the probe added in Task 4:

```m4
AM_CONDITIONAL(ENABLE_NG_OMP,
  test "x$enable_threads" = "xyes" \
    -a "x$nfft_have_fftw_planner_nthreads" = "xyes" \
    -a "x$nfft_fftw3_have_lib_omp" = "xyes")
```

and add `kernel/threads/Makefile` to the `AC_CONFIG_FILES` list beside the
other `kernel/*/Makefile` entries.

- [ ] **Step 6: Verify it builds and links both ways**

```bash
./bootstrap.sh
./configure --enable-all --enable-tests --enable-openmp
make -j
ls .libs/libnfft3_ng_omp.so
nm -D .libs/libnfft3.so | grep -c "plan_with_nthreads"      # expect 0
nm -D .libs/libnfft3_ng_omp.so | grep -c "T nfft_plan_with_nthreads"  # expect 1
nm -D .libs/libnfft3.so | grep -c "fftw_planner_nthreads"   # expect 0
```

Expected: the add-on exports the four entry points, the serial library exports
none of them and references no FFTW threads symbol. If the last check is
non-zero, the hook of Task 4 was bypassed somewhere.

Then confirm the serial build is unaffected:

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

### Task 6: The add-on library's test binary

Implements the R7 and R8 half of the acceptance clause.

**Files:**
- Create: `tests/ngomp.c`, `tests/ngomp.h`, `tests/check_ngomp.c`
- Modify: `tests/Makefile.am`

**Interfaces:**
- Consumes: the four entry points from Task 5.
- Produces: the `checkall_ngomp` binary, built only under `ENABLE_NG_OMP`.

- [ ] **Step 1: Write the failing test**

Create `tests/ngomp.c` with the standard GPL header, then:

```c
#include <stdio.h>
#include <stdlib.h>
#include <CUnit/CUnit.h>

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "ngomp.h"

/* The add-on library's contract while its roster is empty. init_threads
 * reports failure, so plan_with_nthreads leaves the count at 1 and no caller
 * can reach a state with more than one thread and nothing to run on them.
 * When the first threaded solver is added to kernel/threads/conf.c these
 * expectations invert, deliberately and visibly. */
void Y(check_ngomp_empty_roster)(void)
{
  CU_ASSERT_EQUAL(Y(nfft_threads_roster_size)(), 0);
  CU_ASSERT_EQUAL(NFFT(init_threads)(), 0);

  NFFT(plan_with_nthreads)(4);
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1);

  NFFT(cleanup_threads)(); /* must be safe when init never succeeded */
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1);
}

/* Planning still works, at every patience level, through the add-on library. */
void Y(check_ngomp_plans_serially)(void)
{
  const INT N = 64, n = 128, M = 100;
  const unsigned levels[4] = {NFFT_ESTIMATE, NFFT_MEASURE, NFFT_PATIENT,
                              NFFT_EXHAUSTIVE};
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
  for (i = 0; i < 4; i++) {
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

Create `tests/ngomp.h` with the standard GPL header, an include guard
`NGOMP_TEST_H`, `#include "infft.h"`, and the two declarations.

Create `tests/check_ngomp.c` following the shape of `tests/check_ng.c`: CUnit
basic registry, one suite named `ngomp`, both cases added, and the same exit
code and XML-report handling `check_ng.c` uses.

- [ ] **Step 2: Add the binary**

In `tests/Makefile.am`, beside the existing `checkall_ng` rules:

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

and add `$(CHECK_NG_OMP)` to the `CHECK` list and to `check_PROGRAMS`,
following exactly how `checkall_ng_threads` is wired.

- [ ] **Step 3: Run it**

```bash
./bootstrap.sh
./configure --enable-all --enable-tests --enable-openmp
make -j && make check
tests/checkall_ngomp
```

Expected: both `ngomp` cases pass. If `check_ngomp_empty_roster` fails on
`init_threads` returning non-zero, `kernel/threads/conf.c` acquired an entry it
should not have.

- [ ] **Step 4: Commit**

```bash
clang-format -i tests/ngomp.c tests/check_ngomp.c
git add tests/ngomp.c tests/ngomp.h tests/check_ngomp.c tests/Makefile.am
git commit -m "Cover the add-on threading library's contract with its own test binary."
```

---

### Task 7: Wisdom safety

Implements the wisdom half of R1 to R4 and R9 to R10: every one of them changes
what a stored `l`/`u` word or a stored key means.

**Files:**
- Modify: `kernel/planner/planner.c:256-272` (`config_signature`)
- Test: `tests/planner.c`, `tests/planner.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: nothing.
- Produces: a configuration signature that changes whenever the impatience bit
  vocabulary or the key's composition changes.

- [ ] **Step 1: Write the failing test**

Append to `tests/planner.c`:

```c
/* A wisdom file whose flag vocabulary predates PLNR_NO_NONTHREADED must be
 * refused: its l/u words would read as more patient than they were, and its
 * keys were computed without FFTW's thread count. The configuration signature
 * carries the vocabulary tag that does it. */
void Y(check_planner_wisdom_vocabulary)(void)
{
  NFFT(forget_wisdom)();

  /* All-zero signature words: a plausible file from before the change. */
  CU_ASSERT_EQUAL(NFFT(import_wisdom_from_string)(
                       "(" STRINGIZE(Y(wisdom)) "-" PACKAGE_VERSION
                       " #x0 #x0 #x0 #x0)"),
                  0);

  /* Our own export still round-trips. */
  {
    char *w = NFFT(export_wisdom_to_string)();
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

- [ ] **Step 2: Run the test to verify it fails**

Run: `make -j && tests/checkall_ng`
Expected: FAIL — the all-zero signature is accepted today, because nothing in
the signature depends on the flag vocabulary.

- [ ] **Step 3: Tag the vocabulary**

In `kernel/planner/planner.c`, inside `config_signature` (line 256), after the
`sizeof(R)` word and before the registrar loop:

```c
  /* The impatience bit vocabulary and the key's composition. Bump this string
   * whenever a PLNR_* bit is added, removed or renumbered, or whenever the
   * problem hash gains or loses a field: stored words are meaningless under a
   * different vocabulary, and a mismatch must be a rejected import rather than
   * a silently misread entry. */
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

### Task 8: End-to-end behaviour and documentation

**Files:**
- Create: `docs/adr/0005-fftw-patience-lattice-and-threading.md`
- Modify: `include/nfft3.h` (guru doc comment), `tests/nplan.c`,
  `tests/nplan.h`, `tests/check_ng.c`, and the four planner skill documents
- Test: `tests/nplan.c`

**Interfaces:**
- Consumes: everything from Tasks 1 to 7.
- Produces: no new symbols.

- [ ] **Step 1: Write the failing test**

Append to `tests/nplan.c`:

```c
/* Every patience level plans and executes through the serial library, and each
 * level derives its own child FFTW flags. The last part is what stops two
 * levels sharing a wisdom key while planning different child FFTs. */
void Y(check_nplan_patience_levels)(void)
{
  const INT N = 64, n = 128, M = 100;
  const unsigned levels[4] = {NFFT_ESTIMATE, NFFT_MEASURE, NFFT_PATIENT,
                              NFFT_EXHAUSTIVE};
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  md5sig sig[4];
  planner *pl = Y(the_planner)();
  INT j;
  int i, k;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);
  for (j = 0; j < N; j++)
    f_hat[j] = K(0.0);

  for (i = 0; i < 4; i++) {
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

  /* ESTIMATE, MEASURE, PATIENT and EXHAUSTIVE derive four different child
   * words, so four different keys. */
  for (i = 0; i < 4; i++)
    for (k = i + 1; k < 4; k++)
      CU_ASSERT_FALSE(sig[i][0] == sig[k][0] && sig[i][1] == sig[k][1]
                      && sig[i][2] == sig[k][2] && sig[i][3] == sig[k][3]);

  Y(free)(f);
  Y(free)(f_hat);
  Y(free)(x);
}
```

Declare in `tests/nplan.h`; register as
`CU_add_test(nplan_suite, "patience_levels", Y(check_nplan_patience_levels));`.

- [ ] **Step 2: Run the test to verify it fails, then passes**

Run: `make -j && tests/checkall_ng`
Expected: after Tasks 1 to 7 it passes. If two signatures collide, the
derivation in Task 1 step 6 is returning the same word for two levels.

- [ ] **Step 3: Document the levels in the public header**

In `include/nfft3.h`, extend the `plan_ng_guru` doc comment:

```c
 * Patience. NFFT_MEASURE (the default, 0) times candidates on the caller's
 * nodes. NFFT_ESTIMATE skips timing and picks by the analytic cost model.
 * NFFT_PATIENT widens the search and, when more than one thread was
 * requested, lets serial and threaded solvers compete instead of preferring
 * the threaded one; NFFT_EXHAUSTIVE widens it further and implies
 * NFFT_PATIENT. Planning cost rises with patience. ESTIMATE overrides
 * PATIENT if both are given.
 *
 * fftw_flags. Zero means derive the child FFTW plans' patience from the NFFT
 * level: ESTIMATE -> FFTW_ESTIMATE, MEASURE -> FFTW_MEASURE, PATIENT ->
 * FFTW_PATIENT, EXHAUSTIVE -> FFTW_EXHAUSTIVE. A non-zero word is used as
 * given. Either way input preservation is stripped and destruction forced,
 * because the scratch grids belong to the plan. The derived word is part of
 * the wisdom key.
 *
 * Threads. X(plan_with_nthreads) sets the maximum number of threads a plan
 * may use; it lives in the add-on library libnfft3<suffix>_ng_omp, so a
 * program linked against libnfft3<suffix> alone always plans at one thread.
 * The count is part of the wisdom key, as is FFTW's own count when the add-on
 * library is linked. Below NFFT_PATIENT a serial solver declines whenever more
 * than one thread was requested, so the threaded plan is chosen without being
 * timed against the serial one -- FFTW behaves the same way, and NFFT_PATIENT
 * is how the comparison is made to happen.
```

- [ ] **Step 4: Write the ADR**

Create `docs/adr/0005-fftw-patience-lattice-and-threading.md` following the
shape of `docs/adr/0004-in-tree-html-accuracy-reports.md`. It must record:

- Context: `plan_ng` already had the lattice mechanics and a hashed thread
  count, but no patience levels above `MEASURE`, no policy for choosing between
  serial and threaded plans, child FFTW flags inconsistent with the requested
  patience, and a wisdom key blind to FFTW's own thread count. Threaded solvers
  are coming, the convolution first.
- Decision: mirror FFTW. Patience is the absence of restriction; the bits come
  from a mapping stage; serial solvers decline below `PATIENT`; the thread count
  is an explicit parameter that is never searched; the thread-count API and the
  threaded roster ship in an add-on library linked in addition to the serial
  one, so asking for threads is what installs them.
- Why the packaging matters: it is the whole safety argument.
  `NO_NONTHREADEDP` has only FFTW's two conditions, with no marker on
  `solver_adt` and no per-kind counter, because a program that has not linked
  the add-on cannot raise the count. `X(init_threads)` returning 0 on an empty
  roster closes the remaining window.
- Naming: `libnfft3<suffix>_ng_omp` is the FFTW-shaped add-on;
  `libnfft3<suffix>_omp` remains the legacy whole-kernel OpenMP rebuild linked
  instead of the serial library. Two libraries, two models, one letter apart —
  say so plainly wherever either is mentioned.
- The FFTW version floor: `fftw_planner_nthreads` is 3.3.9 and later, so the
  add-on is not built against older FFTW, and the hook then stays null.
- Consequences: no behaviour change for a program linked against `libnfft3`
  alone, apart from the child FFTW flags now following patience and the key now
  admitting FFTW's count. Both invalidate stored wisdom, which the
  configuration signature turns into a clean miss. On OpenMP builds the key
  changes again because `nthr` no longer follows the OpenMP thread count.
- Not decided here: racing thread counts, a pthreads variant, and consumers for
  `PLNR_NO_UGLY`, `PLNR_NO_SLOW`, `PLNR_ALLOW_PRUNING` and
  `PLNR_BELIEVE_PCOST`, which the mapping sets and nothing reads.

- [ ] **Step 5: Update the planner skill**

- `SKILL.md`: extend the planning-flags table with `NFFT_PATIENT (1<<5)`,
  `NFFT_EXHAUSTIVE (1<<6)`, `NFFT_NO_NONTHREADED (1<<7)`; add the four add-on
  entry points to the API listing with a note that they need
  `-lnfft3_ng_omp`; correct the `fftw_flags` paragraph, which currently says
  `0` means `FFTW_MEASURE`, to the derive rule.
- `reference/planning-modes-and-flags.md`: add the patience lattice and the
  decline rule, stating that below `PATIENT` a threaded plan is chosen without
  being timed against the serial one.
- `reference/wisdom.md`: the configuration signature now carries the impatience
  vocabulary tag; the key now carries FFTW's thread count via the hook; `nthr`
  was already there.
- `reference/solvers-problems-windows.md`: a new threaded solver is registered
  in `kernel/threads/conf.c` and nowhere else, and adding the first entry is
  what makes `X(init_threads)` succeed and the decline rule take effect.
- `reference/building-testing-examples.md`: `libnfft3<suffix>_ng_omp`, how it
  differs from `libnfft3<suffix>_omp`, and the `checkall_ngomp` binary.

Also correct the SKILL.md claim that the internal FFTW child plan "defaults to
`FFTW_ESTIMATE`". It never did; it now derives from patience.

- [ ] **Step 6: Run the full matrix**

```bash
./bootstrap.sh
./configure --enable-all --enable-tests --enable-exhaustive-unit-tests && make -j && make check && make distclean
./configure --enable-all --enable-tests --enable-openmp && make -j && make check && make distclean
./configure --enable-all --enable-tests --enable-float && make -j && make check && make distclean
./configure --enable-all --enable-tests --enable-long-double && make -j && make check
```

Expected: PASS in all four. The OpenMP configuration is the one that matters
most: confirm `checkall`, `checkall_threads`, `checkall_ng`,
`checkall_ng_threads` and `checkall_ngomp` are all green.

- [ ] **Step 7: Format and commit**

```bash
clang-format -i include/nfft3.h tests/nplan.c
git add include/nfft3.h docs/adr/0005-fftw-patience-lattice-and-threading.md .claude/skills/understanding-the-planner-api tests/nplan.c tests/nplan.h tests/check_ng.c
git commit -m "Document the patience lattice and pin the guru's behaviour at every level."
```

---

## Self-Review

**Spec coverage.** R1, R2, R3 and R9 are Task 1. R4 is Task 2. R5 is Task 3.
R10 is Task 4. R6, R7 and R8 are Task 5, tested in Task 6. The wisdom clause
common to all of them is Task 7. The acceptance clause's documentation and
end-to-end cases are Task 8.

**Placeholder scan.** No TBD, no "handle errors appropriately", no "similar to
Task N". Every code step carries its code. Task 8 steps 4 and 5 specify prose
deliverables by required content, which is the right granularity for documents.

**Type consistency.** `Y(nfft_map_planning_flags)(unsigned) -> unsigned` is
declared in Task 1 step 3, defined in step 6, used in step 8 and in the tests
of Tasks 2, 3 and 8. `Y(nfft_derive_fftw_flags)(unsigned, unsigned) -> unsigned`
follows the same path and is used again in Task 8 step 1.
`int (*Y(fftw_nthreads_hook))(void)` is declared in Task 4 step 4, defined in
the same step, read in step 5 and set in Task 5 step 4.
`Y(nfft_threads_roster_size)(void) -> int` is declared and defined in Task 5
step 3 and read in step 4 and in Task 6 step 1. `NO_NONTHREADEDP(pl)` takes one
argument everywhere, matching FFTW. `solver_adt` is unchanged, so the
`scope_probe_adt` initialiser in Task 3 has three members.

**Ordering.** Task 2 introduces a decline that nothing can trigger until Task 5
exists, and Task 4 introduces a hook that nothing sets until Task 5 sets it.
Both are tested in isolation before then, so no task depends on a later one.

**Open points for the grilling.**

1. Four of the six lattice bits stay inert. After this work the mapping sets
   `PLNR_NO_UGLY`, `PLNR_NO_SLOW`, `PLNR_ALLOW_PRUNING` and
   `PLNR_BELIEVE_PCOST` correctly and no solver reads any of them. Task 1's
   test asserts what the mapping produces, which is coverage of the mapping and
   of no behaviour.
2. The `NFFT_PATIENT` and `NFFT_EXHAUSTIVE` levels differ from `NFFT_MEASURE`
   only in bits nothing reads, plus the child FFTW word. So today
   `NFFT_EXHAUSTIVE` costs strictly more planning time in the FFTW child for no
   NFFT-side benefit. Whether to ship a level that only slows things down is a
   fair question.
3. Task 4 changes the wisdom key for every existing user, including those who
   never touch threads, because the hash gains a field. Task 7 makes that a
   clean miss rather than a wrong hit, but it is a silent full invalidation.
4. `libnfft3_ng_omp` and `libnfft3_omp` differ by one letter and mean opposite
   things — one is linked in addition, the other instead. That is a naming
   hazard no amount of documentation fully removes.
5. Task 5 has `X(init_threads)` destroy the planner, which discards in-memory
   wisdom. FFTW does the same, but FFTW's `plan_with_nthreads` is documented as
   an early call. A caller who plans, then calls `plan_with_nthreads`, loses
   every earlier decision.
