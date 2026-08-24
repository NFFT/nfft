---
name: understanding-the-planner-api
description: Use when running, reading, changing, or extending NFFT3's next-generation planner API (plan_ng) — the FFTW-style self-tuning NFFT transform, nfft_plan_ng_guru lifecycle, its problem/solver/plan trinity, wisdom cache, native NDFT/fast-NFFT solvers, DECONV/CONV children, or window vtable.
---

# Understanding the planner API (`plan_ng`)

## What this is

`plan_ng` is NFFT3's **next-generation, FFTW-shaped, self-tuning NFFT API**. A
*planner* measures the applicable algorithms for a concrete problem, picks the
fastest, caches the decision (*wisdom*), and hands back a *plan* you execute.
It lives **alongside** the legacy `nfft_plan` API (which is untouched and stays
forever); users opt in. It is **clean-room** — architecture inspired by FFTW3,
implementation written without copying FFTW source (see
[solvers-problems-windows.md](reference/solvers-problems-windows.md#clean-room-protocol)).

**Scope today: NDFT + NFFT only.** The real-valued NFCT/NFST kinds were
prototyped and then **removed** from the new API; only the legacy `nfct_*` /
`nfst_*` API exists for those. Do not look for `nfct_plan_ng` — it does not
exist. (History: [history-and-drift.md](reference/history-and-drift.md).)

Ground truth is always the **current-branch code**, not the (removed) design
docs. Key files:
- `include/nfft3.h` — public API (search `NFFT_DEFINE_PLANNER_API`).
- `include/iplanner.h` — internal planner + problem/solver/plan ADTs (never installed).
- `kernel/nfft/plan_ng.c` — the NFFT bundle lifecycle (guru, precompute, execute, destroy).
- `kernel/nfft/api_ng.c` — wisdom + fprint_plan wrappers.
- `kernel/planner/` — kind-agnostic trinity core (search, wisdom, md5, tensor, printer/scanner, timer).
- `kernel/nfft/{nfast,ndft,ndft-nd,nconst}.c`, `kernel/{deconv,conv}/` — the solvers.

## The mental model (read this first)

Three object families, each with a const vtable (`*_adt`) as its first member,
FFTW-style:

- **problem** — a hashable *description of the work* ("what"). Kinds:
  `NFFT_PROBLEM_NFFT`, `NFFT_PROBLEM_DECONV`, `NFFT_PROBLEM_CONV`.
- **solver** — a stateless *recipe*. Its `mkplan()` returns a **plan** (with an
  estimated `pcost`) or `NULL` when it does not apply. Both inapplicability and
  impatience-gating express as `NULL`.
- **plan** — an *executable* "how", with `apply()` (forward) and
  `apply_adjoint()`. Building node-dependent tables (ψ) happens in `awake()`,
  never in `mkplan()`.

The **planner** (process-global, `the_planner()`) drives a search over the
registered solvers of a problem's kind, keyed by a **128-bit MD5 of the
problem's structural *size class*** (never node coordinates). It memoises
decisions as **wisdom**.

An NFFT `plan_ng` is a **bundle** wrapping one internal plan. Every candidate
plan is **coreless planner-native**: it reads/writes through the problem alone —
there is no shared "legacy core" (that architecture was removed; if a doc
mentions a "provisional/shared core", it is stale).

## The public API

```c
/* Construct. `variant` NULL = all type-I; x/f_hat/f are REQUIRED (every mode). */
nfft_plan_ng *nfft_plan_ng_guru(
    int d, const NFFT_INT *N, const int *variant, const NFFT_INT *n,
    NFFT_INT M, int m, int window,
    R *x, C *f_hat, C *f,           /* R=real, C=complex; precision-mangled */
    unsigned fftw_flags, unsigned planning);

void nfft_precompute(nfft_plan_ng *p);         /* MANDATORY before any execute */
void nfft_execute(nfft_plan_ng *p);            /* forward:  f_hat -> f          */
void nfft_execute_adjoint(nfft_plan_ng *p);    /* adjoint:  f     -> f_hat      */
void nfft_execute_on(nfft_plan_ng *p, C *f_hat, C *f);          /* new-array */
void nfft_execute_adjoint_on(nfft_plan_ng *p, C *f_hat, C *f);  /* new-array */
void nfft_fprint_plan(nfft_plan_ng *p, FILE *out);   /* print the chosen plan tree */
void nfft_plan_ng_destroy(nfft_plan_ng *p);          /* NULL-safe */

/* Wisdom + tuning — one global store, shared library-wide (FFTW model). */
int   nfft_export_wisdom_to_filename(const char *filename);
char *nfft_export_wisdom_to_string(void);          /* caller frees with nfft_free */
int   nfft_import_wisdom_from_filename(const char *filename);
int   nfft_import_wisdom_from_string(const char *s);
void  nfft_forget_wisdom(void);
void  nfft_set_timelimit(double seconds);          /* <0 = unlimited (default) */
int   nfft_get_window_id(void);                    /* compile-time window ordinal */
```

The names shown are the double-precision spellings. In precision-agnostic C use
the mangling macros (`Y(...)`, `X(...)`, `FFTW(...)`); see `CONVENTIONS.md`. The
guru returns `NULL` on bad arguments (a real release-safe guard, not just a
debug assert). Bad arguments include **insufficient oversampling**: unless
`NFFT_NO_FAST_NATIVE` is set, every axis with `N[t] > 1` must satisfy
`n[t] > N[t]`. `sigma <= 1` is refused rather than served by a direct
transform, because the guru does not reveal which solver won and a silent
fallback would cost `O(N*M)` instead of `O(n log n)`; excluding the fast solver
explicitly lifts the requirement. Unit axes (`N[t] == 1`) are elided and
exempt.

## The lifecycle contract — the #1 thing to get right

```
plan_ng_guru(...)   ->   precompute(p)   ->   execute(p) / execute_adjoint(p)
   (may race &              (materialises        (assert the winner is
    measure candidates)      ψ; SLEEPY->AWAKE)    PLNR_AWAKE)
```

1. **`precompute()` is mandatory before every `execute*`**, for *every* plan —
   fast and direct alike. `execute*` asserts the winner is `PLNR_AWAKE`.
   Skipping it is a checked error under `--enable-debug`; in a **release build
   `A(...)` asserts are no-ops**, so skipping precompute is **undefined
   behavior** (you run a `SLEEPY` plan). Almost every "planner crash" is a test
   or caller that forgot `precompute()`.
2. **`x` is copied** into a plan-owned buffer at construction; you keep
   ownership of your `x` and may free/reuse it after the guru returns. The plan
   never writes or frees your `x`.
3. **`f_hat` and `f` are aliased** (borrowed, never freed by the plan).
   **Measured planning may clobber them** during the race — so fill `f_hat`
   *after* constructing, or use the `_on` new-array variants to run on scratch.
   (Estimate mode runs no race, so it never touches your arrays — but filling
   *after* the guru is the safe habit for both modes.)
4. **Measured is the default.** `NFFT_MEASURE == 0`. Passing `NFFT_ESTIMATE`
   skips measurement and picks by an analytic cost model instead (instant, never
   time-bounded).

Full detail: [lifecycle-and-contracts.md](reference/lifecycle-and-contracts.md).

## Planning flags (current values — quote these)

| Flag | Value | Meaning |
|------|-------|---------|
| `NFFT_MEASURE` | `0` | Default: race candidates on your nodes, bless the winner. |
| `NFFT_ESTIMATE` | `1<<0` | Skip the race; pick by analytic cost model. |
| `NFFT_NO_DIRECT` | `1<<1` | Forbid the O(N·M) direct/NDFT solvers. |
| `NFFT_NO_FAST_NATIVE` | `1<<4` | Forbid the native fast NFFT (DECONV+FFT+CONV). It is the *only* fast solver, so this effectively forces a direct NDFT. |

Per-axis NDFT variant (even `N`): `NFFT_NDFT_TYPE_I` (`k=-N/2..N/2-1`),
`NFFT_NDFT_TYPE_II` (`k=-N/2+1..N/2`); odd `N` has only type-I. Windows:
`NFFT_WINDOW_{KAISER_BESSEL,GAUSSIAN,B_SPLINE,SINC_POWER}` work in the fast path;
`NFFT_WINDOW_DIRAC_DELTA` is declined by every fast solver (direct-NDFT only).

## Minimal example (1D forward + adjoint)

```c
Y(nfft_ensure_registered)();          /* implicit inside the guru; safe to call */
int window = NFFT(get_window_id)();   /* or pass a runtime NFFT_WINDOW_* ordinal */
NFFT(plan_ng) *p = NFFT(plan_ng_guru)(
    1, &N, /*variant=*/NULL, &n, M, m, window,
    x, f_hat, f, /*fftw_flags=*/0u, NFFT_ESTIMATE /* or 0 for MEASURE */);
if (!p) { /* bad args: guru returned NULL */ }
/* Fill f_hat HERE — after the guru (a measured race may clobber it), before
 * precompute. x is already copied, so you may free/reuse your x now. */
NFFT(precompute)(p);                  /* MANDATORY — moves winner SLEEPY->AWAKE */
NFFT(execute)(p);                     /* f_hat -> f  */
NFFT(fprint_plan)(p, stdout);         /* see which solver won */
NFFT(execute_adjoint_on)(p, f_hat_adj, f_in);   /* adjoint on scratch arrays */
NFFT(plan_ng_destroy)(p);
```

Working, runnable examples on the current branch:
`examples/nfft/nfast_native.c` (five-way legacy-vs-planner check, prints plan
trees) and `examples/nfft/ndft_fast.c`. Tests: `tests/checkall_ng`
(`plan.c`, `window.c`, `deconv.c`, `conv.c`, `fast_native.c`, `nfft_ng.c`). See
[building-testing-examples.md](reference/building-testing-examples.md).

## Where to go next (progressive disclosure)

- **[lifecycle-and-contracts.md](reference/lifecycle-and-contracts.md)** — awake
  states, precompute-before-execute, array ownership, `execute` vs `_on`,
  forward-only race & adjoint reuse, thread-safety, release-vs-debug asserts.
- **[planning-modes-and-flags.md](reference/planning-modes-and-flags.md)** —
  estimate vs measured, the measured race, impatience/subsumption lattice,
  evidence grades, timelimit, full flag semantics.
- **[wisdom.md](reference/wisdom.md)** — the size-class MD5 key, the wisdom
  store, export/import + file grammar, configuration signature, versioning,
  cache-safe semantics.
- **[planner-internals.md](reference/planner-internals.md)** — the trinity
  ADTs, the search (`mkplan`/`candidates`/`bless`/`hlookup`), coreless-native
  architecture, tensor/mvdim geometry, printer/scanner, md5, the timer,
  generation/ABA guard.
- **[solvers-problems-windows.md](reference/solvers-problems-windows.md)** — the
  problem kinds, the 5-solver roster + gating, the native fast decomposition,
  unit-axis elision, the window vtable + KB self-normalization, **how to add a
  solver/kind** (eager registration!), and the clean-room protocol.
- **[building-testing-examples.md](reference/building-testing-examples.md)** —
  build/configure flags, `checkall_ng`, examples, `--enable-debug` (md5 guard +
  sanitizers), CMake vs Autotools sources.
- **[history-and-drift.md](reference/history-and-drift.md)** — **read before
  trusting any older mental model.** A catalog of superseded decisions and dead
  terms (legacy/shared core, wrapper solvers, ψ-strategy `PRE_FULL_PSI` flags,
  NFCT/NFST new API, `nfft_optimize`, patience ladder, `nfft_flags` param, `x`
  aliasing) with the current correction for each.

## Red flags — you are working from a stale model if you see these

- A "shared/provisional/legacy core" per bundle → **gone**; plans are coreless native.
- Solvers named `nfft_solver_fast_1d/2d/3d`, `fast_nd`, or `nfft_solver_direct`
  → **replaced** by `fast_native` + `ndft_1d`/`ndft_nd` + `const_0d`.
- Flags `NFFT_NO_DIM_SPECIAL`, `NFFT_CONSERVE_MEMORY`, `NFFT_NO_FULL_PSI`,
  `NFFT_SHARED_CORE`, `NFFT_FORWARD_ONLY`, `NFFT_NO_FAST_WRAPPER`, or an
  `nfft_flags` guru parameter → **all removed**.
- `plan_ng_x` / `plan_ng_f_hat` / `plan_ng_f` accessors, or `x` allowed to be
  `NULL` in estimate mode → **removed**; arrays are required, passed to the guru.
- Two awake states (`PLNR_SLEEPY`/`PLNR_AWAKE`) → there are **three**
  (`SLEEPY=0 < AWAKE_ZERO=1 < AWAKE=2`).
- A separate `nfft_optimize()` call → never existed; measurement is at plan time.
