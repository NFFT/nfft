# Fast NFFT type-II and odd-N Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the planner-native fast NFFT (`nfft_solver_fast_native`) accept every
geometry the native direct NDFT solvers already accept -- per-axis type-II, odd axis
lengths, and mixed axes -- and require oversampling `sigma = n/N > 1`.

**Architecture:** The type-II frequency shift and odd-`N` parity enter the fast
decomposition `DECONV -> FFT -> CONV` only through DECONV. CONV, the window vtable and
the FFTW child plan read the index set through no path at all. So each DECONV solver
(1d, 2d, 3d, nd) gets one per-axis slot split `Nneg/Npos` replacing the hard-coded
`N/2` halves, and `guards_ok` in the fast NFFT solver drops its odd, type-II and
unit-axis declines and gains a `n > N` decline.

**Tech Stack:** C99, precision-agnostic (`Y()`/`X()`/`FFTW()` mangling, `R`/`E`/`C`
aliases, `A(...)` asserts). GNU Autotools. CUnit tests in `tests/checkall_ng`.
FFTW3. CodSpeed/CMake for benchmarks.

**Spec:** `docs/plans/nfft-fast-typeii-odd-n.md`

## Global Constraints

- **Slot map (the one rule).** Per axis, with `s = 1` for even-`N` type-II and `s = 0`
  otherwise:
  ```
  ks   = 0 .. N-1                      f_hat slot index
  Nneg = N / 2 - s                      slots  0 .. Nneg-1  hold negative frequencies
  Npos = N - Nneg                       slots  Nneg .. N-1  hold frequency >= 0
  k(ks)   = ks - Nneg                   the frequency at slot ks
  pos(ks) = (ks < Nneg) ? n - Nneg + ks : ks - Nneg     the grid cell
  ```
  Even type-I gives `Nneg == Npos == N/2`, so every existing code path is the special
  case. `mkproblem_deconv` normalizes odd `N` to `NFFT_NDFT_TYPE_I`, so
  `variant == NFFT_NDFT_TYPE_II` always implies even `N` and `Nneg >= 1`.
- **Window table start.** `Y(window_phi_hut_apply)(window, n, N, m, k0, out, count)`
  fills `out[i] = phi_hut(k0 + i)`. Every axis uses `k0 = -Nneg`.
- **Oversampling.** The fast NFFT requires `n > N` strictly on every axis. The DECONV
  solvers themselves only require `n >= N` and decline below that.
- **Name mangling.** Never write a literal `nfft_` / `nfftf_` / `nfftl_` prefix. Use
  `Y(foo)` for library-wide names, `X(foo)` for module-local ones, `FFTW(foo)` for
  FFTW.
- **Formatting.** 2-space indent, BSD braces. Run `clang-format -i <file>` on every
  touched C file before committing (uses the repo `.clang-format`).
- **Comments.** Run the `deslop` skill over every C file you touch before committing.
  No task/item/review numbers, no porting provenance, no restating the code.
- **No public API change.** `include/nfft3.h` is not edited by any task. No new source
  file, no `Makefile.am` or `CMakeLists.txt` change.
- **Build once, then rebuild per task:**
  ```bash
  ./bootstrap.sh                       # only if not already configured
  ./configure --enable-all --enable-openmp --enable-tests
  make -j
  ```
- **Run the tests:**
  ```bash
  make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
  grep -B6 "<FAILURE" CUnitAutomated_ng-Results.xml | head -60
  ```
  `tests/checkall_ng` prints nothing on success and exits 0. On failure it exits 1 and
  the XML names the failing suite/test/assertion.
- **Precision matrix.** After Task 8, `make check` must pass in double and float. For
  long double, build only -- the test run is too slow on this host.

---

## File Structure

| File | Responsibility | Task |
|---|---|---|
| `kernel/deconv/deconv-1d.c` | rank-1 DECONV: slot split, contiguous zero-pad run | 1 |
| `kernel/deconv/deconv-2d.c` | rank-2 DECONV: 2x2 block runs, folded forward/adjoint | 2 |
| `kernel/deconv/deconv-3d.c` | rank-3 DECONV: 2x2x2 block runs, folded forward/adjoint | 3 |
| `kernel/deconv/deconv-nd.c` | rank >= 4 DECONV: generalized carry odometer | 4 |
| `kernel/nfft/nfft-nd.c` | `guards_ok`: the fast NFFT applicability gate | 5 |
| `tests/nfast.c` + `tests/nfast.h` | DECONV per-rank index tests; 2D reference roster | 1-4, 6 |
| `tests/nplan.c` + `tests/nplan.h` | guard declines; odd/type-II acceptance through the fast path | 5, 7 |
| `tests/nplan_data.c` | reference-data mirror pass through the fast path | 7 |
| `tests/check_ng.c` | registers the new test functions | 1-4 |
| `docs/plans/nfft-fast-typeii-odd-n.md`, `CONTEXT.md` | close out the spec | 8 |

Tasks 1-4 are independent of each other: each rewrites one DECONV solver and is
verified directly through `Y(planner_mkplan)` on a DECONV problem, without going
through the NFFT guard. Task 5 flips the guard and must land after 1-4, because it is
the step that first routes odd and type-II NFFT problems into DECONV.

---

### Task 1: rank-1 DECONV -- slot split, odd zero-pad, `n >= N` decline

**Files:**
- Modify: `kernel/deconv/deconv-1d.c`
- Test: `tests/nfast.c` (new function), `tests/nfast.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: `Y(mkproblem_deconv)(d, N, variant, n, m, window, sign, f_hat, g)`,
  `Y(planner_mkplan)(planner *, problem *)`, `Y(plan_awake)(plan *, int)`,
  `Y(window_phi_hut)(window, n, N, m, k)`, `Y(problem_deconv_variant)(p, t)`,
  `Y(problem_deconv_n)(p, t)`, `Y(problem_deconv_N)(p, t)`.
- Produces: three file-static test helpers in `tests/nfast.c` reused by Tasks 2-4:
  `static plan *deconv_awake_plan(problem *p)`,
  `static R deconv_scale(int w, int d, const INT *N, const INT *n, int m, const INT *k)`,
  `static INT count_above(const C *v, INT len, R eps)`.
  Public test entry `void Y(check_nfast_deconv_1d_general)(void)`.

- [ ] **Step 1: Write the failing test**

Add the three shared helpers plus the rank-1 test to `tests/nfast.c`, directly after
the existing `Y(check_nfast_deconv_solver)` function:

```c
/* Plan a DECONV problem through the process-global planner and wake it. */
static plan *deconv_awake_plan(problem *p) {
  plan *pln = Y(planner_mkplan)(Y(the_planner)(), p);
  if (pln)
    Y(plan_awake)
  (pln, PLNR_AWAKE);
  return pln;
}

/* The value a unit spike at frequency k must produce: prod_t 1/phi_hut(k[t]). */
static R deconv_scale(int w, int d, const INT *N, const INT *n, int m,
                      const INT *k) {
  R s = K(1.0);
  int t;
  for (t = 0; t < d; t++)
    s /= Y(window_phi_hut)(w, n[t], N[t], m, k[t]);
  return s;
}

/* Cells whose magnitude exceeds eps. A spike must produce exactly one. */
static INT count_above(const C *v, INT len, R eps) {
  INT j, c = 0;
  for (j = 0; j < len; j++)
    if (CABS(v[j]) > eps)
      c++;
  return c;
}

/* Rank-1 DECONV for odd N and for type-II: the slot -> grid map is
 * pos(ks) = ks < Nneg ? n - Nneg + ks : ks - Nneg, with Nneg = N/2 - type_ii.
 * The grid is pre-dirtied with the expected peak so a stale zero-pad cell shows
 * up as a second above-threshold cell, not as a silent leftover. */
void Y(check_nfast_deconv_1d_general)(void) {
  const int m = 6, w = Y(get_window_id)();
  Y(deconv_ensure_registered)
  ();

  /* (a) odd N = 15, n = 32: Nneg = 7, Npos = 8. */
  {
    INT N = 15, n = 32, ks;
    INT spikes[3] = {0, 7, 14};   /* slots */
    INT freqs[3] = {-7, 0, 7};    /* k(ks) = ks - 7 */
    INT poss[3] = {25, 0, 7};     /* pos(ks) */
    int i;
    C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *g = (C *)Y(malloc)((size_t)n * sizeof(C));
    for (i = 0; i < 3; i++) {
      R sc = deconv_scale(w, 1, &N, &n, m, &freqs[i]);
      problem *p;
      plan *pln;
      for (ks = 0; ks < N; ks++)
        f_hat[ks] = K(0.0);
      f_hat[spikes[i]] = K(1.0);
      for (ks = 0; ks < n; ks++)
        g[ks] = sc; /* pre-dirty at the peak magnitude */
      p = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, 1, f_hat, g);
      pln = deconv_awake_plan(p);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
      pln->adt->apply(pln, p);
      CU_ASSERT(CABS(g[poss[i]] - sc) < K(1e-9) * sc);
      CU_ASSERT_EQUAL(count_above(g, n, K(1e-9) * sc), (INT)1);
      Y(plan_destroy)
      (pln);
      Y(problem_destroy)
      (p);
    }
    /* adjoint: g[25] = 1 gathers into slot 0 only. */
    {
      INT k = -7;
      R sc = deconv_scale(w, 1, &N, &n, m, &k);
      problem *p;
      plan *pln;
      for (ks = 0; ks < N; ks++)
        f_hat[ks] = K(0.0);
      for (ks = 0; ks < n; ks++)
        g[ks] = K(0.0);
      g[25] = K(1.0);
      p = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, -1, f_hat, g);
      pln = deconv_awake_plan(p);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
      pln->adt->apply_adjoint(pln, p);
      CU_ASSERT(CABS(f_hat[0] - sc) < K(1e-9) * sc);
      CU_ASSERT_EQUAL(count_above(f_hat, N, K(1e-9) * sc), (INT)1);
      Y(plan_destroy)
      (pln);
      Y(problem_destroy)
      (p);
    }
    Y(free)
    (f_hat);
    Y(free)
    (g);
  }

  /* (b) even N = 16 type-II, n = 32: Nneg = 7, Npos = 9. */
  {
    INT N = 16, n = 32, ks;
    int v = NFFT_NDFT_TYPE_II;
    INT spikes[2] = {0, 15};
    INT freqs[2] = {-7, 8};
    INT poss[2] = {25, 8};
    int i;
    C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *g = (C *)Y(malloc)((size_t)n * sizeof(C));
    for (i = 0; i < 2; i++) {
      R sc = deconv_scale(w, 1, &N, &n, m, &freqs[i]);
      problem *p;
      plan *pln;
      for (ks = 0; ks < N; ks++)
        f_hat[ks] = K(0.0);
      f_hat[spikes[i]] = K(1.0);
      for (ks = 0; ks < n; ks++)
        g[ks] = sc;
      p = Y(mkproblem_deconv)(1, &N, &v, &n, m, w, 1, f_hat, g);
      pln = deconv_awake_plan(p);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
      pln->adt->apply(pln, p);
      CU_ASSERT(CABS(g[poss[i]] - sc) < K(1e-9) * sc);
      CU_ASSERT_EQUAL(count_above(g, n, K(1e-9) * sc), (INT)1);
      Y(plan_destroy)
      (pln);
      Y(problem_destroy)
      (p);
    }
    Y(free)
    (f_hat);
    Y(free)
    (g);
  }

  /* (c) n < N is declined, not planned: the zero-pad length would wrap. */
  {
    INT N = 32, n = 16;
    C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *g = (C *)Y(malloc)((size_t)n * sizeof(C));
    problem *p = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, 1, f_hat, g);
    CU_ASSERT_PTR_NULL(Y(planner_mkplan)(Y(the_planner)(), p));
    Y(problem_destroy)
    (p);
    Y(free)
    (f_hat);
    Y(free)
    (g);
  }

  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}
```

Declare it in `tests/nfast.h`, immediately after the existing
`void Y(check_nfast_deconv_solver)(void);`:

```c
/* rank-1 DECONV across odd N, type-II and the n < N decline */
void Y(check_nfast_deconv_1d_general)(void);
```

Register it in `tests/check_ng.c`, directly after the `"deconv_solver"` line in the
`nfast_suite` block:

```c
  CU_add_test(nfast_suite, "deconv_1d_general", Y(check_nfast_deconv_1d_general));
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
grep -B6 "<FAILURE" CUnitAutomated_ng-Results.xml | head -40
```

Expected: `exit=1`. The `deconv_1d_general` test fails on case (a) -- the odd-`N`
zero-pad leaves grid cell 24 holding the pre-dirty value, so `count_above` returns 2
instead of 1 -- and on case (c), where `planner_mkplan` returns a plan instead of NULL.

- [ ] **Step 3: Rewrite `kernel/deconv/deconv-1d.c`**

Replace the plan struct (lines 30-38 of the current file):

```c
typedef struct
{
  plan super;
  INT n, N;       /* geometry captured at mkplan */
  INT Nneg, Npos; /* slot split: k(ks) = ks - Nneg, Npos = N - Nneg */
  int m, window;
  R *phi_hut_inv; /* length N: 1/phi_hut(ks - Nneg); at awake */
  int precomputed;
} deconv_plan;
```

Replace `deconv_awake`:

```c
/* precompute 1/phi_hut */
static void deconv_awake(plan *ego_, int wakefulness) {
  deconv_plan *pln = (deconv_plan *)ego_;
  if (wakefulness >= PLNR_AWAKE_ZERO) {
    if (!pln->precomputed) {
      INT N = pln->N, ks;
      /* Slot ks carries frequency ks - Nneg, so the table starts at -Nneg. */
      Y(window_phi_hut_apply)
      (pln->window, pln->n, N, pln->m, -pln->Nneg, pln->phi_hut_inv, N);
      for (ks = 0; ks < N; ks++) /* in-place invert */
        pln->phi_hut_inv[ks] = K(1.0) / pln->phi_hut_inv[ks];
      pln->precomputed = 1;
    }
  } else
    pln->precomputed = 0;
}
```

Replace `deconv_apply` and `deconv_apply_adjoint`:

```c
/* apply the real diagonal scale-and-pad map (f_hat -> g). */
static void deconv_apply(const plan *ego_, const problem *p) {
  const deconv_plan *pln = (const deconv_plan *)ego_;
  const problem_deconv *pd = (const problem_deconv *)p;
  INT N = pln->N, n = pln->n, Nneg = pln->Nneg, Npos = pln->Npos;
  const C *f_hat = pd->f_hat;
  C *g = pd->g;
  INT ks;
  /* The scatter below writes every touched cell exactly once, so they don't need
   * pre-zeroing. The untouched cells are the single contiguous run
   * [Npos .. n - Nneg), of length n - N. */
  memset(g + Npos, 0, (size_t)(n - N) * sizeof(C));
  for (ks = 0; ks < N; ks++) {
    INT pos = (ks < Nneg) ? n - Nneg + ks : ks - Nneg;
    g[pos] = f_hat[ks] * pln->phi_hut_inv[ks]; /* / phi_hut(ks - Nneg) */
  }
}

/* apply the adjoint real diagonal scale-and-pad map (g -> f_hat). The adjoint only swaps scatter->gather. */
static void deconv_apply_adjoint(const plan *ego_, const problem *p) {
  const deconv_plan *pln = (const deconv_plan *)ego_;
  const problem_deconv *pd = (const problem_deconv *)p;
  INT N = pln->N, n = pln->n, Nneg = pln->Nneg;
  const C *g = pd->g;
  C *f_hat = pd->f_hat;
  INT ks;
  for (ks = 0; ks < N; ks++) {
    INT pos = (ks < Nneg) ? n - Nneg + ks : ks - Nneg;
    f_hat[ks] = g[pos] * pln->phi_hut_inv[ks]; /* D^H: same 1/phi_hut as forward */
  }
}
```

In `mkplan_deconv_1d`, add the decline after the window check and set the split
instead of `type_ii`:

```c
  if (Y(problem_deconv_n)(p, 0) < Y(problem_deconv_N)(p, 0))
    return 0; /* n < N aliases grid cells and wraps the zero-pad length */

  pln = (deconv_plan *)Y(plan_create)(sizeof(deconv_plan), &deconv_plan_adt);
  pln->n = Y(problem_deconv_ntot)(p);
  pln->N = Y(problem_deconv_Ntot)(p);
  pln->m = pd->m;
  pln->window = pd->window;
  /* Odd N normalizes to type-I in mkproblem_deconv, so type-II implies even N. */
  pln->Nneg = pln->N / 2 -
              (Y(problem_deconv_variant)(p, 0) == NFFT_NDFT_TYPE_II ? (INT)1 : (INT)0);
  pln->Npos = pln->N - pln->Nneg;
  pln->phi_hut_inv = (R *)Y(malloc)((size_t)pln->N * sizeof(R));
  pln->precomputed = 0;
  pln->super.pcost = 2.0 * (double)pln->N;
  return &pln->super;
```

Update the file header comment: the table is `1/phi_hut(k)` for `k = ks - Nneg`, not
`ks - N/2(+type_ii)`.

- [ ] **Step 4: Run the test to verify it passes**

```bash
clang-format -i kernel/deconv/deconv-1d.c tests/nfast.c
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
```

Expected: `exit=0`. `check_nfast_deconv_solver` (the pre-existing rank-1 test,
including its type-II spike) must still pass unchanged.

- [ ] **Step 5: Commit**

```bash
git add kernel/deconv/deconv-1d.c tests/nfast.c tests/nfast.h tests/check_ng.c
git commit -m "fix(deconv): correct odd-N zero-pad and decline n < N in rank-1 DECONV"
```

---

### Task 2: rank-2 DECONV -- block runs, variant-aware

**Files:**
- Modify: `kernel/deconv/deconv-2d.c`
- Test: `tests/nfast.c` (new function), `tests/nfast.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: `deconv_awake_plan`, `deconv_scale`, `count_above` from Task 1.
- Produces: `void Y(check_nfast_deconv_2d_general)(void)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/nfast.c`, after `Y(check_nfast_deconv_1d_general)`:

```c
/* Rank-2 DECONV with a type-II axis and an odd axis. The rank >= 2 solvers
 * hash `variant` into their wisdom key but ignored it in the index map, so a
 * type-II problem got a type-I plan. Slots are checked as single spikes, so a
 * one-cell index slip shows up as a miss plus a stray. */
void Y(check_nfast_deconv_2d_general)(void) {
  const int m = 6, w = Y(get_window_id)();
  INT N[2] = {16, 15}, n[2] = {32, 32};
  int variant[2] = {NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I};
  /* axis 0: Nneg = 7, Npos = 9.  axis 1: Nneg = 7, Npos = 8. */
  INT Ntot = N[0] * N[1], ntot = n[0] * n[1], j;
  C *f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *g = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  INT slots[3][2] = {{0, 0}, {7, 7}, {15, 14}};
  INT freqs[3][2] = {{-7, -7}, {0, 0}, {8, 7}};
  INT cells[3][2] = {{25, 25}, {0, 0}, {8, 7}};
  int i;

  Y(deconv_ensure_registered)
  ();

  for (i = 0; i < 3; i++) {
    R sc = deconv_scale(w, 2, N, n, m, freqs[i]);
    INT fs = slots[i][0] * N[1] + slots[i][1];
    INT gs = cells[i][0] * n[1] + cells[i][1];
    problem *p;
    plan *pln;
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    f_hat[fs] = K(1.0);
    for (j = 0; j < ntot; j++)
      g[j] = sc;
    p = Y(mkproblem_deconv)(2, N, variant, n, m, w, 1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply(pln, p);
    CU_ASSERT(CABS(g[gs] - sc) < K(1e-9) * sc);
    CU_ASSERT_EQUAL(count_above(g, ntot, K(1e-9) * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  /* adjoint: the corner cell gathers into slot (0,0) only. */
  {
    R sc = deconv_scale(w, 2, N, n, m, freqs[0]);
    problem *p;
    plan *pln;
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    for (j = 0; j < ntot; j++)
      g[j] = K(0.0);
    g[cells[0][0] * n[1] + cells[0][1]] = K(1.0);
    p = Y(mkproblem_deconv)(2, N, variant, n, m, w, -1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply_adjoint(pln, p);
    CU_ASSERT(CABS(f_hat[0] - sc) < K(1e-9) * sc);
    CU_ASSERT_EQUAL(count_above(f_hat, Ntot, K(1e-9) * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  Y(free)
  (f_hat);
  Y(free)
  (g);
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}
```

Declare in `tests/nfast.h` and register in `tests/check_ng.c`:

```c
/* rank-2 DECONV across a type-II axis and an odd axis */
void Y(check_nfast_deconv_2d_general)(void);
```
```c
  CU_add_test(nfast_suite, "deconv_2d_general", Y(check_nfast_deconv_2d_general));
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
grep -B6 "<FAILURE" CUnitAutomated_ng-Results.xml | head -40
```

Expected: `exit=1` in `deconv_2d_general`. In a `--enable-debug` tree it aborts
earlier on `A(pln->N1 % 2 == 0)`, since `N[1] = 15` is odd -- that assert is removed
in Step 3.

- [ ] **Step 3: Rewrite `kernel/deconv/deconv-2d.c`**

Replace the plan struct:

```c
typedef struct
{
  plan super;
  INT N0, N1, n0, n1;         /* geometry captured at mkplan */
  INT Nneg0, Npos0, Nneg1, Npos1; /* per-axis slot split */
  int m, window;
  R *phi_hut_inv0; /* length N0: 1/phi_hut(n0,N0,m, ks0 - Nneg0), at awake */
  R *phi_hut_inv1; /* length N1: 1/phi_hut(n1,N1,m, ks1 - Nneg1), at awake */
  int precomputed;
} deconv_2d_plan;
```

Replace `deconv_2d_awake`:

```c
/* precompute 1/phi_hut */
static void deconv_2d_awake(plan *ego_, int wakefulness) {
  deconv_2d_plan *pln = (deconv_2d_plan *)ego_;
  if (wakefulness >= PLNR_AWAKE_ZERO) {
    if (!pln->precomputed) {
      INT ks;
      Y(window_phi_hut_apply)
      (pln->window, pln->n0, pln->N0, pln->m, -pln->Nneg0, pln->phi_hut_inv0,
       pln->N0);
      for (ks = 0; ks < pln->N0; ks++)
        pln->phi_hut_inv0[ks] = K(1.0) / pln->phi_hut_inv0[ks];
      Y(window_phi_hut_apply)
      (pln->window, pln->n1, pln->N1, pln->m, -pln->Nneg1, pln->phi_hut_inv1,
       pln->N1);
      for (ks = 0; ks < pln->N1; ks++)
        pln->phi_hut_inv1[ks] = K(1.0) / pln->phi_hut_inv1[ks];
      pln->precomputed = 1;
    }
  } else
    pln->precomputed = 0;
}
```

Replace both `deconv_2d_apply` and `deconv_2d_apply_adjoint` with one folded runner
plus two thin entry points:

```c
/* Each axis splits into two contiguous runs of slots: the negative half maps to
 * the grid tail, the non-negative half to the grid head. Even type-I makes the
 * two runs equal (Nneg == Npos == N/2), which is the old quadrant unroll. */
static void deconv_2d_run(const deconv_2d_plan *pln, const problem_deconv *pd,
                          int forward) {
  const INT N1 = pln->N1, n0 = pln->n0, n1 = pln->n1;
  INT len0[2], sof0[2], gof0[2]; /* run length, f_hat slot offset, grid offset */
  INT len1[2], sof1[2], gof1[2];
  C *f_hat = pd->f_hat;
  C *g_hat = pd->g;
  int a, b;

  len0[0] = pln->Nneg0;
  sof0[0] = 0;
  gof0[0] = n0 - pln->Nneg0;
  len0[1] = pln->Npos0;
  sof0[1] = pln->Nneg0;
  gof0[1] = 0;
  len1[0] = pln->Nneg1;
  sof1[0] = 0;
  gof1[0] = n1 - pln->Nneg1;
  len1[1] = pln->Npos1;
  sof1[1] = pln->Nneg1;
  gof1[1] = 0;

  /* In nD the touched frequencies form a box per grid corner, so the zero-pad is
   * fragmented; clearing the whole grid is typically more efficient. The adjoint
   * writes every f_hat slot, so it needs no clear. */
  if (forward)
    memset(g_hat, 0, (size_t)(n0 * n1) * sizeof(C));

  for (a = 0; a < 2; a++)
    for (b = 0; b < 2; b++) {
      const R *p0 = pln->phi_hut_inv0 + sof0[a];
      const R *p1 = pln->phi_hut_inv1 + sof1[b];
      INT i, j;
      for (i = 0; i < len0[a]; i++) {
        R ck0 = p0[i];
        C *gh = g_hat + (gof0[a] + i) * n1 + gof1[b];
        C *fh = f_hat + (sof0[a] + i) * N1 + sof1[b];
        /* Group the two real factors first: complex-by-real costs two real
         * multiplies, so (C * R) * R is four and C * (R * R) is three. */
        if (forward)
          for (j = 0; j < len1[b]; j++)
            gh[j] = fh[j] * (ck0 * p1[j]);
        else
          for (j = 0; j < len1[b]; j++)
            fh[j] = gh[j] * (ck0 * p1[j]);
      }
    }
}

/* apply the real diagonal scale-and-pad map (f_hat -> g). */
static void deconv_2d_apply(const plan *ego_, const problem *p) {
  deconv_2d_run((const deconv_2d_plan *)ego_, (const problem_deconv *)p, 1);
}

/* apply the adjoint real diagonal scale-and-pad map (g -> f_hat). The adjoint only swaps scatter->gather. */
static void deconv_2d_apply_adjoint(const plan *ego_, const problem *p) {
  deconv_2d_run((const deconv_2d_plan *)ego_, (const problem_deconv *)p, 0);
}
```

In `mkplan_deconv_2d`, delete both `A(N % 2 == 0)` asserts and the comment above them,
add the `n < N` decline, and set the split:

```c
  if (Y(problem_deconv_n)(p, 0) < Y(problem_deconv_N)(p, 0) ||
      Y(problem_deconv_n)(p, 1) < Y(problem_deconv_N)(p, 1))
    return 0; /* n < N aliases grid cells */

  pln = (deconv_2d_plan *)Y(plan_create)(sizeof(deconv_2d_plan), &deconv_2d_plan_adt);
  pln->N0 = Y(problem_deconv_N)(p, 0);
  pln->N1 = Y(problem_deconv_N)(p, 1);
  pln->n0 = Y(problem_deconv_n)(p, 0);
  pln->n1 = Y(problem_deconv_n)(p, 1);
  pln->Nneg0 = pln->N0 / 2 -
               (Y(problem_deconv_variant)(p, 0) == NFFT_NDFT_TYPE_II ? (INT)1 : (INT)0);
  pln->Npos0 = pln->N0 - pln->Nneg0;
  pln->Nneg1 = pln->N1 / 2 -
               (Y(problem_deconv_variant)(p, 1) == NFFT_NDFT_TYPE_II ? (INT)1 : (INT)0);
  pln->Npos1 = pln->N1 - pln->Nneg1;
  pln->m = pd->m;
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
clang-format -i kernel/deconv/deconv-2d.c tests/nfast.c
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
```

Expected: `exit=0`. The pre-existing `native_fast_2d` and `native_fast_2d_adjoint`
suites exercise the even type-I path through the same code and must still pass.

- [ ] **Step 5: Commit**

```bash
git add kernel/deconv/deconv-2d.c tests/nfast.c tests/nfast.h tests/check_ng.c
git commit -m "feat(deconv): per-axis slot split in rank-2 DECONV"
```

---

### Task 3: rank-3 DECONV -- block runs, variant-aware

**Files:**
- Modify: `kernel/deconv/deconv-3d.c`
- Test: `tests/nfast.c` (new function), `tests/nfast.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: `deconv_awake_plan`, `deconv_scale`, `count_above` from Task 1.
- Produces: `void Y(check_nfast_deconv_3d_general)(void)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/nfast.c`:

```c
/* Rank-3 DECONV with two type-II axes and one odd axis. */
void Y(check_nfast_deconv_3d_general)(void) {
  const int m = 6, w = Y(get_window_id)();
  INT N[3] = {16, 15, 10}, n[3] = {32, 32, 20};
  int variant[3] = {NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I, NFFT_NDFT_TYPE_II};
  /* axis 0: Nneg 7 / Npos 9.  axis 1: Nneg 7 / Npos 8.  axis 2: Nneg 4 / Npos 6. */
  INT Ntot = N[0] * N[1] * N[2], ntot = n[0] * n[1] * n[2], j;
  C *f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *g = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  INT slots[3][3] = {{0, 0, 0}, {7, 7, 4}, {15, 14, 9}};
  INT freqs[3][3] = {{-7, -7, -4}, {0, 0, 0}, {8, 7, 5}};
  INT cells[3][3] = {{25, 25, 16}, {0, 0, 0}, {8, 7, 5}};
  int i;

  Y(deconv_ensure_registered)
  ();

  for (i = 0; i < 3; i++) {
    R sc = deconv_scale(w, 3, N, n, m, freqs[i]);
    INT fs = (slots[i][0] * N[1] + slots[i][1]) * N[2] + slots[i][2];
    INT gs = (cells[i][0] * n[1] + cells[i][1]) * n[2] + cells[i][2];
    problem *p;
    plan *pln;
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    f_hat[fs] = K(1.0);
    for (j = 0; j < ntot; j++)
      g[j] = sc;
    p = Y(mkproblem_deconv)(3, N, variant, n, m, w, 1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply(pln, p);
    CU_ASSERT(CABS(g[gs] - sc) < K(1e-9) * sc);
    CU_ASSERT_EQUAL(count_above(g, ntot, K(1e-9) * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  /* adjoint: the corner cell gathers into slot (0,0,0) only. */
  {
    R sc = deconv_scale(w, 3, N, n, m, freqs[0]);
    INT gs = (cells[0][0] * n[1] + cells[0][1]) * n[2] + cells[0][2];
    problem *p;
    plan *pln;
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    for (j = 0; j < ntot; j++)
      g[j] = K(0.0);
    g[gs] = K(1.0);
    p = Y(mkproblem_deconv)(3, N, variant, n, m, w, -1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply_adjoint(pln, p);
    CU_ASSERT(CABS(f_hat[0] - sc) < K(1e-9) * sc);
    CU_ASSERT_EQUAL(count_above(f_hat, Ntot, K(1e-9) * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  Y(free)
  (f_hat);
  Y(free)
  (g);
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}
```

Declare in `tests/nfast.h` and register in `tests/check_ng.c`:

```c
/* rank-3 DECONV across type-II and odd axes */
void Y(check_nfast_deconv_3d_general)(void);
```
```c
  CU_add_test(nfast_suite, "deconv_3d_general", Y(check_nfast_deconv_3d_general));
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
grep -B6 "<FAILURE" CUnitAutomated_ng-Results.xml | head -40
```

Expected: `exit=1` in `deconv_3d_general`.

- [ ] **Step 3: Rewrite `kernel/deconv/deconv-3d.c`**

Same shape as Task 2, one axis deeper. Plan struct:

```c
typedef struct
{
  plan super;
  INT N0, N1, N2, n0, n1, n2; /* geometry captured at mkplan */
  INT Nneg0, Npos0, Nneg1, Npos1, Nneg2, Npos2; /* per-axis slot split */
  int m, window;
  R *phi_hut_inv0; /* length N0: 1/phi_hut(n0,N0,m, ks0 - Nneg0), at awake */
  R *phi_hut_inv1; /* length N1: 1/phi_hut(n1,N1,m, ks1 - Nneg1), at awake */
  R *phi_hut_inv2; /* length N2: 1/phi_hut(n2,N2,m, ks2 - Nneg2), at awake */
  int precomputed;
} deconv_3d_plan;
```

`deconv_3d_awake` passes `-pln->Nneg0`, `-pln->Nneg1`, `-pln->Nneg2` as the three `k0`
arguments in place of `-N/2`.

Replace both apply bodies with one folded runner:

```c
/* Each axis splits into two contiguous runs of slots: the negative half maps to
 * the grid tail, the non-negative half to the grid head. Even type-I makes the
 * two runs equal (Nneg == Npos == N/2), which is the old octant unroll. */
static void deconv_3d_run(const deconv_3d_plan *pln, const problem_deconv *pd,
                          int forward) {
  const INT N1 = pln->N1, N2 = pln->N2;
  const INT n0 = pln->n0, n1 = pln->n1, n2 = pln->n2;
  INT len0[2], sof0[2], gof0[2];
  INT len1[2], sof1[2], gof1[2];
  INT len2[2], sof2[2], gof2[2];
  C *f_hat = pd->f_hat;
  C *g_hat = pd->g;
  int a, b, c;

  len0[0] = pln->Nneg0;
  sof0[0] = 0;
  gof0[0] = n0 - pln->Nneg0;
  len0[1] = pln->Npos0;
  sof0[1] = pln->Nneg0;
  gof0[1] = 0;
  len1[0] = pln->Nneg1;
  sof1[0] = 0;
  gof1[0] = n1 - pln->Nneg1;
  len1[1] = pln->Npos1;
  sof1[1] = pln->Nneg1;
  gof1[1] = 0;
  len2[0] = pln->Nneg2;
  sof2[0] = 0;
  gof2[0] = n2 - pln->Nneg2;
  len2[1] = pln->Npos2;
  sof2[1] = pln->Nneg2;
  gof2[1] = 0;

  if (forward)
    memset(g_hat, 0, (size_t)(n0 * n1 * n2) * sizeof(C));

  for (a = 0; a < 2; a++)
    for (b = 0; b < 2; b++)
      for (c = 0; c < 2; c++) {
        const R *p0 = pln->phi_hut_inv0 + sof0[a];
        const R *p1 = pln->phi_hut_inv1 + sof1[b];
        const R *p2 = pln->phi_hut_inv2 + sof2[c];
        INT i, j, k;
        for (i = 0; i < len0[a]; i++) {
          R ck0 = p0[i];
          for (j = 0; j < len1[b]; j++) {
            R ck01 = ck0 * p1[j];
            C *gh = g_hat + ((gof0[a] + i) * n1 + gof1[b] + j) * n2 + gof2[c];
            C *fh = f_hat + ((sof0[a] + i) * N1 + sof1[b] + j) * N2 + sof2[c];
            /* Group the real factors first: three real multiplies per element
             * against six for the left-associated ck01 * ck11 * ck21 chain. */
            if (forward)
              for (k = 0; k < len2[c]; k++)
                gh[k] = fh[k] * (ck01 * p2[k]);
            else
              for (k = 0; k < len2[c]; k++)
                fh[k] = gh[k] * (ck01 * p2[k]);
          }
        }
      }
}

/* apply the real diagonal scale-and-pad map (f_hat -> g). */
static void deconv_3d_apply(const plan *ego_, const problem *p) {
  deconv_3d_run((const deconv_3d_plan *)ego_, (const problem_deconv *)p, 1);
}

/* apply the adjoint real diagonal scale-and-pad map (g -> f_hat). The adjoint only swaps scatter->gather. */
static void deconv_3d_apply_adjoint(const plan *ego_, const problem *p) {
  deconv_3d_run((const deconv_3d_plan *)ego_, (const problem_deconv *)p, 0);
}
```

In `mkplan_deconv_3d`, delete the three `A(N % 2 == 0)` asserts and their comment, add
the decline, and set the split:

```c
  {
    int t;
    for (t = 0; t < 3; t++)
      if (Y(problem_deconv_n)(p, t) < Y(problem_deconv_N)(p, t))
        return 0; /* n < N aliases grid cells */
  }
```
```c
  pln->Nneg0 = pln->N0 / 2 -
               (Y(problem_deconv_variant)(p, 0) == NFFT_NDFT_TYPE_II ? (INT)1 : (INT)0);
  pln->Npos0 = pln->N0 - pln->Nneg0;
  pln->Nneg1 = pln->N1 / 2 -
               (Y(problem_deconv_variant)(p, 1) == NFFT_NDFT_TYPE_II ? (INT)1 : (INT)0);
  pln->Npos1 = pln->N1 - pln->Nneg1;
  pln->Nneg2 = pln->N2 / 2 -
               (Y(problem_deconv_variant)(p, 2) == NFFT_NDFT_TYPE_II ? (INT)1 : (INT)0);
  pln->Npos2 = pln->N2 - pln->Nneg2;
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
clang-format -i kernel/deconv/deconv-3d.c tests/nfast.c
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
```

Expected: `exit=0`. The pre-existing 3D fast suites cover the even type-I path.

- [ ] **Step 5: Commit**

```bash
git add kernel/deconv/deconv-3d.c tests/nfast.c tests/nfast.h tests/check_ng.c
git commit -m "feat(deconv): per-axis slot split in rank-3 DECONV"
```

---

### Task 4: rank >= 4 DECONV -- generalized carry odometer

**Files:**
- Modify: `kernel/deconv/deconv-nd.c`
- Test: `tests/nfast.c` (new function), `tests/nfast.h`, `tests/check_ng.c`

**Interfaces:**
- Consumes: `deconv_awake_plan`, `deconv_scale`, `count_above` from Task 1.
- Produces: `void Y(check_nfast_deconv_nd_general)(void)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/nfast.c`:

```c
/* Rank-4 DECONV with mixed type-II and odd axes, driving the carry odometer. */
void Y(check_nfast_deconv_nd_general)(void) {
  const int m = 3, w = Y(get_window_id)();
  INT N[4] = {8, 7, 6, 5}, n[4] = {16, 16, 16, 16};
  int variant[4] = {NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I, NFFT_NDFT_TYPE_II,
                    NFFT_NDFT_TYPE_I};
  /* Nneg / Npos per axis: 3/5, 3/4, 2/4, 2/3. */
  INT Ntot = N[0] * N[1] * N[2] * N[3];
  INT ntot = n[0] * n[1] * n[2] * n[3], j;
  C *f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *g = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  INT slots[3][4] = {{0, 0, 0, 0}, {3, 3, 2, 2}, {7, 6, 5, 4}};
  INT freqs[3][4] = {{-3, -3, -2, -2}, {0, 0, 0, 0}, {4, 3, 3, 2}};
  INT cells[3][4] = {{13, 13, 14, 14}, {0, 0, 0, 0}, {4, 3, 3, 2}};
  int i, t;

  Y(deconv_ensure_registered)
  ();

  for (i = 0; i < 3; i++) {
    R sc = deconv_scale(w, 4, N, n, m, freqs[i]);
    INT fs = 0, gs = 0;
    problem *p;
    plan *pln;
    for (t = 0; t < 4; t++) {
      fs = fs * N[t] + slots[i][t];
      gs = gs * n[t] + cells[i][t];
    }
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    f_hat[fs] = K(1.0);
    for (j = 0; j < ntot; j++)
      g[j] = sc;
    p = Y(mkproblem_deconv)(4, N, variant, n, m, w, 1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply(pln, p);
    CU_ASSERT(CABS(g[gs] - sc) < K(1e-9) * sc);
    CU_ASSERT_EQUAL(count_above(g, ntot, K(1e-9) * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  /* adjoint: the corner cell gathers into slot (0,0,0,0) only. */
  {
    R sc = deconv_scale(w, 4, N, n, m, freqs[0]);
    INT gs = 0;
    problem *p;
    plan *pln;
    for (t = 0; t < 4; t++)
      gs = gs * n[t] + cells[0][t];
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    for (j = 0; j < ntot; j++)
      g[j] = K(0.0);
    g[gs] = K(1.0);
    p = Y(mkproblem_deconv)(4, N, variant, n, m, w, -1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply_adjoint(pln, p);
    CU_ASSERT(CABS(f_hat[0] - sc) < K(1e-9) * sc);
    CU_ASSERT_EQUAL(count_above(f_hat, Ntot, K(1e-9) * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  Y(free)
  (f_hat);
  Y(free)
  (g);
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}
```

Declare in `tests/nfast.h` and register in `tests/check_ng.c`:

```c
/* rank-4 DECONV across mixed type-II and odd axes */
void Y(check_nfast_deconv_nd_general)(void);
```
```c
  CU_add_test(nfast_suite, "deconv_nd_general", Y(check_nfast_deconv_nd_general));
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
grep -B6 "<FAILURE" CUnitAutomated_ng-Results.xml | head -40
```

Expected: `exit=1` in `deconv_nd_general`.

- [ ] **Step 3: Rewrite the geometry in `kernel/deconv/deconv-nd.c`**

Add the split array to the plan struct:

```c
typedef struct
{
  plan super;
  int d;          /* rank, >= 4 */
  INT *N, *n;     /* owned, length d: geometry captured at mkplan */
  INT *Nneg;      /* owned, length d: slot split, k(ks) = ks - Nneg[t] */
  INT Ntot, ntot; /* owned products, captured at mkplan */
  int m, window;
  R **phi_hut_inv; /* owned array of d owned tables; phi_hut_inv[t] has
                    * length N[t]: 1/phi_hut(n[t],N[t],m,ks-Nneg[t]), at awake */
  int precomputed;
} deconv_nd_plan;
```

In `deconv_nd_awake`, start each table at `-pln->Nneg[t]`:

```c
        Y(window_phi_hut_apply)
        (pln->window, pln->n[t], Nt, pln->m, -pln->Nneg[t], pln->phi_hut_inv[t],
         Nt);
```

In `deconv_nd_run`, bind the split and substitute it in the three odometer sites. The
carry condition `kp[t] == N[t] - 1` is unchanged -- `kp` still counts `0 .. N-1`.

```c
  const INT *Nneg = pln->Nneg;
```

Initialization:

```c
  for (t = d - 1; 0 <= t; t--) {
    kp[t] = k[t] = 0;
    ks[t] = Nneg[t];
  }
```

Carry-reset inside the count loop:

```c
    for (t = d - 1; (t > 0) && (kp[t] == N[t] - 1); t--) {
      kp[t] = k[t] = 0;
      ks[t] = Nneg[t];
    }
    kp[t]++;
    k[t]++;
    ks[t]++;
    /* The non-negative run has Npos = N - Nneg slots; after it, wrap to the
     * grid tail and restart at f_hat slot 0. */
    if (kp[t] == N[t] - Nneg[t]) {
      k[t] = n[t] - Nneg[t];
      ks[t] = 0;
    }
```

In `mkplan_deconv_nd`, allocate `Nneg`, drop the parity assert, and add the decline:

```c
  pln->N = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->n = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->Nneg = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->phi_hut_inv = (R **)Y(malloc)((size_t)d * sizeof(R *));
  pln->Ntot = Y(problem_deconv_Ntot)(p);
  pln->ntot = Y(problem_deconv_ntot)(p);
  for (t = 0; t < d; t++) {
    pln->N[t] = Y(problem_deconv_N)(p, t);
    pln->n[t] = Y(problem_deconv_n)(p, t);
    pln->Nneg[t] =
        pln->N[t] / 2 -
        (Y(problem_deconv_variant)(p, t) == NFFT_NDFT_TYPE_II ? (INT)1 : (INT)0);
    pln->phi_hut_inv[t] = (R *)Y(malloc)((size_t)pln->N[t] * sizeof(R));
  }
```

The `n < N` decline must run before `plan_create`, so it reads from the problem:

```c
  d = pd->sz->rnk;
  for (t = 0; t < d; t++)
    if (Y(problem_deconv_n)(p, t) < Y(problem_deconv_N)(p, t))
      return 0; /* n < N aliases grid cells */
```

Free `Nneg` in `deconv_nd_destroy`, next to `pln->n` and `pln->N`:

```c
  Y(free)
  (pln->Nneg);
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
clang-format -i kernel/deconv/deconv-nd.c tests/nfast.c
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
```

Expected: `exit=0`. The pre-existing rank >= 4 fast suite covers even type-I.

- [ ] **Step 5: Commit**

```bash
git add kernel/deconv/deconv-nd.c tests/nfast.c tests/nfast.h tests/check_ng.c
git commit -m "feat(deconv): per-axis slot split in the rank >= 4 DECONV odometer"
```

---

### Task 5: flip `guards_ok` -- accept odd and type-II, require sigma > 1

**Files:**
- Modify: `kernel/nfft/nfft-nd.c:33-56` (`guards_ok` and its call site)
- Modify: `kernel/deconv/deconv-2d.c`, `kernel/deconv/deconv-3d.c` (stale comments)
- Test: `tests/nplan.c` (`Y(check_nplan_guard_declines)`)

**Interfaces:**
- Consumes: `Y(problem_nfft_N)(p, t)`, `Y(problem_nfft_n)(p, t)`,
  `Y(mkproblem_nfft)(d, N, variant, n, M, m, window, sign, fftw_flags, x, copy_x, f_hat, f)`,
  `expect_winner(planner *, problem *, const char *)` in `tests/nplan.c`.
- Produces: no new symbol. Behaviour change: `nfft_solver_fast_native` now applies to
  odd `N` and type-II axes, and declines any axis with `n <= N`.

- [ ] **Step 1: Rewrite the test to state the new contract**

Replace `Y(check_nplan_guard_declines)` in `tests/nplan.c` (currently lines 344-395)
and its leading comment, in full:

```c
/* The fast solver serves odd N and type-II. It declines any axis with
 * sigma = n/N <= 1, where the deconvolution has no zero-pad to write into and
 * sinc-power's phi_hut vanishes at the band edge; a direct native takes those
 * (cheapest is blocked for 1d), and with NO_DIRECT too the guru returns NULL. */
void Y(check_nplan_guard_declines)(void) {
  planner *pl;
  problem *p;
  Y(nfft_ensure_registered)
  ();
  pl = Y(the_planner)();
  Y(planner_forget)
  (pl, PLNR_FORGET_ALL);
  {
    R x = K(0.1);
    C fh = K(0.0), f = K(0.0);
    INT N = 63, n = 128; /* odd, big M */
    /* x is a scalar, not an M-element array: copy_x=0 so mkproblem_nfft never
     * reads past it. */
    p = Y(mkproblem_nfft)(1, &N, 0, &n, 200000, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, &x, 0, &fh, &f);
    expect_winner(pl, p, "nfft_solver_fast_native");
    Y(problem_destroy)
    (p);
  }
  {
    R x = K(0.1);
    C fh = K(0.0), f = K(0.0);
    INT N = 64, n = 128;
    int v = NFFT_NDFT_TYPE_II; /* type-II, big M */
    p = Y(mkproblem_nfft)(1, &N, &v, &n, 200000, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, &x, 0, &fh, &f);
    expect_winner(pl, p, "nfft_solver_fast_native");
    Y(problem_destroy)
    (p);
  }
  {
    R x = K(0.1);
    C fh = K(0.0), f = K(0.0);
    INT N = 64, n = 64; /* sigma == 1: declined even at big M */
    p = Y(mkproblem_nfft)(1, &N, 0, &n, 200000, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, &x, 0, &fh, &f);
    expect_winner(pl, p, "nfft_solver_ndft_1d_blocked");
    Y(problem_destroy)
    (p);
  }
  /* NO_DIRECT + sigma == 1 -> no solver -> guru returns NULL */
  {
    INT N = 64, n = 64, M = 1000;
    R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
    C *fh = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
    INT j;
    for (j = 0; j < M; j++)
      x[j] = K(0.0);
    CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, fh, f, 0u,
                                       NFFT_ESTIMATE | NFFT_NO_DIRECT));
    Y(free)
    (x);
    Y(free)
    (fh);
    Y(free)
    (f);
  }
  /* NO_DIRECT + odd now plans: the fast serves it. */
  {
    INT N = 63, n = 128, M = 1000;
    R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
    C *fh = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
    Y(plan_ng) * pg;
    INT j;
    for (j = 0; j < M; j++)
      x[j] = K(0.0);
    pg = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, fh, f, 0u,
                         NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL(pg);
    Y(plan_ng_destroy)
    (pg);
    Y(free)
    (x);
    Y(free)
    (fh);
    Y(free)
    (f);
  }
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
grep -B6 "<FAILURE" CUnitAutomated_ng-Results.xml | head -40
```

Expected: `exit=1` in `guard_declines`. The two `expect_winner(..., "nfft_solver_fast_native")`
cases still see `nfft_solver_ndft_1d_blocked`, and the NO_DIRECT-plus-odd case still
gets NULL.

- [ ] **Step 3: Rewrite `guards_ok` in `kernel/nfft/nfft-nd.c`**

Replace the function and its leading comment (lines 33-56):

```c
/* Applicability guard: every axis of sz must satisfy N_t > m, n_t > 2m + 2 and
 * n_t > N_t. The last one is oversampling sigma = n/N strictly above 1; at
 * sigma == 1 there is no zero-pad to deconvolve into and sinc-power's phi_hut
 * is exactly zero at the band edge. Odd N and type-II are both served: the
 * DECONV children carry the per-axis slot split. A surviving unit axis needs no
 * separate decline, since N_t == 1 fails N_t > m for every m >= 1. */
static int guards_ok(const problem *p, int m) {
  const problem_nfft *ego = (const problem_nfft *)p;
  int t;
  for (t = 0; t < ego->sz->rnk; ++t) {
    INT Nt = Y(problem_nfft_N)(p, t);
    INT nt = Y(problem_nfft_n)(p, t);
    if (!(Nt > (INT)m))
      return 0;
    if (!(nt > (INT)(2 * m + 2)))
      return 0;
    if (!(nt > Nt)) /* sigma = n/N > 1 */
      return 0;
  }
  return 1;
}
```

At the call site in `mkplan_native_fast`, replace the trailing comment:

```c
  if (!guards_ok(p, pn->m))
    return 0; /* declines undersized axes and sigma <= 1 */
```

Sweep the comments that name a decline this guard no longer makes, or a kernel file
that no longer exists -- `kernel/nfft/nfast.c` was renamed to `kernel/nfft/nfft-nd.c`.
Find them with:

```bash
grep -rn "nfast\.c" kernel/nfft/*.c kernel/deconv/*.c kernel/conv/*.c tests/*.c
```

Expected hits and their fixes:

- `kernel/nfft/conf.c:51` -- `mkplan_native_fast (nfast.c)` becomes
  `mkplan_native_fast (nfft-nd.c)`.
- `kernel/deconv/deconv-2d.c:195`, `kernel/deconv/deconv-3d.c:241` -- these are the
  comment blocks above the parity asserts. Tasks 2 and 3 delete both, so nothing should
  remain; if either survived, delete it now.
- `tests/nfast.c:1127` -- Task 6 replaces this comment.
- `tests/nplan.c` -- any `guards_ok (kernel/nfft/nfast.c)` becomes
  `guards_ok (kernel/nfft/nfft-nd.c)`.

Also correct the unit-axis comment in `Y(check_nplan_elides_unit_axes)`, which claims a
decline that no longer exists:

```c
  /* A surviving unit axis fails guards_ok (N_t > m), so under NFFT_NO_DIRECT a
   * non-NULL plan proves the middle axis was elided. */
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
clang-format -i kernel/nfft/nfft-nd.c tests/nplan.c
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
```

Expected: `exit=0`. Every other suite must still pass -- this is the step where odd and
type-II problems first reach the DECONV code that Tasks 1-4 generalized, so a failure
elsewhere means one of those tasks is incomplete. Suites to watch:
`nplan/odd_n`, `nplan/type_ii_1d`, `nplan/type_ii_nd`, `nplan_data/native_reference`.

- [ ] **Step 5: Commit**

```bash
git add kernel/nfft/nfft-nd.c kernel/deconv/deconv-2d.c kernel/deconv/deconv-3d.c tests/nplan.c
git commit -m "feat(nfft): fast native accepts odd N and type-II, requires sigma > 1"
```

---

### Task 6: 2D type-II reference cases through the fast path

**Files:**
- Modify: `tests/nfast.c:1124-1148` (the `files_2d` / `adjoint_files_2d` rosters) and
  the two loop bodies that read them

**Interfaces:**
- Consumes: existing `read_nd_case`, `rel_max_errC`, `NFAST_LEGACY_REL_BOUND`.
- Produces: file-scope type `nfast_2d_case` used by both 2D loops.

The reference files `tests/data/nfft_2d_10_20_50_t210.txt`,
`nfft_2d_20_10_50_t210.txt` and their `nfft_adjoint_*` twins already exist -- generated
from `tests/refgen/grids.py`, variant `[TYPE_II, TYPE_I]`. No data regeneration.

- [ ] **Step 1: Extend the roster and pass the variant**

Replace the roster block and its comment:

```c
/* The 2D reference cases, the same geometries tests/nfft.c's 2D suite runs,
 * plus the _t210 type-II variants the composed solver now serves. */
typedef struct
{
  const char *file;
  const int *variant; /* NULL = all type-I */
} nfast_2d_case;

static const int nfast_v_ii_i[2] = {NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I};

static const nfast_2d_case files_2d[] = {
    {"data/nfft_2d_10_10_20.txt", 0},
    {"data/nfft_2d_10_10_50.txt", 0},
    {"data/nfft_2d_10_20_20.txt", 0},
    {"data/nfft_2d_10_20_50.txt", 0},
    {"data/nfft_2d_20_10_20.txt", 0},
    {"data/nfft_2d_20_10_50.txt", 0},
    {"data/nfft_2d_20_20_20.txt", 0},
    {"data/nfft_2d_20_20_50.txt", 0},
    {"data/nfft_2d_10_20_50_t210.txt", nfast_v_ii_i},
    {"data/nfft_2d_20_10_50_t210.txt", nfast_v_ii_i},
};
static const nfast_2d_case adjoint_files_2d[] = {
    {"data/nfft_adjoint_2d_10_10_20.txt", 0},
    {"data/nfft_adjoint_2d_10_10_50.txt", 0},
    {"data/nfft_adjoint_2d_10_20_20.txt", 0},
    {"data/nfft_adjoint_2d_10_20_50.txt", 0},
    {"data/nfft_adjoint_2d_20_10_20.txt", 0},
    {"data/nfft_adjoint_2d_20_10_50.txt", 0},
    {"data/nfft_adjoint_2d_20_20_20.txt", 0},
    {"data/nfft_adjoint_2d_20_20_50.txt", 0},
    {"data/nfft_adjoint_2d_10_20_50_t210.txt", nfast_v_ii_i},
    {"data/nfft_adjoint_2d_20_10_50_t210.txt", nfast_v_ii_i},
};
#define NFAST_2D_NFILES (sizeof(files_2d) / sizeof(files_2d[0]))
```

In `Y(check_nfast_native_fast_2d)`, read through the struct and pass the variant:

```c
    CU_ASSERT_TRUE_FATAL(
        read_nd_case(files_2d[fi].file, &d, &N, &M, &x, &f_hat, &f));
```
```c
      p = Y(plan_ng_guru)(2, N, files_2d[fi].variant, n, M, m,
                          Y(get_window_id)(), x, f_hat, got, 0u,
                          NFFT_ESTIMATE | NFFT_NO_DIRECT);
```

The in-test legacy cross-check calls `NFFT(init_guru)`, which has no type-II concept.
Its block currently opens with the comment line

```c
      /* in-test legacy reference: the same fast algorithm through the old
       * X(plan) API. */
      {
        NFFT(plan)
        lp;
```

and closes with

```c
        NFFT(finalize)
        (&lp);
      }
```

Change only the opening brace line, leaving every statement between the two anchors
byte for byte as it is:

```c
      /* in-test legacy reference: the same fast algorithm through the old
       * X(plan) API, which has no type-II concept -- type-I cases only. */
      if (!files_2d[fi].variant) {
        NFFT(plan)
        lp;
```

Apply the same three edits in `Y(check_nfast_native_fast_2d_adjoint)`, reading
`adjoint_files_2d[fi].file` and `adjoint_files_2d[fi].variant`, and guarding its
legacy `NFFT(adjoint)` block on `!adjoint_files_2d[fi].variant`.

- [ ] **Step 2: Run the test to verify it exercises the new files**

```bash
clang-format -i tests/nfast.c
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
```

Expected: `exit=0`, with 10 cases per 2D loop instead of 8. `read_nd_case` is wrapped
in `CU_ASSERT_TRUE_FATAL`, so a case the loop never reaches cannot fail quietly.
Confirm the new files really are read by hiding one and watching the suite fail -- move
the data file, never edit the source you just wrote:

```bash
mv tests/data/nfft_2d_10_20_50_t210.txt /tmp/
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"   # expect exit=1
mv /tmp/nfft_2d_10_20_50_t210.txt tests/data/
tests/checkall_ng; echo "exit=$?"                                   # back to exit=0
```

- [ ] **Step 3: Commit**

```bash
git add tests/nfast.c
git commit -m "test(nfast): run the 2D type-II reference cases through the fast path"
```

---

### Task 7: acceptance coverage -- fast path forced across odd, type-II and the reference grid

**Files:**
- Modify: `tests/nplan.c` (`odd_case`, `Y(check_nplan_odd_n)`, `type_ii_1d_once`,
  `Y(check_nplan_type_ii_1d)`, `Y(check_nplan_type_ii_nd)`)
- Modify: `tests/nplan_data.c` (`Y(check_nplan_data)` plus a new `run_native_fast`)

**Interfaces:**
- Consumes: `ndft_ref_trafo(d, N, variant, Ntot, M, x, f_hat, f)`,
  `ndft_ref_adjoint(d, N, variant, Ntot, M, x, f, f_hat)`, `fill_nodes`, `fill_fhat`,
  `rel_max_err` in `tests/nplan.c`; `read_case`, `rel_max_errC`, `native_testcases`,
  `native_testcases_count` in `tests/nplan_data.c`.
- Produces: `static void fast_case(int d, const INT *N, const int *variant, const INT *n, INT M, unsigned seed)`
  and `static void assert_plan_names(Y(plan_ng) *p, const char *needle)` in
  `tests/nplan.c`; `static void run_native_fast(const native_testcase_t *tc)` in
  `tests/nplan_data.c`.

The existing tests all planned under `NFFT_ESTIMATE` while the fast declined, so a
direct solver served them. After Task 5 the fast can win those, silently swapping what
is under test. This task pins both sides: the existing calls keep testing direct by
adding `NFFT_NO_FAST_NATIVE`, and new calls force the fast with `NFFT_NO_DIRECT`.

- [ ] **Step 1: Add the fast-path helpers and calls in `tests/nplan.c`**

Insert `assert_plan_names` and `fast_case` immediately before `Y(check_nplan_odd_n)`:

```c
/* The plan tree must name the given solver, so a test that means to exercise
 * the fast path fails loudly if a direct solver served it instead. */
static void assert_plan_names(Y(plan_ng) * p, const char *needle) {
  FILE *tmp = tmpfile();
  char buf[1024];
  size_t got;
  CU_ASSERT_PTR_NOT_NULL_FATAL(tmp);
  Y(fprint_plan)
  (p, tmp);
  rewind(tmp);
  got = fread(buf, 1, sizeof(buf) - 1, tmp);
  buf[got] = '\0';
  fclose(tmp);
  CU_ASSERT_PTR_NOT_NULL(strstr(buf, needle));
}

/* One geometry through the composed fast solver, forced with NFFT_NO_DIRECT and
 * checked against the exact NDFT reference. The bound is looser than odd_case's
 * because the fast path is an approximation: m = 6 at sigma = 2 with the
 * Kaiser-Bessel window lands near 1e-10, and the assertion leaves headroom for
 * -ffast-math reassociation. */
static void fast_case(int d, const INT *N, const int *variant, const INT *n,
                      INT M, unsigned seed) {
  INT Ntot = 1, j;
  int t;
  R *xin;
  C *f_hat, *f, *ref_f, *ref_fhat;
  Y(plan_ng) * p;
  R tol = (R)1.0e-6;
  R err;
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;
  for (t = 0; t < d; t++)
    Ntot *= N[t];
  xin = (R *)Y(malloc)((size_t)(d * M) * sizeof(R));
  f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  f = (C *)Y(malloc)((size_t)M * sizeof(C));
  ref_f = (C *)Y(malloc)((size_t)M * sizeof(C));
  ref_fhat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  fill_nodes(xin, d, M, seed);
  p = Y(plan_ng_guru)(d, N, variant, n, M, 6, NFFT_WINDOW_KAISER_BESSEL, xin,
                      f_hat, f, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  assert_plan_names(p, "nfft_solver_fast_native");
  Y(precompute)
  (p);
  fill_fhat(f_hat, Ntot, seed + 1000u);
  Y(execute)
  (p);
  ndft_ref_trafo(d, N, variant, Ntot, M, xin, f_hat, ref_f);
  err = rel_max_err(f, ref_f, M);
  CU_ASSERT(err < tol);
  for (j = 0; j < M; j++)
    f[j] = ref_f[j];
  Y(execute_adjoint)
  (p);
  ndft_ref_adjoint(d, N, variant, Ntot, M, xin, ref_f, ref_fhat);
  err = rel_max_err(f_hat, ref_fhat, Ntot);
  CU_ASSERT(err < tol);
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (xin);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(free)
  (ref_f);
  Y(free)
  (ref_fhat);
}
```

Pin `odd_case` to the direct path by changing its one guru call:

```c
  p = Y(plan_ng_guru)(d, N, 0, n, M, 6, NFFT_WINDOW_KAISER_BESSEL, xin, f_hat, f, 0u,
                      NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
```

and update its leading comment to say it forces the direct native.

Append fast-path cases to `Y(check_nplan_odd_n)`, before the
`Y(the_planner_destroy)()` call. Every geometry here satisfies `N_t > 6`,
`n_t > 14` and `n_t > N_t`:

```c
  /* The same shapes through the composed fast solver. */
  {
    INT N[1] = {15}, n[1] = {32};
    fast_case(1, N, 0, n, 500, 41u);
  }
  {
    INT N[1] = {25}, n[1] = {64};
    fast_case(1, N, 0, n, 500, 42u);
  }
  {
    INT N[2] = {15, 10}, n[2] = {32, 24};
    fast_case(2, N, 0, n, 400, 43u);
  }
  {
    INT N[2] = {9, 9}, n[2] = {20, 20};
    fast_case(2, N, 0, n, 400, 44u);
  }
  {
    INT N[3] = {15, 10, 7}, n[3] = {32, 24, 16};
    fast_case(3, N, 0, n, 200, 45u);
  }
  {
    INT N[4] = {7, 7, 7, 7}, n[4] = {16, 16, 16, 16};
    fast_case(4, N, 0, n, 100, 46u); /* rank >= 4 odometer */
  }
  {
    INT N[3] = {15, 16, 9}, n[3] = {32, 32, 20};
    int variant[3] = {NFFT_NDFT_TYPE_I, NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I};
    fast_case(3, N, variant, n, 200, 47u); /* odd and type-II on one plan */
  }
```

Pin `type_ii_1d_once` to whatever the caller steers by leaving its guru call as is, and
change the driver to cover both sides:

```c
void Y(check_nplan_type_ii_1d)(void) {
  type_ii_1d_once(NFFT_NO_FAST_NATIVE | NFFT_NO_NDFT_BLOCKED); /* plain direct */
  type_ii_1d_once(NFFT_NO_FAST_NATIVE | NFFT_NO_NDFT_PLAIN);   /* blocked direct */
  type_ii_1d_once(NFFT_NO_DIRECT);                             /* composed fast */
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}
```

`type_ii_1d_once` compares against the exact reference at tol 1e-10, which the fast
path cannot meet, so widen its bound when the fast serves the plan:

```c
  R tol = (steer & NFFT_NO_DIRECT) ? (R)1.0e-6 : (R)1.0e-10;
```

and leave the `1.0e4 * EPSILON` floor and the `> 1.0e-3` type-I-differs check alone.

Convert `Y(check_nplan_type_ii_nd)` the same way: rename its body to
`static void type_ii_nd_once(unsigned steer)`, thread `steer` into the guru flags
(`NFFT_ESTIMATE | steer`) and the tolerance exactly as above, drop the trailing
`Y(the_planner_destroy)()` from the body, and add the driver:

```c
/* per-axis type-II in the nD native (mixed type-I/type-II axes), through the
 * direct native and through the composed fast. */
void Y(check_nplan_type_ii_nd)(void) {
  type_ii_nd_once(NFFT_NO_FAST_NATIVE);
  type_ii_nd_once(NFFT_NO_DIRECT);
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}
```

`N = {16, 16}`, `n = {32, 32}`, `m = 6` passes the guard, so the `NO_DIRECT` pass
plans.

- [ ] **Step 2: Mirror the reference grid through the fast path in `tests/nplan_data.c`**

Add `run_native_fast` after `run_native`:

```c
/* The same reference cases through the composed fast solver, forced with
 * NFFT_NO_DIRECT. Cases whose geometry fails guards_ok (N <= m, n <= 2m+2, or
 * sigma <= 1 at n = 2N) have no fast plan and are skipped; with m = 6 and
 * n = 2N that leaves every axis with N > 6. */
static void run_native_fast(const native_testcase_t *tc) {
  const int m = 6;
  int d, t, ok = 1;
  INT *N, NN, M, *n;
  R *x;
  C *f_hat, *f;
  Y(plan_ng) * p;
  R tol = (R)1.0e-6;
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;
  if (!read_case(tc->filename, &d, &N, &NN, &M, &x, &f_hat, &f)) {
    CU_FAIL("read_case failed");
    Y(free)
    (N);
    Y(free)
    (x);
    Y(free)
    (f_hat);
    Y(free)
    (f);
    return;
  }
  n = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  for (t = 0; t < d; t++) {
    n[t] = 2 * N[t];
    if (!(N[t] > (INT)m) || !(n[t] > (INT)(2 * m + 2)) || !(n[t] > N[t]))
      ok = 0;
  }
  if (ok) {
    if (tc->kind == 0) /* trafo: input f_hat, expected f */
    {
      C *save_fhat = (C *)Y(malloc)((size_t)NN * sizeof(C));
      C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
      INT j;
      for (j = 0; j < NN; j++)
        save_fhat[j] = f_hat[j];
      p = Y(plan_ng_guru)(d, N, tc->variant, n, M, m, NFFT_WINDOW_KAISER_BESSEL,
                          x, f_hat, got, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL_FATAL(p);
      Y(precompute)
      (p);
      for (j = 0; j < NN; j++)
        f_hat[j] = save_fhat[j];
      Y(execute)
      (p);
      CU_ASSERT(rel_max_errC(got, f, M) < tol);
      Y(plan_ng_destroy)
      (p);
      Y(free)
      (save_fhat);
      Y(free)
      (got);
    } else /* adjoint: input f, expected f_hat */
    {
      C *save_f = (C *)Y(malloc)((size_t)M * sizeof(C));
      C *got = (C *)Y(malloc)((size_t)NN * sizeof(C));
      INT j;
      for (j = 0; j < M; j++)
        save_f[j] = f[j];
      p = Y(plan_ng_guru)(d, N, tc->variant, n, M, m, NFFT_WINDOW_KAISER_BESSEL,
                          x, got, f, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL_FATAL(p);
      Y(precompute)
      (p);
      for (j = 0; j < M; j++)
        f[j] = save_f[j];
      Y(execute_adjoint)
      (p);
      CU_ASSERT(rel_max_errC(got, f_hat, NN) < tol);
      Y(plan_ng_destroy)
      (p);
      Y(free)
      (save_f);
      Y(free)
      (got);
    }
  }
  Y(free)
  (N);
  Y(free)
  (n);
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
}
```

Call it from `Y(check_nplan_data)`:

```c
void Y(check_nplan_data)(void) {
  int i;
  for (i = 0; i < native_testcases_count; i++) {
    const native_testcase_t *tc = &native_testcases[i];
    if (tc->d == 1) {
      run_native(tc, NFFT_NO_NDFT_BLOCKED); /* plain */
      run_native(tc, NFFT_NO_NDFT_PLAIN);   /* blocked */
    } else
      run_native(tc, 0u);
    run_native_fast(tc); /* the composed fast, where the guard admits it */
  }
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}
```

- [ ] **Step 3: Run the suite**

```bash
clang-format -i tests/nplan.c tests/nplan_data.c
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"
grep -B6 "<FAILURE" CUnitAutomated_ng-Results.xml | head -60
```

Expected: `exit=0`.

- [ ] **Step 4: Prove the new coverage is live, not skipped**

`run_native_fast` skips silently by design, so confirm it runs something. Temporarily
make the fast path wrong and check the suite notices:

```bash
sed -i 's#pln->Nneg = pln->N / 2 -#pln->Nneg = pln->N / 2 + 1 -#' kernel/deconv/deconv-1d.c
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"   # expect exit=1
git checkout kernel/deconv/deconv-1d.c
make -j -C tests checkall_ng && tests/checkall_ng; echo "exit=$?"   # back to exit=0
```

- [ ] **Step 5: Commit**

```bash
git add tests/nplan.c tests/nplan_data.c
git commit -m "test(nplan): force both the direct and the fast path across odd and type-II cases"
```

---

### Task 8: full matrix, benchmarks, and close out the spec

**Files:**
- Modify: `docs/plans/nfft-fast-typeii-odd-n.md`
- Modify: `CONTEXT.md`
- Modify: `docs/reviews/feature-planner-review.md` (mark the three DECONV items fixed)

**Interfaces:**
- Consumes: nothing. This task adds no symbol.

- [ ] **Step 1: Run the full test suite in double precision**

```bash
make -j && make check
```

Expected: PASS summary with no failures, for `tests/checkall`, `tests/checkall_ng`
and their `_threads` twins.

- [ ] **Step 2: Run the full suite in float precision**

Configure a separate tree -- precisions must not share a build directory:

```bash
B=/tmp/claude-2000/-workspaces-nfft/5aa7a749-fcb3-4b7c-8ae9-099f2d3f9056/scratchpad/nfft-float
mkdir -p "$B" && cd "$B"
/workspaces/nfft/configure --enable-all --enable-openmp --enable-tests --enable-float
make -j && make check
cd /workspaces/nfft
```

Expected: PASS. If a fast-path tolerance fails only in float, widen that one bound and
say why in the comment -- do not widen the double bounds to match.

- [ ] **Step 3: Build long double, do not run**

```bash
B=/tmp/claude-2000/-workspaces-nfft/5aa7a749-fcb3-4b7c-8ae9-099f2d3f9056/scratchpad/nfft-ld
mkdir -p "$B" && cd "$B"
/workspaces/nfft/configure --enable-all --enable-openmp --enable-tests --enable-long-double
make -j
cd /workspaces/nfft
```

Expected: clean build. The long-double test run is too slow on this host; build only.

- [ ] **Step 4: Confirm the even type-I fast path**

The 2D and 3D DECONV rewrites replaced a hand-unrolled quadrant/octant loop with block
runs. The expectation is flat-to-faster, not a regression:

- Per output element the loads and stores are identical; only the trip count changes,
  and the `phi_hut_inv` tables are small enough to stay in L1 either way.
- Forward DECONV at rank >= 2 opens with `memset(g_hat, 0, ntot)` while the scatter
  writes `Ntot` cells, so at sigma = 2 the zeroing moves 4x more memory in 2D and 8x in
  3D than the multiply loop does. The zeroing dominates the step.
- Grouping the real factors cuts the multiply count from 4 to 3 per element in 2D and
  from 6 to 3 in 3D.
- `pcost` puts DECONV at `2*Ntot` against `5*ntot*log2(ntot)` for the FFT alone, so
  DECONV is a sub-1% slice of the transform at any interesting size.

Measure anyway:

```bash
cmake -S . -B build-cmake -DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON \
      -DCMAKE_C_FLAGS="-O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math"
cmake --build build-cmake -j
CODSPEED_PROFILE_FOLDER=/tmp/wt-after build-cmake/benchmarks/bench_nfft_direct
```

Compare against the same run on `HEAD~7` (the commit before Task 1). Report the delta.
Do not add a second code path for the equal-halves case: keeping the unroll alongside
the general form needs per-axis remainder handling (`Nneg` and `Npos` differ by 1 for
odd N and by 2 for type-II), which in 2D means a main block, two edge strips and a
corner, per direction. If the benchmark does show a real regression, bring the number
back and decide then -- do not pre-emptively branch.

- [ ] **Step 5: Update the docs**

In `docs/plans/nfft-fast-typeii-odd-n.md`, mark the three defects fixed and the work
items landed, naming the commits.

In `CONTEXT.md`, update whatever states the fast NFFT's supported geometry: it now
serves odd `N` and per-axis type-II, and requires `sigma > 1`.

In `docs/reviews/feature-planner-review.md` section 6, mark the odd-`N` zero-pad, the
ignored `variant` in rank >= 2, and the missing oversampling guard as fixed.

Record the wisdom consequence in all three docs. The wisdom key already hashes
`variant[]` and carries `N` parity through the tensor dims, and the configuration
signature covers only `sizeof(R)` plus the solver roster -- which this work does not
change. So no key changes and no wisdom file is invalidated. What does change is that
entries blessed before this work, for odd or type-II problems, name a direct solver
that the planner would no longer pick. They stay correct, just slower. Users who want
the fast path for a shape they already measured must call `nfft_forget_wisdom()` or
delete their wisdom file. State that plainly; do not bump `PACKAGE_VERSION` for it.

- [ ] **Step 6: Deslop and commit**

```bash
# Run the deslop skill over every C file touched by Tasks 1-7:
#   kernel/deconv/deconv-{1d,2d,3d,nd}.c kernel/nfft/nfft-nd.c
#   tests/nfast.c tests/nplan.c tests/nplan_data.c
git add -A
git commit -m "docs: record fast-NFFT type-II and odd-N support"
```

---

## Verification checklist

Before calling this done, confirm each by running the command, not by reading the code:

- [ ] `make check` passes in double (`tests/checkall`, `tests/checkall_ng`, both
      `_threads` twins).
- [ ] `make check` passes in a separate float tree.
- [ ] The long-double tree builds.
- [ ] `tests/checkall_ng` includes the four new DECONV suites: `deconv_1d_general`,
      `deconv_2d_general`, `deconv_3d_general`, `deconv_nd_general`.
- [ ] `grep -n "N % 2 == 0" kernel/deconv/` returns nothing.
- [ ] `grep -rn "nfast\.c" kernel/nfft/*.c kernel/deconv/*.c tests/*.c` returns nothing.
      `tests/nfast.c` is a real file, so grep the sources only -- every remaining hit is
      a comment naming the deleted kernel file `kernel/nfft/nfast.c`.
- [ ] The even type-I 2D/3D benchmark is flat against the pre-Task-1 commit.
