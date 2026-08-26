# Lifecycle & contracts

The rules a caller (or a test, or a new example) must obey. Ground truth:
`kernel/nfft/plan_ng.c`, `include/iplanner.h`, `include/nfft3.h`.

## The three awake states

`include/iplanner.h`:

```c
enum {
  PLNR_SLEEPY     = 0,  /* not ready */
  PLNR_AWAKE_ZERO = 1,  /* may execute, results may be INCORRECT (planning state) */
  PLNR_AWAKE      = 2   /* may execute, results correct */
};
```

`Y(plan_awake)(plan, wakefulness)` moves a plan between states and is idempotent
per state; on a transition it calls the plan's `awake()` hook (may be `NULL` for
stateless plans such as the direct NDFT solvers). A fast-native plan's `awake`
hook builds the ψ (window) tables via its DECONV/CONV children, once per awake
period.

- `AWAKE_ZERO` is the **planning-time measurement state**: it builds a
  cost-representative ψ so a candidate can be timed, but its numerical output is
  meaningless. It is **never caller-visible** — only the internal race uses it.
- `AWAKE` is the real, correct state that `precompute()` moves the winner into.
- `plan_destroy` first awakes a plan back to `SLEEPY` (releasing node tables),
  then runs its `destroy` hook, then frees.

## The lifecycle, step by step

```
nfft_plan_ng_guru(...)     construct the bundle; if MEASURE, race candidates
                           SLEEPY -> AWAKE_ZERO (build ψ) -> measure -> SLEEPY.
                           The winner is returned SLEEPY.
        |
nfft_precompute(p)         winner SLEEPY -> AWAKE. Rebuilds ψ once (the accepted
                           one-extra build). MANDATORY. Idempotent.
        |
nfft_execute(p)            asserts winner == PLNR_AWAKE, then apply().
nfft_execute_adjoint(p)    asserts winner == PLNR_AWAKE, then apply_adjoint().
        |
nfft_plan_ng_destroy(p)    SLEEPY-first teardown of the plan(s) + problem(s).
```

### precompute() is mandatory — for every plan

`Y(execute)`, `Y(execute_adjoint)`, `Y(execute_on)`, `Y(execute_adjoint_on)` all
assert `p->dir[FWD]->awake_state == PLNR_AWAKE` (`plan_ng.c`). This is a
**uniform rule** — direct NDFT plans need no ψ, but they still must be walked
`SLEEPY -> AWAKE` by `precompute()` so the state assertion holds.

- Under `--enable-debug`: skipping `precompute()` fires the assert (abort).
- In a **release build `A(...)` is a no-op**: skipping `precompute()` executes a
  `SLEEPY` plan — undefined behavior / crash with no diagnostic.

> The single most common "planner bug" is a test or example that calls
> `execute()` without `precompute()`. If you see a crash in `apply`, check this
> first.

`precompute()` is deferrable and idempotent: build the plan once, fill nodes,
`precompute()`, then `execute()` as many times as you like.

## Array ownership

Passed to the guru: `R *x` (length `d*M`), `C *f_hat` (length `∏N_t`),
`C *f` (length `M`). **All three are required in every mode** (no deferred-fill,
no setter API — FFTW's guru contract).

| Array | Ownership | Notes |
|-------|-----------|-------|
| `x` (nodes) | **Copied** into a plan-owned buffer at construction. | You keep and free your own `x`; the plan never writes or frees it. Under unit-axis elision the copy also *gathers* live-axis columns (see below). The debug md5 guard protects this internal copy. |
| `f_hat` | **Aliased** (borrowed), never freed by the plan. | Forward `apply` reads it; adjoint writes it. |
| `f` | **Aliased** (borrowed), never freed by the plan. | Forward `apply` writes it; adjoint reads it. |

Consequences:

- **Fill `f_hat`/`f` *after* the guru returns**, because a measured race may
  clobber them (the race is destructive on values — it times the real
  node-driven access pattern, values are irrelevant). If you must preserve
  existing contents during planning, plan on scratch and use the `_on` variants.
- There are **no** `plan_ng_x` / `plan_ng_f_hat` / `plan_ng_f` accessors. They
  were removed. The arrays are yours; you already hold the pointers.
- Because `x` is copied, unit axes (`N_t == 1`) are elided at construction and
  the `x` columns are gathered to live axes in the same pass — so a top-level
  `{64,64,1}` problem behaves and hashes identically to `{64,64}`. See
  [solvers-problems-windows.md](solvers-problems-windows.md#unit-axis-elision).

## `execute` vs `execute_on` (new-array execute)

`execute` / `execute_adjoint` run on the plan's bound `f_hat`/`f`.

`execute_on(p, f_hat, f)` / `execute_adjoint_on(p, f_hat, f)` are FFTW
`execute_dft`-style **new-array** variants: they save the problem's `f_hat`/`f`
pointers, swap in your arrays, `apply`, then **restore unconditionally**. `x`
and the ψ tables are untouched (they depend on `x` alone).

- The buffer order is always `(f_hat, f)` for both forward and adjoint (the
  direction is chosen by which function you call, not by argument order).
- No alignment constraint.
- **Not reentrant on the same plan**: two threads calling `_on` on one plan race
  on the swapped pointers — same caveat as `fftw_execute_dft`. Use one plan per
  thread, or serialize.

## Directionality: the forward-only race

A bundle owns a two-slot substrate (`prob[FWD]`/`dir[FWD]` with `sign=+1`,
`prob[ADJ]`/`dir[ADJ]` with `sign=-1`), but **only `[FWD]` is built and raced**.
The single winning forward plan serves the adjoint too, via its
`apply_adjoint()`. The adjoint is **not**
independently optimized or blessed.

- `prob[ADJ]`/`dir[ADJ]` stay `NULL` (dormant substrate for a possible future
  two-internal-plan option). `problem_nfft.sign` is retained as a
  direction-dependent-wisdom hook but is effectively always `+1`.
- `nfft_fprint_plan` therefore prints `(nfft-plan-ng (fwd <plan>) (adj (null)))`
  — the `(null)` adjoint is expected, not a bug.
- Flags `NFFT_FORWARD_ONLY`/`NFFT_ADJOINT_ONLY` do **not** exist (they were
  dropped when the race became forward-only).

## Input validation (release-safe)

`nfft_plan_ng_guru` validates explicitly and returns `NULL` (FFTW's contract)
rather than relying on debug-only `A(...)`:

```c
if (d <= 0 || N == 0 || n == 0 || x == 0 || f_hat == 0 || f == 0) return 0;
for (t = 0; t < d; t++) if (N[t] <= 0 || n[t] <= 0) return 0;  /* N[t]==1 is valid */
```

Checked before any allocation, so a `NULL` return leaks nothing. **Everything
else** — execute-before-precompute, double-destroy, the x-restore guard — is
`A(...)`-gated and therefore **only enforced under `--enable-debug`**. Treat
those as caller discipline in release.

## Thread-safety

- **Planning is not thread-safe** (global `the_planner()`, mirroring FFTW). Do
  not construct/destroy plans or import/export wisdom concurrently.
- **Executing** a plan is fine; **`_on` on the same plan is not reentrant**
  (pointer swap). Thread count is process-global (`Y(get_num_threads)()`, read
  once per planning entry) and folded into the wisdom key; there is no per-plan
  thread knob.

## Destroy

`nfft_plan_ng_destroy(p)` is **NULL-safe**. It destroys each present direction
plan (SLEEPY-first, releasing ψ), then the problem(s) (they outlive plan
selection because their tensors back the plans' `pcost`), then frees the bundle.
