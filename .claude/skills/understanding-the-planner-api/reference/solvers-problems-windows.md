# Solvers, problems, windows — and how to extend them

Ground truth: `include/iplanner.h` (problem structs + register prototypes),
`kernel/nfft/{conf,nfast,ndft,ndft-nd,nconst}.c`, `kernel/deconv/*`,
`kernel/conv/*`, `kernel/util/window.c`.

## The three problem kinds

The fast NFFT is *decomposed* into node-independent and node-dependent stages,
each its own problem kind so the planner can select a solver per stage.

```c
NFFT_PROBLEM_NFFT    /* the full NDFT problem (top level) */
NFFT_PROBLEM_DECONV  /* forward step A: deconvolve + zero-pad f_hat -> oversampled grid g */
NFFT_PROBLEM_CONV    /* forward step C: convolve grid g -> node samples f  (node-dependent) */
```

`problem_nfft` (the top-level problem, `iplanner.h`) carries: `tensor *sz`
(frequency geometry, oriented in dataflow direction — forward `N_t -> n_t`,
adjoint is its `tensor_adjoint`), `tensor *vecsz` (batch loops), `INT M`,
`int m`, `int window`, `int sign`, `unsigned fftw_flags`, `R *x`, `C *f_hat`,
`C *f`, `int *variant` (per-axis NDFT type), `int x_owned` (1 = private copy,
0 = borrowed slice). Constructor:

```c
problem *Y(mkproblem_nfft)(int d, const INT *N, const int *variant,
                           const INT *n, INT M, int m, int window, int sign,
                           unsigned fftw_flags, R *x, int copy_x, C *f_hat, C *f);
```

`copy_x = 1` (top-level) copies `x` into a private buffer *and* gathers live-axis
columns during unit-axis elision; `copy_x = 0` (a child) borrows a slice.

`problem_deconv` and `problem_conv` are the children (borrowed `f_hat`/`g`/`f`,
own `sz`/`vecsz`/`N`). See `iplanner.h` for their exact fields and constructors
(`Y(mkproblem_deconv)`, `Y(mkproblem_conv)`).

## The NFFT solver roster (current)

Registered in `kernel/nfft/conf.c` (`the_roster`), **registration order shown**
(iteration is the reverse — see
[planner-internals.md](planner-internals.md#registration-order-determines-tie-breaks)):

| # | Solver (registrar) | Kind / rank | Gated by | Notes |
|---|--------------------|-------------|----------|-------|
| 1 | `nfft_solver_fast_native` | any rank | `NFFT_NO_FAST_NATIVE` | The composed fast NFFT. Registered *first* → iterated *last* → loses exact ties. |
| 2 | `nfft_solver_ndft_1d` | rnk 1 | `NFFT_NO_NDFT_PLAIN` (+ `NFFT_NO_DIRECT`) | Plain per-term direct NDFT. |
| 3 | `nfft_solver_ndft_1d_blocked` | rnk 1 | `NFFT_NO_NDFT_BLOCKED` (+ `NFFT_NO_DIRECT`) | Blocked-recurrence 1D NDFT (fewer transcendentals). |
| 4 | `nfft_solver_ndft_nd` | rnk ≥ 2 | `NFFT_NO_NDFT_PLAIN` (+ `NFFT_NO_DIRECT`) | One generic direct NDFT (no 2D/3D specialization — direct has none). `NFFT_NO_NDFT_PLAIN` is "plain per-term, **any rank**", so it forbids both `ndft_1d` and `ndft_nd`. |
| 5 | `nfft_solver_const_0d` | rnk 0 | ungated | Rank-0 base case: forward = broadcast `f_hat[0]`, adjoint = reduce `Σ f_j`. The sole rank-0 solver, so it never ties. |

All direct/NDFT solvers are **serial** (OpenMP was stripped from the native
direct path). There are no `nfft_solver_direct` or `fast_1d/2d/3d/nd` solvers —
those were the removed wrapper roster.

## The native fast NFFT (`kernel/nfft/nfft-nd.c`)

`nfft_solver_fast_native` is **coreless** and **rank-general**. Its `mkplan`
recurses into `Y(planner_mkplan)` to build children:

```
forward  apply:  DECONV(f_hat -> g1) -> FFTW forward (g1 -> g2) -> CONV(g2 -> f)
adjoint  apply:  CONV^H(f -> g2)    -> FFTW backward (g2 -> g1) -> DECONV^H(g1 -> f_hat)
```

- **Step B is a direct FFTW plan**, not a problem kind (there is no
  `NFFT_PROBLEM_FFT`). Adjoint uses an **unnormalized `FFTW_BACKWARD`**.
- DECONV is a node-independent diagonal (divide by the window's Fourier
  coefficients); its adjoint multiplies by the **same** `1/phi_hut`
  (self-adjoint real diagonal).
- CONV is the node-dependent window convolution; it owns the ψ table, built in
  `awake()` from the borrowed `x`.
- The internal FFTW child plan defaults to `FFTW_ESTIMATE`.
- The DECONV/CONV leaves are **rank-disjoint** serial transplants: 2D
  (`rnk==2`), 3D (`rnk==3`), generic (`rnk≥4`) — rank alone selects the child.

**Applicability guard** (`guards_ok` in `kernel/nfft/nfft-nd.c`): per axis,
requires `N_t > m`, `n_t > 2m + 2`, and oversampling `sigma = n_t/N_t > 1`
strictly. Odd `N` and per-axis type-I/type-II (mixed across axes) are both
served — the DECONV children carry the per-axis slot split for this. A
geometry that fails the guard for a size reason (returns `NULL`) falls back to
the direct NDFT. **`sigma <= 1` normally never reaches the guard through the
public API**: `plan_ng_guru` rejects it up front and returns `NULL` for the
whole plan, rather than quietly serving a direct transform the caller did not
ask for and cannot see. That rejection is lifted by `NFFT_NO_FAST_NATIVE`,
which takes the fast path out of the running so nothing is lost
unintentionally. The guard keeps its own `sigma > 1` check for callers who
build a problem through `mkproblem_nfft` directly. It is currently serial +
all four real windows (below).

Every rank >= 1 solver additionally declines any axis with `N_t == 1`
(`Y(problem_nfft_has_unit_axis)`). Unit axes are elided at construction, so one
reaching a solver means compression was bypassed. The rank-0 base case is what
serves a problem whose axes were all elided away.

## NDFT variants (type-I / type-II)

Per-axis, via the guru's `variant` array (`NULL` = all type-I):

- `NFFT_NDFT_TYPE_I` (`0`): for even `N`, `k = -N/2 .. N/2-1`.
- `NFFT_NDFT_TYPE_II` (`1`): for even `N`, `k = -N/2+1 .. N/2`.
- Odd `N` has only one type (defined as type-I), range `-(N-1)/2 .. (N-1)/2`;
  the constructor normalizes odd axes to type-I.

Type-II and odd `N` are supported by both the **direct NDFT** solvers and the
**fast** path, subject to the applicability guard above (the size guards; on
the guru path `sigma > 1` has already been enforced on every non-unit axis
unless `NFFT_NO_FAST_NATIVE` lifted it).

## Unit-axis elision

Unit axes (`N_t == 1`) are elided **at construction**, inside `mkproblem_nfft`'s
top-level (`copy_x = 1`) path — FFTW parity. The `x` columns for dropped axes are
gathered to live axes in the same copy pass (`xc[j*dl + t] = x[j*d + live[t]]`).
It is **drop-only, no reordering** (positional `x` pairing forbids a canonical
stride sort). `f_hat`/`f` need no repack (unit axis = stride 1). All axes unit →
**rank 0**, served by `nfft_solver_const_0d`. Consequence:
`{64,64,1}` hashes and behaves identically to `{64,64}`. (Borrowed child
problems keep full rank.)

Because elision is unconditional on the guru path, **no rank >= 1 NFFT solver
accepts a surviving unit axis** — each declines via
`Y(problem_nfft_has_unit_axis)`. A unit axis reaching a solver means
compression was bypassed, so it is surfaced rather than served. Only the
borrowed (`copy_x = 0`) path can build such a problem, and only rank 0 (every
axis elided) is servable, by `nfft_solver_const_0d`.

## The window vtable (`kernel/util/window.c`)

Windows are dispatched at **runtime** by an `NFFT_WINDOW_*` ordinal (the window
is a per-plan guru parameter, hashed in the wisdom key). The vtable:

```c
R Y(window_phi_hut)(int window, INT n, INT N, int m, INT k);  /* Fourier coeff */
R Y(window_phi)(int window, INT n, INT N, int m, R x);        /* window value  */
/* range-apply helpers used by the DECONV/CONV children: */
void Y(window_phi_hut_apply)(int window, INT n, INT N, int m, INT k0, R *out, INT count);
void Y(window_phi_precompute)(int window, INT n, INT N, int m,
                              const R *x, INT x_stride, INT num_nodes, R *out, INT out_stride);
```

Implemented in the fast path: **Kaiser-Bessel** (default), **Gaussian**,
**B-spline**, **sinc-power**. `NFFT_WINDOW_DIRAC_DELTA` is **declined by every
fast solver** (the `default:` returns 0.0) — reachable only via direct NDFT.

**KB self-normalization (single-precision fix).** The Kaiser-Bessel `phi_hut`/ψ
**self-scale** by `exp(-log I0(m·b))` (peak → 1), computed in `window.c` with the
`kernel/util/bessel_i0.c` helpers. Without this, the deconvolved grid
`g1 = f_hat / phi_hut^d` underflows to float32 denormals at `d = 3` (KB-specific;
the composed `CONV·FFT·DECONV` cancels the scaling exactly). It is exact in
double and recovers float accuracy. **Native path only** — the legacy `nfft.c`
window is untouched (legacy float32 3D stays broken by design). The per-axis
`log`-peak is hoisted once per axis by the range-apply helpers, so precompute
cost is at or below baseline.

## Adding a solver or a new kind

1. **Write the `solver_adt`** with `problem_kind` and a `mkplan` that returns a
   plan with `pcost` set or `NULL` (cheap; no node-dependent state — build ψ in
   `awake`). Provide a `Y(..._register)(planner *)` that
   `REGISTER_SOLVER`s a `Y(solver_create)`'d instance.
2. **Add it to the kind's roster** (`the_roster` in `conf.c`) via `SOLVTAB(...)`.
   Mind **registration order**: register a candidate *earlier* to make it lose
   exact-`pcost` ties (iteration is reverse order).
3. **Register cross-kind children EAGERLY.** If your solver's `mkplan` recurses
   into `Y(planner_mkplan)` for another kind (as fast-native does for
   DECONV/CONV), those child solvers **must already be registered** before the
   recursion. Do it in `Y(nfft_ensure_registered)` (see `conf.c`), *before* the
   parent roster — never lazily from inside a `mkplan`. Lazy registration
   reallocs `slvdescs` while the outer `FORALL_SOLVERS_OF_KIND` search still
   holds a raw pointer into it: a use-after-free that reliably crashes on the
   **second** `plan_ng_guru` call in a process, not the first.
   "Reentrancy-safe ≠ registration-safe."
4. **A brand-new problem kind** needs an entry in the `NFFT_PROBLEM_*` enum
   (`iplanner.h`) — note `kind_head[NFFT_PROBLEM_LAST]` hardcodes the kind
   universe, so the enum is the single point of truth. The NFFT `problem_*`
   structs and window decls currently live in the *generic* `iplanner.h`; a
   first port of another transform (NFSFT/NSFFT/FPT) should split those out.
5. **Roster changes change the configuration signature**, so old wisdom becomes
   a safe cache miss (see [wisdom.md](wisdom.md#the-configuration-signature-the-roster-fingerprint)).

## Provenance

`kernel/planner/` (and the per-kind planner glue) follows the architecture of
FFTW3's planner; the implementation was written for this repository. Both
projects are GPL-2.0-or-later. `cycle.h` is FFTW's MIT-licensed file, kept
under FFTW's copyright header for the timing path. Do not paste FFTW code or
comments into planner sources; see `CONTEXT.md` ("provenance").
