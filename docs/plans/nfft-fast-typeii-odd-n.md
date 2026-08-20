# Plan: type-II and odd `N` in the planner-native fast NFFT

Status: landed. Target branch `feature/planner`.

## 1. Goal

`nfft_solver_fast_native` (`kernel/nfft/nfft-nd.c`) must accept every geometry
the native direct NDFT solvers already accept in the new planner API:

| Case | direct NDFT | fast native (today) | fast native (target) |
|---|---|---|---|
| even `N`, type-I | yes | yes | yes |
| even `N`, type-II (per axis, mixed) | yes | **declines** | yes |
| odd `N` (type-I only by definition) | yes | **declines** | yes |
| unit axis `N=1` | yes (elided, rank-0) | elided upstream | unchanged |
| rank 0 | `nfft_solver_const_0d` | n/a | unchanged |

Out of scope: `NFFT_WINDOW_DIRAC_DELTA` (no fast path by design), NFCT/NFST
(not in the new API), OpenMP.

## 2. Where the variant and the parity actually enter

Frequency convention, fixed by the direct solvers
(`kernel/nfft/ndft-1d.c`, `ndft-nd.c`) and by `mkproblem_nfft`:

```
slot ks = 0 .. N-1      (row-major index into f_hat, per axis)
k(ks)   = ks - N/2 + s  (C integer division; s = 1 for even-N type-II, else 0)
```

* even, type-I:  `k = -N/2 .. N/2-1`
* even, type-II: `k = -N/2+1 .. N/2`
* odd:           `k = -(N-1)/2 .. (N-1)/2`, `s` forced to 0 by the constructors

The composed algorithm is
`DECONV (f_hat -> g1) -> FFTW -> CONV (g2 -> f)`.

Write `g1[p] = f_hat[ks] / phi_hut(k(ks))` at `p = k(ks) mod n`. The FFT then
gives `G_l = sum_k (f_hat_k / phi_hut(k)) exp(-2 pi i k l / n)` for **any**
index set, because `exp(-2 pi i (k mod n) l / n) = exp(-2 pi i k l / n)`.
CONV then evaluates the window sum over the grid.

**Consequence: the variant and the parity are confined to DECONV.**

* CONV (`kernel/conv/*`) reads `N` only to derive the window shape
  (`sigma = n/N` inside `window_phi`), never the index set. No change.
* `kernel/util/window.c` is parity-blind and variant-blind (`phi_hut` is even
  in `k`, so type-II's `k=+N/2` is exactly as well conditioned as type-I's
  `k=-N/2`). No change.
* The FFTW child plan depends on `n` only. No change.

Derived quantities used throughout this plan, per axis:

```
Nneg = N/2 - s     /* count of negative frequencies; slots ks = 0 .. Nneg-1  */
Npos = N - Nneg    /* count of k >= 0;               slots ks = Nneg .. N-1  */
```

Grid position of slot `ks`: `ks < Nneg ? n - Nneg + ks : ks - Nneg`.
For even type-I this is `Nneg = Npos = N/2`, so every existing code path is a
special case and stays bit-for-bit unchanged.

## 3. Defects found while scoping (must be fixed by this work)

1. **Fixed (`0f1bcb00`).** `deconv-1d.c` zero-pad was wrong for odd `N`. It
   cleared `[N/2, N/2 + (n-N))`, but the last non-negative frequency lands at
   `Npos-1`, so cell `n - N/2 - 1` was never written and kept stale data.
   Latent only because `guards_ok` declined odd `N`; the 1D DECONV solver
   itself had no `A(N % 2 == 0)`, so it looked supported.
   Fix: `memset(g + Npos, 0, (n - N) * sizeof(C))`.
2. **Fixed (`eb9d0851`, `d3c679cd`, `024f4fa3`).** `deconv-2d.c` /
   `deconv-3d.c` / `deconv-nd.c` ignored `variant` entirely. A rank >= 2
   type-II DECONV problem hashed to its own wisdom key and then got a type-I
   plan: silent wrong result for a direct `mkproblem_deconv` caller.
3. **Fixed (`706c1ad7`).** No oversampling guard anywhere. `guards_ok` checked
   `N > m` and `n > 2m+2` only. With `n < N`, `(size_t)(n - N)` in the 1D
   `memset` wrapped, and the deconv scatter aliased grid cells. **Decision:
   the fast NFFT requires `sigma = n/N > 1`, i.e. `n > N` strictly, on every
   axis.** This also removes the `sigma == 1` band-edge case where, for even
   `N`, sinc-power `phi_hut` is exactly 0 at `|k| = N/2` and `1/phi_hut` is
   `inf`. The DECONV solvers themselves only need `n >= N` and decline below
   that.

These overlapped with items already recorded in
`docs/reviews/feature-planner-review.md` (section 6); that section is now
marked fixed alongside these.

## 4. Work items

### W1 — `kernel/deconv/deconv-1d.c`

* Already applies the type-II shift in `awake`, `apply`, `apply_adjoint`.
  Keep.
* Replace the zero-pad start `N/2 + type_ii` with `Npos = N - N/2 + type_ii`.
  Identical for even `N` (both variants), correct for odd `N`.
* Hoist `N`, `n`, `type_ii`, `Nneg`, `Npos` into `deconv_plan` at `mkplan`
  instead of re-deriving them from the problem on every `apply`.
* `mkplan` returns 0 when `n < N` (a `return 0` decline, not an `A(...)`,
  so release builds are covered too).

### W2 — `kernel/deconv/deconv-2d.c`, `deconv-3d.c`

The current bodies hard-code a 2x2 (2x2x2) quadrant split of equal halves
`N/2`. Generalize to two runs per axis of lengths `Nneg` and `Npos`:

```c
/* per axis t: run 0 = negative k, run 1 = non-negative k */
len[t][0] = Nneg[t];            len[t][1] = Npos[t];
gof[t][0] = n[t] - Nneg[t];     gof[t][1] = 0;        /* grid offset  */
sof[t][0] = 0;                  sof[t][1] = Nneg[t];  /* f_hat / phi_hut offset */
```

and loop `for (b0 = 0; b0 < 2; b0++) for (b1 = 0; b1 < 2; b1++) { ... }`
around a contiguous inner loop of length `len[d-1][b_{d-1}]`. A zero-length
run (unit axis, `Nneg == 0`) drops out on its own.

Also:
* `awake` must start the `phi_hut` table at `k0 = -N/2 + s` per axis, not
  `-N/2` (the table stays indexed by `ks`).
* Delete the `A(N % 2 == 0)` asserts and the "declined upstream" comments.
* Fold the forward and adjoint mirror bodies into one
  `run(pln, pd, forward)` helper, as `deconv-nd.c` already does. The two
  bodies differ only in the innermost assignment, and doubling the number of
  hand-written index expressions doubles the odd/type-II risk surface.

### W3 — `kernel/deconv/deconv-nd.c`

The odometer generalizes by substituting three constants; the carry logic is
untouched:

| current | replacement |
|---|---|
| `ks[t] = N[t] / 2` (init and carry reset) | `ks[t] = Nneg[t]` |
| `if (kp[t] == N[t] / 2)` (switch to negative block) | `if (kp[t] == Npos[t])` |
| `k[t] = n[t] - N[t] / 2` | `k[t] = n[t] - Nneg[t]` |

* Store `Nneg[]`, `Npos[]` in the plan at `mkplan`.
* `awake` starts each axis table at `-N[t]/2 + s[t]`.
* Drop `A(N[t] % 2 == 0)`; add the `n[t] >= N[t]` decline.

### W4 — `kernel/nfft/nfft-nd.c` (`guards_ok`)

New guard set, per axis:

```c
if (!(Nt > (INT)m))         return 0;
if (!(nt > (INT)(2*m + 2))) return 0;
if (!(nt > Nt))             return 0;   /* new: sigma = n/N > 1 */
```

Removed: the odd-`N` decline, the type-II decline. The explicit unit-axis
decline becomes redundant (`Nt == 1` fails `Nt > m` for every `m >= 1`, and
`nt == 1` fails `nt > 2m+2`); delete it and say so in the comment.

`mkproblem_deconv` is already called with `pn->variant`, in the same
compressed axis order as `pln->Nc` / `pln->nc`. No wiring change.

### W5 — cost model

`pcost` in `nfft-nd.c` is parity- and variant-independent
(`Ntot`, `ntot`, `M`, `m`, `d`). No change. Note that this changes *which*
solver wins for odd / type-II geometries that previously had only the direct
NDFT as a candidate — that is the point of the work, and the estimate
comparison is the same one already used for even type-I.

## 5. Test plan

Existing tests that assert the decline must be **inverted, not deleted**:

| File | Test | Change |
|---|---|---|
| `tests/nplan.c` | `check_nplan_guard_declines` | The odd-`N` and type-II cases must now select `nfft_solver_fast_native` at large `M`. Keep a genuine decline case (`n < N`, or `N <= m`) so the guard is still exercised. The `NO_DIRECT + odd -> guru returns NULL` case must become `NO_DIRECT + odd -> a valid fast plan`; move the NULL assertion to a geometry that fails the new guard. |
| `tests/nfast.c` | `files_2d` comment at ~line 1126 | Drop the "`_t210` variants excluded" exclusion; add `nfft_2d_10_20_50_t210.txt`, `nfft_2d_20_10_50_t210.txt` and the two adjoint files to the fast-forced 2D roster. |
| `tests/nplan.c` | `check_nplan_odd_n`, `check_nplan_type_ii_1d`, `check_nplan_type_ii_nd` | Currently run under `NFFT_ESTIMATE` and are served by the direct solvers. Parameterize each on a steering flag and run twice: `NFFT_NO_FAST_NATIVE` (direct, present behaviour) and `NFFT_NO_DIRECT` (fast). Fast needs a looser tolerance than `1e-10` — use the window bound the other `nfast` accuracy tests use, not the direct-NDFT tolerance. |

New tests:

1. **`nplan_data` fast mirror.** `tests/nplan_data.c:check_nplan_data` drives
   every generated reference case with `NFFT_NO_FAST_NATIVE`. Add a second
   pass with `NFFT_NO_DIRECT`, skipping cases that fail `guards_ok`
   (`N <= m` or `n <= 2m+2`); the grid in `tests/refgen/grids.py` already
   contains odd-`N` and `_t210` entries, so no data regeneration is needed.
   Log the skipped cases so the coverage gap is visible.
2. **DECONV unit tests, rank >= 2.** `tests/nfast.c:check_nfast_deconv_solver`
   only covers 1D type-II. Add a spike test per rank (2, 3, 4): a single
   `f_hat` slot set to 1, assert it lands at the expected wrapped grid cell
   and is divided by `phi_hut` of the expected frequency. Cover
   (even type-I, even type-II, odd) x (rank 2, 3, 4) and mixed per-axis
   variants. This is the test that would have caught defect 2 above.
3. **Odd + type-II zero-pad completeness.** Pre-fill `g` with a sentinel,
   run forward DECONV, assert every untouched cell is exactly zero. This is
   the direct regression test for defect 1.
4. **Adjoint consistency.** For each new geometry, check
   `<A f_hat, f> == <f_hat, A^H f>` to the window tolerance, so a forward-only
   index fix cannot pass.

Precision matrix: run `make check` for double and float. Skip long-double
runs (build only) — they are too slow on this host.

**Result:** double `make check` passes clean. Float `make check` initially
failed three of the four new DECONV general spike tests
(`deconv_2d_general`, `deconv_3d_general`, `deconv_nd_general`) because their
tolerance was a fixed literal (`K(1e-9) * sc`) rather than scaled to the
working precision; observed float relative error peaked around 1.1e-7 (about
1 `eps`), comfortably past that literal. Fixed by replacing the literal with
a helper, `deconv_general_tol()`, returning `64 * eps` (`tests/nfast.c`);
float now passes clean with the same margin reasoning the file's
`err_bound` helper already uses elsewhere (a multiple-of-`eps` floor for
round-off accumulation). The double bounds were not touched. Long-double
built clean; the test run was skipped as planned.

## 6. Wisdom and cache

* No new problem field, no roster change, so the configuration signature
  (`kernel/planner/planner.c:config_signature`, over `sizeof(R)` plus every
  registered solver's `reg_id` and registrar name) is unchanged, and existing
  wisdom files stay importable.
* `variant[]` is already hashed by both `problem_nfft::hash` and
  `problem_deconv::dv_hash`; parity is carried by the tensor dims. Keys are
  correct as-is.
* **Stale-but-safe wisdom:** an entry blessed before this work for an
  odd-`N` or type-II key names a direct NDFT solver and keeps being honoured.
  The result stays correct, only slower. No key changed and no wisdom file is
  invalidated. A user who wants the fast path for a shape they already
  measured before this work must call `nfft_forget_wisdom()` or delete their
  wisdom file to force a fresh search. `PACKAGE_VERSION` is not bumped for
  this change.

## 7. Risks

| Risk | Mitigation |
|---|---|
| Rewriting the 2D/3D quadrant unroll regresses throughput on the even type-I fast path (the common case). | **Resolved:** CodSpeed walltime run (`bench_nfft_direct`, forward and adjoint, 2D and 3D) comparing this branch against `112275bc` (the commit before this work) showed every case within about -3% to +0.7%, i.e. flat to slightly faster, no regression. The unroll was not kept. |
| `sigma = 1` (`n == N`) puts `\|k\| = N/2` at the window band edge: sinc-power `phi_hut` is exactly 0 there, so `1/phi_hut` is `inf`. | **Resolved, then strengthened.** `guards_ok` requires `n > N` strictly, so the fast NFFT never reaches `sigma = 1`. Falling back to a direct NDFT was the first answer and proved wrong: the guru does not reveal which solver won, so a caller passing `n == N` silently lost the fast path and paid `O(N*M)`. `plan_ng_guru` now **rejects** `sigma <= 1` and returns `NULL` on every axis with `N > 1`, whenever the fast solver is in play; `NFFT_NO_FAST_NATIVE` lifts the requirement, since nothing is then lost unintentionally. Unit axes are elided and exempt. The guard stays as the check for direct `mkproblem_nfft` callers, and the DECONV solvers keep the weaker `n >= N` decline, since a DECONV problem alone is well defined at `sigma = 1` for every window except sinc-power. |
| Forward/adjoint mirror bodies drift while being generalized. | W2 folds them into one `run(..., forward)` helper first, then generalizes once. |
| Rank >= 2 DECONV has no per-rank index test today, so an index slip lands silently. | Test item 2 lands before W2/W3. |

## 8. Commit order (landed)

1. `0f1bcb00` — W1 (`deconv-1d`): odd zero-pad fix, `n >= N` decline, hoisted
   geometry.
2. `eb9d0851` — W2 (`deconv-2d`): fold mirrors, generalize to the per-axis
   slot split.
3. `d3c679cd` — W2 (`deconv-3d`): same.
4. `024f4fa3` — W3 (`deconv-nd`): odometer generalization.
5. `706c1ad7` — W4 (`guards_ok`): opens the fast path to odd `N` and
   per-axis type-II, adds the `sigma > 1` guard. The behaviour-flipping
   commit; lands with the `tests/nplan.c` decline-test inversions.
6. `0d52d0ae` — the `_t210` 2D reference cases run through the fast path.
7. `5e8ae5de` — direct and fast paths pinned separately across odd,
   type-II, and the reference grid.
8. This closeout: full precision matrix, the CodSpeed walltime comparison,
   and these doc updates, plus the `guards_ok` header comment correction
   (scoping the sinc-power band-edge clause to even `N`) and a float-only
   DECONV spike-test tolerance fix (`tests/nfast.c`, see the test plan note
   below).

`.claude/skills/understanding-the-planner-api/reference/solvers-problems-windows.md`
carried the stale "declines type-II NDFT, odd `N`" / "the fast path declines
them" wording and stale `nfast.c` path references from before this work.
Fixed in `ee22da38` (guard-set wording and the two "declines" sentences) and
`6e74a52b` (the remaining `kernel/nfft/nfast.c` -> `nfft-nd.c` path
references in `solvers-problems-windows.md` and
`building-testing-examples.md`).

## 9. Files touched

```
kernel/deconv/deconv-1d.c              W1
kernel/deconv/deconv-2d.c              W2
kernel/deconv/deconv-3d.c              W2
kernel/deconv/deconv-nd.c              W3
kernel/nfft/nfft-nd.c                  W4 (guards_ok + header comment)
tests/nfast.c                          new DECONV rank tests, 2D _t210 roster,
                                        float-precision spike tolerance fix
tests/nplan.c                          decline-test inversion, odd/type-II fast passes
tests/nplan_data.c                     NO_DIRECT mirror pass
docs/plans/nfft-fast-typeii-odd-n.md   this file, closeout
CONTEXT.md                             fast-path geometry glossary entry
docs/reviews/feature-planner-review.md section 6 items marked fixed
```

No new source files, so no `Makefile.am` / `CMakeLists.txt` change.
No public API change: `nfft3.h`, `iplanner.h` untouched.
