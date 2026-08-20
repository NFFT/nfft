# History & drift — what changed, and dead terms to ignore

The planner API went through several designs before this one. The original
design docs were **removed** from the tree. This page is the
distilled catalog of **superseded decisions and dead terminology**, so a future
agent who has read an old commit, an external note, or a stale mental model can
map it to the current reality. **When a doc and this page disagree, the
current-branch code wins.**

## The two biggest architectural shifts

1. **Wrapper solvers → coreless native.** Early designs *wrapped* the legacy
   `nfft.c` `trafo`/`adjoint` behind solvers sharing a "legacy core" per bundle.
   That entire model — the shared/provisional core, `set_core`, `psi_valid`
   coordination, and the `nfft_solver_direct` / `fast_1d/2d/3d/nd` wrapper roster
   — was **removed**. Today every solver is **planner-native and coreless**, and
   the fast NFFT is *decomposed* into DECONV + FFTW + CONV children.

2. **NFCT/NFST new API added, then removed.** The real-valued kinds got a full
   `plan_ng` surface, which was **rolled back before review**.
   The new API ships **NDFT + NFFT only**. The legacy `nfct_*`/`nfst_*` API is
   untouched and still the way to use those transforms.

## Superseded facts → current correction

| Stale belief (from old docs) | Current reality |
|-----------------------------------|-----------------|
| `nfct_plan_ng` / `nfst_plan_ng` exist; two-solver `{direct, fast}` roster per kind; `X(check)` guard mirror; `R*` accessors | **Removed.** New API is NDFT/NFFT only. Use legacy `nfct_*`/`nfst_*`. |
| Bundle owns a shared/provisional/legacy core; `Y(nfft_solver_plan_set_core)`; `psi_valid` | **Gone.** Coreless native plans; each reads/writes through its problem. Fast plan builds ψ in its own `awake` via CONV child. |
| Roster is `direct + fast_{1,2,3}d + fast_nd` | **Replaced** by `fast_native`, `ndft_1d`, `ndft_1d_blocked`, `ndft_nd`, `const_0d`. |
| `NFFT_NO_DIM_SPECIAL`, `NFFT_NO_FAST_WRAPPER` flags | **Removed.** Fast is gated by `NFFT_NO_FAST_NATIVE` (`1<<4`); direct by `NFFT_NO_DIRECT`/`NFFT_NO_NDFT_PLAIN`/`NFFT_NO_NDFT_BLOCKED`. |
| ψ-strategy flags `NFFT_NO_FULL_PSI`, `NFFT_SHARED_CORE`, `NFFT_PREFER_FORWARD`, `NFFT_PREFER_ADJOINT`, `NFFT_CONSERVE_MEMORY`; `fast_*_fullpsi` solvers; share-when-agree cores | **Reverted.** `PRE_PSI` (sparse) is the only fast ψ-strategy in the planner. `PRE_FULL_PSI` survives only in the *legacy* kernel, untouched. |
| `NFFT_MEASURE`/`NFFT_PATIENT`/`NFFT_EXHAUSTIVE` patience ladder; `BELIEVE_PCOST` | **Dropped.** Only `NFFT_MEASURE` (=0) and `NFFT_ESTIMATE` exist. Planning-effort control is measured-vs-estimate plus the timelimit. |
| `nfft_flags` parameter on `plan_ng_guru`; legacy `PRE_*` flags are caller inputs; `nfft_flags` in the wisdom key | **Removed.** Guru takes only `fftw_flags` + `planning`. Each solver derives its own child flags. |
| `x` is aliased/borrowed by the plan (or plan-owned and freed) | **Copied.** The plan holds a private copy; the caller keeps and frees their own `x`; the plan never writes it. |
| `plan_ng_x` / `plan_ng_f_hat` / `plan_ng_f` accessors; `x` may be `NULL` in estimate mode | **Removed.** All three arrays are required and passed to the guru in every mode; no accessors. |
| Two awake states `PLNR_SLEEPY`/`PLNR_AWAKE` | **Three:** `SLEEPY=0 < AWAKE_ZERO=1 < AWAKE=2`. `AWAKE_ZERO` is the race's measurement state. |
| Per-direction race; forward and adjoint each own a wisdom entry (`NFFT_FORWARD_ONLY`/`NFFT_ADJOINT_ONLY`) | **Forward-only race.** One plan serves both directions via `apply_adjoint`; the adjoint is not separately blessed. Those flags were dropped. |
| A separate `nfft_optimize()` verb; measurement happens after `precompute` | **Never shipped.** Measurement is at plan time (nodes arrive at the guru). There is no `nfft_optimize`. |
| The measured race zeroes `f_hat`/`f` (and/or the ψ index tables) | **No.** The race is value-blind but does *not* touch `f_hat`/`f`; ψ is *precomputed, not zeroed* (zeroing the node-derived tables would destroy the access pattern being measured). |
| Cost model uses `TICKS_PER_SECOND` / wall seconds; window is a compile-time key term | **Cycle-counter ticks** (arbitrary units) for measurement, wall clock only for the budget. Window is a **runtime ordinal**, hashed in the key. |
| Accuracy floor is a wisdom-key term / an `applicable()` gate | **Reverted.** Redesigned as a (still-unbuilt) construction-time `digits → (m, σ)` helper; not in the key, never in the cost scalar. |
| The native fast solver is 1D-only / KB-only | **Rank-general** and **all four real windows (KB, Gaussian, B-spline, sinc)**; DIRAC declined. Serial. |
| Constructor has 8–10 positional args | **12 args** including per-axis `variant` and a runtime `window` ordinal (see SKILL.md). |

## Known sharp edges still open (warn future work)

These are real, unresolved limitations — not drift:

- **`A(...)` is a no-op in release.** Every `A()`-gated invariant
  (execute-before-precompute, double-destroy, the x-restore guard) is absent in
  production builds. There is no death-test coverage. Treat as caller discipline.
- **Planning is not thread-safe**, and `_on` is non-reentrant on one plan.
  Thread count is process-global and baked into the wisdom key; no per-plan
  thread knob. OpenMP determinism is untested.
- **`table_lookup` returns the greedy min-`u` subsuming entry**, not the true
  lattice minimum for incomparable `u` sets. Dormant today (one `u` per key);
  becomes reachable if a patience ladder is ever reintroduced.
- **`uses_x` seam** (`plan_s.uses_x`) is written but not yet consumed —
  speculative generality for a future data-view/permuting solver.
- **Extensibility debt:** the NFFT `problem_*` structs and window declarations
  live in the *generic* `iplanner.h`; `kind_head[NFFT_PROBLEM_LAST]` hardcodes
  the kind universe. A first port of another transform (NFSFT/NSFFT/FPT) must do
  a header split.
- **Exported test symbols:** five `Y(plan_ng_test_*)` accessors ship as `T`
  symbols in the shared library (no `-export-symbols-regex`). Reviewed and
  accepted.
