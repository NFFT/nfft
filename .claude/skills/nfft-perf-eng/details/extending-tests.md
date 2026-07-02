# Extending the net (proving or disproving a risk)

*[← Overview & map](../REFERENCE.md) — cross-cutting reference; the response to a risk
([risk-assessment.md](risk-assessment.md)) the existing net can't settle.*

When the [risk assessment](risk-assessment.md) surfaces a **material** risk the current
net can't see — most often *"accuracy may have dropped at sizes larger than any tested
case"* — the right move is often to **add a test that would expose it**, run it, and
convert the risk from `residual` to `retired` (or to `proven`, which means fix/revert).
This is a judgement call, not a mandatory gate: extend the net when the risk is plausible
*and* matters *and* a probing check is cheap relative to the doubt it removes.

This doc is the *how*. The test harness and the two test classes it builds on are
documented in [`test-methodology.md`](../../../docs/agents/test-methodology.md) — read it
first; this doc only adds the perf-loop-specific framing.

## The two levers (glossary class names)

`CONTEXT.md` fixes the two CUnit test classes as the **file-based check** and the
**online check** (see also `test-methodology.md`). They give two ways to extend coverage,
with very different costs:

1. **Online check (fast-vs-direct) — cheap, no data file.** `setup_online` /
   `setup_adjoint_online` generate random input at a chosen size, build the reference by
   running the C **direct** transform, and compare the **fast** transform against it. To
   add one you register another case in the test harness (a new size in the online list
   for the suite) — **no `tests/data/*.txt` file is needed**. This is the primary lever
   for settling a **fast-path** optimization's risk at a larger `N`/`M`: the trusted
   direct transform *is* the oracle — **but only when that direct transform is itself
   sound** (its own error already at the floor, per the [Phase-D analysis](phase-d-error-analysis.md)).
   If you are optimizing the direct transform, it cannot also be its own oracle — use lever 2.
2. **File-based check (reference-data generator) — adds a high-precision oracle.** A
   `tests/data/*.txt` reference is produced offline by the **reference-data generator**
   (`tests/refgen/`, "refgen"), and both the direct and fast transforms are checked
   against it. Use this when you need an *independent* oracle — chiefly when the **direct**
   transform itself is the target (an online check would compare the direct transform
   against itself and prove nothing), or for a permanent regression pin at a specific
   grid. Generate with:

   ```bash
   uv run --with mpmath==1.3.0 python -m tests.refgen.generate --module <module> --precision 64
   ```

   (NFFT/NFCT/NFST only; see `test-methodology.md` for grids, options, and the
   committed-artifact workflow. Per [[prefer-uv-not-pip]], always `uv run` — never pip.)

## Which lever — by what you're optimizing

The lever differs by optimization type because the *oracle* differs: a fast-path change
can be pinned by the direct transform (an online check), but a direct-path change needs an
external (refgen) oracle.

| Target of the optimization | Cheapest settling lever |
|----------------------------|-------------------------|
| **fast path** (`trafo`, `adjoint`, `precompute_one_psi`) | online check at a larger size — direct is the oracle, no data file |
| **direct transform** (the worked example) | file-based check at a new grid (refgen) — an online check would self-compare |
| **OpenMP-only** change | the same case in `checkall_threads` (links the `_omp` lib); a serial-only case won't see it |
| **non-NFFT module / no benchmark** | the module's own suite; coverage and even the benchmark may be missing first (see [caveats.md](caveats.md)) |

## Temporary vs permanent

Decide this explicitly and record it in the risk table:

- **Temporary probe** — added solely to settle one risk for this run; **prefer an online
  check** (no data files to manage). Run it on the *baseline* (unoptimized) code to
  confirm it passes, then on the change; the delta is the evidence. **Remove it before
  close-out** so Phase F's "`git diff` is only the intended optimization" exit condition
  holds — and record in the risk table that it was a temporary probe and what it showed.
- **Permanent check** — kept in the suite and CI; this is usually where a *file-based*
  (refgen) addition lands, since generated data is awkward to add and remove for one run.
  Justified when the gap it closes is real and its run cost is acceptable (online checks at
  moderate sizes are cheap; large refgen grids add to every CI run — weigh it). A permanent
  addition is a *legitimate part of the change's `git diff`*: harness/registration edits,
  and any committed `tests/data/*.txt` + generated headers. **Confirm it green on the
  baseline commit** (it pins existing behaviour at a new condition, so it must pass
  unoptimized) — that makes it a known-green case Phase F compares against alongside the
  Phase-A baseline. Note it in the Phase-F deliverable so the reviewer expects test files.

## The settle loop (proving or disproving)

1. Build the settling check (online or refgen) at the size/distribution the risk names.
2. Confirm it **passes on the baseline** code — `git stash` the change (or check out the
   baseline commit), build, run the new case. If it fails on baseline, the check is
   mis-specified (or the bound is wrong), not evidence about your change — fix the check.
3. Restore the change, run the case again **in all three precisions**:
   - passes → risk **retired** (record it; remove if temporary).
   - fails → risk **proven** → for a material accuracy drop this is a **hard no**: fix the
     optimization or revert (Phase E / F). Never ship it.
4. Record the outcome in the risk table ([deliverables.md](deliverables.md#canonical-formats)).
   If no settling check is constructible after honest attempts, the risk stays `residual`
   → the run lands as `partial` with the gap surfaced (see [risk-assessment.md](risk-assessment.md)).

## Differential trend analysis (the strongest settle for accuracy-for-speed)

A single pass/fail at one size only shows the change is *within the bound there* — it cannot
tell you whether the optimization changed the error's **order of growth**. An optimization can
sit under the bound at every committed size yet still be heading for trouble at larger N (e.g.
a reassociation or recurrence that turns O(√N) round-off into O(N)). The check the existing
tests can't give: **the trend with respect to a parameter, compared against the unmodified code.**

So for a size-dependent / accuracy-for-speed risk, prefer to *also* run a **differential trend
study**:

1. Pick a geometric sweep of the parameter (e.g. `N = 32, 64, …, 4096`), with **one
   higher-precision oracle per point** (the [file-based check](#the-two-levers-glossary-class-names)
   — a direct-transform target needs an *external* oracle, since it cannot validate itself).
   > **The oracle must be strictly more precise than the build under test — use refgen
   > (mpmath, arbitrary `dps`), never a same-precision C reference.** Computing the
   > reference in C `long double` works for the *double* and *float* builds (quad > double >
   > float), but it is **not an oracle for the long-double build itself** — a quad reference
   > cannot out-resolve a quad transform, so the long-double error it reports is meaningless
   > (it measures the recurrence's *added* error vs an exact-trig long-double sum, not the
   > true error). Drive `tests.refgen.transforms` directly at the swept `N`/`M` with
   > `mpmath.mp.dps` set comfortably above the build's decimal digits (e.g. `dps ≥ 50` covers
   > double; `dps ≥ 40` covers quad) — the *same* arbitrary-precision oracle serves all three
   > precisions uniformly. mpmath is slow, so shrink `M` or thin the sweep if a point is
   > expensive; that is cheaper than maintaining a parallel C oracle and is correct in every
   > precision.
2. Measure the suite's error metric at each point for **both** the optimized code and the
   **baseline** (`git stash` the change between runs), **in each precision** (float and long
   double can diverge — float cancels/overflows sooner, so a double-only sweep hides a
   precision-specific change). Keep the other parameters fixed (e.g. `M`) so only the swept one
   varies. (The harness `MAX_SECONDS` cost guard skips large direct cases — raise it temporarily
   for the sweep, or shrink `M`; revert after.)
3. Fit `error ~ C·Nᵖ` (log-log least-squares) for each precision with
   [`scripts/perf-trend.py`](../scripts/perf-trend.py) and **compare the exponents** per
   precision. **The verdict is Δp (optimized − baseline), not the absolute `p`.** The absolute
   exponent often sits *above* the textbook `√N` (`p≈0.5`): a sub-dominant working-precision
   term — e.g. forming the phase `2π(k−N/2)x_j` carries absolute error `~N·u` — lifts the
   measured `p` to ≈0.6–0.7 even for the *baseline* code, which is fine and expected. What
   matters is that the optimization does not *steepen* it: a Δslope near zero (both arms on the
   same curve) ⇒ order of growth preserved ⇒ risk **retired**, with far more confidence than any
   single-size pass. A clearly steeper optimized slope (in *any* precision) ⇒ a **proven**
   order-of-growth regression (hard no), even if every individual point is still under the flat
   bound.
4. The sweep grids/oracles are a *temporary probe* (revert the grid/cost edits at close-out);
   keep the per-precision `(N, error)` data as `artifacts/trend-{baseline,optimized}-{d,f,l}.dat`
   ([`perf-summary.py charts`](../scripts/perf-summary.py) charts them tabbed in `summary.html`)
   and cite the exponents per precision in the risk table. (Worked example: the `trafo-direct`
   reference run — blocked
   recurrence vs baseline: Δp = +0.034 (double), +0.0012 (float), +0.005 (long double), all
   ≪ the ±0.15 tolerance ⇒ retired. Absolute `p` ≈ 0.62–0.68, above √N, in *both* arms — the
   working-precision phase term, not a recurrence artifact; the small Δp is what settles it.)

See [`test-methodology.md`](../../../docs/agents/test-methodology.md) for the harness
mechanics (case registration, `NN` formula, bound tables, adding a whole new transform).
