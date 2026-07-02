# Phase D — rounding-error analysis (set the accuracy objective)

*[← Overview & map](../REFERENCE.md) · Prev: [Phase C — performance metric](phase-c-performance-metric.md) · Next: [Phase E — inner loop](phase-e-inner-loop.md)*

The two gates are green: the target is **covered** (B) and **measurable** (C). Before
changing a single line **for speed**, scrutinise the target's *numerical behaviour* with a
meticulous **rounding-error analysis under the standard model**. The optimization is then
held to a sharper standard than "tests still pass": a faster transform that is also *less
accurate* is usually not progress, and — read the other way — the analysis often reveals an
**avoidable** error source whose removal makes the code both **more accurate and faster**.

This phase makes no optimization edits (the tree is clean again after C's revert). It only
*reads, derives, and measures*, and emits a verdict that **sets the objective** of the
inner loop:

- **`clean`** — the target's error is at (or provably near) the best achievable floor for
  what it computes. Objective for Phase E: **pure performance**, with *no accuracy
  regression* — and now the derived bound, not just the pass/fail tests, is the yardstick.
- **`improve-first`** — the analysis and/or the baseline test errors expose **error beyond
  the floor**: an *avoidable* dependence on the transform size `N`/`M`, catastrophic
  cancellation, a poorly conditioned intermediate, an unstable summation/recurrence.
  Objective for Phase E: **first root out the avoidable error, then tune speed.** Often the
  same restructuring does both ("often enough you win on both metrics").

## D1. The standard model

State and use the standard model of floating-point arithmetic (unit roundoff `u`):

> for `∘ ∈ {+,−,×,÷}`, `fl(a ∘ b) = (a ∘ b)(1 + δ)`, `|δ| ≤ u`, and likewise
> `fl(√a) = √a (1 + δ)`. For a length-`n` accumulation the errors compose into the
> Higham constant `γₙ = n·u / (1 − n·u) ≈ n·u`.

`u` per build precision (the library's `MANT_DIG` branches): float `2⁻²⁴`, double `2⁻⁵³`,
long double `2⁻⁶⁴` (80-bit) / `2⁻¹¹³` (quad). The analysis is **precision-agnostic** in
form (a `γ`-bound in `n` and `u`) and instantiates per precision — the same reason the loop
measures all three (see [precision-matrix](precision-matrix.md)).

## D2. Derive the target's error bound (theory)

Write the target as a sequence of floating-point operations and bound the forward error,
**keeping the dependence on `N` and `M` explicit** — that dependence is the whole point.

Worked example — the direct NDFT `f_j = Σ_{k} f̂_k e^{-2πi k x_j}` (`X(trafo_direct)`):

1. Each summand is a product `f̂_k · w` with `w = e^{-2πi k x_j}` formed by `COS`/`SIN`
   (or, after an optimization, a recurrence). Bound the **per-summand** error first: the
   transcendental's relative error, plus argument-reduction error in `k·x_j`, plus the
   complex multiply.
2. The length-`N` **sequential (recursive) summation** contributes the Higham factor
   `γ_{N-1}`. Give **both** bounds, because they differ sharply for this transform's data:
   - **Worst case (guarantee):** `E∞/‖f̂‖₁ ≤ γ_{N-1} + c·u ≈ (N + c)·u` (measured as the suite
     does — see [`test-methodology.md`](../../../../docs/agents/test-methodology.md); the
     normalisation works because `|e^{-2πi k x_j}| = 1`, so `Σ|f̂_k·w_k| = ‖f̂‖₁`).
   - **Expected (random signs):** the test inputs have random signs, so the rounding errors
     partially cancel — the RMS error grows like `~√N·u`, the full `N·u` being a pessimistic
     envelope rarely realised. `√N·u` is the *dominant* term random data exhibits. (In a
     measured trend the fitted exponent often lands a little **above** 0.5 — ≈0.6–0.7 — because
     a sub-dominant working-precision term, e.g. forming the phase `2π(k−N/2)x_j` to absolute
     error `~N·u`, lifts it; this shows up in the *baseline* too and is not an alarm. See D3.)
3. Compare with the suite's **bound**: the direct transforms use a *flat* round-off floor
   `bound = C·ε` (`C` = 48 NFFT / 120 NFCT / 130 NFST). Under the **√N** regime the test data
   exhibits, a flat constant `C` is robust across a wide size range — *not* a fragile
   coincidence. Only the **worst-case** `N·u` envelope would eventually exceed `C·ε`; flag
   that as a *worst-case-only* `size-dependent` note (D3 decides from the measured trend
   whether it is a real risk) — never raise the alarm on the `N·u` envelope alone when the
   data sits on `√N` ([risk-assessment](risk-assessment.md)).

Separate **inherent** from **avoidable**: summing `N` terms is fundamentally `O(N·u)`
(mitigable by compensated summation, never zero); an implementation artifact — a recurrence
that compounds to *worse* than `N·u`, a cancelling reformulation, an un-reduced large
argument — is **avoidable** and is exactly what `improve-first` targets.

## D2′. Calibrate to where the error matters (don't minimise blindly)

The objective is **not** "drive the target's error to the smallest possible". The right
amount of accuracy depends on the target's role in the larger computation — situate it:

- **Downstream-limited.** If a *later* stage deliberately approximates (a truncated window,
  a fixed cutoff `m`), the end-to-end accuracy is capped there. Pushing the target below that
  cap with heavy machinery (compensated, or **doubly**-compensated summation, exact products)
  buys nothing the user sees — and may cost speed. Match the target's error to the cap, no
  finer.
- **Leveraged / amortised.** If the target's result is *reused* — above all a **precomputation**
  feeding many transforms — its error propagates into every consumer, and the cost of computing
  it more accurately is paid **once**. Here a smaller error is highly desirable even at a real
  one-off cost: weigh per-call vs one-off explicitly.
- **Terminal.** If the target *is* (part of) the returned result, its error is the user's
  error; the derived bound is the thing to protect.

State which case the target is in — it sets *how hard* `improve-first` should push, and it is
what makes an accuracy/speed trade-off a **reasoned** call rather than a reflex.

## D3. Cross-check against the baseline (measurement)

Theory predicts an *order of growth*; the Phase-A baseline already measured the *actual*
error. Reconcile them:

- The Phase-A `ctest` logs (`artifacts/baseline-tests-{d,f,l}.log`) print each case's
  measured error and its bound. Tabulate **measured error vs `N` (and `M`)** for the
  target's cases, per precision. Expect the trend dominated by **`√N`** (random-sign
  cancellation), well below the worst-case `N·u` envelope — that is the *healthy* baseline.
  A fitted exponent a little above 0.5 (≈0.6–0.7) is normal — the sub-dominant
  working-precision phase term (above), present in the baseline too; the `trafo-direct`
  reference run measured
  ≈0.62–0.68 in *both* arms. What would be alarming is the **optimized** trend growing
  *faster than the baseline* (an avoidable source the theory pinned, or the worst case being
  realised) — i.e. a non-trivial Δp, not the absolute `p`.
- When the regime is in question, run the **differential trend study in all three precisions**
  (float and long double can diverge — float cancels/overflows sooner, so a double-only trend
  hides a precision-specific change): a geometric sweep of `N` against an oracle **strictly more
  precise than the build** — `tests.refgen.transforms` at high `mpmath` `dps`, *not* a
  same-precision C reference (a quad C oracle cannot validate the long-double build; see
  [extending-tests.md](extending-tests.md#differential-trend-analysis-the-strongest-settle-for-accuracy-for-speed)).
  Fit `error ≈ C·Nᵖ` with [`scripts/perf-trend.py`](../scripts/perf-trend.py) and **compare the
  exponents** per precision — the verdict is `Δp = p_opt − p_base` (≈0 ⇒ preserved), not the
  absolute `p`. The data
  lands as `artifacts/trend-{baseline,optimized}-{d,f,l}.dat` and is charted tabbed (one per
  precision) in `summary.html` by [`scripts/perf-summary.py`](../scripts/perf-summary.py). This is the empirical half of the analysis and the same
  tool Phase E/the risk loop use to settle accuracy-for-speed risks — see
  [extending-tests.md](extending-tests.md#differential-trend-analysis-the-strongest-settle-for-accuracy-for-speed).

## D4. Verdict — classify and set the objective

Conclude with one of `clean` / `improve-first`, the derived bound (per precision), the
measured-vs-derived reconciliation, the **role classification** from D2′ (downstream-limited
/ leveraged / terminal), and — for `improve-first` — the **named avoidable source** and how
the optimization should remove it. This verdict is what [Phase E](phase-e-inner-loop.md)
optimizes toward and what [Phase F](phase-f-exit-gate.md) checks was honoured. It also
**seeds the risk table** with the size/precision-dependence the analysis exposed, so the
inner loop carries it forward ([risk-assessment](risk-assessment.md)).

> The honest outcome may be **`clean`: the current error behaviour is the best achievable**
> for what the target computes — that is a valid, common result, not a failure to find work.
> Say so plainly and optimize for speed alone (without regressing accuracy).

**When accuracy and speed genuinely conflict, escalate — don't decide silently.** The happy
path is one restructuring that improves both. But a *real* trade-off exists — removing the
avoidable error **costs** speed with no joint win (e.g. doubly-compensated summation on a
target whose accuracy is anyway downstream-limited, or a more accurate but slower
precomputation amortised differently than you'd guess). That is a **human-judgement call**,
not the agent's to make unilaterally: the loop is autonomous, but this decision is the
reviewer's. Pick a documented default (the role classification from D2′ is the argument —
protect leveraged/terminal accuracy, don't over-invest in downstream-limited), record the
trade-off and the alternative explicitly in `error-analysis.html`, and carry it to
[Phase F](phase-f-exit-gate.md) so the run surfaces it (a `partial` outcome with the trade-off
laid out) rather than burying a chosen point on the accuracy/speed curve.

## Deliverables (exit criteria)

Two deliverables, both in `.perfeng/`:

- **`error-analysis.html`** — the **primary** deliverable: a self-contained, MathJax-rendered
  document carrying the full standard-model derivation (the math the reviewer must follow),
  the measured-vs-derived reconciliation, and the verdict. Copy
  [`../templates/error-analysis.html`](../templates/error-analysis.html) and fill it; it is
  self-contained (MathJax via CDN `<script>`, inline CSS) so it opens in any browser. The
  math belongs here, where it renders — not in markdown.
- **`phase-d-error-analysis.md`** — a short markdown companion (fill
  [`../templates/phase-d-error-analysis.md`](../templates/phase-d-error-analysis.md)) that
  records, for the next agent and the exit gate *without rendering HTML*: the verdict
  (`clean` | `improve-first`), the derived bound per precision, the named avoidable source
  (if any) and the objective it sets, and the seeded risk-table rows. It links to the HTML.

If the trend study produced data, keep it under `artifacts/` (e.g.
`artifacts/error-trend-vs-N.md`) as in [extending-tests.md](extending-tests.md); the sweep's
grid/cost edits are a temporary probe, reverted before close-out.

*Deliverable = exit gate:* Phase D is not exitable until `error-analysis.html` **and**
`phase-d-error-analysis.md` exist with a stated verdict, the risk table is seeded, and the
tracker Phase D row is `✅` with its **Exit signal** = the verdict (`clean` or
`improve-first: <source>`). The verdict is mandatory — there is no "skip the analysis" path;
`clean` is a conclusion you reach, not the absence of one.

*[← Overview & map](../REFERENCE.md) · Prev: [Phase C — performance metric](phase-c-performance-metric.md) · Next: [Phase E — inner loop](phase-e-inner-loop.md)*
