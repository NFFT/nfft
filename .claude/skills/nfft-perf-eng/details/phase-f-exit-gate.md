# Phase F — exit gate (full baseline re-check)

*[← Overview & map](../REFERENCE.md) · Prev: [Phase E — inner loop](phase-e-inner-loop.md)*

A scoped optimization can have effects outside its scope — a shared helper, a header
change, a compiler-visibility shift, an aliasing assumption, **or a precision-specific
break**. So the task is **not done** until the *entire* Phase-A baseline is re-run and
compared **in all three precisions** ([precision-matrix](precision-matrix.md)). The capture
is the deterministic [`scripts/perf-capture.sh`](../scripts/perf-capture.sh) (the same one
Phase A ran, now with the `final` phase), and the noise-rule comparison is
[`scripts/perf-bench.py compare`](../scripts/perf-bench.py) — so the pass/fail call is
reproducible, not eyeballed:

```bash
SCR=.claude/skills/nfft-perf-eng/scripts
$SCR/perf-build.sh walltime                                  # rebuild all 3 trees (or incremental if untouched)
$SCR/perf-capture.sh final docs/perfeng/0001-trafo-direct    # FULL ctest + ALL bench cases → artifacts/final-*
uv run python $SCR/perf-bench.py compare --taskdir docs/perfeng/0001-trafo-direct
#   → the Comparison table (all cases × d/f/l) + the noise-rule verdict; exit 1 if any case regressed
$SCR/perf-confirm.sh docs/perfeng/0001-trafo-direct
#   → if compare flagged any case, RE-MEASURE the affected precision(s) and report which survive
#     (a whole capture can be inflated by background load; re-measure on a quiet box). Exit 1 if
#     a flagged regression survives — then attribute it (real coupling vs code-layout artifact).
```

The full walltime pass is **slow** — ×3 for the precision matrix (the 2d/3d cases run for
seconds each, here and in Phase A). That cost is the price of catching out-of-scope and
precision-specific regressions; if you must trim it, run the affected families at full fidelity
and the rest as a coarse sanity pass, and **say so** in the report rather than silently skipping
cases or precisions.

**Exit condition (all must hold, in float · double · long double):**

1. The full test suite passes exactly as in Phase A — no new failures in **any** precision, in
   either the single-threaded or the OpenMP (`checkall_threads`) library.
2. **No benchmark regresses beyond noise** versus the Phase-A baseline — *every* case in *every*
   precision, not just the target's. `perf-bench.py compare` applies the [noise
   rule](measurement-modes.md#the-noise-rule-this-is-the-metric) (a case counts as regressed
   only if its `median_ns` rises past `max(3·stdev, 2% of the median)`) and exits non-zero if
   any case trips it; **re-run** the tripped cases (`perf-confirm.sh`) before believing it —
   noise rarely survives a second run, a real regression does. Do **not** fail the gate on raw
   walltime jitter. The target's own metric should improve (or be equal). Walltime is the metric
   CI gates on too (CodSpeed Macro Runners), so this *is* the CI check, just noisier.
3. **Any surviving regression on an untouched case is attributed** — real coupling vs a
   code-layout/cache artifact. A regression that survives the quiet re-run (condition 2) on a
   case your change *did not touch* must be explained before exit: build a one-off `simulation`
   tree and compare the case's callgrind `I refs` (the [deterministic
   tie-breaker](measurement-modes.md#deterministic-tie-breaker-optional--is-a-stuck-control-case-real-work-or-layout)).
   **`I refs` rose ⇒ real coupling — fails the gate** (fix or revert). **Same `I refs`, higher
   walltime ⇒ layout/cache artifact:** pin layout with the CI-matching alignment flags so it
   does not appear (in your run or CI — CI is walltime too, so it is *not* invisible to the
   gate); if it still persists with aligned flags, attribute and disclose it. Note the
   attribution in the report; `N/A` when no untouched case survived re-run.
4. `git diff` contains only the intended optimization.
5. **The Phase-D accuracy objective is honoured.** For an `improve-first` verdict: the named
   avoidable error source is removed, and the error behaviour is **no worse — ideally better**
   than the Phase-A baseline (a differential trend study with
   [`scripts/perf-trend.py`](../scripts/perf-trend.py) is the strongest evidence the order of
   growth improved or held). For a `clean` verdict: the **derived bound is not regressed** —
   the optimization stayed within the error model Phase D established (not merely "tests still
   pass"). A speed win that *worsened* the error behaviour Phase D flagged as fixable, or
   regressed a `clean` bound, fails this condition. **A genuine accuracy/speed trade-off Phase
   D flagged** (removing avoidable error costs speed, no joint win) is *not* the agent's to
   resolve silently: it lands only as a **`partial`** run with the trade-off, the chosen
   default, and the alternative laid out in `summary.html` for the reviewer.
6. **The risk assessment is complete and honest.** Every *material* risk in the risk table
   ([risk-assessment.md](risk-assessment.md)) — including the rows **seeded in Phase D** —
   ends as `proven`, `retired`, `accepted`, or `residual`, never unexamined or hidden. A
   `proven` **material accuracy drop is a hard no**: it does not land — fix the optimization,
   or revert (a faster target that quietly lost accuracy is not a success, and the agent never
   makes a case for shipping one). `residual` risks do *not* block landing, **but** an
   unsettled `residual` material drop forces the run to exit **`partial`** (landed with
   caveats), with the open risk and the permanent test that would settle it surfaced in
   `summary.html` for the reviewer — every change is human-reviewed before merge, so the
   summary is the case, and a surfaced caveat beats a suppressed one. If a temporary probe was
   added to settle a risk, it is removed here so condition 4 holds; its finding stays in the
   risk table. A *permanent* check addition is a legitimate part of the `git diff` — call it
   out so the reviewer expects test files. **Outcome class:** `ok` (no unsettled material
   risk) · `partial` (lands with an unsettled `residual` material risk) · `fail`/`reverted` (a
   `proven` drop that couldn't be fixed, or a hard-gate block).

> **On condition 4 and added tests:** a *temporary* probe must be gone by close-out (so
> the diff is only the optimization). A *permanent* test case added to retire a risk is
> expected *in* the diff — harness/registration edits and any committed `tests/data/*.txt`
> + generated headers. Either way, [extending-tests.md](extending-tests.md) governs how
> the case was built and the temporary-vs-permanent call.

If any check fails, the optimization is **not complete**. The agent may loop back to
Phase E with further changes and re-evaluate this gate, or — if it cannot satisfy all
of them — **give up and revert**, reporting why. A faster target bought with a
regression or a broken test elsewhere is not a success.

## Deliverables (exit criteria)

This phase **owns ending the run** — fill
[`../templates/phase-f-exit-gate.md`](../templates/phase-f-exit-gate.md), capture the
raw artifacts, and do the close-out. See [deliverables.md](deliverables.md) for layout
and canonical formats.

`phase-f-exit-gate.md` records:

- the final full re-run results (`ctest` + ALL benchmark cases) **for each precision**;
- a **Comparison table** (canonical format with the `prec` column — see
  [deliverables.md](deliverables.md#canonical-formats)) of baseline-vs-final over
  *every* benchmark case in *every* precision, with the noise rule applied and a per-case
  verdict (the `perf-bench.py compare` output drops straight in);
- the six-point exit-condition checklist above, annotated pass/fail per precision —
  including condition 5, **the Phase-D accuracy objective**;
- a **Risk assessment** — the consolidated **Risk table** (canonical format in
  [deliverables.md](deliverables.md#canonical-formats)): every material risk with its
  category, final state, and evidence/disposition. This is the source for the
  always-present Risk assessment section of `summary.html`;
- the final verdict.

Raw, verbatim under `artifacts/` (per precision): `final-tests-{d,f,l}.log`,
`final-bench-{d,f,l}.json` (collated by `perf-capture.sh` exactly as in Phase A),
and `change.diff` (the landed change). `perf-bench.py compare` diffs each
`final-bench-<p>.json` against Phase A's `baseline-bench-<p>.json` to build the per-precision
**Comparison table**. (If a precision was correctness-only at baseline, compare its tests and
say so.)

**Close-out** (part of this deliverable — the run ends here): in the tracker
(`README.md`) set header **Status** = `complete` (gate passed) or `reverted` (gave
up), fill the **Outcome** one-liner, and flip the Phase F row. Then update the
matching row in the [`docs/perfeng/README.md`](../../../../docs/perfeng/) index to the
same status + outcome. Finally, write the **human report** `summary.html` (from
[`../templates/summary.html`](../templates/summary.html), `<body class>` = `ok` or
`partial`) — the reviewer-facing walkthrough of the whole run. It must present **every
phase's result with numbers, the required charts, and links to every deliverable and raw
artifact** (see [deliverables.md](deliverables.md#required-visualizations)). Generate the
charts and verify completeness deterministically:

```bash
SCR=.claude/skills/nfft-perf-eng/scripts; T=docs/perfeng/0001-trafo-direct
uv run python $SCR/perf-summary.py charts --taskdir $T   # → artifacts/chart-speedup-{d,f,l}.svg (+ chart-trend-{d,f,l}.svg if a trend study ran)
# ... fill summary.html: embed the charts, link all deliverables + artifacts ...
uv run python $SCR/perf-summary.py check  --taskdir $T   # exit 0 only when nothing is orphaned and required charts are linked
```

*Deliverable = exit gate:* the run is not closed — Phase F not exitable — until all
six conditions hold **in every precision** (or `reverted`-with-reason is recorded) —
including the honoured Phase-D accuracy objective and a completed risk assessment with no
`proven` regression left unaddressed — `phase-f-exit-gate.md`, the per-precision artifacts +
`change.diff`, **and** `summary.html` exist, `summary.html` **embeds the required charts and
links every asset** (`perf-summary.py check` passes), **and** both the tracker header/row and
the index row are flipped. A green verdict with the tracker still `in-progress` — or green in
double but untested in float/long double, or a `summary.html` that omits charts or orphans an
artifact — is not done.
