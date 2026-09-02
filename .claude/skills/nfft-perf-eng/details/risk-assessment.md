# Risk assessment (what a green run can still hide)

*[← Overview & map](../REFERENCE.md) — cross-cutting reference; applies to Phases B, D, E, F.*

The two signals this loop pins — the Phase-B correctness **net** and the Phase-C
performance **metric** — are deliberately *narrow*, and even the full Phase-A suite is
*bounded*: it checks specific transforms at specific sizes and grids (the refgen file
cases) plus the online fast-vs-direct cases. **Passing every test is authoritative only
for what the tests cover.** An optimization can be faster, leave the whole suite green,
and still have **silently dropped accuracy on a path the suite never exercises** — the
classic case is a change that only loses precision at transform sizes larger than any
tested case, or for node/coefficient distributions no case spans.

So an optimization is not characterised by its speedup alone. **Every run also produces a
risk assessment**: a deliberate survey of the negative side effects the signals in hand
do *not* rule out — proven or not — recorded in a **risk table** and surfaced in
`summary.html` on *every* exit (completed, reverted, or hard-gate blocked). A faster
target with an un-surveyed risk is not a finished run.

## What counts as a risk

A **risk** is a plausible negative side effect of the change that the correctness /
performance signals in hand do not, or may not, detect. Common categories here:

- **accuracy-for-speed** — the change trades numerical precision for speed (reassociated
  sums, a cheaper transcendental, fused multiply-adds, `-ffast-math`-sensitive
  reorderings) in a way the tested tolerances and sizes don't expose.
- **size-dependent** — the effect only appears above the largest tested `N`/`M` (error
  that grows with the transform length; an accumulator that loses bits only after many
  terms). The decisive settle here is rarely a single larger case — it is the **error
  *trend*** vs the parameter, fit and compared against the unmodified code, so a changed
  *order of growth* (e.g. O(√N) → O(N)) is caught even when every individual size still
  passes the flat bound. See the differential trend study in
  [extending-tests.md](extending-tests.md#differential-trend-analysis-the-strongest-settle-for-accuracy-for-speed).
- **precision-specific** — residual beyond what the float·double·long-double matrix
  already catches (e.g. a path taken only under one `MANT_DIG` branch).
- **input-distribution** — the change assumes node positions or coefficient magnitudes
  the test cases don't span (clustered nodes, large dynamic range, near-boundary `x`).
- **out-of-scope coupling** — a shared helper, header, or aliasing/visibility change the
  target touches is used by callers no benchmark or test in the net exercises.

## The risk table (how it threads through the loop)

The risk table is built up as the run proceeds and consolidated at the exit gate. It is *not*
a new phase or a separate file — it lives in the Phase-D journal while iterating and in
the Phase-E deliverable at close-out, in the **Risk table** canonical format
([deliverables.md](deliverables.md#canonical-formats)): `| risk | category | state |
evidence / disposition |`.

- **Phase B** — when you pin the net, note what it visibly does *not* cover near the
  target (sizes, dims, distributions). That is the seed of the risk table.
- **Phase D** — the [rounding-error analysis](phase-d-error-analysis.md) adds the
  size/precision dependence it *derives* (e.g. a flat suite bound the standard model says
  must break above some `N`) as further seed rows, before any change is made. If it finds a
  genuine **accuracy/speed trade-off** (removing avoidable error costs speed, no joint win),
  that is a reviewer decision, not the agent's: record it with a documented default and
  surface it — it forces a `partial` outcome, never a silent choice of operating point.
- **Phase E** — for each kept change, ask what it could break that the narrow net can't
  see, and add a risk-table row. When a risk is **material** (plausible *and* would matter)
  **and cheaply testable**, extend the net to prove or disprove it
  ([extending-tests.md](extending-tests.md)); otherwise carry it forward as `residual`.
- **Phase F** — consolidate. Every material risk must end in one of:
  - `proven` → an extended check now exposes the side effect. For a **material accuracy
    drop this is a hard no**: fix the optimization, or revert it. The agent never ships a
    known material drop, and never "makes a case" for one — a faster transform that lost
    accuracy is not a success.
  - `retired` → an extended check that *would* have exposed it now passes (record the
    check; if it was a temporary probe, record that it was removed and what it showed).
  - `accepted` → shown to stay **within the documented error bound** the suite already
    enforces (so not a meaningful degradation), or a non-accuracy concern reasoned
    through — with the evidence. *Never* a material accuracy drop.
  - `residual` → a *suspected* material drop the agent could **not construct a check to
    settle**, after exhausting reasonable extension attempts. Surfaced loudly in
    `summary.html` as an open risk **and** a test-coverage gap, with the permanent test
    that would settle it proposed for the reviewer.

## Why this discipline exists

A material accuracy drop can only reach this point **because the existing tests were
insufficient** — had a case exercised the path, Phase B/D/E would already show it as a
`-> FAIL`. So the risk assessment's real job is to *suspect* the uncaught drop, and
[extending-tests.md](extending-tests.md)'s job is to *settle* it. The two are one loop:
suspect → extend → expose (`proven` → fix/revert) or disprove (`retired`).

## Disposition rules

- Every change is **human-reviewed before merge**; the agent's deliverables (above all
  `summary.html`) are the *case for the reviewer*, not an autonomous merge. The one
  unforgivable outcome is a **hidden or unsurveyed** risk — a drop presented as a pure
  win. Surfacing always beats suppressing.
- A `proven` material accuracy drop is a **hard no** — fix or revert; it never lands.
- `residual` risks **do not block** landing, but they are mandatory in the summary and
  they change the outcome class: a run that lands with an unsettled `residual`
  material-accuracy risk exits **`partial`** (landed with caveats — `summary.html`
  `<body class="partial">`), never `ok`. "We didn't look" and "we looked and couldn't
  settle it" are different; only the second is acceptable, and only as `partial`.
- Prefer **retiring** a material risk by extending the net over carrying it as
  `residual`. `residual` is the option of last resort, reached only after honest attempts
  to build a settling check failed.
- Outcome classes: `ok` (every material risk retired/within-bound, no unsettled
  residual) · `partial` (lands, but an unsettled residual material risk the reviewer must
  weigh) · `fail` (reverted — incl. a `proven` drop that couldn't be fixed — or a Phase
  B/C hard-gate block).

See also: [extending-tests.md](extending-tests.md) (how to build the settling check),
[precision-matrix.md](precision-matrix.md) (the one risk class already systematised),
[caveats.md](caveats.md) (the narrow-net and narrow-benchmark gotchas this generalises).
