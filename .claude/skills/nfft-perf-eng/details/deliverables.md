# Deliverables & the `.perfeng/` directory

*[← Overview & map](../REFERENCE.md) — cross-cutting reference; read before Phase A, consult from every phase.*

This loop is **front-to-back auditable**: every phase produces a concrete
deliverable, and a phase is not exitable until its deliverable exists in the
directory and the tracker row is updated. The tracker is the live front page; the phase
docs are the durable record of what was measured and decided.

**Deliverables are never committed.** They live in a fixed, **gitignored** directory,
`.perfeng/`. Git only ever sees your *source* changes — so commit each self-contained unit of
work as you go (a kept optimization step, a permanent test addition), and let
[Phase G](phase-g-conclude.md) squash them into one commit and carry the deliverables to the
reviewer as a **zip attached to the PR**. The run's permanent record is that commit + PR +
archive, not files in the tree.

**The rule:** *deliverable = exit gate.* You may not mark a phase `done` (or move to
the next phase) until its deliverable file is written and the tracker reflects it.
For the hard-gate phases (B, C) the deliverable also records the gate verdict — and
if the gate fails, the deliverable is the **stop report** (status `blocked`), not a
green light.

**Don't hand-write the deliverables — fill the templates.** Every deliverable below
has a fill-in skeleton under [`../templates/`](../templates/) (one per phase, plus the
tracker). At each phase, copy the matching template into `.perfeng/` and fill its
`<…>` placeholders; the canonical tables below are pre-embedded in them. This is what
keeps runs comparable. The phase docs name the exact template to copy.

## Where it lives

One fixed, **gitignored** directory — `.perfeng/` — holds the **current** run. It is a
*static* path (no `NNNN` sequence, no shared index): each run overwrites it, and Step 0's
[`scripts/perf-init.sh`](../scripts/perf-init.sh) moves any prior `.perfeng/` aside to
`.perfeng.bak/` (a one-level undo) before recreating it. Nothing here is committed. At
[conclude](phase-g-conclude.md), once the PR number `N` is known, `.perfeng/` is renamed to
`.perfeng-pr-<N>/` (no leading zeros) and zipped, so a downloaded results archive unpacks into a
PR-unique directory. (`.perfeng-pr-*` is gitignored too.)

```
.perfeng/                           # gitignored; the current optimization, front-to-back
  BASE                              # squash base — the commit + branch HEAD was at, at Step 0 (for Phase G)
  README.md                         # THE TRACKER — status per phase, links, outcome (front page)
  phase-a-baseline.md               # Phase A deliverable — baseline snapshot (the exit reference)
  phase-b-correctness-net.md        # Phase B deliverable — the pinned net (or blocked report)
  phase-c-performance-metric.md     # Phase C deliverable — the pinned metric (or blocked report)
  phase-d-error-analysis.md         # Phase D deliverable — rounding-error verdict & objective (companion to the HTML)
  error-analysis.html               # Phase D deliverable — full standard-model analysis (MathJax, self-contained)
  phase-e-inner-loop.md             # Phase E deliverable — iteration journal
  phase-f-exit-gate.md              # Phase F deliverable — final verdict
  summary.html                      # Close-out — human-facing report (any outcome)
  artifacts/                        # raw captured data, kept verbatim for exact diffing
    baseline-tests-{d,f,l}.log      # Phase A  (tee'd ctest output, one per precision)
    baseline-bench-{d,f,l}.json     # Phase A  (per-case stats, collated by perf-capture.sh)
    fault.diff                      # Phase B  (the injected correctness fault)
    fault.log                       # Phase B  (checkall stdout under the fault — parsed by perf-net.py)
    net-cases.txt                   # Phase B  (failing-case names = the net; revert + inner-loop checks)
    slowdown.diff                   # Phase C  (the injected slowdown)
    trend-{baseline,optimized}-{d,f,l}.dat  # Phase D  (differential trend study: N→error pairs, per precision, if run)
    change.diff                     # Phase E/F (the actual optimization)
    final-tests-{d,f,l}.log         # Phase F
    final-bench-{d,f,l}.json        # Phase F  (collated like baseline-bench-*.json)
    chart-speedup-{d,f,l}.svg       # Phase F  (perf-summary.py — embedded tabbed in summary.html)
    chart-trend-{d,f,l}.svg         # Phase F  (perf-summary.py — when a trend study ran, per precision)
```

Captures are **per precision** — `-d` (double), `-f` (float), `-l` (long double) — because
Phases A, D, E run the whole float·double·long-double matrix; see
[precision-matrix.md](precision-matrix.md). (`fault.diff`/`slowdown.diff`/`change.diff` are
source edits, precision-independent, so they are single files.)

**Test additions.** When the [risk assessment](risk-assessment.md) leads you to extend
the net ([extending-tests.md](extending-tests.md)), a *permanent* case lands as ordinary
source: harness/registration edits under `tests/` and any committed `tests/data/*.txt` +
generated headers — part of the change's `git diff` (commit it), not under `.perfeng/artifacts/`.
A *temporary* probe is removed before close-out; only its finding survives, in the risk table.

The benchmark binary writes one file per process to
`$CODSPEED_PROFILE_FOLDER/results/<pid>.json` (a transient scratch dir).
[`scripts/perf-capture.sh`](../scripts/perf-capture.sh) collates these into the **kept**
deliverable — one flat-array JSON per precision under `artifacts/`
(`baseline-bench-{d,f,l}.json` in Phase A, `final-bench-{d,f,l}.json` in Phase F).
[`scripts/perf-bench.py compare`](../scripts/perf-bench.py) diffs the two, per precision,
applying the noise rule.

- The phase docs mirror the skill's own `details/` phase names 1:1, so the mapping
  from instruction to deliverable is obvious.
- Keep `artifacts/` raw and verbatim (logs, JSON, diffs) — the narrative docs embed
  *summaries* (tables) of these; the raw files are what Phase F diffs against.
- Redirect commands that the phase docs write to `/tmp/...` into `.perfeng/artifacts/`
  instead, so the capture is durable. Example: `... | tee
  .perfeng/artifacts/baseline-tests.log`.

## Step 0 — open the `.perfeng/` directory (gated)

Bootstrap the directory first. This is deterministic — run
[`scripts/perf-init.sh <target-slug>`](../scripts/perf-init.sh): it creates the fixed
`.perfeng/` (with `artifacts/`), copies **every** template in (tracker→`README.md`, the
phase docs, and `error-analysis.html`), stamps the baseline commit, and records the **squash
base** — the commit + branch HEAD points at right now — into `.perfeng/BASE` so
[Phase G](phase-g-conclude.md) can collapse the run's commits deterministically. If `.perfeng/`
already exists it is moved aside to `.perfeng.bak/` and recreated fresh. What remains manual:

1. **Set the tracker's Target line** — open `.perfeng/README.md` (copied from
   [`../templates/tracker.md`](../templates/tracker.md)) and set the **Target** line; leave
   **Outcome** as `—` and **Status** as `in-progress`.

*Deliverable = exit gate:* Step 0 is not done until `.perfeng/README.md` exists with all six
phase rows present and `.perfeng/BASE` records the squash base. Every later "flip the tracker
row" instruction assumes this infrastructure already exists.

## The tracker (`README.md`)

The front page of the directory. Created in Step 0 by copying
[`../templates/tracker.md`](../templates/tracker.md), updated at every phase boundary:
flip the row's status, fill its **Deliverable** link and **Exit signal**, and append a
**Log** line. The template ships with **Phase A** `🔄` (everything below `⬜`, `Outcome`
undecided) so a fresh copy already reflects Step 0's end state — Phase A is the first phase.

Status legend: `⬜` todo · `🔄` in-progress · `✅` done · `⛔` blocked (hard gate
failed — see that phase's deliverable for why). **Lifecycle:** a row is `⬜` until
entered, `🔄` on entry, `✅` on exit-gate pass, `⛔` on a hard-gate stop or abandon.
**Header mapping:** a `⛔` on any row ⇒ header **Status** = `reverted`. When the loop
ends, set the header **Status** to `complete` (Phase F passed) or `reverted` (gave up
/ hard gate) and fill **Outcome**. **Reopen:** if Phase F bounces the work back
([phase-f-exit-gate.md](phase-f-exit-gate.md)), reset the Phase E row from `✅` to `🔄` and add a
fresh journal row — `✅` means "passed *its own* gate", not "final".

## The human summary (`summary.html`)

The markdown deliverables are the audit trail for the next agent; `summary.html` is the
report a **human reviewer** reads to approve the work. Produce it from
[`../templates/summary.html`](../templates/summary.html) at **close-out — on every
exit**: a completed run (Phase F passed), a reverted run (gave up at Phase F), or a
hard-gate block (Phase B or C). It walks the reviewer concisely from the idea (what to
change and why), through the process, to the results — and on a failure documents the
known and hypothetical causes. It also banks **patterns used/discovered** and
**carry-forward notes** so knowledge accumulates across optimization targets. Set its
`<body class>` to `ok` / `partial` / `fail` to match the outcome. In [Phase
G](phase-g-conclude.md) `summary.html` (with the rest of `.perfeng/`) is zipped and attached
to the PR — it is the reviewer's entry point into the evidence bundle.

The run's *other* human-facing HTML is the Phase-D **`error-analysis.html`** — the
self-contained, MathJax-rendered rounding-error analysis (the math the reviewer must follow
to trust the accuracy verdict). It is produced mid-run (before any speed edit), not at
close-out, and `summary.html` links to it rather than repeating the derivation.

`summary.html` must present the **entire** run to the reviewer — every phase's result, with
the numbers, **and links to every deliverable and raw artifact** (the per-phase docs,
`error-analysis.html`, and each file under `artifacts/`: tests logs, benchmark JSON, the
diffs, the net, the charts). Nothing in the directory may be orphaned;
[`scripts/perf-summary.py check`](../scripts/perf-summary.py) enforces this at close-out
(exit 1 on any unlinked file or missing required chart).

### Required visualizations (the lower bar)

`summary.html` must *show*, not only tabulate, the results a reviewer judges the run on.
Generate these deterministically with
[`scripts/perf-summary.py charts`](../scripts/perf-summary.py) — self-contained inline SVG,
no JS or network — and embed them (`<img src="artifacts/chart-*.svg">`). This is the
**minimum**; richer/extra charts are welcome.

1. **Performance — % faster per case, one chart per precision**
   (`artifacts/chart-speedup-{d,f,l}.svg`, shown **tabbed** in `summary.html` so the full case
   list stays readable). **Required on every completed run.** Base→final per case; green
   faster · red regressed · grey within-noise, with the noise threshold drawn. Tab across
   float·double·long double to see whether the win holds everywhere and whether any control
   case regressed. (A correctness-only precision has no chart — drop its tab.)
2. **Accuracy trend — error vs `N`, log-log, baseline vs optimized, one chart per precision**
   (`artifacts/chart-trend-{d,f,l}.svg`, **tabbed** in `summary.html`). When a differential
   trend study is run at all (any `improve-first` verdict, or a size-dependent accuracy risk),
   **run it in all three precisions** — float and long double can behave differently (float
   cancels/overflows sooner), so a double-only trend would hide a precision-specific accuracy
   change. Each chart plots both fitted exponents `p` (expect ≈0.5, √N) side by side — the
   clearest single picture of an order-of-growth change.
3. **Accuracy margin — measured error vs the documented bound**, per case × precision (table
   always; a bar chart optional). Show the headroom to `C·ε`; if the change touched accuracy,
   show error before→after.

On a **Phase B/C hard-gate block** there is no final benchmark/trend data, so the charts are
`N/A` — `summary.html` then carries the §"Why it stalled" analysis instead, and
`perf-summary.py check` does not demand the charts.

## Canonical formats

Use these exact shapes so snapshots are comparable across phases (Phase F diffs the
Phase-A and final benchmark snapshots; the B net is re-checked in E and F).

**Benchmark snapshot** — Phase A, C, E, F. Embed a table in the narrative doc; keep
the raw per-case JSON in `artifacts/`. The `prec` column (`d`/`f`/`l`) keeps all three
precisions in one comparable table ([precision-matrix](precision-matrix.md)).

```markdown
| prec | case                       | median_ns | stdev_ns | rounds |
|------|----------------------------|-----------|----------|--------|
| d    | nfft_forward_direct_1d/…   |    123456 |      789 |     50 |
| f    | nfft_forward_direct_1d/…   |     98765 |      654 |     50 |
| l    | nfft_forward_direct_1d/…   |    210987 |      912 |     50 |
```

**Correctness net** — Phase B. The set of cases that flip to `-> FAIL` under the
injected fault, plus the suite to run in the inner loop and the net size (the latter
two as header bullets above the table, as in the template):

```markdown
- **Suite to run in inner loop:** `nfft`
- **Net size:** 149 cases

| suite | case                                   | error   | bound    |
|-------|----------------------------------------|---------|----------|
| nfft  | nfft_1d_50_50.txt … trafo_direct       | 5.7e+14 | 1.07e-14 |
```

**Comparison table** — Phase F. Baseline vs final per case **per precision**, with the
noise rule applied and a per-case verdict. A regression in *any* precision fails the gate.

```markdown
| prec | case                     | base median_ns | final median_ns | Δ%   | threshold | verdict |
|------|--------------------------|-----------------|-----------------|------|-----------|---------|
| d    | nfft_forward_direct_1d/… |          123456 |           95012 | −23% | 2%/3σ     | ✅ faster |
| f    | nfft_forward_direct_1d/… |           98765 |           80120 | −19% | 2%/3σ     | ✅ faster |
| l    | nfft_forward_direct_1d/… |          210987 |          165430 | −22% | 2%/3σ     | ✅ faster |
```

**Iteration journal** — Phase E. One row per change attempt, appended as you go;
`net` is `green` or the B-net case that flipped to `-> FAIL`; the median columns are
the Phase-C metric case(s), before→after (medians, never single runs).

```markdown
| iter | change                         | net                              | metric median (ns) before→after |
|------|--------------------------------|----------------------------------|----------------------------------|
| 1    | hoist `K[j]` out of inner loop | green                            | 123456 → 118900                  |
| 2    | unroll ×4                      | FAIL: nfft_1d_50_50 trafo_direct | — (reverted)                     |
```

**Risk table** — Phase B (seed), D (analysis seeds), E (per-iteration rows), F (consolidated). The negative
side effects the narrow net can't see, each with its category, state, and disposition.
States: `proven` (fix/revert) · `retired` (a check that would expose it now passes) ·
`accepted` (within the documented error bound / non-accuracy only — never a material drop) ·
`residual` (suspected material drop, no constructible check; surfaced in `summary.html`).
See [risk-assessment.md](risk-assessment.md).

| risk                                        | category          | state    | evidence / disposition                          |
|---------------------------------------------|-------------------|----------|-------------------------------------------------|
| accuracy may drop for N≫ tested 1d sizes    | size-dependent    | retired  | added online check N=4096, green d/f/l (temp probe)    |
| reassociated sum loses bits in long double  | accuracy-for-speed| residual | plausible; no case isolates it — see summary    |
| shared `K[j]` helper reused by `nnfft`      | out-of-scope coupling | accepted | `nnfft` suite green in Phase F; within bound |
