# Deliverables & the task directory

*[← Overview & map](../REFERENCE.md) — cross-cutting reference; read before Phase A, consult from every phase.*

This loop is **front-to-back auditable**: every phase produces a concrete
deliverable, and a phase is not exitable until its deliverable exists in the task
directory and the tracker row is updated. Deliverables accumulate in one directory
per optimization — the same way decisions accumulate as ADRs under
[`docs/adr/`](../../../../docs/adr/). The tracker is the live front page; the phase
docs are the durable record of what was measured and decided.

**The rule:** *deliverable = exit gate.* You may not mark a phase `done` (or move to
the next phase) until its deliverable file is written and the tracker reflects it.
For the hard-gate phases (B, C) the deliverable also records the gate verdict — and
if the gate fails, the deliverable is the **stop report** (status `blocked`), not a
green light.

**Don't hand-write the deliverables — fill the templates.** Every deliverable below
has a fill-in skeleton under [`../templates/`](../templates/) (one per phase, plus the
tracker). At each phase, copy the matching template into the task dir and fill its
`<…>` placeholders; the canonical tables below are pre-embedded in them. This is what
keeps runs comparable. The phase docs name the exact template to copy.

## Where it lives

One directory per optimization, created at the very start (a Step 0, before
Preflight), under `docs/perfeng/`:

```
docs/perfeng/
  README.md                         # convention + index of all optimization tasks
  NNNN-<target-slug>/               # one optimization, front-to-back (e.g. 0001-trafo-direct)
    README.md                       # THE TRACKER — status per phase, links, outcome (front page)
    preflight.md                    # Preflight deliverable — chosen measurement track + evidence
    phase-a-baseline.md             # Phase A deliverable — baseline snapshot (the exit reference)
    phase-b-correctness-net.md      # Phase B deliverable — the pinned net (or blocked report)
    phase-c-performance-metric.md   # Phase C deliverable — the pinned metric (or blocked report)
    phase-d-inner-loop.md           # Phase D deliverable — iteration journal
    phase-e-exit-gate.md            # Phase E deliverable — final verdict
    summary.html                    # Close-out — human-facing report (any outcome)
    artifacts/                      # raw captured data, kept verbatim for exact diffing
      baseline-tests-{d,f,l}.log    # Phase A  (tee'd ctest output, one per precision)
      baseline-bench-{d,f,l}.json   # Phase A  (per-case stats, collated from the codspeed scratch dir)
      fault.diff                    # Phase B  (the injected correctness fault)
      slowdown.diff                 # Phase C  (the injected slowdown)
      change.diff                   # Phase D/E (the actual optimization)
      final-tests-{d,f,l}.log       # Phase E
      final-bench-{d,f,l}.json      # Phase E  (collated like baseline-bench-*.json)
```

Captures are **per precision** — `-d` (double), `-f` (float), `-l` (long double) — because
Phases A, D, E run the whole float·double·long-double matrix; see
[precision-matrix.md](precision-matrix.md). (`fault.diff`/`slowdown.diff`/`change.diff` are
source edits, precision-independent, so they are single files.)

The benchmark binary writes one file per process to
`$CODSPEED_PROFILE_FOLDER/results/<pid>.json`; that folder is a transient scratch dir
(point it at `/tmp`). The **kept** deliverable is a single collated JSON under
`artifacts/` — `jq -s '.' /tmp/<scratch>/results/*.json > artifacts/baseline-bench.json`
(Phase A) and likewise `final-bench.json` (Phase E). Phase E diffs these two files.

- `NNNN` is a zero-padded sequence (`0001`, `0002`, …), like ADRs; `<target-slug>`
  names the region, e.g. `0001-trafo-direct`.
- The phase docs mirror the skill's own `details/` phase names 1:1, so the mapping
  from instruction to deliverable is obvious.
- Keep `artifacts/` raw and verbatim (logs, JSON, diffs) — the narrative docs embed
  *summaries* (tables) of these; the raw files are what Phase E diffs against.
- Redirect commands that the phase docs write to `/tmp/...` into `artifacts/`
  instead, so the capture is durable. Example: `... | tee
  docs/perfeng/0001-trafo-direct/artifacts/baseline-tests.log`.

## Step 0 — open the task directory (gated)

Before Preflight, bootstrap the directory. Concretely:

1. **Pick `NNNN`** — highest existing number under `docs/perfeng/` plus one,
   zero-padded (`ls docs/perfeng/` → next is `0001`, `0002`, …).
2. **Derive `<target-slug>`** from the region (e.g. `trafo-direct`), giving
   `docs/perfeng/NNNN-<target-slug>/`.
3. **Create the tracker** — copy [`../templates/tracker.md`](../templates/tracker.md)
   into the new dir as `README.md`, set the **Target** line; leave **Track**/**Outcome**
   as `—` and the **Status** as `in-progress`. Create the empty `artifacts/` dir.
4. **Add the index row** to `docs/perfeng/README.md` (see [below](#index-in-docsperfengreadmemd)).

*Deliverable = exit gate:* Step 0 is not done until
`docs/perfeng/NNNN-<slug>/README.md` exists with all six phase rows present **and**
the `docs/perfeng/README.md` index row is added. Every later "flip the tracker row"
instruction assumes this infrastructure already exists.

## The tracker (`README.md`)

The front page of the task directory. Created in Step 0 by copying
[`../templates/tracker.md`](../templates/tracker.md), updated at every phase boundary:
flip the row's status, fill its **Deliverable** link and **Exit signal**, and append a
**Log** line. The template ships *mid-Preflight* (Preflight `🔄`, everything below `⬜`,
`Track`/`Outcome` undecided) so a fresh copy already reflects Step 0's end state.

Status legend: `⬜` todo · `🔄` in-progress · `✅` done · `⛔` blocked (hard gate
failed — see that phase's deliverable for why). **Lifecycle:** a row is `⬜` until
entered, `🔄` on entry, `✅` on exit-gate pass, `⛔` on a hard-gate stop or abandon.
**Header mapping:** a `⛔` on any row ⇒ header **Status** = `reverted`. When the loop
ends, set the header **Status** to `complete` (Phase E passed) or `reverted` (gave up
/ hard gate) and fill **Outcome**. **Reopen:** if Phase E bounces the work back
(§[phase-d](phase-e-exit-gate.md)), reset the Phase D row from `✅` to `🔄` and add a
fresh journal row — `✅` means "passed *its own* gate", not "final".

## The human summary (`summary.html`)

The markdown deliverables are the audit trail for the next agent; `summary.html` is the
report a **human reviewer** reads to approve the work. Produce it from
[`../templates/summary.html`](../templates/summary.html) at **close-out — on every
exit**: a completed run (Phase E passed), a reverted run (gave up at Phase E), or a
hard-gate block (Phase B or C). It walks the reviewer concisely from the idea (what to
change and why), through the process, to the results — and on a failure documents the
known and hypothetical causes. It also banks **patterns used/discovered** and
**carry-forward notes** so knowledge accumulates across optimization targets. Set its
`<body class>` to `ok` / `partial` / `fail` to match the outcome.

## Canonical formats

Use these exact shapes so snapshots are comparable across phases (Phase E diffs the
Phase-A and final benchmark snapshots; the B net is re-checked in D and E).

**Benchmark snapshot** — Phase A, C, D, E. Embed a table in the narrative doc; keep
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

**Comparison table** — Phase E. Baseline vs final per case **per precision**, with the
noise rule applied and a per-case verdict. A regression in *any* precision fails the gate.

```markdown
| prec | case                     | base median_ns | final median_ns | Δ%   | threshold | verdict |
|------|--------------------------|-----------------|-----------------|------|-----------|---------|
| d    | nfft_forward_direct_1d/… |          123456 |           95012 | −23% | 2%/3σ     | ✅ faster |
| f    | nfft_forward_direct_1d/… |           98765 |           80120 | −19% | 2%/3σ     | ✅ faster |
| l    | nfft_forward_direct_1d/… |          210987 |          165430 | −22% | 2%/3σ     | ✅ faster |
```

**Iteration journal** — Phase D. One row per change attempt, appended as you go;
`net` is `green` or the B-net case that flipped to `-> FAIL`; the median columns are
the Phase-C metric case(s), before→after (medians, never single runs).

```markdown
| iter | change                         | net                              | metric median (ns) before→after |
|------|--------------------------------|----------------------------------|----------------------------------|
| 1    | hoist `K[j]` out of inner loop | green                            | 123456 → 118900                  |
| 2    | unroll ×4                      | FAIL: nfft_1d_50_50 trafo_direct | — (reverted)                     |
```

## Index in `docs/perfeng/README.md`

When you create a task directory (Step 0), add one row to the index table in
`docs/perfeng/README.md`, matching its columns exactly — `Task` (a link to the task
dir), `Target`, `Track`, `Status`, `Outcome`. At Step 0 the `Track`/`Outcome` cells
are `—` (not yet known); fill `Track` after Preflight and `Status`/`Outcome` at the
close-out (Phase E, **or** a Phase B/C hard-gate `reverted` exit). This is the
at-a-glance registry of every optimization run, the way `docs/adr/` is the registry
of decisions.
