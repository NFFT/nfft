# Phase E — exit gate (full baseline re-check)

*[← Overview & map](../REFERENCE.md) · Prev: [Phase D — inner loop](phase-d-inner-loop.md)*

A scoped optimization can have effects outside its scope — a shared helper, a header
change, a compiler-visibility shift, an aliasing assumption. So the task is **not
done** until the *entire* Phase-A baseline is re-run and compared:

```bash
cmake --build build-cmake -j
ctest --test-dir build-cmake                            # FULL suite — must be all-pass
CODSPEED_PROFILE_FOLDER=/tmp/bench-final build-cmake/benchmarks/bench_nfft_direct \
    2>/tmp/final-bench.log                              # ALL cases (walltime) vs Phase-A
```

The full walltime pass is **slow** (the 2d/3d cases run for seconds each, twice — here
and in Phase A). That cost is the price of catching out-of-scope regressions; if you
must trim it, run the affected families at full fidelity and the rest as a coarse
sanity pass, and **say so** in the report rather than silently skipping cases.

**Exit condition (all must hold):**

1. The full test suite passes exactly as in Phase A — no new failures, in either the
   single-threaded or the OpenMP (`checkall_threads`) library.
2. **No benchmark regresses beyond noise** versus the Phase-A baseline — *every* case,
   not just the target's. Apply the noise rule from
   [Working without CodSpeed](measurement-modes.md#working-without-codspeed): a case
   counts as regressed only if its `median_ns` rises past `max(3·stdev, 2% of the
   median)`, **and** the rise survives a re-run. Do **not** fail the gate on raw
   walltime jitter — untouched code routinely swings a few percent. The target's own
   metric should improve (or be equal). *(With CodSpeed/simulation available, judge this
   on the deterministic instruction count instead — no noise rule needed.)*
3. **Deterministic cross-check (simulation), if available:** the instruction count
   should not regress either. *With a CodSpeed account:* `codspeed run` for clean
   per-case counts, or the MCP server for the base branch's CI numbers. *Account-free:*
   there is no simulation baseline from Phase A (it captured walltime), so reconstruct
   a rough one — `git stash` the change, build a `simulation` tree, `valgrind` the
   affected cases *one at a time* (process-total `I refs`), `git stash pop`, rebuild,
   re-measure, compare deltas. *No simulation at all:* **skip** this step and note in
   the report that the verdict rests on walltime only and CI is the final authority.
4. `git diff` contains only the intended optimization.

If any check fails, the optimization is **not complete**. The agent may loop back to
Phase D with further changes and re-evaluate this gate, or — if it cannot satisfy all
of them — **give up and revert**, reporting why. A faster target bought with a
regression or a broken test elsewhere is not a success.

## Deliverables (exit criteria)

This phase **owns ending the run** — fill
[`../templates/phase-e-exit-gate.md`](../templates/phase-e-exit-gate.md), capture the
raw artifacts, and do the close-out. See [deliverables.md](deliverables.md) for layout
and canonical formats.

`phase-e-exit-gate.md` records:

- the final full re-run results (`ctest` + ALL benchmark cases);
- a **Comparison table** (canonical format — see
  [deliverables.md](deliverables.md#canonical-formats)) of baseline-vs-final over
  *every* benchmark case, with the noise rule applied and a per-case verdict;
- the four-point exit-condition checklist above, annotated pass/fail each;
- the final verdict.

Raw, verbatim under `artifacts/`: `final-tests.log`, `final-bench.json` (collated
from the codspeed scratch dir exactly as in Phase A —
`jq -s '.' /tmp/bench-final/results/*.json > …/artifacts/final-bench.json`), and
`change.diff` (the landed change). Phase E diffs `final-bench.json` against Phase A's
`baseline-bench.json` to build the **Comparison table**.

**Close-out** (part of this deliverable — the run ends here): in the tracker
(`README.md`) set header **Status** = `complete` (gate passed) or `reverted` (gave
up), fill the **Outcome** one-liner, and flip the Phase E row. Then update the
matching row in the [`docs/perfeng/README.md`](../../../../docs/perfeng/) index to the
same status + outcome. Finally, write the **human report** `summary.html` (from
[`../templates/summary.html`](../templates/summary.html), `<body class>` = `ok` or
`partial`) — the reviewer-facing walkthrough of the whole run.

*Deliverable = exit gate:* the run is not closed — Phase E not exitable — until all
four conditions hold (or `reverted`-with-reason is recorded), `phase-e-exit-gate.md`,
the three artifacts, **and** `summary.html` exist, **and** both the tracker header/row
and the index row are flipped. A green verdict with the tracker still `in-progress` is
not done.
