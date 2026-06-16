# Phase D — exit gate (full baseline re-check)

*[← Overview & map](../REFERENCE.md) · Prev: [Phase C — inner loop](phase-c-inner-loop.md)*

A scoped optimization can have effects outside its scope — a shared helper, a header
change, a compiler-visibility shift, an aliasing assumption. So the task is **not
done** until the *entire* Phase-0 baseline is re-run and compared:

```bash
cmake --build build-cmake -j
ctest --test-dir build-cmake                            # FULL suite — must be all-pass
CODSPEED_PROFILE_FOLDER=/tmp/bench-final build-cmake/benchmarks/bench_nfft_direct \
    2>/tmp/final-bench.log                              # ALL cases (walltime) vs Phase-0
```

The full walltime pass is **slow** (the 2d/3d cases run for seconds each, twice — here
and in Phase 0). That cost is the price of catching out-of-scope regressions; if you
must trim it, run the affected families at full fidelity and the rest as a coarse
sanity pass, and **say so** in the report rather than silently skipping cases.

**Exit condition (all must hold):**

1. The full test suite passes exactly as in Phase 0 — no new failures, in either the
   single-threaded or the OpenMP (`checkall_threads`) library.
2. **No benchmark regresses beyond noise** versus the Phase-0 baseline — *every* case,
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
   there is no simulation baseline from Phase 0 (it captured walltime), so reconstruct
   a rough one — `git stash` the change, build a `simulation` tree, `valgrind` the
   affected cases *one at a time* (process-total `I refs`), `git stash pop`, rebuild,
   re-measure, compare deltas. *No simulation at all:* **skip** this step and note in
   the report that the verdict rests on walltime only and CI is the final authority.
4. `git diff` contains only the intended optimization.

If any check fails, the optimization is **not complete**. The agent may loop back to
Phase C with further changes and re-evaluate this gate, or — if it cannot satisfy all
of them — **give up and revert**, reporting why. A faster target bought with a
regression or a broken test elsewhere is not a success.
