# Phase E — inner loop (optimize against the scoped net + metric)

*[← Overview & map](../REFERENCE.md) · Prev: [Phase D — rounding-error analysis](phase-d-error-analysis.md) · Next: [Phase F — exit gate](phase-f-exit-gate.md)*

**Optimize toward the [Phase-D objective](phase-d-error-analysis.md), not speed alone.** If
Phase D's verdict was `improve-first`, the inner loop's *first* job is to remove the named
avoidable error source — *then* tune speed (often the same restructuring does both). If it
was `clean`, optimize for speed while **not regressing the derived bound** — the
standard-model bound, not just test pass/fail, is now the accuracy yardstick.

This is the fast loop, run against the *narrow* subset identified in B and C (not the
full suite — that is Phase F) — but **in all three precisions** (a change that holds in
double can break float or long double; see [precision-matrix](precision-matrix.md)). The
net is narrow, so 3× is cheap and catches a precision-specific break *early*. After every
change, for each tree `build-cmake{,-f,-l}`:

1. `cmake --build <tree> -j && <tree>/tests/checkall` (add `<tree>/tests/checkall_threads`
   if the change touches OpenMP) — the Phase-B net must stay green in **every** precision.
2. Re-run the Phase-C benchmark case(s) (walltime) in each — `median_ns` should drop, and
   **must not** rise beyond noise, versus that precision's saved baseline.

```bash
SCR=.claude/skills/nfft-perf-eng/scripts
for t in build-cmake build-cmake-f build-cmake-l; do
  cmake --build $t -j >/dev/null && \
    $t/tests/checkall 2>&1 | uv run python $SCR/perf-net.py check -   # GREEN in every precision (exit 0)
  echo "$t net check -> $?"
done
```

A kept change must leave the B-net **green** (`perf-net.py check` exit 0) in every
precision; if it reports failures, that change broke the net — `perf-net.py names` tells you
which case(s) flipped, then revert that attempt. (Incremental rebuilds here — `cmake --build`
only. The full reconfigure-and-rebuild
[`scripts/perf-build.sh`](../scripts/perf-build.sh) and the full capture
[`scripts/perf-capture.sh`](../scripts/perf-capture.sh) belong to Phases A and F, not this
scoped loop.)

3. **Update the risk table** (seeded in [Phase D](phase-d-error-analysis.md)). For each kept
   change, ask what it could break that this narrow net can't see — reassociated arithmetic,
   a cheaper transcendental, a size-dependent error, an input-range assumption (the
   categories in [risk-assessment.md](risk-assessment.md)). Add a risk-table row. When a risk
   is *material* (plausible and would matter) **and** cheaply testable, **extend the net** to
   settle it ([extending-tests.md](extending-tests.md)) — typically an online fast-vs-direct
   case at a larger size, fit with [`scripts/perf-trend.py`](../scripts/perf-trend.py) when
   the question is order-of-growth — and record the outcome (`retired`, `proven`, or
   `accepted`). Otherwise carry the risk forward as `residual` for the Phase-F summary.

4. **Commit the kept change.** When an attempt sticks — net green in all three precisions and
   the metric improved — commit it as its own self-contained commit with a clear message (e.g.
   `git commit -am "trafo_direct: hoist K[j] out of the inner loop"`); a permanent test addition
   is its own commit too. `.perfeng/` is gitignored, so the commit is source-only. A reverted
   attempt leaves nothing to commit. These per-unit commits are what
   [Phase G](phase-g-conclude.md) squashes into one at the end.

Iterate until the metric is satisfactory in all three **and** the Phase-D accuracy objective
is met. The scoped checks here are *necessary but not sufficient*: they are fast feedback,
but they only see the narrow slice. The authoritative verdict is Phase F.

## Deliverables (exit criteria)

This phase produces a **living** deliverable — `phase-e-inner-loop.md`, an *iteration
journal* (fill [`../templates/phase-e-inner-loop.md`](../templates/phase-e-inner-loop.md))
— plus the current `artifacts/change.diff`. Both update **every iteration**, not once at
the end. Write into `.perfeng/` (gitignored — these deliverables are never committed).

- **`phase-e-inner-loop.md`** — one row per change attempt, appended as you go, in the
  **Iteration journal** canonical format
  ([deliverables.md](deliverables.md#canonical-formats)): `| iter | change | net |
  metric median (ns) before→after |`. `net` = green, or the B-net case that flipped to
  `-> FAIL` (then revert that attempt) — and it must report **all three precisions** (e.g.
  `green d/f/l`, or name the precision that broke); the median figures are the Phase-C
  metric case(s), before→after — never single runs (note the precision if they differ).
  - The journal also carries a **risk note** per iteration (the side effects the change
    could have that the net can't see) and, when the net was extended, the probe and its
    result. These roll up into the Phase-F risk assessment — see
    [risk-assessment.md](risk-assessment.md) and the **Risk table** format in
    [deliverables.md](deliverables.md#canonical-formats).
- **`artifacts/change.diff`** — the current optimization (`git diff` of the target),
  overwritten as it evolves. It always reflects the *latest kept* state, not the journal.

*deliverable = exit gate:* Phase E is not exitable until — at the latest kept state — the
B-net is green **in float, double, and long double**, the Phase-C metric median has dropped
clearly beyond noise (the goal of the loop), **the Phase-D accuracy objective is met** (the
avoidable source removed for `improve-first`, or the derived bound not regressed for
`clean`), the journal is current through the last iteration, and `artifacts/change.diff`
matches that state. Then flip the tracker Phase E row to ✅ (`🔄` while iterating; `⛔` only
if you abandon and revert). The ✅ means "passed *its own* gate", not "final" — if
[Phase F](phase-f-exit-gate.md) bounces the work back, reopen the row to `🔄` and add a
fresh journal row. Any criterion missing ⇒ stay in E.
