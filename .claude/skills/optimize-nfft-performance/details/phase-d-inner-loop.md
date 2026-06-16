# Phase D — inner loop (optimize against the scoped net + metric)

*[← Overview & map](../REFERENCE.md) · Prev: [Phase C — performance metric](phase-c-performance-metric.md) · Next: [Phase E — exit gate](phase-e-exit-gate.md)*

This is the fast loop, run against the *narrow* subset identified in B and C (not the
full suite — that is Phase E). After every change:

1. `cmake --build build-cmake -j && build-cmake/tests/checkall` (add
   `build-cmake/tests/checkall_threads` if the change touches OpenMP) — the Phase-B
   net must stay green.
2. Re-run the Phase-C benchmark case(s) (walltime) — `median_ns` should drop, and
   **must not** rise beyond noise, versus the saved baseline.

Iterate until the metric is satisfactory. The scoped checks here are *necessary but
not sufficient*: they are fast feedback, but they only see the narrow slice. The
authoritative verdict is Phase E.

## Deliverables (exit criteria)

This phase produces a **living** deliverable — `phase-d-inner-loop.md`, an *iteration
journal* (fill [`../templates/phase-d-inner-loop.md`](../templates/phase-d-inner-loop.md))
— plus the current `artifacts/change.diff`. Both update **every iteration**, not once at
the end. Write into `docs/perfeng/NNNN-<target-slug>/` (worked example:
`docs/perfeng/0001-trafo-direct/`).

- **`phase-d-inner-loop.md`** — one row per change attempt, appended as you go, in the
  **Iteration journal** canonical format
  ([deliverables.md](deliverables.md#canonical-formats)): `| iter | change | net |
  metric median (ns) before→after |`. `net` = green, or the B-net case that flipped to
  `-> FAIL` (then revert that attempt); the median figures are the Phase-C metric
  case(s), before→after — never single runs.
- **`artifacts/change.diff`** — the current optimization (`git diff` of the target),
  overwritten as it evolves. It always reflects the *latest kept* state, not the journal.

*deliverable = exit gate:* Phase D is not exitable until — at the latest kept state —
the B-net is green **and** the Phase-C metric median has dropped clearly beyond noise
(the goal of the loop), the journal is current through the last iteration, and
`artifacts/change.diff` matches that state. Then flip the tracker Phase D row to ✅
(`🔄` while iterating; `⛔` only if you abandon and revert). The ✅ means "passed *its
own* gate", not "final" — if [Phase E](phase-e-exit-gate.md) bounces the work back,
reopen the row to `🔄` and add a fresh journal row. Any criterion missing ⇒ stay in C.
