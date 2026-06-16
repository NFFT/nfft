# Phase C — inner loop (optimize against the scoped net + metric)

*[← Overview & map](../REFERENCE.md) · Prev: [Phase B — performance metric](phase-b-performance-metric.md) · Next: [Phase D — exit gate](phase-d-exit-gate.md)*

This is the fast loop, run against the *narrow* subset identified in A and B (not the
full suite — that is Phase D). After every change:

1. `cmake --build build-cmake -j && build-cmake/tests/checkall` (add
   `build-cmake/tests/checkall_threads` if the change touches OpenMP) — the Phase-A
   net must stay green.
2. Re-run the Phase-B benchmark case(s) (walltime) — `median_ns` should drop, and
   **must not** rise beyond noise, versus the saved baseline.

Iterate until the metric is satisfactory. The scoped checks here are *necessary but
not sufficient*: they are fast feedback, but they only see the narrow slice. The
authoritative verdict is Phase D.
