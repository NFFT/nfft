# Performance engineering log (`docs/perfeng/`)

One directory per optimization, front-to-back — the way [`docs/adr/`](../adr/) keeps
one file per decision. Each directory records a single scoped optimization (a
function, a hot loop) as it moves through the test-pinned, benchmark-measured loop,
with a concrete **deliverable** captured at every phase.

These directories are **produced by the `optimize-nfft-performance` skill**
(`.claude/skills/optimize-nfft-performance/`). Don't hand-author them; the skill
creates the directory, the tracker, and the per-phase deliverables as it runs. The
deliverable format and the exit-gate-per-phase rule live in the skill's
[`details/deliverables.md`](../../.claude/skills/optimize-nfft-performance/details/deliverables.md).

## Layout

```
docs/perfeng/
  README.md                         # this file — convention + index
  NNNN-<target-slug>/               # one optimization (e.g. 0001-trafo-direct)
    README.md                       # the tracker: phase status, links, outcome
    preflight.md  phase-a-baseline.md  phase-b-correctness-net.md
    phase-c-performance-metric.md   phase-d-inner-loop.md  phase-e-exit-gate.md
    summary.html                    # human-facing close-out report (any outcome)
    artifacts/                      # raw logs, benchmark JSON, diffs
```

The tracker (`README.md`) is the agent-facing front page; `summary.html` is the
**human report** — open it to read the run end-to-end (idea → process → results, plus
failure analysis and reusable lessons) and approve the work.

## Index

| Task | Target | Track | Status | Outcome |
|------|--------|-------|--------|---------|
| [0001-trafo-direct](0001-trafo-direct/) | `X(trafo_direct)` (kernel/nfft/nfft.c:145) | walltime | complete (caveat) | recurrence: ~1.3×/1.6×/3.9× (d/f/l), net green all 3; one confirmed +8% float layout regression on an untouched control case → reviewer decides |
