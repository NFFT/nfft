# Performance engineering log (`docs/perfeng/`)

One directory per optimization, front-to-back — the way [`docs/adr/`](../adr/) keeps
one file per decision. Each directory records a single scoped optimization (a
function, a hot loop) as it moves through the test-pinned, benchmark-measured loop,
with a concrete **deliverable** captured at every phase.

These directories are **produced by the `nfft-perf-eng` skill**
(`.claude/skills/nfft-perf-eng/`). Don't hand-author them; the skill
creates the directory, the tracker, and the per-phase deliverables as it runs. The
deliverable format and the exit-gate-per-phase rule live in the skill's
[`details/deliverables.md`](../../.claude/skills/nfft-perf-eng/details/deliverables.md).

## Layout

```
docs/perfeng/
  README.md                         # this file — convention + index
  NNNN-<target-slug>/               # one optimization (e.g. the `trafo-direct` task)
    README.md                       # the tracker: phase status, links, outcome
    phase-a-baseline.md  phase-b-correctness-net.md
    phase-c-performance-metric.md   phase-d-error-analysis.md
    phase-e-inner-loop.md  phase-f-exit-gate.md
    error-analysis.html             # Phase D — rounding-error analysis (MathJax, human-facing)
    summary.html                    # human-facing close-out report (charts + links to all assets)
    artifacts/                      # raw logs, benchmark JSON, diffs, inline-SVG charts
```

The tracker (`README.md`) is the agent-facing front page; `summary.html` is the
**human report** — open it to read the run end-to-end (idea → process → results, plus
failure analysis and reusable lessons) and approve the work.

## Index

| Task | Target | Status | Outcome |
|------|--------|--------|---------|
