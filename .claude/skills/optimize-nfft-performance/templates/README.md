# Deliverable templates

Fill-in skeletons for the per-optimization deliverables, so every run under
`docs/perfeng/` records the same fields in the same shape and runs stay comparable
(no reinventing the structure each time). The authoritative format spec — layout,
canonical tables, the exit-gate rule — is
[`../details/deliverables.md`](../details/deliverables.md); these files are its
copy-paste counterparts.

## How to use

At each phase, **copy** the matching file into the task directory
`docs/perfeng/NNNN-<target-slug>/`, keeping the filename, and fill it in:

| When | Copy | → becomes |
|------|------|-----------|
| Step 0    | `tracker.md`                    | the task dir `README.md` (front page) |
| Preflight | `preflight.md`                  | `preflight.md` |
| Phase A   | `phase-a-baseline.md`           | `phase-a-baseline.md` |
| Phase B   | `phase-b-correctness-net.md`    | `phase-b-correctness-net.md` |
| Phase C   | `phase-c-performance-metric.md` | `phase-c-performance-metric.md` |
| Phase D   | `phase-d-inner-loop.md`         | `phase-d-inner-loop.md` |
| Phase E   | `phase-e-exit-gate.md`          | `phase-e-exit-gate.md` |
| Close-out | `summary.html`                  | `summary.html` |

(`tracker.md` is the one renamed on copy — to `README.md`; the rest keep their name.)

`summary.html` is the **human-facing report** — a reviewer reads it to approve the work.
Produce it at close-out, whether the run **completed** (Phase E passed) or was
**reverted/blocked** (a Phase B/C hard gate, or give-up). It is the one deliverable
meant for a person, not the next agent; the markdown docs are its evidence base. Set
`<body class>` to `ok` / `partial` / `fail` to match the outcome.

## Conventions

- `<angle brackets>` — replace with your value, **brackets and all** (they mark the
  slot, they are not part of the value).
- `<!-- HTML comment -->` — guidance; **delete** it once the section is filled. A few
  comments hold a *blocked-path* section to **un-comment** instead — they say so.
- Tables — keep the header row, replace the example row(s), add as many rows as needed.
  Where two example rows show distinct outcomes (e.g. faster vs within-noise), they
  illustrate the range — replace both with real rows.
- Don't pad: a field with nothing to say gets `—`, not invented prose.
- The embedded tables are the canonical formats (Benchmark snapshot, Correctness net,
  Iteration journal, Comparison table) — keep their columns exactly so runs compare.
