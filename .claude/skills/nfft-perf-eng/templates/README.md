# Deliverable templates

Fill-in skeletons for the per-optimization deliverables, so every run under the gitignored
`.perfeng/` records the same fields in the same shape and runs stay comparable
(no reinventing the structure each time). The authoritative format spec — layout,
canonical tables, the exit-gate rule — is
[`../details/deliverables.md`](../details/deliverables.md); these files are its
copy-paste counterparts.

## How to use

At each phase, **copy** the matching file into the gitignored `.perfeng/` directory,
keeping the filename, and fill it in. (Step 0's
[`scripts/perf-init.sh`](../scripts/perf-init.sh) copies them all in one shot.)

| When | Copy | → becomes |
|------|------|-----------|
| Step 0    | `tracker.md`                    | the task dir `README.md` (front page) |
| Phase A   | `phase-a-baseline.md`           | `phase-a-baseline.md` |
| Phase B   | `phase-b-correctness-net.md`    | `phase-b-correctness-net.md` |
| Phase C   | `phase-c-performance-metric.md` | `phase-c-performance-metric.md` |
| Phase D   | `phase-d-error-analysis.md` + `error-analysis.html` | `phase-d-error-analysis.md` + `error-analysis.html` |
| Phase E   | `phase-e-inner-loop.md`         | `phase-e-inner-loop.md` |
| Phase F   | `phase-f-exit-gate.md`          | `phase-f-exit-gate.md` |
| Close-out | `summary.html`                  | `summary.html` |

(`tracker.md` is the one renamed on copy — to `README.md`; the rest keep their name.)

Two deliverables are **human-facing HTML**, not markdown. `error-analysis.html` (Phase D)
is the MathJax-rendered rounding-error derivation, produced *mid-run* before any speed edit;
`summary.html` is the close-out report a reviewer reads to approve the work, produced whether
the run **completed** (Phase F passed) or was **reverted/blocked** (a Phase B/C hard gate, or
give-up). Set each one's `<body class>` to match its outcome (`clean`/`improve-first` for the
analysis; `ok`/`partial`/`fail` for the summary). `summary.html` must **embed the required
charts** (`perf-summary.py charts` — speedup; error-vs-`N` trend) and **link every deliverable
and raw artifact**; `perf-summary.py check` verifies nothing in the task dir is orphaned. See
[`../details/deliverables.md`](../details/deliverables.md#required-visualizations).

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
