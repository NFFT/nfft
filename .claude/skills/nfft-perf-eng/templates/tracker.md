# Perf: <target> — <one-line goal>

- **Target:** `<X(func)>` — <file:line>
- **Baseline commit:** <sha>
- **Status:** in-progress         <!-- in-progress | complete | reverted -->
- **Outcome:** —                  <!-- one line, filled at close-out -->

## Phase status

| Phase | Status | Deliverable | Exit signal |
|-------|--------|-------------|-------------|
| Phase A — baseline       | 🔄 | [phase-a-baseline.md](phase-a-baseline.md) | — |
| Phase B — net            | ⬜ | —                                    | — |
| Phase C — metric         | ⬜ | —                                    | — |
| Phase D — error analysis | ⬜ | —                                    | — |
| Phase E — inner          | ⬜ | —                                    | — |
| Phase F — exit           | ⬜ | —                                    | — |

<!-- Status: ⬜ todo · 🔄 in-progress · ✅ done · ⛔ blocked. A fresh copy ships with Phase A
     🔄 (Step 0 just finished; Phase A is the first phase).
     Lifecycle: ⬜ until entered, 🔄 on entry, ✅ on exit-gate pass, ⛔ on hard-gate/abandon.
     A ⛔ on any row ⇒ header Status = reverted.
     Updating a row: on ENTERING a phase set Status 🔄 and fill its Deliverable link
     (e.g. [phase-a-baseline.md](phase-a-baseline.md)); on EXIT set ✅/⛔ and fill the
     Exit signal (each phase doc names the exact string, e.g. "full suite green; N cases"). -->

## Log

- <YYYY-MM-DD> — created; baseline commit <sha>.
