# Perf: <target> — <one-line goal>

- **Target:** `<X(func)>` — <file:line>
- **Track:** —                    <!-- set in Preflight: walltime | simulation -->
- **Baseline commit:** <sha>
- **Status:** in-progress         <!-- in-progress | complete | reverted -->
- **Outcome:** —                  <!-- one line, filled at close-out -->

## Phase status

| Phase | Status | Deliverable | Exit signal |
|-------|--------|-------------|-------------|
| Preflight          | 🔄 | [preflight.md](preflight.md)                         | — |
| Phase A — baseline | ⬜ | —                                                    | — |
| Phase B — net      | ⬜ | —                                                    | — |
| Phase C — metric   | ⬜ | —                                                    | — |
| Phase D — inner    | ⬜ | —                                                    | — |
| Phase E — exit     | ⬜ | —                                                    | — |

<!-- Status: ⬜ todo · 🔄 in-progress · ✅ done · ⛔ blocked.
     Lifecycle: ⬜ until entered, 🔄 on entry, ✅ on exit-gate pass, ⛔ on hard-gate/abandon.
     A ⛔ on any row ⇒ header Status = reverted.
     Updating a row: on ENTERING a phase set Status 🔄 and fill its Deliverable link
     (e.g. [phase-a-baseline.md](phase-a-baseline.md)); on EXIT set ✅/⛔ and fill the
     Exit signal (each phase doc names the exact string, e.g. "full suite green; N cases"). -->

## Log

- <YYYY-MM-DD> — created; baseline commit <sha>.
