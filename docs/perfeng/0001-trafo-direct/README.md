# Perf: `X(trafo_direct)` — reduce per-call cost of the direct NDFT

- **Target:** `X(trafo_direct)()` — kernel/nfft/nfft.c:145
- **Track:** walltime
- **Baseline commit:** a692216c
- **Status:** complete (with caveat)   <!-- in-progress | complete | reverted -->
- **Outcome:** Phase recurrence (cos/sin → one complex multiply) on the 1d forward path: **~1.3× double, ~1.6× float, ~3.9× long double**; 116-case net green in all three precisions. **Caveat:** a reproducible **+8% float regression** on the untouched `adjoint_direct_1d[32/100]` micro-benchmark (code-layout/i-cache effect) trips the exit gate — reviewer decides ship/fix/revert. Double-only would have missed it.

## Phase status

| Phase | Status | Deliverable | Exit signal |
|-------|--------|-------------|-------------|
| Preflight          | ✅ | [preflight.md](preflight.md)                         | track = walltime |
| Phase A — baseline | ✅ | [phase-a-baseline.md](phase-a-baseline.md)           | green d/f/l (1854/1843/1854); 22 cases × 3 |
| Phase B — net      | ✅ | [phase-b-correctness-net.md](phase-b-correctness-net.md) | 116-case net pinned (suite `nfft`) |
| Phase C — metric   | ✅ | [phase-c-performance-metric.md](phase-c-performance-metric.md) | metric = `nfft_forward_direct_1d[*]` (10× under slowdown) |
| Phase D — inner    | ✅ | [phase-d-inner-loop.md](phase-d-inner-loop.md)       | recurrence; net green d/f/l; 1.3×/1.6×/3.9× |
| Phase E — exit     | ✅ | [phase-e-exit-gate.md](phase-e-exit-gate.md)         | landed-with-caveat: target 1.3×/1.6×/3.9×; 1 float layout regression |

<!-- Status: ⬜ todo · 🔄 in-progress · ✅ done · ⛔ blocked. A ⛔ on any row ⇒ header Status = reverted. -->

## Log

- 2026-06-16 — created; baseline commit a692216c; Preflight: walltime (no CodSpeed MCP/CLI).
- 2026-06-16 — precision matrix enabled: build float·double·long-double (stale in-source config.h cleared first).
- 2026-06-16 — Phase A green d/f/l (1854/1843/1854, 22 bench cases ×3). Long double ~250× slower than double (COSL/SINL).
- 2026-06-16 — Phase B: 116-case net. Phase C: metric = forward_direct_1d[*].
- 2026-06-16 — Phase D: recurrence; net green d/f/l; speedup 1.3×/1.6×/3.9× (precision-dependent).
- 2026-06-16 — Phase E: tests green d/f/l; **float adjoint_1d[32/100] +8% layout regression (confirmed)**; landed-with-caveat. summary.html written.
