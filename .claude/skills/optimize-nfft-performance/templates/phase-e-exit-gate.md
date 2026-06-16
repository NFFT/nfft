# Phase E — exit gate

- **Verdict:** <complete ✅ | reverted ⛔>
- **Target metric:** <case>: <base> → <final> ns (<Δ%>)
- **Landed change:** `artifacts/change.diff`

## Comparison (baseline → final, all cases)

<!-- Comparison table format — keep these columns; one row per benchmark case (ALL cases).
     The two example rows are illustrative verdicts (faster vs within-noise) — replace
     them with real per-case rows. `threshold` is the constant shorthand `2%/3σ` =
     max(3·stdev_ns, 2% of base median). Noise rule: a case counts as regressed only if
     its median rises past that threshold AND the rise survives a re-run. -->
| case                       | base median_ns | final median_ns | Δ%   | threshold | verdict   |
|----------------------------|----------------|-----------------|------|-----------|-----------|
| <nfft_forward_direct_1d/…> | <123456>       | <95012>         | <−23%> | 2%/3σ   | <✅ faster> |
| <nfft_adjoint_direct_2d/…> | <45000>        | <45600>         | <+1%>  | 2%/3σ   | <✅ noise>  |

## Exit conditions (all must hold)

1. Full suite passes as in Phase A (incl. `checkall_threads`): <PASS> — `artifacts/final-tests.log`
2. No benchmark regresses beyond noise (every case): <PASS>
3. Deterministic cross-check (simulation), if available: <PASS | N/A — no CodSpeed, walltime-only>
4. `git diff` is only the intended change: <PASS>

## Verdict

<!-- complete: what shipped + the speedup, in one or two lines.
     reverted: why — what was tried and what blocked the gate. -->
<...>

## Close-out

<!-- The run ends here — do ALL of these, or the exit gate is not met: -->
- tracker header **Status** = <complete | reverted>
- tracker **Outcome** one-liner filled
- tracker **Phase E** row flipped (`✅` complete / `⛔` reverted) + Exit signal set
- `docs/perfeng/README.md` index row updated to the same **Status** + **Outcome**
- **`summary.html` written** (from `../templates/summary.html`) — the human report;
  set `<body class>` to `ok` (completed) or `fail`/`partial` (reverted)
