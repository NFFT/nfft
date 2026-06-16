# Phase E — exit gate

- **Verdict:** <complete ✅ | reverted ⛔>
- **Target metric:** <case>: <base> → <final> ns (<Δ%>)
- **Landed change:** `artifacts/change.diff`

## Comparison (baseline → final, all cases × all precisions)

<!-- Comparison table format — keep these columns; one row per benchmark case PER PRECISION
     (prec = d/f/l). The example rows illustrate verdicts (faster vs within-noise) — replace
     with real rows. `threshold` is the shorthand `2%/3σ` = max(3·stdev_ns, 2% of base median).
     Noise rule: regressed only if median rises past threshold AND survives a re-run. A
     regression in ANY precision fails the gate. -->
| prec | case                     | base median_ns | final median_ns | Δ%   | threshold | verdict   |
|------|--------------------------|----------------|-----------------|------|-----------|-----------|
| d    | <nfft_forward_direct_1d/…> | <123456>     | <95012>         | <−23%> | 2%/3σ   | <✅ faster> |
| f    | <nfft_forward_direct_1d/…> | <98765>      | <80120>         | <−19%> | 2%/3σ   | <✅ faster> |
| l    | <nfft_forward_direct_1d/…> | <210987>     | <165430>        | <−22%> | 2%/3σ   | <✅ faster> |

## Exit conditions (all must hold, in float · double · long double)

1. Full suite passes as in Phase A (incl. `checkall_threads`), every precision: <PASS d/f/l> — `artifacts/final-tests-{d,f,l}.log`
2. No benchmark regresses beyond noise (every case, every precision): <PASS d/f/l>
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
