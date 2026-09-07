# Phase F — exit gate

- **Verdict:** <complete ✅ | reverted ⛔>
- **Target metric:** <case>: <base> → <final> ns (<Δ%>)
- **Phase-D objective:** <clean | improve-first: <source>> — honoured: <yes>
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
2. No benchmark regresses beyond noise (every case, every precision) — `perf-bench.py compare --taskdir .`: <PASS d/f/l>
3. Surviving regressions on untouched cases attributed (callgrind `I refs`: layout artifact vs real coupling): <PASS | N/A — none survived re-run>
4. `git diff` is only the intended change: <PASS>
5. Phase-D accuracy objective honoured — `improve-first`: avoidable source removed, error no worse (trend study if order-of-growth); `clean`: derived bound not regressed: <PASS>
6. Risk assessment complete & honest — every material risk (incl. Phase-D-seeded rows) `proven`/`retired`/`accepted`/`residual`; no `proven` material drop landed (fixed or reverted); any unsettled `residual` → outcome `partial`: <PASS>

## Verdict

<!-- complete: what shipped + the speedup, in one or two lines.
     reverted: why — what was tried and what blocked the gate. -->
<...>

## Risk assessment (consolidated risk table)

<!-- Risk table format — every material risk this optimization carries, with its final
     state. proven = fix/revert (must not ship). retired = a check that would expose it now
     passes (name it; note if it was a temporary probe and what it showed). accepted =
     tolerable, with rationale. residual = open — these MUST also appear in summary.html.
     See ../details/risk-assessment.md. -->
| risk                                    | category       | state    | evidence / disposition                  |
|-----------------------------------------|----------------|----------|-----------------------------------------|
| <accuracy may drop for N≫ tested sizes> | <size-dependent> | <retired> | <online check N=4096 green d/f/l (temp probe)> |
| <reassociated sum in long double>       | <accuracy-for-speed> | <residual> | <suspected; no constructible check — proposes refgen long-double case at N≥…; forces outcome `partial`> |

Note: `proven` material accuracy drop is a hard no — it never reaches this consolidated table as a *landed* state (it was fixed or the run reverted). An unsettled `residual` material drop ⇒ this run's outcome is `partial`.

## Close-out

<!-- The run ends here — do ALL of these, or the exit gate is not met: -->
- temporary probing tests removed (their findings kept in the risk table); permanent test additions noted as expected in the `git diff`
- tracker header **Status** = <complete | reverted>
- tracker **Outcome** one-liner filled
- tracker **Phase F** row flipped (`✅` complete / `⛔` reverted) + Exit signal set
- **charts generated** — `perf-summary.py charts --taskdir .` → `artifacts/chart-speedup-{d,f,l}.svg` (+ `chart-trend-{d,f,l}.svg` if a trend study ran)
- **`summary.html` written** (from `../templates/summary.html`) — the human report;
  presents every phase's result + numbers, embeds the required charts, links every deliverable + artifact;
  set `<body class>` to `ok` (completed clean) · `partial` (landed with caveats — an unsettled `residual` material risk) · `fail` (reverted / hard-gate block)
- **completeness verified** — `perf-summary.py check --taskdir .` exits 0 (nothing orphaned; required charts linked)
- **concluded (Phase G)** — `perf-conclude.sh squash -m "…"` squashed the run to one commit; then (opt-in) pushed, opened a `perf-eng`-labelled PR to `develop`, and `perf-conclude.sh package <N>` renamed `.perfeng/` → `.perfeng-pr-<N>/` + zipped it, attached to the PR
