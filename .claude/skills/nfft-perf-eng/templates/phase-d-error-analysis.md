# Phase D — rounding-error analysis (verdict & objective)

<!-- Short markdown companion to error-analysis.html (the full, MathJax-rendered derivation).
     This file is what the next agent and the Phase-F exit gate read WITHOUT rendering HTML:
     the verdict, the derived bound, the avoidable source (if any), and the seeded risks.
     The math itself lives in error-analysis.html — link, don't duplicate it here. -->

- **Verdict:** <clean | improve-first>
- **Target role:** <downstream-limited | leveraged | terminal> — <how hard to push accuracy, and why>
- **Full analysis:** [error-analysis.html](error-analysis.html)
- **Standard model `u`:** float `2⁻²⁴` · double `2⁻⁵³` · long double `2⁻⁶⁴`/`2⁻¹¹³`

## Derived bound

<!-- The forward-error bound from D2, with the N/M dependence explicit, per precision.
     Give BOTH the worst-case envelope and the expected (random-sign) behaviour. -->
- **Worst case (guarantee):** <e.g. `E∞/‖f̂‖₁ ≲ γ_{N-1} + c·u ≈ N·u` for the direct NDFT>
- **Expected (random signs):** <`~√N·u` — what random data exhibits; the N·u envelope is loose>
- **Suite bound vs derived:** <flat `C·ε` (C=<48>) robust under the √N regime; only the worst-case `N·u` envelope breaches, at N≈<…> — worst-case-only>

## Measured vs derived (baseline reconciliation)

<!-- From artifacts/baseline-tests-{d,f,l}.log: measured error vs N (and M) for the target's
     cases. Does the measured trend track the derived order of growth? If a trend study was
     run, cite the fitted exponent p (scripts/perf-trend.py) and artifacts/error-trend-vs-N.md. -->
| prec | N (M) | measured error | bound | trend p (expect ≈0.5) |
|------|-------|----------------|-------|-----------------------|
| d    | <…>   | <…>            | <…>   | <p≈0.5 √N / →1 worst-case or avoidable> |

## Objective for Phase E

<!-- clean       -> "pure performance, no accuracy regression; derived bound is the yardstick"
     improve-first -> name the AVOIDABLE source and how to remove it, THEN tune speed
     Calibrate to the target role (don't minimise blindly): protect leveraged/terminal
     accuracy; don't over-invest where downstream-limited. -->
<...>

<!-- Accuracy/speed trade-off, if any (no joint win): state it, the chosen default, and the
     alternative. This is a REVIEWER decision → forces a `partial` outcome. Delete if none. -->
- **Trade-off:** <none | removing <source> costs <X>% speed — default <chosen>, alternative <other>; reviewer to weigh>

## Seeded risks (roll into the Phase-E/F risk table)

<!-- The size/precision dependence the analysis exposes — the first risk-table rows.
     category ∈ accuracy-for-speed | size-dependent | precision-specific | input-distribution
     | out-of-scope coupling. See ../details/risk-assessment.md. -->
| risk                                   | category       | state    | evidence / disposition                |
|----------------------------------------|----------------|----------|---------------------------------------|
| <worst-case N·u envelope exceeds C·ε for N≫ tested>| <size-dependent (worst-case only)> | <accepted> | <measured p≈0.5 (√N) stays under C·ε; only the loose N·u envelope breaches, at N≈…> |

<!-- After filling the two deliverables, flip the tracker Phase D row to ✅ with
     Exit signal = the verdict ("clean" or "improve-first: <source>"). -->
