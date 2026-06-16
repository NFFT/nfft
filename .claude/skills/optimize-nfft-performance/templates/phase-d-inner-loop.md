# Phase D — inner loop (iteration journal)

- **B-net:** `<suite>` (<N> cases) · **C-metric:** <metric case(s)>
- **Current change:** `artifacts/change.diff`

## Iterations

<!-- Iteration journal format — append one row per change attempt as you go (living doc).
     net = green, or the B-net case that flipped to -> FAIL (then revert that attempt).
     medians = the Phase-C metric case(s), before→after — never single runs. -->
| iter | change                         | net                                | metric median (ns) before→after |
|------|--------------------------------|------------------------------------|----------------------------------|
| 1    | <hoist `K[j]` out of inner loop> | <green>                          | <123456 → 118900>                |
| 2    | <unroll ×4>                    | <FAIL: nfft_1d_50_50 trafo_direct> | <— (reverted)>                   |

## Exit state

- B-net green at latest kept state: <yes>
- metric median dropped beyond noise: <base → final ns (Δ%)>
- `artifacts/change.diff` matches latest kept state: <yes>

<!-- ✅ means "passed its own gate", not "final" — if Phase E bounces the work back,
     reopen the tracker Phase D row to 🔄 and add a fresh iteration row above. -->
