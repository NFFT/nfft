# Phase E — inner loop (iteration journal)

- **B-net:** `<suite>` (<N> cases) · **C-metric:** <metric case(s)>
- **Phase-D objective:** <clean — speed only, don't regress the derived bound | improve-first — remove <avoidable source> first, then tune speed>
- **Current change:** `artifacts/change.diff`

## Iterations

<!-- Iteration journal format — append one row per change attempt as you go (living doc).
     net = green in ALL THREE precisions (write "green d/f/l"), or name the precision+case
     that flipped to -> FAIL (then revert that attempt). medians = the Phase-C metric case(s),
     before→after — never single runs; note the precision if it differs. -->
| iter | change                         | net                                  | metric median (ns) before→after |
|------|--------------------------------|--------------------------------------|----------------------------------|
| 1    | <hoist `K[j]` out of inner loop> | <green d/f/l>                      | <d 123456→118900; f/l similar>   |
| 2    | <unroll ×4>                    | <FAIL (f): nfft_1d_50_50 trafo_direct> | <— (reverted)>                 |

## Risks raised (rolls up to Phase F)

<!-- Risk table format — one row per side effect this change could have that the narrow
     net can't see. category ∈ accuracy-for-speed | size-dependent | precision-specific |
     input-distribution | out-of-scope coupling. state ∈ proven | retired | accepted |
     residual. When material AND cheap, extend the net (../details/extending-tests.md) and
     record the probe + result here; else carry as residual for the summary. -->
| risk                                     | category       | state    | evidence / disposition                   |
|------------------------------------------|----------------|----------|------------------------------------------|
| <accuracy may drop for N≫ tested sizes>  | <size-dependent> | <retired> | <online check N=4096 green d/f/l (temp probe)> |

## Exit state

- B-net green at latest kept state, **all precisions**: <yes — d/f/l>
- metric median dropped beyond noise: <d: base→final (Δ%); f/l similar>
- Phase-D accuracy objective met: <yes — avoidable source removed / derived bound not regressed>
- `artifacts/change.diff` matches latest kept state: <yes>
- risk table current; every material risk `proven`/`retired`/`accepted`/`residual`: <yes>

<!-- ✅ means "passed its own gate", not "final" — if Phase F bounces the work back,
     reopen the tracker Phase E row to 🔄 and add a fresh iteration row above. -->
