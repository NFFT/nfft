# Phase B — correctness net

- **Outcome:** <net pinned ✅ | blocked ⛔>
- **Target:** `<X(func)>` — <file:line>
- **Fault injected:** <one-line description, e.g. imaginary-sign flip in the 1d kernel> — `artifacts/fault.diff`
- **Suite to run in inner loop:** `<suite, e.g. nfft>`
- **Net size:** <N> cases

## Correctness net

<!-- Correctness net format — keep these columns; one row per case that flipped to -> FAIL.
     In the example, `…` marks the elided middle of the real case label — expand it. -->
| suite  | case                               | error     | bound      |
|--------|------------------------------------|-----------|------------|
| <nfft> | <nfft_1d_50_50.txt … trafo_direct> | <5.7e+14> | <1.07e-14> |

## Revert confirmation

- suite green again: <yes> · `git diff` empty: <yes>

<!-- ✅ path: after filling the two sections above, flip the tracker Phase B row to ✅. -->

<!-- ⛔ BLOCKED path: if no test failed even under a more destructive fault (region
     uncovered), set Outcome = blocked ⛔, DELETE the Correctness net + Revert sections
     above, and UN-COMMENT the block below. The loop STOPS — no Phase C/D/E.

## Coverage gap (blocked ⛔)

- **Most destructive fault tried:** <description> — `artifacts/fault.diff`
- **Why no test caught it:** <which tests sit near the region and why they stayed green>
- **Resolution:** add a test that fails on this fault, then resume — or escalate.

Close-out: set the tracker header **Status** = `reverted`, update the
`docs/perfeng/README.md` index row (status `reverted`, one-line blocked outcome), and
write `summary.html` (from `../templates/summary.html`, `<body class>` = `fail`) — the
human report documents the coverage gap as the failure. -->
