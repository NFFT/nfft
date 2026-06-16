# Phase B — correctness net

- **Outcome:** net pinned ✅
- **Target:** `X(trafo_direct)()` — kernel/nfft/nfft.c:145 (1d branch)
- **Fault injected:** imaginary-kernel sign flip in the 1d inner loop (line 163),
  `(COS(omega) - II * SIN(omega))` → `(COS(omega) + II * SIN(omega))` — `artifacts/fault.diff`
- **Suite to run in inner loop:** `nfft` (run in all three precisions in Phase D)
- **Net size:** 116 cases (21 `nfft` test-data files: `nfft_1d_{2,4,10,20,50}_{1,10,20,50}` + `nfft_online`)

## Correctness net

<!-- Correctness net (canonical format). The fault flips the 1d `trafo_direct` cases of the
     `nfft` suite; 2d/3d and all adjoint cases stay green (fault is 1d-only). Pinned in double;
     the same suite/cases are the net in float and long double. -->
| suite | case | error | bound |
|-------|------|-------|-------|
| nfft  | nfft_1d_2_1.txt … m=8, trafo_direct   | 7.21e-02 | 1.07e-14 |
| nfft  | nfft_1d_50_50.txt … m=8, trafo_direct | (FAIL)   | 1.07e-14 |
| nfft  | nfft_online … trafo_direct (32 cases) | (FAIL)   | 1.07e-14 |

(116 `-> FAIL` lines total across 21 files; small N=2,4 files contribute 9 sub-cases each.)

## Revert confirmation

- suite green again: **yes** (`checkall` 1854/1854 `-> OK`, exit 0) · `git diff` empty: **yes**

The net to run in the Phase D inner loop is the **`nfft`** suite (`checkall` / `checkall_threads`),
run in **all three precisions** — the recurrence must keep it green in float, double, and long double.
