# Phase B — pin the correctness net

*[← Overview & map](../REFERENCE.md) · Prev: [Phase A — baseline](phase-a-baseline.md) · Next: [Phase C — performance metric](phase-c-performance-metric.md)*

**A1. Identify the target.** A specific function / region, e.g. `X(trafo_direct)()` —
the direct, O(N·M) NDFT in `kernel/nfft/nfft.c:145`.

**A2. Inject a fault.** Make the *smallest* edit that changes the target's behaviour
— flip an operator, drop a term. The goal is to make dependent tests fail, not to be
realistic. Example: in `trafo_direct`'s 1d branch, flip the sign of the imaginary
kernel — `v += f_hat[k_L] * (COS(omega) - II * SIN(omega));` →
`v += f_hat[k_L] * (COS(omega) + II * SIN(omega));`.

**A3. Rebuild and see what flips to FAIL.** This set is your correctness net.

```bash
cmake --build build-cmake -j >/dev/null            # rebuild the library
build-cmake/tests/checkall >/tmp/buggy.log 2>&1; echo "exit=$?"
grep -E '\-> (FAIL|ERROR)' /tmp/buggy.log          # granular failing cases
```

Each line names the suite/case, the measured error and the bound, e.g.
`nfft_1d_50_50.txt … trafo_direct -> FAIL 5.7e+14 ( 1.07e-14)`. The affected
**suite** (`nfft` — `trafo_direct` is NFFT-specific, so not `nfct`/`nfst`/`util`) is
the one to run in the inner loop. The detailed machine-readable report is
`tests/CUnitAutomated-Results.xml` (and `…_threads-Results.xml`).

> **HARD GATE — no failing test ⇒ no coverage ⇒ stop.** If the injected fault leaves
> *every* test green (exit 0, zero `-> FAIL`), the target is **not covered** by the
> suite. This is a blocking condition, not a green light: you cannot safely optimize
> code whose behaviour no test pins. Either add a test that fails on the fault first
> (then resume), or stop and report the coverage gap. Never optimize an uncovered
> region. (Try a *more destructive* fault before concluding there is no coverage — a
> too-subtle edit may stay within tolerance.)

**A4. Revert and re-confirm green.** Restore the exact original, rebuild, re-run;
`build-cmake/tests/checkall` must exit 0 with zero FAIL lines and `git diff` must be
empty.

You now know the precise tests that guard this region.

## Deliverables (exit criteria)

Fill [`../templates/phase-b-correctness-net.md`](../templates/phase-b-correctness-net.md)
in the task dir (`docs/perfeng/NNNN-<target-slug>/`, e.g. `0001-trafo-direct`). It has
two outcomes — record exactly one:

- **Net pinned ✅** — the fault flipped ≥1 case. The doc records: the injected fault
  (saved verbatim as `artifacts/fault.diff`); the resulting net as a **Correctness
  net** table (canonical format — see [deliverables.md](deliverables.md#canonical-formats))
  *with the suite to run in the inner loop and the net size*; and the revert
  confirmation (`git diff` empty, suite green again).
- **Blocked ⛔** — no test failed even under a more destructive fault ⇒ region
  uncovered. The doc is a **blocked report**: it documents the coverage gap and states
  that the loop STOPS here — no Phase C/D/E. `artifacts/fault.diff` still records the
  fault tried.

*Deliverable = exit gate:* Phase B is not exitable until `phase-b-correctness-net.md`
and `artifacts/fault.diff` exist AND the tracker Phase B row is flipped — `✅` (net
pinned, proceed to C) or `⛔` (blocked, stop). On `⛔` the run ends here: set the
tracker header **Status** = `reverted`, update the `docs/perfeng/README.md` index
row (status `reverted`, one-line blocked outcome), and write the human report
`summary.html` (`<body class>` = `fail`, documenting the coverage gap) — the index must
not be left showing `in-progress`, and a blocked run still gets a reviewer-facing report.
