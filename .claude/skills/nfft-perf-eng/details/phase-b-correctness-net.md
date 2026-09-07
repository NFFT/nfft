# Phase B — pin the correctness net

*[← Overview & map](../REFERENCE.md) · Prev: [Phase A — baseline](phase-a-baseline.md) · Next: [Phase C — performance metric](phase-c-performance-metric.md)*

**A1. Identify the target.** A specific function / region, e.g. `X(trafo_direct)()` —
the direct, O(N·M) NDFT in `kernel/nfft/nfft.c:145`.

**A2. Inject a fault.** Make the *smallest* edit that changes the target's behaviour
— flip an operator, drop a term. The goal is to make dependent tests fail, not to be
realistic. Example: in `trafo_direct`'s 1d branch, flip the sign of the imaginary
kernel — `v += f_hat[k_L] * (COS(omega) - II * SIN(omega));` →
`v += f_hat[k_L] * (COS(omega) + II * SIN(omega));`.

**A3. Rebuild and see what flips to FAIL.** This set is your correctness net. The fault was
yours to design (it depends on the code); turning the resulting `checkall` stdout *into the
net* is deterministic, so [`scripts/perf-net.py`](../scripts/perf-net.py) does it — the
canonical table and the failing-case name set, no hand-grepping:

```bash
SCR=.claude/skills/nfft-perf-eng/scripts
D=.perfeng/artifacts
cmake --build build-cmake -j >/dev/null                                  # rebuild the library
build-cmake/tests/checkall > $D/fault.log 2>&1; echo "exit=$?"
uv run python $SCR/perf-net.py table $D/fault.log                        # the Correctness-net table (suite|case|error|bound)
uv run python $SCR/perf-net.py names $D/fault.log > $D/net-cases.txt     # the failing-case set (revert + inner-loop checks)
```

Each FAIL line names the case (`<file> … <trafo delegate>`), the measured error and the
bound, e.g. `nfft_1d_50_50.txt … trafo_direct -> FAIL 5.7e+14 ( 1.07e-14)`. `perf-net.py`
guesses the **suite** from the filename prefix (`nfft` — `trafo_direct` is NFFT-specific, so
not `nfct`/`nfst`/`util`); that suite is the one to run in the inner loop. The detailed
machine-readable report is `tests/CUnitAutomated-Results.xml` (and `…_threads-Results.xml`).

> **HARD GATE — no failing test ⇒ no coverage ⇒ stop.** If the injected fault leaves
> *every* test green (`perf-net.py check` reports GREEN, zero `-> FAIL`), the target is **not covered** by the
> suite. This is a blocking condition, not a green light: you cannot safely optimize
> code whose behaviour no test pins. Either add a test that fails on the fault first
> (then resume) — see [extending-tests.md](extending-tests.md) for how to add an online
> or refgen case — or stop and report the coverage gap. Never optimize an uncovered
> region. (Try a *more destructive* fault before concluding there is no coverage — a
> too-subtle edit may stay within tolerance.)

**A4. Revert and re-confirm green.** Restore the exact original, rebuild, re-run;
`perf-net.py check` must report **GREEN** (exit 0) and `git diff` must be empty:

```bash
cmake --build build-cmake -j >/dev/null
build-cmake/tests/checkall 2>&1 | uv run python $SCR/perf-net.py check -   # GREEN, exit 0
git diff --quiet && echo "diff empty"
```

You now know the precise tests that guard this region.

**Seed the risk table.** While the net is fresh, note what it visibly does *not* cover
near the target — sizes, dims, node/coefficient distributions. Those gaps are the first
entries the [risk assessment](risk-assessment.md) carries forward into Phase D (the
[rounding-error analysis](phase-d-error-analysis.md)) and the inner loop.

## Deliverables (exit criteria)

Fill [`../templates/phase-b-correctness-net.md`](../templates/phase-b-correctness-net.md)
in `.perfeng/`. It has two outcomes — record exactly one:

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
tracker header **Status** = `reverted` and write the human report
`summary.html` (`<body class>` = `fail`, documenting the coverage gap) — a blocked run
still gets a reviewer-facing report, then [conclude](phase-g-conclude.md) (there is no
optimization to land, so a PR is optional — offer one only if a permanent test closed the gap).
