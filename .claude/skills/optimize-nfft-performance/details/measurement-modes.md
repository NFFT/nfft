# Measurement modes (walltime locally, simulation in CI)

*[← Overview & map](../REFERENCE.md) — cross-cutting reference, consult from any phase.*

Mode is **baked in at build time** via `-DNFFT_BENCHMARK_MODE=` (`off`, the default,
builds no benchmarks). The two measuring modes play different roles:

| | **walltime** — local workhorse | **simulation** — CI gate |
|---|---|---|
| Build | `-DNFFT_BENCHMARK_MODE=walltime` | `-DNFFT_BENCHMARK_MODE=simulation` |
| Run | binary directly, offline | under `valgrind --tool=callgrind` (or in CI) |
| Metric | wall-clock ns: min/mean/median/stdev/IQR | deterministic instruction count |
| Local result | JSON at `$CODSPEED_PROFILE_FOLDER/results/<pid>.json` | callgrind `I refs` on stderr |
| Authoritative source | the local JSON | **CodSpeed** (uploaded from CI) |
| Noise | real timing noise — compare medians | none (deterministic) |

**Use walltime for the local loop.** It is the only mode that writes results locally
and offline, so Phases 0/B/C/D read its JSON. Switching mode is just a
reconfigure (`-DNFFT_BENCHMARK_MODE=…`); you do **not** need two build trees unless you
want both modes at once. Because walltime is noisy in a container, compare **medians**
over the (auto-tuned) rounds, not single runs.

**Simulation is the deterministic enhancement** (and what CI gates on). It is
*optional* — see [Working without CodSpeed](#working-without-codspeed). How to get
simulation numbers, in increasing fidelity:

1. **Raw `valgrind --tool=callgrind` on a `simulation` build** — no account needed, but
   the `I refs` it prints is a **process total** that *includes fixed startup*
   (dynamic-linker/regex setup, ~3M Ir). Filter to a *single* case and treat it as a
   rough before/after — the delta is meaningful, the ratio is diluted by startup.
2. **CodSpeed CLI (`codspeed run` / `exec`)** — gives clean per-benchmark counts
   (its instrumentation excludes startup), **but requires a CodSpeed token**: the
   runner validates auth before it will run (verified — there is no offline /
   `--skip-upload` mode in 4.17.5). So this needs a CodSpeed account.
3. **CodSpeed MCP / cloud** — the base branch's CI numbers, for cross-commit parity.
   Needs a CodSpeed account and the repo onboarded; see
   [Tooling status](tooling-status.md). Use it in
   [Phase D](phase-d-exit-gate.md) to confirm against CI.

Instruction count and wall-clock can diverge a lot (e.g. a change measured here was
2.3× fewer instructions but only 1.3× faster wall-clock — the loop became
memory/arithmetic-bound once the transcendentals went). The wall-clock figure is the
user-facing speedup; instruction count is the low-noise proxy.

## Working without CodSpeed

A coding agent may have **no CodSpeed access at all** (no account, repo not onboarded,
no MCP). That is a supported — if non-ideal — configuration: **the entire loop
runs on walltime alone.** What you lose is the deterministic metric and the CI-history
comparison; what you keep is a complete, offline measure → change → re-measure cycle.
The cost is **timing noise**, which you must manage explicitly:

- Compare **medians**, never single runs.
- Treat a case as regressed only when its median rises **beyond noise** — past, say,
  `max(3·stdev, 2% of the median)`. Identical, untouched code routinely swings several
  percent here (worst on the few-iteration 2d/3d cases), so a hard "no case may be
  slower at all" rule produces **false regressions**.
- **Re-run** any case that trips the threshold before believing it; noise rarely
  survives a second run, a real regression does.

A raw `valgrind` single-case process total (option 1 above) is a cheap deterministic
sanity check when a walltime result is ambiguous, but clean per-case instruction
counts need a CodSpeed account (the CLI is auth-gated), so on a no-CodSpeed setup
**walltime with the noise rule is the metric** — accept that and lean on re-runs.
