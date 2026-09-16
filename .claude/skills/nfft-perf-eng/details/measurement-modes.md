# Measurement: walltime everywhere (local and CI)

*[← Overview & map](../REFERENCE.md) — cross-cutting reference, consult from any phase.*

**There is one performance metric: wall-clock time.** CI measures it with **CodSpeed
Macro Runners** in walltime mode; you measure it locally with the same `walltime` build,
run directly and offline. Local and CI agree on *what* is measured, so a local result is
a faithful (if noisier) preview of the CI gate. There is no separate simulation /
instruction-count gate to satisfy, and no CodSpeed account, token, MCP, or cloud upload
in the loop.

Mode is **baked in at build time** via `-DNFFT_BENCHMARK_MODE=` (`off`, the default,
builds no benchmarks):

| | **walltime** — the metric, local and CI |
|---|---|
| Build | `-DNFFT_BENCHMARK_MODE=walltime` |
| Run | the binary directly, offline (CI runs it on a CodSpeed Macro Runner) |
| Metric | wall-clock ns: min/mean/median/stdev/IQR |
| Local result | JSON at `$CODSPEED_PROFILE_FOLDER/results/<pid>.json` |
| Noise | real timing noise — compare **medians**, re-run to confirm |

The benchmark build still uses codspeed-cpp (FetchContent fetches it; submodules
auto-recurse) — that is just the harness library and needs **no account**. The whole
loop, Phases A/C/E/F, reads the local walltime JSON.

## The noise rule (this is the metric)

Because walltime is noisy in a container, **manage the noise explicitly** — this is the
metric, so the discipline below *is* the gate, not a workaround:

- Compare **medians**, never single runs.
- Treat a case as regressed only when its median rises **beyond noise** — past, say,
  `max(3·stdev, 2% of the median)`. Identical, untouched code routinely swings several
  percent here (worst on the few-iteration 2d/3d cases), so a hard "no case may be
  slower at all" rule produces **false regressions**.
- **Re-run** any case that trips the threshold before believing it; noise rarely
  survives a second run, a real regression does. This is scripted:
  [`perf-confirm.sh <task-dir>`](../scripts/perf-confirm.sh) re-measures only the flagged
  precision(s) and reports which regressions survive. Re-measure on a **quiet** machine —
  `perf-capture.sh` logs the 1-min loadavg and WARNs when it is high, because concurrent
  builds/sweeps inflate medians (a whole capture can read +10% under load yet sit within
  noise when re-run idle). Use `perf-capture.sh … --bench-only` to re-measure without the
  slow deterministic `ctest`.

## Deterministic tie-breaker (optional) — is a stuck control case real work or layout?

When an **untouched** control case stays flagged even after a quiet re-run, the question
is *why*: real coupling to your change (extra work executed) or a **code-layout/cache
artifact** (your edit shifted the alignment of an untouched neighbour, so the same
instructions now run on a worse cache/branch boundary). Wall-clock alone cannot tell
these apart; the deterministic instruction count can:

- Build a one-off `simulation` tree (`perf-build.sh simulation`, or just one precision)
  and run the affected case under `valgrind --tool=callgrind …
  --benchmark_filter='<one case>'`. Read `I refs` from the callgrind output. (This is a
  process **total** that includes fixed startup ~3M Ir, so use it as a *before/after
  delta* on a single case, not an absolute.)
- **Same `I refs` before and after, but a higher wall-clock ⇒ layout/cache, not real
  work.** The fix is to **pin layout** with the CI-matching alignment flags
  (`-falign-functions=64 -falign-loops=32`, already in the Phase-A flag list and in
  `bench-linux.yml`) so the artifact does not appear — *in your run or in CI*. Note: CI
  is **also** walltime now, so a layout artifact is **not** invisible to the gate the way
  it would have been under an instruction-count gate — it must be mitigated (aligned
  flags) or attributed, not waved through. See [caveats.md](caveats.md) (layout regressions).
- **`I refs` rose ⇒ real extra work** — a genuine regression to attribute and fix,
  regardless of layout.

This callgrind cross-check is the **only** use of the `simulation` build mode in the
loop — a diagnostic for an ambiguous control case, not a routine metric and not a gate.
Clean per-case instruction counts (`codspeed run`) need a CodSpeed token and are **not**
part of this loop; the raw single-case `I refs` delta above needs no account.
