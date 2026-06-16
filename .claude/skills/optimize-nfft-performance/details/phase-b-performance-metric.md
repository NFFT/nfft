# Phase B — pin the performance metric

*[← Overview & map](../REFERENCE.md) · Prev: [Phase A — correctness net](phase-a-correctness-net.md) · Next: [Phase C — inner loop](phase-c-inner-loop.md)*

**B1. Measure a baseline for the target's case(s)** with the walltime binary from
Phase 0 (run directly; it writes a per-case JSON):

```bash
CODSPEED_PROFILE_FOLDER=/tmp/b build-cmake/benchmarks/bench_nfft_direct \
    --benchmark_filter='nfft_forward_direct_1d.*'
# /tmp/b/results/<pid>.json → per case: stats.median_ns (use the median, it is robust)
```

**B2. Inject a slowdown** in the target — e.g. wrap its body in a `for` loop that
repeats the work N times. Keep results correct (so you isolate cost, not correctness).

**B3. Re-run and compare medians against the baseline.** Whichever benchmark cases'
`median_ns` rise clearly (well beyond `stdev_ns`) are the ones that exercise the
target — your performance metric. A case whose median is unchanged does **not** touch
the target.

> **HARD GATE — no obtainable metric ⇒ no progress.** The agent must be able to
> produce a *concrete, comparable* number for the target — in either **simulation**
> (callgrind instruction count) **or wall-clock** mode. Two ways this gate fails:
> (a) the measurement tooling cannot run at all (see [Tooling status](tooling-status.md)),
> or (b) the tooling runs but the deliberate slowdown moves *no* benchmark — meaning
> no benchmark exercises the target, so a real improvement would be equally
> invisible. Either way, **stop**: add a benchmark that covers the target (see
> [caveats](caveats.md)) or pick a different target. Optimization without a metric is
> unverifiable and must not proceed.

**B4. Revert** the slowdown; `git diff` must be empty.
