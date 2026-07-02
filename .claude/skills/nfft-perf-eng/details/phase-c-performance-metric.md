# Phase C — pin the performance metric

*[← Overview & map](../REFERENCE.md) · Prev: [Phase B — correctness net](phase-b-correctness-net.md) · Next: [Phase D — rounding-error analysis](phase-d-error-analysis.md)*

**B1. Measure a baseline for the target's case(s)** with the walltime binary from
Phase A (run directly; it writes a per-case JSON):

```bash
CODSPEED_PROFILE_FOLDER=/tmp/b build-cmake/benchmarks/bench_nfft_direct \
    --benchmark_filter='nfft_forward_direct_1d.*'
# /tmp/b/results/<pid>.json → per case: stats.median_ns (use the median, it is robust)
```

**B2. Inject a slowdown** in the target — e.g. wrap its body in a `for` loop that
repeats the work N times. Keep results correct (so you isolate cost, not correctness).

> **Make the repeated work un-foldable — the build is `-O3 -ffast-math`.** A naive ×N
> wrapper *inside* the innermost body (recomputing the same value with the same arguments,
> then discarding all but the last) is silently removed by the optimizer: pure functions
> like `COS`/`SIN` with identical arguments are **common-subexpression-eliminated** to one
> call, and an overwritten accumulator is **dead-store-eliminated** — so the benchmark
> barely moves and the gate looks (wrongly) un-pinnable. Force work the optimizer cannot
> fold: wrap a **whole loop pass** (two full passes over the k-loop, which the compiler does
> not fuse), or make every repetition *observable* (e.g. accumulate all N copies into the
> live result and divide by N — `acc += term; … result = acc/N`). Verify the slowdown
> actually landed (the metric ≈doubles for a ×2) before trusting a "no movement" reading.
> A consistent few-percent rise across *all* sizes is the tell-tale of a partly-folded
> slowdown — not a pinned metric.

**B3. Re-run and compare medians against the baseline.** Whichever benchmark cases'
`median_ns` rise clearly (well beyond `stdev_ns`) are the ones that exercise the
target — your performance metric. A case whose median is unchanged does **not** touch
the target. Collate both runs (`jq -s '[.[].benchmarks[]]' …/results/*.json > a.json`)
and let [`scripts/perf-bench.py compare`](../scripts/perf-bench.py) do the comparison
deterministically — under the *deliberate* slowdown the cases it flags `❌ regressed`
(`--base` = pre-slowdown, `--final` = slowed) are exactly your metric; here a "regression"
is the goal, so a non-zero exit confirms the metric moved.

> **HARD GATE — no obtainable metric ⇒ no progress.** The agent must be able to
> produce a *concrete, comparable* **walltime** number (`median_ns`) for the target.
> Two ways this gate fails:
> (a) the benchmark tooling cannot run at all (see [Tooling status](tooling-status.md)),
> or (b) the tooling runs but the deliberate slowdown moves *no* benchmark — meaning
> no benchmark exercises the target, so a real improvement would be equally
> invisible. Either way, **stop**: add a benchmark that covers the target (see
> [caveats](caveats.md)) or pick a different target. Optimization without a metric is
> unverifiable and must not proceed.

**B4. Revert** the slowdown; `git diff` must be empty.

## Deliverables (exit criteria)

Fill [`../templates/phase-c-performance-metric.md`](../templates/phase-c-performance-metric.md)
in `.perfeng/` and save
the injected slowdown verbatim as `artifacts/slowdown.diff`. Two outcomes — record
exactly one:

- **Metric pinned ✅** — the doc records: the target's baseline measurement as a
  [**Benchmark snapshot**](deliverables.md#canonical-formats) table (canonical format
  — do not re-define columns); the injected slowdown (saved as `artifacts/slowdown.diff`);
  the metric case(s) whose `median_ns` rose clearly, with before/after numbers; and
  the revert confirmation (`git diff` empty). Flip the tracker Phase C row to ✅.
- **Blocked ⛔** — no benchmark moves *or* tooling can't measure (the HARD GATE,
  branches (a)/(b) above). The deliverable is a **blocked report**: the doc explains
  the no-metric verdict (uncovered target / tooling can't run), `artifacts/slowdown.diff`
  is the slowdown that moved nothing. Set the tracker Phase C row to ⛔ — **the loop
  STOPS**, do not proceed (no Phase D/E/F).

*Deliverable = exit gate:* Phase C is not exitable until `phase-c-performance-metric.md`
exists, `artifacts/slowdown.diff` is saved, and the tracker row is flipped — ✅ (metric
pinned) ⇒ continue to D, or ⛔ (blocked) ⇒ the run ends here: like a Phase B block, set
the tracker header **Status** = `reverted`, write the human report
`summary.html` (`<body class>` = `fail`), then [conclude](phase-g-conclude.md) and report.
