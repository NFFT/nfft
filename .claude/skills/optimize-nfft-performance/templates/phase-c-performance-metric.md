# Phase C — performance metric

- **Outcome:** <metric pinned ✅ | blocked ⛔>
- **Benchmark filter:** `<--benchmark_filter=…>`
- **Slowdown injected:** <one-line description, e.g. body wrapped in an ×N loop> — `artifacts/slowdown.diff`

## Target baseline

<!-- Benchmark snapshot format — keep these columns. -->
| case     | median_ns | stdev_ns | rounds |
|----------|-----------|----------|--------|
| <case>   | <n>       | <n>      | <n>    |

## Metric — cases that moved under the slowdown

<!-- The cases whose median rose clearly (well beyond stdev) ARE the metric. -->
| case   | median_ns base | median_ns slowed | moves the target? |
|--------|----------------|------------------|-------------------|
| <case> | <n>            | <n>              | <yes — beyond stdev> |

## Revert confirmation

- `git diff` empty: <yes>

<!-- ✅ path: after filling the sections above, flip the tracker Phase C row to ✅. -->

<!-- ⛔ BLOCKED path: if no benchmark moves OR the tooling can't measure (no obtainable
     metric), set Outcome = blocked ⛔, DELETE the Metric table (and the Target-baseline
     table too if the tooling couldn't measure at all), and UN-COMMENT the block below.
     The loop STOPS — no Phase D.

## No metric (blocked ⛔)

- **Why no metric:** <uncovered target — slowdown moved nothing | which tool failed and how>
- **Slowdown tried:** <description> — `artifacts/slowdown.diff`
- **Resolution:** add a benchmark that covers the target, then resume — or pick another target.

Close-out: set the tracker header **Status** = `reverted`, update the
`docs/perfeng/README.md` index row (status `reverted`, one-line blocked outcome), and
write `summary.html` (from `../templates/summary.html`, `<body class>` = `fail`) — the
human report documents the no-metric block as the failure. -->
