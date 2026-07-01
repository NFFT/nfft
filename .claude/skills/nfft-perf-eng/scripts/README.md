# Helper scripts (the deterministic spine of the loop)

Many steps in this loop are **deterministic** — building the three precision trees,
capturing the full test + benchmark state, collating the per-process JSON, applying the
noise rule, fitting an error-growth trend. Those steps are driven by the scripts here so
they are **reproducible and not re-derived by hand** each run. The judgement steps
(injecting a fault, deciding what to optimize, writing the rounding-error analysis,
reading the verdict) stay with the agent; the scripts do the mechanical work around them.

Run the shell scripts from the **repo root**; run the Python helpers with `uv` (per
[[prefer-uv-not-pip]] — they are stdlib-only, no third-party deps):

| Script | What it does (deterministic) | Used in |
|--------|------------------------------|---------|
| `perf-init.sh <slug>` | create `docs/perfeng/NNNN-<slug>/`, copy every template in, stamp the commit, print the index row | Step 0 |
| `perf-build.sh [walltime\|simulation]` | clear stale `config.h`, configure + build `build-cmake{,-f,-l}` with CI-aligned flags (use `walltime` — `simulation` only for the callgrind tie-breaker) | Phase A, E, and rebuilds in D |
| `perf-capture.sh <baseline\|final> <task-dir> [--bench-only\|--tests-only] [--prec P] [--filter RE]` | run full `ctest` + benchmark per precision into `<task-dir>/artifacts/`, collate JSON; logs loadavg + WARNs if loaded; exit 1 if any precision isn't fully green/captured. `--bench-only` re-measures benchmarks without the slow ctest | Phase A (`baseline`), Phase F (`final`), noisy-capture re-measure |
| `perf-confirm.sh <task-dir>` | re-measure the case(s) the noise rule flagged (affected precisions only) and report which survive — the scripted "re-run before believing it"; exit 1 if any survives | Phase F, when `compare` flags a regression |
| `perf-net.py table\|names\|check <log>` | parse `checkall` stdout into the **correctness-net** table / failing-case names / a green-or-diff check (exit 1 if not green) | Phase B (enumerate net + revert), Phase E (net unchanged) |
| `perf-bench.py snapshot <bench.json> [--prec d]` | emit a **Benchmark snapshot** markdown table | A / C / E deliverables |
| `perf-bench.py compare (--base B --final F \| --taskdir DIR)` | emit a **Comparison table** + noise-rule verdict; exit 1 if any case regressed beyond `max(3σ, 2%)` | Phase F exit gate |
| `perf-trend.py fit <data>` / `compare <base> <opt>` | log-log least-squares error-growth fit; `compare` flags an order-of-growth regression | Phase D analysis, differential trend study |
| `perf-summary.py charts --taskdir DIR` | generate the required inline-SVG charts — speedup + error-vs-N trend, **one per precision** (tabbed in summary.html) — from the artifacts | Phase F close-out |
| `perf-summary.py check --taskdir DIR` | verify `summary.html` links every deliverable + artifact and embeds the required charts (exit 1 on any gap) | Phase F close-out |

```bash
# typical front-to-back wiring (walltime — the one metric)
SCR=.claude/skills/nfft-perf-eng/scripts
$SCR/perf-init.sh trafo-direct                       # Step 0 -> docs/perfeng/0001-trafo-direct
$SCR/perf-build.sh walltime                          # Phase A: build all 3 precisions
$SCR/perf-capture.sh baseline docs/perfeng/0001-trafo-direct   # Phase A: capture the exit reference
# ... Phases B, C (fault / slowdown — agent-driven), D (rounding-error analysis) ...
$SCR/perf-capture.sh final docs/perfeng/0001-trafo-direct      # Phase F: re-capture
uv run python $SCR/perf-bench.py compare --taskdir docs/perfeng/0001-trafo-direct  # Phase F: verdict
```

These scripts encode the rules the prose docs explain — read the matching phase doc for
the *why*; run the script for the *what*. They never inject faults, edit the target, or
mutate the shared `docs/perfeng/README.md` index (that row is printed for you to add), so
nothing the reviewer must judge is done silently.
