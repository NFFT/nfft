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
| `perf-init.sh <slug>` | create the fixed, gitignored `.perfeng/`, copy every template in, stamp the baseline commit, record the squash base in `.perfeng/BASE` | Step 0 |
| `perf-build.sh [walltime\|simulation]` | clear stale `config.h`, configure + build `build-cmake{,-f,-l}` with CI-aligned flags (use `walltime` — `simulation` only for the callgrind tie-breaker) | Phase A, E, and rebuilds in D |
| `perf-capture.sh <baseline\|final> <task-dir> [--bench-only\|--tests-only] [--prec P] [--filter RE]` | run full `ctest` + benchmark per precision into `<task-dir>/artifacts/`, collate JSON; logs loadavg + WARNs if loaded; exit 1 if any precision isn't fully green/captured. `--bench-only` re-measures benchmarks without the slow ctest | Phase A (`baseline`), Phase F (`final`), noisy-capture re-measure |
| `perf-confirm.sh <task-dir>` | re-measure the case(s) the noise rule flagged (affected precisions only) and report which survive — the scripted "re-run before believing it"; exit 1 if any survives | Phase F, when `compare` flags a regression |
| `perf-net.py table\|names\|check <log>` | parse `checkall` stdout into the **correctness-net** table / failing-case names / a green-or-diff check (exit 1 if not green) | Phase B (enumerate net + revert), Phase E (net unchanged) |
| `perf-bench.py snapshot <bench.json> [--prec d]` | emit a **Benchmark snapshot** markdown table | A / C / E deliverables |
| `perf-bench.py compare (--base B --final F \| --taskdir DIR)` | emit a **Comparison table** + noise-rule verdict; exit 1 if any case regressed beyond `max(3σ, 2%)` | Phase F exit gate |
| `perf-trend.py fit <data>` / `compare <base> <opt>` | log-log least-squares error-growth fit; `compare` flags an order-of-growth regression | Phase D analysis, differential trend study |
| `perf-summary.py charts --taskdir DIR` | generate the required inline-SVG charts — speedup + error-vs-N trend, **one per precision** (tabbed in summary.html) — from the artifacts | Phase F close-out |
| `perf-summary.py check --taskdir DIR` | verify `summary.html` links every deliverable + artifact and embeds the required charts (exit 1 on any gap) | Phase F close-out |
| `perf-conclude.sh squash -m "msg"` | preflight (clean tree, feature branch, base recorded), squash the run's commits to one (`reset --soft <base>` + commit), print push/label/PR/package commands | Phase G |
| `perf-conclude.sh package <pr-number>` | rename `.perfeng/` → `.perfeng-pr-<N>/` (no leading zeros) and zip it to `perfeng-pr-<N>.zip` outside the tree (standard ZIP), ready to attach to the PR; re-run after follow-up work to refresh the archive (then overwrite the zip in the existing PR comment) | Phase G (after the PR exists; and on any post-conclude follow-up) |

```bash
# typical front-to-back wiring (walltime — the one metric)
SCR=.claude/skills/nfft-perf-eng/scripts
$SCR/perf-init.sh trafo-direct                       # Step 0 -> .perfeng/ (gitignored)
$SCR/perf-build.sh walltime                          # Phase A: build all 3 precisions
$SCR/perf-capture.sh baseline .perfeng               # Phase A: capture the exit reference
# ... Phases B, C (fault / slowdown — agent-driven), D (rounding-error analysis) ...
# ... Phase E: optimize; commit each kept change ...
$SCR/perf-capture.sh final .perfeng                  # Phase F: re-capture
uv run python $SCR/perf-bench.py compare --taskdir .perfeng     # Phase F: verdict
$SCR/perf-conclude.sh squash -m "Optimize trafo_direct: …" # Phase G: squash; then push/PR (opt-in)
$SCR/perf-conclude.sh package 231                    # Phase G: after the PR — rename+zip → attach
```

These scripts encode the rules the prose docs explain — read the matching phase doc for
the *why*; run the script for the *what*. They never inject faults or edit the target, so
nothing the reviewer must judge is done silently. `perf-conclude.sh` does rewrite history
(the squash) and creates the archive, but pushing, opening the PR, and attaching the zip stay
opt-in — it only prints those commands.
