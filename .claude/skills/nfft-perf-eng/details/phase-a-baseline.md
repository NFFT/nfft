# Phase A — build tree + full baseline (the exit reference)

*[← Overview & map](../REFERENCE.md) · Prev: [Step 0 — .perfeng directory](deliverables.md#step-0--open-the-perfeng-directory-gated) · Next: [Phase B — correctness net](phase-b-correctness-net.md)*

**One self-contained CMake tree per precision drives the loop** — float, double, and long
double (the sources are precision-agnostic, so a change can pass in one precision and break
another; see [precision-matrix](precision-matrix.md)). Each tree is built in **walltime**
mode (`-DNFFT_BENCHMARK_MODE=walltime`) — the one metric, the same mode CI runs on its
CodSpeed Macro Runners ([measurement-modes](measurement-modes.md)). Each tree produces the
library, the CUnit tests, **and** the benchmark binaries. (Don't use the Autotools
`--with-codspeed` path — it is legacy.) This build + the capture below is also what
**confirms the walltime harness is operable** end to end (codspeed-cpp fetched, the binary
writes its per-case JSON) — no separate preflight; if the harness can't run, it fails here.

The build is **deterministic**, so it is driven by [`scripts/perf-build.sh`](../scripts/perf-build.sh)
— it clears a stale in-source `config.h` (a leftover Autotools double config otherwise
poisons the float/long-double trees, mis-linking `nfft_*` vs `nfftf_*`; CI does the same),
then configures + builds `build-cmake{,-f,-l}` with the **CI-aligned flags**
(`-O3 … -falign-functions=64 -falign-loops=32` — the alignment pair must match
`.github/workflows/bench-linux.yml` so an untouched neighbour can't show a phantom alignment
regression; see [caveats](caveats.md)):

```bash
SCR=.claude/skills/nfft-perf-eng/scripts
$SCR/perf-build.sh walltime          # configures + builds all three precision trees
#   configures + builds build-cmake (double), build-cmake-f (float), build-cmake-l (long double)
```

(The flags, the `config.h` clearing, and the per-precision `-DNFFT_ENABLE_{FLOAT,LONG_DOUBLE}`
toggles live *inside* the script, the single source of truth — read it for the exact recipe.
Reuse one fetched codspeed-cpp checkout across trees with `PERF_CODSPEED_SRC=<path>`.)

Both signals come from each tree:

- **correctness** — `ctest --test-dir <tree>`, or run `<tree>/tests/checkall` directly for the
  granular `-> OK/FAIL` stdout (same CUnit suite as `make check`, at precision-appropriate
  tolerances — the very thing that catches a precision-specific break).
- **performance** — `<tree>/benchmarks/bench_nfft_direct[_omp]`, run directly (walltime) — the
  same metric CI gates on ([Measurement modes](measurement-modes.md)); compare medians under the
  noise rule.

(Autotools `make check` is a valid, CI-canonical correctness path — see
[`test-methodology.md`](../../../../docs/agents/test-methodology.md) — but it does not build the
benchmarks, so the loop stays in the CMake trees.)

Now record the complete state of **both** signals on each unmodified tree — the contract the
finished work is judged against in [Phase F](phase-f-exit-gate.md). It must be the *full* suite,
not the scoped subset used later, **in every precision**. The capture is deterministic too —
[`scripts/perf-capture.sh`](../scripts/perf-capture.sh) runs the full `ctest` and the full
benchmark per precision straight into the task dir's `artifacts/`, collating the per-process
codspeed JSON into one flat array per precision:

```bash
$SCR/perf-capture.sh baseline .perfeng
#   → artifacts/baseline-tests-{d,f,l}.log  and  artifacts/baseline-bench-{d,f,l}.json
#   exit 0 = every precision fully green and captured; exit 1 = a precision was not (see WARNs)
```

If **any** precision's baseline is not fully green (`perf-capture.sh` exits non-zero), **stop** —
optimization starts only from a clean tree. Keep these artifacts for the whole task.

## Deliverables (exit criteria)

`phase-a-baseline.md` is **the exit reference** — [Phase F](phase-f-exit-gate.md) is
judged against it, so it must be complete and durable (the `artifacts/` captures, not
vanishing `/tmp`). `perf-capture.sh` already wrote the durable artifacts; build the snapshot
table from them with [`scripts/perf-bench.py snapshot`](../scripts/perf-bench.py):

```bash
for p in d f l; do uv run python $SCR/perf-bench.py snapshot \
  .perfeng/artifacts/baseline-bench-$p.json --prec $p; done
```

Fill [`../templates/phase-a-baseline.md`](../templates/phase-a-baseline.md), which
records: the build config (the `perf-build.sh` walltime mode + the CI-aligned flags), the
baseline commit SHA, the full-suite test result summary **for each precision**,
and the baseline benchmark snapshot as a **Benchmark snapshot** table (canonical format with
the `prec` column — [deliverables.md](deliverables.md#canonical-formats)) covering **all** cases
in **all three precisions** (the `perf-bench.py snapshot` rows drop straight in).

`artifacts/` (verbatim, per precision): `baseline-tests-{d,f,l}.log` (tee'd `ctest`) and
`baseline-bench-{d,f,l}.json` (per-case stats, collated by the script). The narrative doc embeds
the *summary* table; these raw files are what Phase F diffs against. (If a precision can't
build/run the **benchmark** in your environment, `perf-capture.sh` captures its tests only and
WARNs — record that, never silently drop a precision; see [precision-matrix](precision-matrix.md).)

**Exit gate** — *deliverable = exit gate*: Phase A is not exitable until the baseline is
**fully green in float, double, and long double**, `phase-a-baseline.md` + the per-precision
raw artifacts exist, and the tracker Phase A row reads `✅` with exit signal `full suite green
(d/f/l); N cases captured`. A non-green baseline in *any* precision ⇒ stop; do not proceed to Phase B.

*[← Overview & map](../REFERENCE.md) · Prev: [Step 0 — .perfeng directory](deliverables.md#step-0--open-the-perfeng-directory-gated) · Next: [Phase B — correctness net](phase-b-correctness-net.md)*
