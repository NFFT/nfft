# Phase A — build tree + full baseline (the exit reference)

*[← Overview & map](../REFERENCE.md) · Prev: [Preflight](preflight.md) · Next: [Phase B — correctness net](phase-b-correctness-net.md)*

**One self-contained CMake tree per precision drives the loop** — float, double, and long
double (the sources are precision-agnostic, so a change can pass in one precision and break
another; see [precision-matrix](precision-matrix.md)). Each tree is built in the mode your
[preflight](preflight.md) picked — `walltime` (local track, shown here) or `simulation`
(CodSpeed track). `FetchContent` fetches/builds codspeed-cpp once; reuse it across trees with
`-DFETCHCONTENT_SOURCE_DIR_CODSPEED=<path>`. Each tree produces the library, the CUnit tests,
**and** the benchmark binaries. (Don't use the Autotools `--with-codspeed` path — it is legacy.)

```bash
# Clear a stale in-source config.h first — a leftover Autotools (double) config poisons the
# float/long-double trees (mis-link as nfft_* vs nfftf_*). CI does this; see precision-matrix.md.
[ -f include/config.h ] && mv include/config.h include/config.h.autotools.bak
# bash array — keeps the quoted CMAKE_C_FLAGS a SINGLE argument (an unquoted $FLAGS would
# split it on spaces and break configure).
FLAGS=(-DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON \
       -DCMAKE_C_FLAGS="-O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math")
cmake -S . -B build-cmake   "${FLAGS[@]}"                              # double      (suffix d)
cmake -S . -B build-cmake-f "${FLAGS[@]}" -DNFFT_ENABLE_FLOAT=ON       # float       (suffix f)
cmake -S . -B build-cmake-l "${FLAGS[@]}" -DNFFT_ENABLE_LONG_DOUBLE=ON # long double (suffix l)
for t in build-cmake build-cmake-f build-cmake-l; do cmake --build $t -j; done
```

Both signals come from each tree:

- **correctness** — `ctest --test-dir <tree>`, or run `<tree>/tests/checkall` directly for the
  granular `-> OK/FAIL` stdout (same CUnit suite as `make check`, at precision-appropriate
  tolerances — the very thing that catches a precision-specific break).
- **performance** — `<tree>/benchmarks/bench_nfft_direct[_omp]`, run directly (walltime). See
  [Measurement modes](measurement-modes.md) for when to *also* consult the CI **simulation** metric.

(Autotools `make check` is a valid, CI-canonical correctness path — see
[`test-methodology.md`](../../../../docs/agents/test-methodology.md) — but it does not build the
benchmarks, so the loop stays in the CMake trees.)

Now record the complete state of **both** signals on each unmodified tree — the contract the
finished work is judged against in [Phase E](phase-e-exit-gate.md). It must be the *full* suite,
not the scoped subset used later, **in every precision**.

```bash
for p in d:build-cmake f:build-cmake-f l:build-cmake-l; do s=${p%%:*}; t=${p#*:}
  ctest --test-dir $t 2>&1 | tee /tmp/baseline-tests-$s.log                     # full suite
  CODSPEED_PROFILE_FOLDER=/tmp/bench-base-$s $t/benchmarks/bench_nfft_direct \
      2>/tmp/baseline-bench-$s.log                       # ALL cases → /tmp/bench-base-$s/results/*.json
done
```

(*CodSpeed track:* the benchmark baseline is the base branch's CI run, read via the
MCP in Phase E — no local capture needed; still record the test baseline per precision.)

If **any** precision's baseline is not fully green, **stop** — optimization starts only from a
clean tree. Keep these artifacts for the whole task.

## Deliverables (exit criteria)

`phase-a-baseline.md` is **the exit reference** — [Phase E](phase-e-exit-gate.md) is
judged against it, so it must be complete and durable (no `/tmp` captures that vanish).
Capture the kept artifacts under `artifacts/`, collating the codspeed scratch results
into one JSON with `jq`:

```bash
D=docs/perfeng/0001-trafo-direct/artifacts
for p in d:build-cmake f:build-cmake-f l:build-cmake-l; do s=${p%%:*}; t=${p#*:}
  ctest --test-dir $t 2>&1 | tee $D/baseline-tests-$s.log
  CODSPEED_PROFILE_FOLDER=/tmp/bench-base-$s $t/benchmarks/bench_nfft_direct        # → results/*.json
  jq -s '.' /tmp/bench-base-$s/results/*.json > $D/baseline-bench-$s.json           # collate the kept artifact
done
```

Fill [`../templates/phase-a-baseline.md`](../templates/phase-a-baseline.md), which
records: the build config (exact `cmake` flags), the baseline commit SHA, the
measurement track, the full-suite test result summary **for each precision**, and the
baseline benchmark snapshot as a **Benchmark snapshot** table (canonical format with the
`prec` column — [deliverables.md](deliverables.md#canonical-formats)) covering **all** cases
in **all three precisions**.

`artifacts/` (verbatim, per precision): `baseline-tests-{d,f,l}.log` (tee'd `ctest`) and
`baseline-bench-{d,f,l}.json` (per-case stats, collated from the codspeed scratch dir). The
narrative doc embeds the *summary* table; these raw files are what Phase E diffs against.
(If a precision can't build/run the **benchmark** in your environment, capture its tests only
and **say so** — never silently drop a precision; see [precision-matrix](precision-matrix.md).)

**Exit gate** — *deliverable = exit gate*: Phase A is not exitable until the baseline is
**fully green in float, double, and long double**, `phase-a-baseline.md` + the per-precision
raw artifacts exist, and the tracker Phase A row reads `✅` with exit signal `full suite green
(d/f/l); N cases captured`. A non-green baseline in *any* precision ⇒ stop; do not proceed to Phase B.

*[← Overview & map](../REFERENCE.md) · Prev: [Preflight](preflight.md) · Next: [Phase B — correctness net](phase-b-correctness-net.md)*
