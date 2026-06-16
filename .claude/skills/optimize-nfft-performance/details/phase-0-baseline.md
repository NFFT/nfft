# Phase 0 — build tree + full baseline (the exit reference)

*[← Overview & map](../REFERENCE.md) · Prev: [Preflight](preflight.md) · Next: [Phase A — correctness net](phase-a-correctness-net.md)*

**One self-contained CMake tree drives the whole loop**, built in the mode your
[preflight](preflight.md) picked — `walltime` (local track, shown here) or `simulation`
(CodSpeed track). `FetchContent` fetches/builds codspeed-cpp (submodules and all) and
the tree produces the library, the CUnit tests, **and** the benchmark binaries together.
(Don't use the Autotools `--with-codspeed` path — it is legacy.)

```bash
cmake -S . -B build-cmake -DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON \
  -DCMAKE_C_FLAGS="-O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math"
cmake --build build-cmake -j   # libnfft3, checkall(_threads), bench_nfft_direct(_omp)
```

Both signals come from this tree:

- **correctness** — `ctest --test-dir build-cmake`, or run `build-cmake/tests/checkall`
  directly for the granular `-> OK/FAIL` stdout (same CUnit suite as `make check`).
- **performance** — `build-cmake/benchmarks/bench_nfft_direct[_omp]`, run directly
  (walltime). See [Measurement modes](measurement-modes.md) for when to *also* consult
  the CI **simulation** metric.

For network-free rebuilds, reuse a checkout with
`-DFETCHCONTENT_SOURCE_DIR_CODSPEED=<path>`. (Autotools `make check` is a valid,
CI-canonical correctness path — see
[`test-methodology.md`](../../../../docs/agents/test-methodology.md) — but it does not
build the benchmarks, so the loop stays in the CMake tree.)

Now record the complete state of **both** signals on the unmodified tree — the
contract the finished work is judged against in
[Phase D](phase-d-exit-gate.md). It must be the *full* suite, not the scoped subset
used later.

```bash
ctest --test-dir build-cmake 2>&1 | tee /tmp/baseline-tests.log              # full suite
CODSPEED_PROFILE_FOLDER=/tmp/bench-base build-cmake/benchmarks/bench_nfft_direct \
    2>/tmp/baseline-bench.log                          # ALL cases → /tmp/bench-base/results/*.json
```

(*CodSpeed track:* the benchmark baseline is the base branch's CI run, read via the
MCP in Phase D — no local capture needed; still record the test baseline.)

If the baseline is not fully green, **stop** — optimization starts only from a clean
tree. Keep these artifacts for the whole task.
