# Caveats

*[← Overview & map](../REFERENCE.md) — cross-cutting reference, consult from any phase.*

- **Crash-class faults lose granularity.** A Phase-B fault that corrupts buffer sizes
  or indices (rather than just values) can make a case under-allocate and `SIGABRT`
  (exit 134) partway through the suite. The `-> FAIL` lines printed *before* the abort
  are still a reliable net, but the XML may be truncated. Prefer a fault that yields
  wrong-but-finite results — e.g. the imaginary-sign flip in `trafo_direct`'s kernel
  (Phase B) — over one that aborts when you want a complete failure list.
- **OpenMP-only changes show only in `checkall_threads`.** That binary links the
  `_omp` library; `checkall` does not. A change guarded by `#ifdef _OPENMP` (or
  affecting only parallel scheduling) will leave `checkall` green and must be judged
  by `checkall_threads`.
- **A wrong target is possible — but Phase C catches it, so you needn't pre-judge.**
  Some code is simply not on any timed path: `init`-only helpers like `intprod()` run
  only in `*_init*` (plan setup), and the benchmarks call `init_*` **outside** the
  timed `for (auto _ : state)` loop (`benchmarks/bench_nfft_direct.cpp:65`) — they time
  only `trafo_direct` / `adjoint_direct`. You don't have to know that up front: the
  Phase C hard gate makes it self-evident — the deliberate slowdown moves *no*
  benchmark, so there is no metric, so the loop **stops**. The remedy is to target code
  inside the timed region (`trafo*` / `adjoint*`) or add a benchmark that covers it.
- **Benchmark coverage is narrow.** Only the **direct** (slow, O(N·M)) transforms
  are benchmarked, forward and adjoint, in 1d/2d/3d (`bench_nfft_direct.cpp`). The
  **fast** NFFT path (`trafo`, `adjoint`, the `precompute_one_psi` strategies) and
  every non-NFFT module have no benchmark yet. Optimizing those requires adding a
  benchmark first (Phase C with no coverage is meaningless).
- **`[ERROR] instrument-hooks: failed to write environment.json` is benign.** A
  walltime run prints it when the **current working directory** (where the binary
  tries to drop `environment.json`) is not writable; the per-case results JSON is
  still written to `$CODSPEED_PROFILE_FOLDER/results/` and the exit code is 0. Run
  from a writable dir to silence it, or ignore it.
- **Walltime regressions on untouched code are usually noise — but not always random
  noise.** A change confined to one function (e.g. `trafo_direct`'s 1d branch) will still
  show neighbouring, *byte-identical* cases (`adjoint_direct`, other dims) swinging ±several
  percent run to run — worst on the few-iteration 2d/3d cases. This is why Phase F uses a
  noise threshold and a re-run, not a bare "nothing got slower" check.
- **A *reproducible* regression on untouched code can still be a code-layout artifact, not a
  real slowdown.** Editing the target changes its size, which shifts the alignment of whatever
  the linker places after it; an untouched neighbour's hot loop can land on a worse boundary
  and read a *stable* +N% that survives re-runs (so the noise-rule re-run does **not** clear
  it). The fix is to pin layout: build with `-falign-functions=64 -falign-loops=32` (CI's
  `bench-linux.yml` does — the Phase-A flag list now includes them). If a control case still
  regresses reproducibly, confirm attribution by re-measuring the **baseline** under the same
  aligned flags (`git stash` the change, rebuild, measure): a layout artifact vanishes; a true
  coupling persists. Observed concretely: a float `adjoint_direct_1d[32/100]` control read +7.5%
  without the flags and −0.9% (flat) with them.
- **A green suite is narrow, not a guarantee.** The Phase-B net is a deliberate subset,
  and even the full Phase-A suite only checks specific transforms at specific sizes and
  grids. An optimization can be faster, pass every test, and still have dropped accuracy
  on a path no case exercises — most often at transform sizes larger than any tested
  case, or for node/coefficient distributions no case spans. Passing tests is
  authoritative *for what the tests cover*; the gap is what the **risk assessment**
  ([risk-assessment.md](risk-assessment.md)) surveys and [extending-tests.md](extending-tests.md)
  can close.
