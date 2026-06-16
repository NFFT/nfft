# Caveats

*[← Overview & map](../REFERENCE.md) — cross-cutting reference, consult from any phase.*

- **Crash-class faults lose granularity.** The Phase-A `intprod` fault corrupts the
  allocation sizes `N_total`/`n_total` (`kernel/nfft/nfft.c:5962`), so multi-dim
  cases under-allocate and the binary `SIGABRT`s (exit 134) partway through. The
  `-> FAIL` lines printed *before* the abort are still a reliable net, but the XML
  may be truncated. Prefer a fault that yields wrong-but-finite results over one that
  aborts when you want a complete failure list.
- **OpenMP-only changes show only in `checkall_threads`.** That binary links the
  `_omp` library; `checkall` does not. A change guarded by `#ifdef _OPENMP` (or
  affecting only parallel scheduling) will leave `checkall` green and must be judged
  by `checkall_threads`.
- **`intprod` is a deliberately bad example.** It runs only in `*_init*` (plan
  setup), and the current benchmarks call `init_*` **outside** the timed
  `for (auto _ : state)` loop (`benchmarks/bench_nfft_direct.cpp:65`) — they time
  only `trafo_direct` / `adjoint_direct`. So Phase B will show *no* benchmark moving
  for an `intprod` change: it is not a meaningful performance target under the
  current suite. Real targets live inside `trafo*` / `adjoint*`.
- **Benchmark coverage is narrow.** Only the **direct** (slow, O(N·M)) transforms
  are benchmarked, forward and adjoint, in 1d/2d/3d (`bench_nfft_direct.cpp`). The
  **fast** NFFT path (`trafo`, `adjoint`, the `precompute_one_psi` strategies) and
  every non-NFFT module have no benchmark yet. Optimizing those requires adding a
  benchmark first (Phase B with no coverage is meaningless).
- **`[ERROR] instrument-hooks: failed to write environment.json` is benign.** A
  walltime run prints it when the **current working directory** (where the binary
  tries to drop `environment.json`) is not writable; the per-case results JSON is
  still written to `$CODSPEED_PROFILE_FOLDER/results/` and the exit code is 0. Run
  from a writable dir to silence it, or ignore it.
- **Walltime regressions on untouched code are noise, not bugs.** A change confined to
  one function (e.g. `trafo_direct`'s 1d branch) will still show neighbouring,
  *byte-identical* cases (`adjoint_direct`, other dims) swinging ±several percent run
  to run — worst on the few-iteration 2d/3d cases. This is why Phase D uses a noise
  threshold and a re-run, not a bare "nothing got slower" check.
