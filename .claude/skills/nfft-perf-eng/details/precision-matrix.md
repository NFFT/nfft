# Precision matrix (float · double · long double)

*[← Overview & map](../REFERENCE.md) — cross-cutting reference; applies to Phases A, E, F.*

The library is **precision-agnostic**: the same C sources compile in float, double, and
long-double precision (`Y()`/`X()`/`R`/`C` mangling — see `CONVENTIONS.md`). A change can
therefore pass in **double** yet regress in another precision — float has ~7 significant
digits (overflow/underflow and cancellation bite sooner), long double exercises a different
codepath/ABI and `COSL/SINL`. A double-only loop **misses single-precision-only and
long-double-only regressions** — both correctness (a fault that stays within the double
tolerance but breaks float) and performance. So **correctness *and* performance are measured
in all three precisions** at the baseline (A), the inner loop (D), and the exit gate (E).

## Three build trees — one precision each (the FFTW model)

Precision is fixed at configure time; `NFFT_ENABLE_FLOAT` and `NFFT_ENABLE_LONG_DOUBLE` are
mutually exclusive, default is double. Use a **separate tree per precision**, suffixed like
the library (`""`/`f`/`l`):

```bash
# CLEAR a stale in-source config.h FIRST. A prior Autotools `./configure` leaves a
# double-precision include/config.h in the tree; the CMake compile picks it up over the
# per-tree generated float/long-double config, so the library compiles double-mangled
# (nfft_*) while the float-compiled tests want nfftf_* → "undefined reference to nfftf_init_1d".
# CI does exactly this (build-linux.yml). Skip it and the float/long-double trees mis-link.
[ -f include/config.h ] && mv include/config.h include/config.h.autotools.bak

# bash array — keeps the quoted CMAKE_C_FLAGS a SINGLE argument (unquoted $FLAGS would
# split it on spaces and break configure with stray -g/-f… args).
# -falign-functions=64 -falign-loops=32 match CI (bench-linux.yml) — deterministic code
# layout, so an untouched neighbour can't show a phantom alignment regression (see caveats.md).
FLAGS=(-DNFFT_BENCHMARK_MODE=walltime -DNFFT_ENABLE_OPENMP=ON \
       -DCMAKE_C_FLAGS="-O3 -g -fomit-frame-pointer -fstrict-aliasing -ffast-math -falign-functions=64 -falign-loops=32")
cmake -S . -B build-cmake   "${FLAGS[@]}"                              # double
cmake -S . -B build-cmake-f "${FLAGS[@]}" -DNFFT_ENABLE_FLOAT=ON       # float  (suffix f)
cmake -S . -B build-cmake-l "${FLAGS[@]}" -DNFFT_ENABLE_LONG_DOUBLE=ON # long double (suffix l)
for t in build-cmake build-cmake-f build-cmake-l; do cmake --build $t -j; done
```

> **Symptom of a stale `config.h`:** the non-double library exports unsuffixed `nfft_*`
> names (e.g. `nfft_trafo_direct` in `libnfft3f.so`) and is missing `nfftf_*`, so tests,
> examples, and the benchmark fail to link with `undefined reference to nfftf_init_1d`. If
> you see that, you skipped the clear-config step — it is **not** a repo build defect.

Each tree builds its own `checkall`(`_threads`) **and** `bench_nfft_direct`(`_omp`): the
CUnit suite always compiles per precision (with precision-appropriate tolerances — the very
thing that catches a precision-specific break), and the benchmark is precision-agnostic too
(`BENCHMARK_AGNOSTIC_PRECISION` defaults ON), so it measures in each tree. Reuse the fetched
codspeed checkout across trees with `-DFETCHCONTENT_SOURCE_DIR_CODSPEED=<path>` to avoid
re-downloading.

> If a precision's tree cannot build or run the **benchmark** in your environment (e.g. the
> float/long-double FFTW variant is missing, or the benchmark is disabled), fall back to
> **correctness-only** for that precision, and **say so** in the deliverable — never silently
> drop a precision.

## What each phase runs across the matrix

- **Phase A** — build all three trees; capture the full test + benchmark baseline in each.
- **Phase E** — after every kept change, run the **B-net** in all three (`checkall` per tree,
  `checkall_threads` if OpenMP touched) and re-measure the **C-metric** in all three. The net
  is the scoped subset, so 3× is cheap; this is what catches a float-only break *early*.
- **Phase F** — re-run the ENTIRE baseline (full `ctest` + ALL benchmark cases) in all three;
  a regression in *any* precision fails the gate.

## Artifacts & deliverable formats

Suffix per-precision captures `-d` / `-f` / `-l`: `baseline-tests-d.log`,
`baseline-bench-f.json`, `final-bench-l.json`, … . The **Benchmark snapshot** and
**Comparison table** carry a `prec` column (`d`/`f`/`l`) so all three precisions sit in one
comparable table (see [deliverables.md](deliverables.md#canonical-formats)). The cost is real
— 3× build/test/benchmark time — but it is the price of catching precision-specific
regressions, and the inner loop pays it only on the narrow net/metric.
