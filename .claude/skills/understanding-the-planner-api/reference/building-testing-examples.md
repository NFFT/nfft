# Building, testing & running the planner API

See the repo `CLAUDE.md`/`AGENTS.md` for the general build. This page is the
**planner-API-specific** delta. Ground truth: `configure.ac`, `kernel/*/Makefile.am`,
`kernel/CMakeLists.txt`, `tests/Makefile.am`, `examples/nfft/Makefile.am`.

## The library always contains the planner

`kernel/planner/` builds **unconditionally** (like `kernel/util`) in both build
systems, all three precisions, serial + OpenMP — there is no configure gate for
the planner core or the NFFT `plan_ng` surface. So the API is in
`libnfft3<suffix>` whenever you build the NFFT module.

## Configure / build

```bash
./bootstrap.sh                                   # after editing any Makefile.am / configure.ac
./configure --enable-all --enable-openmp --enable-tests
make -j
```

- **`--enable-tests` is required to build the CUnit binaries** (including
  `checkall_ng`) and is **not** implied by `--enable-all`. Without it, `make
  check` runs nothing for the planner suite.
- `--enable-openmp` additionally produces `libnfft3<suffix>_omp` and the
  threaded test binaries (`checkall_ng_threads`).
- Precision: default double; `--enable-float` / `--enable-long-double` (separate
  configured trees). Long-double `make check` is very slow on some hosts — prefer
  build-only for long-double, and run tests in double + float.
- Window: `--with-window=kaiserbessel|gaussian|bspline|sinc|dirac` sets the
  **compile-time** window (`Y(get_window_id)()`). The native fast path can also
  take *any* implemented window at runtime via the guru's `window` argument.

## `--enable-debug` (the planner's extra checks)

`--enable-debug` does two things that matter here:

1. Compiles the `A(...)` runtime assertions — including the **execute-asserts-AWAKE**
   check and the **md5 byte-identity guard** (`kernel/nfft/xcheck.c`,
   `Y(nfft_x_md5)`/`Y(nfft_x_verify)`) that verifies the internal `x` copy is
   restored byte-for-byte between racing candidates. The guard is sized from the 
   **compressed** rank (`sz->rnk * M`, not `d * M`) — it must be, or it over-reads.
2. Switches `CFLAGS` to an Address/UBSan build
   (`-fsanitize=address,undefined -fno-sanitize-recover=all`). Debug binaries are
   slower and need the sanitizer runtimes.

In **release** these asserts are no-ops. Run the debug build to actually catch
lifecycle violations (missing `precompute`, double-destroy, x-restore breakage).

## Tests

The new-planner suite is a **separate** CUnit binary, `checkall_ng` (plus
`checkall_ng_threads` under OpenMP), distinct from the legacy `checkall`.

```bash
make check                # builds + runs both suites; PASS/FAIL summary
tests/checkall_ng         # run the planner suite directly, full output
```

`checkall_ng` sources (`tests/Makefile.am`): `check_ng.c` + shared `util.c` +
`planner.c` (wisdom-store unit tests) + `nplan.c` / `nplan_data.c` /
`nplan_perm.c` (the bundle lifecycle, file-based reference checks, and the
permuting-solver x-restore guard tests) + `nfast.c` (native-fast vs legacy/direct
checks). Reference data: `tests/data/nfft_{,adjoint_}{1d,2d}_*.txt`, listed in
`tests/data/generated/nfft_native_testcases.h`.

> `make check` can no-op on an already-built tree; run `tests/checkall_ng`
> directly to be sure it actually executed (verify `exit 0`).

Notes: the measured-race timing tests are known-skips (visible `CU_PASS`) on
hosts without a cycle counter, and hard-fail only under
`--enable-exhaustive-unit-tests`. Negative-argument tests
(`check_nplan_guru_rejects_null_args`, `..._bad_geometry`) cover the release-safe
guard; the x-restore guard tests are no-ops unless `--enable-debug`.

## CMake

CMake builds **the same** new-API sources:
`kernel/CMakeLists.txt` lists `nfft/nfft-nd.c`, `nfft/nconst.c`, `nfft/xcheck.c`,
the whole `conv/` and `deconv/` modules, plus `checkall_ng` / `checkall_ng_threads`
`ctest` targets. `nfct`/`nfst` are legacy-only (their planner glue was removed,
not merely un-wired). If you add a planner source file, add it to **both**
`Makefile.am` and `kernel/CMakeLists.txt`.

## Runnable examples (`examples/nfft/`)

- **`nfast_native.c`** — the best worked example. A six-way check over one 1D
  problem (`data/nfft_1d_8192_128.txt`): legacy direct NDFT, legacy fast NFFT,
  planner direct NDFT, and planner native fast (each fast one twice, FFTW
  estimate vs measure). It forces specific solvers via flags (e.g.
  `NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE` to isolate the direct NDFT), prints each plan tree with
  `nfft_fprint_plan`, times precompute/forward/adjoint separately, and checks the
  forward error against the file reference. Shows correct guru → `precompute` →
  `execute`/`execute_adjoint_on` → `plan_ng_destroy` usage.
- **`ndft_fast.c`** — direct-NDFT-focused example.

Build them with `--enable-all`; both are `noinst_PROGRAMS` linking
`libnfft3<suffix>`.

## Reading a plan tree

`nfft_fprint_plan(p, stdout)` prints an S-expression naming the chosen solver, e.g.:

```
(nfft-plan-ng (fwd (nfft_solver_fast_native pcost=... (deconv (deconv_solver_nd pcost=...)) (conv (conv_solver_nd pcost=...)))) (adj (null)))
```

The `(adj (null))` is expected (forward-only race — the forward plan serves the
adjoint). A direct winner prints `(<solver-name> pcost=...)` — e.g.
`nfft_solver_ndft_1d`, `nfft_solver_ndft_nd`,
`nfft_solver_const_0d`.
