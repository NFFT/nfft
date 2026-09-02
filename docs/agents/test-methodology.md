# Test methodology (CUnit reference + accuracy tests)

How the `nfft` / `nfct` / `nfst` CUnit suites validate the transforms, and how to
extend the same pattern to a new transform. Uses the vocabulary fixed in
[`CONTEXT.md`](../../CONTEXT.md) (*ctest suite*, *precision suffix*, *orphaned test*).

## Two classes of test

1. **File-based check.** A high-precision reference is read from `tests/data/*.txt`.
   Both the **direct** (slow, `*_direct_file`) and the **fast** (`*_fast_file`) C
   transform are run on the file's input and compared against the file's reference
   output. The reference is computed offline at arbitrary precision by the
   *reference-data generator* (`tests/refgen/`, see below). This is what proves the
   **direct transform itself is correct** — it is the trusted oracle for class 2.

2. **Online check.** `setup_online` / `setup_adjoint_online` generate random input
   with the fixed `SEED`, build the reference by running the C **direct** transform,
   then compare the **fast** transform against it at larger sizes. No file is read.

The split: file-based checks pin the *direct* transform to an independent
high-precision oracle on small sizes; online checks pin the *fast* transform to the
(now-trusted) direct transform on larger sizes.

## Error model

`compare_trafo` (e.g. `tests/nfft.c:792`) returns `E∞ / ‖f̂‖₁` =
`max_j |f[j]-p.f[j]| / Σ_k |f_hat[k]|`; the adjoint variant swaps the roles.
A case passes when `err < bound`:

- Direct transforms: `bound = C · float_property(NFFT_EPSILON)` — a pure round-off
  floor. `C` = 48 (NFFT), 120 (NFCT), 130 (NFST).
- Fast transforms: `err_trafo()` returns a window- and precision-dependent analytic
  truncation bound, branched on `MANT_DIG` (113/64/53/24) and the window macro,
  floored by the direct round-off bound.

Because bounds are `multiplier · ε(MANT_DIG)` and `float_property` derives `ε` from
`MANT_DIG`, tolerances auto-scale across precisions. The `MANT_DIG == 113` (quad)
branch carries `// TODO` notes — its `a`/`b` constants are not yet tuned. The
`MANT_DIG == 64` (Intel 80-bit extended) branch is tuned and passes
(see [ADR-0003](../adr/0003-quad-precision-readiness.md)).

## Reference-data file format

Blank-line-separated sections, one scalar per line:

```
<d> / <N[0..d-1]> / <M> / <x: M*d node-major> / <f_hat: NN> / <f: M>
```

NFFT values are complex (`re im` per line); NFCT/NFST are real (one value).
NN = `∏N` (NFFT/NFCT) or `∏(N-1)` (NFST). Frequency index order is the tensor
product with dimension 0 lowest; per-dim ranges are `[⌈-N/2⌉..⌊(N-1)/2⌋]` (NFFT),
`[0..N-1]` (NFCT), `[1..N-1]` (NFST). Filenames:
`<module>[_adjoint]_<d>d_<N0>[_<N1>…]_<M>.txt`.

The two **input** sections — the nodes `x` and the input coefficients (`f_hat` for a
trafo case, `f` for an adjoint case) — are drawn as single-precision floats, so they
are exactly representable in every build precision (float through quad) and are written
**compactly**, with trailing zeros stripped. The remaining **output** section is the
high-precision reference and carries ≥64 significant digits — more than long-double
(≈19) or quad (≈34) — so the same files serve every precision; the C reader (`__FI__`)
reads each value and rounds it to the build's type on input.

## The reference-data generator

`tests/refgen/` (Python + mpmath) is the **single source of truth**. From the grids
in `tests/refgen/grids.py` it emits, in one run:

- the `tests/data/*.txt` files,
- the committed C headers `tests/data/generated/<module>_testcases.h` (the
  `testcase_delegate_file_t` declarations and `testcases_*_file[]` arrays the `.c`
  files `#include`), and
- the `tests/data/Makefile.am` `EXTRA_DIST` list.

The C build never runs Python; the artifacts are committed. The generator is a
dev-only tool run via [uv](https://docs.astral.sh/uv/); its single dependency
(`mpmath`, pinned for byte-stable output) is declared inline on the run command, so
there is no `pyproject.toml` or lockfile. Regenerate with:

```bash
uv run --with mpmath==1.3.0 python -m tests.refgen.generate --module all --precision 64
```

See [`tests/refgen/README.md`](../../tests/refgen/README.md) for all CLI options and
the generator self-tests.

The committed `tests/data/*.txt` are now produced by this generator (regenerated with
`--module all --precision 64 --seed 1`), replacing the original Mathematica-era data.

## Window evaluation reference tests

`tests/window.c` also bounds the window evaluations themselves against
reference values, separately from the transform tests above. The Kaiser-Bessel
tables were taken from mpmath by hand, with the recipe in the comment over each
one. The B-spline and sinc-power tables come from `tests/windowref`:

```bash
uv run --with mpmath==1.3.0 python -m tests.windowref.generate
```

It prints C literals for pasting into `tests/window.c`; there is no generated
header, so the values stay next to the test. Sinc powers come from mpmath at 60
digits, cardinal B-splines from exact rational arithmetic. See
`tests/windowref/README.md` for the parameter choices.

These tests call the kernels in `include/infft.h` directly rather than through
the `PHI`/`PHI_HUT` macros, so they run in every window build, not only the one
whose window they cover.

## How to add a new transform's tests

1. **Generator:** add the transform's direct forward/adjoint and frequency index
   set to `tests/refgen/transforms.py`; add its grids to `tests/refgen/grids.py`;
   teach `tests/refgen/registration.py` the module name. Add pytest self-tests.
2. **Test file:** clone `tests/nfct.c` (real) or `tests/nfft.c` (complex). Adjust
   the trailing `#include "<module>.h"`, the `NN` formula, the `err_trafo_direct`
   constant, the window×precision `a`/`b` bound tables, and the registered grids.
   Replace the inline file-testcase block with
   `#include "data/generated/<module>_testcases.h"`.
3. **Registration:** in `tests/check.c` add `#include`, an `#ifdef HAVE_<MODULE>`
   guard, `#define X(name) <MODULE>(name)`, a `CU_add_suite`, and `CU_add_test`
   lines mirroring the 1d/2d/3d (+ exhaustive online) pattern.
4. **Build wiring:** add the sources to `tests/Makefile.am` (`checkall_SOURCES`) and
   the generated header to `EXTRA_DIST`.
5. Run `uv run --with mpmath==1.3.0 python -m tests.refgen.generate --module <module>`
   then `make check`.

The harness requires the module to expose `X(init*)`, `X(precompute_one_psi)`,
`X(trafo)`, `X(trafo_direct)`, `X(adjoint)`, `X(adjoint_direct)`, `X(check)`,
`X(finalize)`, and a plan struct with `.x/.f/.f_hat/.m/.sigma/.d/.N/.M_total/.flags`.
Transforms without an exact reference (e.g. **solver**, iterative) need a different
assertion (convergence / round-trip recovery), not file-based reference data.
