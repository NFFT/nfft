# Reference-Data Generator & Quad-Precision Readiness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the Mathematica reference-data notebooks with an agent-modifiable, arbitrary-precision Python (mpmath) generator that is the single source of truth for the NFFT/NFCT/NFST CUnit reference data, document the whole test methodology so it can be applied to new transforms, and make the data pipeline ready for 128-bit quadruple precision.

**Architecture:** A self-contained Python package `tests/refgen/` computes the *direct* (slow) NDFT/NDCT/NDST at arbitrary precision with `mpmath`, and emits three coupled artifacts from one set of grid definitions: (1) the `tests/data/*.txt` reference files in the exact byte format the existing `setup_file()` parsers consume, (2) committed C headers `tests/data/generated/<module>_testcases.h` holding the `testcase_delegate_file_t` declarations + `testcases_*_file[]` arrays (eliminating the manual-paste drift that left `check_nfsft.nb` stale), and (3) the `tests/data/Makefile.am` `EXTRA_DIST` list. The C build never invokes Python — the generated artifacts are committed; the generator is a dev/CI tool, exactly as the notebooks were. Precision is a CLI parameter (default 64 decimal digits, ≥ quad's ~34), so the same tool serves double, long-double, and future quad/octuple reference data. The legacy notebooks are removed; full `__float128` *build* support is scoped in an ADR as follow-on work, not implemented here.

**Tech Stack:** Python 3 + `mpmath` (dev dependency, managed with uv via `pyproject.toml`/`uv.lock`; not a build or runtime dependency), `pytest` for generator self-tests; C99 + CUnit + Autotools for the consumed tests; Markdown for the methodology doc and ADRs.

---

## Background facts this plan relies on (verified against the tree)

- **Two test classes** in `tests/nfft.c` / `nfct.c` / `nfst.c`:
  - **File-based check** — loads a `tests/data/*.txt` reference and compares both the *direct* (`*_direct_file`) and *fast* (`*_fast_file`) C transforms against it. The `.txt` holds high-precision reference values.
  - **Online check** — `setup_online`/`setup_adjoint_online` generate random input with a fixed `SEED`, build the reference by running the C *direct* transform, then compare the *fast* transform against it. No file involved.
  This plan only changes the **file-based** reference data and its generation; the online path is untouched.
- **Exact `.txt` format** (parsed by `setup_file`, e.g. `tests/nfft.c:426`): sections separated by a blank line, each scalar on its own line:
  ```
  <d>
  <blank>
  <N[0]> ... <N[d-1]>            # one integer per line
  <blank>
  <M>
  <blank>
  <x>                           # M*d lines, node-major: node j's d coords contiguous
  <blank>
  <f_hat>                       # NN lines (NFFT: "re im"; NFCT/NFST: one real)
  <blank>
  <f>                           # M lines (NFFT: "re im"; NFCT/NFST: one real)
  <blank>
  ```
  The same section order (x, then f_hat, then f) is used for **both** forward and adjoint files. In an `*_adjoint_*` file, `f` is the random input and `f_hat` is the computed reference.
- **NN (coefficient count) per module:** NFFT `∏ N[i]`; NFCT `∏ N[i]`; NFST `∏ (N[i]-1)`.
- **Frequency index sets (must match the C `f_hat` ordering byte-for-byte):**
  - NFFT: per dim `k ∈ [⌈-N/2⌉ … ⌊(N-1)/2⌋]`.
  - NFCT: per dim `k ∈ [0 … N-1]`.
  - NFST: per dim `k ∈ [1 … N-1]`.
  - Multi-dim ordering is the tensor product with **dimension 0 varying slowest** (row-major), i.e. `itertools.product(freqs[0], freqs[1], …)`.
- **Transform definitions:**
  - NFFT forward `f[j] = Σ_k f_hat[k]·exp(-2πi·k·x[j])`; adjoint `f_hat[k] = Σ_j f[j]·exp(+2πi·k·x[j])`.
  - NFCT forward `f[j] = Σ_k f_hat[k]·∏_i cos(2π k_i x_{j,i})`; adjoint swaps the summation index. Real.
  - NFST as NFCT with `sin`. Real.
- **Random input ranges:** NFFT nodes `x∈[-1/2,1/2)`, coeffs complex with re,im `∈[-1,1]`; NFCT/NFST nodes `x∈[0,1/2]`, coeffs real `∈[-1,1]`.
- **Reproducibility:** Mathematica's `SeedRandom[1]` cannot be byte-reproduced outside Mathematica, so the regenerated files carry *new* (but equally valid) random inputs. The acceptance test is therefore "regenerate, build, `make check` passes within the existing error bounds", **not** byte-equality with the old files. The error bounds (`tests/nfft.c:213` `48·eps`, `nfct.c:215` `120·eps`, `nfst.c:215` `130·eps`) already cover the direct-vs-reference round-off and do not depend on the specific random draw.
- **Filename scheme:** `<module>[_adjoint]_<d>d_<N0>[_<N1>…]_<M>.txt` (e.g. `nfft_1d_10_1.txt`, `nfct_2d_10_25_50.txt`, `nfft_adjoint_3d_10_10_10_10.txt`).
- **Grids currently shipped** (the generator must reproduce these exact `(d, N, M)` sets so existing filenames are preserved):
  - NFFT 1d: N∈{1,2,4,10,20,50} × M∈{1,10,20,50}; 2d: N∈{(10,10),(10,20),(20,10),(20,20)} × M∈{20,50}; 3d: (10,10,10) × M∈{10}.
  - NFCT 1d: N∈{1,2,4,10,25,50} × M∈{1,10,25,50}; 2d: N∈{(10,10),(10,25),(25,10),(25,25)} × M∈{25,50}; 3d: (10,10,10) × M∈{10}.
  - NFST 1d: N∈{2,4,10,25,50} × M∈{1,10,25,50}; 2d: same N as NFCT × M∈{25,50}; 3d: (10,10,10) × M∈{10}.
    (NFST starts at N=2 because it has N−1 frequencies; N=1 → zero coefficients.)
  - Forward **and** adjoint files for every grid entry.
- **Mangling discipline for the methodology doc:** `Y(name)` library-wide, `X(name)` module-local, `FFTW(name)` for FFTW; `R`/`C` types; `K(x)` literals; `__FI__`/`__FE__` scan/print formats (`include/infft.h`).
- **Domain vocabulary** already fixed in `CONTEXT.md`: *precision suffix*, *ctest suite* (= util/nfft/nfct/nfst), *orphaned test* (`tests/check_nfsft.c`). New terms introduced here (e.g. *reference-data generator*, *file-based check*, *online check*) are candidates for `/grill-with-docs` to add to the glossary — do not unilaterally rename existing terms.

---

## File Structure

**New — the generator (dev/CI tool, not built by `make`):**
- `tests/refgen/__init__.py` — package marker.
- `tests/refgen/transforms.py` — arbitrary-precision NDFT/NDCT/NDST forward+adjoint and frequency index sets. One responsibility: the math.
- `tests/refgen/io_format.py` — number formatting and `.txt` writing/parsing in the exact harness format. One responsibility: serialization.
- `tests/refgen/grids.py` — the `(d, N, M)` test-case grids per module. Single source of truth for *which* cases exist.
- `tests/refgen/registration.py` — emits the C `<module>_testcases.h` header and the `EXTRA_DIST` list. One responsibility: codegen.
- `tests/refgen/generate.py` — CLI entry point tying the above together (`--module`, `--precision`, `--seed`, `--data-dir`, `--header-dir`).
- `tests/refgen/tests/test_transforms.py` — pytest self-tests for the math.
- `tests/refgen/tests/test_io_format.py` — pytest round-trip tests for the format.
- `tests/refgen/README.md` — how to install deps (uv) and regenerate.

**New — generated, committed C artifacts (consumed by the build):**
- `tests/data/generated/nfft_testcases.h`
- `tests/data/generated/nfct_testcases.h`
- `tests/data/generated/nfst_testcases.h`

**New — documentation:**
- `docs/agents/test-methodology.md` — the documented methodology + "how to add a transform" recipe.
- `docs/adr/0002-python-reference-data-generator.md` — decision to replace Mathematica.
- `docs/adr/0003-quad-precision-readiness.md` — scoped `__float128` build gap.

**Modified:**
- `tests/nfft.c`, `tests/nfct.c`, `tests/nfst.c` — replace the inline `testcase_delegate_file_t` declarations + `testcases_*_file[]` arrays with `#include "data/generated/<module>_testcases.h"`.
- `tests/data/Makefile.am` — `EXTRA_DIST` regenerated to match.
- `tests/Makefile.am` — add the generated headers to `EXTRA_DIST` / `noinst_HEADERS` so `make dist` ships them.
- `CONTEXT.md` — new glossary terms (left for `/grill-with-docs`, stubbed here).

**Removed (legacy transform-data notebooks, superseded):**
- `tests/check_nfft.nb`, `tests/check_nfft.m`
- `tests/check_nfct.nb`, `tests/check_nfct.m`
- `tests/check_nfst.nb`, `tests/check_nfst.m`
- `tests/check_nfsft.nb`, `tests/check_nfsft.m` (already-stale NFST clone)

**Kept (out of scope for this plan — deliberately not removed):**
- `tests/PrintVector.m` — still used by `check_bspline.m` (`<<PrintVector``); removing it would break that generator.
- `tests/check_bspline.{nb,m}` — *also* a test-data generator (it samples the B-spline window and emits the `b0…b30` arrays in `tests/bspline.c`), but a **different artifact** (a sampled-window C array, not a transform `.txt`). Deferred to a follow-up that exercises the methodology doc's "add a new generator target" recipe.
- `kernel/util/bessel_i0.nb`, `kernel/util/bspline.nb` — genuine **library codegen** (minimax `I₀` coefficients; the library's B-spline evaluator), not test data; out of scope.

---

## Phase 0 — Document the current methodology

This phase is pure documentation. It de-risks the rest by forcing a written, reviewable model of the harness, and it is itself a deliverable ("document the approach so it can be applied to other transforms").

### Task 0.1: Write the test-methodology document

**Files:**
- Create: `docs/agents/test-methodology.md`

- [ ] **Step 1: Write the document**

Create `docs/agents/test-methodology.md` with exactly this content:

````markdown
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
and `64` (extended) branches carry `// TODO` notes — their `a`/`b` constants are not
yet tuned (see [ADR-0003](../adr/0003-quad-precision-readiness.md)).

## Reference-data file format

Blank-line-separated sections, one scalar per line:

```
<d> / <N[0..d-1]> / <M> / <x: M*d node-major> / <f_hat: NN> / <f: M>
```

NFFT values are complex (`re im` per line); NFCT/NFST are real (one value).
NN = `∏N` (NFFT/NFCT) or `∏(N-1)` (NFST). Frequency index order is the tensor
product with dimension 0 slowest; per-dim ranges are `[⌈-N/2⌉..⌊(N-1)/2⌋]` (NFFT),
`[0..N-1]` (NFCT), `[1..N-1]` (NFST). Filenames:
`<module>[_adjoint]_<d>d_<N0>[_<N1>…]_<M>.txt`.

The values carry ≥64 significant digits — more than long-double (≈19) or quad (≈34),
so the same files serve every precision; the C reader (`__FI__`) rounds on input.

## The reference-data generator

`tests/refgen/` (Python + mpmath) is the **single source of truth**. From the grids
in `tests/refgen/grids.py` it emits, in one run:

- the `tests/data/*.txt` files,
- the committed C headers `tests/data/generated/<module>_testcases.h` (the
  `testcase_delegate_file_t` declarations and `testcases_*_file[]` arrays the `.c`
  files `#include`), and
- the `tests/data/Makefile.am` `EXTRA_DIST` list.

The C build never runs Python; the artifacts are committed. The generator is a
uv-managed dev tool (deps pinned in `pyproject.toml` / `uv.lock`). Regenerate with:

```bash
uv sync                                  # one-time: create .venv from the lockfile
uv run python -m tests.refgen.generate --module all --precision 64
```

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
5. Run `python -m tests.refgen.generate --module <module>` then `make check`.

The harness requires the module to expose `X(init*)`, `X(precompute_one_psi)`,
`X(trafo)`, `X(trafo_direct)`, `X(adjoint)`, `X(adjoint_direct)`, `X(check)`,
`X(finalize)`, and a plan struct with `.x/.f/.f_hat/.m/.sigma/.d/.N/.M_total/.flags`.
Transforms without an exact reference (e.g. **solver**, iterative) need a different
assertion (convergence / round-trip recovery), not file-based reference data.
````

- [ ] **Step 2: Verify the doc renders and links resolve**

Run: `test -f docs/agents/test-methodology.md && grep -c '##' docs/agents/test-methodology.md`
Expected: prints a count ≥ 5 (the five `## ` section headers exist; the Step 1 content defines exactly five H2 sections).

- [ ] **Step 3: Commit**

```bash
git add docs/agents/test-methodology.md
git commit -m "docs: document CUnit reference/accuracy test methodology"
```

---

## Phase 1 — Build the Python reference-data generator (with self-tests)

TDD here means: the generator's *math* is verified by independent pytest checks
(equispaced NDFT equals the DFT; adjoint equals conjugate-transpose; hand-computed
small cases; format round-trip) **before** it is wired into the C suite. The C
`make check` in Phase 2 is the integration test.

### Task 1.1: Package skeleton and dependency note

**Files:**
- Create: `tests/refgen/__init__.py`
- Create: `tests/refgen/tests/__init__.py`
- Create: `pyproject.toml` (repo root)
- Create: `uv.lock` (repo root, generated by `uv lock`)
- Create: `tests/refgen/README.md`

- [ ] **Step 0: Pin the dependency (byte-reproducibility contract)**

The generator is managed with [uv](https://docs.astral.sh/uv/) via a repo-root
`pyproject.toml` (a `[tool.uv] package = false` virtual project — uv manages the
venv + locked deps but never tries to build the C repo as a Python package).

`pyproject.toml` (repo root):
```toml
# Dev tooling only — NOT part of the C library build (autotools/CMake ignore this).
# Declares the uv-managed environment for the reference-data generator in
# tests/refgen/ (see docs/agents/test-methodology.md and docs/adr/0002-*.md).
# `package = false` makes this a virtual project: uv manages the venv + locked
# deps but never tries to build the C repo as a Python package.
[project]
name = "nfft-refgen"
version = "0"
description = "Arbitrary-precision reference-data generator for the NFFT3 CUnit test suite"
requires-python = ">=3.9"
dependencies = [
    "mpmath==1.3.0", # pinned: the drift guard regenerates and diffs, so output must be byte-stable
]

[dependency-groups]
dev = [
    "pytest>=7",
]

[tool.uv]
package = false
```

Then produce the lockfile:
```bash
cd /workspaces/nfft && uv lock && uv sync   # writes uv.lock and creates .venv
```

> The CI drift guard (Phase 4) regenerates and `git diff --exit-code`s, so the
> *formatting* must be byte-stable. `mpmath` is pinned (in `pyproject.toml`) and the
> exact resolution is captured in `uv.lock`; anyone regenerating locally runs
> `uv sync --locked` to get that same version, or the guard will flag spurious diffs.
> (`random`/Mersenne-Twister and the `nstr` algorithm are otherwise stable across all
> CPython 3.x.) If `mpmath==1.3.0` is unavailable at execution time, pin whatever
> version `uv lock` resolves to and update `pyproject.toml` + `uv.lock` + the README.

- [ ] **Step 1: Create the package markers**

`tests/refgen/__init__.py`:
```python
"""Arbitrary-precision reference-data generator for the NFFT3 CUnit suite.

Replaces the legacy Mathematica notebooks (check_nfft/nfct/nfst.nb). See
docs/agents/test-methodology.md.
"""
```

`tests/refgen/tests/__init__.py`:
```python
```
(empty file)

- [ ] **Step 2: Write the README**

`tests/refgen/README.md`:
````markdown
# Reference-data generator

Computes the direct NDFT/NDCT/NDST at arbitrary precision (mpmath) and emits the
CUnit reference data + generated C registration headers + the data `Makefile.am`
list. Single source of truth — see `docs/agents/test-methodology.md`.

## Install (dev only; not a build dependency)

Managed by [uv](https://docs.astral.sh/uv/). From the repository root:

```bash
uv sync   # creates .venv from the locked deps (mpmath, pytest); byte-stable output
```

## Regenerate everything

```bash
# from the repository root
uv run python -m tests.refgen.generate --module all --precision 64
```

## Run the generator's own tests

```bash
uv run python -m pytest tests/refgen/tests -q
```

## Options

- `--module {nfft,nfct,nfst,all}`
- `--precision N`   working/printed significant decimal digits (default 64; ≥34 for quad)
- `--seed N`        PRNG seed (default 1)
- `--data-dir P`    output dir for *.txt (default tests/data)
- `--header-dir P`  output dir for generated headers (default tests/data/generated)
````

- [ ] **Step 3: Verify the package imports**

Run: `cd /workspaces/nfft && uv run python -c "import tests.refgen; print('ok')"`
Expected: `ok`

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml uv.lock tests/refgen/__init__.py \
        tests/refgen/tests/__init__.py tests/refgen/README.md
git commit -m "test(refgen): add generator package skeleton"
```

### Task 1.2: Frequency index sets + transforms (TDD)

**Files:**
- Create: `tests/refgen/transforms.py`
- Test: `tests/refgen/tests/test_transforms.py`

- [ ] **Step 1: Write the failing test**

`tests/refgen/tests/test_transforms.py`:
```python
import itertools
import cmath

import mpmath
import pytest

from tests.refgen import transforms as T


def test_freqs_match_c_ordering():
    # NFFT even/odd, NFCT, NFST per-dim ranges
    assert T.freqs("nfft", [10]) == [list(range(-5, 5))]
    assert T.freqs("nfft", [5]) == [list(range(-2, 3))]
    assert T.freqs("nfft", [1]) == [[0]]
    assert T.freqs("nfct", [4]) == [[0, 1, 2, 3]]
    assert T.freqs("nfst", [4]) == [[1, 2, 3]]


def test_nn_counts():
    assert T.nn("nfft", [10, 20]) == 200
    assert T.nn("nfct", [10, 25]) == 250
    assert T.nn("nfst", [10, 25]) == 9 * 24


def test_ndft_equispaced_matches_dft():
    # On equispaced nodes x_j = j/N - 1/2 the NDFT reduces to a DFT we can
    # cross-check with cmath in double precision.
    mpmath.mp.dps = 30
    N, M = 8, 8
    x = [(mpmath.mpf(j) / N - mpmath.mpf(1) / 2,) for j in range(M)]
    f_hat = [mpmath.mpc(j + 1, -j) for j in range(N)]
    f = T.trafo("nfft", [N], M, x, f_hat)
    ks = list(range(-(N // 2), N - (N // 2)))
    for j in range(M):
        xj = float(x[j][0])
        ref = sum((j_c.real + 1j * j_c.imag) * cmath.exp(-2j * cmath.pi * k * xj)
                  for k, j_c in zip(ks, [complex(float(c.real), float(c.imag)) for c in f_hat]))
        got = complex(float(f[j].real), float(f[j].imag))
        assert abs(got - ref) < 1e-10


def test_adjoint_is_conjugate_transpose_nfft():
    mpmath.mp.dps = 30
    N, M = 5, 4
    x = [(mpmath.mpf(2 * j + 1) / (4 * M) - mpmath.mpf(1) / 2,) for j in range(M)]
    f = [mpmath.mpc(j + 1, j) for j in range(M)]
    f_hat = T.adjoint("nfft", [N], M, x, f)
    ks = list(range(-(N // 2), N - (N // 2)))
    for idx, k in enumerate(ks):
        ref = sum(f[j] * mpmath.expj(2 * mpmath.pi * k * x[j][0]) for j in range(M))
        assert abs(f_hat[idx] - ref) < mpmath.mpf(10) ** -25


def test_ndct_single_term_handcomputed():
    mpmath.mp.dps = 40
    # one coefficient f_hat[k=2]=3, one node x=1/8 -> f = 3*cos(2*pi*2*1/8)=3*cos(pi/2)=0
    x = [(mpmath.mpf(1) / 8,)]
    f_hat = [mpmath.mpf(0), mpmath.mpf(0), mpmath.mpf(3), mpmath.mpf(0)]
    f = T.trafo("nfct", [4], 1, x, f_hat)
    assert abs(f[0]) < mpmath.mpf(10) ** -30


def test_ndst_single_term_handcomputed():
    mpmath.mp.dps = 40
    # f_hat indexes k=1..3; coeff for k=2 is 5, node x=1/8 -> 5*sin(2*pi*2/8)=5*sin(pi/2)=5
    x = [(mpmath.mpf(1) / 8,)]
    f_hat = [mpmath.mpf(0), mpmath.mpf(5), mpmath.mpf(0)]  # k=1,2,3
    f = T.trafo("nfst", [4], 1, x, f_hat)
    assert abs(f[0] - 5) < mpmath.mpf(10) ** -30
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /workspaces/nfft && uv run python -m pytest tests/refgen/tests/test_transforms.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'tests.refgen.transforms'` (or `AttributeError`).

- [ ] **Step 3: Write the implementation**

`tests/refgen/transforms.py`:
```python
"""Arbitrary-precision direct NDFT / NDCT / NDST (forward + adjoint).

All functions operate at the current mpmath working precision (mpmath.mp.dps),
which the caller sets. Frequency index sets and ordering reproduce the C f_hat
layout byte-for-byte (see docs/agents/test-methodology.md).
"""
import itertools

import mpmath

_REAL_MODULES = ("nfct", "nfst")


def freqs(module, N):
    """Per-dimension frequency ranges (lists of ints), one list per dimension."""
    if module == "nfft":
        return [list(range(-(n // 2), n - (n // 2))) for n in N]
    if module == "nfct":
        return [list(range(0, n)) for n in N]
    if module == "nfst":
        return [list(range(1, n)) for n in N]
    raise ValueError("unknown module: %r" % (module,))


def index_set(module, N):
    """All frequency tuples, tensor product with dimension 0 varying slowest."""
    return list(itertools.product(*freqs(module, N)))


def nn(module, N):
    """Number of Fourier coefficients."""
    count = 1
    for r in freqs(module, N):
        count *= len(r)
    return count


def _basis(module, k, xj, adjoint):
    """basis(k, x_j). For NFFT the adjoint conjugates (flips the exponent sign)."""
    twopi = 2 * mpmath.pi
    if module == "nfft":
        dot = mpmath.fsum([ki * xji for ki, xji in zip(k, xj)])
        ang = twopi * dot
        return mpmath.expj(ang) if adjoint else mpmath.expj(-ang)
    if module == "nfct":
        prod = mpmath.mpf(1)
        for ki, xji in zip(k, xj):
            prod *= mpmath.cos(twopi * ki * xji)
        return prod
    if module == "nfst":
        prod = mpmath.mpf(1)
        for ki, xji in zip(k, xj):
            prod *= mpmath.sin(twopi * ki * xji)
        return prod
    raise ValueError("unknown module: %r" % (module,))


def _zero(module):
    return mpmath.mpf(0) if module in _REAL_MODULES else mpmath.mpc(0)


def trafo(module, N, M, x, f_hat):
    """Forward: f[j] = sum_k f_hat[k] * basis(k, x[j]).  Returns list of length M."""
    K = index_set(module, N)
    out = []
    for j in range(M):
        acc = _zero(module)
        for idx, k in enumerate(K):
            acc += f_hat[idx] * _basis(module, k, x[j], adjoint=False)
        out.append(acc)
    return out


def adjoint(module, N, M, x, f):
    """Adjoint: f_hat[k] = sum_j f[j] * conj(basis)(k, x[j]).  Length = nn()."""
    K = index_set(module, N)
    out = []
    for k in K:
        acc = _zero(module)
        for j in range(M):
            acc += f[j] * _basis(module, k, x[j], adjoint=True)
        out.append(acc)
    return out
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /workspaces/nfft && uv run python -m pytest tests/refgen/tests/test_transforms.py -q`
Expected: PASS (6 passed).

- [ ] **Step 5: Commit**

```bash
git add tests/refgen/transforms.py tests/refgen/tests/test_transforms.py
git commit -m "test(refgen): arbitrary-precision NDFT/NDCT/NDST with self-tests"
```

### Task 1.3: Number formatting + .txt write/parse (TDD)

**Files:**
- Create: `tests/refgen/io_format.py`
- Test: `tests/refgen/tests/test_io_format.py`

- [ ] **Step 1: Write the failing test**

`tests/refgen/tests/test_io_format.py`:
```python
import mpmath

from tests.refgen import io_format as IO


def test_format_real_has_enough_digits():
    mpmath.mp.dps = 64
    s = IO.fmt_scalar(mpmath.mpf(1) / 3, ndig=50, is_complex=False)
    # at least 40 significant digits, no exponent for an O(1) value
    digits = s.replace("-", "").replace(".", "").lstrip("0")
    assert len(digits) >= 40
    assert "e" not in s and "E" not in s


def test_format_complex_two_columns():
    mpmath.mp.dps = 64
    s = IO.fmt_scalar(mpmath.mpc("0.5", "-0.25"), ndig=40, is_complex=True)
    parts = s.split()
    assert len(parts) == 2
    assert parts[0].startswith("0.5")
    assert parts[1].startswith("-0.25")


def test_write_then_parse_roundtrip_real(tmp_path):
    mpmath.mp.dps = 64
    d, N, M = 1, [4], 2
    x = [(mpmath.mpf("0.1"),), (mpmath.mpf("0.2"),)]
    f_hat = [mpmath.mpf("0.3"), mpmath.mpf("0.4"), mpmath.mpf("0.5")]  # nfst N-1=3
    f = [mpmath.mpf("0.6"), mpmath.mpf("0.7")]
    p = tmp_path / "nfst_1d_4_2.txt"
    IO.write_testcase(str(p), d, N, M, x, f_hat, f, is_complex=False, ndig=50)
    rd, rN, rM, rx, rfh, rf = IO.parse_testcase(str(p), is_complex=False)
    assert (rd, rN, rM) == (d, N, M)
    assert abs(rx[0][0] - x[0][0]) < mpmath.mpf(10) ** -40
    assert abs(rfh[2] - f_hat[2]) < mpmath.mpf(10) ** -40
    assert abs(rf[1] - f[1]) < mpmath.mpf(10) ** -40


def test_write_then_parse_roundtrip_complex(tmp_path):
    mpmath.mp.dps = 64
    d, N, M = 1, [1], 2
    x = [(mpmath.mpf("-0.1"),), (mpmath.mpf("0.3"),)]
    f_hat = [mpmath.mpc("0.2", "0.1")]
    f = [mpmath.mpc("0.4", "-0.5"), mpmath.mpc("-0.6", "0.7")]
    p = tmp_path / "nfft_1d_1_2.txt"
    IO.write_testcase(str(p), d, N, M, x, f_hat, f, is_complex=True, ndig=50)
    rd, rN, rM, rx, rfh, rf = IO.parse_testcase(str(p), is_complex=True)
    assert (rd, rN, rM) == (d, N, M)
    assert abs(rf[0] - f[0]) < mpmath.mpf(10) ** -40
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /workspaces/nfft && uv run python -m pytest tests/refgen/tests/test_io_format.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'tests.refgen.io_format'`.

- [ ] **Step 3: Write the implementation**

`tests/refgen/io_format.py`:
```python
"""Serialization of testcases into the exact CUnit reference-data format.

File layout (blank line between sections, one scalar per line):
    d / N[0..d-1] / M / x (M*d node-major) / f_hat (NN) / f (M)
NFFT scalars are complex ("re im"); NFCT/NFST are real (one value).
"""
import mpmath


def fmt_scalar(value, ndig, is_complex):
    """Format one scalar. Complex -> 're im'; real -> 're'. Fixed-point, no exponent."""
    if is_complex:
        return "%s %s" % (_fmt_real(value.real, ndig), _fmt_real(value.imag, ndig))
    return _fmt_real(value, ndig)


def _fmt_real(v, ndig):
    # mpmath.nstr with no exponent: values here are O(1)..O(M), safe as fixed-point.
    s = mpmath.nstr(mpmath.mpf(v), ndig, strip_zeros=False, min_fixed=-mpmath.inf,
                    max_fixed=mpmath.inf)
    return s


def write_testcase(path, d, N, M, x, f_hat, f, is_complex, ndig):
    lines = []
    lines.append(str(d))
    lines.append("")
    for n in N:
        lines.append(str(n))
    lines.append("")
    lines.append(str(M))
    lines.append("")
    for xj in x:                       # node-major: node j's d coords contiguous
        for coord in xj:
            lines.append(_fmt_real(coord, ndig))
    lines.append("")
    for v in f_hat:
        lines.append(fmt_scalar(v, ndig, is_complex))
    lines.append("")
    for v in f:
        lines.append(fmt_scalar(v, ndig, is_complex))
    lines.append("")
    with open(path, "w") as fp:
        fp.write("\n".join(lines) + "\n")


def parse_testcase(path, is_complex):
    """Inverse of write_testcase (used by tests; mirrors C setup_file token order)."""
    with open(path) as fp:
        toks = fp.read().split()
    it = iter(toks)
    d = int(next(it))
    N = [int(next(it)) for _ in range(d)]
    M = int(next(it))
    x = []
    for _ in range(M):
        x.append(tuple(mpmath.mpf(next(it)) for _ in range(d)))

    def read_vec(count):
        out = []
        for _ in range(count):
            if is_complex:
                re = mpmath.mpf(next(it))
                im = mpmath.mpf(next(it))
                out.append(mpmath.mpc(re, im))
            else:
                out.append(mpmath.mpf(next(it)))
        return out

    # NN is whatever remains minus M values; compute from remaining token count.
    per = 2 if is_complex else 1
    remaining = sum(1 for _ in iter(lambda: next(it, None), None))
    # re-tokenize because the generator above exhausted `it`; simpler: re-split.
    toks2 = open(path).read().split()
    base = 1 + d + 1 + M * d
    tail = toks2[base:]
    total_scalars = len(tail) // per
    nn = total_scalars - M
    it2 = iter(tail)

    def read_vec2(count):
        out = []
        for _ in range(count):
            if is_complex:
                out.append(mpmath.mpc(mpmath.mpf(next(it2)), mpmath.mpf(next(it2))))
            else:
                out.append(mpmath.mpf(next(it2)))
        return out

    f_hat = read_vec2(nn)
    f = read_vec2(M)
    return d, N, M, x, f_hat, f
```

> Note: `parse_testcase` is test-only scaffolding; keep it simple over elegant. The
> double-tokenization above is intentional — it recomputes `NN` from the token count
> so the parser needs no module knowledge.

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /workspaces/nfft && uv run python -m pytest tests/refgen/tests/test_io_format.py -q`
Expected: PASS (4 passed).

- [ ] **Step 5: Commit**

```bash
git add tests/refgen/io_format.py tests/refgen/tests/test_io_format.py
git commit -m "test(refgen): exact reference-data file format with round-trip tests"
```

### Task 1.4: Grids (single source of truth)

**Files:**
- Create: `tests/refgen/grids.py`

- [ ] **Step 1: Write the grids module**

`tests/refgen/grids.py`:
```python
"""The (d, N, M) test-case grids per module. Single source of truth for which
reference files exist. Reproduces the currently-shipped grids so filenames are
preserved. Edit here to add/remove cases, then regenerate."""

# Each entry: (d, N-as-list, M)
GRIDS = {
    "nfft": (
        [(1, [n], m) for n in (1, 2, 4, 10, 20, 50) for m in (1, 10, 20, 50)]
        + [(2, list(N), m) for N in ((10, 10), (10, 20), (20, 10), (20, 20)) for m in (20, 50)]
        + [(3, [10, 10, 10], 10)]
    ),
    "nfct": (
        [(1, [n], m) for n in (1, 2, 4, 10, 25, 50) for m in (1, 10, 25, 50)]
        + [(2, list(N), m) for N in ((10, 10), (10, 25), (25, 10), (25, 25)) for m in (25, 50)]
        + [(3, [10, 10, 10], 10)]
    ),
    "nfst": (
        [(1, [n], m) for n in (2, 4, 10, 25, 50) for m in (1, 10, 25, 50)]
        + [(2, list(N), m) for N in ((10, 10), (10, 25), (25, 10), (25, 25)) for m in (25, 50)]
        + [(3, [10, 10, 10], 10)]
    ),
}

KINDS = ("trafo", "adjoint")  # both forward and adjoint files for every grid entry


def basename(module, kind, d, N, M):
    adj = "_adjoint" if kind == "adjoint" else ""
    return "%s%s_%dd_%s_%d" % (module, adj, d, "_".join(str(n) for n in N), M)
```

- [ ] **Step 2: Verify the grids reproduce the existing committed filename set EXACTLY**

This is the guard against silently changing coverage: the generated filename set (all
modules, forward **and** adjoint, all dims) must equal the committed `tests/data/*.txt`
set. A mismatch means a dropped or extra case. (The CI diff-guard cannot catch an
*orphaned* tracked `.txt`, so this set-equality check is the real safety net.)

Run:
```bash
cd /workspaces/nfft && python -c "
import os
from tests.refgen import grids as G
gen = set()
for m in ('nfft','nfct','nfst'):
    for kind in G.KINDS:
        for (d,N,M) in G.GRIDS[m]:
            gen.add(G.basename(m,kind,d,N,M)+'.txt')
have = {f for f in os.listdir('tests/data') if f.endswith('.txt')}
print('in committed but NOT generated (would be orphaned):', sorted(have - gen))
print('in generator but NOT committed (new cases):', sorted(gen - have))
"
```
Expected: both lists empty. If `have - gen` is non-empty, either add those cases to the
grids or `git rm` the orphans deliberately (and note it). If `gen - have` is non-empty,
the grids added cases — intentional only if you meant to expand coverage.

- [ ] **Step 3: Commit**

```bash
git add tests/refgen/grids.py
git commit -m "test(refgen): per-module test-case grids reproducing shipped filenames"
```

### Task 1.5: C registration + Makefile codegen (TDD)

**Files:**
- Create: `tests/refgen/registration.py`
- Test: add to `tests/refgen/tests/test_io_format.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/refgen/tests/test_io_format.py`:
```python
from tests.refgen import registration as REG


def test_header_contains_decls_and_arrays():
    h = REG.render_header("nfft")
    assert "data/nfft_1d_1_1.txt" in h
    assert "testcase_delegate_file_t nfft_1d_1_1 =" in h
    assert "*testcases_1d_file[]" in h
    assert "*testcases_adjoint_1d_file[]" in h
    assert "*testcases_3d_file[]" in h
    # header guard
    assert "#ifndef" in h and "#endif" in h


def test_extra_dist_lists_all_modules():
    txt = REG.render_extra_dist()
    assert txt.startswith("EXTRA_DIST =")
    assert "nfft_1d_1_1.txt" in txt
    assert "nfct_2d_10_25_50.txt" in txt
    assert "nfst_adjoint_3d_10_10_10_10.txt" in txt
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd /workspaces/nfft && uv run python -m pytest tests/refgen/tests/test_io_format.py -q -k registration`
Expected: FAIL — `No module named 'tests.refgen.registration'`.

- [ ] **Step 3: Write the implementation**

`tests/refgen/registration.py`:
```python
"""Generate the committed C registration header and the data Makefile.am list.

The header replaces the hand-maintained testcase_delegate_file_t declarations and
testcases_*_file[] arrays inside tests/<module>.c. The .c #includes it. This is the
mechanism that eliminates the manual-paste drift that left check_nfsft.nb stale."""

from tests.refgen import grids as G

_MODULES = ("nfft", "nfct", "nfst")


def _decl(module, kind, d, N, M):
    base = G.basename(module, kind, d, N, M)
    return ('static const testcase_delegate_file_t %s = '
            '{setup_file, destroy_file, ABSPATH("data/%s.txt")};' % (base, base))


def _array(module, kind, d):
    # Match the existing C symbols exactly: "adjoint" precedes the dimension.
    # Forward:  testcases_<d>d_file[]          (alias testcases_<d>d_file_)
    # Adjoint:  testcases_adjoint_<d>d_file[]  (alias testcases_adjoint_<d>d_file_)
    prefix = "adjoint_" if kind == "adjoint" else ""
    arr = "testcases_%s%dd_file" % (prefix, d)
    entries = [G.basename(module, kind, d, N, M)
               for (dd, N, M) in G.GRIDS[module] if dd == d]
    lines = ["static const testcase_delegate_file_t *%s[] =" % arr, "{"]
    lines += ["  &%s," % e for e in entries]
    lines += ["};",
              "static const testcase_delegate_t **%s_ = "
              "(const testcase_delegate_t**)%s;" % (arr, arr)]
    return "\n".join(lines)


def render_header(module):
    guard = "NFFT3_TESTS_DATA_GENERATED_%s_TESTCASES_H" % module.upper()
    dims = sorted({d for (d, _N, _M) in G.GRIDS[module]})
    out = [
        "/* GENERATED by tests/refgen — do not edit by hand.",
        " * Regenerate: python -m tests.refgen.generate --module %s" % module,
        " */",
        "#ifndef %s" % guard,
        "#define %s" % guard,
        "",
    ]
    for kind in G.KINDS:
        for d in dims:
            for (dd, N, M) in G.GRIDS[module]:
                if dd == d:
                    out.append(_decl(module, kind, d, N, M))
            out.append("")
            out.append(_array(module, kind, d))
            out.append("")
    out.append("#endif /* %s */" % guard)
    return "\n".join(out) + "\n"


def render_extra_dist():
    files = []
    for module in _MODULES:
        for kind in G.KINDS:
            for (d, N, M) in G.GRIDS[module]:
                files.append(G.basename(module, kind, d, N, M) + ".txt")
    files = sorted(set(files))
    body = "  \\\n".join(files)
    return "EXTRA_DIST = " + body + "\n"
```

- [ ] **Step 4: Run to verify it passes**

Run: `cd /workspaces/nfft && uv run python -m pytest tests/refgen/tests/test_io_format.py -q -k registration`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add tests/refgen/registration.py tests/refgen/tests/test_io_format.py
git commit -m "test(refgen): C registration header + EXTRA_DIST codegen"
```

### Task 1.6: CLI entry point

**Files:**
- Create: `tests/refgen/generate.py`

- [ ] **Step 1: Write the CLI**

`tests/refgen/generate.py`:
```python
"""CLI: regenerate reference data + generated C headers + the data EXTRA_DIST list.

    python -m tests.refgen.generate --module all --precision 64
"""
import argparse
import os
import random

import mpmath

from tests.refgen import grids as G
from tests.refgen import io_format as IO
from tests.refgen import registration as REG
from tests.refgen import transforms as T

_MODULES = ("nfft", "nfct", "nfst")


def _is_complex(module):
    return module == "nfft"


def _draw_inputs(module, kind, d, N, M, rng):
    """Reproducible inputs. Nodes drawn as doubles (exactly representable in any R,
    so the C reader rounds losslessly even at quad); summation is high precision."""
    lo, hi = (-0.5, 0.5) if module == "nfft" else (0.0, 0.5)
    x = [tuple(mpmath.mpf(rng.uniform(lo, hi)) for _ in range(d)) for _ in range(M)]
    count = T.nn(module, N) if kind == "trafo" else M
    if _is_complex(module):
        coeff = [mpmath.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(count)]
    else:
        coeff = [mpmath.mpf(rng.uniform(-1, 1)) for _ in range(count)]
    return x, coeff


def _gen_module(module, precision, seed, data_dir):
    mpmath.mp.dps = precision
    is_c = _is_complex(module)
    for kind in G.KINDS:
        for (d, N, M) in G.GRIDS[module]:
            # Deterministic per-file seed so a single file regenerates identically
            # regardless of iteration order. random.Random(str) uses a stable
            # internal hash (NOT PYTHONHASHSEED), so this is reproducible across
            # processes — unlike the builtin hash() of a tuple.
            rng = random.Random("%d:%s" % (seed, G.basename(module, kind, d, N, M)))
            x, coeff = _draw_inputs(module, kind, d, N, M, rng)
            if kind == "trafo":
                f_hat = coeff
                f = T.trafo(module, N, M, x, f_hat)
            else:
                f = coeff
                f_hat = T.adjoint(module, N, M, x, f)
            path = os.path.join(data_dir, G.basename(module, kind, d, N, M) + ".txt")
            IO.write_testcase(path, d, N, M, x, f_hat, f,
                              is_complex=is_c, ndig=precision)
            print("wrote", path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--module", choices=("nfft", "nfct", "nfst", "all"), default="all")
    ap.add_argument("--precision", type=int, default=64,
                    help="working/printed significant decimal digits (>=34 for quad)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--data-dir", default=os.path.join("tests", "data"))
    ap.add_argument("--header-dir", default=os.path.join("tests", "data", "generated"))
    args = ap.parse_args(argv)

    modules = _MODULES if args.module == "all" else (args.module,)
    os.makedirs(args.header_dir, exist_ok=True)
    for m in modules:
        _gen_module(m, args.precision, args.seed, args.data_dir)
        hpath = os.path.join(args.header_dir, "%s_testcases.h" % m)
        with open(hpath, "w") as fp:
            fp.write(REG.render_header(m))
        print("wrote", hpath)

    # EXTRA_DIST always reflects all modules (the file lists every module's data).
    mk = os.path.join(args.data_dir, "Makefile.am")
    with open(mk, "w") as fp:
        fp.write(REG.render_extra_dist())
    print("wrote", mk)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke-test the CLI into a scratch dir (do not touch tests/data yet)**

Run:
```bash
cd /workspaces/nfft && rm -rf /tmp/refgen_smoke && mkdir -p /tmp/refgen_smoke/generated && \
uv run python -m tests.refgen.generate --module nfct --precision 40 \
  --data-dir /tmp/refgen_smoke --header-dir /tmp/refgen_smoke/generated && \
head -12 /tmp/refgen_smoke/nfct_1d_1_1.txt && \
grep -c "testcase_delegate_file_t" /tmp/refgen_smoke/generated/nfct_testcases.h
```
Expected: prints a well-formed `nfct_1d_1_1.txt` (d=1, N=1, M=1, one node, one f_hat, one f) and a count > 0 for the header decls.

- [ ] **Step 3: Run the full generator test suite**

Run: `cd /workspaces/nfft && uv run python -m pytest tests/refgen/tests -q`
Expected: PASS (all tests).

- [ ] **Step 4: Commit**

```bash
git add tests/refgen/generate.py
git commit -m "test(refgen): CLI to emit reference data, headers, and EXTRA_DIST"
```

---

## Phase 2 — Regenerate & replace; wire generated headers into the C suite

### Task 2.1: Refactor the three .c files to include the generated headers

**Files:**
- Modify: `tests/nfft.c` (remove inline file-testcase decls + `testcases_*_file[]` arrays; add include)
- Modify: `tests/nfct.c` (same)
- Modify: `tests/nfst.c` (same)

- [ ] **Step 1: Generate the headers in place first (so the include target exists)**

Run:
```bash
cd /workspaces/nfft && uv run python -m tests.refgen.generate --module all --precision 64
git status --short tests/data | head
```
Expected: `tests/data/generated/{nfft,nfct,nfst}_testcases.h` created; many `tests/data/*.txt` modified; `tests/data/Makefile.am` modified.

- [ ] **Step 2: In `tests/nfft.c`, delete the generated regions and add the include**

Locate the block beginning at the first `static const testcase_delegate_file_t nfft_1d_1_1 = {setup_file,destroy_file,ABSPATH("data/nfft_1d_1_1.txt")};` and the matching `testcases_*_file[]` / `testcases_*_adjoint_file[]` arrays (for 1d, 2d, 3d, forward and adjoint). Remove **only** those `static const testcase_delegate_file_t <name> = {...};` declarations and the `testcases_<d>d[_adjoint]_file[]` array definitions plus their `_ = (const testcase_delegate_t**)...` aliases. **Keep** the entry functions (`X(check_1d_fast_file)` etc.), the `trafos_*` arrays, and everything else.

Immediately before the first entry function that references `testcases_1d_file_`, add:
```c
#include "data/generated/nfft_testcases.h"
```

> The generated header defines, for each `d`, both `testcases_<d>d_file[]` (+ `_`)
> and `testcases_<d>d_adjoint_file[]` (+ `_`). The entry functions already reference
> those exact names, so no call sites change.

- [ ] **Step 3: Repeat for `tests/nfct.c` and `tests/nfst.c`** with `#include "data/generated/nfct_testcases.h"` and `#include "data/generated/nfst_testcases.h"` respectively.

- [ ] **Step 4: Reconcile array names between header and .c**

Run:
```bash
cd /workspaces/nfft && for m in nfft nfct nfst; do
  echo "== $m =="
  echo "header arrays:"; grep -o 'testcases_[0-9]d[a-z_]*\[\]' tests/data/generated/${m}_testcases.h | sort -u
  echo "used in .c:"; grep -o 'testcases_[0-9]d[a-z_]*_' tests/${m}.c | sort -u
done
```
Expected: every `testcases_*_` symbol used in the `.c` has a matching `testcases_*[]` defined in the header. If the `.c` uses a name the header doesn't emit (e.g. a 4d array, or a differently-suffixed adjoint), fix `render_header`/`grids.py` to emit it and regenerate, OR adjust the `.c` call site — they must agree. Resolve before building.

- [ ] **Step 5: Build the tests**

Run:
```bash
cd /workspaces/nfft && ./bootstrap.sh >/tmp/boot.log 2>&1 && \
./configure --enable-all --enable-openmp >/tmp/conf.log 2>&1 && \
make -j >/tmp/make.log 2>&1; tail -5 /tmp/make.log
```
Expected: build succeeds (no errors). If `setup_file`/`destroy_file`/`ABSPATH` are reported undefined in the header context, move the `#include` to *after* their definitions in the `.c`.

- [ ] **Step 6: Commit**

```bash
git add tests/nfft.c tests/nfct.c tests/nfst.c tests/data/generated/
git commit -m "test: include generated file-testcase headers in nfft/nfct/nfst"
```

### Task 2.2: Verify the regenerated reference data passes the suite

**Files:** none (verification + the regenerated data committed in 2.1/here)

- [ ] **Step 1: Run the full CUnit suite**

Run: `cd /workspaces/nfft && make check 2>&1 | tee /tmp/check.log | tail -30`
Expected: overall `PASSED` / zero failures for the `nfft`, `nfct`, `nfst`, `util` suites.

- [ ] **Step 2: If any file-based case fails, diagnose with the systematic-debugging skill**

The most likely causes and checks:
- **Frequency ordering mismatch** (the f_hat layout differs from C). Symptom: large error on every `*_file` case for one module, but `*_online` passes. Fix: re-examine `transforms.freqs`/`index_set` ordering against the C direct transform's k-loop for that module.
- **NN miscount** (NFST `∏(N-1)`). Symptom: parse/length error or systematic failure on 2d/3d NFST. Fix: confirm `grids`/`transforms.nn`.
- **Sign convention** (forward vs adjoint exponent). Symptom: adjoint files fail, forward pass. Fix: `_basis(adjoint=...)`.
Do **not** loosen the error bounds to make a case pass — the bounds are the contract.

- [ ] **Step 3: Confirm the regenerated data is the only data lineage**

Run: `cd /workspaces/nfft && git status --short tests/data/*.txt | wc -l`
Expected: a nonzero count of modified `.txt` files (the regenerated reference set). Commit them.

- [ ] **Step 4: Commit the regenerated data + Makefile list**

```bash
git add tests/data/*.txt tests/data/Makefile.am
git commit -m "test(data): regenerate nfft/nfct/nfst reference data via tests/refgen"
```

### Task 2.3: Ship generated headers in the dist tarball

**Files:**
- Modify: `tests/Makefile.am`

- [ ] **Step 1: Inspect the current header handling**

Run: `cd /workspaces/nfft && grep -nE 'EXTRA_DIST|noinst_HEADERS|\.h' tests/Makefile.am`
Expected: shows how existing `*.h` test headers are listed.

- [ ] **Step 2: Add the generated headers to EXTRA_DIST**

In `tests/Makefile.am`, add (matching the file's existing style):
```makefile
EXTRA_DIST += data/generated/nfft_testcases.h \
              data/generated/nfct_testcases.h \
              data/generated/nfst_testcases.h
```
(If `EXTRA_DIST` is not yet defined in this file, use `EXTRA_DIST =` for the first occurrence.)

- [ ] **Step 3: Verify dist includes them**

Run: `cd /workspaces/nfft && make distdir 2>/dev/null >/dev/null; ls */tests/data/generated/ 2>/dev/null || find . -path '*tests/data/generated/*_testcases.h' -newer tests/Makefile.am | head`
Expected: the three headers appear under the dist dir (or are found). If `make distdir` is unavailable in this tree, at minimum `make -j` must still pass (Step 4).

- [ ] **Step 4: Rebuild to confirm nothing broke**

Run: `cd /workspaces/nfft && make -j >/tmp/make2.log 2>&1 && make check 2>&1 | tail -5`
Expected: build + `make check` pass.

- [ ] **Step 5: Commit**

```bash
git add tests/Makefile.am
git commit -m "build(tests): distribute generated file-testcase headers"
```

### Task 2.4: Remove the legacy Mathematica transform-data notebooks

**Files:**
- Remove: `tests/check_nfft.{nb,m}`, `tests/check_nfct.{nb,m}`, `tests/check_nfst.{nb,m}`, `tests/check_nfsft.{nb,m}` (8 files)
- **Keep** `tests/PrintVector.m` and `tests/check_bspline.{nb,m}` — see "Kept" in File Structure.

- [ ] **Step 1: Confirm `PrintVector.m` is only needed by the kept bspline notebook**

Run: `cd /workspaces/nfft && grep -rln 'PrintVector' tests/*.m`
Expected: only `tests/check_bspline.m` (and the transform notebooks being removed). This
is **why `PrintVector.m` must be kept** — `check_bspline.m` still loads it.

- [ ] **Step 2: Confirm nothing in the build references the transform notebooks**

Run: `cd /workspaces/nfft && grep -rnE 'check_nf(ft|ct|st|sft)\.(nb|m)' --include=Makefile.am --include=configure.ac --include=*.in . | grep -v matlab/`
Expected: no build-system references (the notebooks were always run manually). If anything appears, stop and reassess.

- [ ] **Step 3: Remove the 8 transform-notebook files (git retains history)**

Run:
```bash
cd /workspaces/nfft && git rm tests/check_nfft.nb tests/check_nfft.m \
  tests/check_nfct.nb tests/check_nfct.m tests/check_nfst.nb tests/check_nfst.m \
  tests/check_nfsft.nb tests/check_nfsft.m
```
Expected: 8 files staged for removal. (`PrintVector.m` and `check_bspline.{nb,m}` are
**not** removed.)

- [ ] **Step 4: Confirm the build + tests still pass without them**

Run: `cd /workspaces/nfft && make check 2>&1 | tail -5`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git commit -m "test: remove legacy Mathematica transform-data notebooks (superseded by tests/refgen)"
```

> Note: `tests/data/Makefile.am` previously listed only `.txt` files in `EXTRA_DIST`,
> so removing the notebooks does not affect `make dist`. `tests/check_nfsft.c`
> (the *orphaned test*, per `CONTEXT.md`) is intentionally left untouched — wiring a
> real NFSFT suite is future work enabled by `docs/agents/test-methodology.md`.
> `tests/check_bspline.{nb,m}` + `tests/PrintVector.m` are intentionally retained:
> the B-spline window-sample data (`tests/bspline.c`) is a separate artifact, slated
> as the first follow-up exercise of the generator's "add a new target" recipe.

---

## Phase 3 — Quad-precision readiness

### Task 3.1: Verify the regenerated data across the full window × precision matrix

The **fast** file-based bound (`err_trafo()`) is a hand-tuned `a`/`b` constant table
indexed by **window** × **`MANT_DIG`**. New random inputs change the *actual* error, so
a marginally-tuned bound that the old draw cleared could now trip. The test compiles
for exactly four windows — `kaiserbessel`(`KAISER_BESSEL`), `gaussian`(`GAUSSIAN`),
`bspline`(`B_SPLINE`), `sinc`(`SINC_POWER`); `dirac` hits `#error Unsupported window
function` (`tests/nfft.c:327`), so it has no suite. Each precision selects a distinct
`MANT_DIG` branch (24/53/64). The full local gate is therefore **4 × 3 = 12 cells**.

The reference `.txt` files are window-independent (pure math); the sweep varies only
the *fast algorithm* and its bound. Out-of-source build trees share the one committed
data set via `ABS_SRCDIR`, so no data is regenerated per cell.

**Files:** none (verification only)

- [ ] **Step 1: Run the 12-cell matrix (serial, out-of-source)**

Run:
```bash
cd /workspaces/nfft
mkdir -p /tmp/sweep && : > /tmp/sweep/summary.txt
for win in kaiserbessel gaussian bspline sinc; do
  for prec in float double long-double; do
    case $prec in
      float)       pflag=--enable-float ;;
      double)      pflag= ;;
      long-double) pflag=--enable-long-double ;;
    esac
    bd="/tmp/sweep/${win}-${prec}"; rm -rf "$bd"; mkdir -p "$bd"
    ( cd "$bd" && \
      /workspaces/nfft/configure --enable-all --with-window="$win" $pflag \
        >configure.log 2>&1 && make -j >make.log 2>&1 && \
      make check >check.log 2>&1 ) \
      && echo "PASS  $win  $prec" >> /tmp/sweep/summary.txt \
      || echo "FAIL  $win  $prec  (see /tmp/sweep/${win}-${prec}/*.log)" >> /tmp/sweep/summary.txt
  done
done
cat /tmp/sweep/summary.txt
```
Expected: 12 lines, all `PASS`.

- [ ] **Step 2: Triage any FAIL**

For a failing cell, inspect `/tmp/sweep/<win>-<prec>/check.log`:
- **Direct file case fails** (`*_direct_file`, bound `48/120/130·ε`): a real generator
  defect (ordering / formula / NN). Fix the generator, regenerate, re-run. Do **not**
  touch bounds.
- **Fast file case fails** (`*_fast_file`, `err_trafo()` bound) by a *small* margin in a
  window/precision whose `a`/`b` carries a `// TODO` or is known-marginal: this is a
  legitimate **bound re-tune** for that `(window, MANT_DIG)` cell in
  `tests/{nfft,nfct,nfst}.c` — adjust the constant minimally and record it in
  ADR-0003. Re-tuning a documented tuning constant is maintenance, **not** a generator
  defect; **never** alter generator output to satisfy a tuned bound.
- **Online case fails:** unrelated to this plan (the online path is untouched) — flag
  as pre-existing.

- [ ] **Step 3: Restore the primary double build tree**

Run:
```bash
cd /workspaces/nfft && ./configure --enable-all --enable-openmp >/tmp/conf-d.log 2>&1 && \
make -j >/tmp/make-d.log 2>&1; tail -3 /tmp/make-d.log
```
Expected: clean double + OpenMP build.

- [ ] **Step 4: Record the matrix result in ADR-0003** (the `/tmp/sweep/summary.txt`
  table and any bound re-tunes). No source commit here unless Step 2 required a bound
  re-tune — in that case commit those `.c` edits:

```bash
git add tests/nfft.c tests/nfct.c tests/nfst.c   # only if bounds were re-tuned
git commit -m "test: re-tune fast-bound constants for regenerated reference data"
```

### Task 3.2: Confirm the generator emits quad-sufficient digits

**Files:**
- Test: add to `tests/refgen/tests/test_io_format.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/refgen/tests/test_io_format.py`:
```python
def test_quad_precision_digits(tmp_path):
    import mpmath
    from tests.refgen import generate as GEN
    # Generate one nfft file at quad-target precision and check digit count.
    GEN.main(["--module", "nfft", "--precision", "40",
              "--data-dir", str(tmp_path), "--header-dir", str(tmp_path)])
    # pick the smallest file
    p = tmp_path / "nfft_1d_1_1.txt"
    toks = p.read_text().split()
    # node value should carry >=34 significant digits (quad)
    node = toks[3]  # d, N, M, then first node
    sig = node.replace("-", "").replace(".", "").lstrip("0")
    assert len(sig) >= 34
```

- [ ] **Step 2: Run to verify it passes (the generator already supports `--precision`)**

Run: `cd /workspaces/nfft && uv run python -m pytest tests/refgen/tests/test_io_format.py -q -k quad`
Expected: PASS. (If it fails because `_fmt_real` truncates below 34 digits, raise the working `mpmath.mp.dps` / `ndig` handling so `--precision 40` yields ≥34 significant figures.)

- [ ] **Step 3: Commit**

```bash
git add tests/refgen/tests/test_io_format.py
git commit -m "test(refgen): assert generator emits quad-sufficient digit counts"
```

### Task 3.3: ADR — scoped `__float128` build gap

**Files:**
- Create: `docs/adr/0003-quad-precision-readiness.md`

- [ ] **Step 1: Write the ADR**

`docs/adr/0003-quad-precision-readiness.md`:
````markdown
# ADR-0003: 128-bit quadruple precision readiness

## Status
Accepted (data-pipeline readiness); the `__float128` build is scoped, not implemented.

## Context
We want the test reference-data pipeline ready for 128-bit quadruple precision and,
later, even wider types. The reference-data generator (`tests/refgen/`, see
[ADR-0002](0002-python-reference-data-generator.md)) now emits ≥64-digit values and
takes `--precision`, so the **data** is quad-ready. Wiring an actual `__float128`
*build* is a separate, larger effort with these touch points:

- `configure.ac`: a new precision (`NFFT_QUAD`, `PRECISION=q`, `PREC_SUFFIX=q`,
  `NFFT_PRECISION_MACRO=NFFT_PRECISION_QUAD`), conflict checks vs single/long-double,
  and `-lquadmath` detection. `MANT_DIG==113` is already accepted at `configure.ac`.
- FFTW: `m4/nfft_lib_fftw3.m4` is fully parameterized by `$PREC_SUFFIX`, so it finds
  `fftw3q`/`fftwq_*` automatically once `PREC_SUFFIX=q` — plus libquadmath in the link.
- `include/infft.h`: an `#elif defined(NFFT_QUAD)` arm for `R`/`C`, the `K(x)`
  `Q`-suffix branch, the `EPSILON`/`MANT_DIG` arm (`FLT128_*`), and the full quad
  math-function macro set (`fabsq`/`sqrtq`/`expq`/`cosq`/`csqrtq`/…). This macro set
  is the largest mechanical change.
- I/O: `printf`/`scanf` cannot handle `__float128`. The `__FI__`-based `fscanf`
  (`tests/*.c`) and `__FE__`-based prints must switch to `strtoflt128` /
  `quadmath_snprintf` (`%Qe`). The current `MANT_DIG==113 → "...LE"` format arms are
  written for 113-bit *long double*, not `__float128`, and would misbehave for true
  quad.
- Public API: `include/nfft3.h` needs `NFFT_MANGLE_QUAD` + a 4th `*_DEFINE_API`
  instantiation per module; `include/nfft3mp.h` needs an `NFFT_PRECISION_QUAD` arm.
- Test bounds: the `MANT_DIG==113`/`64` fast-bound `a`/`b` tables in
  `tests/{nfft,nfct,nfst}.c` carry `// TODO`s and need tuning once a quad build runs.

## Decision
Make the reference-data pipeline quad-ready now (done: generator `--precision`,
≥64-digit files, precision-portable format). Defer the `__float128` build to a
dedicated follow-on, tracked against the touch-point list above.

## Verification evidence
The regenerated data was validated across the full **4 windows
(kaiserbessel/gaussian/bspline/sinc) × 3 precisions (float/double/long-double)**
matrix — the three `MANT_DIG` bound branches (24/53/64) every supported window can
hit. (`dirac` has no CUnit suite: `#error Unsupported window function`.) Any fast-bound
`a`/`b` re-tunes made for the new draw are listed here:

<!-- paste /tmp/sweep/summary.txt and note any (window, MANT_DIG) constants changed -->

The `MANT_DIG==113` (quad) branch cannot be exercised until the `__float128` build
above exists; its `a`/`b` constants remain `// TODO` and must be tuned then.

## Consequences
- Existing double/long-double trees are unaffected.
- A future quad build reuses the same `tests/data/*.txt` files (the C reader rounds
  on input); only the build plumbing and I/O path above must be added.
- The `__FI__`/`__FE__` `__float128` I/O work is the genuine blocker for *running*
  the suite at quad and must be done before the bound tables can be tuned.
````

- [ ] **Step 2: Commit**

```bash
git add docs/adr/0003-quad-precision-readiness.md
git commit -m "docs(adr): scope 128-bit quad-precision readiness and build gap"
```

### Task 3.4: ADR — replace Mathematica with the Python generator

**Files:**
- Create: `docs/adr/0002-python-reference-data-generator.md`

- [ ] **Step 1: Write the ADR**

`docs/adr/0002-python-reference-data-generator.md`:
````markdown
# ADR-0002: Python (mpmath) reference-data generator replaces Mathematica notebooks

## Status
Accepted.

## Context
The high-precision reference data for the `nfft`/`nfct`/`nfst` *file-based checks*
was produced by Mathematica notebooks (`tests/check_*.nb`/`.m`) at 64-digit working
precision, then **manually pasted** into the `.c` files and `tests/data/Makefile.am`.
Problems: (1) requires a Mathematica licence and the GUI; (2) coding agents cannot
read/modify `.nb` blobs; (3) the manual paste drifts — `check_nfsft.nb`/`.m` were a
stale verbatim copy of the NFST generator (Sin kernel, `nfst_` filenames) and never
produced spherical-harmonic data.

## Decision
Replace the test-data notebooks with `tests/refgen/`, a Python + mpmath generator:
- Arbitrary precision via `mpmath` (CLI `--precision`); plain-text, agent-editable.
- Single source of truth: one run emits the `tests/data/*.txt` files, the committed
  C registration headers `tests/data/generated/<module>_testcases.h`, and the
  `tests/data/Makefile.am` `EXTRA_DIST` list — no manual paste, no drift.
- The C build never invokes Python; artifacts are committed (as the notebook outputs
  were). `mpmath`/`pytest` are dev/CI dependencies only. Generation is byte-stable
  (pinned `mpmath`, stable PRNG), so a CI job (`.github/workflows/refgen.yml`)
  regenerates and `git diff --exit-code`s to enforce that committed artifacts never
  drift from the generator — the guard the old manual paste lacked.
- The generator has its own pytest self-tests (equispaced NDFT = DFT, adjoint =
  conjugate transpose, hand-computed cases, format round-trip).

### Input values are drawn as doubles, summed at high precision
Nodes and input coefficients are drawn with a 53-bit (double) PRNG and only the
*summation* is done at the configured arbitrary precision. Rationale: the stored
decimal of a double is representable **exactly** in `long double` and `__float128`,
so the C reader reproduces the generator's input bit-for-bit at *every* precision —
the cross-precision comparison carries **no input-rounding term**, only summation
error, which the existing `48·ε`/`120·ε`/`130·ε` bounds cover. This keeps a future
quad-accuracy measurement (for tuning the `MANT_DIG==113` bound constants) free of an
input-rounding artifact. Drawing full-precision random inputs (as Mathematica did)
also passes, but injects a tiny bounded input perturbation at every precision; we
prefer the cleaner property. Node entropy is irrelevant for these deterministic
accuracy tests, and the transforms are continuous in `x`, so 53-bit nodes lose
nothing.

Mathematica's `SeedRandom[1]` is not byte-reproducible outside Mathematica, so the
regenerated files carry new (equally valid) random inputs; correctness is verified by
`make check` against the unchanged error bounds, not by byte-equality.

The legacy **transform-data** notebooks are removed (git history retains them). Kept,
out of scope for this plan: `tests/check_bspline.{nb,m}` + `tests/PrintVector.m` — also
test-data generators, but for the B-spline window samples in `tests/bspline.c` (a
different artifact: a sampled C array, not a transform `.txt`), slated as the first
follow-up of the "add a new generator target" recipe; and
`kernel/util/{bspline,bessel_i0}.nb`, which are genuine library codegen (the B-spline
evaluator and the minimax `I₀` coefficients), not test data.

## Consequences
- Adding a transform's reference tests is now a code change in `tests/refgen/`
  (transforms + grids + module name) plus the C clone described in
  `docs/agents/test-methodology.md`.
- Quad/octuple readiness is a `--precision` flag (see
  [ADR-0003](0003-quad-precision-readiness.md)).
````

- [ ] **Step 2: Verify ADRs cross-link**

Run: `cd /workspaces/nfft && grep -l '0003-quad' docs/adr/0002-*.md && grep -l '0002-python' docs/adr/0003-*.md`
Expected: both filenames print (cross-links resolve).

- [ ] **Step 3: Commit**

```bash
git add docs/adr/0002-python-reference-data-generator.md
git commit -m "docs(adr): record Python reference-data generator decision"
```

### Task 3.5: Glossary stubs for `/grill-with-docs`

**Files:**
- Modify: `CONTEXT.md`

- [ ] **Step 1: Add a Tests-section stub for the new terms**

In `CONTEXT.md`, under the existing `### Tests` section, append:
```markdown
**Reference-data generator**:
The Python + mpmath tool `tests/refgen/` that is the single source of truth for the
file-based check data: it emits `tests/data/*.txt`, the generated C headers
`tests/data/generated/<module>_testcases.h`, and the data `EXTRA_DIST` list. Not a
build dependency. See ADR-0002.
_Avoid_: data scripts, Mathematica notebooks (removed).

**File-based check / Online check**:
The two CUnit test classes. A **file-based check** reads a high-precision reference
from `tests/data/*.txt` and validates both the direct and fast transform against it.
An **online check** generates random input, builds the reference with the direct
transform, and validates the fast transform against it. See
`docs/agents/test-methodology.md`.
_Avoid_: reference test / accuracy test (ambiguous — name the class).
```

> These are provisional; `/grill-with-docs` will sharpen the wording and placement
> during the adversarial review the user runs next.

- [ ] **Step 2: Commit**

```bash
git add CONTEXT.md
git commit -m "docs(context): add reference-data generator and check-class terms"
```

---

## Phase 4 — CI enforcement (the teeth behind "single source of truth")

### Task 4.1: CI job — generator self-tests + regenerate-and-diff drift guard

**Files:**
- Create: `.github/workflows/refgen.yml`

- [ ] **Step 1: Confirm regeneration is a no-op on the committed tree (local pre-check)**

Run:
```bash
cd /workspaces/nfft && uv sync --locked && \
uv run python -m tests.refgen.generate --module all --precision 64 --seed 1 && \
git diff --stat -- tests/data tests/data/generated
```
Expected: **no diff** — the committed artifacts already equal a fresh run with the
default precision/seed. If there is a diff, the committed data was produced with a
different mpmath version or options; reconcile (re-commit the fresh output, update the
pin in `pyproject.toml`/`uv.lock`) before adding the guard, or the CI job will fail on
its first run.

- [ ] **Step 2: Write the workflow**

`.github/workflows/refgen.yml`:
```yaml
name: refgen
on:
  push:
    paths:
      - 'tests/refgen/**'
      - 'tests/data/**'
      - 'pyproject.toml'
      - 'uv.lock'
      - '.github/workflows/refgen.yml'
  pull_request:
    paths:
      - 'tests/refgen/**'
      - 'tests/data/**'
      - 'pyproject.toml'
      - 'uv.lock'
      - '.github/workflows/refgen.yml'
jobs:
  refgen:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Install uv
        uses: astral-sh/setup-uv@v5
        with:
          enable-cache: true
      - name: Sync locked deps
        run: uv sync --locked
      - name: Generator self-tests
        run: uv run python -m pytest tests/refgen/tests -q
      - name: Regenerate reference data + headers
        run: uv run python -m tests.refgen.generate --module all --precision 64 --seed 1
      - name: Fail if committed artifacts are stale
        run: |
          if ! git diff --exit-code -- tests/data tests/data/generated; then
            echo "::error::Committed reference data/headers are stale. Run:"
            echo "::error::  uv run python -m tests.refgen.generate --module all --precision 64 --seed 1"
            exit 1
          fi
```

- [ ] **Step 3: Validate the workflow locally (lint + dry behaviour)**

Run:
```bash
cd /workspaces/nfft && uv run --with pyyaml python -c "import yaml; yaml.safe_load(open('.github/workflows/refgen.yml')); print('yaml ok')"
```
Expected: `yaml ok`. (`--with pyyaml` pulls PyYAML ephemerally; it is not a project dependency.)

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/refgen.yml
git commit -m "ci: add refgen self-tests and reference-data drift guard"
```

> The guard does not run Python during the C build matrix — it is an independent job,
> consistent with `CONTEXT.md`'s rule that the C build never depends on the generator.

---

## Self-Review

**Spec coverage:**
- "Analyse current test methodology and document the approach" → Phase 0 (`docs/agents/test-methodology.md`) + the verified Background section. ✓
- "Analyse how high-precision reference data was calculated" → captured in Background + ADR-0002 context (Mathematica, 64-digit, manual paste, stale nfsft). ✓
- "Refactor … into arbitrary precision in C or Python so an agent can modify" → Phase 1 (Python + mpmath, decision recorded), single source of truth. ✓
- "128-bit quadruple precision readiness" → Phase 3 (generator `--precision`, ≥34-digit assertion, full 4-window × 3-precision verification matrix proving wider-type consumption, ADR-0003 build gap). ✓
- Apply methodology to other transforms → enabled (not implemented) via the doc's recipe + the `transforms.py`/`grids.py` extension points; explicitly out of scope per the chosen quad/scope answers. ✓

**Placeholder scan:** No "TBD"/"add error handling"/"similar to Task N" — every code step contains full content. The only deliberately-deferred item is the `__float128` build (ADR-0003, by user decision). Glossary wording was finalized during the `/grill-with-docs` review and is captured verbatim in Task 3.5.

**Type/name consistency:** Generator symbols are consistent across modules: `transforms.freqs/index_set/nn/trafo/adjoint`, `io_format.fmt_scalar/write_testcase/parse_testcase`, `grids.GRIDS/KINDS/basename`, `registration.render_header/render_extra_dist`, `generate.main`. The C side reuses the **existing** `testcases_<d>d_file[]` (forward) and `testcases_adjoint_<d>d_file[]` (adjoint) names — verified against `tests/nfft.c:903` and `:976` (note: "adjoint" precedes the dimension; `registration._array` was corrected to match). Task 2.1 Step 4 reconciles header-emitted names against `.c` call sites before building.

**Known risk flagged for execution:** the generated `testcases_*[]` array set must exactly match what each `.c` already references (including any module-specific adjoint suffix or absence of a 4d file array). Task 2.1 Step 4 is the gate; resolve there before `make`.

---

## Execution Handoff

**Plan complete, reviewed via `/grill-with-docs`, and saved to
`docs/superpowers/plans/2026-06-13-reference-data-generator-and-quad-readiness.md`.**

Decisions resolved during the adversarial review (all reflected above):
- **Generator:** Python + mpmath; single source of truth (data + C headers + `EXTRA_DIST`).
- **Existing data:** regenerate & replace; correctness gate is `make check`, not byte-equality.
- **Inputs:** drawn as doubles, summed at high precision (lossless cross-precision read; rationale in ADR-0002).
- **Verification:** full 4-window × 3-precision (12-cell) matrix locally, not just the default.
- **Enforcement:** CI drift guard (`refgen.yml`) + pinned `mpmath`.
- **Quad:** data-pipeline readiness + ADR-0003 scoping the `__float128` build gap (not built here).
- **B-spline:** out of scope (Option A); `PrintVector.m` + `check_bspline.{nb,m}` kept; queued as the first follow-up.
- Two generator bugs caught & fixed in the plan: adjoint array naming (`testcases_adjoint_<d>d_file`) and reproducible per-file seeding.

Choose an execution approach:

1. **Subagent-Driven (recommended)** — fresh subagent per task, review between tasks.
2. **Inline Execution** — execute tasks in this session with checkpoints.
