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
  were). `mpmath`/`pytest` are dev-only dependencies, declared inline on the
  `uv run --with …` command (pinned `mpmath==1.3.0`) — no `pyproject.toml` or lockfile
  to maintain. Generation is byte-stable (pinned `mpmath`, stable PRNG), and the
  generator carries its own pytest self-tests (run locally with
  `uv run --with mpmath==1.3.0 --with pytest python -m pytest tests/refgen/tests`).
  (No CI workflow is wired up for now; a regenerate-and-`git diff --exit-code` drift
  guard could be added later to enforce that the committed artifacts stay in step
  with the generator.)
- The generator has its own pytest self-tests (equispaced NDFT = DFT, adjoint =
  conjugate transpose, hand-computed cases, format round-trip).

### Input values are drawn as single-precision floats, summed at high precision
Nodes and input coefficients are drawn with the PRNG and **rounded to 24-bit single
precision (`float`)**; only the *summation* is done at the configured arbitrary
precision. Rationale: a `float` value is representable **exactly** in `float`,
`double`, `long double`, and `__float128`, so the C reader reproduces the generator's
input bit-for-bit at *every* precision the library builds — **including the `float`
build**. The cross-precision comparison therefore carries **no input-rounding term**,
only summation error, which the existing `48·ε`/`120·ε`/`130·ε` bounds cover.
(Drawing as `double`, as we did initially, has this property only for `long double`
and `__float128`; the `float` build then has to round the stored double input down to
single precision, injecting an input-rounding term that the `float` tests would carry.
Drawing as `float` is the smallest draw precision that is exact in all four types, so
it cleanly covers the `float` suite too.) This also keeps a future quad-accuracy
measurement (for tuning the `MANT_DIG==113` bound constants) free of an input-rounding
artifact. Node entropy is irrelevant for these deterministic accuracy tests, and the
transforms are continuous in `x`, so 24-bit nodes lose nothing.

Mathematica's `SeedRandom[1]` is not byte-reproducible outside Mathematica, so the
regenerated `.txt` carry new (equally valid) random inputs; correctness is verified by
`make check` against the unchanged error bounds, not by byte-equality. The generator,
its headers, and the build wiring landed first; the shipped `tests/data/*.txt` were
then regenerated from it in a follow-up change (filenames unchanged), so the data is
now produced from the Python source rather than the retired notebooks.

The legacy **transform-data** notebooks are moved to `tests/legacy/` (slated for
removal in a later change; git history also retains them). Kept,
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
