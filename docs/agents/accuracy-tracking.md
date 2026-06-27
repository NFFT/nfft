# Accuracy tracking with Bencher

Every CUnit case in the `nfft` / `nfct` / `nfst` suites computes a relative error
`err` and compares it to an analytic `bound` (`tests/nfft.c:393` and the clones).
The pass/fail gate stays in C (`err < bound`); Bencher additionally tracks the
*value* over time so slow drift is visible before a case crosses its bound. See
[ADR-0004](../adr/0004-accuracy-tracking-with-bencher.md) and the
*Accuracy tracking* vocabulary in [`CONTEXT.md`](../../CONTEXT.md).

## Pipeline

1. **Emit (C).** `bench_emit_accuracy` (`tests/bench_emit.c`) appends one raw
   NDJSON record per case to `$NFFT_BENCH_OUT`. No-op when unset, so `make check`
   is unaffected.
2. **Aggregate (Python).** `tests/bench/ndjson_to_bmf.py` groups by the
   *error-shaping parameters* and collapses the *bound-absorbed* `N`/`M` via `max`,
   emitting BMF with two measures per metric: `tightness-ratio` = `max(err/bound)`
   (primary) and `max-error` = `max(err)` (secondary). Metric name:
   `<module>/<oracle>/<speed>/<direction>/<dim>d/<init-slug>`.
3. **Upload (CI).** `.github/workflows/bench-accuracy-linux.yml` runs the serial
   `tests/checkall` across the window×precision matrix (each cell a testbed) and
   calls `bencher run --adapter json --file ...` track-only (no thresholds, no
   `--err`), `develop` baseline, start-point recipe on PRs.

## Run it locally

```bash
NFFT_BENCH_OUT="$PWD/tests/accuracy.ndjson" tests/checkall > /dev/null
uv run python tests/bench/ndjson_to_bmf.py tests/accuracy.ndjson tests/accuracy.bmf.json
bencher run --dry-run --project nfft-accuracy --branch local \
  --testbed local --adapter json --file tests/accuracy.bmf.json
```

## Conventions / gotchas

- **`checkall` is a `check_PROGRAM`** — `make all` does not build it; use
  `make -C tests checkall` (or `make check`).
- **Track-only, phased.** No thresholds yet; never `--err`. `tightness-ratio` is
  lower-is-better → add an *upper* boundary when alerts come later.
- **Determinism.** Serial `tests/checkall`, fixed `SEED`. Not `checkall_threads`.
- **Metric-name stability.** Changing the grouping key breaks Bencher history.
- **Speed axis granularity.** The `fast` speed axis does not distinguish the
  dimension-specialized kernels (`trafo_1d/2d/3d`, `adjoint_1d/2d/3d`) from the
  generic guru transform (`trafo`/`adjoint`); within a group their errors are
  `max`-merged. The `init` variant usually separates them in practice (the
  dedicated 1D path pairs with `init_1d`, the guru path with `init_guru …`).
- **Scope.** Only nfft/nfct/nfst emit today. To add a module, give its harness the
  same `#include "bench_emit.h"` + one `bench_emit_accuracy(...)` call.
- **Quota dial.** Merge file+online in `group_key` (180 → 108 metrics/cell) first.
