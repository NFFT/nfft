# Bencher Accuracy Tracking Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **Revision (2026-06-29) — implemented, with two changes from the original Task 4/§measures below.** The authoritative description is now [ADR-0004](../../adr/0004-accuracy-tracking-with-bencher.md) and [`docs/agents/accuracy-tracking.md`](../../agents/accuracy-tracking.md):
> 1. **No separate workflow.** Instead of a dedicated `bench-accuracy-linux.yml` that re-runs the 12-cell test matrix, emission piggybacks on the **existing `make check`** in `build-linux.yml` (set `NFFT_BENCH_OUT`); each gcc cell publishes an `accuracy-bmf-<BUILD_CONFIG>` artifact with no secret, and a separate **environment-gated `bencher-upload` job** is the only step that touches the API key (PRs gated, `develop` pushes automatic). Tests run once.
> 2. **Serial + OpenMP both tracked** as a `runtime` (serial/omp) axis — `checkall_threads` routes to `<file>.threads` and tags `openmp:1`; the metric name gains a `/<runtime>/` segment. Bencher project is **`nfft`** (not `nfft-accuracy`).
>
> Tasks 2–3 below are accurate; Task 1 (Bencher project/secret) and Task 5 (docs) stand; Task 4's standalone-workflow steps are superseded by the above.

**Goal:** Capture the per-case accuracy figures (`err` vs `bound`) the NFFT/NFCT/NFST CUnit suites already compute, aggregate them by error-shaping parameters, and track them over time in Bencher Cloud — track-only, never changing the existing pass/fail gate.

**Architecture:** `check_single` (`tests/nfft.c:393`, `tests/nfct.c:392`, `tests/nfst.c:398`) computes a relative error `err` and an analytic `bound` per case. We add an *opt-in* emitter (`tests/bench_emit.{c,h}`) that, when `NFFT_BENCH_OUT` is set, appends **one raw, structured NDJSON record per case**. A pure-stdlib Python converter then does **all grouping and aggregation**: it collapses the *bound-absorbed* parameters (`N`, `M`) via `max` within each combination of *error-shaping* parameters (oracle file/online, direct/fast, forward/adjoint, dimension, init variant) and emits [Bencher Metric Format (BMF)](https://bencher.dev/docs/reference/bencher-metric-format/) JSON with two measures per group. A new GitHub Actions workflow runs the serial suite across the window×precision matrix (each cell an *accuracy testbed*) and uploads via `bencher run` — **track-only** (no thresholds, no `--err`).

**Tech Stack:** C99 (CUnit harness), Autotools + CMake (both updated), Python 3.13 stdlib (converter, no third-party deps), GitHub Actions, Bencher Cloud + `bencher run`.

**Design decisions (resolved in the grilling session — see [ADR-0004](../../adr/0004-accuracy-tracking-with-bencher.md)):**
- **Aggregation, not per-case.** Per-case is ~1680 metrics/cell (measured) → ~20k/run; illegible and quota-heavy. Instead one **accuracy metric** per error-shaping combination, `N`/`M` collapsed via `max`.
- **Two measures per group:** `tightness-ratio` = `max(err/bound)` (primary; distance-to-failure), `max-error` = `max(err)` (secondary; absolute).
- **Split file vs online** (different oracle → error-shaping): preserves the large-`N` online signal from being masked. ~**180 metrics/cell** × 2 measures × 12 testbeds ≈ **4,300 metrics/run**.
- **All 12 gcc cells** report on **push** (baseline) and on **gated PRs** (comparison comment); `develop` is the Bencher baseline; PRs use the standard start-point recipe.
- **Track-only, phased:** pure recording now (no thresholds); non-blocking alerts later once baselines settle; hard gating only by explicit later decision.
- **Aggregation lives in the converter**, not C — the emitter stays dumb and emits raw per-case rows.

## Global Constraints

- **Precision-agnostic.** The emitter compiles and behaves identically in float/double/long-double. It takes `long double` (cast at the call site) and never touches `R`/`C`. The build matrix must stay green.
- **Zero behavior change when disabled.** When `NFFT_BENCH_OUT` is unset/empty, `bench_emit_accuracy` does nothing beyond one `getenv`. `make check` output, exit codes, and the `CU_ASSERT` gate are unchanged. Hard requirement — verify it.
- **`checkall` is a `check_PROGRAM`.** `make -j` / `make all` does **not** build it; only `make check` (which also runs it) or `make -C tests checkall` does. Build it explicitly before running it with the env var.
- **Track-only.** `bencher run` is invoked **without** `--err` and **without** thresholds. CUnit's `err < bound` remains the only gate.
- **Metric-name stability.** The grouping key fixes the Bencher series identity (like CodSpeed benchmark names). Do not change the error-shaping / bound-absorbed split casually — it breaks historical continuity.
- **Names exclude window/precision** (those are the *accuracy testbed* = `BUILD_CONFIG`), so the same metric-name set appears in every testbed.
- **Pin all GitHub Actions to a commit SHA** (repo convention). Reuse SHAs already in `.github/workflows/build-linux.yml`.
- **Security model mirrors `bench-linux.yml`:** `pull_request` (never `pull_request_target`); PR runs gated behind the protected `benchmarks` environment; pushes to default branches run automatically; fork PRs (empty secret) fall back to `--dry-run`.
- **Scope = NFFT/NFCT/NFST** only. Other modules are deliberate follow-up.

## Risks & Notes (grounded in a real exhaustive run)

- **Measured counts (double/kaiserbessel, `--enable-exhaustive-unit-tests`):** 1680 emit lines/cell → after aggregation **180 metrics/cell**, ×2 measures ×12 testbeds ≈ **4,320 metrics/run**. Confirm Bencher Cloud quota covers this against expected runs/month (**Task 1 Step 2**). First dial if quota bites: collapse file+online back into one group (one-line converter change → 108/cell), then reduce PR-scope.
- **Low-bit noise.** `-ffast-math`, runner CPU differences jitter the least-significant bits of `err`. We run the **serial** `tests/checkall` with the fixed test `SEED` for determinism, and stay track-only (no thresholds) so noise can't cause false alerts or failures.
- **`tightness-ratio` is lower-is-better.** Not needed now (no thresholds); when alerts are added later, configure an **upper** boundary on `tightness-ratio` (a regression = ratio creeps toward 1.0).

---

### Task 1: Bencher project prerequisites (manual, one-time)

**Files:** none (account/secret setup).

**Interfaces:**
- Produces: a Bencher project slug `nfft-accuracy`, a repository secret `BENCHER_API_TOKEN`, and confirmation the protected `benchmarks` environment exists (already used by `bench-linux.yml`).

- [ ] **Step 1: Create the Bencher org/project and API token**

Sign in at https://bencher.dev, create/reuse an org, create a project with slug `nfft-accuracy`, generate a user API token. Verify with the CLI:
```bash
# Install per https://bencher.dev/docs/how-to/install-cli/ , then:
export BENCHER_API_TOKEN=<token>
bencher project view nfft-accuracy
```
Expected: JSON describing `nfft-accuracy` (not 404/401).

- [ ] **Step 2: Confirm the Cloud quota is adequate**

Check the org's monthly metric quota against the estimate in *Risks* (~4,320 metrics/run × runs/month). If inadequate, decide the fallback (merge file+online; reduce PR scope; paid/self-hosted) before proceeding, and record it in the PR description.

- [ ] **Step 3: Add the repository secret**
```bash
gh secret set BENCHER_API_TOKEN --body "<token>"
gh secret list | grep BENCHER_API_TOKEN
```
Expected: `BENCHER_API_TOKEN` listed.

- [ ] **Step 4: Confirm the `benchmarks` environment exists**
```bash
gh api "repos/{owner}/{repo}/environments" --jq '.environments[].name'
```
Expected: output includes `benchmarks`. If absent, create it (repo Settings → Environments) with the same required reviewer(s) as the CodSpeed gate.

---

### Task 2: C accuracy emitter (raw structured NDJSON) and harness wiring

**Files:**
- Create: `tests/bench_emit.h`
- Create: `tests/bench_emit.c`
- Modify: `tests/nfft.c` (include near line 13; emit call after line 396)
- Modify: `tests/nfct.c` (include near line 13; emit call after line 395)
- Modify: `tests/nfst.c` (include near line 13; emit call after line 401)
- Modify: `tests/Makefile.am:37` (`checkall_SOURCES`)
- Modify: `tests/CMakeLists.txt` (`NFFT_TEST_SOURCES`)

**Interfaces:**
- Produces: `void bench_emit_accuracy(const char *module, const char *oracle, int d, const int *N, int M, const char *init_name, const char *trafo_name, long double accuracy, long double bound, int ok);` — appends one NDJSON record per case to `$NFFT_BENCH_OUT`, or does nothing if unset/empty. Record shape (consumed by Task 3):
  ```json
  {"module":"nfft","oracle":"file","dim":1,"N":[16],"M":8,"init":"init_guru (PRE PSI)","trafo":"trafo_direct","accuracy":1.5e-14,"bound":4.8e-13,"ok":1}
  ```
- Consumes (existing locals in `check_single`): `d`, `int *N`, `M`, the `testcase` pointer, `init_delegate->name`, `trafo_delegate->name`, `err`, `bound`, `ok`, and the static `setup_online` symbol (file vs online by pointer identity). `init`/`trafo` names contain only `[A-Za-z0-9_ ()]`, so they need no JSON escaping.

- [ ] **Step 1: Write the emitter header**

Create `tests/bench_emit.h`:
```c
/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef NFFT_TESTS_BENCH_EMIT_H
#define NFFT_TESTS_BENCH_EMIT_H

/* Append one raw accuracy datapoint as an NDJSON record to the file named by
 * the NFFT_BENCH_OUT environment variable. A no-op (beyond a getenv) when that
 * variable is unset or empty, so ordinary `make check` runs are unaffected.
 *
 * Grouping and aggregation are NOT done here -- the converter
 * (tests/bench/ndjson_to_bmf.py) collapses N/M and builds the metric names.
 * `oracle` is "file" or "online". Values are emitted as long double with ample
 * precision and are valid JSON numbers. */
void bench_emit_accuracy(const char *module, const char *oracle, int d,
                         const int *N, int M, const char *init_name,
                         const char *trafo_name, long double accuracy,
                         long double bound, int ok);

#endif /* NFFT_TESTS_BENCH_EMIT_H */
```

- [ ] **Step 2: Write the emitter implementation**

Create `tests/bench_emit.c`:
```c
/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>

#include "bench_emit.h"

void bench_emit_accuracy(const char *module, const char *oracle, int d,
                         const int *N, int M, const char *init_name,
                         const char *trafo_name, long double accuracy,
                         long double bound, int ok)
{
  const char *path = getenv("NFFT_BENCH_OUT");
  FILE *fp;
  int j;

  if (path == NULL || path[0] == '\0')
    return;

  fp = fopen(path, "a");
  if (fp == NULL)
    return;

  fprintf(fp, "{\"module\": \"%s\", \"oracle\": \"%s\", \"dim\": %d, \"N\": [",
          module, oracle, d);
  for (j = 0; j < d; j++)
    fprintf(fp, "%s%d", j ? "," : "", N[j]);
  /* %.20Le yields a valid JSON number (no leading '+'/spaces; err,bound >= 0). */
  fprintf(fp, "], \"M\": %d, \"init\": \"%s\", \"trafo\": \"%s\", "
              "\"accuracy\": %.20Le, \"bound\": %.20Le, \"ok\": %d}\n",
          M, init_name, trafo_name, accuracy, bound, ok);

  fclose(fp);
}
```

- [ ] **Step 3: Wire the emitter into all three harnesses (include + call)**

In `tests/nfft.c`, add after line 13 (`#include "nfft.h"`):
```c
#include "bench_emit.h"
```
Then in `check_single`, immediately after the existing `printf(" -> %-4s " __FE__ ...` line (`tests/nfft.c:396`), before the block's closing `}`, insert:
```c
    bench_emit_accuracy("nfft",
      testcase->setup == setup_online ? "online" : "file",
      d, N, M, init_delegate->name, trafo_delegate->name,
      (long double)err, (long double)bound, ok);
```
Repeat **identically** in `tests/nfct.c` (include after line 13 `#include "nfct.h"`; call after the printf at `tests/nfct.c:395`) with module `"nfct"`, and in `tests/nfst.c` (include after line 13 `#include "nfst.h"`; call after the printf at `tests/nfst.c:401`) with module `"nfst"`. All referenced locals (`testcase`, `setup_online`, `d`, `N`, `M`, `init_delegate`, `trafo_delegate`, `err`, `bound`, `ok`) exist identically in each file.

- [ ] **Step 4: Add the new module to both build systems**

In `tests/Makefile.am:37`, change:
```make
checkall_SOURCES = check.c util.c util.h reflect.c reflect.h bspline.c bspline.h bessel.c bessel.h nfft.c nfft.h $(NFCT_SOURCES) $(NFST_SOURCES)
```
to:
```make
checkall_SOURCES = check.c util.c util.h bench_emit.c bench_emit.h reflect.c reflect.h bspline.c bspline.h bessel.c bessel.h nfft.c nfft.h $(NFCT_SOURCES) $(NFST_SOURCES)
```

In `tests/CMakeLists.txt`, add `bench_emit.c bench_emit.h` to `NFFT_TEST_SOURCES`, right after the `util.c util.h` line:
```cmake
  util.c    util.h
  bench_emit.c bench_emit.h
  reflect.c reflect.h
```

- [ ] **Step 5: Build the test program (Autotools)**

```bash
./bootstrap.sh
./configure --enable-all --enable-tests --enable-exhaustive-unit-tests
make -j                  # builds the library
make -C tests checkall   # checkall is a check_PROGRAM; `make all` does NOT build it
```
Expected: clean build; `tests/checkall` exists.

- [ ] **Step 6: Verify the disabled path is a no-op (gate unchanged)**
```bash
make check
```
Expected: same PASS summary as before this task; **no** `tests/accuracy.ndjson` created (env unset).

- [ ] **Step 7: Verify enabled emission produces valid structured NDJSON**
```bash
rm -f /tmp/acc.ndjson
NFFT_BENCH_OUT=/tmp/acc.ndjson tests/checkall > /dev/null
echo "lines: $(wc -l < /tmp/acc.ndjson)"   # expect ~1680 for the exhaustive double/kaiserbessel cell
python3 - <<'PY'
import json
recs=[json.loads(l) for l in open("/tmp/acc.ndjson")]
need={"module","oracle","dim","N","M","init","trafo","accuracy","bound","ok"}
assert all(need<=set(r) for r in recs), "missing fields"
assert any(r["oracle"]=="online" for r in recs), "no online rows"
assert any("direct" in r["trafo"] for r in recs), "no direct rows"
print("ok:", len(recs), "records, all fields present")
PY
```
Expected: `~1680` lines; `ok: <n> records, all fields present`.

- [ ] **Step 8: Commit**
```bash
git add tests/bench_emit.c tests/bench_emit.h tests/nfft.c tests/nfct.c tests/nfst.c tests/Makefile.am tests/CMakeLists.txt
git commit -m "test: emit per-case accuracy NDJSON behind NFFT_BENCH_OUT"
```

---

### Task 3: NDJSON → aggregated BMF converter (Python, TDD)

**Files:**
- Create: `tests/bench/ndjson_to_bmf.py`
- Test: `tests/bench/test_ndjson_to_bmf.py`

**Interfaces:**
- Consumes: the structured NDJSON records from Task 2.
- Produces: `convert(records) -> dict` returning BMF where each key is a metric name `"<module>/<oracle>/<speed>/<direction>/<dim>d/<init-slug>"` and each value is `{"tightness-ratio": {"value": max(err/bound)}, "max-error": {"value": max(err)}}`. Helpers `read_ndjson(text) -> list[dict]`, `group_key(rec) -> tuple`, `metric_name(key) -> str`, `slug(s) -> str`. CLI: `python tests/bench/ndjson_to_bmf.py <in.ndjson> <out.bmf.json>`. Derivations: `speed = "direct" if "direct" in trafo else "fast"`; `direction = "adjoint" if trafo.startswith("adjoint") else "forward"`. A record with `bound <= 0` raises `ValueError` (the analytic bound is always `> 0`; this catches emitter bugs).

- [ ] **Step 1: Write the failing tests**

Create `tests/bench/test_ndjson_to_bmf.py`:
```python
import pytest

from ndjson_to_bmf import convert, read_ndjson, group_key, metric_name, slug


def _rec(**kw):
    base = {"module": "nfft", "oracle": "file", "dim": 1, "N": [16], "M": 8,
            "init": "init_guru ()", "trafo": "trafo_direct",
            "accuracy": 1e-14, "bound": 1e-13, "ok": 1}
    base.update(kw)
    return base


def test_slug_normalizes_init_name():
    assert slug("init_guru (PRE PSI)") == "init_guru_PRE_PSI"


def test_group_key_collapses_N_and_M():
    a = group_key(_rec(N=[16], M=8))
    b = group_key(_rec(N=[64], M=99))
    assert a == b  # N and M are bound-absorbed -> same group


def test_speed_and_direction_derived_from_trafo():
    assert metric_name(group_key(_rec(trafo="trafo"))) == \
        "nfft/file/fast/forward/1d/init_guru"
    assert metric_name(group_key(_rec(trafo="adjoint_direct"))) == \
        "nfft/file/direct/adjoint/1d/init_guru"


def test_convert_takes_max_ratio_and_max_error():
    recs = [_rec(N=[16], accuracy=1e-14, bound=1e-13),   # ratio 0.1
            _rec(N=[64], accuracy=5e-13, bound=1e-12)]   # ratio 0.5, err 5e-13
    bmf = convert(recs)
    (name, measures), = bmf.items()
    assert measures["tightness-ratio"]["value"] == pytest.approx(0.5)
    assert measures["max-error"]["value"] == pytest.approx(5e-13)


def test_file_and_online_are_separate_metrics():
    bmf = convert([_rec(oracle="file"), _rec(oracle="online")])
    assert len(bmf) == 2


def test_nonpositive_bound_raises():
    with pytest.raises(ValueError, match="bound"):
        convert([_rec(bound=0.0)])


def test_read_ndjson_skips_blank_and_reports_bad_line():
    assert len(read_ndjson('{"module":"nfft"}\n\n')) == 1
    with pytest.raises(ValueError, match="line 1"):
        read_ndjson("not json\n")
```

- [ ] **Step 2: Run the tests to verify they fail**
```bash
cd tests/bench && uv run --with pytest python -m pytest -v
```
Expected: FAIL — `ModuleNotFoundError: No module named 'ndjson_to_bmf'`.

- [ ] **Step 3: Write the converter**

Create `tests/bench/ndjson_to_bmf.py`:
```python
"""Convert per-case accuracy NDJSON (emitted by the CUnit harness via
NFFT_BENCH_OUT) into aggregated Bencher Metric Format (BMF) JSON.

Each input line is one raw case:
    {"module","oracle","dim","N","M","init","trafo","accuracy","bound","ok"}

Output is one BMF object per *accuracy metric* -- a combination of the
error-shaping parameters (module, oracle file/online, speed direct/fast,
direction forward/adjoint, dimension, init variant) -- with the bound-absorbed
parameters (N, M) collapsed via max:

    {"<module>/<oracle>/<speed>/<direction>/<dim>d/<init-slug>": {
        "tightness-ratio": {"value": max(err/bound)},   # primary
        "max-error":       {"value": max(err)}}, ...}    # secondary
"""
from __future__ import annotations

import argparse
import json
import re
import sys


def slug(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")


def group_key(rec):
    trafo = rec["trafo"]
    speed = "direct" if "direct" in trafo else "fast"
    direction = "adjoint" if trafo.startswith("adjoint") else "forward"
    return (rec["module"], rec["oracle"], speed, direction,
            int(rec["dim"]), slug(rec["init"]))


def metric_name(key):
    module, oracle, speed, direction, dim, init = key
    return f"{module}/{oracle}/{speed}/{direction}/{dim}d/{init}"


def convert(records):
    agg = {}  # key -> [max_ratio, max_err]
    for rec in records:
        bound = float(rec["bound"])
        if bound <= 0.0:
            raise ValueError(f"non-positive bound in record: {rec!r}")
        err = float(rec["accuracy"])
        ratio = err / bound
        key = group_key(rec)
        if key not in agg:
            agg[key] = [ratio, err]
        else:
            agg[key][0] = max(agg[key][0], ratio)
            agg[key][1] = max(agg[key][1], err)
    return {
        metric_name(key): {
            "tightness-ratio": {"value": r},
            "max-error": {"value": e},
        }
        for key, (r, e) in agg.items()
    }


def read_ndjson(text):
    records = []
    for lineno, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid JSON on line {lineno}: {exc}") from exc
    return records


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", help="NDJSON file written via NFFT_BENCH_OUT")
    ap.add_argument("output", help="BMF JSON output path")
    args = ap.parse_args(argv)

    with open(args.input, encoding="utf-8") as f:
        records = read_ndjson(f.read())
    bmf = convert(records)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(bmf, f, indent=2, sort_keys=True)
    print(f"wrote {len(bmf)} accuracy metrics to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run the tests to verify they pass**
```bash
cd tests/bench && uv run --with pytest python -m pytest -v
```
Expected: PASS (7 passed).

- [ ] **Step 5: Smoke-test against real harness output and confirm the count**
```bash
[ -s /tmp/acc.ndjson ] || NFFT_BENCH_OUT=/tmp/acc.ndjson tests/checkall > /dev/null
uv run python tests/bench/ndjson_to_bmf.py /tmp/acc.ndjson /tmp/acc.bmf.json
python3 -c 'import json; d=json.load(open("/tmp/acc.bmf.json")); print(len(d),"metrics"); k=sorted(d)[0]; print(k,"->",d[k])'
```
Expected: **~180 metrics** (double/kaiserbessel exhaustive cell); a sample entry shaped `{'tightness-ratio': {'value': ...}, 'max-error': {'value': ...}}`.

- [ ] **Step 6: Commit**
```bash
git add tests/bench/ndjson_to_bmf.py tests/bench/test_ndjson_to_bmf.py
git commit -m "test: aggregate accuracy NDJSON into BMF (ratio + raw, max over N/M)"
```

---

### Task 4: CI workflow — build matrix, convert, upload (track-only)

**Files:**
- Create: `.github/workflows/bench-accuracy-linux.yml`

**Interfaces:**
- Consumes: `tests/checkall` (Task 2), `tests/bench/ndjson_to_bmf.py` (Task 3), the `BENCHER_API_TOKEN` secret and `benchmarks` environment (Task 1).
- Produces: per-cell uploads to project `nfft-accuracy`, testbed = `BUILD_CONFIG`, two measures per metric. Baseline branch `develop`; PRs compared via start-point.

- [ ] **Step 1: Resolve the Bencher action SHA to pin**
```bash
# Reuse from build-linux.yml (already pinned there):
#   actions/checkout   -> de0fac2e4500dabe0009e67214ff5f5447ce83dd  (v6)
#   astral-sh/setup-uv -> 08807647e7069bb48b6ef5acd8ec9567f424441b  (v7)
# Resolve the Bencher CLI action's latest release SHA:
gh api repos/bencherdev/bencher/tags --jq '.[0] | "\(.name) \(.commit.sha)"'
```
Expected: a `vX.Y.Z <40-hex-sha>` line; use that SHA for `__BENCHER_SHA__` in Step 2.

- [ ] **Step 2: Write the workflow**

Create `.github/workflows/bench-accuracy-linux.yml`:
```yaml
name: Bencher Accuracy

on:
  # SECURITY INVARIANT: `pull_request`, NEVER `pull_request_target`. Fork PRs run
  # with a read-only token and an empty BENCHER_API_TOKEN -> the upload step falls
  # back to --dry-run. Same-repo PRs are gated by the `benchmarks` environment;
  # pushes to default branches keep the develop baseline fresh.
  push:
    branches: [ main, develop, master ]
  pull_request:
    branches: [ main, develop, master ]
  workflow_dispatch:

permissions:
  contents: read
  pull-requests: write

concurrency:
  group: bench-accuracy-${{ github.workflow }}-${{ github.ref }}
  cancel-in-progress: true

env:
  BENCHER_PROJECT: nfft-accuracy

jobs:
  accuracy:
    # OPT-IN GATE on pull_request (same protected environment as bench-linux);
    # empty name on push/dispatch -> no gate.
    environment: ${{ github.event_name == 'pull_request' && 'benchmarks' || '' }}
    runs-on: ubuntu-latest
    strategy:
      fail-fast: false
      matrix:
        window: ["kaiserbessel", "gaussian", "bspline", "sinc"]
        precision_opt: ["", "--enable-float", "--enable-long-double"]
    env:
      BUILD_CONFIG: "ubuntu-latest_gcc_${{ matrix.window }}_${{ matrix.precision_opt == '' && 'double' || matrix.precision_opt == '--enable-float' && 'float' || matrix.precision_opt == '--enable-long-double' && 'long-double' || 'unknown' }}"
    steps:
      - uses: actions/checkout@de0fac2e4500dabe0009e67214ff5f5447ce83dd # v6
        with:
          persist-credentials: false
          fetch-depth: 0   # start-point-hash needs base-branch history

      - name: Install dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y build-essential autoconf automake libtool \
            libfftw3-dev libcunit1-dev

      - name: Set up uv
        uses: astral-sh/setup-uv@08807647e7069bb48b6ef5acd8ec9567f424441b # v7
        with:
          enable-cache: true

      - name: Bootstrap
        run: ./bootstrap.sh

      - name: Configure
        env:
          WINDOW: ${{ matrix.window }}
          PRECISION_OPT: ${{ matrix.precision_opt }}
        run: |
          ./configure --with-window="${WINDOW}" ${PRECISION_OPT} \
            --enable-all --enable-tests --enable-exhaustive-unit-tests

      - name: Build library and test program
        run: |
          make -j                  # library
          make -C tests checkall   # check_PROGRAM, not built by `make all`

      - name: Run accuracy suite (serial, emit NDJSON)
        run: |
          rm -f tests/accuracy.ndjson
          NFFT_BENCH_OUT="$PWD/tests/accuracy.ndjson" tests/checkall > /dev/null
          test -s tests/accuracy.ndjson

      - name: Aggregate NDJSON -> BMF JSON
        run: |
          uv run python tests/bench/ndjson_to_bmf.py \
            tests/accuracy.ndjson tests/accuracy.bmf.json
          head -c 400 tests/accuracy.bmf.json; echo

      - name: Set up Bencher CLI
        uses: bencherdev/bencher@__BENCHER_SHA__ # pin via Step 1

      - name: Upload to Bencher (track-only)
        env:
          BENCHER_API_TOKEN: ${{ secrets.BENCHER_API_TOKEN }}
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        run: |
          # Track-only: no --err and no thresholds, so nothing ever fails CI.
          # Fork PRs have an empty token -> --dry-run validates the pipeline.
          ARGS=()
          if [ -z "${BENCHER_API_TOKEN}" ]; then ARGS+=(--dry-run); fi
          if [ "${GITHUB_EVENT_NAME}" = "pull_request" ]; then
            # Standard start-point recipe: clone+reset the base branch's baseline
            # so the PR is compared against the line it targets (usually develop).
            ARGS+=(--branch "${GITHUB_HEAD_REF}"
                   --start-point "${GITHUB_BASE_REF}"
                   --start-point-hash "$(git rev-parse "origin/${GITHUB_BASE_REF}")"
                   --start-point-reset)
          else
            ARGS+=(--branch "${GITHUB_REF_NAME}")   # push/dispatch: develop/main/...
          fi
          bencher run \
            "${ARGS[@]}" \
            --project "${BENCHER_PROJECT}" \
            --token "${BENCHER_API_TOKEN}" \
            --testbed "${BUILD_CONFIG}" \
            --adapter json \
            --file tests/accuracy.bmf.json \
            --github-actions "${GH_TOKEN}"
```

> **Quota fallback (Task 1 Step 2):** if quota is too small, first merge file+online in the converter (`group_key` drops `oracle` → 108/cell), then reduce the PR matrix.

> **`bencher run --file` note:** `--file` ingests an existing results file. If the pinned CLI predates `--file`, replace it with a command form: `bencher run "${ARGS[@]}" ... "cat tests/accuracy.bmf.json"`. Verify against `bencher run --help` in Step 4.

- [ ] **Step 3: Lint the workflow YAML**
```bash
uv run --with yamllint yamllint -d relaxed .github/workflows/bench-accuracy-linux.yml
command -v actionlint >/dev/null && actionlint .github/workflows/bench-accuracy-linux.yml || echo "actionlint not installed; skipping"
```
Expected: no errors (warnings acceptable).

- [ ] **Step 4: Dry-run the upload locally (no token needed)**
```bash
bencher run --help | grep -- '--file' || echo "CLI lacks --file; use the cat fallback"
bencher run --dry-run --project nfft-accuracy --branch local-test \
  --testbed ubuntu-latest_gcc_kaiserbessel_double \
  --adapter json --file /tmp/acc.bmf.json
```
Expected: the CLI parses and reports it would upload ~180 benchmarks (dry run), exit 0.

- [ ] **Step 5: Commit and trigger the workflow**
```bash
git add .github/workflows/bench-accuracy-linux.yml
git commit -m "ci: track per-case test accuracy in Bencher (track-only)"
git push
gh workflow run "Bencher Accuracy" --ref "$(git rev-parse --abbrev-ref HEAD)" 2>/dev/null || true
gh run list --workflow "Bencher Accuracy" --limit 1
gh run watch "$(gh run list --workflow 'Bencher Accuracy' --limit 1 --json databaseId --jq '.[0].databaseId')"
```
Expected: job succeeds; "Upload to Bencher" logs an upload (or `--dry-run` on a fork). Confirm the `nfft-accuracy` project shows testbeds like `ubuntu-latest_gcc_kaiserbessel_double` with `tightness-ratio` + `max-error` measures on metrics named `nfft/online/fast/forward/3d/...`.

---

### Task 5: Ignore artifacts and document the pipeline

**Files:**
- Create: `tests/.gitignore`
- Create: `docs/agents/accuracy-tracking.md`
- Modify: `CLAUDE.md` (pointer under "Running the tests")

**Interfaces:**
- Consumes: the full pipeline (Tasks 2–4). Produces: docs + ignored local artifacts.

- [ ] **Step 1: Ignore generated artifacts**

Create `tests/.gitignore`:
```gitignore
# Local Bencher accuracy artifacts (see docs/agents/accuracy-tracking.md).
accuracy.ndjson
accuracy.bmf.json
```

- [ ] **Step 2: Write the pipeline doc**

Create `docs/agents/accuracy-tracking.md`:
````markdown
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
- **Scope.** Only nfft/nfct/nfst emit today. To add a module, give its harness the
  same `#include "bench_emit.h"` + one `bench_emit_accuracy(...)` call.
- **Quota dial.** Merge file+online in `group_key` (180 → 108 metrics/cell) first.
````

- [ ] **Step 3: Add a pointer from CLAUDE.md**

In `CLAUDE.md`, under "3. Running the tests" (after the reference-data paragraph), add:
```markdown
Test **accuracy figures** (each case's `err` vs `bound`) are additionally tracked
over time in Bencher; set `NFFT_BENCH_OUT=<file>` when running `tests/checkall` to
emit them. See [`docs/agents/accuracy-tracking.md`](docs/agents/accuracy-tracking.md)
and [ADR-0004](docs/adr/0004-accuracy-tracking-with-bencher.md).
```

- [ ] **Step 4: Verify docs reference real paths**
```bash
grep -q "bench_emit.c" docs/agents/accuracy-tracking.md && \
test -f tests/bench_emit.c && test -f tests/bench/ndjson_to_bmf.py && \
test -f .github/workflows/bench-accuracy-linux.yml && \
test -f docs/adr/0004-accuracy-tracking-with-bencher.md && echo "all referenced paths exist"
```
Expected: `all referenced paths exist`.

- [ ] **Step 5: Commit**
```bash
git add tests/.gitignore docs/agents/accuracy-tracking.md CLAUDE.md
git commit -m "docs: document Bencher accuracy-tracking pipeline"
```

---

## Self-Review

**Spec coverage:**
- "Capture accuracy figures in Bencher" → Tasks 2 (emit), 3 (aggregate), 4 (upload). ✓
- "Many/all tests return an accuracy figure vs a threshold" → emitter captures `err`+`bound` at the exact `check_single` site for nfft/nfct/nfst; other modules scoped out. ✓
- "Bencher for tests, CodSpeed keeps benchmarks" → separate workflow + project `nfft-accuracy`, `json` adapter, no CodSpeed overlap; vocabulary kept distinct in CONTEXT.md. ✓
- "Integrate into CI" → Task 4 mirrors `bench-linux.yml` security/gating. ✓
- Grilling decisions captured: aggregate by error-shaping params (max over N/M), ratio + raw measures, file/online split, 12-cell push + gated PR, develop baseline + start-point, track-only/phased, aggregation in converter. ✓ (ADR-0004)

**Placeholder scan:** Only intentional, instruction-backed tokens remain: `__BENCHER_SHA__` (resolved in Task 4 Step 1) and `{owner}/{repo}` in a `gh api` call. No "TBD"/"add error handling"/"write tests for the above". ✓

**Type/name consistency:** `bench_emit_accuracy` signature identical in `.h`, `.c`, and all three call sites (now `oracle` param, structured emission). NDJSON keys (`module/oracle/dim/N/M/init/trafo/accuracy/bound/ok`) match between the C `fprintf` and the Python `convert`/tests. `convert`/`read_ndjson`/`group_key`/`metric_name`/`slug` names match converter↔tests. Measure slugs `tightness-ratio`/`max-error` consistent across converter, docs, and Task 4 verification. `BUILD_CONFIG`/testbed expression mirrors `build-linux.yml`. ✓
