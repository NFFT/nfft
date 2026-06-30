# Accuracy Reporting Layer (P1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **Revision (2026-06-30) — implemented, with three post-build refinements** (authoritative: the spec + the actual code):
> 1. **No-baseline path.** The first PR (or any PR before `develop` has published a baseline) now shows absolute accuracy with a clear "no baseline yet" message — not a misleading "unchanged". New `report.comment_body_no_baseline` / `check_summary_no_baseline` and `pr_report.py --no-baseline`.
> 2. **Permanent PR heatmap links.** PR images are linked via `raw.githubusercontent.com/<owner>/<repo>/gh-pages/pr/<n>/…` (stable, inline-rendering, work before Pages is enabled) rather than the Pages site URL.
> 3. **gh-pages branch confirmed** as the archive (permanent, accumulating per-PR; `actions/deploy-pages` replaces the site each deploy and artifacts expire). Fork PRs cannot write `gh-pages` (read-only token) → emoji/text only.

**Goal:** Replace Bencher's noisy PR comments and unscalable dashboard with a GitHub Pages heatmap dashboard of absolute accuracy plus a quiet, bidirectional, always-present PR signal (Check + diff comment), keeping Bencher purely as the data archive.

**Architecture:** New pure-Python modules under `tests/bench/` (`diff.py`, `report.py`, `heatmap.py`) plus two thin CLIs (`dashboard.py`, `pr_report.py`) transform the existing per-testbed BMF JSON into heatmap PNGs, an emoji grid, a PR comment body, and a Check summary. A new `accuracy-report` job in `build-linux.yml` runs them: on `develop` push it publishes the absolute dashboard + baseline BMFs to a `gh-pages` branch; on a PR it fetches the baseline, diffs, and posts the Check + comment. The diff is a pure function of two sets of BMF files — no Bencher API in the PR path.

**Tech Stack:** Python 3.13 (stdlib + `matplotlib`, run via `uv run --with matplotlib`), pytest, GitHub Actions, bash.

## Global Constraints

- **Spec:** [`docs/superpowers/specs/2026-06-30-accuracy-reporting-layer-design.md`](../specs/2026-06-30-accuracy-reporting-layer-design.md). P2 (sweep tests + scaling curves) is OUT OF SCOPE.
- **BMF input shape (existing):** one file per testbed, `{ "<metric>": {"accuracy-digits": {"value": <float>}, "max-error": {"value": <float>}}, ... }`. Metric name = `<module>/<runtime>/<oracle>/<speed>/<direction>/<dim>d/<init-slug>`. The reporting layer reads only `accuracy-digits`.
- **Significance gate:** a metric is *changed* when `|Δ accuracy-digits| ≥ 0.5` (default; configurable). Improvements and regressions are symmetric.
- **Track-only:** the GitHub Check NEVER fails CI (`conclusion` is `neutral`). Nothing here may fail a build.
- **Comment is always posted and upserted** (one comment, updated in place), compact when nothing changed.
- **Pure modules are import-by-module-name** and tested from `tests/bench/` (matching the existing `ndjson_to_bmf.py` convention): run `cd tests/bench && uv run --with pytest[,matplotlib] python -m pytest`.
- **Pin every GitHub Action to a commit SHA** (repo convention); reuse SHAs already in `build-linux.yml`.
- **Fork-gating mirrors the existing `bencher-upload` job:** fork PRs gated behind the reviewer environment; everything else unattended. Fork PRs get no `gh-pages` write → comment omits PNG links and uses the emoji grid only.
- **build-linux concurrency** is `group: build-linux, cancel-in-progress: false` → report-job `gh-pages` writes are already serialized (no concurrent-push races).

## Data model (used across Tasks 1–4)

```python
# tests/bench/diff.py
from dataclasses import dataclass, field

@dataclass
class Change:
    testbed: str
    name: str            # metric name within the testbed
    base_digits: float
    pr_digits: float
    delta_digits: float  # pr - base (positive = improvement)
    pct: float           # 100 * delta_digits / base_digits

@dataclass
class DiffResult:
    improvements: list           # list[Change], |Δ|>=gate, delta>0, sorted by |Δ| desc
    regressions: list            # list[Change], |Δ|>=gate, delta<0, sorted by |Δ| desc
    unchanged_count: int
    added: list = field(default_factory=list)    # list[str] "testbed::name" only in PR
    removed: list = field(default_factory=list)  # list[str] "testbed::name" only in baseline
    by_testbed: dict = field(default_factory=dict)  # testbed -> {"improved","regressed","unchanged"}
```

---

### Task 1: `diff.py` — pure BMF diff

**Files:**
- Create: `tests/bench/diff.py`
- Test: `tests/bench/test_diff.py`

**Interfaces:**
- Produces: `Change`, `DiffResult` (above); `digits(bmf, name) -> float`; `load_bmf_tree(dir) -> dict[str, dict]` (testbed → BMF, testbed = filename minus `.bmf.json`); `diff(pr_tree, base_tree, gate=0.5) -> DiffResult`.
- Consumes: nothing (first task).

- [ ] **Step 1: Write the failing tests**

Create `tests/bench/test_diff.py`:
```python
import json

import pytest

from diff import Change, DiffResult, diff, digits, load_bmf_tree


def _bmf(**digit_by_name):
    return {n: {"accuracy-digits": {"value": v}, "max-error": {"value": 10 ** -v}}
            for n, v in digit_by_name.items()}


def test_digits_reads_accuracy_digits():
    assert digits(_bmf(a=13.5), "a") == 13.5


def test_improvement_and_regression_classified_and_pct():
    pr = {"t1": _bmf(a=14.3, b=12.0)}     # a improved +0.8, b regressed -1.0
    base = {"t1": _bmf(a=13.5, b=13.0)}
    r = diff(pr, base, gate=0.5)
    assert len(r.improvements) == 1 and len(r.regressions) == 1
    imp = r.improvements[0]
    assert (imp.testbed, imp.name) == ("t1", "a")
    assert imp.delta_digits == pytest.approx(0.8)
    assert imp.pct == pytest.approx(100 * 0.8 / 13.5)
    assert r.regressions[0].name == "b"
    assert r.unchanged_count == 0


def test_below_gate_is_unchanged():
    pr = {"t1": _bmf(a=13.9)}
    base = {"t1": _bmf(a=13.5)}            # +0.4 < 0.5 gate
    r = diff(pr, base, gate=0.5)
    assert r.improvements == [] and r.regressions == []
    assert r.unchanged_count == 1
    assert r.by_testbed["t1"] == {"improved": 0, "regressed": 0, "unchanged": 1}


def test_groups_sorted_by_abs_delta_desc():
    pr = {"t1": _bmf(a=15.0, b=16.0)}     # a +1.0, b +3.0
    base = {"t1": _bmf(a=14.0, b=13.0)}
    r = diff(pr, base, gate=0.5)
    assert [c.name for c in r.improvements] == ["b", "a"]


def test_added_and_removed_tracked_not_counted():
    pr = {"t1": _bmf(a=13.5, c=13.5)}
    base = {"t1": _bmf(a=13.5, d=13.5)}
    r = diff(pr, base, gate=0.5)
    assert r.added == ["t1::c"] and r.removed == ["t1::d"]
    assert r.improvements == [] and r.regressions == []


def test_load_bmf_tree_keys_by_testbed(tmp_path):
    (tmp_path / "ubuntu_gcc_kb_double.bmf.json").write_text(json.dumps(_bmf(a=13.5)))
    tree = load_bmf_tree(str(tmp_path))
    assert set(tree) == {"ubuntu_gcc_kb_double"}
    assert tree["ubuntu_gcc_kb_double"]["a"]["accuracy-digits"]["value"] == 13.5
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd tests/bench && uv run --with pytest python -m pytest test_diff.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'diff'`.

- [ ] **Step 3: Write `diff.py`**

Create `tests/bench/diff.py`:
```python
"""Pure diff of two sets of per-testbed BMF files on the accuracy-digits measure."""
from __future__ import annotations

import glob
import json
import os
from dataclasses import dataclass, field


@dataclass
class Change:
    testbed: str
    name: str
    base_digits: float
    pr_digits: float
    delta_digits: float
    pct: float


@dataclass
class DiffResult:
    improvements: list
    regressions: list
    unchanged_count: int
    added: list = field(default_factory=list)
    removed: list = field(default_factory=list)
    by_testbed: dict = field(default_factory=dict)


def digits(bmf, name):
    return float(bmf[name]["accuracy-digits"]["value"])


def load_bmf_tree(directory):
    tree = {}
    for path in sorted(glob.glob(os.path.join(directory, "*.bmf.json"))):
        testbed = os.path.basename(path)[: -len(".bmf.json")]
        with open(path, encoding="utf-8") as f:
            tree[testbed] = json.load(f)
    return tree


def diff(pr_tree, base_tree, gate=0.5):
    improvements, regressions, added, removed = [], [], [], []
    unchanged = 0
    by_testbed = {}
    for testbed in sorted(set(pr_tree) | set(base_tree)):
        pr_bmf = pr_tree.get(testbed, {})
        base_bmf = base_tree.get(testbed, {})
        counts = {"improved": 0, "regressed": 0, "unchanged": 0}
        for name in sorted(set(pr_bmf) | set(base_bmf)):
            in_pr, in_base = name in pr_bmf, name in base_bmf
            if in_pr and not in_base:
                added.append(f"{testbed}::{name}")
                continue
            if in_base and not in_pr:
                removed.append(f"{testbed}::{name}")
                continue
            b, p = digits(base_bmf, name), digits(pr_bmf, name)
            delta = p - b
            pct = (100.0 * delta / b) if b != 0 else 0.0
            if abs(delta) >= gate:
                change = Change(testbed, name, b, p, delta, pct)
                (improvements if delta > 0 else regressions).append(change)
                counts["improved" if delta > 0 else "regressed"] += 1
            else:
                unchanged += 1
                counts["unchanged"] += 1
        by_testbed[testbed] = counts
    improvements.sort(key=lambda c: abs(c.delta_digits), reverse=True)
    regressions.sort(key=lambda c: abs(c.delta_digits), reverse=True)
    return DiffResult(improvements, regressions, unchanged, added, removed, by_testbed)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd tests/bench && uv run --with pytest python -m pytest test_diff.py -q`
Expected: PASS (6 passed).

- [ ] **Step 5: Commit**
```bash
git add tests/bench/diff.py tests/bench/test_diff.py
git commit -m "test: add pure BMF accuracy-digits diff (improvements/regressions)"
```

---

### Task 2: `report.py` — Check summary + comment markdown

**Files:**
- Create: `tests/bench/report.py`
- Test: `tests/bench/test_report.py`

**Interfaces:**
- Consumes: `DiffResult`, `Change` from `diff.py`.
- Produces: `check_summary(result) -> (conclusion: str, title: str, summary_md: str)` (`conclusion` always `"neutral"`); `comment_body(result, png_urls: dict | None) -> str`. `png_urls` is `{"absolute": url, "relative": url}` or `None`/`{}` for fork PRs. Marker line `<!-- nfft-accuracy-report -->` is the FIRST line of every comment (used by the workflow to find-and-update).

- [ ] **Step 1: Write the failing tests**

Create `tests/bench/test_report.py`:
```python
from diff import Change, DiffResult
from report import MARKER, check_summary, comment_body


def _result(impr=0, regr=0, unchanged=5):
    imps = [Change("t1", f"i{k}", 13.0, 13.0 + 1 + k, 1 + k, 100 * (1 + k) / 13.0)
            for k in range(impr)]
    regs = [Change("t1", f"r{k}", 14.0, 14.0 - 1 - k, -(1 + k), -100 * (1 + k) / 14.0)
            for k in range(regr)]
    return DiffResult(imps, regs, unchanged, by_testbed={})


def test_check_never_fails():
    conclusion, _, _ = check_summary(_result())
    assert conclusion == "neutral"


def test_check_title_flat_vs_changed():
    assert "unchanged" in check_summary(_result())[1].lower()
    assert check_summary(_result(impr=2, regr=1))[1] == "1 regressed, 2 improved"


def test_comment_has_marker_first_line():
    assert comment_body(_result(), None).splitlines()[0] == MARKER


def test_flat_comment_is_one_line_summary():
    body = comment_body(_result(unchanged=42), None)
    assert "42 cases unchanged" in body
    assert "Improvements" not in body and "Regressions" not in body


def test_changed_comment_lists_groups_and_links():
    body = comment_body(_result(impr=1, regr=1),
                        {"absolute": "http://x/abs.png", "relative": "http://x/rel.png"})
    assert "1 unchanged" not in body  # uses the aggregate line format below
    assert "improved" in body and "regressed" in body
    assert "Improvements" in body and "Regressions" in body
    assert "13.00 → 14.00 digits" in body  # i0: base 13 -> pr 14
    assert "http://x/abs.png" in body and "http://x/rel.png" in body


def test_group_capped_at_10_with_more_note():
    body = comment_body(_result(impr=13), None)
    assert body.count("i") >= 10
    assert "+3 more" in body
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd tests/bench && uv run --with pytest python -m pytest test_report.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'report'`.

- [ ] **Step 3: Write `report.py`**

Create `tests/bench/report.py`:
```python
"""Render the PR Check summary and the upserted comment body from a DiffResult."""
from __future__ import annotations

MARKER = "<!-- nfft-accuracy-report -->"
CAP = 10


def check_summary(result):
    y, z = len(result.improvements), len(result.regressions)
    if y == 0 and z == 0:
        title = "accuracy unchanged"
    else:
        title = f"{z} regressed, {y} improved"
    summary = (f"{result.unchanged_count} unchanged · {y} improved · "
               f"{z} regressed")
    return "neutral", title, summary


def _line(c):
    return (f"- `{c.testbed} {c.name}`: {c.base_digits:.2f} → "
            f"{c.pr_digits:.2f} digits ({c.pct:+.0f}%)")


def _group(heading, changes):
    out = [f"### {heading}"]
    for c in changes[:CAP]:
        out.append(_line(c))
    if len(changes) > CAP:
        out.append(f"- … +{len(changes) - CAP} more")
    return "\n".join(out)


def comment_body(result, png_urls):
    y, z = len(result.improvements), len(result.regressions)
    if y == 0 and z == 0:
        return f"{MARKER}\nAccuracy: {result.unchanged_count} cases unchanged."
    lines = [MARKER, "## Accuracy report",
             f"{result.unchanged_count} unchanged · {y} improved · "
             f"{z} regressed", ""]
    if png_urls:
        lines.append(f"[absolute heatmap]({png_urls['absolute']}) · "
                     f"[relative heatmap]({png_urls['relative']})")
        lines.append("")
    if result.improvements:
        lines.append(_group("Improvements", result.improvements))
        lines.append("")
    if result.regressions:
        lines.append(_group("Regressions", result.regressions))
    return "\n".join(lines).rstrip() + "\n"
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd tests/bench && uv run --with pytest python -m pytest test_report.py -q`
Expected: PASS (6 passed).

- [ ] **Step 5: Commit**
```bash
git add tests/bench/report.py tests/bench/test_report.py
git commit -m "test: render accuracy PR Check summary and upserted comment body"
```

---

### Task 3: `heatmap.py` — PNG renderer + emoji grid

**Files:**
- Create: `tests/bench/heatmap.py`
- Test: `tests/bench/test_heatmap.py`

**Interfaces:**
- Consumes: a BMF tree (`dict[testbed, bmf]`) for absolute; a `DiffResult` for relative/emoji.
- Produces: `emoji_grid(result) -> str` (markdown table; 💚/🟩 improve by magnitude, 🟥/🟧 regress by magnitude, `·` unchanged); `render_absolute(tree, out_png) -> None`; `render_relative(result, out_png) -> None`. `BIG = 1.0` (digit magnitude splitting the two emoji tiers).

- [ ] **Step 1: Write the failing tests**

Create `tests/bench/test_heatmap.py`:
```python
import os

from diff import Change, DiffResult
from heatmap import emoji_grid, render_absolute, render_relative


def _result():
    return DiffResult(
        improvements=[Change("t1", "m/serial/file/fast/forward/1d/a", 13, 14.5, 1.5, 11)],
        regressions=[Change("t1", "m/serial/file/fast/forward/2d/a", 14, 13.7, -0.3, -2)]
        and [Change("t1", "m/serial/file/fast/forward/2d/a", 14, 12.5, -1.5, -11)],
        unchanged_count=3,
        by_testbed={"t1": {"improved": 1, "regressed": 1, "unchanged": 3}},
    )


def test_emoji_grid_uses_magnitude_tiers():
    grid = emoji_grid(_result())
    assert "\U0001F49A" in grid  # 💚 large improvement (>=1.0)
    assert "\U0001F7E5" in grid  # 🟥 large regression (>=1.0)
    assert "|" in grid           # markdown table


def test_render_absolute_writes_png(tmp_path):
    tree = {"t1": {"m/serial/file/fast/forward/1d/a":
                   {"accuracy-digits": {"value": 13.5}, "max-error": {"value": 1e-13}}}}
    out = str(tmp_path / "abs.png")
    render_absolute(tree, out)
    assert os.path.getsize(out) > 0


def test_render_relative_writes_png(tmp_path):
    out = str(tmp_path / "rel.png")
    render_relative(_result(), out)
    assert os.path.getsize(out) > 0
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd tests/bench && uv run --with pytest,matplotlib python -m pytest test_heatmap.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'heatmap'`.

- [ ] **Step 3: Write `heatmap.py`**

Create `tests/bench/heatmap.py`:
```python
"""Render absolute / relative accuracy heatmaps (PNG) and an inline emoji grid.

matplotlib uses the non-interactive Agg backend so it runs headless in CI.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BIG = 1.0


def _emoji(delta, gate=0.5):
    if delta >= BIG:
        return "\U0001F49A"   # 💚
    if delta >= gate:
        return "\U0001F7E9"   # 🟩
    if delta <= -BIG:
        return "\U0001F7E5"   # 🟥
    if delta <= -gate:
        return "\U0001F7E7"   # 🟧
    return "·"           # ·


def emoji_grid(result):
    rows = {c.name: c for c in result.improvements + result.regressions}
    if not rows:
        return "_No significant changes._"
    testbeds = sorted(result.by_testbed)
    out = ["| case | " + " | ".join(testbeds) + " |",
           "|" + "---|" * (len(testbeds) + 1)]
    by_name = {}
    for c in result.improvements + result.regressions:
        by_name.setdefault(c.name, {})[c.testbed] = c.delta_digits
    for name in sorted(by_name):
        cells = [_emoji(by_name[name].get(t, 0.0)) for t in testbeds]
        out.append(f"| `{name}` | " + " | ".join(cells) + " |")
    return "\n".join(out)


def _heatmap(matrix, rows, cols, out_png, cmap, vcenter=None):
    fig, ax = plt.subplots(figsize=(max(6, len(cols) * 0.6),
                                    max(4, len(rows) * 0.25)))
    im = ax.imshow(matrix, aspect="auto", cmap=cmap)
    ax.set_xticks(range(len(cols)), cols, rotation=90, fontsize=6)
    ax.set_yticks(range(len(rows)), rows, fontsize=5)
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def render_absolute(tree, out_png):
    testbeds = sorted(tree)
    names = sorted({n for bmf in tree.values() for n in bmf})
    matrix = [[tree[t][n]["accuracy-digits"]["value"] if n in tree[t] else float("nan")
               for t in testbeds] for n in names]
    _heatmap(matrix, names, testbeds, out_png, cmap="RdYlGn")


def render_relative(result, out_png):
    changes = result.improvements + result.regressions
    testbeds = sorted(result.by_testbed)
    names = sorted({c.name for c in changes})
    delta = {(c.name, c.testbed): c.delta_digits for c in changes}
    matrix = [[delta.get((n, t), 0.0) for t in testbeds] for n in names] or [[0.0]]
    _heatmap(matrix, names or ["(none)"], testbeds or ["(none)"], out_png, cmap="RdYlGn")
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd tests/bench && uv run --with pytest,matplotlib python -m pytest test_heatmap.py -q`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**
```bash
git add tests/bench/heatmap.py tests/bench/test_heatmap.py
git commit -m "test: render accuracy heatmap PNGs and inline emoji grid"
```

---

### Task 4: CLIs — `dashboard.py` (develop) and `pr_report.py` (PR)

**Files:**
- Create: `tests/bench/dashboard.py`
- Create: `tests/bench/pr_report.py`
- Test: `tests/bench/test_clis.py`

**Interfaces:**
- Consumes: `diff.load_bmf_tree`, `diff.diff`, `heatmap.render_absolute/render_relative`, `report.check_summary/comment_body`.
- Produces:
  - `dashboard.py main(argv)`: `dashboard.py <bmf_dir> <out_dir>` → writes `<out_dir>/absolute.png`, `<out_dir>/index.html`, and copies each `<bmf_dir>/<testbed>.bmf.json` to `<out_dir>/baseline/<testbed>.bmf.json`.
  - `pr_report.py main(argv)`: `pr_report.py <pr_bmf_dir> <base_bmf_dir> <out_dir> [--abs-url U] [--rel-url U]` → writes `<out_dir>/absolute.png`, `<out_dir>/relative.png`, `<out_dir>/comment.md`, `<out_dir>/check.json` (`{"conclusion","title","summary"}`). PNG URLs are embedded in `comment.md` only when both `--abs-url` and `--rel-url` are given.

- [ ] **Step 1: Write the failing tests**

Create `tests/bench/test_clis.py`:
```python
import json
import os

import dashboard
import pr_report


def _write_bmf(path, **digit_by_name):
    with open(path, "w", encoding="utf-8") as f:
        json.dump({n: {"accuracy-digits": {"value": v}, "max-error": {"value": 10 ** -v}}
                   for n, v in digit_by_name.items()}, f)


def test_dashboard_emits_png_index_and_baseline(tmp_path):
    src = tmp_path / "bmf"; src.mkdir()
    _write_bmf(src / "tb1.bmf.json", **{"m/serial/file/fast/forward/1d/a": 13.5})
    out = tmp_path / "site"
    dashboard.main([str(src), str(out)])
    assert (out / "absolute.png").stat().st_size > 0
    assert (out / "index.html").exists()
    assert (out / "baseline" / "tb1.bmf.json").exists()


def test_pr_report_changed_writes_all_artifacts(tmp_path):
    pr = tmp_path / "pr"; pr.mkdir()
    base = tmp_path / "base"; base.mkdir()
    _write_bmf(pr / "tb1.bmf.json", **{"m/serial/file/fast/forward/1d/a": 15.0})
    _write_bmf(base / "tb1.bmf.json", **{"m/serial/file/fast/forward/1d/a": 13.0})
    out = tmp_path / "out"
    pr_report.main([str(pr), str(base), str(out),
                    "--abs-url", "http://x/a.png", "--rel-url", "http://x/r.png"])
    check = json.load(open(out / "check.json"))
    assert check["conclusion"] == "neutral" and "improved" in check["title"]
    body = (out / "comment.md").read_text()
    assert "Improvements" in body and "http://x/a.png" in body
    assert (out / "relative.png").stat().st_size > 0


def test_pr_report_flat_has_no_links_and_one_liner(tmp_path):
    pr = tmp_path / "pr"; pr.mkdir()
    base = tmp_path / "base"; base.mkdir()
    _write_bmf(pr / "tb1.bmf.json", **{"m/serial/file/fast/forward/1d/a": 13.5})
    _write_bmf(base / "tb1.bmf.json", **{"m/serial/file/fast/forward/1d/a": 13.5})
    out = tmp_path / "out"
    pr_report.main([str(pr), str(base), str(out)])
    body = (out / "comment.md").read_text()
    assert "cases unchanged" in body and "http" not in body
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd tests/bench && uv run --with pytest,matplotlib python -m pytest test_clis.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'dashboard'`.

- [ ] **Step 3: Write `dashboard.py`**

Create `tests/bench/dashboard.py`:
```python
"""CLI: render the absolute accuracy dashboard + copy baseline BMFs (develop)."""
from __future__ import annotations

import argparse
import os
import shutil
import sys

from diff import load_bmf_tree
from heatmap import render_absolute

INDEX = """<!doctype html><meta charset=utf-8>
<title>NFFT accuracy dashboard</title>
<h1>NFFT accuracy (develop baseline)</h1>
<p>Worst-case accurate digits per case &times; testbed. Green = more digits.</p>
<img src="absolute.png" style="max-width:100%">
"""


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("bmf_dir")
    ap.add_argument("out_dir")
    args = ap.parse_args(argv)
    os.makedirs(os.path.join(args.out_dir, "baseline"), exist_ok=True)
    tree = load_bmf_tree(args.bmf_dir)
    render_absolute(tree, os.path.join(args.out_dir, "absolute.png"))
    with open(os.path.join(args.out_dir, "index.html"), "w", encoding="utf-8") as f:
        f.write(INDEX)
    for testbed in tree:
        shutil.copyfile(os.path.join(args.bmf_dir, f"{testbed}.bmf.json"),
                        os.path.join(args.out_dir, "baseline", f"{testbed}.bmf.json"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Write `pr_report.py`**

Create `tests/bench/pr_report.py`:
```python
"""CLI: diff PR BMFs vs baseline, render heatmaps, write comment.md + check.json."""
from __future__ import annotations

import argparse
import json
import os
import sys

from diff import diff, load_bmf_tree
from heatmap import emoji_grid, render_absolute, render_relative
from report import check_summary, comment_body


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("pr_bmf_dir")
    ap.add_argument("base_bmf_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--abs-url")
    ap.add_argument("--rel-url")
    ap.add_argument("--gate", type=float, default=0.5)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    pr_tree = load_bmf_tree(args.pr_bmf_dir)
    base_tree = load_bmf_tree(args.base_bmf_dir)
    result = diff(pr_tree, base_tree, gate=args.gate)

    render_absolute(pr_tree, os.path.join(args.out_dir, "absolute.png"))
    render_relative(result, os.path.join(args.out_dir, "relative.png"))

    png_urls = None
    if args.abs_url and args.rel_url:
        png_urls = {"absolute": args.abs_url, "relative": args.rel_url}

    # The emoji grid is embedded in the comment ABOVE the itemized groups.
    body = comment_body(result, png_urls)
    if result.improvements or result.regressions:
        body = body.replace("## Accuracy report\n",
                            "## Accuracy report\n\n" + emoji_grid(result) + "\n\n", 1)
    with open(os.path.join(args.out_dir, "comment.md"), "w", encoding="utf-8") as f:
        f.write(body)

    conclusion, title, summary = check_summary(result)
    with open(os.path.join(args.out_dir, "check.json"), "w", encoding="utf-8") as f:
        json.dump({"conclusion": conclusion, "title": title, "summary": summary}, f)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `cd tests/bench && uv run --with pytest,matplotlib python -m pytest -q`
Expected: PASS (all of test_diff, test_report, test_heatmap, test_clis, and the existing test_ndjson_to_bmf).

- [ ] **Step 6: Commit**
```bash
git add tests/bench/dashboard.py tests/bench/pr_report.py tests/bench/test_clis.py
git commit -m "test: dashboard and PR-report CLIs over BMF directories"
```

---

### Task 5: CI wiring — `accuracy-report` job + suppress Bencher comment

**Files:**
- Modify: `.github/workflows/build-linux.yml` (add `accuracy-report` job; remove `--github-actions` from the `bencher-upload` job)

**Interfaces:**
- Consumes: the `accuracy-bmf-*` artifacts (existing), `tests/bench/dashboard.py` + `pr_report.py` (Task 4).
- Produces: a `gh-pages` branch (dashboard + `baseline/*.bmf.json` + `pr/<n>/*.png`), a GitHub Check, and an upserted PR comment.

- [ ] **Step 1: Remove Bencher's PR comment**

In `.github/workflows/build-linux.yml`, in the `bencher-upload` job's `bencher run` invocation, delete the line:
```yaml
              --github-actions "${GH_TOKEN}"
```
(and the trailing `\` on the now-last argument line `--file "${file}"`). Bencher keeps uploading; it no longer comments.

- [ ] **Step 2: Add the `accuracy-report` job**

In `.github/workflows/build-linux.yml`, add this job at the `jobs:` level (e.g. after `bencher-upload`, before `distcheck`). Replace `__CHECKOUT_SHA__` with the same checkout SHA used elsewhere in the file (`de0fac2e4500dabe0009e67214ff5f5447ce83dd`), `__UV_SHA__` with the setup-uv SHA (`08807647e7069bb48b6ef5acd8ec9567f424441b`), and `__DLA_SHA__` with the download-artifact SHA already used by `bencher-upload` (`37930b1c2abaa49bbe596cd826c3c89aef350131`):

```yaml
  accuracy-report:
    needs: build
    runs-on: ubuntu-latest
    # Fork PRs gated (no gh-pages write); everything else unattended. Mirrors
    # bencher-upload. Never fails CI (the Check is informational).
    environment: ${{ github.event_name == 'pull_request' && github.event.pull_request.head.repo.fork && 'benchmarks' || 'bencher-baseline' }}
    permissions:
      contents: write          # push to gh-pages
      checks: write            # post the Check
      pull-requests: write      # upsert the comment
    steps:
      - uses: actions/checkout@de0fac2e4500dabe0009e67214ff5f5447ce83dd # v6
        with:
          persist-credentials: true   # need to push gh-pages
          fetch-depth: 0

      - name: Set up uv
        uses: astral-sh/setup-uv@08807647e7069bb48b6ef5acd8ec9567f424441b # v7
        with:
          enable-cache: true

      - name: Download accuracy BMF artifacts
        uses: actions/download-artifact@37930b1c2abaa49bbe596cd826c3c89aef350131 # v7.0.0
        with:
          pattern: accuracy-bmf-*
          path: bmf-artifacts

      - name: Arrange BMFs as <testbed>.bmf.json
        shell: bash
        run: |
          mkdir -p pr-bmf
          shopt -s nullglob
          for d in bmf-artifacts/accuracy-bmf-*/; do
            tb="$(basename "$d")"; tb="${tb#accuracy-bmf-}"
            [ -f "${d}accuracy.bmf.json" ] && cp "${d}accuracy.bmf.json" "pr-bmf/${tb}.bmf.json"
          done
          ls pr-bmf

      - name: Develop — publish dashboard + baseline to gh-pages
        if: github.event_name == 'push'
        shell: bash
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        run: |
          uv run --with matplotlib python tests/bench/dashboard.py pr-bmf site
          bash .github/scripts/gh-pages-publish.sh site "develop dashboard ${GITHUB_SHA::7}"

      - name: PR — fetch baseline, diff, post Check + comment
        if: github.event_name == 'pull_request'
        shell: bash
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        run: bash .github/scripts/accuracy-pr-report.sh
```

- [ ] **Step 3: Add the gh-pages publish helper**

Create `.github/scripts/gh-pages-publish.sh` (publishes a directory's contents into `gh-pages`, keeping existing files; serialized by the build-linux concurrency group):
```bash
#!/usr/bin/env bash
# Usage: gh-pages-publish.sh <src_dir> <commit_msg> [<dest_subdir>]
set -euo pipefail
src="$1"; msg="$2"; dest="${3:-}"
repo_url="https://x-access-token:${GH_TOKEN}@github.com/${GITHUB_REPOSITORY}.git"
tmp="$(mktemp -d)"
git clone --depth 1 --branch gh-pages "$repo_url" "$tmp" 2>/dev/null || {
  git clone --depth 1 "$repo_url" "$tmp"; git -C "$tmp" checkout --orphan gh-pages
  git -C "$tmp" rm -rf . >/dev/null 2>&1 || true
}
mkdir -p "$tmp/$dest"
cp -r "$src"/. "$tmp/$dest/"
git -C "$tmp" add -A
git -C "$tmp" -c user.name=ci -c user.email=ci@nfft commit -m "$msg" || { echo "no changes"; exit 0; }
git -C "$tmp" push origin gh-pages
```

- [ ] **Step 4: Add the PR report helper**

Create `.github/scripts/accuracy-pr-report.sh`:
```bash
#!/usr/bin/env bash
set -euo pipefail
owner="${GITHUB_REPOSITORY%%/*}"; repo="${GITHUB_REPOSITORY##*/}"
pages_raw="https://raw.githubusercontent.com/${owner}/${repo}/gh-pages/baseline"
pr="${GITHUB_REF_NAME%%/*}"   # refs/pull/<n>/merge -> GITHUB_REF gives the number via event
pr="$(jq -r .number "$GITHUB_EVENT_PATH")"

# Fetch baseline BMFs (skip cleanly if the dashboard hasn't published yet).
mkdir -p base-bmf
for f in pr-bmf/*.bmf.json; do
  tb="$(basename "$f")"
  curl -fsSL "${pages_raw}/${tb}" -o "base-bmf/${tb}" || echo "no baseline yet: ${tb}"
done

is_fork='${{ github.event.pull_request.head.repo.fork }}'
if [ -d base-bmf ] && ls base-bmf/*.bmf.json >/dev/null 2>&1; then
  if [ "$is_fork" = "true" ]; then
    uv run --with matplotlib python tests/bench/pr_report.py pr-bmf base-bmf out
  else
    abs="https://${owner}.github.io/${repo}/pr/${pr}/absolute.png"
    rel="https://${owner}.github.io/${repo}/pr/${pr}/relative.png"
    uv run --with matplotlib python tests/bench/pr_report.py pr-bmf base-bmf out \
      --abs-url "$abs" --rel-url "$rel"
    mkdir -p site/pr/${pr}; cp out/absolute.png out/relative.png site/pr/${pr}/
    bash .github/scripts/gh-pages-publish.sh site/pr/${pr} "pr ${pr} heatmaps" "pr/${pr}"
  fi
else
  # No baseline published yet: still render PR-only heatmaps + a flat-ish comment.
  uv run --with matplotlib python tests/bench/pr_report.py pr-bmf pr-bmf out
fi

# Post the Check (never fails CI).
title="$(jq -r .title out/check.json)"; summary="$(jq -r .summary out/check.json)"
gh api -X POST "repos/${GITHUB_REPOSITORY}/check-runs" \
  -f name="accuracy" -f head_sha="${GITHUB_SHA}" -f status="completed" \
  -f conclusion="neutral" -f output[title]="$title" -f output[summary]="$summary"

# Upsert the comment (find by marker, else create).
body="$(cat out/comment.md)"
existing="$(gh api "repos/${GITHUB_REPOSITORY}/issues/${pr}/comments" \
  --jq '.[] | select(.body | startswith("<!-- nfft-accuracy-report -->")) | .id' | head -1)"
if [ -n "$existing" ]; then
  gh api -X PATCH "repos/${GITHUB_REPOSITORY}/issues/comments/${existing}" -f body="$body"
else
  gh api -X POST "repos/${GITHUB_REPOSITORY}/issues/${pr}/comments" -f body="$body"
fi
```

- [ ] **Step 5: Make helpers executable and validate**

Run:
```bash
chmod +x .github/scripts/gh-pages-publish.sh .github/scripts/accuracy-pr-report.sh
bash -n .github/scripts/gh-pages-publish.sh .github/scripts/accuracy-pr-report.sh && echo "bash syntax OK"
uv run --with yamllint yamllint -d relaxed .github/workflows/build-linux.yml | grep -v "line too long\|trailing spaces" || echo "no new yamllint errors"
uv run --with pyyaml python -c "import yaml; d=yaml.safe_load(open('.github/workflows/build-linux.yml')); print('accuracy-report' in d['jobs'])"
```
Expected: `bash syntax OK`; no new yamllint errors beyond the pre-existing ones; prints `True`.

- [ ] **Step 6: Local dry-run of the PR pipeline (no posting)**

Run (uses the converter to make real BMFs, perturbs one to force a change):
```bash
cd /workspaces/nfft2
NFFT_BENCH_OUT="$PWD/tests/accuracy.ndjson" tests/checkall > /dev/null
uv run python tests/bench/ndjson_to_bmf.py tests/accuracy.ndjson /tmp/dry/pr/ubuntu-latest_gcc_kaiserbessel_double.bmf.json 2>/dev/null || \
  { mkdir -p /tmp/dry/pr /tmp/dry/base; uv run python tests/bench/ndjson_to_bmf.py tests/accuracy.ndjson /tmp/dry/pr/ubuntu-latest_gcc_kaiserbessel_double.bmf.json; }
cp -r /tmp/dry/pr/. /tmp/dry/base/
python3 -c "import json,glob;f=glob.glob('/tmp/dry/base/*.json')[0];d=json.load(open(f));k=sorted(d)[0];d[k]['accuracy-digits']['value']-=2.0;json.dump(d,open(f,'w'))"
uv run --with matplotlib python tests/bench/pr_report.py /tmp/dry/pr /tmp/dry/base /tmp/dry/out --abs-url http://x/a.png --rel-url http://x/r.png
echo "--- comment.md ---"; cat /tmp/dry/out/comment.md; echo "--- check.json ---"; cat /tmp/dry/out/check.json
```
Expected: a comment with the marker, an emoji grid, an Improvements group (the perturbed case), the PNG links, and `check.json` with `"conclusion":"neutral"`.

- [ ] **Step 7: Commit**
```bash
git add .github/workflows/build-linux.yml .github/scripts/gh-pages-publish.sh .github/scripts/accuracy-pr-report.sh
git commit -m "ci: add accuracy-report job (Pages dashboard + PR Check/comment), drop Bencher comment"
```

---

### Task 6: Documentation

**Files:**
- Modify: `docs/agents/accuracy-tracking.md`
- Modify: `tests/.gitignore` (ignore local report artifacts)

**Interfaces:** Consumes the whole pipeline. Produces docs.

- [ ] **Step 1: Ignore local report artifacts**

Append to `tests/.gitignore`:
```gitignore
# Local accuracy-report scratch
/bench/out/
/bench/site/
```

- [ ] **Step 2: Document the reporting layer**

In `docs/agents/accuracy-tracking.md`, add a `## Reporting` section after the Pipeline section:
```markdown
## Reporting

Bencher is now only the long-term archive. Human-facing reporting is the
`accuracy-report` job in `build-linux.yml`, built from pure modules in
`tests/bench/`:

- **`diff.py`** — compares two sets of per-testbed BMFs on `accuracy-digits`;
  a case is *changed* when `|Δ digits| ≥ 0.5` (configurable).
- **`heatmap.py`** — absolute and relative heatmap PNGs + an inline emoji grid.
- **`report.py`** — the Check summary and the upserted PR comment body.
- **`dashboard.py` / `pr_report.py`** — CLIs the workflow runs.

On **develop** push: the absolute heatmap + `baseline/*.bmf.json` are published
to the `gh-pages` branch (the standing dashboard). On a **PR**: the baseline is
fetched from `gh-pages`, diffed against the PR's BMFs, and a non-failing
**Check** plus an always-upserted **comment** (emoji grid + itemized
improvements/regressions, capped at 10/group, links to the absolute + relative
heatmap PNGs) are posted. Fork PRs get the emoji grid only (no `gh-pages` write).

Scope: P1. The convergence-curve view (err vs N) is P2 — see the spec.
```

- [ ] **Step 3: Verify references resolve**

Run:
```bash
cd /workspaces/nfft2
for f in tests/bench/diff.py tests/bench/heatmap.py tests/bench/report.py \
         tests/bench/dashboard.py tests/bench/pr_report.py \
         .github/scripts/gh-pages-publish.sh .github/scripts/accuracy-pr-report.sh; do
  test -f "$f" || { echo "MISSING $f"; exit 1; }
done
echo "all referenced files exist"
```
Expected: `all referenced files exist`.

- [ ] **Step 4: Commit**
```bash
git add docs/agents/accuracy-tracking.md tests/.gitignore
git commit -m "docs: document the accuracy reporting layer (dashboard + PR signal)"
```

---

## Self-Review

**Spec coverage:**
- Heatmap dashboard on Pages (absolute) → Tasks 3, 4 (dashboard.py), 5 (develop publish). ✓
- PR Check (non-failing) + always-on bidirectional comment → Tasks 2, 4 (pr_report.py), 5 (post). ✓
- Two PR heatmaps (absolute + relative) → Task 3 + pr_report.py. ✓
- Inline emoji grid + linked PNGs; fork PRs emoji-only → `emoji_grid`, `comment_body(png_urls)`, Task 5 fork branch. ✓
- Self-contained diff (Approach 1), no Bencher API in PR path → `diff.py` + baseline fetched from gh-pages raw. ✓
- Significance gate 0.5, cap 10, sorted by magnitude → `diff.diff(gate)`, `report.CAP`, sort in `diff`. ✓
- Suppress Bencher comment → Task 5 Step 1. ✓
- Bencher stays as archive → `bencher-upload` job untouched except `--github-actions`. ✓
- P2 out of scope → stated in Global Constraints + docs. ✓

**Placeholder scan:** Only intentional, instruction-backed tokens: `__CHECKOUT_SHA__`/`__UV_SHA__`/`__DLA_SHA__` (resolved inline in Task 5 Step 2 with exact values). No "TBD"/"add error handling"/"similar to".

**Type consistency:** `Change`/`DiffResult` fields identical across `diff.py`, `report.py`, `heatmap.py`, and tests. `comment_body(result, png_urls)`, `check_summary(result)`, `render_absolute(tree,out)`, `render_relative(result,out)`, `emoji_grid(result)`, `load_bmf_tree(dir)`, `diff(pr_tree,base_tree,gate)` signatures consistent between definitions, CLIs, and tests. BMF access path `["accuracy-digits"]["value"]` matches the converter's output from `ndjson_to_bmf.py`. Artifact naming `accuracy-bmf-<testbed>` matches the existing upload job and Task 5's arranger.
