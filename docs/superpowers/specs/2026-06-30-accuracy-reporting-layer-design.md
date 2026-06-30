# Accuracy Reporting Layer (P1) — Design

**Status:** approved (brainstorming), ready for implementation planning
**Date:** 2026-06-30
**Related:** [ADR-0004](../../adr/0004-accuracy-tracking-with-bencher.md), [`docs/agents/accuracy-tracking.md`](../../agents/accuracy-tracking.md)

## 1. Background

The accuracy-tracking pipeline already works: `build-linux.yml`'s gcc cells run
`make check` with `NFFT_BENCH_OUT` set, emit per-case NDJSON, aggregate it with
`tests/bench/ndjson_to_bmf.py` into one BMF file per testbed, and upload to the
Bencher project `nfft`. Each metric is named
`<module>/<runtime>/<oracle>/<speed>/<direction>/<dim>d/<init-slug>` with two
measures: **`accuracy-digits`** = `-log10(worst err)` (primary, higher = better)
and **`max-error`** (raw worst error). There are ~150–290 metrics per testbed
across 12 testbeds (4 windows × 3 precisions).

**The problem this design solves:** Bencher's presentation does not fit this
project. Its PR bot comment dumps every metric (hundreds of rows) on every PR —
overwhelming, and the wrong model for a numerical library where most PRs don't
touch accuracy at all. Its per-metric line-plot dashboard does not scale to this
many high-cardinality series, so there is no good "scan for weak spots" view.

## 2. Goals / non-goals

**Goals (P1):**
- A readable **standing dashboard** (GitHub Pages) showing the absolute accuracy
  state of the `develop` baseline as a heatmap.
- A quiet, **bidirectional PR signal**: a one-line GitHub Check plus a single
  always-present comment that is compact when accuracy is unchanged and itemizes
  improvements *and* regressions when it moves.
- Keep **Bencher purely as the long-term data archive** (and future statistical
  threshold engine); stop using its PR comment.

**Non-goals (explicitly deferred):**
- **P2 — sweep-structured tests + scaling curves.** Re-designing the refgen grids
  so each metric is a controlled single-parameter sweep (e.g. err vs `N` at fixed
  everything else), and adding "err vs parameter" convergence-curve plots. That
  touches the test oracle (refgen, reference-data regeneration, analytic bounds)
  and gets its own brainstorm + spec. P1's heatmap works on the current grids and
  keeps working after P2; P2 only *adds* the curve view.
- Failing CI on accuracy changes (stays track-only; the Check never fails).
- Changing the metric structure, the emitter, or the converter measures.

## 3. Key decisions (from brainstorming, with rationale)

| Decision | Choice | Why |
|---|---|---|
| Keep Bencher? | **Yes — as archive only** | Detection/storage are fine; only its *presentation* is wrong. |
| Diff data source | **Self-contained** (Approach 1) | The PR diff is a pure function of two BMF files → testable, robust, fork-safe; no Bencher API/token in the PR path. |
| Dashboard home | **GitHub Pages** | Stable URL, scannable, room to grow (P2 curves slot in). User will enable Pages. |
| Dashboard view | **Heatmap** | Fastest "scan for hotspots" view for high-cardinality data; line plots don't scale. |
| PR surfaces | **Check + always-on comment** | Check = at-a-glance status; comment = detail, compact when flat. |
| PR comment policy | **Always posted, upserted** | Consistent presence; one line when unchanged, itemized when moved. |
| PR heatmaps | **Two: absolute + relative-change** | Absolute = where the PR stands; relative = where it moved things. |
| Inline figure | **Emoji heatmap (relative) + linked hi-res PNGs** | Self-contained inline (fork-safe), real figures linked from Pages. |
| Bencher PR comment | **Suppressed** | Replaced by our Check + comment. |
| "Changed" gate | **\|Δ accuracy-digits\| ≥ 0.5** (configurable) | ~3× error change; above `-ffast-math`/compiler noise. |

## 4. Architecture & data flow

```
build-linux gcc cells (existing, unchanged)
  make check + NFFT_BENCH_OUT  →  per-testbed BMF (existing converter)
        │                                  │
        ├──► bencher-upload job (existing): archive to Bencher (drop --github-actions)
        │
        └──► accuracy-report job (NEW, needs: build): downloads all accuracy-bmf-* artifacts
               ├─ develop push:  render ABSOLUTE heatmaps + copy baseline BMFs → publish to gh-pages
               └─ pull_request:  GET baseline BMFs from Pages → diff → render ABSOLUTE + RELATIVE heatmaps
                                 → post Check (counts) + upsert comment (emoji + itemized + PNG links)
```

The entire report is a pure transform of BMF JSON; rendering and diffing never
call the Bencher API.

## 5. Components

All new Python lives under `tests/bench/` next to the existing converter, and is
pure and unit-testable. Rendering uses `matplotlib` (run via `uv`, consistent
with how the repo already invokes Python tooling in CI).

### 5.1 `tests/bench/diff.py`
Pure comparison of two BMF objects.

- `load_bmf(path) -> dict` — read one BMF file.
- `diff(pr_bmf, baseline_bmf, gate=0.5) -> DiffResult` where `DiffResult` carries:
  - `improvements: list[Change]`, `regressions: list[Change]` (each a metric whose
    `|Δ accuracy-digits| ≥ gate`), `unchanged_count: int`,
  - each `Change` = `{name, base_digits, pr_digits, delta_digits, pct}` with
    `pct = (pr_digits - base_digits) / base_digits * 100`,
  - both lists **sorted by `|delta_digits|` descending**,
  - `by_testbed`: per-testbed counts for the heatmap.
- Metrics present in only one side (added/removed cases) are reported in a
  separate `added`/`removed` list, not counted as improve/regress.

### 5.2 `tests/bench/heatmap.py`
One renderer, two color modes; emits both a raster figure and an inline markdown
grid.

- `render_absolute(bmf_by_testbed, out_png)` — cell = worst-case `accuracy-digits`,
  colormap green(high)→red(low).
- `render_relative(diff_result, out_png)` — cell = signed `delta_digits`,
  diverging colormap green(+)/neutral(0)/red(−).
- `emoji_grid(diff_result) -> str` — a compact markdown table of 🟩/🟥/· (shaded by
  magnitude, e.g. 🟩/💚 tiers) for inline use in the comment.
- **Default layout** (tunable in implementation): faceted per module
  (`nfft`/`nfct`/`nfst`); rows = cases (`oracle/speed/direction/dim/init`);
  columns = the 12 testbeds × runtime (`serial`/`omp`). Missing cells rendered
  blank.

### 5.3 `tests/bench/report.py`
Assembles the PR surfaces from a `DiffResult`.

- `check_summary(diff) -> (conclusion, title, summary_md)` — `conclusion` is always
  non-failing (`neutral`/`success`); `title` e.g. `accuracy unchanged` or
  `3 regressed, 2 improved`.
- `comment_body(diff, png_urls) -> str` — the full markdown:
  - **Flat:** single line `Accuracy: N cases unchanged.`
  - **Changed:** the inline emoji relative heatmap; an aggregate line
    `X unchanged · Y improved · Z regressed`; an **Improvements** group then a
    **Regressions** group, each listing cases as
    `‹name›: ‹base›→‹pr› digits (‹±pct›%)`, capped at **10** with `… +N more`;
    and links to the absolute + relative detailed heatmap PNGs.
  - When `png_urls` is empty (fork PR, no archive write), omit the links; the emoji
    heatmap and itemized lists remain.
- `check_summary_no_baseline()` / `comment_body_no_baseline(abs_url=None)` — the
  **no-baseline** state (the first PR, or any PR opened before `develop` has
  published a baseline). The Check title is `baseline pending` (still `neutral`);
  the comment says plainly "No `develop` baseline yet … showing this PR's absolute
  accuracy", links the absolute heatmap (if `abs_url`), and shows **no diff** — it
  must NOT report a misleading "unchanged".

### 5.4 Archiving & publishing (gh-pages branch)
The heatmaps are committed to a **`gh-pages` branch** — chosen deliberately over
`actions/deploy-pages` (which replaces the whole site each deploy, so per-PR
history doesn't accumulate) and over workflow artifacts (which expire ~90 days →
stale links). The branch gives **permanent, accumulating per-PR archives**.

- **develop push:** render absolute heatmaps + a minimal `index.html`, and copy
  the per-testbed baseline BMFs to `baseline/<testbed>.bmf.json` (so PR runs can
  fetch them). The human dashboard is served by **GitHub Pages** from this branch.
- **same-repo PR:** publish the PR's heatmap PNGs to `pr/<number>/…` on `gh-pages`.
  The comment links them via **`raw.githubusercontent.com/<owner>/<repo>/gh-pages/
  pr/<number>/…`** — permanent, render inline, and resolve **without Pages enabled
  or built** (also why the first PR works before Pages is configured).
- **no baseline yet:** same-repo PR still archives + links the **absolute** heatmap
  only (no relative diff).
- **fork PR:** GitHub forces a read-only `GITHUB_TOKEN` for fork `pull_request`
  runs regardless of `permissions:`/environment approval — so forks can write
  neither `gh-pages` nor a Check/comment. The `accuracy-report` job is therefore
  **skipped entirely for fork PRs** (job-level `if`), so it can't fail trying.
  Fork PRs get no accuracy report. (A future `workflow_run`-triggered companion,
  running with the base-repo token after the untrusted build, could add it — out
  of scope for P1.)

### 5.5 Workflow wiring (`accuracy-report` job in `build-linux.yml`)
- `needs: build`; downloads all `accuracy-bmf-*` artifacts (same as the upload job).
- Uses **only `GITHUB_TOKEN`** (no external secret), so **no environment gating**.
  Permissions: `contents: write` (push `gh-pages`), `checks: write`, `pull-requests:
  write`. `persist-credentials: false` (the publish script authenticates via its own
  clone URL).
- **Job-level `if` skips fork PRs** (read-only token can't write); runs for pushes,
  same-repo PRs, and `workflow_dispatch`.
- `push` → publish dashboard + baseline **only on `refs/heads/develop`** (main/master
  pushes must not clobber the develop baseline). `pull_request` → diff + Check +
  comment (+ PR PNG archive).
- **Both report steps are `continue-on-error: true`** — a publish/post/network
  failure is informational and must never red CI (the Check *conclusion* is also
  always `neutral`). The Check is posted on the **PR head SHA**
  (`github.event.pull_request.head.sha`), not the ephemeral merge SHA.
- The existing `bencher-upload` job stays; remove its `--github-actions` flag so
  Bencher no longer comments.

## 6. Testing

- `diff.py`: golden pytest over crafted BMF pairs — improvement/regression
  classification, the 0.5 gate boundary, `pct` math, sort order, added/removed
  handling, all-unchanged case.
- `report.py`: golden pytest of `comment_body` for the flat case and a changed
  case (assert aggregate line, grouping, sort, 10-cap `… +N more`, link presence
  vs absence), and `check_summary` conclusions.
- `heatmap.py`: smoke test — renders a valid PNG of expected size from sample
  data; `emoji_grid` returns a well-formed markdown table.
- Workflow: `yamllint`; a local end-to-end dry-run of the report on real BMFs
  (generate via `tests/checkall`, diff a perturbed copy, render, build the comment
  markdown) without posting.

## 7. Configuration defaults

- Significance gate: `0.5` accuracy-digits (env/flag configurable).
- Comment cap: `10` cases per group.
- Check: never fails CI (track-only).
- Heatmap: faceted per module.

## 8. Out of scope

P2 (sweep-structured tests + scaling curves), any change to the emitter/converter
measures, statistical thresholds (left to Bencher for later), and failing CI on
accuracy regressions.
