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


def check_summary_no_baseline():
    """Check output when no develop baseline exists yet (e.g. the first PR)."""
    return "neutral", "baseline pending", "no develop baseline yet to compare against"


def comment_body_no_baseline(abs_url=None):
    """Comment when there is no baseline: show absolute accuracy only, no diff."""
    lines = [MARKER, "## Accuracy report", "",
             "No `develop` baseline yet to compare against — showing this PR's "
             "absolute accuracy. A baseline appears once changes land on `develop`."]
    if abs_url:
        lines += ["", f"[absolute heatmap]({abs_url})"]
    return "\n".join(lines).rstrip() + "\n"


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
