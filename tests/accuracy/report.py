"""Render the PR Check summary and the upserted comment body from a DiffResult.

The comment is a lightweight per-module rollup (improved/regressed counts) plus a
link to the full HTML accuracy report on GitHub Pages. Per-case detail lives in
the report, not the comment.
"""

from __future__ import annotations

from tests.accuracy.htmlreport import aggregate_by_module

MARKER = "<!-- nfft-accuracy-report -->"
CAP = 10


def _counts(result):
    y, z = len(result.improvements), len(result.regressions)
    na, nr = len(result.added), len(result.removed)
    line = f"{result.unchanged_count} unchanged · {y} improved · {z} regressed"
    if na or nr:
        line += f" · {na} added · {nr} removed"
    return y, z, na, nr, line


# A renamed or regrouped metric shows up as one added plus one removed name, and
# the max-merged neighbours it left change too, so those counts are not accuracy.
def _regroup_note(na, nr):
    if not (na or nr):
        return ""
    return (
        f"The metric set changed ({na} added, {nr} removed), so improved and "
        "regressed counts may reflect regrouping rather than accuracy changes."
    )


def check_summary(result):
    y, z, na, nr, summary = _counts(result)
    if y == 0 and z == 0:
        title = "accuracy unchanged"
    else:
        title = f"{z} regressed, {y} improved"
    if na or nr:
        title += f", {na} added, {nr} removed"
    return "neutral", title, summary


def check_summary_no_baseline():
    return "neutral", "baseline pending", "no develop baseline yet to compare against"


def _report_link(report_url):
    return f"📊 [Full accuracy report]({report_url})" if report_url else ""


def _module_table(result):
    agg = aggregate_by_module(result)
    lines = ["| module | improved | regressed |", "|---|---|---|"]
    for module in sorted(agg):
        a = agg[module]
        lines.append(f"| `{module}` | {a['improved']} | {a['regressed']} |")
    return "\n".join(lines)


def comment_body(result, report_url):
    y, z, na, nr, summary = _counts(result)
    note = _regroup_note(na, nr)
    if y == 0 and z == 0:
        body = f"{MARKER}\nAccuracy: {result.unchanged_count} cases unchanged."
        if note:
            body += f" {na} metrics added, {nr} removed."
        link = _report_link(report_url)
        return (body + (f"\n\n{link}" if link else "")).rstrip() + "\n"
    lines = [MARKER, "## Accuracy report", summary, ""]
    if note:
        lines += [note, ""]
    lines += [_module_table(result), ""]
    link = _report_link(report_url)
    if link:
        lines.append(link)
    return "\n".join(lines).rstrip() + "\n"


def comment_body_no_baseline(report_url):
    lines = [
        MARKER,
        "## Accuracy report",
        "",
        "No `develop` baseline yet to compare against — showing this PR's "
        + "absolute accuracy in the full report. A baseline appears once "
        + "changes land on `develop`.",
    ]
    link = _report_link(report_url)
    if link:
        lines += ["", link]
    return "\n".join(lines).rstrip() + "\n"
