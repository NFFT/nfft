"""Render accuracy reports as a self-contained HTML document (stdlib only).

CSS-only window tabs (hidden radio inputs), one <table> per module per tab,
testbeds as columns, submetrics as rows. Cells are colored by their margin
(accuracy-digits minus bound-digits). No JavaScript.
"""

from __future__ import annotations

import html

# Window functions whose token may appear in a testbed name
# (`<os>_<compiler>_<window>_<precision>`). Used to split the window axis out.
KNOWN_WINDOWS = ("kaiserbessel", "gaussian", "bspline", "sinc", "dirac")
# Precision tokens of a testbed name, in mantissa order; columns sort by it.
PRECISIONS = ("float", "double", "long-double")

# Margin color bands. Hex equals round(255*c) of the old matplotlib RGBA floats.
RED = "#d63026"  # margin < 0      (error exceeds bound -- should never happen)
YELLOW = "#fcd94c"  # 0 <= margin <0.5 (barely passing)
GREY = "#d9d9d9"  # missing cell    (metric absent for this testbed)
# Green ramp endpoints (light -> dark), from the old LinearSegmentedColormap.
# The dark end is lightened from the original (0.0, 0.27, 0.106) = #00451b so the
# black cell text stays legible: #29944c gives ~5.4:1 contrast on black vs ~2.3:1
# before (WCAG AA needs 4.5:1). It is still clearly the dark end of the ramp.
_GREEN_LO_RGB = (0.78, 0.91, 0.55)
_GREEN_HI_RGB = (0.16, 0.58, 0.30)
GREEN_LO, GREEN_HI = 0.5, 3.0  # margin range the ramp spans, then clamps

CAP = 10  # max itemized rows per group in the Changes section


def _hex(rgb):
    return "#%02x%02x%02x" % tuple(max(0, min(255, round(c * 255))) for c in rgb)


def _lerp(a, b, t):
    return tuple(a[i] + (b[i] - a[i]) * t for i in range(3))


def parse_window(testbed):
    """Split a testbed name into (window, label-without-window)."""
    tokens = testbed.split("_")
    for i, tok in enumerate(tokens):
        if tok in KNOWN_WINDOWS:
            rest = tokens[:i] + tokens[i + 1 :]
            return tok, "_".join(rest)
    return "other", testbed


def precision_rank(label):
    """Position of the label's precision token in PRECISIONS; unknown last."""
    tokens = label.split("_")
    for i, p in enumerate(PRECISIONS):
        if p in tokens:
            return i
    return len(PRECISIONS)


def metric_module(name):
    """Split a metric name into (module, rest) on the first '/'."""
    module, _, rest = name.partition("/")
    return module, rest


def margin_color(margin):
    """Map a margin (accurate digits beyond the bound) to a CSS hex color."""
    if margin < 0.0:
        return RED
    if margin < GREEN_LO:
        return YELLOW
    t = (min(margin, GREEN_HI) - GREEN_LO) / (GREEN_HI - GREEN_LO)
    return _hex(_lerp(_GREEN_LO_RGB, _GREEN_HI_RGB, t))


def cell_text(digits, bound_digits, max_error):
    """Cell label: accurate digits with the bound (required digits) in parens.

    A perfect result (max_error == 0) shows the infinity glyph for the digits.
    """
    head = "∞" if max_error == 0.0 else f"{digits:.1f}"
    return f"{head} ({bound_digits:.1f})"


_STYLE = """
body{font:13px/1.4 system-ui,sans-serif;margin:1.5rem;color:#111}
h1{font-size:1.1rem}
p.legend{font-size:11px;color:#444}
.legend span{display:inline-block;width:1.1em;height:1.1em;vertical-align:-2px;
  margin:0 .2em 0 .8em;border:1px solid #0002}
.tabs input{position:absolute;opacity:0;pointer-events:none}
.tabs label{display:inline-block;padding:.3em .7em;margin:0 .2em;cursor:pointer;
  border:1px solid #ccc;border-bottom:none;border-radius:4px 4px 0 0;background:#f4f4f4}
.tabs input:checked + label{background:#fff;font-weight:600;border-color:#888}
.panel{display:none;border:1px solid #888;padding:.8rem;border-radius:0 4px 4px 4px}
table.hm{border-collapse:collapse;margin:.4rem 0 1.2rem}
table.hm caption{text-align:left;font-weight:600;margin:.2em 0}
table.hm th{font:11px monospace;font-weight:400;text-align:right;padding:1px 4px;
  white-space:nowrap}
table.hm thead th{text-align:center}
table.hm td{font:10px monospace;text-align:center;padding:1px 5px;
  border:1px solid #fff;white-space:nowrap}
td.miss{color:#0000}
td.chg{outline:2px solid #1a1a1a;outline-offset:-2px;font-weight:700}
"""

# Each window radio: `:checked + label` styles the tab; `:checked ~ #panel-<win>`
# reveals the panel. All inputs+labels are emitted first, then all panels, so the
# general-sibling selector reaches every panel.
_PANEL_SEL = "\n".join(
    f"#tab-{w}:checked ~ #panel-{w}{{display:block}}"
    for w in list(KNOWN_WINDOWS) + ["other"]
)


def _esc(s):
    return html.escape(str(s), quote=True)


def _module_subs(tree, testbeds):
    """module -> sorted union of within-module submetric names across testbeds."""
    per_module = {}
    for tb, _label in testbeds:
        for name in tree[tb]:
            module, sub = metric_module(name)
            per_module.setdefault(module, set()).add(sub)
    return {m: sorted(subs) for m, subs in per_module.items()}


def _cell_html(cell, change=None):
    """One <td>. `cell` is the BMF entry or None (missing). `change` marks a diff."""
    if cell is None:
        return '<td class="miss" style="background:%s"></td>' % GREY
    digits = cell["accuracy-digits"]["value"]
    bound_digits = cell["bound-digits"]["value"]
    max_error = cell["max-error"]["value"]
    color = margin_color(digits - bound_digits)
    text = cell_text(digits, bound_digits, max_error)
    cls, tip = "", f"err-digits={digits:.2f} bound-digits={bound_digits:.2f}"
    if change is not None:
        arrow = "▲" if change.delta_digits > 0 else "▼"
        cls = ' class="chg"'
        pr_head = "∞" if max_error == 0.0 else f"{digits:.1f}"
        text = f"{change.base_digits:.1f} → {pr_head} ({bound_digits:.1f}) {arrow}"
        tip = (
            f"{change.base_digits:.2f} -> {change.pr_digits:.2f} digits "
            f"(delta {change.delta_digits:+.2f})"
        )
    return f'<td{cls} style="background:{color}" title="{_esc(tip)}">{_esc(text)}</td>'


def _module_table(module, window, subs, testbeds, tree, changes):
    head = "".join(f"<th>{_esc(lbl)}</th>" for _tb, lbl in testbeds)
    rows = [f"<tr><th></th>{head}</tr>"]
    for sub in subs:
        name = f"{module}/{sub}"
        cells = []
        for tb, _lbl in testbeds:
            cell = tree[tb].get(name)
            change = changes.get((tb, name)) if changes else None
            cells.append(_cell_html(cell, change))
        rows.append(f"<tr><th>{_esc(sub)}</th>{''.join(cells)}</tr>")
    return (
        f'<table class="hm"><caption>{_esc(module)} · {_esc(window)}'
        f"</caption><thead>{rows[0]}</thead><tbody>"
        f"{''.join(rows[1:])}</tbody></table>"
    )


def _changes_index(diff_result):
    """{(testbed, metric_name): Change} over improvements + regressions."""
    if diff_result is None:
        return {}
    out = {}
    for c in diff_result.improvements + diff_result.regressions:
        out[(c.testbed, c.name)] = c
    return out


def aggregate_by_module(diff_result):
    """{module: {"improved","regressed"}} from a DiffResult."""
    agg = {}
    for c in diff_result.improvements:
        agg.setdefault(metric_module(c.name)[0], {"improved": 0, "regressed": 0})[
            "improved"
        ] += 1
    for c in diff_result.regressions:
        agg.setdefault(metric_module(c.name)[0], {"improved": 0, "regressed": 0})[
            "regressed"
        ] += 1
    return agg


def _change_items(changes):
    out = []
    for c in changes[:CAP]:
        out.append(
            f"<li><code>{_esc(c.testbed)} {_esc(c.name)}</code>: "
            f"{c.base_digits:.2f} → {c.pr_digits:.2f} digits "
            f"({c.pct:+.0f}%)</li>"
        )
    if len(changes) > CAP:
        out.append(f"<li>… +{len(changes) - CAP} more</li>")
    return "".join(out)


def _render_changes(diff_result):
    y, z = len(diff_result.improvements), len(diff_result.regressions)
    na, nr = len(diff_result.added), len(diff_result.removed)
    agg = aggregate_by_module(diff_result)
    rows = "".join(
        f"<tr><td><code>{_esc(m)}</code></td><td>{agg[m]['improved']}</td>"
        f"<td>{agg[m]['regressed']}</td></tr>"
        for m in sorted(agg)
    )
    parts = [
        '<section class="changes"><h2>Changes</h2>',
        f"<p>{diff_result.unchanged_count} unchanged &middot; {y} improved "
        f"&middot; {z} regressed",
    ]
    if na or nr:
        parts.append(f" &middot; {na} added &middot; {nr} removed")
    parts.append("</p>")
    if rows:
        parts.append(
            '<table class="hm"><thead><tr><th>module</th>'
            "<th>improved</th><th>regressed</th></tr></thead>"
            f"<tbody>{rows}</tbody></table>"
        )
    if diff_result.improvements:
        parts.append(
            f"<h3>Improvements</h3><ul>{_change_items(diff_result.improvements)}</ul>"
        )
    if diff_result.regressions:
        parts.append(
            f"<h3>Regressions</h3><ul>{_change_items(diff_result.regressions)}</ul>"
        )
    parts.append("</section>")
    return "".join(parts)


_LEGEND = (
    '<p class="legend">Correct decimal digits per case (bound in parens), colored '
    "by margin beyond each case’s bound: "
    f'<span style="background:{RED}"></span>over bound '
    f'<span style="background:{YELLOW}"></span>barely passing '
    f'<span style="background:{margin_color(0.5)}"></span>→'
    f'<span style="background:{margin_color(3.0)}"></span>healthy '
    f'<span style="background:{GREY}"></span>missing &middot; ∞ = exact.</p>'
)


def render_report(tree, diff_result=None, *, title):
    """Return a self-contained HTML document for `tree`.

    Absolute (develop) view when diff_result is None; PR view (changed cells
    marked + a Changes section on top) when a DiffResult is given.
    """
    changes = _changes_index(diff_result)
    windows = {}
    for tb in tree:
        win, label = parse_window(tb)
        windows.setdefault(win, []).append((tb, label))

    order = [w for w in (list(KNOWN_WINDOWS) + ["other"]) if w in windows]
    inputs, panels = [], []
    for idx, win in enumerate(order):
        testbeds = sorted(windows[win], key=lambda x: (precision_rank(x[1]), x[1]))
        checked = " checked" if idx == 0 else ""
        inputs.append(
            f'<input type="radio" name="wintab" id="tab-{win}"{checked}>'
            f'<label for="tab-{win}">{_esc(win)}</label>'
        )
        subs_by_module = _module_subs(tree, testbeds)
        tables = "".join(
            _module_table(m, win, subs_by_module[m], testbeds, tree, changes)
            for m in sorted(subs_by_module)
        )
        panels.append(f'<section class="panel" id="panel-{win}">{tables}</section>')

    changes_section = _render_changes(diff_result) if diff_result is not None else ""
    return (
        "<!doctype html>\n<html><head><meta charset=utf-8>"
        f"<title>{_esc(title)}</title><style>{_STYLE}{_PANEL_SEL}</style></head>"
        f"<body><h1>{_esc(title)}</h1>{_LEGEND}{changes_section}"
        f'<div class="tabs">{"".join(inputs)}{"".join(panels)}</div>'
        "</body></html>\n"
    )
