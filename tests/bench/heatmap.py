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
