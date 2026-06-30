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
