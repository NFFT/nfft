import re

from tests.accuracy.diff import Change, DiffResult
from tests.accuracy.htmlreport import (
    GREY,
    RED,
    YELLOW,
    aggregate_by_module,
    cell_text,
    margin_color,
    metric_module,
    parse_window,
    render_report,
)


def test_parse_window_extracts_and_strips_known_window():
    assert parse_window("ubuntu-latest_gcc_kaiserbessel_double") == (
        "kaiserbessel",
        "ubuntu-latest_gcc_double",
    )


def test_parse_window_unknown_is_other():
    assert parse_window("ubuntu-latest_gcc_double") == (
        "other",
        "ubuntu-latest_gcc_double",
    )


def test_metric_module_splits_first_segment():
    assert metric_module("nfft/serial/file/fast/forward/1d/a") == (
        "nfft",
        "serial/file/fast/forward/1d/a",
    )


def test_margin_color_bands_are_hex():
    assert margin_color(-1.0) == RED
    assert margin_color(-5.0) == RED
    assert margin_color(0.0) == YELLOW
    assert margin_color(0.49) == YELLOW
    light = margin_color(0.5)
    dark = margin_color(3.0)
    assert light != dark
    assert margin_color(10.0) == dark  # clamped beyond 3.0
    for c in (RED, YELLOW, GREY, light, dark):
        assert len(c) == 7 and c[0] == "#"


def test_green_ramp_darkens():
    def chan_sum(hexc):
        return int(hexc[1:3], 16) + int(hexc[3:5], 16) + int(hexc[5:7], 16)

    assert chan_sum(margin_color(3.0)) < chan_sum(margin_color(0.5))


def test_cell_text_digits_and_bound():
    assert cell_text(13.72, 11.6, 1e-13) == "13.7 (11.6)"


def test_cell_text_perfect_is_infinity():
    assert cell_text(300.0, 11.6, 0.0) == "∞ (11.6)"


# --- render --------------------------------------------------------------


def _cell(digits, bound_digits, max_error=None):
    if max_error is None:
        max_error = 10.0 ** (-digits)
    return {
        "accuracy-digits": {"value": digits},
        "max-error": {"value": max_error},
        "bound-digits": {"value": bound_digits},
    }


def _tree():
    return {
        "ci_gcc_kaiserbessel_double": {
            "nfft/serial/file/fast/forward/1d/a": _cell(13.7, 11.0),
            "nfct/serial/file/fast/forward/1d/a": _cell(12.0, 11.0),
        },
        "ci_gcc_kaiserbessel_float": {
            # 1d metric missing here -> grey cell; only 2d present
            "nfft/serial/file/fast/forward/2d/a": _cell(6.0, 5.0),
        },
        "ci_gcc_gaussian_double": {
            "nfft/serial/file/fast/forward/1d/a": _cell(9.0, 8.0),
        },
    }


def test_render_absolute_has_a_tab_per_window():
    doc = render_report(_tree(), title="T")
    assert 'id="tab-kaiserbessel"' in doc
    assert 'id="tab-gaussian"' in doc
    assert 'id="panel-kaiserbessel"' in doc
    assert "<script" not in doc  # no JavaScript
    assert doc.lstrip().startswith("<!doctype html>")


def test_render_absolute_tab_ids_are_unique():
    doc = render_report(_tree(), title="T")
    ids = re.findall(r'id="(tab-[^"]+)"', doc)
    assert len(ids) == len(set(ids))


def test_render_absolute_has_table_per_module():
    doc = render_report(_tree(), title="T")
    assert "nfft" in doc and "nfct" in doc
    assert "1d/a" in doc and "2d/a" in doc  # union of submetrics


def test_render_absolute_missing_cell_is_grey():
    doc = render_report(_tree(), title="T")
    assert GREY in doc


def test_render_absolute_colors_a_passing_cell_green():
    doc = render_report(_tree(), title="T")
    assert margin_color(2.7) in doc  # margin 13.7 - 11.0


# --- PR diff view --------------------------------------------------------


def _diff():
    return DiffResult(
        improvements=[
            Change(
                "ci_gcc_kaiserbessel_double",
                "nfft/serial/file/fast/forward/1d/a",
                13.0,
                14.5,
                1.5,
                11.0,
            )
        ],
        regressions=[
            Change(
                "ci_gcc_kaiserbessel_double",
                "nfct/serial/file/fast/forward/1d/a",
                14.0,
                12.5,
                -1.5,
                -11.0,
            )
        ],
        unchanged_count=3,
        added=["ci_gcc_kaiserbessel_double::nfft/serial/file/fast/forward/3d/a"],
        removed=[],
        by_testbed={},
    )


def test_aggregate_by_module_counts_per_module():
    agg = aggregate_by_module(_diff())
    assert agg["nfft"] == {"improved": 1, "regressed": 0}
    assert agg["nfct"] == {"improved": 0, "regressed": 1}


def test_pr_view_has_changes_section_with_itemized_list():
    doc = render_report(_tree(), _diff(), title="PR")
    assert "Changes" in doc
    assert "13.00" in doc and "14.50" in doc
    assert "added" in doc.lower()


def test_pr_view_marks_changed_cell():
    doc = render_report(_tree(), _diff(), title="PR")
    assert 'class="chg"' in doc
    assert "▲" in doc
    # Changed cell shows baseline → new (bound) using cell's own digits (13.7 from tree).
    # In production change.pr_digits == cell digits; the mock has them mismatched on purpose
    # so we test with the cell value that actually ends up in the HTML.
    assert "13.0 → 13.7 (11.0) ▲" in doc
    # Regression: base=14.0, cell digits=12.0 (bound=11.0)
    assert "14.0 → 12.0 (11.0) ▼" in doc


def test_changes_itemized_capped_at_10():
    big = DiffResult(
        improvements=[
            Change("t", f"nfft/serial/file/fast/forward/{k}d/a", 13.0, 14.0, 1.0, 7.0)
            for k in range(13)
        ],
        regressions=[],
        unchanged_count=0,
        added=[],
        removed=[],
        by_testbed={},
    )
    doc = render_report(_tree(), big, title="PR")
    assert "+3 more" in doc
