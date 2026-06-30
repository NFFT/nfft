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
    # Exactly 10 itemized rows are rendered (the cap), plus a "+N more" note.
    item_lines = [ln for ln in body.splitlines() if ln.startswith("- `t1 i")]
    assert len(item_lines) == 10
    assert "+3 more" in body
