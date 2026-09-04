from tests.accuracy.diff import Change, DiffResult
from tests.accuracy.report import (
    MARKER,
    check_summary,
    check_summary_no_baseline,
    comment_body,
    comment_body_no_baseline,
)


def _result(impr=0, regr=0, unchanged=5, added=0, removed=0):
    imps = [
        Change("t1", f"nfft/i{k}", 13.0, 14.0 + k, 1.0 + k, 7.0) for k in range(impr)
    ]
    regs = [
        Change("t1", f"nfct/r{k}", 14.0, 13.0 - k, -(1.0 + k), -7.0)
        for k in range(regr)
    ]
    return DiffResult(
        imps,
        regs,
        unchanged,
        added=[f"t1::nfft/a{k}" for k in range(added)],
        removed=[f"t1::nfft/d{k}" for k in range(removed)],
        by_testbed={},
    )


def test_check_never_fails():
    assert check_summary(_result())[0] == "neutral"


def test_check_title_flat_vs_changed():
    assert "unchanged" in check_summary(_result())[1].lower()
    assert check_summary(_result(impr=2, regr=1))[1] == "1 regressed, 2 improved"


def test_comment_marker_first_line():
    assert comment_body(_result(), None).splitlines()[0] == MARKER


def test_flat_comment_is_one_line_summary():
    body = comment_body(_result(unchanged=42), None)
    assert "42 cases unchanged" in body
    assert "| module |" not in body  # no table when nothing changed


def test_changed_comment_has_per_module_table_and_link():
    body = comment_body(_result(impr=1, regr=1), "https://x.github.io/y/pr/9/")
    assert "5 unchanged" in body and "1 improved" in body and "1 regressed" in body
    assert "| module | improved | regressed |" in body
    assert "`nfft`" in body and "`nfct`" in body
    assert "https://x.github.io/y/pr/9/" in body
    # the comment NO LONGER carries itemized per-case lists
    assert "Improvements" not in body and "→" not in body


def test_no_baseline_comment_links_and_notes():
    body = comment_body_no_baseline("https://x.github.io/y/pr/9/")
    assert body.splitlines()[0] == MARKER
    assert "No `develop` baseline" in body
    assert "https://x.github.io/y/pr/9/" in body


def test_no_baseline_check_is_neutral_and_pending():
    conclusion, title, _ = check_summary_no_baseline()
    assert conclusion == "neutral" and "pending" in title


def test_added_removed_absent_when_metric_set_is_stable():
    _, title, summary = check_summary(_result(impr=2, regr=1))
    assert "added" not in title and "added" not in summary
    assert "regrouping" not in comment_body(_result(impr=2, regr=1), "u")


def test_added_removed_counts_in_check_and_comment():
    r = _result(impr=2, regr=1, added=3, removed=4)
    _, title, summary = check_summary(r)
    assert title == "1 regressed, 2 improved, 3 added, 4 removed"
    assert summary.endswith("· 3 added · 4 removed")
    body = comment_body(r, "u")
    assert "5 unchanged · 2 improved · 1 regressed · 3 added · 4 removed" in body
    assert "metric set changed (3 added, 4 removed)" in body
    assert "regrouping" in body


def test_flat_comment_mentions_added_removed():
    body = comment_body(_result(added=3, removed=4), "u")
    assert "5 cases unchanged. 3 metrics added, 4 removed." in body
