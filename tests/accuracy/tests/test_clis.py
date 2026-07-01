import json

from tests.accuracy import dashboard
from tests.accuracy import pr_report


def _write_bmf(path, **digit_by_name):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(
            {
                n: {
                    "accuracy-digits": {"value": v},
                    "max-error": {"value": 10**-v},
                    "bound-digits": {"value": v - 2.0},
                }
                for n, v in digit_by_name.items()
            },
            f,
        )


def test_dashboard_emits_html_and_baseline(tmp_path):
    src = tmp_path / "bmf"
    src.mkdir()
    _write_bmf(
        src / "ci_gcc_kaiserbessel_double.bmf.json",
        **{"nfft/serial/file/fast/forward/1d/a": 13.5},
    )
    out = tmp_path / "site"
    dashboard.main([str(src), str(out)])
    html = (out / "index.html").read_text(encoding="utf-8")
    assert html.lstrip().startswith("<!doctype html>")
    assert 'id="tab-kaiserbessel"' in html
    assert not list(out.glob("*.png"))
    assert (out / "baseline" / "ci_gcc_kaiserbessel_double.bmf.json").exists()


def test_pr_report_changed_writes_html_comment_check(tmp_path):
    pr = tmp_path / "pr"
    pr.mkdir()
    base = tmp_path / "base"
    base.mkdir()
    _write_bmf(
        pr / "ci_gcc_kaiserbessel_double.bmf.json",
        **{"nfft/serial/file/fast/forward/1d/a": 15.0},
    )
    _write_bmf(
        base / "ci_gcc_kaiserbessel_double.bmf.json",
        **{"nfft/serial/file/fast/forward/1d/a": 13.0},
    )
    out = tmp_path / "out"
    pr_report.main(
        [str(pr), str(base), str(out), "--report-url", "https://x.github.io/y/pr/9/"]
    )
    check = json.loads((out / "check.json").read_text())
    assert check["conclusion"] == "neutral" and "improved" in check["title"]
    body = (out / "comment.md").read_text()
    assert "| module | improved | regressed |" in body
    assert "https://x.github.io/y/pr/9/" in body
    html = (out / "index.html").read_text()
    assert 'class="chg"' in html  # changed cell marked
    assert not list(out.glob("*.png"))


def test_pr_report_no_baseline(tmp_path):
    pr = tmp_path / "pr"
    pr.mkdir()
    _write_bmf(
        pr / "ci_gcc_kaiserbessel_double.bmf.json",
        **{"nfft/serial/file/fast/forward/1d/a": 13.5},
    )
    out = tmp_path / "out"
    pr_report.main(
        [
            str(pr),
            str(pr),
            str(out),
            "--no-baseline",
            "--report-url",
            "https://x.github.io/y/pr/9/",
        ]
    )
    body = (out / "comment.md").read_text()
    assert "No `develop` baseline" in body
    assert "https://x.github.io/y/pr/9/" in body
    check = json.loads((out / "check.json").read_text())
    assert check["conclusion"] == "neutral" and "pending" in check["title"]
    assert (out / "index.html").exists()


def test_pr_report_flat_one_liner(tmp_path):
    pr = tmp_path / "pr"
    pr.mkdir()
    base = tmp_path / "base"
    base.mkdir()
    _write_bmf(
        pr / "ci_gcc_kaiserbessel_double.bmf.json",
        **{"nfft/serial/file/fast/forward/1d/a": 13.5},
    )
    _write_bmf(
        base / "ci_gcc_kaiserbessel_double.bmf.json",
        **{"nfft/serial/file/fast/forward/1d/a": 13.5},
    )
    out = tmp_path / "out"
    pr_report.main([str(pr), str(base), str(out)])
    body = (out / "comment.md").read_text()
    assert "cases unchanged" in body and "http" not in body
