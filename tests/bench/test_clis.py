import json
import os

import dashboard
import pr_report


def _write_bmf(path, **digit_by_name):
    with open(path, "w", encoding="utf-8") as f:
        json.dump({n: {"accuracy-digits": {"value": v}, "max-error": {"value": 10 ** -v}}
                   for n, v in digit_by_name.items()}, f)


def test_dashboard_emits_png_index_and_baseline(tmp_path):
    src = tmp_path / "bmf"; src.mkdir()
    _write_bmf(src / "tb1.bmf.json", **{"m/serial/file/fast/forward/1d/a": 13.5})
    out = tmp_path / "site"
    dashboard.main([str(src), str(out)])
    assert (out / "absolute.png").stat().st_size > 0
    assert (out / "index.html").exists()
    assert (out / "baseline" / "tb1.bmf.json").exists()


def test_pr_report_changed_writes_all_artifacts(tmp_path):
    pr = tmp_path / "pr"; pr.mkdir()
    base = tmp_path / "base"; base.mkdir()
    _write_bmf(pr / "tb1.bmf.json", **{"m/serial/file/fast/forward/1d/a": 15.0})
    _write_bmf(base / "tb1.bmf.json", **{"m/serial/file/fast/forward/1d/a": 13.0})
    out = tmp_path / "out"
    pr_report.main([str(pr), str(base), str(out),
                    "--abs-url", "http://x/a.png", "--rel-url", "http://x/r.png"])
    check = json.load(open(out / "check.json"))
    assert check["conclusion"] == "neutral" and "improved" in check["title"]
    body = (out / "comment.md").read_text()
    assert "Improvements" in body and "http://x/a.png" in body
    assert (out / "absolute.png").stat().st_size > 0
    assert (out / "relative.png").stat().st_size > 0


def test_pr_report_no_baseline(tmp_path):
    pr = tmp_path / "pr"; pr.mkdir()
    _write_bmf(pr / "tb1.bmf.json", **{"m/serial/file/fast/forward/1d/a": 13.5})
    out = tmp_path / "out"
    pr_report.main([str(pr), str(pr), str(out), "--no-baseline",
                    "--abs-url", "http://x/a.png"])
    body = (out / "comment.md").read_text()
    assert "No `develop` baseline" in body and "http://x/a.png" in body
    check = json.load(open(out / "check.json"))
    assert check["conclusion"] == "neutral" and "pending" in check["title"]
    assert (out / "absolute.png").stat().st_size > 0
    assert not (out / "relative.png").exists()  # no relative heatmap without baseline


def test_pr_report_flat_has_no_links_and_one_liner(tmp_path):
    pr = tmp_path / "pr"; pr.mkdir()
    base = tmp_path / "base"; base.mkdir()
    _write_bmf(pr / "tb1.bmf.json", **{"m/serial/file/fast/forward/1d/a": 13.5})
    _write_bmf(base / "tb1.bmf.json", **{"m/serial/file/fast/forward/1d/a": 13.5})
    out = tmp_path / "out"
    pr_report.main([str(pr), str(base), str(out)])
    body = (out / "comment.md").read_text()
    assert "cases unchanged" in body and "http" not in body
