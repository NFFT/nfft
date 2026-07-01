"""CLI: diff PR BMFs vs baseline, render the HTML report, write comment + check."""

from __future__ import annotations

import argparse
import json
import os

from tests.accuracy.diff import diff, load_bmf_tree
from tests.accuracy.htmlreport import render_report
from tests.accuracy.report import (
    check_summary,
    check_summary_no_baseline,
    comment_body,
    comment_body_no_baseline,
)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("pr_bmf_dir")
    ap.add_argument("base_bmf_dir")
    ap.add_argument("out_dir")
    ap.add_argument(
        "--report-url",
        help="published Pages URL of this report, embedded in the comment",
    )
    ap.add_argument("--gate", type=float, default=0.5)
    ap.add_argument(
        "--no-baseline",
        action="store_true",
        help="no develop baseline yet: absolute report only, no diff",
    )
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    pr_tree = load_bmf_tree(args.pr_bmf_dir)
    if args.no_baseline:
        doc = render_report(pr_tree, title="NFFT accuracy (this PR)")
        body = comment_body_no_baseline(args.report_url)
        conclusion, title, summary = check_summary_no_baseline()
    else:
        base_tree = load_bmf_tree(args.base_bmf_dir)
        result = diff(pr_tree, base_tree, gate=args.gate)
        doc = render_report(pr_tree, result, title="NFFT accuracy (this PR)")
        body = comment_body(result, args.report_url)
        conclusion, title, summary = check_summary(result)

    with open(os.path.join(args.out_dir, "index.html"), "w", encoding="utf-8") as f:
        f.write(doc)
    with open(os.path.join(args.out_dir, "comment.md"), "w", encoding="utf-8") as f:
        f.write(body)
    with open(os.path.join(args.out_dir, "check.json"), "w", encoding="utf-8") as f:
        json.dump({"conclusion": conclusion, "title": title, "summary": summary}, f)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
