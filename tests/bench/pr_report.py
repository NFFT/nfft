"""CLI: diff PR BMFs vs baseline, render heatmaps, write comment.md + check.json."""
from __future__ import annotations

import argparse
import json
import os
import sys

from diff import diff, load_bmf_tree
from heatmap import emoji_grid, render_absolute, render_relative
from report import check_summary, comment_body


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("pr_bmf_dir")
    ap.add_argument("base_bmf_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--abs-url")
    ap.add_argument("--rel-url")
    ap.add_argument("--gate", type=float, default=0.5)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    pr_tree = load_bmf_tree(args.pr_bmf_dir)
    base_tree = load_bmf_tree(args.base_bmf_dir)
    result = diff(pr_tree, base_tree, gate=args.gate)

    render_absolute(pr_tree, os.path.join(args.out_dir, "absolute.png"))
    render_relative(result, os.path.join(args.out_dir, "relative.png"))

    png_urls = None
    if args.abs_url and args.rel_url:
        png_urls = {"absolute": args.abs_url, "relative": args.rel_url}

    # The emoji grid is embedded in the comment ABOVE the itemized groups.
    body = comment_body(result, png_urls)
    if result.improvements or result.regressions:
        body = body.replace("## Accuracy report\n",
                            "## Accuracy report\n\n" + emoji_grid(result) + "\n\n", 1)
    with open(os.path.join(args.out_dir, "comment.md"), "w", encoding="utf-8") as f:
        f.write(body)

    conclusion, title, summary = check_summary(result)
    with open(os.path.join(args.out_dir, "check.json"), "w", encoding="utf-8") as f:
        json.dump({"conclusion": conclusion, "title": title, "summary": summary}, f)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
