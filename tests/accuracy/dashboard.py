"""CLI: render the absolute accuracy dashboard (HTML) + copy baseline BMFs.

Run on develop pushes; the output is published to the gh-pages branch root.
"""

from __future__ import annotations

import argparse
import os
import shutil

from tests.accuracy.diff import load_bmf_tree
from tests.accuracy.htmlreport import render_report


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("bmf_dir")
    ap.add_argument("out_dir")
    args = ap.parse_args(argv)
    os.makedirs(os.path.join(args.out_dir, "baseline"), exist_ok=True)
    tree = load_bmf_tree(args.bmf_dir)
    doc = render_report(tree, title="NFFT accuracy (develop baseline)")
    with open(os.path.join(args.out_dir, "index.html"), "w", encoding="utf-8") as f:
        f.write(doc)
    for testbed in tree:
        shutil.copyfile(
            os.path.join(args.bmf_dir, f"{testbed}.bmf.json"),
            os.path.join(args.out_dir, "baseline", f"{testbed}.bmf.json"),
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
