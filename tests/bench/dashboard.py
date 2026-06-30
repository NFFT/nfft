"""CLI: render the absolute accuracy dashboard + copy baseline BMFs (develop)."""
from __future__ import annotations

import argparse
import os
import shutil
import sys

from diff import load_bmf_tree
from heatmap import render_absolute

INDEX = """<!doctype html><meta charset=utf-8>
<title>NFFT accuracy dashboard</title>
<h1>NFFT accuracy (develop baseline)</h1>
<p>Worst-case accurate digits per case &times; testbed. Green = more digits.</p>
<img src="absolute.png" style="max-width:100%">
"""


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("bmf_dir")
    ap.add_argument("out_dir")
    args = ap.parse_args(argv)
    os.makedirs(os.path.join(args.out_dir, "baseline"), exist_ok=True)
    tree = load_bmf_tree(args.bmf_dir)
    render_absolute(tree, os.path.join(args.out_dir, "absolute.png"))
    with open(os.path.join(args.out_dir, "index.html"), "w", encoding="utf-8") as f:
        f.write(INDEX)
    for testbed in tree:
        shutil.copyfile(os.path.join(args.bmf_dir, f"{testbed}.bmf.json"),
                        os.path.join(args.out_dir, "baseline", f"{testbed}.bmf.json"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
