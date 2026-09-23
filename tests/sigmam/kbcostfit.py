#!/usr/bin/env python3
"""Refit the cost constants in kbmodel.json from measured transform times.

usage: uv run --with numpy python kbcostfit.py [--write]

Input: data/cost-trafo-d.csv.gz, costprobe.c over data/cost-configs.txt:

    gcc -O2 -std=gnu99 -I<src>/include costprobe.c -L<lib> -lnfft3 -lfftw3 -lm
    ./a.out < data/cost-configs.txt | gzip > data/cost-trafo-d.csv.gz

Model, seconds per trafo (PRE_PHI_HUT | PRE_PSI, FFTW_ESTIMATE):

    t = fft * n log2 n + M_ * M + tap * M (2m + 2)

FFT, and the convolution's per-node and per-tap work. The O(N) deconvolution
cannot be told apart from the FFT, n being 1.25 to 2.5 N; the FFT term carries it. The FFT
term is deliberately smooth in n: power-of-two sizes, where FFTW's cost jumps
on this host, are left out of the fit: every point whose FFT alone takes more
than IRREGULAR times the median per n log2 n. That is not only powers of two
but any n with a large power-of-two factor, 5 * 2^k among them.
"""

import csv
import gzip
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = ("fft", "M", "tap")
IRREGULAR = 1.5


def features(N, M, n, m):
    return [n * math.log2(n), M, M * (2 * m + 2)]


def main():
    with gzip.open(os.path.join(HERE, "data", "cost-trafo-d.csv.gz"), "rt") as f:
        rows = list(csv.DictReader(f))
    per = [float(r["fft"]) / (int(r["n"]) * math.log2(int(r["n"]))) for r in rows]
    cut = IRREGULAR * float(np.median(per))
    X, y, skipped = [], [], set()
    for r, q in zip(rows, per):
        N, M, n, m = (int(r[k]) for k in ("N", "M", "n", "m"))
        if q > cut:
            skipped.add(n)
            continue
        X.append(features(N, M, n, m))
        y.append(float(r["trafo"]))
    print("irregular FFT sizes left out:", sorted(skipped))
    X, y = np.array(X), np.array(y)
    w = 1 / y  # relative least squares
    c, *_ = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)
    rel = y / (X @ c)
    cost = {k: float("%.3g" % v) for k, v in zip(NAMES, c)}
    print("points %d  measured/model p5 %.2f p50 %.2f p95 %.2f  (min %.2f max %.2f)"
          % (len(y), *np.percentile(rel, [5, 50, 95]), rel.min(), rel.max()))
    print("cost", cost, " tap/fft %.2f" % (cost["tap"] / cost["fft"]))
    if "--write" in sys.argv:
        p = os.path.join(HERE, "kbmodel.json")
        k = json.load(open(p))
        k["cost"] = cost
        with open(p, "w") as f:
            json.dump(k, f, indent=2)
            f.write("\n")
        print("wrote kbmodel.json")


if __name__ == "__main__":
    main()
