"""Fit the cost weight w in cost(n, m) = n*log2(n) + w*M*(2m+2).

Reads the CSVs written by costfit.c. Two answers, both per precision:

  fit      least squares on t ~ a*n*log2(n) + b*M*(2m+2) with relative
           residuals, giving w = b/a
  ranking  the w that orders the most pairs of candidates correctly within a
           group of equal (N, M, direction) -- the only thing the tuner ever
           asks of the model

Usage: python3 costfit.py <costfit-*.csv>...
"""
import csv
import math
import sys
from collections import defaultdict


def load(paths):
    rows = []
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh):
                rows.append({
                    "prec": r["prec"],
                    "N": int(r["N"]),
                    "M": int(r["M"]),
                    "n": int(r["n"]),
                    "m": int(r["m"]),
                    "dir": r["dir"],
                    "t": float(r["secs"]),
                })
    return [r for r in rows if r["t"] > 0.0]


def terms(r):
    return r["n"] * math.log2(r["n"]), r["M"] * (2 * r["m"] + 2)


def fit(rows):
    """Minimise sum (a*F/t + b*G/t - 1)^2, so every point weighs the same."""
    saa = sab = sbb = sa = sb = 0.0
    for r in rows:
        f, g = terms(r)
        f /= r["t"]
        g /= r["t"]
        saa += f * f
        sab += f * g
        sbb += g * g
        sa += f
        sb += g
    det = saa * sbb - sab * sab
    if det == 0.0:
        return None
    a = (sa * sbb - sb * sab) / det
    b = (sb * saa - sa * sab) / det
    return a, b


def rank_score(rows, w):
    """Fraction of same-group candidate pairs the model orders as measured.

    A group is one (precision, N, M, direction): the tuner only ever ranks
    candidates within one of those.
    """
    groups = defaultdict(list)
    for r in rows:
        groups[(r["prec"], r["N"], r["M"], r["dir"])].append(r)
    ok = tot = 0
    for g in groups.values():
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                ri, rj = g[i], g[j]
                if ri["t"] == rj["t"]:
                    continue
                fi, gi = terms(ri)
                fj, gj = terms(rj)
                ci = fi + w * gi
                cj = fj + w * gj
                if ci == cj:
                    continue
                tot += 1
                if (ci < cj) == (ri["t"] < rj["t"]):
                    ok += 1
    return ok / tot if tot else 0.0


def main(paths):
    rows = load(paths)
    precs = sorted({r["prec"] for r in rows})
    grid = [round(0.1 * i, 2) for i in range(1, 61)]

    print("%-12s %6s  %8s  %8s   %s" % ("precision", "rows", "w (fit)",
                                        "w (rank)", "rank score @ w"))
    for p in precs + ["ALL"]:
        sub = rows if p == "ALL" else [r for r in rows if r["prec"] == p]
        ab = fit(sub)
        w_fit = ab[1] / ab[0] if ab and ab[0] else float("nan")
        best_w, best_s = max(((w, rank_score(sub, w)) for w in grid),
                             key=lambda t: t[1])
        print("%-12s %6d  %8.3f  %8.2f   %.3f" %
              (p, len(sub), w_fit, best_w, best_s))

    print()
    print("rank score against w, all precisions pooled:")
    for w in [0.4, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 4.0]:
        print("  w = %-4s  %.3f" % (w, rank_score(rows, w)))

    print()
    print("per precision and direction:")
    for p in precs:
        for d in ("forward", "adjoint"):
            sub = [r for r in rows if r["prec"] == p and r["dir"] == d]
            if not sub:
                continue
            ab = fit(sub)
            w_fit = ab[1] / ab[0] if ab and ab[0] else float("nan")
            best_w, best_s = max(((w, rank_score(sub, w)) for w in grid),
                                 key=lambda t: t[1])
            print("  %-12s %-8s w_fit=%7.3f  w_rank=%5.2f  score=%.3f  "
                  "score@1.5=%.3f" %
                  (p, d, w_fit, best_w, best_s, rank_score(sub, 1.5)))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
