#!/usr/bin/env python3
"""Measured attainable floor per (sigma, precision), expressed in units of eps.

If the floor is a fixed multiple of eps, oversampling has bought back the whole
mantissa and the transform reaches the machine-epsilon region. If it grows with
the precision, the eps**(D/b) law is still biting.

Usage: floors.py sweep-*.csv
"""

import math
import sys
from collections import defaultdict

from fit import D, A, EPS, load

ORDER = ["float", "double", "long_double"]


def main():
    pts = load(sys.argv[1:])
    best = defaultdict(lambda: float("inf"))
    for r in pts:
        for key in ("fwd", "adj"):
            if r[key] > 0:
                k = (key, round(r["sigma"], 3), r["prec"], r["N"])
                best[k] = min(best[k], r[key])
    # worst N per (key, sigma, prec)
    agg = defaultdict(lambda: 0.0)
    for (key, s, prec, N), v in best.items():
        agg[(key, s, prec)] = max(agg[(key, s, prec)], v)

    for key in ("fwd", "adj"):
        print(f"\n=== {key}: measured floor / eps ===")
        print(
            f"{'sigma':>7} {'gamma':>7} | "
            + " ".join(f"{p:>13}" for p in ORDER)
            + "   verdict"
        )
        sigmas = sorted({s for (k, s, p) in agg if k == key})
        for s in sigmas:
            d = D(s)
            b = d + A(s)
            g = d / b
            vals = []
            for p in ORDER:
                v = agg.get((key, s, p))
                vals.append(v / EPS[p] if v else None)
            cells = " ".join(
                f"{v:13.1f}" if v is not None else f"{'-':>13}" for v in vals
            )
            known = [v for v in vals if v is not None]
            if len(known) >= 2:
                spread = max(known) / min(known)
                verdict = "flat (eps region)" if spread < 8 else f"grows x{spread:.0f}"
            else:
                verdict = ""
            print(f"{s:7.3f} {g:7.4f} | {cells}   {verdict}")


if __name__ == "__main__":
    main()
