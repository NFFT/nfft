#!/usr/bin/env python3
"""The tie map behind the regime bounds in kbmodel.json "choice".

Per type, goal, N and M/N: how many distinct m lie within the tie tolerance of
the cheapest pair, and a * where the cheapest pair sits more than 0.1 in sigma
above the smallest-sigma pair within the tolerance. The balanced band is where
the * appear; outside it one pair wins or all near-ties share the smallest
sigma. Set "balanced" to enclose the band.

usage: python3 kbties.py
"""

import kbmodel as km

RATIOS = [2 ** (e / 2) for e in range(-8, 11)]
CASES = (("double", (1e-5, 1e-9, 1e-13)), ("float", (1e-3, 1e-5)),
         ("quad", (1e-12, 1e-25)))


def main():
    tie = km.K["choice"]["tie"]
    print("M/N" + " " * 22 + " ".join("%5.3g" % r for r in RATIOS))
    for typ, goals in CASES:
        eps = km.EPS[typ]
        for g in goals:
            for N in (1024, 65536):
                cells = []
                for r in RATIOS:
                    M = max(1, round(r * N))
                    c = km.candidates(N, M, g, eps)
                    b = min(c, key=lambda x: x[2])
                    near = [x for x in c if x[2] <= (1 + tie) * b[2]]
                    smin = min(x[0] for x in near) / N
                    cells.append("%d%s" % (len({x[1] for x in near}),
                                           "*" if b[0] / N - smin > 0.1 else " "))
                print("%-6s %.0e N=%-6d  " % (typ, g, N) + " ".join("%5s" % x for x in cells))
    print("\ncurrent balanced band:", km.K["choice"]["balanced"])


if __name__ == "__main__":
    main()
