#!/usr/bin/env python3
"""Summarise hostcheck results: helper grid against the legacy grid, same m.

usage: python3 summarize.py results.csv [more.csv ...]

A case is a win when the legacy arm takes more than 2 % longer (trafo plus
adjoint), a loss when the helper does, a tie otherwise or when both arms
chose the same n. Any helper error above the goal is a failure of the error
model and is listed first.
"""

import csv
import sys

TIE = 0.02
SMALL_GRID = 4096  # n of a 64 KB complex double grid


def load(paths):
    info, pairs = [], {}
    for path in paths:
        rows = []
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    info.append(line[1:].strip())
                else:
                    rows.append(line)
        for r in csv.DictReader(rows):
            key = (int(r["N"]), int(r["M"]), float(r["goal"]))
            pairs.setdefault(key, {})[r["arm"]] = r
    return info, {k: v for k, v in pairs.items() if len(v) == 2}


def secs(r):
    return float(r["t_trafo"]) + float(r["t_adjoint"])


def classify(h, l):
    s = secs(l) / secs(h)
    if h["n"] == l["n"] or abs(s - 1) <= TIE:
        return "tie", s
    return ("win" if s > 1 else "loss"), s


def table(title, groups):
    print("\n%s" % title)
    print("  %-26s %5s %5s %5s %5s   %s" % ("", "cases", "win", "tie", "loss", "speedup median / max"))
    for name, cases in groups:
        if not cases:
            continue
        c = {"win": 0, "tie": 0, "loss": 0}
        s = sorted(x[1] for x in cases)
        for k, _ in cases:
            c[k] += 1
        print("  %-26s %5d %5d %5d %5d   %.2f / %.2f" % (name, len(cases), c["win"], c["tie"], c["loss"], s[len(s) // 2], s[-1]))


def main():
    info, pairs = load(sys.argv[1:])
    for line in info:
        print(line)
    rows = []
    for (N, M, g), a in sorted(pairs.items()):
        h, l = a["helper"], a["legacy"]
        k, s = classify(h, l)
        rows.append((N, M, g, h, l, k, s))

    miss = [r for r in rows if float(r[3]["err"]) > r[2]]
    print("\nhelper error above the goal: %d of %d cases" % (len(miss), len(rows)))
    for N, M, g, h, l, k, s in miss:
        print("  N=%d M=%d goal %.0e: err %s, n=%s m=%s" % (N, M, g, h["err"], h["n"], h["m"]))

    table("all cases", [("all", [(r[5], r[6]) for r in rows])])
    table("by helper grid", [
        ("small (n <= %d)" % SMALL_GRID, [(r[5], r[6]) for r in rows if int(r[3]["n"]) <= SMALL_GRID]),
        ("large (n > %d)" % SMALL_GRID, [(r[5], r[6]) for r in rows if int(r[3]["n"]) > SMALL_GRID])])
    table("by nodes", [
        ("M < N/4", [(r[5], r[6]) for r in rows if r[1] < r[0] / 4]),
        ("N/4 <= M <= 8N", [(r[5], r[6]) for r in rows if r[0] / 4 <= r[1] <= 8 * r[0]]),
        ("M > 8N", [(r[5], r[6]) for r in rows if r[1] > 8 * r[0]])])
    table("by goal", [("goal %.0e" % g, [(r[5], r[6]) for r in rows if r[2] == g])
                      for g in sorted({r[2] for r in rows}, reverse=True)])

    losses = sorted((r for r in rows if r[5] == "loss"), key=lambda r: r[6])
    print("\nlosses, slowest first (speedup = legacy time / helper time):")
    for N, M, g, h, l, k, s in losses:
        print("  %.2f  N=%-6d M=%-7d goal %.0e  n %s vs %s, m=%s, helper %.1f us" % (
            s, N, M, g, h["n"], l["n"], h["m"], 1e6 * secs(h)))
    if not losses:
        print("  none")


if __name__ == "__main__":
    main()
