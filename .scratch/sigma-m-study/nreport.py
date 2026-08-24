#!/usr/bin/env python3
"""Emit the n-sweep data as compact JSON for the HTML summary.

Usage: nreport.py nsweep.csv > nsweep.json
"""

import json
import sys

from nanalyze import load, factorise


def bucket(lpf):
    if lpf <= 2:
        return 0
    if lpf <= 5:
        return 1
    if lpf <= 13:
        return 2
    return 3


def smallest_smooth(rowmap, N, limit):
    for n in range(2 * N, 4 * N + 1):
        r = rowmap.get(n)
        if r and r["lpf"] <= limit:
            return r
    return None


def next_power_of_2(x):
    p = 1
    while p < x:
        p *= 2
    return p


def main():
    rows = load(sys.argv[1])
    out = {"panels": [], "lpf": {}, "summary": []}

    rel_by_bucket = {0: [], 1: [], 2: [], 3: []}

    for N in sorted(rows):
        v = sorted(rows[N], key=lambda r: r["n"])
        rowmap = {r["n"]: r for r in v}
        base = rowmap[2 * N]
        best = min(v, key=lambda r: r["t_fwd"])
        smooth5 = smallest_smooth(rowmap, N, 5)
        legacy_n = 2 * next_power_of_2(N)
        legacy = rowmap.get(legacy_n)

        for r in v:
            rel_by_bucket[bucket(r["lpf"])].append(r["t_fwd"] / base["t_fwd"])

        out["panels"].append(
            {
                "N": N,
                "factors": "x".join(str(f) for f in factorise(N)),
                "n": [r["n"] for r in v],
                "sigma": [round(r["sigma"], 4) for r in v],
                "t": [round(r["t_fwd"] * 1e6, 3) for r in v],
                "err": [float(f"{r['err_fwd']:.3e}") for r in v],
                "b": [bucket(r["lpf"]) for r in v],
                "lpf": [r["lpf"] for r in v],
                "mark": {
                    "base": base["n"],
                    "best": best["n"],
                    "smooth5": smooth5["n"] if smooth5 else None,
                    "legacy": legacy_n if legacy else None,
                },
            }
        )

        out["summary"].append(
            {
                "N": N,
                "factors": "x".join(str(f) for f in factorise(N)),
                "base_n": base["n"],
                "base_t": round(base["t_fwd"] * 1e6, 2),
                "legacy_n": legacy_n,
                "legacy_sigma": round(legacy_n / N, 3) if legacy else None,
                "legacy_t": round(legacy["t_fwd"] * 1e6, 2) if legacy else None,
                "smooth_n": smooth5["n"] if smooth5 else None,
                "smooth_sigma": round(smooth5["sigma"], 3) if smooth5 else None,
                "smooth_t": round(smooth5["t_fwd"] * 1e6, 2) if smooth5 else None,
                "best_n": best["n"],
                "best_sigma": round(best["sigma"], 3),
                "best_t": round(best["t_fwd"] * 1e6, 2),
                "err_2N": float(f"{base['err_fwd']:.3e}"),
                "err_4N": float(f"{v[-1]['err_fwd']:.3e}"),
            }
        )

    labels = ["2 (power of two)", "3 or 5", "7 to 13", "17 and above"]
    for b in range(4):
        s = sorted(rel_by_bucket[b])
        out["lpf"][b] = {
            "label": labels[b],
            "count": len(s),
            "median": round(s[len(s) // 2], 3),
            "p90": round(s[int(0.9 * len(s))], 3),
        }

    json.dump(out, sys.stdout, separators=(",", ":"))


if __name__ == "__main__":
    main()
