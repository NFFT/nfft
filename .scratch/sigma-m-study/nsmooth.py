#!/usr/bin/env python3
"""Compare picking n by smoothness against n = 2N and against the power of two.

Usage: nsmooth.py nsweep.csv
"""

import sys

from nanalyze import load


def smooth_rule(rows, limit):
    print(f"\nSmallest {limit}-smooth n >= 2N, versus n = 2N and the power of two:")
    print(
        f"{'N':>6} {'n_smooth':>9} {'sigma':>7} {'t':>9} {'vs 2N':>7} {'vs pow2':>8}"
    )
    for N in sorted(rows):
        v = {r["n"]: r for r in rows[N]}
        base = v[2 * N]
        cand = None
        for n in range(2 * N, 4 * N + 1):
            if n in v and v[n]["lpf"] <= limit:
                cand = v[n]
                break
        if not cand:
            continue
        p2 = [r for r in rows[N] if r["pow2"]]
        b2 = min(p2, key=lambda r: r["t_fwd"]) if p2 else None
        vs2 = base["t_fwd"] / cand["t_fwd"]
        vsp = (b2["t_fwd"] / cand["t_fwd"]) if b2 else float("nan")
        print(
            f"{N:>6} {cand['n']:>9} {cand['sigma']:>7.3f} "
            f"{cand['t_fwd']*1e6:>8.2f}u {vs2:>6.2f}x {vsp:>7.2f}x"
        )


if __name__ == "__main__":
    rows = load(sys.argv[1])
    smooth_rule(rows, 5)
    smooth_rule(rows, 7)
