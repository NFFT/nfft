#!/usr/bin/env python3
"""Summarise the fine-grained n sweep: is sigma > 2 worth it to reach a good n?

Usage: nanalyze.py nsweep.csv
"""

import csv
import sys
from collections import defaultdict


def load(path):
    rows = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            rows[int(r["N"])].append(
                dict(
                    n=int(r["n"]),
                    sigma=float(r["sigma"]),
                    lpf=int(r["lpf"]),
                    pow2=int(r["is_pow2"]),
                    err_fwd=float(r["err_fwd"]),
                    err_adj=float(r["err_adj"]),
                    t_fwd=float(r["t_fwd"]),
                    t_adj=float(r["t_adj"]),
                    t_pre=float(r["t_pre"]),
                )
            )
    return rows


def factorise(v):
    out, f = [], 2
    while f * f <= v:
        while v % f == 0:
            out.append(f)
            v //= f
        f += 1
    if v > 1:
        out.append(v)
    return out


def main():
    rows = load(sys.argv[1])
    print(
        f"{'N':>6} {'factors':>14} | {'t(n=2N)':>9} {'best t':>9} {'best n':>7} "
        f"{'sigma*':>7} {'speedup':>8} | {'pow2 n':>7} {'t(pow2)':>9} {'vs 2N':>7}"
    )
    for N in sorted(rows):
        v = sorted(rows[N], key=lambda r: r["n"])
        base = next(r for r in v if r["n"] == 2 * N)
        best = min(v, key=lambda r: r["t_fwd"])
        p2 = [r for r in v if r["pow2"]]
        fac = "*".join(str(f) for f in factorise(N))
        cells = (
            f"{N:>6} {fac:>14} | {base['t_fwd']*1e6:>8.2f}u "
            f"{best['t_fwd']*1e6:>8.2f}u {best['n']:>7} {best['sigma']:>7.3f} "
            f"{base['t_fwd']/best['t_fwd']:>7.2f}x |"
        )
        if p2:
            b2 = min(p2, key=lambda r: r["t_fwd"])
            cells += (
                f" {b2['n']:>7} {b2['t_fwd']*1e6:>8.2f}u "
                f"{base['t_fwd']/b2['t_fwd']:>6.2f}x"
            )
        else:
            cells += f" {'-':>7} {'-':>9} {'-':>7}"
        print(cells)

    print("\nTime by largest prime factor of n (median over all N, forward):")
    by = defaultdict(list)
    for N in rows:
        # normalise per N so different N are comparable
        base = next(r for r in rows[N] if r["n"] == 2 * N)["t_fwd"]
        for r in rows[N]:
            by[min(r["lpf"], 17)].append(r["t_fwd"] / base)
    print(f"{'lpf':>5} {'count':>6} {'median t / t(2N)':>18}")
    for lpf in sorted(by):
        v = sorted(by[lpf])
        lab = f"{lpf}" if lpf < 17 else ">=17"
        print(f"{lab:>5} {len(v):>6} {v[len(v)//2]:>18.2f}")

    print("\nAccuracy across the band (forward), at fixed m:")
    print(f"{'N':>6} {'err(2N)':>11} {'err(4N)':>11} {'ratio':>8}")
    for N in sorted(rows):
        v = sorted(rows[N], key=lambda r: r["n"])
        a = next(r for r in v if r["n"] == 2 * N)["err_fwd"]
        b = v[-1]["err_fwd"]
        print(f"{N:>6} {a:>11.3e} {b:>11.3e} {a/b:>8.2f}")


if __name__ == "__main__":
    main()
