#!/usr/bin/env python3
"""Choose the error-model form by validation, not by in-sample fit quality.

Candidate forms, all keeping the two exact rates D and A:

    T(m,s)   = a * u^p * m^r * exp(-D(s)*m)          u = 1 - 1/sigma
    U(m,s,e) = c * e * u^q * m^w * exp(alpha*A(s)*m)

with r, q, w and alpha optionally pinned to 0 / 1. Each candidate is fitted,
its (a,c) raised to an upper envelope over every measurement, then scored on
what a tuner actually has to get right:

    misses      the returned m fails the goal in measurement   (must be 0)
    exact       the returned m equals the smallest sufficient m
    tail        the returned m is more than 2 above it         (must be ~0)
    unreach     goals declared unattainable that measurement reaches
    floor       median pessimism of the reported attainable accuracy

Usage: select.py sweep-*.csv
"""

import math
import sys
from collections import defaultdict

from fit import D, A, EPS, load, lstsq, split_regimes

MARGIN = 1.25


def make_model(f):
    def model(rec):
        s, m = rec["sigma"], rec["m"]
        u = 1.0 - 1.0 / s
        t = f["a"] * u ** f["p"] * m ** f["r"] * math.exp(-D(s) * m)
        uu = (
            f["c"]
            * rec["eps"]
            * u ** f["q"]
            * m ** f["w"]
            * math.exp(f["alpha"] * A(s) * m)
        )
        return t + uu

    return model


def fit_form(pts, key, use_r, use_qw, use_alpha, use_w=True):
    trunc, round_ = split_regimes(pts, key)

    cols = [lambda r: 1.0, lambda r: math.log(1.0 - 1.0 / r["sigma"])]
    if use_r:
        cols.append(lambda r: math.log(r["m"]))
    rowsA = [[c(r) for c in cols] for r in trunc]
    rhs = [math.log(r[key]) + D(r["sigma"]) * r["m"] for r in trunc]
    sol = lstsq(rowsA, rhs)
    a, p = math.exp(sol[0]), sol[1]
    rr = sol[2] if use_r else 0.0

    cols = [lambda r: 1.0]
    if use_alpha:
        cols.append(lambda r: A(r["sigma"]) * r["m"])
    if use_qw:
        cols.append(lambda r: math.log(1.0 - 1.0 / r["sigma"]))
        if use_w:
            cols.append(lambda r: math.log(r["m"]))
    rowsA = [[c(r) for c in cols] for r in round_]
    rhs = [
        math.log(r[key] / r["eps"]) - (0.0 if use_alpha else A(r["sigma"]) * r["m"])
        for r in round_
    ]
    sol = lstsq(rowsA, rhs)
    c = math.exp(sol[0])
    i = 1
    alpha = 1.0
    if use_alpha:
        alpha = sol[i]
        i += 1
    q = w = 0.0
    if use_qw:
        q = sol[i]
        if use_w:
            w = sol[i + 1]

    f = dict(a=a, p=p, r=rr, c=c, alpha=alpha, q=q, w=w)

    # Raise (a,c) to the tightest dominating pair.
    live = [r for r in pts if r[key] > 0.0]
    best = None
    for step in range(-40, 81):
        cc = c * (1.15**step)
        g = dict(f, a=1.0, c=cc)
        mg = make_model(g)
        need = 0.0
        for rec in live:
            u = 1.0 - 1.0 / rec["sigma"]
            t = u ** f["p"] * rec["m"] ** f["r"] * math.exp(-D(rec["sigma"]) * rec["m"])
            if t <= 0:
                continue
            need = max(need, (rec[key] - (mg(rec) - t)) / t)
        aa = max(need, 1e-12)
        h = dict(f, a=aa, c=cc)
        mh = make_model(h)
        slack = sum(math.log(mh(rec) / rec[key]) for rec in live)
        if best is None or slack < best[0]:
            best = (slack, aa, cc)
    f["a"] = best[1] * MARGIN
    f["c"] = best[2] * MARGIN
    return f


def score(pts, key, f):
    model = make_model(f)
    over = sum(
        1 for r in pts if r[key] > 0 and r[key] > model(r)
    )
    series = defaultdict(list)
    for r in pts:
        if r[key] > 0.0:
            series[(r["prec"], r["sigma"], r["N"])].append(r)
    goals = [10.0**-e for e in range(2, 18)]
    misses = exact = tail = unreach = total = 0
    floor_ratio = []
    for (prec, s, N), v in series.items():
        v.sort(key=lambda r: r["m"])
        m_max = v[-1]["m"]
        eps = EPS[prec]
        meas = {r["m"]: r[key] for r in v}
        curve = [
            (m, model(dict(sigma=s, m=m, eps=eps, N=N))) for m in range(1, m_max + 1)
        ]
        floor_ratio.append(min(e for _, e in curve) / min(meas.values()))
        for goal in goals:
            mm = next((m for m, e in curve if e <= goal), None)
            ok = [m for m in sorted(meas) if meas[m] <= goal]
            if mm is None:
                if ok:
                    unreach += 1
                continue
            total += 1
            if meas.get(mm, float("inf")) > goal:
                misses += 1
            if ok:
                d = mm - ok[0]
                if d == 0:
                    exact += 1
                elif d > 2:
                    tail += 1
    floor_ratio.sort()
    return dict(
        over=over,
        misses=misses,
        exact=100.0 * exact / max(total, 1),
        tail=100.0 * tail / max(total, 1),
        unreach=unreach,
        floor=floor_ratio[len(floor_ratio) // 2],
        total=total,
    )


FORMS = [
    ("T:u^p,m^r  U:c,alpha        ", True, False, True, False),
    ("T:u^p,m^r  U:c,alpha,u^q    ", True, True, True, False),
    ("T:u^p,m^r  U:c,alpha,u^q,m^w", True, True, True, True),
]


def main():
    pts = load(sys.argv[1:])
    print(f"{len(pts)} aggregated points\n")
    for key in ("fwd", "adj"):
        print(f"=== {key} ===")
        print(
            f"{'form':<34} {'over':>5} {'miss':>5} {'exact%':>7} "
            f"{'tail%':>6} {'unreach':>8} {'floor':>7}"
        )
        for name, ur, uqw, ua, uw in FORMS:
            f = fit_form(pts, key, ur, uqw, ua, uw)
            s = score(pts, key, f)
            print(
                f"{name:<34} {s['over']:>5} {s['misses']:>5} {s['exact']:>7.0f} "
                f"{s['tail']:>6.1f} {s['unreach']:>8} {s['floor']:>7.0f}"
            )
        print()


if __name__ == "__main__":
    main()
