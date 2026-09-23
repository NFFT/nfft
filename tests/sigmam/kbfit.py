#!/usr/bin/env python3
"""Refit kbmodel.json from the worst-case sweeps and validate both building blocks.

usage: uv run --with scipy python kbfit.py [--write] [margin]

Inputs (run_worst.sh regenerates them):
  data/trunc-exact-width.csv.gz        exact-arithmetic truncation, max over theta, u
  data/floor-worst-width-{f,d,l}.csv.gz measured worst case, one row per (sigma, N, m)

The whole model carries the margin: model * (1 - margin) stays above every
measured point from m = 1 up to the model's best m, the stretch smallest_m can
return. The model is linear in (CT, C1, C0), so the fit is a linear program;
the envelope of T * I0(w D) / R over the exact scan, divided by (1 - margin), is
a lower bound on CT. G, the trough level of the ripple term, is not linear
in the model; it is scanned, and the value giving the most exact m is kept. Without --write the fit is only printed.
"""

import csv
import gzip
import json
import math
import os
import sys

import kbmodel as km

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SWEEPS = {"float": "f", "double": "d", "quad": "l"}
QUAD_TRUNC_LIMIT = 1e-19  # below this the quad scan reads its own round-off


def rows(name):
    with gzip.open(os.path.join(DATA, name), "rt") as f:
        yield from csv.DictReader(f)


def fit_ct(k):
    ct = 0.0
    for r in rows("trunc-exact-width.csv.gz"):
        s, m, t = float(r["sigma"]), int(r["m"]), float(r["trunc"])
        if t > QUAD_TRUNC_LIMIT:
            D = 2 * math.pi * math.sqrt(1 - 1 / s)
            ct = max(ct, t * math.exp(km.log_i0((m + 1) * D)) / km.ripple(s, m + 1, k))
    return ct


def load_sweeps():
    """{(type, sigma, N): {m: err}}, cut where the grid caps m."""
    out, cap = {}, {}
    for typ, p in SWEEPS.items():
        path = "floor-worst-width-%s.csv.gz" % p
        if not os.path.exists(os.path.join(DATA, path)):
            continue
        for r in rows(path):
            e = float(r["err"])
            if not (e > 0 and math.isfinite(e)):
                continue
            key = (typ, float(r["sigma"]), int(r["N"]))
            out.setdefault(key, {})[int(r["m"])] = e
            cap[key] = min(int(r["N"]) - 1, (int(r["n"]) - 3) // 2)
        for key in list(out):
            # a curve that never turns up inside the swept range has no floor
            ms = out[key]
            if min(ms, key=ms.get) >= min(cap[key], max(ms)):
                del out[key]
    return out


def fit_lp(sw, g, margin):
    """Tightest (CT, C1, C0) keeping model * (1 - margin) above every point up
    to the model's best m. The model is linear in the constants, so this is a
    linear program; its constraint set depends on best m, hence the loop."""
    from scipy.optimize import linprog

    k = dict(km.K, G=g)
    ct_lo = fit_ct(k) / (1 - margin)
    for _ in range(20):
        A, rhs, obj = [], [], [0.0, 0.0, 0.0]
        for (typ, s, N), ms in sw.items():
            eps = km.EPS[typ]
            b = math.pi * (2 - 1 / s)
            D = 2 * math.pi * math.sqrt(1 - 1 / s)
            mb, _ = km.best_m(s, eps, k)
            for m, e in ms.items():
                if m > mb:
                    continue
                ld = km.log_i0((m + 1) * D)
                row = [km.ripple(s, m + 1, k) * math.exp(-ld) / e, eps * math.exp(km.log_i0((m + 1) * b) - ld) / e, eps / e]
                A.append([-v for v in row])
                rhs.append(-1 / (1 - margin))
                obj = [o + v for o, v in zip(obj, row)]
        res = linprog(obj, A_ub=A, b_ub=rhs, bounds=[(ct_lo, None), (0, None), (0, None)])
        assert res.success, res.message
        new = dict(zip(("CT", "C1", "C0"), res.x), G=g)
        if all(abs(new[n] - k[n]) <= 1e-9 * max(1, k[n]) for n in new):
            break
        k = new
    return new


def block2(sw, k):
    """{type: ({m_model - m_measured or "MISS": count}, [misses])}"""
    out = {}
    for typ in SWEEPS:
        hist, worst = {}, []
        for (t, s, N), ms in sw.items():
            if t != typ:
                continue
            eps = km.EPS[typ]
            fl = km.floor_error(s, eps, k)
            for j in range(0, 400):
                g = 10 ** (-1 - j / 8)
                if g < fl:
                    break
                m = km.smallest_m(s, g, eps, k)
                mm = min((q for q, e in ms.items() if e <= g), default=None)
                if mm is None:
                    continue
                d = m - mm
                if ms.get(m, math.inf) > g:
                    d = "MISS"
                    worst.append((s, N, g, m))
                hist[d] = hist.get(d, 0) + 1
        if hist:
            out[typ] = (hist, worst)
    return out


def exact_share(sw, k):
    tot = ex = 0
    for hist, worst in block2(sw, k).values():
        if worst:
            return -1.0
        tot += sum(hist.values())
        ex += hist.get(0, 0)
    return ex / tot


def validate(sw, k):
    print("\nbuilding block 2, per (type, sigma, N) curve and goal:")
    print("  model m minus measured smallest m meeting the goal")
    for typ, (hist, worst) in block2(sw, k).items():
        if hist:
            tot = sum(hist.values())
            print("  %-6s %s" % (typ, "  ".join(
                "%s: %d (%.0f%%)" % (d, n, 100 * n / tot)
                for d, n in sorted(hist.items(), key=lambda x: str(x[0])))))
            for w in worst[:5]:
                print("         miss sigma %.4f N %d goal %.1e m %d" % w)

    print("\nbuilding block 1: worst measured floor over N, per requested sigma")
    for typ in SWEEPS:
        eps = km.EPS[typ]
        fl = {}
        for (t, s, N), ms in sw.items():
            if t == typ:
                sr = round(s * 8) / 8 if abs(s * 8 - round(s * 8)) < 0.05 else round(s, 2)
                fl[sr] = max(fl.get(sr, 0.0), min(ms.values()))
        if not fl:
            continue
        fails = 0
        for j in range(0, 400):
            g = 10 ** (-1 - j / 4)
            if g < km.floor_error(2.5, eps, k) / 100:
                break
            s = km.smallest_sigma(g, eps, k=k)
            gg = max(g, km.floor_error(2.5, eps, k))
            fails += sum(1 for sr, e in fl.items() if sr >= s - 1e-3 and e > gg)
        r = [e / km.floor_error(sr, eps, k) for sr, e in fl.items()]
        print("  %-6s measured/model floor %.2f .. %.2f   failures %d" % (typ, min(r), max(r), fails))


def main():
    margin = float(next((a for a in sys.argv[1:] if a[0].isdigit()), 0.05))
    sw = load_sweeps()
    # G enters R non-linearly: scan it, the LP fits the rest
    fits = [fit_lp(sw, 1.0 + 0.05 * i, margin) for i in range(9)]
    for f in fits:
        print("  G %.2f  exact m %.1f%%" % (f["G"], 100 * exact_share(sw, f)))
    k = max(fits, key=lambda f: exact_share(sw, f))
    # rounding up keeps the envelope
    k = {n: math.ceil(v * 1000) / 1000 for n, v in k.items()}
    print("curves %d, margin %.0f%%" % (len(sw), 100 * margin))
    print("fit     ", k)
    print("current ", km.K)
    validate(sw, k)
    if "--write" in sys.argv:
        path = os.path.join(HERE, "kbmodel.json")
        with open(path) as f:
            out = json.load(f)
        out.update(k)  # cost and choice belong to kbcostfit.py
        with open(path, "w") as f:
            json.dump(out, f, indent=2)
            f.write("\n")
        print("\nwrote kbmodel.json")


if __name__ == "__main__":
    main()
