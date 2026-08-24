#!/usr/bin/env python3
"""Fit the Kaiser-Bessel error model to the sigma/m sweep.

The measured error of the 1-D NFFT with the Kaiser-Bessel window splits into
a window-truncation term that falls with m and a roundoff term that the
deconvolution step amplifies with m:

    D(s) = 2*pi*sqrt(1 - 1/s)                 truncation decay rate
    A(s) = 2*pi*(1 - 1/(2s)) - D(s)           amplification rate

    T(m,s)   = a * (1-1/s)**p * m**r * exp(-D(s)*m)
    U(m,s,e) = c * e * (1-1/s)**q * m**w * exp(alpha*A(s)*m)
    E(m,s,e) = T + U

Both rates come from the window's Fourier transform and are used as given;
only the prefactors a, p, r, c, alpha are fitted. a and c are then raised to
the smallest pair that dominates every measurement (upper envelope), chosen
to keep the envelope as tight as possible.

Usage: fit.py sweep-*.csv
"""

import csv
import math
import sys
from collections import defaultdict

EPS = {"float": 2.0**-23, "double": 2.0**-52, "long_double": 2.0**-112}

# Include an N power in the roundoff branch (the 1/||f_hat||_1 normalisation).
USE_N = False
N_REF = 1024.0
DIRECTIONS = ("fwd", "adj")


def D(s):
    return 2.0 * math.pi * math.sqrt(1.0 - 1.0 / s)


def A(s):
    return 2.0 * math.pi * (1.0 - 1.0 / (2.0 * s)) - D(s)


def lstsq(rowsA, rhs):
    """Solve the normal equations for a small dense least-squares problem."""
    k = len(rowsA[0])
    ata = [[sum(r[i] * r[j] for r in rowsA) for j in range(k)] for i in range(k)]
    atb = [sum(r[i] * y for r, y in zip(rowsA, rhs)) for i in range(k)]
    for i in range(k):  # Gauss-Jordan
        piv = max(range(i, k), key=lambda t: abs(ata[t][i]))
        ata[i], ata[piv] = ata[piv], ata[i]
        atb[i], atb[piv] = atb[piv], atb[i]
        d = ata[i][i]
        if abs(d) < 1e-300:
            return None
        for j in range(i, k):
            ata[i][j] /= d
        atb[i] /= d
        for t in range(k):
            if t == i:
                continue
            f = ata[t][i]
            for j in range(i, k):
                ata[t][j] -= f * ata[i][j]
            atb[t] -= f * atb[i]
    return atb


def load(paths):
    """Aggregate the sweep to the max error over trials per (prec,s,N,m)."""
    agg = defaultdict(lambda: [0.0, 0.0])
    for path in paths:
        with open(path) as fh:
            for r in csv.DictReader(fh):
                key = (r["prec"], float(r["sigma"]), int(r["N"]), int(r["m"]))
                cur = agg[key]
                cur[0] = max(cur[0], float(r["err_fwd"]))
                cur[1] = max(cur[1], float(r["err_adj"]))
    pts = []
    for (prec, s, N, m), (ef, ea) in agg.items():
        pts.append(dict(prec=prec, sigma=s, N=N, m=m, eps=EPS[prec], fwd=ef, adj=ea))
    return pts


def split_regimes(pts, key):
    """Per (prec,sigma,N) series, label points truncation- or roundoff-side."""
    series = defaultdict(list)
    for r in pts:
        if r[key] > 0.0:
            series[(r["prec"], r["sigma"], r["N"])].append(r)
    trunc, round_ = [], []
    for _, v in series.items():
        v.sort(key=lambda r: r["m"])
        i_min = min(range(len(v)), key=lambda i: v[i][key])
        e_min = v[i_min][key]
        for i, r in enumerate(v):
            if i <= i_min - 2 and r[key] > 100.0 * e_min:
                trunc.append(r)
            elif i >= i_min + 2:
                round_.append(r)
    return trunc, round_


def fit_direction(pts, key, margin):
    trunc, round_ = split_regimes(pts, key)

    # --- truncation branch: log E + D*m = log a + p*log(1-1/s) + r*log m ---
    rowsA = [
        [1.0, math.log(1.0 - 1.0 / r["sigma"]), math.log(r["m"])] for r in trunc
    ]
    rhs = [math.log(r[key]) + D(r["sigma"]) * r["m"] for r in trunc]
    sol = lstsq(rowsA, rhs)
    a, p, rr = math.exp(sol[0]), sol[1], sol[2]

    # --- roundoff branch: log(E/eps) = log c + alpha*A*m (+ g*log N) ------
    rowsA = [
        [
            1.0,
            A(r["sigma"]) * r["m"],
            math.log(1.0 - 1.0 / r["sigma"]),
            math.log(r["m"]),
        ]
        for r in round_
    ]
    rhs = [math.log(r[key] / r["eps"]) for r in round_]
    sol = lstsq(rowsA, rhs)
    c, alpha, q, w = math.exp(sol[0]), sol[1], sol[2], sol[3]
    gN = 0.0

    def raw(rec, aa, cc):
        s, m = rec["sigma"], rec["m"]
        t = aa * (1.0 - 1.0 / s) ** p * m**rr * math.exp(-D(s) * m)
        u = (
            cc
            * rec["eps"]
            * (1.0 - 1.0 / s) ** q
            * m**w
            * math.exp(alpha * A(s) * m)
        )
        return t, u

    # --- raise (a,c) to the tightest dominating pair ----------------------
    live = [r for r in pts if r[key] > 0.0]
    best = None
    for step in range(-30, 61):
        cc = c * (1.15**step)
        need = 0.0
        for rec in live:
            t, u = raw(rec, 1.0, cc)
            if t <= 0:
                continue
            need = max(need, (rec[key] - u) / t)
        aa = max(need, 1e-12)
        slack = 0.0
        for rec in live:
            t, u = raw(rec, aa, cc)
            slack += math.log((t + u) / rec[key])
        if best is None or slack < best[0]:
            best = (slack, aa, cc)
    _, a_env, c_env = best
    return dict(
        a=a_env * margin,
        p=p,
        r=rr,
        c=c_env * margin,
        alpha=alpha,
        q=q,
        w=w,
        gN=gN,
        n_trunc=len(trunc),
        n_round=len(round_),
        a_ls=a,
        c_ls=c,
    )


def model(rec, f):
    s, m = rec["sigma"], rec["m"]
    t = f["a"] * (1.0 - 1.0 / s) ** f["p"] * m ** f["r"] * math.exp(-D(s) * m)
    u = (
        f["c"]
        * rec["eps"]
        * (1.0 - 1.0 / s) ** f["q"]
        * m ** f["w"]
        * math.exp(f["alpha"] * A(s) * m)
    )
    return t + u


def report(pts, key, f):
    ratios = []
    worst = (0.0, None)
    for rec in pts:
        if rec[key] <= 0:
            continue
        q = rec[key] / model(rec, f)
        ratios.append(q)
        if q > worst[0]:
            worst = (q, rec)
    ratios.sort()
    over = sum(1 for q in ratios if q > 1.0)
    med = ratios[len(ratios) // 2]
    print(
        f"  a={f['a']:.5g} p={f['p']:.4f} r={f['r']:.4f}"
    )
    print(
        f"  c={f['c']:.5g} alpha={f['alpha']:.4f} q={f['q']:.4f} w={f['w']:.4f}"
    )
    print(f"  fit points: trunc={f['n_trunc']} round={f['n_round']}")
    print(
        f"  err/E: median={med:.3g} p90={ratios[int(0.9*len(ratios))]:.3g} "
        f"max={worst[0]:.3g}  over-envelope={over}/{len(ratios)}"
    )
    if worst[1]:
        w = worst[1]
        print(
            f"  worst: prec={w['prec']} sigma={w['sigma']:.4f} N={w['N']} m={w['m']}"
        )


def solve_m(f, s, eps, m_max, N=N_REF):
    """Mirror of the C tune(): smallest m meeting goal, else the argmin."""
    best_m, best_e = 1, float("inf")
    curve = []
    for m in range(1, m_max + 1):
        e = model(dict(sigma=s, m=m, eps=eps, N=N), f)
        curve.append((m, e))
        if e < best_e:
            best_e, best_m = e, m
    return curve, best_m, best_e


def validate(pts, key, f):
    """Compare the m the model picks with the m the measurement allows."""
    series = defaultdict(list)
    for r in pts:
        if r[key] > 0.0:
            series[(r["prec"], r["sigma"], r["N"])].append(r)
    goals = [10.0**-e for e in range(2, 18)]
    deltas = defaultdict(int)
    misses = 0
    checked = 0
    for (prec, s, N), v in series.items():
        v.sort(key=lambda r: r["m"])
        m_max = v[-1]["m"]
        eps = EPS[prec]
        meas = {r["m"]: r[key] for r in v}
        curve, argmin_m, _ = solve_m(f, s, eps, m_max, N)
        for goal in goals:
            mm = next((m for m, e in curve if e <= goal), None)
            reachable_meas = [m for m in sorted(meas) if meas[m] <= goal]
            if mm is None:
                # model says unreachable: report the argmin instead
                if reachable_meas:
                    deltas["model-says-unreachable-but-measurable"] += 1
                continue
            checked += 1
            if meas.get(mm, float("inf")) > goal:
                misses += 1
            if reachable_meas:
                deltas[mm - reachable_meas[0]] += 1
    print(f"  goal ladder: {checked} (sigma,N,goal) cases where the model returns an m")
    print(f"  cases where that m MISSES the goal in measurement: {misses}")
    # How pessimistic is the reported attainable floor?
    ratios = []
    for (prec, s2, N), v in series.items():
        v.sort(key=lambda r: r["m"])
        _, _, e_min_model = solve_m(f, s2, EPS[prec], v[-1]["m"], N)
        e_min_meas = min(r[key] for r in v)
        if e_min_meas > 0:
            ratios.append(e_min_model / e_min_meas)
    ratios.sort()
    print(
        f"  attainable floor, model/measured: median={ratios[len(ratios)//2]:.3g}"
        f" p90={ratios[int(0.9*len(ratios))]:.3g} max={ratios[-1]:.3g}"
    )
    num = {k: v for k, v in deltas.items() if isinstance(k, int)}
    tot = sum(num.values())
    if tot:
        line = "  extra m vs the smallest measured-sufficient m: " + ", ".join(
            f"+{k}:{v*100.0/tot:.0f}%" for k, v in sorted(num.items()))
        print(line)
    for k, v in deltas.items():
        if not isinstance(k, int):
            print(f"  {k}: {v}")


def main():
    global USE_N
    argv = list(sys.argv[1:])
    if "--no-n" in argv:
        USE_N = False
        argv.remove("--no-n")
    pts = load(argv)
    counts = defaultdict(int)
    for r in pts:
        counts[r["prec"]] += 1
    print(f"loaded {len(pts)} aggregated points: " + ", ".join(
        f"{k}={v}" for k, v in sorted(counts.items())))
    for key in DIRECTIONS:
        print(f"\n=== {key} ===")
        f = fit_direction(pts, key, margin=1.25)
        report(pts, key, f)
        print(
            f"  CONSTANTS {key}: a={f['a']:.6g} p={f['p']:.6g} r={f['r']:.6g} "
            f"c={f['c']:.6g} alpha={f['alpha']:.6g} q={f['q']:.6g} w={f['w']:.6g}"
        )
        validate(pts, key, f)


if __name__ == "__main__":
    main()
