"""Fit the geometry-aware error model to the gsweep data.

fit.py fits E as a function of (sigma, m, eps) alone, because sweep.c ties
M to 2N and cannot separate the two. gsweep.c crosses them, so the same two
branches gain a geometry factor:

    T = a * u^p * m^r * exp(-D*m)      * N^tn * M^tm      u = 1 - 1/sigma
    U = c * eps * u^q * exp(alpha*A*m) * N^un * M^um
    E = T + U

The exponents are not free-form curve fitting -- they are the measure's own
normalisation. The forward error is a max over M nodes of a sum of N terms
with random phases, divided by an l1 norm over N terms, so the truncation
branch falls like N^-1/2 and rises slowly with M; the adjoint swaps the two.
The roundoff branch accumulates rather than cancels, so its N exponent has
the opposite sign.

Same envelope discipline as fit.py: fit the branches by least squares in log
space, then raise a and c to the tightest pair that dominates every measured
point, and apply a 1.25 margin.

Usage: python3 gfit.py <gsweep-*.csv>... [--holdout <csv>...]

Files after --holdout are excluded from the fit and validated separately.
The geometry exponents are the part most at risk of not extrapolating, so
the holdout should sit outside the fitted box -- larger N above all.
"""
import csv
import math
import sys
from collections import defaultdict

EPS = {"float": 2.0 ** -24, "double": 2.0 ** -53, "long_double": 2.0 ** -113}
DIRECTIONS = ("fwd", "adj")
MARGIN = 1.25


def D(s):
    return 2.0 * math.pi * math.sqrt(1.0 - 1.0 / s)


def A(s):
    return 2.0 * math.pi * (1.0 - 1.0 / (2.0 * s)) - D(s)


def lstsq(rowsA, rhs):
    k = len(rowsA[0])
    ata = [[sum(r[i] * r[j] for r in rowsA) for j in range(k)] for i in range(k)]
    atb = [sum(r[i] * y for r, y in zip(rowsA, rhs)) for i in range(k)]
    for i in range(k):
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
    """Max over trials per (prec, sigma, N, M, m) -- the envelope's target."""
    agg = defaultdict(lambda: [0.0, 0.0])
    for path in paths:
        with open(path) as fh:
            for r in csv.DictReader(fh):
                if r.get("err_fwd") is None or r.get("err_adj") is None:
                    continue
                key = (r["prec"], float(r["sigma"]), int(r["N"]), int(r["M"]),
                       int(r["m"]))
                cur = agg[key]
                cur[0] = max(cur[0], float(r["err_fwd"]))
                cur[1] = max(cur[1], float(r["err_adj"]))
    return [dict(prec=prec, sigma=s, N=N, M=M, m=m, eps=EPS[prec], fwd=ef,
                 adj=ea)
            for (prec, s, N, M, m), (ef, ea) in agg.items()]


def split_regimes(pts, key):
    series = defaultdict(list)
    for r in pts:
        if r[key] > 0.0:
            series[(r["prec"], r["sigma"], r["N"], r["M"])].append(r)
    trunc, round_ = [], []
    for v in series.values():
        v.sort(key=lambda r: r["m"])
        i_min = min(range(len(v)), key=lambda i: v[i][key])
        e_min = v[i_min][key]
        for i, r in enumerate(v):
            if i <= i_min - 2 and r[key] > 100.0 * e_min:
                trunc.append(r)
            elif i >= i_min + 2:
                round_.append(r)
    return trunc, round_


def fit_direction(pts, key, alpha_fixed=None):
    """Fit both branches and raise them to an envelope.

    alpha_fixed pins the roundoff branch's exponent instead of fitting it.
    That exponent is a correction to the exactly derived rate A = b - D, and
    where A*m stays small over the fitted range there is nothing for it to
    correct: the branch is flat in m and the least squares reads noise. Pass
    1.0 -- the derivation's own value -- for such a range.
    """
    trunc, round_ = split_regimes(pts, key)

    rows = [[1.0, math.log(1.0 - 1.0 / r["sigma"]), math.log(r["m"]),
             math.log(r["N"]), math.log(r["M"])] for r in trunc]
    rhs = [math.log(r[key]) + D(r["sigma"]) * r["m"] for r in trunc]
    sol = lstsq(rows, rhs)
    a, p, rr, tn, tm = math.exp(sol[0]), sol[1], sol[2], sol[3], sol[4]

    if alpha_fixed is None:
        rows = [[1.0, A(r["sigma"]) * r["m"], math.log(1.0 - 1.0 / r["sigma"]),
                 math.log(r["N"]), math.log(r["M"])] for r in round_]
        rhs = [math.log(r[key] / r["eps"]) for r in round_]
        sol = lstsq(rows, rhs)
        c, alpha, q, un, um = math.exp(sol[0]), sol[1], sol[2], sol[3], sol[4]
    else:
        alpha = alpha_fixed
        rows = [[1.0, math.log(1.0 - 1.0 / r["sigma"]), math.log(r["N"]),
                 math.log(r["M"])] for r in round_]
        rhs = [math.log(r[key] / r["eps"]) - alpha * A(r["sigma"]) * r["m"]
               for r in round_]
        sol = lstsq(rows, rhs)
        c, q, un, um = math.exp(sol[0]), sol[1], sol[2], sol[3]

    def raw(rec, aa, cc):
        s, m = rec["sigma"], rec["m"]
        u = 1.0 - 1.0 / s
        t = (aa * u ** p * m ** rr * math.exp(-D(s) * m)
             * rec["N"] ** tn * rec["M"] ** tm)
        try:
            v = (cc * rec["eps"] * u ** q * math.exp(alpha * A(s) * m)
                 * rec["N"] ** un * rec["M"] ** um)
        except OverflowError:
            v = float("inf")
        return t, v

    live = [r for r in pts if r[key] > 0.0]
    best = None
    for step in range(-30, 61):
        cc = c * (1.15 ** step)
        need = 0.0
        for rec in live:
            t, v = raw(rec, 1.0, cc)
            if t <= 0 or not math.isfinite(v):
                continue
            need = max(need, (rec[key] - v) / t)
        aa = max(need, 1e-12)
        slack = 0.0
        ok = True
        for rec in live:
            t, v = raw(rec, aa, cc)
            if not math.isfinite(t + v) or t + v <= 0:
                ok = False
                break
            slack += math.log((t + v) / rec[key])
        if ok and (best is None or slack < best[0]):
            best = (slack, aa, cc)
    _, a_env, c_env = best
    return dict(a=a_env * MARGIN, p=p, r=rr, tn=tn, tm=tm,
                c=c_env * MARGIN, alpha=alpha, q=q, un=un, um=um,
                n_trunc=len(trunc), n_round=len(round_))


def model(f, s, m, eps, N, M):
    u = 1.0 - 1.0 / s
    t = (f["a"] * u ** f["p"] * m ** f["r"] * math.exp(-D(s) * m)
         * N ** f["tn"] * M ** f["tm"])
    try:
        v = (f["c"] * eps * u ** f["q"] * math.exp(f["alpha"] * A(s) * m)
             * N ** f["un"] * M ** f["um"])
    except OverflowError:
        v = float("inf")
    return t + v


def validate(pts, key, f, label):
    """The only check that matters: the m the model picks against the m the
    measurement shows to be enough, over a ladder of goals."""
    series = defaultdict(list)
    for r in pts:
        if r[key] > 0.0:
            series[(r["prec"], r["sigma"], r["N"], r["M"])].append(r)
    goals = [10.0 ** -e for e in range(2, 30)]
    miss = 0
    extra = []
    cases = 0
    for (prec, s, N, M), v in series.items():
        v.sort(key=lambda r: r["m"])
        eps = EPS[prec]
        for g in goals:
            need = next((r["m"] for r in v if r[key] <= g), None)
            if need is None:
                continue
            got = next((r["m"] for r in v
                        if model(f, s, r["m"], eps, N, M) <= g), None)
            if got is None:
                continue
            cases += 1
            err_at_got = next(r[key] for r in v if r["m"] == got)
            if err_at_got > g:
                miss += 1
            extra.append(got - need)
    extra.sort()
    print("  %-8s cases=%-6d misses=%-4d extra m: mean %+.2f  median %+d  "
          "p90 %+d  max %+d  same-m %.0f %%"
          % (label, cases, miss, sum(extra) / len(extra), extra[len(extra)//2],
             extra[int(0.9 * len(extra))], extra[-1],
             100.0 * sum(1 for e in extra if e <= 0) / len(extra)))
    return miss


# What ships today, for a like-for-like baseline: same validation, same data,
# no geometry term.
SHIPPED = {
    "fwd": dict(a=3.0431, p=0.902205, r=-0.0183106, tn=0.0, tm=0.0,
                c=68.9787, alpha=0.967705, q=1.67263, un=0.0, um=0.0),
    "adj": dict(a=2.03698, p=0.234342, r=0.401585, tn=0.0, tm=0.0,
                c=33.5633, alpha=0.994013, q=1.15639, un=0.0, um=0.0),
}


def main(paths, holdout):
    pts = load(paths)
    hpts = load(holdout) if holdout else []
    print("points: %d fitted" % len(pts)
          + (", %d held out" % len(hpts) if hpts else "")
          + "  (N x M x sigma x m, worst trial)")
    if hpts:
        print("held-out N: %s" % sorted({r["N"] for r in hpts}))
    print()
    out = {}
    for key, name in (("fwd", "forward"), ("adj", "adjoint")):
        f = fit_direction(pts, key)
        out[key] = f
        print("%s: %d truncation points, %d roundoff points" %
              (name, f["n_trunc"], f["n_round"]))
        print("  T: a=%.6g p=%.6g r=%.6g  N^%.4f M^%.4f"
              % (f["a"], f["p"], f["r"], f["tn"], f["tm"]))
        print("  U: c=%.6g alpha=%.6g q=%.6g  N^%.4f M^%.4f"
              % (f["c"], f["alpha"], f["q"], f["un"], f["um"]))
        validate(pts, key, SHIPPED[key], "shipped")
        validate(pts, key, f, "geometry")
        if hpts:
            validate(hpts, key, SHIPPED[key], "shipped/holdout")
            validate(hpts, key, f, "geom/holdout")
        print()

    print("C constants for kernel/nfft/tune.c:")
    for key, name in (("fwd", "tune_forward"), ("adj", "tune_adjoint")):
        f = out[key]
        print("static const tune_coeff %s = {" % name)
        print("    K(%.6g), K(%.6g), K(%.6g), K(%.6g), K(%.6g),"
              % (f["a"], f["p"], f["r"], f["tn"], f["tm"]))
        print("    K(%.6g), K(%.6g), K(%.6g), K(%.6g), K(%.6g)};"
              % (f["c"], f["alpha"], f["q"], f["un"], f["um"]))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    argv = sys.argv[1:]
    if "--holdout" in argv:
        i = argv.index("--holdout")
        main(argv[:i], argv[i + 1:])
    else:
        main(argv, [])
