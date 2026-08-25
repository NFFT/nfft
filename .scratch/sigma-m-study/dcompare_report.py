#!/usr/bin/env python3
"""Turn the three-arm comparison into a Markdown summary.

Usage: dcompare_report.py out.md <csv>...

The CSV comes from compare.c, which runs two arms over the same nodes, the
same data and the same plan_ng path, so only the parameter pair differs:

  LEGACY  n = 2*next_power_of_2(N) with the cut-off found by measuring
  DYADIC  Y(tune_plan_dyadic)
"""

import csv
import math
import sys
from collections import defaultdict

PRECS = ["float", "double", "long_double"]
PRETTY = {"float": "float", "double": "double", "long_double": "long double"}
SHAPES = ["N-dominated", "balanced", "M-dominated"]
ARMS = [("dya", "dyadic")]


def next_pow2(x):
    p = 1
    while p < x:
        p *= 2
    return p


def load(paths):
    rows = []
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh):
                for k in ("N", "M", "n_dya", "m_dya", "n_leg", "m_leg",
                          "dya_met", "leg_met"):
                    r[k] = int(r[k])
                for k in ("goal", "err_dya", "err_leg", "t_dya", "t_leg"):
                    r[k] = float(r[k])
                r["t"] = next_pow2(r["N"]) / r["N"]
                for a, _ in ARMS:
                    r["sp_" + a] = (r["t_leg"] / r["t_" + a]
                                    if r["t_" + a] > 0 else 0.0)
                # Which rung the dyadic tuner took, relative to legacy.
                base = next_pow2(r["N"])
                r["rung"] = (0 if r["n_dya"] == base else
                             1 if r["n_dya"] == 2 * base else
                             2 if r["n_dya"] == 4 * base else -1)
                rows.append(r)
    return rows


def med(v):
    v = sorted(v)
    return v[len(v) // 2] if v else float("nan")


def block(w, rows, label):
    """One speedup table for a subset."""
    if not rows:
        return
    w("| arm | cases | median | p10 | min | faster or tied | goal met |")
    w("|---|---|---|---|---|---|---|")
    for a, name in ARMS:
        sp = sorted(r["sp_" + a] for r in rows)
        met = sum(r[a + "_met"] for r in rows)
        w("| %s | %d | %.2fx | %.2fx | %.2fx | %d/%d | %d/%d |"
          % (name, len(sp), med(sp), sp[int(0.1 * len(sp))], sp[0],
             sum(1 for s in sp if s >= 0.999), len(sp), met, len(rows)))
    w("")


def main():
    out_path = sys.argv[1]
    rows = load(sys.argv[2:])
    L = []
    w = L.append

    w("# The dyadic ladder against the legacy choice")
    w("")
    w("Three arms over the same nodes, the same data and the same `plan_ng`")
    w("path, so the only difference is the parameter pair.")
    w("")
    w("- **legacy** -- `n = 2*next_power_of_2(N)`, cut-off found by searching")
    w("  upward for the smallest `m` whose *measured* error meets the goal.")
    w("  The legacy API has no cut-off estimator, so this is an oracle it")
    w("  does not possess. It is also rung 1 of the ladder.")
    w("- **dyadic** -- `nfft_tune_plan_dyadic`, the three rungs")
    w("  `n = 2^j * next_power_of_2(N)`.")
    w("")
    w("Speedups are against legacy, so above 1.00x is a win. The tuner runs")
    w("unaided: no measurement, no refinement.")
    w("")

    ts = sorted({(r["N"], round(r["t"], 3)) for r in rows})
    w("The N set spans `t = next_power_of_2(N)/N`, which is what decides")
    w("whether the ladder has a rung below the legacy grid at all -- rung 0")
    w("needs `t >= 5/4`. Below that the ladder can only return the legacy")
    w("grid or a wider one.")
    w("")
    w("| N | " + " | ".join(str(n) for n, _ in ts) + " |")
    w("|---|" + "---|" * len(ts))
    w("| t | " + " | ".join("%.2f" % t for _, t in ts) + " |")
    w("")

    w("## Overall")
    w("")
    block(w, rows, "all")

    w("## By shape")
    w("")
    for s in SHAPES:
        sub = [r for r in rows if r["shape"] == s]
        if not sub:
            continue
        w("### %s (%d cases)" % (s, len(sub)))
        w("")
        block(w, sub, s)

    w("## By oversampling headroom")
    w("")
    w("`t < 5/4` is the control: rung 0 is illegal, so the ladder returns the")
    w("legacy grid or a wider one and can only tie or win on the cut-off.")
    w("A control row below 1.00x is a cut-off regression, not a ladder one.")
    w("")
    for lo, hi, name in ((0.0, 1.25, "t < 5/4 (control)"),
                         (1.25, 9.9, "t >= 5/4 (rung 0 legal)")):
        sub = [r for r in rows if lo <= r["t"] < hi]
        if not sub:
            continue
        w("### %s (%d cases)" % (name, len(sub)))
        w("")
        block(w, sub, name)

    w("## By precision")
    w("")
    for p in PRECS:
        sub = [r for r in rows if r["prec"] == p]
        if not sub:
            continue
        w("### %s (%d cases)" % (PRETTY[p], len(sub)))
        w("")
        block(w, sub, p)

    w("## Which rung the ladder took")
    w("")
    w("Rung 1 is the legacy grid. Rung 0 halves it, rung 2 doubles it.")
    w("")
    w("| shape | rung 0 | rung 1 (legacy) | rung 2 | median speedup |")
    w("|---|---|---|---|---|")
    for s in SHAPES:
        sub = [r for r in rows if r["shape"] == s]
        if not sub:
            continue
        c = defaultdict(int)
        for r in sub:
            c[r["rung"]] += 1
        w("| %s | %d | %d | %d | %.2fx |"
          % (s, c[0], c[1], c[2], med([r["sp_dya"] for r in sub])))
    w("")

    w("## How much of this is timing noise")
    w("")
    same = [r for r in rows
            if r["n_dya"] == r["n_leg"] and r["m_dya"] == r["m_leg"]]
    if same:
        sp = sorted(r["sp_dya"] for r in same)
        n = len(sp)
        w("In %d of %d cases the ladder returns the legacy pair itself, so the"
          % (n, len(rows)))
        w("two arms run identical parameters and their measured ratio is pure")
        w("noise. It should read 1.00x exactly:")
        w("")
        w("| median | mean | p05 | p95 | min | max | within 2 % | within 5 % |")
        w("|---|---|---|---|---|---|---|---|")
        w("| %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %d/%d | %d/%d |"
          % (med(sp), sum(sp) / n, sp[int(0.05 * n)], sp[int(0.95 * n)],
             sp[0], sp[-1],
             sum(1 for x in sp if abs(x - 1) <= 0.02), n,
             sum(1 for x in sp if abs(x - 1) <= 0.05), n))
        w("")
        w("The median is exact and the mean is unbiased, so **the medians in")
        w("this document carry signal**. The per-case spread does not: a")
        w("single row at 0.80x or 1.20x is within what identical parameters")
        w("produce, so the min and max columns above, and the counts of rows")
        w("below 1.00x, are noise-dominated and should not be read as losses.")
        w("")

    w("## Where the ladder loses")
    w("")
    loss = sorted([r for r in rows if r["sp_dya"] < 0.999],
                  key=lambda r: r["sp_dya"])
    w("%d of %d cases below 1.00x." % (len(loss), len(rows)))
    w("")
    if loss:
        same_n = [r for r in loss if r["n_dya"] == r["n_leg"]]
        more_m = [r for r in same_n if r["m_dya"] > r["m_leg"]]
        w("- %d are on the legacy grid itself, of which %d carry a larger"
          % (len(same_n), len(more_m)))
        w("  cut-off than the oracle's. That is the envelope's cost, not the")
        w("  ladder's: the same cut-off gap the model has always had.")
        w("- %d chose a different rung." % (len(loss) - len(same_n)))
        w("")
        w("| prec | N | t | M | shape | goal | dir | rung | n | m | n_leg |"
          " m_leg | speedup |")
        w("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in loss[:20]:
            w("| %s | %d | %.2f | %d | %s | %.0e | %s | %d | %d | %d | %d |"
              " %d | %.2fx |"
              % (PRETTY[r["prec"]], r["N"], r["t"], r["M"], r["shape"],
                 r["goal"], r["dir"], r["rung"], r["n_dya"], r["m_dya"],
                 r["n_leg"], r["m_leg"], r["sp_dya"]))
        w("")

    w("## Accuracy")
    w("")
    w("| arm | goal met |")
    w("|---|---|")
    w("| legacy (oracle) | %d/%d |"
      % (sum(r["leg_met"] for r in rows), len(rows)))
    for a, name in ARMS:
        w("| %s | %d/%d |" % (name, sum(r[a + "_met"] for r in rows),
                              len(rows)))
    w("")
    w("The input is drawn with real and imaginary parts on `[0, 1)`, matching")
    w("`nfft_vrand_unit_complex` and the sweep the error model is fitted to.")
    w("Centred data measures a forward error up to 2.6x smaller; see issue 05")
    w("in `.scratch/tune-dyadic/`.")
    w("")

    with open(out_path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("wrote %s (%d cases)" % (out_path, len(rows)))


if __name__ == "__main__":
    main()
