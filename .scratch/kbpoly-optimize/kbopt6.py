"""Degree reduction at equal accuracy, develop layout.

For each (sigma, m) and format: develop's stored degree 2 h_dev from its table
(kbpoly_data.h of the develop worktree), the accuracy of the plain Chebyshev
window at that degree (J_inf, J_2), and the smallest h at which

  plain   : the plain window itself is within 5% of its h_dev accuracy (both measures)
  LS      : the least-squares window is at least as accurate in BOTH measures
  LS-inf  : ... in the worst case only
  Lawson  : the Lawson window is at least as accurate in the worst case

Horner work per tap in develop's evaluator is proportional to h, so the
saving is 1 - h_min / h_dev.

Usage: python kbopt6.py MANT_DIG sigma [m ...]
"""
import math
import re
import sys
import numpy as np
from kbopt5 import SymModel

HEADER = "develop/kernel/util/kbpoly_data.h"
DEG_MAX = {113: 30, 64: 20, 53: 16, 24: 9}


def parse_table(mant):
    txt = open(HEADER).read()
    blocks = re.split(r"#(?:el)?if MANT_DIG == (\d+)|#else", txt)
    # blocks: [pre, '24', body, '64', body, '113', body, None, body(else), ...]
    tabs = {}
    i = 1
    while i < len(blocks) - 1:
        key = blocks[i]
        body = blocks[i + 1]
        fmt = int(key) if key else 53
        g = lambda name: re.search(name + r"[^{]*\{([^}]*)\}", body).group(1)
        tabs[fmt] = dict(
            mmax=int(re.search(r"KB_POLY_TAB_MMAX (\d+)", body).group(1)),
            br_opt=float(re.search(r"KB_POLY_BR_OPT K\(([\d.]+)\)", body).group(1)),
            deg0=[int(v) for v in g("kb_poly_tab_deg0").replace(",", " ").split()],
            off=[int(v) for v in g("kb_poly_tab_off").replace(",", " ").split()],
            brk=[int(v) for v in g("kb_poly_tab_brk").replace(",", " ").split()],
        )
        i += 2
    return tabs[mant]


def dev_h(tab, mant, sigma, m):
    b = math.pi * (2 - 1 / sigma)
    mt = math.ceil(tab["br_opt"] / b) - 1
    mt = min(m, mt, tab["mmax"])
    mt = max(mt, 1)
    deg = tab["deg0"][mt]
    for k in range(tab["off"][mt], tab["off"][mt + 1]):
        deg += 1 if sigma >= 1e-4 * tab["brk"][k] else 0
    deg = min(deg, DEG_MAX[mant])
    return (deg + 1) // 2


def main():
    mant = int(sys.argv[1]); sigma = float(sys.argv[2])
    ms = [int(v) for v in sys.argv[3:]]
    tab = parse_table(mant)
    eps = 2.0 ** (1 - mant)
    xi_d = np.linspace(0, 1 / (2 * sigma), 120); t_d = (np.arange(96) + 0.5) / 96
    print("MANT_DIG=%d sigma=%g   (h = half the stored degree; work per tap ~ h)" % (mant, sigma))
    print("%3s %5s %9s %8s | %6s %6s %6s %6s | %s" % (
        "m", "h_dev", "J_inf", "J_inf/eps", "plain", "LS", "LS-inf", "Lawson", "saving LS / Lawson"))
    for m in ms:
        h_dev = dev_h(tab, mant, sigma, m)
        ref = SymModel(sigma, m, h_dev)
        c_ref = ref.cheb_kb()
        Jref_inf, Jref_2 = ref.metrics(c_ref)
        found = {"plain": None, "LS": None, "LS-inf": None, "Lawson": None}
        rows = []
        for h in range(1, h_dev + 1):
            M = SymModel(sigma, m, h)
            c0 = M.cheb_kb()
            jp = M.metrics(c0)
            jl = M.metrics(M.fit_ls(c0))
            jw = M.metrics(M.fit_lawson(c0, xi_d, t_d, iters=20))
            rows.append((h, jp, jl, jw))
            if found["plain"] is None and jp[0] <= 1.05 * Jref_inf and jp[1] <= 1.05 * Jref_2:
                found["plain"] = h
            if found["LS"] is None and jl[0] <= Jref_inf and jl[1] <= Jref_2:
                found["LS"] = h
            if found["LS-inf"] is None and jl[0] <= Jref_inf:
                found["LS-inf"] = h
            if found["Lawson"] is None and jw[0] <= Jref_inf:
                found["Lawson"] = h
            if all(v is not None for v in found.values()):
                break
        f = lambda k: ("%d" % found[k]) if found[k] else "-"
        sv = lambda k: ("%3.0f%%" % (100 * (1 - found[k] / h_dev))) if found[k] else "  - "
        print("%3d %5d %9.2e %8.0f | %6s %6s %6s %6s | %s %s" % (
            m, h_dev, Jref_inf, Jref_inf / eps, f("plain"), f("LS"), f("LS-inf"), f("Lawson"),
            sv("LS"), sv("Lawson")))
        for h, jp, jl, jw in rows:
            print("      h=%2d  plain %.2e/%.2e   LS %.2e/%.2e   Lawson %.2e/%.2e" % (
                h, jp[0], jp[1], jl[0], jl[1], jw[0], jw[1]))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
