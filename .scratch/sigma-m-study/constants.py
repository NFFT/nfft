#!/usr/bin/env python3
"""Print the chosen model's constants, ready to paste into kernel/nfft/tune.c.

The form is the one selectform.py picks: the two exact rates, a truncation
prefactor with a sigma and an m power, and a roundoff prefactor with a fitted
amplification rate.

    T(m,s)   = a * u^p * m^r * exp(-D(s)*m)      u = 1 - 1/sigma
    U(m,s,e) = c * e * u^q * exp(alpha*A(s)*m)

Usage: constants.py sweep-*.csv
"""

import sys

from fit import load
from selectform import fit_form, score

FORM = dict(use_r=True, use_qw=True, use_alpha=True, use_w=False)


def main():
    pts = load(sys.argv[1:])
    print(f"{len(pts)} aggregated points\n")
    out = {}
    for key in ("fwd", "adj"):
        f = fit_form(pts, key, FORM["use_r"], FORM["use_qw"], FORM["use_alpha"],
                      FORM["use_w"])
        s = score(pts, key, f)
        out[key] = f
        print(f"=== {key} ===")
        print(f"  a     = {f['a']:.6g}")
        print(f"  p     = {f['p']:.6g}")
        print(f"  r     = {f['r']:.6g}")
        print(f"  c     = {f['c']:.6g}")
        print(f"  alpha = {f['alpha']:.6g}")
        print(f"  q     = {f['q']:.6g}")
        print(
            f"  over-envelope={s['over']}  misses={s['misses']}  "
            f"exact={s['exact']:.0f}%  tail(>2)={s['tail']:.1f}%  "
            f"unreach={s['unreach']}  floor={s['floor']:.0f}x"
        )
        print()

    print("/* paste into kernel/nfft/tune.c */")
    for key, name in (("fwd", "tune_forward"), ("adj", "tune_adjoint")):
        f = out[key]
        print(f"static const tune_coeff {name} = {{")
        print(
            f"    K({f['a']:.6g}), K({f['p']:.6g}), K({f['r']:.6g}),"
        )
        print(f"    K({f['c']:.6g}), K({f['alpha']:.6g})}};")


if __name__ == "__main__":
    main()
