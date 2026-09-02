"""Print the reference tables for tests/window.c.

    uv run --with mpmath==1.3.0 python -m tests.windowref.generate

The output goes to stdout for pasting into tests/window.c; nothing in the tree
is written. Each table names the window evaluation it bounds and the parameters
it was taken at.

Sinc powers come from mpmath at 60 digits. Cardinal B-spline values come from
the truncated power basis in exact rational arithmetic: that formula cancels
catastrophically in floating point, but over the rationals it is exact, and it
shares nothing with the de Boor scheme it is used to check.
"""
import argparse
import sys
from fractions import Fraction
from math import comb, factorial

import mpmath as mp

DPS = 60
DIGITS = 45

# n = 512, N = 256, so sigma = 2 throughout.
N_GRID = 512
N_FREQ = 256
# One dyadic fraction of a grid cell, exact in every precision, so the run
# offset itself contributes no rounding.
FRAC = Fraction(5, 16)
# Every eighth frequency out to the band edge N/2.
K_STEP = 8

# The B-spline tables run at the cutoff nfft_init picks. The sinc-power tables
# run at 8 instead, where w = (2 sigma - 1)/(2 m sigma) is 3/32 and so exactly
# representable; at 11 it is 3/44, and its rounding would be charged to the
# window rather than to the constant.
M_BSPLINE = 11
M_SINCPOW = 8


def bspline(k, x):
    """The cardinal B-spline of order k on [0, k], exactly, at rational x."""
    if x <= 0 or x >= k:
        return Fraction(0)

    s = Fraction(0)

    for j in range(k + 1):
        d = x - j
        if d <= 0:
            break
        s += (-1) ** j * comb(k, j) * d ** (k - 1)

    return s / factorial(k - 1)


def mpf(q):
    return mp.mpf(q.numerator) / mp.mpf(q.denominator)


def cnum(v):
    if v == 0:
        return "0.0"
    s = mp.nstr(v, DIGITS, min_fixed=0, max_fixed=0, strip_zeros=False)
    return s if "e" in s or "E" in s else s + "e0"


def emit(name, values, comment, out):
    out.append(comment)
    out.append("static const R %s[] =" % name)
    out.append("{")
    for v in values:
        out.append("  K(%s)," % cnum(v))
    out.append("};")
    out.append("")


def bspline_phi_hut(out):
    """phi_hut(k) = (sin t / t)^(2m) / n, t = k pi / n."""
    m, n = M_BSPLINE, N_GRID
    ks = range(0, N_FREQ // 2 + 1, K_STEP)
    vals = []

    for k in ks:
        if k == 0:
            vals.append(mp.mpf(1) / n)
            continue
        t = mp.pi * k / n
        vals.append((mp.sin(t) / t) ** (2 * m) / n)

    emit("bspline_phi_hut_ref_512_256_11", vals,
         "/* B-spline phi_hut for n = %d, N = %d, m = %d, k = 0, %d, .. %d, at\n"
         " * %d digits. Regenerate with tests/windowref. */"
         % (n, N_FREQ, m, K_STEP, N_FREQ // 2, DIGITS), out)


def bspline_phi(out):
    """phi over one run: B_{2m}(t - l) / n, t = 2m + FRAC."""
    m, n = M_BSPLINE, N_GRID
    t = Fraction(2 * m) + FRAC
    vals = [mpf(bspline(2 * m, t - l) / n) for l in range(2 * m + 2)]
    peak = max(vals)

    emit("bspline_phi_ref_512_11", vals,
         "/* B-spline phi over one full run for n = %d, m = %d, with the node\n"
         " * %s of a cell past a grid point, exactly. The run straddles the\n"
         " * support, so its ends are zero and its tail is 1e-30 of the peak.\n"
         " * Regenerate with tests/windowref. */" % (n, m, FRAC), out)
    out.append("#define BSPLINE_PHI_PEAK_512_11 K(%s)" % cnum(peak))
    out.append("")


def sincpow_phi(out):
    """phi over one run: w (sin a / a)^(2m), a = pi w (nx0 - l)."""
    m, n = M_SINCPOW, N_GRID
    w = Fraction(2 * (N_GRID // N_FREQ) - 1, 2 * m * (N_GRID // N_FREQ))
    nx0 = Fraction(m) + FRAC
    vals = []

    for l in range(2 * m + 2):
        a = mp.pi * mpf(w) * mpf(nx0 - l)
        vals.append(mpf(w) * (mp.sin(a) / a) ** (2 * m))

    emit("sincpow_phi_ref_512_8", vals,
         "/* Sinc-power phi over one full run for n = %d, m = %d, sigma = 2,\n"
         " * with the node %s of a cell past a grid point, at %d digits.\n"
         " * Regenerate with tests/windowref. */" % (n, m, FRAC, DIGITS), out)
    out.append("#define SINCPOW_PHI_PEAK_512_8 K(%s)" % cnum(max(vals)))
    out.append("")


def sincpow_phi_hut(out):
    """phi_hut(k) = B_{2m}(k/(w n) + m)."""
    m, n = M_SINCPOW, N_GRID
    w = Fraction(2 * (N_GRID // N_FREQ) - 1, 2 * m * (N_GRID // N_FREQ))
    ks = range(0, N_FREQ // 2 + 1, K_STEP)
    vals = [mpf(bspline(2 * m, Fraction(k, 1) / (w * n) + m)) for k in ks]

    emit("sincpow_phi_hut_ref_512_256_8", vals,
         "/* Sinc-power phi_hut for n = %d, N = %d, m = %d, sigma = 2,\n"
         " * k = 0, %d, .. %d, exactly. Regenerate with tests/windowref. */"
         % (n, N_FREQ, m, K_STEP, N_FREQ // 2), out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dps", type=int, default=DPS)
    args = ap.parse_args(argv)

    mp.mp.dps = args.dps

    out = []
    bspline_phi_hut(out)
    bspline_phi(out)
    sincpow_phi(out)
    sincpow_phi_hut(out)
    sys.stdout.write("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
