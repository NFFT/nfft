"""Gate for the dyadic tuner: does a per-band refit tighten the cut-off?

The dyadic tuner offers three grids, n = 2^j * next_power_of_2(N) for j in
{0, 1, 2}. With t = next_power_of_2(N)/N in (1, 2], their oversamplings are
sigma = t, 2t, 4t, so each j occupies a disjoint band:

    j = 0   sigma in [5/4, 2]    (legal only when t >= 5/4)
    j = 1   sigma in (2, 4]      the legacy grid
    j = 2   sigma in (4, 8]

The shipped envelope is raised until it dominates every measurement anywhere
in sigma in [5/4, 4]. A per-band envelope only has to dominate its own band,
so it can only be tighter -- the question this script answers is whether it
is tighter by enough to cost one less cut-off, which is the whole of the
M-dominated shortfall against the legacy choice.

Usage: python3 dfit.py <gsweep-*.csv>...

Reports, per band and direction: misses (must stay 0), the extra cut-off the
envelope costs against the measured smallest sufficient one, and how often it
picks that one exactly. Baseline is the constants in kernel/nfft/tune.c,
validated on the same band-restricted population.
"""
import math
import sys

import gfit

# What kernel/nfft/tune.c carries today, for a like-for-like baseline.
SHIPPED = {
    "fwd": dict(a=47.1246, p=0.93096, r=-0.00821491, tn=-0.577425,
                tm=0.0480177, c=216.943, alpha=0.951476, q=2.25144,
                un=-0.0252356, um=0.11049),
    "adj": dict(a=24.7348, p=0.312453, r=0.434428, tn=0.0374448,
                tm=-0.506233, c=321.955, alpha=0.989474, q=1.48935,
                un=0.415192, um=-0.427671),
}

BANDS = (
    ("j=0  sigma [1.25, 2]", 1.25, 2.0),
    ("j=1  sigma [2, 4]", 2.0, 4.0),
    ("j=2  sigma [4, 8]", 4.0, 8.0),
)

# The roundoff branch rises as exp(alpha*A*m) with A = b - D, and alpha is a
# fitted correction to that exactly derived rate. A shrinks fast as the
# oversampling grows -- 0.96 at sigma = 5/4, 0.27 at 2, 0.056 at 4, 0.012 at
# 8 -- so on the widest band the branch is nearly flat over m <= 32 and the
# least squares has nothing to read. Where the band's largest A*m_max leaves
# the branch varying by less than e^2, alpha is pinned at the derivation's
# own value instead of fitted.
ALPHA_IDENTIFIABLE = 2.0


def band_alpha(lo, hi, m_max=32):
    """None to fit alpha, 1.0 to pin it."""
    a = gfit.A(lo)  # A falls with sigma, so the low edge is the best case
    return None if a * m_max >= ALPHA_IDENTIFIABLE else 1.0


def main(paths):
    pts = gfit.load(paths)
    print("points: %d  (N x M x sigma x m, worst trial)" % len(pts))
    print("sigma grid: %s" % sorted({r["sigma"] for r in pts}))
    print()

    fits = {}
    for j, (label, lo, hi) in enumerate(BANDS):
        band = [r for r in pts if lo <= r["sigma"] <= hi]
        sigmas = sorted({r["sigma"] for r in band})
        print("=== %s : %d points, %d distinct sigma" % (label, len(band),
                                                         len(sigmas)))
        if len(sigmas) < 2:
            # u = 1 - 1/sigma is then constant, so it is collinear with the
            # intercept and the least squares is singular.
            print("  need at least two distinct sigma to fit -- run "
                  "run_gsweep.sh with a sigma list covering this band\n")
            continue
        af = band_alpha(lo, hi)
        if af is not None:
            print("  alpha pinned at %.1f: A(%.2f)*32 = %.2f, below %.1f"
                  % (af, lo, gfit.A(lo) * 32, ALPHA_IDENTIFIABLE))
        for key, name in (("fwd", "forward"), ("adj", "adjoint")):
            f = gfit.fit_direction(band, key, alpha_fixed=af)
            fits[(j, key)] = f
            print("  %s" % name)
            print("    T: a=%.6g p=%.6g r=%.6g  N^%.4f M^%.4f"
                  % (f["a"], f["p"], f["r"], f["tn"], f["tm"]))
            print("    U: c=%.6g alpha=%.6g q=%.6g  N^%.4f M^%.4f"
                  % (f["c"], f["alpha"], f["q"], f["un"], f["um"]))
            gfit.validate(band, key, SHIPPED[key], "shipped")
            gfit.validate(band, key, f, "per-band")
        print()

    emit(fits)


def emit(fits):
    """The tables kernel/nfft/tune.c carries, one row per rung."""
    print("C constants for kernel/nfft/tune.c:")
    for key, name in (("fwd", "tune_dyadic_forward"),
                      ("adj", "tune_dyadic_adjoint")):
        print("static const tune_coeff %s[3] = {" % name)
        for j in range(3):
            f = fits.get((j, key))
            if f is None:
                print("    /* j=%d: NOT FITTED, band has no data */" % j)
                continue
            print("    {K(%.6g), K(%.6g), K(%.6g), K(%.6g), K(%.6g),"
                  % (f["a"], f["p"], f["r"], f["tn"], f["tm"]))
            print("     K(%.6g), K(%.6g), K(%.6g), K(%.6g), K(%.6g)},"
                  % (f["c"], f["alpha"], f["q"], f["un"], f["um"]))
        print("};")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
