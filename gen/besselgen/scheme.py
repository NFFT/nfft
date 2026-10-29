"""The I0 approximation scheme: what is fitted, in which variable, and how the
runtime evaluates it.

Two ranges, split at S:

  x <= S   I0(x) = 1 + y * P1(y),          y = (x/2)^2
           P1 is fitted so that the *relative* error of the whole expression
           is equioscillating. The leading 1 makes the form exact at x = 0 and
           every P1 coefficient is positive, so Horner cannot cancel.

  x > S    I0(x) = exp(x)/sqrt(x) * P2(u), u = 1/x in (0, 1/S]
           P2 is fitted against the exponentially scaled Bessel function
           g(u) = sqrt(1/u) * exp(-1/u) * I0(1/u), smooth and bounded away
           from zero on [0, 1/S], then shipped as monomial coefficients in u.
           The runtime then evaluates a plain sum, which splits into
           independent chains; the Clenshaw recurrence it replaced does not,
           and that cost the evaluation a factor of four.

           The monomial form is well conditioned only while the degree stays
           low, and that is what pins the split per format. See
           `branch2_growth` and the SPEC notes in generate.py.
"""
import mpmath as mp

from . import remez

# Mantissa bits -> (name, C guard, K() literal digits).
FORMATS = {
    24: ("float", "MANT_DIG == 24", 12),
    53: ("double", "MANT_DIG == 53", 22),
    64: ("long double", "MANT_DIG == 64", 25),
    113: ("quad", "MANT_DIG == 113", 40),
}


def i0(x):
    return mp.besseli(0, x)


def f1(y):
    """(I0(2*sqrt(y)) - 1) / y, extended by its limit 1 at y = 0."""
    if y == 0:
        return mp.mpf(1)
    return (i0(2 * mp.sqrt(y)) - 1) / y


def w1(y):
    """Weight turning the fit error of P1 into the relative error of 1+y*P1."""
    if y == 0:
        return mp.mpf(0)
    return y / i0(2 * mp.sqrt(y))


def g2(u):
    """sqrt(x) * exp(-x) * I0(x) at x = 1/u, extended by 1/sqrt(2*pi) at u = 0."""
    if u == 0:
        return 1 / mp.sqrt(2 * mp.pi)
    x = 1 / u
    return mp.sqrt(x) * mp.exp(-x) * i0(x)


def w2(u):
    return 1 / g2(u)


def _horner(c, x):
    r = c[-1]
    for j in range(len(c) - 2, -1, -1):
        r = r * x + c[j]
    return r


def verify_branch1(P1, split, npts=1200):
    """Relative error of the shipped form 1 + y*P1(y), in exact arithmetic.

    The fit runs in the Chebyshev basis but the runtime evaluates a monomial
    Horner chain, so the error that matters is this one; reporting the fit's
    own residual would hide any loss in the basis change.
    """
    worst = mp.mpf(0)
    for i in range(1, npts + 1):
        x = mp.mpf(split) * i / npts
        y = (x / 2) ** 2
        ref = i0(x)
        worst = max(worst, abs((1 + y * _horner(P1, y)) - ref) / ref)
    return worst


def verify_branch2(P2, split, npts=1200):
    """Relative error of the shipped monomial form P2(u), in exact arithmetic.

    The runtime groups the sum into four chains rather than one Horner pass.
    That is the same value in exact arithmetic; what the grouping changes is
    rounding, which `branch2_growth` bounds.
    """
    worst = mp.mpf(0)
    for i in range(npts + 1):
        u = (1 / mp.mpf(split)) * i / npts
        ref = g2(u)
        worst = max(worst, abs(_horner(P2, u) - ref) / ref)
    return worst


def branch2_growth(P2, split, npts=800):
    """max sum|c_j| u^j / |P2(u)| over the branch.

    1 means no ordering of the sum can cancel, so the four-chain grouping is
    as safe as Horner. This pins the split: the monomial coefficients of g2
    grow with the degree, and past roughly degree 24 the factor leaves 1 and
    climbs fast -- 466 at MANT_DIG 64 with a split of 10, 1.2e9 at MANT_DIG
    113 with a split of 15. Raising the split shortens branch 2 until the
    factor returns to 1. Branch 1 absorbs the degree, and its coefficients are
    positive at any degree, so it stays stable however long it gets.
    """
    hi = 1 / mp.mpf(split)
    worst = mp.mpf(0)
    for i in range(npts + 1):
        u = hi * i / npts
        num = sum(abs(c) * u ** j for j, c in enumerate(P2))
        den = abs(_horner(P2, u))
        if den > 0:
            worst = max(worst, num / den)
    return worst


def fit_branch1(split, n, prec_digits=None):
    """Monomial coefficients of P1 for `x <= split`, degree n."""
    with mp.workdps(prec_digits or mp.mp.dps):
        hi = (mp.mpf(split) / 2) ** 2
        c, _ = remez.minimax_poly(f1, w1, 0, hi, n)
        mono = remez.cheb_to_monomial(c, mp.mpf(0), hi)
        mono = [+v for v in mono]
        return mono, +verify_branch1(mono, split)


def fit_branch2(split, n, prec_digits=None):
    """Monomial coefficients of P2 in u = 1/x for `x > split`, degree n."""
    with mp.workdps(prec_digits or mp.mp.dps):
        hi = 1 / mp.mpf(split)
        c, _ = remez.minimax_poly(g2, w2, 0, hi, n)
        mono = [+v for v in remez.cheb_to_monomial(c, mp.mpf(0), hi)]
        return mono, +verify_branch2(mono, split)

