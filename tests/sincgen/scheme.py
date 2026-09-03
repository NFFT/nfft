"""The log|sinc| approximation scheme: what is fitted and how the runtime
evaluates it.

Two ranges, split at S:

  |x| <= S   log|sinc(x)| = y * Q(y),  y = x*x, Horner in y.
             Q is fitted so that the *relative* error of y*Q(y) equioscillates.
             The leading y makes the form exact at x = 0, and the Taylor
             coefficients of log sinc are all negative
             (-1/6, -1/180, -1/2835, ...), so the Horner chain cannot cancel.

  |x| >  S   log|sinc(x)| = LOG(FABS(SIN(x)/x)), evaluated directly.
             Past the split sinc is far enough from 1 that the logarithm keeps
             its argument's relative accuracy; no coefficients are needed.

The quantity callers need is the *absolute* error of the result, because they
form EXP(2m * log_sinc(x)): a relative error there is 2m times the absolute
error here. Branch 1 delivers |L| * eps, which is what makes the whole scheme
worth having -- the direct LOG(SIN(x)/x) can only manage eps, independent of
how small |L| is.
"""
import mpmath as mp

from ..besselgen import remez

# Mantissa bits -> (name, C guard, K() literal digits).
FORMATS = {
    24: ("float", "MANT_DIG == 24", 12),
    53: ("double", "MANT_DIG == 53", 22),
    64: ("long double", "MANT_DIG == 64", 25),
    113: ("quad", "MANT_DIG == 113", 40),
}


def log_sinc(x):
    """log(sin(x)/x), extended by its limit 0 at x = 0."""
    if x == 0:
        return mp.mpf(0)
    return mp.log(mp.sin(x) / x)


def fq(y):
    """log(sinc(sqrt(y))) / y, extended by its limit -1/6 at y = 0."""
    if y == 0:
        return -mp.mpf(1) / 6
    return log_sinc(mp.sqrt(y)) / y


def wq(y):
    """Weight turning the fit error of Q into the relative error of y*Q(y).

    fq is bounded away from zero on the whole branch (it runs from -1/6 at
    y = 0 to about -0.197 at y = 4), so unlike the I0 branch-1 weight this one
    vanishes nowhere and needs no stepped-in search interval.
    """
    return 1 / abs(fq(y))


def _horner(c, x):
    r = c[-1]
    for j in range(len(c) - 2, -1, -1):
        r = r * x + c[j]
    return r


def verify_branch(Q, split, npts=1600):
    """Relative error of the shipped form y*Q(y), in exact arithmetic.

    The fit runs in the Chebyshev basis but the runtime evaluates a monomial
    Horner chain, so this is the error that matters; the fit's own residual
    would hide any loss in the basis change.
    """
    worst = mp.mpf(0)
    for i in range(1, npts + 1):
        x = mp.mpf(split) * i / npts
        y = x * x
        ref = log_sinc(x)
        worst = max(worst, abs(y * _horner(Q, y) - ref) / abs(ref))
    return worst


def horner_growth(Q, split, npts=800):
    """max sum|c_j| y^j / |sum c_j y^j| over the branch.

    1 means the Horner chain cannot cancel. This is what pins the split: at
    S = 2 the factor is 1, at S = 2.5 it is already about 53.
    """
    worst = mp.mpf(0)
    for i in range(npts + 1):
        y = (mp.mpf(split) ** 2) * i / npts
        num = sum(abs(c) * y ** j for j, c in enumerate(Q))
        den = abs(sum(c * y ** j for j, c in enumerate(Q)))
        if den > 0:
            worst = max(worst, num / den)
    return worst


def fit_branch(split, n, prec_digits=None):
    """Monomial coefficients of Q for |x| <= split, degree n."""
    with mp.workdps(prec_digits or mp.mp.dps):
        hi = mp.mpf(split) ** 2
        c, _ = remez.minimax_poly(fq, wq, 0, hi, n)
        mono = remez.cheb_to_monomial(c, mp.mpf(0), hi)
        mono = [+v for v in mono]
        return mono, +verify_branch(mono, split)
