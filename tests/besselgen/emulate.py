"""Faithful fixed-precision emulation of the C runtime evaluation.

mpmath rounds every operation to the working precision with round-to-nearest,
ties-to-even, so evaluating at prec = MANT_DIG reproduces the C arithmetic of
the matching IEEE format (barring the exponent range, which the callers stay
inside). This is what lets the quad scheme be measured on a host with no
__float128 available, and it is cross-checked against the compiled C for
float, double and long double.
"""
import mpmath as mp


def horner(c, x):
    r = c[-1]
    for j in range(len(c) - 2, -1, -1):
        r = r * x + c[j]
    return r


def clenshaw_c(c, x):
    """The `evaluate_chebyshev` recurrence, as the pre-generator I0 wrote it."""
    n = len(c)
    a = c[n - 2]
    b = c[n - 1]
    for j in range(n - 2, 0, -1):
        t = c[j - 1] - b
        b = a + 2 * x * b
        a = t
    return a + x * b


def poly4(c, u):
    """The `poly4` four-chain sum as written in kernel/util/bessel_i0.c."""
    n = len(c)
    assert n >= 4 and n % 4 == 0, "table length must be a multiple of four"
    v = u * u
    w = v * v
    j = n - 4
    a = [c[j], c[j + 1], c[j + 2], c[j + 3]]
    for j in range(n - 8, -1, -4):
        for k in range(4):
            a[k] = a[k] * w + c[j + k]
    return (a[0] + u * a[1]) + v * (a[2] + u * a[3])


def i0_new(x, prec, P1, P2, split, split_exp=True):
    """Candidate runtime evaluation of I0 at `prec` mantissa bits."""
    with mp.workprec(prec):
        x = +abs(x)
        if x == 0:
            return mp.mpf(1)
        if x <= split:
            h = x * mp.mpf(0.5)
            y = h * h
            return 1 + y * poly4(P1, y)
        p = poly4(P2, 1 / x)
        if split_exp:
            # exp(x/2) twice keeps the intermediate below the overflow
            # threshold, so the usable range ends where I0 itself overflows.
            e = mp.exp(x * mp.mpf(0.5))
            return (e / mp.sqrt(x)) * e * p
        return (mp.exp(x) / mp.sqrt(x)) * p


def i0_current_cheb(x, prec, P1, Q1, P2, Q2, split):
    """The existing rational-Clenshaw evaluation (float/double/long double)."""
    with mp.workprec(prec):
        x = +abs(x)
        if x == 0:
            return mp.mpf(1)
        if x <= split:
            y = x * x
            return clenshaw_c(P1, y) / clenshaw_c(Q1, y)
        t = (mp.mpf(30) - x) / x
        return (mp.exp(x) / mp.sqrt(x)) * (clenshaw_c(P2, t) / clenshaw_c(Q2, t))


def ulp_error(got, ref, prec):
    if ref == 0:
        return mp.mpf(0)
    e = mp.frexp(abs(ref))[1]
    return abs(got - ref) / mp.mpf(2) ** (int(e) - prec)
