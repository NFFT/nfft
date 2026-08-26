"""Arbitrary-precision direct NDFT / NDCT / NDST (forward + adjoint).

All functions operate at the current mpmath working precision (mpmath.mp.dps),
which the caller sets.
"""
import itertools

import mpmath

_REAL_MODULES = ("nfct", "nfst")

# NDFT "problem type" per axis: TYPE_I is the legacy symmetric-about-zero range;
# TYPE_II (even N only) is the ascending range with the Nyquist term (+n/2) moved
# to the last slot instead of straddling zero. Odd N has no ascending/symmetric
# split (n//2 and n-n//2 already coincide at the boundary), so TYPE_II collapses
# to TYPE_I for odd n.
TYPE_I = 0
TYPE_II = 1


def _freq_axis_nfft(n, variant):
    if variant == TYPE_II and n % 2 == 0:
        return list(range(-(n // 2) + 1, n - (n // 2) + 1))  # ascending; +n/2 last
    return list(range(-(n // 2), n - (n // 2)))              # type-I; odd-safe


def freqs(module, N, variant=None):
    """Per-dimension frequency ranges (lists of ints), one list per dimension.

    variant, if given, is a list of per-axis TYPE_I/TYPE_II tags (nfft only).
    variant=None is byte-identical to the legacy (type-I-only) behavior."""
    if module == "nfft":
        var = variant if variant is not None else [TYPE_I] * len(N)
        return [_freq_axis_nfft(n, v) for n, v in zip(N, var)]
    if module == "nfct":
        return [list(range(0, n)) for n in N]
    if module == "nfst":
        return [list(range(1, n)) for n in N]
    raise ValueError("unknown module: %r" % (module,))


def index_set(module, N, variant=None):
    """All frequency tuples, tensor product with dimension 0 varying slowest."""
    return list(itertools.product(*freqs(module, N, variant)))


def nn(module, N, variant=None):
    """Number of Fourier coefficients."""
    count = 1
    for r in freqs(module, N, variant):
        count *= len(r)
    return count


def _basis(module, k, xj, adjoint):
    """basis(k, x_j). For NFFT the adjoint conjugates (flips the exponent sign)."""
    twopi = 2 * mpmath.pi
    if module == "nfft":
        dot = mpmath.fsum([ki * xji for ki, xji in zip(k, xj)])
        ang = twopi * dot
        return mpmath.expj(ang) if adjoint else mpmath.expj(-ang)
    if module == "nfct":
        prod = mpmath.mpf(1)
        for ki, xji in zip(k, xj):
            prod *= mpmath.cos(twopi * ki * xji)
        return prod
    if module == "nfst":
        prod = mpmath.mpf(1)
        for ki, xji in zip(k, xj):
            prod *= mpmath.sin(twopi * ki * xji)
        return prod
    raise ValueError("unknown module: %r" % (module,))


def _zero(module):
    return mpmath.mpf(0) if module in _REAL_MODULES else mpmath.mpc(0)


def trafo(module, N, M, x, f_hat, variant=None):
    """Forward: f[j] = sum_k f_hat[k] * basis(k, x[j]).  Returns list of length M."""
    K = index_set(module, N, variant)
    out = []
    for j in range(M):
        acc = _zero(module)
        for idx, k in enumerate(K):
            acc += f_hat[idx] * _basis(module, k, x[j], adjoint=False)
        out.append(acc)
    return out


def adjoint(module, N, M, x, f, variant=None):
    """Adjoint: f_hat[k] = sum_j f[j] * conj(basis)(k, x[j]).  Length = nn()."""
    K = index_set(module, N, variant)
    out = []
    for k in K:
        acc = _zero(module)
        for j in range(M):
            acc += f[j] * _basis(module, k, x[j], adjoint=True)
        out.append(acc)
    return out
