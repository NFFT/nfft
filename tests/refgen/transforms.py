"""Arbitrary-precision direct NDFT / NDCT / NDST (forward + adjoint).

All functions operate at the current mpmath working precision (mpmath.mp.dps),
which the caller sets. Frequency index sets and ordering reproduce the C f_hat
layout byte-for-byte (see docs/agents/test-methodology.md).
"""
import itertools

import mpmath

_REAL_MODULES = ("nfct", "nfst")


def freqs(module, N):
    """Per-dimension frequency ranges (lists of ints), one list per dimension."""
    if module == "nfft":
        return [list(range(-(n // 2), n - (n // 2))) for n in N]
    if module == "nfct":
        return [list(range(0, n)) for n in N]
    if module == "nfst":
        return [list(range(1, n)) for n in N]
    raise ValueError("unknown module: %r" % (module,))


def index_set(module, N):
    """All frequency tuples, tensor product with dimension 0 varying slowest."""
    return list(itertools.product(*freqs(module, N)))


def nn(module, N):
    """Number of Fourier coefficients."""
    count = 1
    for r in freqs(module, N):
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


def trafo(module, N, M, x, f_hat):
    """Forward: f[j] = sum_k f_hat[k] * basis(k, x[j]).  Returns list of length M."""
    K = index_set(module, N)
    out = []
    for j in range(M):
        acc = _zero(module)
        for idx, k in enumerate(K):
            acc += f_hat[idx] * _basis(module, k, x[j], adjoint=False)
        out.append(acc)
    return out


def adjoint(module, N, M, x, f):
    """Adjoint: f_hat[k] = sum_j f[j] * conj(basis)(k, x[j]).  Length = nn()."""
    K = index_set(module, N)
    out = []
    for k in K:
        acc = _zero(module)
        for j in range(M):
            acc += f[j] * _basis(module, k, x[j], adjoint=True)
        out.append(acc)
    return out
