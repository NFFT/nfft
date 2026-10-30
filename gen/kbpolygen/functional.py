"""The transform's error functional of a piecewise polynomial window.

The NFFT applies the window as 2m + 2 taps psi(a + t), a = -(m+1) .. m, at
the node's offset t = frac(n x) into its grid cell, and divides the
coefficients by C(xi) = PHI_HUT, xi = k / n. For an even window stored on the
cells j = 0 .. m of |nx| as polynomials in the offset s = t - 1/2 from the
cell centre (the layout of kernel/util/kbpoly.c), the relative error of the
frequency xi at a node with offset t is

    E(xi, t) = D(xi, t) / C(xi) - 1,
    D(xi, t) = sum_j [ p_j(s) e^(-2 pi i xi (j + t)) + p_j(-s) e^(-2 pi i xi (t - j - 1)) ].

D is linear in the polynomial coefficients. What the degree table needs is
the part of E the polynomial adds on top of the exact window: with the
coefficients of the interpolation error e_{j,k} (Chebyshev basis T_k(2s)),

    poly_part = max_{xi, t} |D_e(xi, t)| / C(xi),

the same sum built from e instead of p. The sum has no cancellation of
large terms, so it is evaluated in double from the mpmath error coefficients
and is accurate for every format, down to the quad floor.
"""
import numpy as np

NXI, NT = 48, 32


def grid(sigma):
    xi = np.linspace(0.0, 1.0 / (2.0 * sigma), NXI)
    t = (np.arange(NT) + 0.5) / NT
    return xi, t


def cheb_values(nk, s):
    """T_k(2 s), k < nk, at the offsets s: array (len(s), nk)."""
    x = 2.0 * s
    T = np.empty((len(s), nk))
    T[:, 0] = 1.0
    if nk > 1:
        T[:, 1] = x
    for k in range(2, nk):
        T[:, k] = 2.0 * x * T[:, k - 1] - T[:, k - 2]
    return T


def poly_part(err_coef, C, sigma, m):
    """max |D_e| / C over the grid; err_coef[j][k] for j = 0 .. m, C[i] at
    xi[i] of grid(sigma)."""
    xi, t = grid(sigma)
    s = t - 0.5
    e = np.asarray(err_coef, dtype=float)              # (m+1, nk)
    nk = e.shape[1]
    T = cheb_values(nk, s)                              # (nt, nk)
    sign = (-1.0) ** np.arange(nk)
    P = T @ e.T                                         # p_j(s)      (nt, m+1)
    Pm = (T * sign) @ e.T                               # p_j(-s)     (nt, m+1)
    j = np.arange(m + 1)
    ph = np.exp(-2j * np.pi * np.outer(xi, j))          # (nxi, m+1)
    phm = np.exp(-2j * np.pi * np.outer(xi, -j - 1))
    D = (ph @ P.T + phm @ Pm.T) * np.exp(-2j * np.pi * np.outer(xi, t))
    return float((np.abs(D) / np.asarray(C, dtype=float)[:, None]).max())
