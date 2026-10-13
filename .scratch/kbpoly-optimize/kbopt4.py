"""Two follow-ups.

A. Why the degree table is conservative: per degree, the table's tail bound
   2 sum_{j>d} max_cells |c_j| / peak  versus the measured contribution of the
   Chebyshev truncation to the transform error,
   max |D_cheb - D_KB| / |C|, and versus the Kaiser-Bessel error itself.

B. Working-precision fit: the least-squares and Lawson fits done in float32
   (complex64), in the residual formulation (fit the correction to the
   Chebyshev interpolant, not the window), compared with the float64 fits.
"""
import sys
import numpy as np
from numpy.polynomial import chebyshev as Ch
from kbopt import Model, kb_window, PI


def tail_bound(M, d, nc=44):
    n = nc
    nodes = (np.cos(PI * (np.arange(n) + 0.5) / n) + 1) / 2
    env = np.zeros(n)
    for a in range(-(M.m + 2), M.m + 2):
        v = kb_window(a + nodes, M.b, M.R)
        env = np.maximum(env, np.abs(Ch.chebfit(2 * nodes - 1, v, n - 1)))
    peak = kb_window(np.array([0.0]), M.b, M.R)[0]
    return 2 * env[d + 1:].sum() / peak


def part_a(sigma, m):
    probe = Model(sigma, m, 1)
    dlut = probe.lut_degree(53)
    P_exact = kb_window(probe.tf[:, None] + probe.cells[None, :], probe.b, probe.R)
    D_kb = probe.D_of_values(P_exact, probe.xif, probe.tf)
    Cf = probe.C_lib(probe.xif)[:, None]
    J_kb = np.abs(D_kb / Cf - 1).max()
    tol = 0.05 * PI * np.exp(-probe.b * probe.R)
    print("\n=== A: sigma=%g m=%d  table d=%d   KB J_inf=%.2e   rule reference pi e^-bR=%.2e (KB is %.1fx it), tol 5%% of it=%.2e"
          % (sigma, m, dlut, J_kb, PI * np.exp(-probe.b * probe.R), J_kb / (PI * np.exp(-probe.b * probe.R)), tol))
    print("  d | tail bound | measured poly part | bound/measured | poly part / KB | chebKB J_inf / KB")
    for d in range(max(2, dlut - 5), dlut + 1):
        M = Model(sigma, m, d)
        c = M.cheb_kb()
        D_c = M.D_of_values(M.values(c, M.tf), M.xif, M.tf)
        poly = np.abs((D_c - D_kb) / Cf).max()
        J_c = np.abs(D_c / Cf - 1).max()
        tb = tail_bound(M, d)
        print(" %2d |   %.2e |          %.2e |        %6.1f |         %6.3f |  %.3f"
              % (d, tb, poly, tb / poly, poly / J_kb, J_c / J_kb))


def fit_residual(M, C, dtype, lawson=0):
    """Fit the correction delta to the Chebyshev interpolant in the given real
    dtype: rows A = Phi/C, rhs = -(A c0 - 1). Returns c0 + delta in float64."""
    cdt = np.complex64 if dtype == np.float32 else np.complex128
    c0 = M.cheb_kb().astype(dtype)
    A = (M.design(M.xi, M.t) / np.repeat(C, len(M.t))[:, None]).astype(cdt)
    r = (A @ c0.ravel().astype(cdt)) - dtype(1)        # residual, in working precision
    rhs = -r
    Ar = np.vstack([A.real, A.imag]).astype(dtype)
    br = np.concatenate([rhs.real, rhs.imag]).astype(dtype)
    w = np.ones(len(rhs), dtype=dtype) / dtype(len(rhs))
    best, bx = np.inf, None
    for it in range(max(1, lawson)):
        sw = np.concatenate([np.sqrt(w), np.sqrt(w)])
        x, *_ = np.linalg.lstsq(Ar * sw[:, None], br * sw, rcond=None)
        res = np.abs((A @ x.astype(cdt)) - rhs)
        if res.max() < best:
            best, bx = res.max(), x
        if lawson:
            w = w * res.astype(dtype)
            w = w / w.sum()
    return (c0.ravel() + bx).astype(np.float64).reshape(M.ncell, M.deg + 1)


def part_b(sigma, m, nxi=40, nt=32):
    probe = Model(sigma, m, 1)
    dlut = probe.lut_degree(53)
    print("\n=== B: sigma=%g m=%d  (float eps 6e-8)" % (sigma, m))
    print("  d |   KB J_inf | LS f64 J_inf/J_2       | LS f32 J_inf/J_2       | Lawson20 f64 J_inf | Lawson20 f32 J_inf")
    for d in (dlut - 2, dlut - 1):
        M = Model(sigma, m, d, nxi=nxi, nt=nt)
        C = M.C_lib(M.xi)
        P_exact = kb_window(M.tf[:, None] + M.cells[None, :], M.b, M.R)
        J_kb = M.metrics(P_exact, M.C_lib, poly=False)[0]
        ls64 = M.metrics(fit_residual(M, C, np.float64), M.C_lib)
        ls32 = M.metrics(fit_residual(M, C, np.float32), M.C_lib)
        mm64 = M.metrics(fit_residual(M, C, np.float64, lawson=20), M.C_lib)[0]
        mm32 = M.metrics(fit_residual(M, C, np.float32, lawson=20), M.C_lib)[0]
        print(" %2d |   %.2e | %.2e / %.2e    | %.2e / %.2e    |          %.2e |          %.2e"
              % (d, J_kb, ls64[0], ls64[1], ls32[0], ls32[1], mm64, mm32))


if __name__ == "__main__":
    s, m = float(sys.argv[1]), int(sys.argv[2])
    if sys.argv[3] == "a":
        part_a(s, m)
    else:
        part_b(s, m)
