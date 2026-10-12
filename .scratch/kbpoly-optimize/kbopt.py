"""Feasibility study: optimise the per-cell polynomial Kaiser-Bessel window
of the NFFT against the transform's own aliasing error, instead of fitting
the KB function per cell.

Model (grid units, y = n x, offset t = frac(y), cells a = -R .. R-1, R = m+1):

  psi(a + t) = p_a(t)                      window actually used (2m+2 taps)
  D(xi, t)   = sum_a p_a(t) e^{-2 pi i xi (a + t)}   one-point rule for psi_hat(xi)
  E(xi, t)   = D(xi, t) / C(xi) - 1        relative error of frequency xi = k/n
                                           at a node with offset t

  C(xi) = I0(R sqrt(b^2 - (2 pi xi)^2)) / I0(R b)   (the deconvolution the
          library uses, the FT of the entire sinh/sin function)

Worst case error  J_inf = max_{|xi| <= 1/(2 sigma), t} |E|
RMS error         J_2   = sqrt(mean |E|^2)        (uniform in xi and t)

Variants per degree d (one Chebyshev basis per cell, (2m+2)(d+1) unknowns):
  chebKB   per-cell Chebyshev interpolant of KB (what the branch does)
  chebKB*  the same, deconvolved by the polynomial's own FT psi_hat
  LS       least squares in (xi, t) for the library's C
  LS*      least squares with C re-optimised (alternating, C = psi_hat_poly)
  MM       Lawson IRLS minimax for the library's C
  MM*      minimax with C = psi_hat_poly
  exactKB  the truncated KB window itself (no polynomial): the target the
           degree table aims 5% above.
"""
import sys
import numpy as np
from numpy.polynomial import chebyshev as Ch
from scipy.special import i0e, roots_legendre

PI = np.pi


def kb_window(y, b, R):
    """sinh(b ra)/(pi ra)/I0(R b), ra = sqrt(R^2 - y^2), for |y| <= R, else the
    sin branch; same normalisation as the library (phi_hut(0) = 1)."""
    y = np.asarray(y, dtype=float)
    lg_tail = np.log(i0e(R * b))
    out = np.empty_like(y)
    inside = np.abs(y) < R
    ra = np.sqrt(np.maximum(R * R - y[inside] ** 2, 0.0))
    # sinh(b ra)/I0(Rb) = exp(b ra - R b - lg_tail) (1 - exp(-2 b ra)) / 2
    e = np.exp(-b * y[inside] ** 2 / (ra + R) - lg_tail)
    val = e * (-np.expm1(-2 * b * ra)) / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(ra > 1e-9, val / (PI * ra), b / PI * np.exp(-R * b - lg_tail))
    out[inside] = v
    ri = np.sqrt(np.maximum(y[~inside] ** 2 - R * R, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        vo = np.where(ri > 1e-9, np.sin(b * ri) / (PI * ri), b / PI)
    out[~inside] = vo * np.exp(-R * b - lg_tail)
    return out


def phi_hut(xi, b, R):
    """I0(R sqrt(b^2 - w^2)) / I0(R b), w = 2 pi xi (the library's PHI_HUT)."""
    w = 2 * PI * np.asarray(xi, dtype=float)
    s = np.sqrt(np.maximum((b - w) * (b + w), 0.0))
    return np.exp(-R * w * w / (s + b)) * i0e(R * s) / i0e(R * b)


class Model:
    def __init__(self, sigma, m, deg, nxi=40, nt=32, nfine=160):
        self.sigma, self.m, self.deg = sigma, m, deg
        self.R = m + 1
        self.b = PI * (2 - 1 / sigma)
        self.cells = np.arange(-self.R, self.R)          # 2m+2 cells
        self.ncell = len(self.cells)
        self.xi_max = 1 / (2 * sigma)
        # coarse design grid (for the fits)
        self.xi = np.linspace(0, self.xi_max, nxi)
        self.t = (np.arange(nt) + 0.5) / nt
        # fine grid (for the metrics)
        self.xif = np.linspace(-self.xi_max, self.xi_max, nfine)
        self.tf = (np.arange(nfine) + 0.5) / nfine
        # Gauss-Legendre on [0,1] for psi_hat of a polynomial window
        xg, wg = roots_legendre(deg + 12)
        self.tg, self.wg = (xg + 1) / 2, wg / 2

    # polynomial values p_a(t) from Chebyshev coefficients c[a, j]
    def values(self, c, t):
        V = Ch.chebvander(2 * t - 1, self.deg)         # (nt, deg+1)
        return V @ c.T                                   # (nt, ncell)

    # D(xi, t) from a table P[t, a] of window values on the cells
    def D_of_values(self, P, xi, t):
        ph = np.exp(-2j * PI * np.outer(xi, self.cells))      # (nxi, ncell)
        return (ph @ P.T) * np.exp(-2j * PI * np.outer(xi, t))  # (nxi, nt)

    def design(self, xi, t):
        """Phi[(xi,t), (a,j)] with D = Phi @ c.ravel()."""
        V = Ch.chebvander(2 * t - 1, self.deg)                # (nt, deg+1)
        ph = np.exp(-2j * PI * np.outer(xi, self.cells))      # (nxi, ncell)
        et = np.exp(-2j * PI * np.outer(xi, t))               # (nxi, nt)
        # Phi[x, s, a, j] = ph[x,a] et[x,s] V[s,j]
        Phi = ph[:, None, :, None] * et[:, :, None, None] * V[None, :, None, :]
        return Phi.reshape(len(xi) * len(t), self.ncell * (self.deg + 1))

    def psi_hat_poly(self, c, xi):
        """Exact FT of the piecewise polynomial window on the grid xi."""
        P = self.values(c, self.tg)                           # (ng, ncell)
        ph = np.exp(-2j * PI * np.outer(xi, self.cells))
        et = np.exp(-2j * PI * np.outer(xi, self.tg)) * self.wg
        return np.real((ph * (et @ P)).sum(axis=1))

    def C_lib(self, xi):
        return phi_hut(xi, self.b, self.R)

    def metrics(self, P_or_c, C_fn, poly=True):
        if poly:
            P = self.values(P_or_c, self.tf)
        else:
            P = P_or_c
        D = self.D_of_values(P, self.xif, self.tf)
        E = D / C_fn(self.xif)[:, None] - 1
        return np.abs(E).max(), np.sqrt(np.mean(np.abs(E) ** 2))

    # ---- fits -------------------------------------------------------------
    def cheb_kb(self):
        """Per-cell Chebyshev interpolant of KB at deg+1 Chebyshev points."""
        n = self.deg + 1
        k = np.arange(n)
        nodes = (np.cos(PI * (k + 0.5) / n) + 1) / 2
        c = np.empty((self.ncell, n))
        for i, a in enumerate(self.cells):
            v = kb_window(a + nodes, self.b, self.R)
            c[i] = Ch.chebfit(2 * nodes - 1, v, self.deg)
        return c

    def lut_degree(self, mant=53, nc=44):
        """The branch's table rule: smallest d with
        2 sum_{j>d} max_cells |c_j| <= max(0.05 pi e^{-bR}, 2^(1-mant)) peak."""
        n = nc
        k = np.arange(n)
        nodes = (np.cos(PI * (k + 0.5) / n) + 1) / 2
        env = np.zeros(n)
        for a in range(-(self.m + 2), self.m + 2):
            v = kb_window(a + nodes, self.b, self.R)
            env = np.maximum(env, np.abs(Ch.chebfit(2 * nodes - 1, v, n - 1)))
        tail = np.cumsum(env[::-1])[::-1]        # tail[j] = sum_{i>=j} env[i]
        peak = kb_window(np.array([0.0]), self.b, self.R)[0]
        tol = max(0.05 * PI * np.exp(-self.b * self.R), 2.0 ** (1 - mant)) * peak
        for d in range(n - 1):
            if 2 * tail[d + 1] <= tol:
                return d
        return n - 1

    def _solve(self, A, rhs, w=None):
        if w is not None:
            A = A * w[:, None]
            rhs = rhs * w
        Ar = np.vstack([A.real, A.imag])
        br = np.concatenate([rhs.real, rhs.imag])
        x, *_ = np.linalg.lstsq(Ar, br, rcond=None)
        return x

    def fit_ls(self, C):
        """Least squares: rows Phi/C, rhs 1."""
        Phi = self.design(self.xi, self.t)
        Cr = np.repeat(C, len(self.t))
        x = self._solve(Phi / Cr[:, None], np.ones(len(Cr)))
        return x.reshape(self.ncell, self.deg + 1)

    def fit_minimax(self, C, iters=120):
        """Lawson's IRLS on |Phi c / C - 1| (discrete minimax)."""
        Phi = self.design(self.xi, self.t)
        Cr = np.repeat(C, len(self.t))
        A = Phi / Cr[:, None]
        rhs = np.ones(len(Cr))
        w = np.ones(len(Cr)) / len(Cr)
        best, bestx = np.inf, None
        for _ in range(iters):
            x = self._solve(A, rhs, np.sqrt(w))
            r = np.abs(A @ x - rhs)
            if r.max() < best:
                best, bestx = r.max(), x
            w = w * r
            w /= w.sum()
        return bestx.reshape(self.ncell, self.deg + 1)

    def fit_alt(self, kind, rounds=4):
        """Alternate the polynomial fit with C := psi_hat of the polynomial."""
        C = self.C_lib(self.xi)
        for _ in range(rounds):
            c = self.fit_ls(C) if kind == "ls" else self.fit_minimax(C)
            C = self.psi_hat_poly(c, self.xi)
        return c

    def monomial_growth(self, c):
        """max |monomial coefficient| / peak, a conditioning indicator for the
        run-time Horner form in t."""
        peak = kb_window(np.array([0.0]), self.b, self.R)[0]
        g = 0.0
        for row in c:
            mono = Ch.cheb2poly(row)        # in s = 2t - 1
            # substitute s = 2t - 1
            pt = np.polynomial.Polynomial([-1, 2])
            q = np.polynomial.Polynomial(mono)(pt)
            g = max(g, np.abs(q.coef).max())
        return g / peak


def run_case(sigma, m, degs=None, out=sys.stdout):
    probe = Model(sigma, m, 1)
    dlut = probe.lut_degree(53)
    dlutf = probe.lut_degree(24)
    R, b = probe.R, probe.b
    if degs is None:
        degs = list(range(max(1, dlut - 4), dlut + 2))
    # exact truncated KB as a value table on the fine grid
    P_exact = kb_window(probe.tf[:, None] + probe.cells[None, :], b, R)
    Jinf0, J20 = probe.metrics(P_exact, probe.C_lib, poly=False)
    out.write("\n=== sigma=%g m=%d R=%d b=%.4f  e^-bR=%.2e  LUT degree double=%d float=%d\n"
              % (sigma, m, R, b, np.exp(-b * R), dlut, dlutf))
    out.write("exact truncated KB, library C:   J_inf=%.3e  J_2=%.3e\n" % (Jinf0, J20))
    out.write("%4s | %9s %9s | %9s %9s | %9s %9s | %9s %9s | %9s %9s | %9s %9s | %6s %6s\n" % (
        "deg", "chebKB", "", "chebKB*", "", "LS", "", "LS*", "", "MM", "", "MM*", "", "grow", "growMM"))
    out.write("%4s | %9s %9s | %9s %9s | %9s %9s | %9s %9s | %9s %9s | %9s %9s |\n" % (
        "", "J_inf", "J_2", "J_inf", "J_2", "J_inf", "J_2", "J_inf", "J_2", "J_inf", "J_2", "J_inf", "J_2"))
    for d in degs:
        M = Model(sigma, m, d)
        res = []
        c0 = M.cheb_kb()
        res += M.metrics(c0, M.C_lib)
        res += M.metrics(c0, lambda xi: M.psi_hat_poly(c0, xi))
        c1 = M.fit_ls(M.C_lib(M.xi))
        res += M.metrics(c1, M.C_lib)
        c2 = M.fit_alt("ls")
        res += M.metrics(c2, lambda xi: M.psi_hat_poly(c2, xi))
        c3 = M.fit_minimax(M.C_lib(M.xi))
        res += M.metrics(c3, M.C_lib)
        res += (float("nan"), float("nan"))
        g0, g3 = M.monomial_growth(c0), M.monomial_growth(c3)
        out.write("%4d | " % d + " | ".join("%9.2e %9.2e" % (res[i], res[i + 1]) for i in range(0, 12, 2))
                  + " | %6.1f %6.1f\n" % (g0, g3))
        out.flush()


if __name__ == "__main__":
    cases = [(2.0, 2), (2.0, 3), (2.0, 5), (1.25, 6), (1.5, 4), (4.0, 2)]
    if len(sys.argv) > 1:
        cases = [(float(sys.argv[1]), int(sys.argv[2]))]
    for s, m in cases:
        run_case(s, m)
