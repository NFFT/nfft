"""Optimised polynomial Kaiser-Bessel window in the develop layout.

Develop stores an even window on columns c = 0 .. m + 1, column c being the
cell [m + 1 - c, m + 2 - c) of |nx| as a polynomial in the offset s in
[-1/2, 1/2] from the cell centre, rows k = 0 .. 2h holding the coefficient
of s^k. The run of 2m + 2 taps for a node with t = frac(n x) uses columns
1 .. m + 1 at s = t - 1/2 and the same columns at -s for the mirror taps.

This script fits, for (sigma, m, h), the Chebyshev interpolant (develop's
fit), a least-squares fit and a Lawson minimax fit of the stencil cells
j = 0 .. m against the transform's error functional, predicts their errors,
and writes coefficient files for kbgain.c:

    line 1: deg (= 2h) m
    then 2h + 1 rows of m + 1 values: coefficient of s^k for columns 1 .. m+1

Usage: python kbopt5.py sigma m h [nxi nt]   (writes coef_<variant>.txt)
"""
import sys
import numpy as np
import mpmath as mp
from numpy.polynomial import chebyshev as Ch
from kbopt import kb_window, phi_hut, PI


class SymModel:
    def __init__(self, sigma, m, h, nxi=40, nt=32, nfine=240):
        self.sigma, self.m, self.h = sigma, m, h
        self.R = m + 1
        self.b = PI * (2 - 1 / sigma)
        self.deg = 2 * h
        self.nk = self.deg + 1                 # T_0 .. T_2h in 2s
        self.cells = np.arange(0, m + 1)       # j: cell [j, j+1), mirror [-j-1, -j)
        self.xi_max = 1 / (2 * sigma)
        self.xi = np.linspace(0, self.xi_max, nxi)
        self.t = (np.arange(nt) + 0.5) / nt
        self.xif = np.linspace(-self.xi_max, self.xi_max, nfine)
        self.tf = (np.arange(nfine) + 0.5) / nfine

    # window values on the 2m+2 taps from coefficients c[j, k] (T_k(2s))
    def values(self, c, t):
        s = t - 0.5
        Vp = Ch.chebvander(2 * s, self.deg)            # (nt, nk)
        Vm = Ch.chebvander(-2 * s, self.deg)
        P = Vp @ c.T                                    # cells j = 0..m at t  (nt, m+1)
        Pm = Vm @ c.T                                   # mirror cells -j-1 at t
        return P, Pm

    def D_of(self, c, xi, t):
        P, Pm = self.values(c, t)
        ph = np.exp(-2j * PI * np.outer(xi, self.cells))          # e^{-2pi i xi j}
        phm = np.exp(-2j * PI * np.outer(xi, -self.cells - 1))    # e^{-2pi i xi (-j-1)}
        return (ph @ P.T + phm @ Pm.T) * np.exp(-2j * PI * np.outer(xi, t))

    def D_exact(self, xi, t):
        y = t[:, None] + self.cells[None, :]
        P = kb_window(y, self.b, self.R)
        Pm = kb_window(t[:, None] - self.cells[None, :] - 1, self.b, self.R)
        ph = np.exp(-2j * PI * np.outer(xi, self.cells))
        phm = np.exp(-2j * PI * np.outer(xi, -self.cells - 1))
        return (ph @ P.T + phm @ Pm.T) * np.exp(-2j * PI * np.outer(xi, t))

    def design(self, xi, t):
        s = t - 0.5
        Vp = Ch.chebvander(2 * s, self.deg)
        Vm = Ch.chebvander(-2 * s, self.deg)
        ph = np.exp(-2j * PI * np.outer(xi, self.cells))
        phm = np.exp(-2j * PI * np.outer(xi, -self.cells - 1))
        et = np.exp(-2j * PI * np.outer(xi, t))
        # Phi[x, s, j, k] = (ph[x,j] Vp[s,k] + phm[x,j] Vm[s,k]) et[x,s]
        Phi = (ph[:, None, :, None] * Vp[None, :, None, :]
               + phm[:, None, :, None] * Vm[None, :, None, :]) * et[:, :, None, None]
        return Phi.reshape(len(xi) * len(t), (self.m + 1) * self.nk)

    def C(self, xi):
        return phi_hut(xi, self.b, self.R)

    def metrics(self, c=None):
        D = self.D_exact(self.xif, self.tf) if c is None else self.D_of(c, self.xif, self.tf)
        E = D / self.C(self.xif)[:, None] - 1
        return np.abs(E).max(), np.sqrt(np.mean(np.abs(E) ** 2))

    def single_freq(self, c, ks, n):
        """Predicted max over t of |E| for the frequencies k on a grid of n."""
        xi = np.array(ks, dtype=float) / n
        D = self.D_exact(xi, self.tf) if c is None else self.D_of(c, xi, self.tf)
        return np.abs(D / self.C(xi)[:, None] - 1).max(axis=1)

    def cheb_kb(self):
        """Develop's fit: 2h+1 Chebyshev points per cell of |nx|."""
        n = self.nk
        node = np.cos(PI * (np.arange(n) + 0.5) / n) / 2       # s in [-1/2, 1/2]
        c = np.empty((self.m + 1, n))
        for j in self.cells:
            v = kb_window(j + 0.5 + node, self.b, self.R)
            c[j] = Ch.chebfit(2 * node, v, self.deg)
        return c

    def _solve(self, A, rhs, w=None):
        if w is not None:
            A = A * w[:, None]; rhs = rhs * w
        Ar = np.vstack([A.real, A.imag]); br = np.concatenate([rhs.real, rhs.imag])
        x, *_ = np.linalg.lstsq(Ar, br, rcond=None)
        return x

    def fit_ls(self, c0, xi=None, t=None):
        xi = self.xi if xi is None else xi; t = self.t if t is None else t
        A = self.design(xi, t) / np.repeat(self.C(xi), len(t))[:, None]
        r = A @ c0.ravel() - 1
        d = self._solve(A, -r)
        return c0 + d.reshape(c0.shape)

    def fit_lawson(self, c0, xi, t, iters=20):
        A = self.design(xi, t) / np.repeat(self.C(xi), len(t))[:, None]
        rhs = -(A @ c0.ravel() - 1)
        w = np.ones(len(rhs)) / len(rhs); best, bx = np.inf, None
        for _ in range(iters):
            x = self._solve(A, rhs, np.sqrt(w))
            r = np.abs(A @ x - rhs)
            if r.max() < best:
                best, bx = r.max(), x
            w = w * r; w /= w.sum()
        return c0 + bx.reshape(c0.shape)

    def to_monomial(self, c):
        """Rows k = 0..2h of the coefficient of s^k per column 1..m+1
        (column c = cell m + 1 - c), converted in 50-digit arithmetic."""
        mp.mp.dps = 50
        rows = np.zeros((self.nk, self.m + 1))
        for j in self.cells:
            # series sum_k c[j,k] T_k(2s): T_k(x) in powers of x, then x = 2s
            poly = [mp.mpf(0)] * self.nk
            tk = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(1)]]       # T_0, T_1 in x
            for k in range(2, self.nk):
                prev, pp = tk[-1], tk[-2]
                nxt = [mp.mpf(0)] * (k + 1)
                for i, v in enumerate(prev):
                    nxt[i + 1] += 2 * v
                for i, v in enumerate(pp):
                    nxt[i] -= v
                tk.append(nxt)
            for k in range(self.nk):
                for i, v in enumerate(tk[k]):
                    poly[i] += mp.mpf(float(c[j, k])) * v * mp.mpf(2) ** i
            col = self.m + 1 - j
            rows[:, col - 1] = [float(v) for v in poly]
        return rows


def write_coef(path, M, rows):
    with open(path, "w") as f:
        f.write("%d %d\n" % (M.deg, M.m))
        for k in range(M.nk):
            f.write(" ".join("%.17e" % v for v in rows[k]) + "\n")


def main():
    sigma, m, h = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    nxi = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    nt = int(sys.argv[5]) if len(sys.argv) > 5 else 32
    tag = sys.argv[6] if len(sys.argv) > 6 else ""
    M = SymModel(sigma, m, h, nxi=nxi, nt=nt)
    N = 2048; n = int(round(sigma * N)); n += n % 2   # as kbgain.c sizes the grid
    ks = [0, N // 8, N // 4, N // 2 - 1]
    print("=== sigma=%g m=%d deg=%d (h=%d)  b=%.4f  KB (truncated, exact): J_inf=%.3e J_2=%.3e"
          % ((sigma, m, M.deg, h, M.b) + M.metrics()))
    c0 = M.cheb_kb()
    variants = [("cheb", c0)]
    variants.append(("ls", M.fit_ls(c0)))
    xi_d = np.linspace(0, M.xi_max, 120); t_d = (np.arange(96) + 0.5) / 96
    variants.append(("lawson", M.fit_lawson(c0, xi_d, t_d, iters=20)))
    print("%-8s %10s %10s   predicted max|E| at k = %s (n=%d)" % ("variant", "J_inf", "J_2", ks, n))
    for name, c in variants:
        ji, j2 = M.metrics(c)
        sf = M.single_freq(c, ks, n)
        print("%-8s %10.3e %10.3e   %s" % (name, ji, j2, " ".join("%.3e" % v for v in sf)))
        write_coef("coef_%s_s%g_m%d_h%d%s.txt" % (name, sigma, m, h, tag), M, M.to_monomial(c))
    sf = M.single_freq(None, ks, n)
    print("%-8s %10s %10s   %s" % ("exactKB", "", "", " ".join("%.3e" % v for v in sf)))


if __name__ == "__main__":
    main()
