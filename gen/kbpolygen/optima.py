# /// script
# requires-python = ">=3.10"
# dependencies = ["python-flint==0.9.0", "numpy", "mpmath==1.3.0"]
# ///
"""Emit kernel/util/kbpoly_opt_data.h, the optimized polynomial windows.

    uv run --with python-flint==0.9.0 --with numpy --with mpmath==1.3.0 \\
        python -m gen.kbpolygen.optima -j 10

Per cell (sigma, m, h) the table of Y(kb_poly_init) at degree 2h with the
correction of Y(kb_poly_correct) for OPT_WINDOW, solved to 0.1% instead of
1%, in 200-bit ball arithmetic (python-flint). The algorithm is the one in
kernel/util/kbpoly.c: the design of 24 x 9 points, one weighted least-squares
pass from the interpolant, safeguarded Dinkelbach steps and bisection on
second-order cone subproblems (socp3.c), exchange rounds on a 401 x 65 fine
grid. mpmath is needed only by generate.py, which gives the cell list.
"""

import argparse
import math
import os
import sys
import time
from multiprocessing import Pool

import flint
import numpy as np
from flint import arb, arb_mat

from gen.kbpolygen import generate

PREC = 200
PREC_QR = 2 * PREC + 50  # Gram matrix of the rank-revealing factorization
SIGMAS = (1.1, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0)
NXI, NT = 24, 16  # design: frequencies, offsets 0 and NT / 2 midpoints
FXI, FT = 401, 64  # fine grid
ADD, ROUNDS, NPTS = 80, 40, 400
DELTA, DELTA_X, GAIN = 0.001, 0.005, 1.3
TOL_DINK, TOL_BISECT = 1e-4, DELTA / 100
SOLVES, MAXIT = 30, 60
GL_EXTRA = 24  # quadrature nodes h + GL_EXTRA: exact to 200 bits
DIGITS = 36
OUT = "kernel/util/kbpoly_opt_data.h"
LOG = ".scratch/kbpoly-optimize/results/optima_log.txt"


def _ufunc(f):
    return np.frompyfunc(f, 1, 1)


MID = _ufunc(lambda a: a.mid())
SQRT = _ufunc(lambda a: a.sqrt())
COS = _ufunc(lambda a: (2 * a).cos_pi())  # cos(2 pi a)
SIN = _ufunc(lambda a: (2 * a).sin_pi())
FIN = _ufunc(lambda a: a.is_finite())
FLT = np.vectorize(lambda a: float(a.mid()), otypes=[float])


def eps():
    return arb(2) ** (1 - PREC)


def amat(a):
    a = np.asarray(a, dtype=object)
    if a.ndim == 1:
        a = a[:, None]
    return arb_mat(a.shape[0], a.shape[1], a.ravel().tolist())


def arr(mat):
    return np.array(mat.entries(), dtype=object).reshape(mat.nrows(), mat.ncols())


def mm(a, b):
    """Product of object arrays through arb_mat; 1-d b gives 1-d."""
    out = arr(amat(a) * amat(b))
    return out[:, 0] if np.ndim(b) == 1 else out


def zeros(*shape):
    z = np.empty(shape, dtype=object)
    z.ravel()[:] = [arb(0)] * z.size
    return z


def cheb(nk, x):
    """T_k(x), k < nk, for the arbs x: (len(x), nk)."""
    T = zeros(len(x), nk)
    T[:, 0] = arb(1)
    if nk > 1:
        T[:, 1] = x
    for k in range(2, nk):
        T[:, k] = 2 * x * T[:, k - 1] - T[:, k - 2]
    return T


class Cell:
    """Window, quadrature and design geometry of one (sigma, m, h)."""

    def __init__(self, sigma, m, h):
        self.sigma, self.m, self.h = sigma, m, h
        self.nc, self.nk, self.w = m + 1, 2 * h + 1, m + 2
        sig = arb(repr(float(sigma)))
        self.b = arb.pi() * (2 - 1 / sig)
        self.ximax = (1 / (2 * sig)).mid()
        ng = h + GL_EXTRA
        nodes = [arb.legendre_p_root(ng, k, weight=True) for k in range(ng)]
        self.tg = np.array([((1 - x) / 2).mid() for x, _ in nodes], dtype=object)
        self.wg = np.array([(w / 2).mid() for _, w in nodes], dtype=object)
        self.sgn = np.array([(-1) ** k for k in range(self.nk)], dtype=object)
        # T_k(2 s) at the quadrature nodes, weighted: (ng, nk).
        self.wT = cheb(self.nk, 2 * self.tg - 1) * self.wg[:, None]

    def phi(self, nu):
        """The Kaiser-Bessel window over I0(R b), R = m + 1, as kb_phi."""
        R = self.m + 1
        z = R * R - nu * nu
        peak = (R * self.b).bessel_i(0)
        if z > 0:
            s = z.sqrt()
            return (self.b * s).sinh() / (arb.pi() * s) / peak
        s = (-z).sqrt()
        return (self.b * s).sin() / (arb.pi() * s) / peak

    def interpolant(self):
        """Chebyshev coefficients in T_j(2s) of the interpolant per column,
        (nk, w), on the nodes of Y(kb_poly_fit)."""
        n = self.nk
        th = [arb(k * 2 + 1) / (2 * n) for k in range(n)]
        c = zeros(n, self.w)
        for col in range(self.w):
            cen = arb(self.m + 1 - col) + arb(1) / 2
            val = [self.phi(cen + th[k].cos_pi() / 2) for k in range(n)]
            for j in range(n):
                s = sum((val[k] * (j * th[k]).cos_pi() for k in range(n)), arb(0))
                c[j, col] = (2 * s / n / (2 if j == 0 else 1)).mid()
        return c

    def mono(self, ch):
        """Chebyshev (nk, w) to monomials in s, the layout of the library."""
        n = self.nk
        tp = [[0] * n for _ in range(n)]
        tp[0][0] = 1
        if n > 1:
            tp[1][1] = 2
        for j in range(2, n):
            for i in range(j + 1):
                tp[j][i] = -tp[j - 2][i] + (4 * tp[j - 1][i - 1] if i else 0)
        tpt = np.array([[arb(tp[j][i]) for j in range(n)] for i in range(n)])
        return MID(mm(tpt, ch))

    def dmat(self, mono, xi, t):
        """Real and imaginary part of the tap sum D(xi, t), (len xi, len t)."""
        m, nk = self.m, self.nk
        s = t - arb(1) / 2
        S = zeros(len(t), nk)
        S[:, 0] = arb(1)
        for k in range(1, nk):
            S[:, k] = S[:, k - 1] * s
        cols = mono[:, 1:]
        lo, hi = mm(S, cols), mm(S * self.sgn, cols)
        j = np.array([arb(m - c) for c in range(m + 1)], dtype=object)
        alo, ahi = -xi[:, None] * j, xi[:, None] * (j + 1)
        LH = np.vstack([lo.T, hi.T])
        Ar = mm(np.hstack([COS(alo), COS(ahi)]), LH)
        Ai = mm(np.hstack([SIN(alo), SIN(ahi)]), LH)
        at = -xi[:, None] * t[None, :]
        ce, se = COS(at), SIN(at)
        return ce * Ar - se * Ai, ce * Ai + se * Ar

    def errors(self, mono, xi, t):
        """|E(xi, t)| = |D - C| / C with C = psi_hat by quadrature."""
        Dr, Di = self.dmat(mono, xi, t)
        C = mm(self.dmat(mono, xi, self.tg)[0], self.wg)
        if not all(C > 0):
            return None
        Er, Ei = (Dr - C[:, None]) / C[:, None], Di / C[:, None]
        return SQRT(Er * Er + Ei * Ei)

    def worst(self, mono):
        """Fine-grid maximum of |E| before refinement and the ADD largest
        local maxima, refined on 11 x 11 patches, as (|E|, xi, t) by
        decreasing |E|; worst_points() of kbpoly.c."""
        dx, dt = self.ximax / (FXI - 1), arb(1) / FT
        fx = np.array([dx * i for i in range(FXI)], dtype=object)
        ft = np.array([arb(0)] + [(q - arb(1) / 2) * dt for q in range(1, FT + 1)])
        Ea = self.errors(mono, fx, ft)
        if Ea is None:
            return None, []
        ev = FLT(Ea)
        nt = FT + 1
        up = np.vstack([ev[1:], np.full((1, nt), -1.0)])
        dn = np.vstack([np.full((1, nt), -1.0), ev[:-1]])
        qd = np.r_[1, np.arange(nt - 1)]
        qu = np.r_[np.arange(1, nt), nt - 2]
        loc = (up <= ev) & (dn <= ev) & (ev[:, qd] <= ev) & (ev[:, qu] <= ev)
        idx = sorted(zip(*np.nonzero(loc)), key=lambda p: -ev[p])[:ADD]
        out = []
        for i, q in idx:
            best = (Ea[i, q], fx[i], ft[q])
            k5 = [arb(k - 5) / 5 for k in range(11)]
            px = [min(max(fx[i] + k * dx, arb(0)), self.ximax) for k in k5]
            pt = [ft[q] + k * dt for k in k5]
            pt = [p + 1 if p < 0 else (p - 1 if p >= 1 else p) for p in pt]
            px, pt = np.array(px, dtype=object), np.array(pt, dtype=object)
            Ep = self.errors(mono, px, pt)
            if Ep is None:
                return None, []
            fp = FLT(Ep)
            a, c = np.unravel_index(np.argmax(fp), fp.shape)
            if fp[a, c] > float(best[0]):
                best = (Ep[a, c], px[a], pt[c])
            out.append(best)
        out.sort(key=lambda p: -float(p[0]))
        return Ea[np.unravel_index(np.argmax(ev), ev.shape)], out

    def rows(self, xi, t):
        """Rows of the points (xi, t) on the full unknowns u = c nk + k (the
        T_k(2s) coefficient of column c + 1): Re and Im of D - C, and C."""
        m, nc, nk = self.m, self.nc, self.nk
        T = cheb(nk, 2 * t - 1)
        j = np.array([arb(m - c) for c in range(nc)], dtype=object)
        alo = -xi[:, None] * (j[None, :] + t[:, None])
        ahi = xi[:, None] * (j[None, :] + 1 - t[:, None])
        Ts = T * self.sgn
        Dr = (
            T[:, None, :] * COS(alo)[:, :, None] + Ts[:, None, :] * COS(ahi)[:, :, None]
        )
        Di = (
            T[:, None, :] * SIN(alo)[:, :, None] + Ts[:, None, :] * SIN(ahi)[:, :, None]
        )
        ag = -xi[:, None] * self.tg[None, :]
        Fr, Fi = mm(COS(ag), self.wT), mm(SIN(ag), self.wT)
        a0 = -xi[:, None] * j[None, :]
        c0, s0 = COS(a0), SIN(a0)
        B = 2 * (Fr[:, None, :] * c0[:, :, None] - Fi[:, None, :] * s0[:, :, None])
        n = len(xi)
        return [MID(a.reshape(n, nc * nk)) for a in (Dr - B, Di, B)]


def pivchol(A, trunc):
    """Column-pivoted Cholesky of A^T A, A real (rows, cols): the R of a
    pivoted QR of A, stopped at R_kk <= trunc R_00. Returns the pivots S and
    R^-1 (k, k) at PREC_QR."""
    with_prec = flint.ctx.prec
    flint.ctx.prec = PREC_QR
    try:
        Am = amat(A)
        Gm = arr(Am.transpose() * Am)
        n = Gm.shape[0]
        d = np.array([Gm[i, i] for i in range(n)], dtype=object)
        L = zeros(n, n)
        S, free = [], np.ones(n, dtype=bool)
        t2 = trunc * trunc
        for j in range(n):
            df = np.where(free, FLT(d), -np.inf)
            p = int(np.argmax(df))
            if j == 0:
                d0 = d[p]
            if not d[p] > t2 * d0:
                break
            col = Gm[:, p] - (L[:, :j] @ L[p, :j] if j else 0)
            col = MID(col / d[p].sqrt())
            col[~free] = arb(0)
            col[p] = d[p].sqrt().mid()
            L[:, j] = col
            free[p] = False
            S.append(p)
            d = MID(d - col * col)
        k = len(S)
        R = zeros(k, k)
        for i in range(k):
            R[: i + 1, i] = L[S[i], : i + 1]
        Rinv = MID(arr(amat(R).inv()))
        return S, Rinv
    finally:
        flint.ctx.prec = with_prec


def socp(h, G, q, y, early, tol, maxit=MAXIT):
    """socp3_solve() of kernel/util/socp3.c: min q^T x subject to
    s = h - G x in (Q^3)^M, from x = 0 and the dual warm start y (M, 3),
    overwritten. h (M, 3), G (M, 3, n). Returns x, status, it, pobj, dobj,
    pres, dres."""
    M, _, n = G.shape
    e = eps()
    Gf = G.reshape(3 * M, n)
    Gmat = amat(Gf)
    GmatT = Gmat.transpose()
    hf = h.ravel()
    hn = sum(hf * hf).sqrt()
    qn = sum(q * q).sqrt()
    x = zeros(n)
    s = h.copy()
    marg = min(s[:, 0] - SQRT(s[:, 1] ** 2 + s[:, 2] ** 2))
    if marg <= 0:
        s[:, 0] = s[:, 0] + 1 - marg
    ymax = max(abs(v) for v in y[:, 0])
    shift = 0.1 * ymax - min(y[:, 0] - SQRT(y[:, 1] ** 2 + y[:, 2] ** 2))
    y[:, 0] = MID(y[:, 0] + max(shift, arb(0)) + 1e-3 * ymax)
    dres_min = None
    k = 0
    while True:
        rx = q + arr(GmatT * amat(y.ravel()))[:, 0]
        rz = arr(Gmat * amat(x))[:, 0].reshape(M, 3) + s - h
        gap = sum((s * y).ravel())
        mu = gap / M
        pobj = sum(q * x)
        dobj = -sum(hf * y.ravel())
        pres = sum((rz * rz).ravel()).sqrt() / (1 + hn)
        dres = sum(rx * rx).sqrt() / (1 + qn)
        slack = tol * (1 + abs(pobj))
        if pres < tol and dres < tol and gap < slack:
            status = "optimal"
            break
        if k == 0:
            dres_min = dres
        if dres > 100 * dres_min and dres_min < 10 * tol and gap < slack:
            status = "stalled"
            break
        dres_min = min(dres_min, dres)
        if early and dres < tol and dobj > 10 * tol:
            status = "infeasible"
            break
        if early and pres < tol and pobj < -10 * tol:
            status = "feasible"
            break
        if k == maxit:
            status = "maxit"
            break
        W, Winv = nt_scaling(s, y)
        lam = mv3(W, y)
        Wrz = mv3(Winv, rz)
        Gt = np.stack(
            [sum(Winv[r][c][:, None] * G[:, c, :] for c in range(3)) for r in range(3)],
            axis=1,
        )
        Gtm = amat(MID(Gt.reshape(3 * M, n)))
        del Gt
        GtmT = Gtm.transpose()
        Mn = GtmT * Gtm
        dmax = max(Mn[i, i] for i in range(n))
        reg = 100 * e * dmax

        for _ in range(8):
            Mreg = Mn + arb_mat(
                n, n, [reg if i == j else 0 for i in range(n) for j in range(n)]
            )
            try:
                nw = (Mreg, Gtm, GtmT, Gmat, Winv, rx, rz, Wrz)
                dx, dy, ds = newton(nw, -lam)
                break
            except ZeroDivisionError:
                reg *= 100
        else:
            status = "numerical"
            break
        if not (all(FIN(dx)) and all(FIN(dy.ravel())) and all(FIN(ds.ravel()))):
            status = "numerical"
            break
        a = min(step_to_boundary(s, ds), step_to_boundary(y, dy), arb(1))
        mu_aff = sum(((s + a * ds) * (y + a * dy)).ravel()) / M
        sig = min(max(mu_aff / mu, arb(0)), arb(1)) ** 3
        l0 = lam[:, 0]
        if not all(l0 > 0) or not all(
            l0 * l0 - lam[:, 1] ** 2 - lam[:, 2] ** 2 > e * l0 * l0
        ):
            status = "numerical"
            break
        ws, wy = mv3(Winv, ds), mv3(W, dy)
        wv = np.stack(
            [
                sig * mu - (ws * wy).sum(axis=1),
                -(ws[:, 0] * wy[:, 1] + wy[:, 0] * ws[:, 1]),
                -(ws[:, 0] * wy[:, 2] + wy[:, 0] * ws[:, 2]),
            ],
            axis=1,
        )
        dx, dy, ds = newton(nw, jsolve(lam, wv) - lam)
        if not (all(FIN(dx)) and all(FIN(dy.ravel())) and all(FIN(ds.ravel()))):
            status = "numerical"
            break
        a = min(0.99 * min(step_to_boundary(s, ds), step_to_boundary(y, dy)), arb(1))
        for _ in range(40):
            if interior(s + a * ds) and interior(y + a * dy):
                break
            a = a / 2
        x = MID(x + a * dx)
        s = MID(s + a * ds)
        y[:] = MID(y + a * dy)
        k += 1
    return x, status, k, pobj, dobj, pres, dres


J = (1, -1, -1)


def newton(nw, t):
    """Newton direction of socp() for the scaled complementarity right-hand
    side t (M, 3)."""
    Mreg, Gtm, GtmT, Gmat, Winv, rx, rz, Wrz = nw
    M = len(t)
    t = t + Wrz
    rhs = -rx - arr(GtmT * amat(t.ravel()))[:, 0]
    dx = MID(arr(Mreg.solve(amat(rhs), algorithm="approx"))[:, 0])
    t = t + arr(Gtm * amat(dx))[:, 0].reshape(M, 3)
    ds = -rz - arr(Gmat * amat(dx))[:, 0].reshape(M, 3)
    return dx, MID(mv3(Winv, t)), MID(ds)


def nt_scaling(s, y):
    """Nesterov-Todd scaling per cone: W, Winv as 3 x 3 nested lists of
    (M,) arrays."""
    sJs = s[:, 0] ** 2 - s[:, 1] ** 2 - s[:, 2] ** 2
    yJy = y[:, 0] ** 2 - y[:, 1] ** 2 - y[:, 2] ** 2
    rs, ry = SQRT(sJs), SQRT(yJy)
    beta = SQRT(SQRT(sJs / yJy))
    gamma = SQRT((1 + (s * y).sum(axis=1) / (rs * ry)) / 2)
    v = [(s[:, r] / rs + J[r] * y[:, r] / ry) / (2 * gamma) for r in range(3)]
    nv = SQRT(2 * (v[0] + 1))
    v[0] = v[0] + 1
    v = [vr / nv for vr in v]
    W = [
        [beta * (2 * v[r] * v[c] - (J[r] if r == c else 0)) for c in range(3)]
        for r in range(3)
    ]
    Wi = [
        [
            (2 * J[r] * v[r] * J[c] * v[c] - (J[r] if r == c else 0)) / beta
            for c in range(3)
        ]
        for r in range(3)
    ]
    return W, Wi


def mv3(A, u):
    return np.stack(
        [sum(A[r][c] * u[:, c] for c in range(3)) for r in range(3)], axis=1
    )


def jsolve(lam, w):
    l0, l1, l2 = lam[:, 0], lam[:, 1], lam[:, 2]
    det = l0 - (l1 * l1 + l2 * l2) / l0
    x0 = (w[:, 0] - (l1 * w[:, 1] + l2 * w[:, 2]) / l0) / det
    return np.stack([x0, (w[:, 1] - l1 * x0) / l0, (w[:, 2] - l2 * x0) / l0], axis=1)


def interior(u):
    return all(u[:, 0] > 0) and all(u[:, 0] ** 2 > u[:, 1] ** 2 + u[:, 2] ** 2)


def step_to_boundary(u, du):
    """Largest alpha, capped at 2, with u + alpha du in every cone."""
    alpha = arb(2)
    for (s0, s1, s2), (d0, d1, d2) in zip(u, du):
        a = (d0 * d0 - d1 * d1 - d2 * d2).mid()
        b = (2 * (s0 * d0 - s1 * d1 - s2 * d2)).mid()
        c = (s0 * s0 - s1 * s1 - s2 * s2).mid()
        D = (b * b - 4 * a * c).mid()
        if D >= 0:
            sq = D.sqrt()
            qq = (-(b + (sq if b >= 0 else -sq)) / 2).mid()
            if a != 0 and 0 < qq / a < alpha:
                alpha = (qq / a).mid()
            if qq != 0 and 0 < c / qq < alpha:
                alpha = (c / qq).mid()
        if a == 0 and b < 0 and -c / b < alpha:
            alpha = (-c / b).mid()
        if d0 < 0 and -s0 / d0 < alpha:
            alpha = (-s0 / d0).mid()
    return alpha


class Design:
    """Points with their rows reduced to nr unknowns by the elimination of
    the unknown piv, and G c, B c of the interpolant (gc, bc)."""

    def __init__(self, cell, x0, piv, el):
        self.cell, self.x0, self.piv, self.el = cell, x0, piv, el
        self.keep = [u for u in range(len(x0)) if u != piv]
        self.xi = np.empty(0, dtype=object)
        self.t = np.empty(0, dtype=object)
        self.Gr = self.Gi = self.B = None
        self.gr = self.gi = self.bc = np.empty(0, dtype=object)

    def add(self, xi, t):
        Gr, Gi, B = self.cell.rows(
            np.array(xi, dtype=object), np.array(t, dtype=object)
        )
        red = [
            MID(a[:, self.keep] + a[:, [self.piv]] * self.el[None, :])
            for a in (Gr, Gi, B)
        ]
        base = [MID(mm(a, self.x0)) for a in (Gr, Gi, B)]
        cat = np.concatenate
        self.xi, self.t = cat([self.xi, xi]), cat([self.t, t])
        if self.Gr is None:
            self.Gr, self.Gi, self.B = red
        else:
            self.Gr, self.Gi, self.B = (
                np.vstack([o, r]) for o, r in zip((self.Gr, self.Gi, self.B), red)
            )
        self.gr, self.gi, self.bc = (
            cat([o, n]) for o, n in zip((self.gr, self.gi, self.bc), base)
        )

    def values(self, d):
        """Re, Im of G c and B c at the window x0 + reduced change d."""
        return (
            MID(self.gr + mm(self.Gr, d)),
            MID(self.gi + mm(self.Gi, d)),
            MID(self.bc + mm(self.B, d)),
        )

    def ratio(self, d):
        gr, gi, bc = self.values(d)
        if not all(bc > 0):
            return None, (gr, gi, bc)
        r = SQRT(gr * gr + gi * gi) / bc
        return max(r).mid(), (gr, gi, bc)

    def full(self, d):
        """The full unknowns x0 + change."""
        x = self.x0.copy()
        x[self.keep] = x[self.keep] + d
        x[self.piv] = x[self.piv] + sum(self.el * d)
        return MID(x)


def driver(ds, dtot, lo, y, trunc, stats):
    """Safeguarded Dinkelbach steps and bisection on the design from dtot,
    driver() of kbpoly.c. Returns hi, dtot, lo."""
    hi, (gr, gi, bc) = ds.ratio(dtot)
    if hi is None:
        return None, dtot, lo
    # Preconditioner from the rows over C at the start; the pivots and R
    # depend on the denominators only through a relative change of the size
    # of the error, so they are kept for the call.
    A = np.vstack([ds.Gr / bc[:, None], ds.Gi / bc[:, None]])
    S, Rinv = pivchol(MID(A), trunc)
    k = len(S)
    flint.ctx.prec = PREC_QR
    Rm = amat(Rinv)
    GR, GI, BR = (arr(amat(a[:, S]) * Rm) for a in (ds.Gr, ds.Gi, ds.B))
    flint.ctx.prec = PREC
    dink, fresh = True, False
    npt = len(ds.xi)
    for _ in range(SOLVES):
        if not hi > (1 + DELTA) * lo:
            break
        tau = hi if dink else ((lo * hi).sqrt().mid() if lo > 0 else hi / 4)
        if not fresh:
            kappa = max(SQRT(gr * gr + gi * gi) / bc).mid()
            E0r, E0i = MID(gr / (bc * kappa)), MID(gi / (bc * kappa))
            Gam_r, Gam_i = MID(GR / bc[:, None]), MID(GI / bc[:, None])
            beta = MID(BR * kappa / bc[:, None])
            fresh = True
        tn = (tau / kappa).mid()
        h = np.stack([np.full(npt, tn, dtype=object), E0r, E0i], axis=1)
        G = zeros(npt, 3, k + 1)
        G[:, 0, :k] = MID(-tn * beta)
        G[:, 0, k] = arb(-1)
        G[:, 1, :k] = -Gam_r
        G[:, 2, :k] = -Gam_i
        q = zeros(k + 1)
        q[k] = arb(1)
        tol = TOL_DINK if dink else TOL_BISECT
        x, st, it, _, dobj, _, dres = socp(h, G, q, y, not dink, tol)
        stats["nit"] += it + 1
        flint.ctx.prec = PREC_QR
        step = arr(Rm * amat(x[:k]))[:, 0]
        dnew = dtot.copy()
        dnew[S] = dnew[S] + kappa * step
        flint.ctx.prec = PREC
        dnew = MID(dnew)
        r, vals = ds.ratio(dnew)
        if r is not None and r < tau:
            dink = r * GAIN < hi
            hi, dtot, fresh = r, dnew, False
            gr, gi, bc = vals
        else:
            if st == "infeasible" or tau < hi or (dres < tol and dobj > 10 * tol):
                lo = max(lo, tau)
            dink = False
    return hi, dtot, lo


def solve_cell(sigma, m, h, verbose=False):
    """Optimized table of one cell: dict with e (fine-grid maximum with
    refinement), coef (monomial, row-major nk x w), pts, nit."""
    flint.ctx.prec = PREC
    t0 = time.time()
    cell = Cell(sigma, m, h)
    nc, nk = cell.nc, cell.nk
    nu = nc * nk
    ch = cell.interpolant()
    x0 = np.array([ch[u % nk, 1 + u // nk] for u in range(nu)], dtype=object)
    trunc = 500 * eps()

    # kbpoly.c's 24 frequencies for m <= 9; above, 3 per column.
    nxi = NXI if m <= 9 else 3 * (m + 1)
    th = [arb(i) / (nxi - 1) for i in range(nxi)]
    xi = [(cell.ximax * (1 - v.cos_pi()) / 2).mid() for v in th]
    xw = [
        (0.5 if i in (0, nxi - 1) else 1) * (v.sin_pi() + 0.02)
        for i, v in enumerate(th)
    ]
    # Offsets 0 and nt / 2 midpoints: kbpoly.c's 17 for h <= 8; above, nt =
    # 2 h keeps nt / 2 + 1 >= h + 1 offsets per T_k(2 s) half-range.
    nt = max(NT, 2 * h)
    ts = [arb(0)] + [(q - arb(1) / 2) / nt for q in range(1, nt // 2 + 1)]

    # The scale: sum_i B_i d = 0 over the initial design eliminates the
    # unknown with the largest |Bsum|.
    xis = np.array([v for v in xi for _ in ts], dtype=object)
    tss = np.array([v for _ in xi for v in ts], dtype=object)
    wts = np.array([v for v in xw for _ in ts], dtype=object)
    Bsum = cell.rows(np.array(xi, dtype=object), np.zeros(nxi, dtype=object) + arb(0))[
        2
    ].sum(axis=0) * len(ts)
    piv = int(np.argmax(np.abs(FLT(Bsum))))
    el = MID(
        np.array([-Bsum[u] / Bsum[piv] for u in range(nu) if u != piv], dtype=object)
    )
    ds = Design(cell, x0, piv, el)
    ds.add(xis, tss)
    nr = nu - 1

    # One weighted LS pass, columns scaled to unit norm.
    sc = np.array(
        [1 / v.sqrt() if v > 0 else arb(1) for v in (ds.Gr**2 + ds.Gi**2).sum(axis=0)],
        dtype=object,
    )
    wa = np.array([w.sqrt() for w in wts], dtype=object) / ds.bc
    A = MID(np.vstack([ds.Gr * wa[:, None], ds.Gi * wa[:, None]]) * sc[None, :])
    yv = MID(np.concatenate([-wa * ds.gr, -wa * ds.gi]))
    S, Rinv = pivchol(A, trunc)
    flint.ctx.prec = PREC_QR
    Rm = amat(Rinv)
    g = arr(Rm * (Rm.transpose() * (amat(A[:, S]).transpose() * amat(yv))))[:, 0]
    flint.ctx.prec = PREC
    dtot = zeros(nr)
    dtot[S] = MID(g * sc[S])
    if ds.ratio(dtot)[0] is None:
        dtot = zeros(nr)

    y = np.array([[arb(1), arb(0), arb(0)] for _ in range(len(ds.xi))], dtype=object)
    stats = {"nit": 0}
    lo, fbest, dbest = arb(0), None, dtot
    floor = 500 * arb(2) ** -PREC
    for r in range(ROUNDS):
        mx, dtot, lo = driver(ds, dtot, lo, y, trunc, stats)
        fine, pts = cell.worst(cell.mono(table(ch, ds.full(dtot), nk)))
        if verbose:
            v = [f"{float(a):.6e}" if a is not None else "-" for a in (mx, lo, fine)]
            print(
                f"  round {r} npt {len(ds.xi)} hi {v[0]} lo {v[1]} fine {v[2]}",
                flush=True,
            )
        if fine is None:
            break
        # The refined maximum: e of the header and the gap e / lo.
        top = pts[0][0]
        if fbest is None or top < fbest:
            fbest, dbest, best = top, dtot, pts
        # Against lo, not hi: failed() tests the same e / lo.
        if mx is None or top <= (1 + DELTA_X) * lo or (r >= 1 and mx < floor):
            break
        ds.add(
            np.array([p[1] for p in pts], dtype=object),
            np.array([p[2] for p in pts], dtype=object),
        )
        y = np.vstack(
            [y, np.array([[arb(1), arb(0), arb(0)] for _ in pts], dtype=object)]
        )

    if fbest is None:
        raise SystemExit(f"cell {sigma} {m} {h}: psi_hat <= 0 off the design")
    mono = cell.mono(table(ch, ds.full(dbest), nk))
    _, (gr, gi, bc) = ds.ratio(dbest)
    ep = FLT(SQRT(gr * gr + gi * gi) / bc)
    return {
        "sigma": sigma,
        "m": m,
        "h": h,
        "e": best[0][0].mid(),
        "coef": list(mono.ravel()),
        "pts": points(ds.xi, ds.t, ep, cell.ximax),
        "nit": stats["nit"],
        "lo": float(lo),
        "rounds": r + 1,
        "seconds": time.time() - t0,
    }


def table(ch, x, nk):
    """Chebyshev table (nk, w) with columns 1 .. m + 1 from the unknowns."""
    out = ch.copy()
    for u, v in enumerate(x):
        out[u % nk, 1 + u // nk] = v
    return out


def points(xi, t, ev, ximax):
    """(round(65535 xi / ximax), round(65535 t)), one per key with its
    largest |E|, the NPTS largest, sorted."""
    best = {}
    for x, s, v in zip(xi, t, ev):
        key = (round(65535 * float(x / ximax)), round(65535 * float(s)))
        best[key] = max(best.get(key, -1.0), v)
    top = sorted(best, key=lambda k: -best[k])[:NPTS]
    return sorted(top)


def fmt(v):
    """Decimal text of an arb, DIGITS significant digits; text passes."""
    if isinstance(v, str):
        return v
    s = v.mid().str(DIGITS, radius=False)
    return s if any(c in s for c in ".e") else s + ".0"


def check(res):
    """Re-evaluate the fine-grid maximum from the emitted digits."""
    flint.ctx.prec = PREC
    cell = Cell(res["sigma"], res["m"], res["h"])
    mono = np.array([arb(fmt(v)) for v in res["coef"]], dtype=object).reshape(
        cell.nk, cell.w
    )
    _, best = cell.worst(mono)
    moved = abs(float(best[0][0] / res["e"]) - 1)
    if moved > 1e-3:
        raise SystemExit(
            f"sigma={res['sigma']} m={res['m']}: e moved by {moved:.2e} from the emitted digits"
        )
    return moved


def cell_text(res, i):
    m, coef, pts = res["m"], res["coef"], res["pts"]
    name = f"{m}_{i}"
    vals = [f"K({fmt(v)})" for v in coef]
    out = [f"static const R kb_poly_opt_coef_{name}[] = {{"]
    out += ["  " + ", ".join(vals[k : k + 2]) + "," for k in range(0, len(vals), 2)]
    out.append("};")
    flat = [str(v) for p in pts for v in p]
    out.append(f"static const unsigned short kb_poly_opt_pts_{name}[] = {{")
    out += ["  " + ", ".join(flat[k : k + 12]) + "," for k in range(0, len(flat), 12)]
    out.append("};")
    return out


def header(results):
    """The header text of the results, sorted by (m, sigma)."""
    rows = {}
    for res in sorted(results, key=lambda r: (r["m"], r["sigma"])):
        rows.setdefault(res["m"], []).append(res)
    # A partial header (--mmax) caps every MMAX at its largest m.
    top = max(rows, default=0)
    mmax = {
        p: min(top, math.floor(generate.br_ls(p) / bval(1.1)) - 1)
        for p in (24, 53, 64, 113)
    }
    out = [
        "/* Generated by gen/kbpolygen/optima.py; do not edit. */",
        "#if MANT_DIG == 24",
        f"  #define KB_POLY_OPT_MMAX {mmax[24]}",
        "#elif MANT_DIG == 53",
        f"  #define KB_POLY_OPT_MMAX {mmax[53]}",
        "#elif MANT_DIG == 64",
        f"  #define KB_POLY_OPT_MMAX {mmax[64]}",
        "#else",
        f"  #define KB_POLY_OPT_MMAX {mmax[113]}",
        "#endif",
    ]
    for m, cells in rows.items():
        out.append(f"#if KB_POLY_OPT_MMAX >= {m}")
        for i, res in enumerate(cells):
            out += cell_text(res, i)
        out.append(f"static const kb_poly_opt_cell kb_poly_opt_row_{m}[] = {{")
        for i, res in enumerate(cells):
            out.append(
                f"  {{ K({float(res['sigma'])!r}), K({fmt(res['e'])}), {m}, {res['h']}, "
                f"{len(res['pts'])}, kb_poly_opt_coef_{m}_{i}, kb_poly_opt_pts_{m}_{i} }},"
            )
        out += ["};", "#endif"]
    out.append("static const kb_poly_opt_row kb_poly_opt_rows[] = {")
    for m, cells in rows.items():
        out += [
            f"#if KB_POLY_OPT_MMAX >= {m}",
            f"  {{ {m}, {len(cells)}, kb_poly_opt_row_{m} }},",
            "#endif",
        ]
    out += [
        "};",
        "#define KB_POLY_OPT_NROWS (sizeof(kb_poly_opt_rows) / sizeof(kb_poly_opt_rows[0]))",
    ]
    return "\n".join(out) + "\n"


def bval(sigma):
    return math.pi * (2 - 1 / sigma)


def cells(mmax=None):
    """(sigma, m) of the table, largest m first."""
    out = []
    for s in SIGMAS:
        top = math.floor(generate.br_ls(113) / bval(s)) - 1
        out += [(s, m) for m in range(1, top + 1) if mmax is None or m <= mmax]
    return sorted(out, key=lambda c: (-c[1], c[0]))


def _init():
    generate._init()
    flint.ctx.prec = PREC


def run(job):
    sigma, m = job
    h = generate.halves(sigma, m)[113]
    res = solve_cell(sigma, m, h)
    check(res)
    # arb does not pickle: send the digits to the parent.
    res["e"], res["coef"] = fmt(res["e"]), [fmt(v) for v in res["coef"]]
    return res


def failed(res):
    """The exchange left a gap e / lo - 1 above DELTA_X."""
    return float(res["e"]) > (1 + DELTA_X) * res["lo"]


def logline(res):
    return (
        f"{res['sigma']:<5} {res['m']:2d} {res['h']:2d} {float(res['e']):.4e} "
        f"{res['nit']:5d} {res['seconds']:8.1f} {res['lo']:.4e} {res['rounds']:2d}"
        + (" FAIL" if failed(res) else "")
    )


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("-j", type=int, default=os.cpu_count())
    ap.add_argument("--cell", nargs=3, metavar=("SIGMA", "M", "H"))
    ap.add_argument("--mmax", type=int)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args(argv)

    if args.cell:
        _init()
        s, m, h = float(args.cell[0]), int(args.cell[1]), int(args.cell[2])
        res = solve_cell(s, m, h, verbose=True)
        print(f"e {float(res['e']):.6e} nit {res['nit']} seconds {res['seconds']:.1f}")
        return

    results = []
    # One cell per worker process: memory of a large cell is returned.
    with Pool(args.j, initializer=_init, maxtasksperchild=1) as pool:
        for res in pool.imap_unordered(run, cells(args.mmax)):
            print(logline(res), flush=True)
            results.append(res)
    good = [r for r in results if not failed(r)]
    with open(args.out, "w") as f:
        f.write(header(good))
    sys.stderr.write(f"{len(good)} cells to {args.out}\n")
    if len(good) < len(results):
        raise SystemExit(f"{len(results) - len(good)} cells above the exchange gap")


if __name__ == "__main__":
    main()
