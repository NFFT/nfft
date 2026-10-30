"""Rounding-aware window design (goal 4 of analysis.md section 13).

Even window in develop's layout (cells j = 0..m of |nx|, Chebyshev basis
T_k(2s) in the offset s from the cell centre, degree 2h). Unknowns c[j, k].
C(xi) = psi_hat(xi) of the polynomial itself (real, linear in c).

For an error level tau the feasibility problem

    |D(xi,t) - C(xi)| + u gamma A(t) <= tau C(xi)   for all (xi, t) on the grid
    C(0) = 1

is a second-order cone program; tau is bisected. A(t) is the rounding
amplification: sum over the 2m+2 taps of |p| (the gather) plus the Horner
terms sum_k |mono_k| |s|^k (both cells and mirrors), gamma a calibration
factor against the library.

Usage: python kbopt7.py sigma m h u gamma [N]
Writes coef_rnd_s<sigma>_m<m>_h<h>.txt (with the C block for N) and
coef_alias_... (the u = 0 optimum with its own C), and prints predictions.
"""
import sys
import numpy as np
import cvxpy as cp
import mpmath as mp
from numpy.polynomial import chebyshev as Ch
from scipy.special import roots_legendre
from kbopt5 import SymModel, write_coef
from kbopt import kb_window, PI


def cheb_to_mono_matrix(nk):
    """T[k, i]: coefficient of s^i in T_k(2s)."""
    mp.mp.dps = 50
    tk = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(1)]]
    for k in range(2, nk):
        prev, pp = tk[-1], tk[-2]
        nxt = [mp.mpf(0)] * (k + 1)
        for i, v in enumerate(prev):
            nxt[i + 1] += 2 * v
        for i, v in enumerate(pp):
            nxt[i] -= v
        tk.append(nxt)
    T = np.zeros((nk, nk))
    for k in range(nk):
        for i, v in enumerate(tk[k]):
            T[k, i] = float(v * mp.mpf(2) ** i)
    return T


class Design:
    def __init__(self, M, nxi=48, nt=24):
        self.M = M
        self.xi = np.linspace(0, M.xi_max, nxi)
        self.t = (np.arange(nt) + 0.5) / nt
        self.s = self.t - 0.5
        # D rows: Phi (nxi*nt, nu) complex, D = Phi @ c.ravel()
        self.Phi = M.design(self.xi, self.t)
        # C(xi) = G @ c.ravel(), real: integral over t of the window's taps
        xg, wg = roots_legendre(M.deg + 12)
        tg, wgt = (xg + 1) / 2, wg / 2
        Phig = M.design(self.xi, tg).reshape(len(self.xi), len(tg), -1)
        self.G = np.real(np.einsum("g,xgu->xu", wgt, Phig))
        self.G0 = np.real(np.einsum("g,gu->u", wgt, M.design(np.array([0.0]), tg)[: len(tg)]))
        # values p_j(t) and p_j(1-t): V (nt, nu) each
        Vp = Ch.chebvander(2 * self.s, M.deg)
        Vm = Ch.chebvander(-2 * self.s, M.deg)
        nc = M.m + 1
        self.Vp = np.zeros((nt, nc, nc * M.nk)); self.Vm = np.zeros((nt, nc, nc * M.nk))
        for j in range(nc):
            self.Vp[:, j, j * M.nk:(j + 1) * M.nk] = Vp
            self.Vm[:, j, j * M.nk:(j + 1) * M.nk] = Vm
        self.T = cheb_to_mono_matrix(M.nk)
        self.spow = np.abs(self.s)[:, None] ** np.arange(M.nk)[None, :]   # (nt, nk)


def total_on_grid(Dn, c, u, gamma):
    """max over the design grid of |D - C|/C + u gamma A/C, in plain arithmetic."""
    M = Dn.M; nt = len(Dn.t)
    cv = c.ravel()
    C = Dn.G @ cv
    Dm = (Dn.Phi @ cv).reshape(len(Dn.xi), nt)
    P = np.einsum("tju,u->tj", Dn.Vp, cv); Pm = np.einsum("tju,u->tj", Dn.Vm, cv)
    mono = np.abs(c @ Dn.T)
    A = np.abs(P).sum(axis=1) + np.abs(Pm).sum(axis=1) + 2 * (Dn.spow @ mono.sum(axis=0))
    tot = np.abs(Dm - C[:, None]) / np.abs(C)[:, None] + u * gamma * A[None, :] / np.abs(C)[:, None]
    return tot.max()


def feasible(Dn, tau, u, gamma):
    """Scaled by 1/tau so the solver sees O(1) quantities."""
    M = Dn.M; nc = M.m + 1; nu = nc * M.nk
    c = cp.Variable(nu)
    q = cp.Variable((len(Dn.t), nc)); qm = cp.Variable((len(Dn.t), nc))
    r = cp.Variable((nc, M.nk))
    cons = [Dn.G0 @ c == 1]
    for j in range(nc):
        mono = Dn.T.T @ c[j * M.nk:(j + 1) * M.nk]
        cons += [r[j, :] >= mono, r[j, :] >= -mono]
    P = cp.vstack([Dn.Vp[i] @ c for i in range(len(Dn.t))])
    Pm = cp.vstack([Dn.Vm[i] @ c for i in range(len(Dn.t))])
    cons += [q >= P, q >= -P, qm >= Pm, qm >= -Pm]
    A = cp.sum(q, axis=1) + cp.sum(qm, axis=1) + 2 * (Dn.spow @ cp.sum(r, axis=0))
    nxi, nt = len(Dn.xi), len(Dn.t)
    Gx = np.repeat(Dn.G, nt, axis=0)
    Re = (np.real(Dn.Phi) - Gx) / tau; Im = np.imag(Dn.Phi) / tau
    rhs = Gx @ c - (u * gamma / tau) * cp.hstack([A] * nxi)
    cons += [cp.norm(cp.vstack([Re @ c, Im @ c]), 2, axis=0) <= rhs]
    prob = cp.Problem(cp.Minimize(0), cons)
    try:
        prob.solve(solver=cp.CLARABEL, max_iter=200)
    except Exception:
        return None
    if prob.status != "optimal" or c.value is None:
        return None
    cc = c.value.reshape(nc, M.nk)
    if total_on_grid(Dn, cc, u, gamma) > 1.05 * tau:
        return None
    return cc


def minimize_tau(Dn, u, gamma, c_start, tau0, outer=5):
    """Minimize tau with C(xi) frozen at the previous iterate (sequential
    convex), scaled by tau0; returns the best verified (c, total)."""
    M = Dn.M; nc = M.m + 1; nu = nc * M.nk
    nxi, nt = len(Dn.xi), len(Dn.t)
    c_prev = c_start
    best_c, best_t = c_start, total_on_grid(Dn, c_start, u, gamma)
    for it in range(outer):
        Cprev = np.repeat(Dn.G @ c_prev.ravel(), nt)
        c = cp.Variable(nu); tau = cp.Variable(nonneg=True)
        q = cp.Variable((nt, nc)); qm = cp.Variable((nt, nc)); r = cp.Variable((nc, M.nk))
        cons = [Dn.G0 @ c == 1]
        for j in range(nc):
            mono = Dn.T.T @ c[j * M.nk:(j + 1) * M.nk]
            cons += [r[j, :] >= mono, r[j, :] >= -mono]
        P = cp.vstack([Dn.Vp[i] @ c for i in range(nt)])
        Pm = cp.vstack([Dn.Vm[i] @ c for i in range(nt)])
        cons += [q >= P, q >= -P, qm >= Pm, qm >= -Pm]
        A = cp.sum(q, axis=1) + cp.sum(qm, axis=1) + 2 * (Dn.spow @ cp.sum(r, axis=0))
        Gx = np.repeat(Dn.G, nt, axis=0)
        Re = (np.real(Dn.Phi) - Gx) / tau0; Im = np.imag(Dn.Phi) / tau0
        rhs = cp.multiply(tau, Cprev) - (u * gamma / tau0) * cp.hstack([A] * nxi)
        cons += [cp.norm(cp.vstack([Re @ c, Im @ c]), 2, axis=0) <= rhs]
        prob = cp.Problem(cp.Minimize(tau), cons)
        try:
            prob.solve(solver=cp.CLARABEL, max_iter=300)
            if c.value is None:
                raise RuntimeError(prob.status)
        except Exception as e:
            print("   clarabel failed (%s), trying SCS" % e)
            try:
                prob.solve(solver=cp.SCS, eps_abs=1e-9, eps_rel=1e-9, max_iters=50000)
            except Exception as e2:
                print("   scs failed:", e2); break
        if c.value is None:
            print("   status", prob.status); break
        cc = c.value.reshape(nc, M.nk)
        tot = total_on_grid(Dn, cc, u, gamma)
        print("   outer %d: status %s  tau=%.3e  verified total=%.3e" % (it, prob.status, tau.value * tau0, tot))
        if tot < best_t:
            best_c, best_t = cc, tot
        if abs(tot - best_t) <= 0.01 * best_t and it > 0 and tot <= best_t * 1.001:
            pass
        c_prev = cc
    return best_c, best_t


def predict(M, c, u, gamma, C_fn=None, nfine=200):
    """(aliasing J_inf, predicted total with rounding, 1/C at band edge)."""
    xi = np.linspace(-M.xi_max, M.xi_max, nfine); t = (np.arange(nfine) + 0.5) / nfine
    D = M.D_of(c, xi, t)
    if C_fn is None:
        xg, wg = roots_legendre(M.deg + 12); tg, wgt = (xg + 1) / 2, wg / 2
        C = np.real(np.einsum("g,xg->x", wgt, M.D_of(c, xi, tg)))
    else:
        C = C_fn(xi)
    E = np.abs(D / C[:, None] - 1)
    P, Pm = M.values(c, t)
    T = cheb_to_mono_matrix(M.nk)
    mono = np.abs(c @ T)                                   # (nc, nk)
    spow = np.abs(t - 0.5)[:, None] ** np.arange(M.nk)[None, :]
    A = np.abs(P).sum(axis=1) + np.abs(Pm).sum(axis=1) + 2 * (spow @ mono.sum(axis=0))
    tot = E + u * gamma * A[None, :] / np.abs(C)[:, None]
    return E.max(), tot.max(), 1 / np.abs(C).min(), C


def main():
    sigma, m, h = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    u, gamma = float(sys.argv[4]), float(sys.argv[5])
    N = int(sys.argv[6]) if len(sys.argv) > 6 else 2048
    n = int(round(sigma * N)); n += n % 2
    M = SymModel(sigma, m, h)
    Dn = Design(M)
    c0 = M.cheb_kb()
    J0, T0, amp0, _ = predict(M, c0, u, gamma, C_fn=M.C)
    print("=== sigma=%g m=%d deg=%d u=%.1e gamma=%g" % (sigma, m, M.deg, u, gamma))
    print("%-10s aliasing J_inf=%.3e  predicted total=%.3e  1/C(edge)=%.2f" % ("stock", J0, T0, amp0))
    T0g = total_on_grid(Dn, c0, u, gamma)
    ca, ta = minimize_tau(Dn, 0.0, gamma, c0, T0g)
    Ja, Ta, ampa, Ca = predict(M, ca, u, gamma)
    print("%-10s aliasing J_inf=%.3e  predicted total=%.3e  1/C(edge)=%.2f  (tau=%.2e)" % ("alias-opt", Ja, Ta, ampa, ta))
    cr, tr = minimize_tau(Dn, u, gamma, ca, T0g)
    Jr, Tr, ampr, Cr = predict(M, cr, u, gamma)
    print("%-10s aliasing J_inf=%.3e  predicted total=%.3e  1/C(edge)=%.2f  (tau=%.2e)" % ("round-opt", Jr, Tr, ampr, tr))
    ks = np.arange(-N // 2, N // 2)
    xg, wg = roots_legendre(M.deg + 12); tg, wgt = (xg + 1) / 2, wg / 2
    for name, c in (("ownc", c0), ("alias", ca), ("rnd", cr)):
        Ck = np.real(np.einsum("g,xg->x", wgt, M.D_of(c, ks / n, tg)))
        path = "coef_%s_s%g_m%d_h%d.txt" % (name, sigma, m, h)
        write_coef(path, M, M.to_monomial(c))
        with open(path, "a") as f:
            f.write("C %d\n" % N)
            f.write("\n".join("%.17e" % v for v in Ck) + "\n")
    # window shape summary
    for name, c in (("stock", c0), ("alias", ca), ("rnd", cr)):
        P, _ = M.values(c, np.array([0.5]))
        print("   %-6s centre values per cell j=0..m: %s" % (name, " ".join("%.3e" % v for v in P[0])))


if __name__ == "__main__":
    main()
