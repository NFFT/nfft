"""Rounding-aware window by a smooth NLP (SLSQP, epigraph form).

minimize tau  s.t.  [ |D(xi,t) - C(xi)| + u gamma A(t) ] / C(xi) <= tau  on the grid,  C(0) = 1,

|.| smoothed as sqrt(x^2 + eps^2). Variables: the (m+1)(2h+1) Chebyshev
coefficients and tau, all scaled by tau0 (the stock total). Start: the
aliasing optimum of kbopt7 if its file exists, else the Chebyshev interpolant.

Usage: python kbopt8.py sigma m h u gamma [N]
"""
import sys, os
import numpy as np
from scipy.optimize import minimize
from scipy.special import roots_legendre
from kbopt5 import SymModel, write_coef
from kbopt7 import Design, cheb_to_mono_matrix, total_on_grid, predict


def read_coef_cheb(path, M):
    """Chebyshev coefficients from a monomial coefficient file (inverse of to_monomial)."""
    with open(path) as f:
        deg, m = map(int, f.readline().split())
        rows = np.array([[float(v) for v in f.readline().split()] for _ in range(deg + 1)])
    T = cheb_to_mono_matrix(M.nk)             # mono = c @ T
    c = np.zeros((M.m + 1, M.nk))
    for j in range(M.m + 1):
        col = M.m + 1 - j
        c[j] = np.linalg.solve(T.T, rows[:, col - 1])
    return c


def main():
    sigma, m, h = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    u, gamma = float(sys.argv[4]), float(sys.argv[5])
    N = int(sys.argv[6]) if len(sys.argv) > 6 else 2048
    n = int(round(sigma * N)); n += n % 2
    M = SymModel(sigma, m, h)
    Dn = Design(M, nxi=40, nt=20)
    nc, nk = M.m + 1, M.nk; nu = nc * nk
    c0 = M.cheb_kb()
    start = "coef_alias_s%g_m%d_h%d.txt" % (sigma, m, h)
    cs = read_coef_cheb(start, M) if os.path.exists(start) else c0
    tau0 = total_on_grid(Dn, c0, u, gamma)
    nxi, nt = len(Dn.xi), len(Dn.t)
    Re = np.real(Dn.Phi); Im = np.imag(Dn.Phi)
    Gx = np.repeat(Dn.G, nt, axis=0)
    ReG = Re - Gx
    Vp = Dn.Vp.reshape(nt * nc, nu); Vm = Dn.Vm.reshape(nt * nc, nu)
    # Horner term: 2 * sum_j sum_k |mono_{j,k}| |s|^k, mono = c @ T  (per cell)
    Tm = np.kron(np.eye(nc), Dn.T.T)          # mono.ravel() = Tm @ c.ravel()   (nc*nk)
    spow = Dn.spow                            # (nt, nk)
    eps_s = 1e-3 * tau0                       # smoothing of |.| on O(tau0) quantities
    tile_t = np.tile(np.arange(nt), nxi)      # row -> t index

    def unpack(z):
        return z[:nu], z[nu]

    def smooth_abs(x, e):
        return np.sqrt(x * x + e * e)

    def cons(z):
        c, tau = unpack(z)
        dr = ReG @ c; di = Im @ c
        mag = np.sqrt(dr * dr + di * di + (eps_s) ** 2)
        P = Vp @ c; Pm = Vm @ c
        mono = Tm @ c
        A = (smooth_abs(P, 1e-3).reshape(nt, nc).sum(axis=1)
             + smooth_abs(Pm, 1e-3).reshape(nt, nc).sum(axis=1)
             + 2 * (spow @ smooth_abs(mono, 1e-3).reshape(nc, nk).sum(axis=0)))
        C = Gx @ c
        g = tau - (mag + u * gamma * A[tile_t]) / C / tau0      # >= 0
        return g

    def cons_jac(z):
        c, tau = unpack(z)
        dr = ReG @ c; di = Im @ c
        mag = np.sqrt(dr * dr + di * di + eps_s ** 2)
        dmag = (dr[:, None] * ReG + di[:, None] * Im) / mag[:, None]
        P = Vp @ c; Pm = Vm @ c; mono = Tm @ c
        sP = smooth_abs(P, 1e-3); sPm = smooth_abs(Pm, 1e-3); sM = smooth_abs(mono, 1e-3)
        A = (sP.reshape(nt, nc).sum(axis=1) + sPm.reshape(nt, nc).sum(axis=1)
             + 2 * (spow @ sM.reshape(nc, nk).sum(axis=0)))
        # dA/dc (nt, nu)
        dA = ((P / sP)[:, None] * Vp).reshape(nt, nc, nu).sum(axis=1) \
           + ((Pm / sPm)[:, None] * Vm).reshape(nt, nc, nu).sum(axis=1)
        dM = ((mono / sM)[:, None] * Tm).reshape(nc, nk, nu).sum(axis=0)   # (nk, nu)
        dA = dA + 2 * (spow @ dM)
        C = Gx @ c
        num = mag + u * gamma * A[tile_t]
        dnum = dmag + u * gamma * dA[tile_t]
        J = np.zeros((len(num), nu + 1))
        J[:, :nu] = -(dnum / C[:, None] - (num / C / C)[:, None] * Gx) / tau0
        J[:, nu] = 1.0
        return J

    z0 = np.concatenate([cs.ravel(), [total_on_grid(Dn, cs, u, gamma) / tau0]])
    res = minimize(lambda z: z[nu], z0, jac=lambda z: np.r_[np.zeros(nu), 1.0], method="SLSQP",
                   constraints=[{"type": "ineq", "fun": cons, "jac": cons_jac},
                                {"type": "eq", "fun": lambda z: Dn.G0 @ z[:nu] - 1,
                                 "jac": lambda z: np.r_[Dn.G0, 0.0]}],
                   options={"maxiter": 400, "ftol": 1e-12})
    c = res.x[:nu].reshape(nc, nk)
    print("=== sigma=%g m=%d deg=%d u=%.1e gamma=%g   SLSQP: %s (%d it)" % (sigma, m, M.deg, u, gamma, res.message, res.nit))
    for name, cc, Cfn in (("stock", c0, M.C), ("start", cs, None), ("nlp", c, None)):
        J, T, amp, _ = predict(M, cc, u, gamma, C_fn=Cfn)
        print("%-6s aliasing J_inf=%.3e  predicted total=%.3e  1/C(edge)=%.2f" % (name, J, T, amp))
    ks = np.arange(-N // 2, N // 2)
    xg, wg = roots_legendre(M.deg + 12); tg, wgt = (xg + 1) / 2, wg / 2
    Ck = np.real(np.einsum("g,xg->x", wgt, M.D_of(c, ks / n, tg)))
    path = "coef_nlp_s%g_m%d_h%d.txt" % (sigma, m, h)
    write_coef(path, M, M.to_monomial(c))
    with open(path, "a") as f:
        f.write("C %d\n" % N); f.write("\n".join("%.17e" % v for v in Ck) + "\n")
    P, _ = M.values(c, np.array([0.5]))
    print("   nlp centre values per cell: %s" % " ".join("%.3e" % v for v in P[0]))


if __name__ == "__main__":
    main()
