"""Rounding-aware shaping in the null-ish space of the aliasing operator.

c = c_start + V beta, V = right singular vectors of (Phi - G) with the
smallest singular values (shape directions that leave the aliasing residual
almost unchanged), beta found by Nelder-Mead on the plainly evaluated total
error max [ |D - C| + u gamma A ] / C, then a least-squares aliasing polish
in the complementary space. Usage: kbopt9.py sigma m h u gamma p [N]
"""
import sys, os
import numpy as np
from scipy.optimize import minimize
from scipy.special import roots_legendre
from kbopt5 import SymModel, write_coef
from kbopt7 import Design, total_on_grid, predict
from kbopt8 import read_coef_cheb

sigma, m, h = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
u, gamma, p = float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
N = int(sys.argv[7]) if len(sys.argv) > 7 else 2048
n = int(round(sigma * N)); n += n % 2
M = SymModel(sigma, m, h)
Dn = Design(M, nxi=48, nt=24)
nc, nk = M.m + 1, M.nk; nu = nc * nk
c0 = M.cheb_kb()
start = "coef_alias_s%g_m%d_h%d.txt" % (sigma, m, h)
cs = read_coef_cheb(start, M) if os.path.exists(start) else M.fit_ls(c0)
Gx = np.repeat(Dn.G, len(Dn.t), axis=0)
Aop = np.vstack([np.real(Dn.Phi) - Gx, np.imag(Dn.Phi)])        # aliasing operator
# remove the normalization direction: work in the space orthogonal to G0
Q, _ = np.linalg.qr(np.c_[Dn.G0, np.eye(nu)])
Qn = Q[:, 1:]                                                   # basis of {G0 . c = 0}
U, S, Vt = np.linalg.svd(Aop @ Qn, full_matrices=False)
V = Qn @ Vt[-p:].T                                              # p softest directions
print("=== sigma=%g m=%d deg=%d u=%.1e gamma=%g  p=%d  singular values (softest): %s"
      % (sigma, m, M.deg, u, gamma, p, " ".join("%.1e" % v for v in S[-p:])))
scale = np.linalg.norm(cs)

def polish(c):
    """least-squares aliasing polish keeping C(0) = 1 (small correction)."""
    r = Aop @ c.ravel()
    d, *_ = np.linalg.lstsq(Aop @ Qn, -r, rcond=None)
    return (c.ravel() + Qn @ d).reshape(nc, nk)

def total(beta):
    c = (cs.ravel() + V @ (beta * scale)).reshape(nc, nk)
    c = polish(c)
    return total_on_grid(Dn, c, u, gamma)

best = (total(np.zeros(p)), np.zeros(p))
for trial in range(4):
    x0 = best[1] + (0.05 if trial else 0.0) * np.random.default_rng(trial).standard_normal(p)
    res = minimize(total, x0, method="Nelder-Mead", options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-12})
    if res.fun < best[0]:
        best = (res.fun, res.x)
c = polish((cs.ravel() + V @ (best[1] * scale)).reshape(nc, nk))
for name, cc, Cfn in (("stock", c0, M.C), ("start", cs, None), ("shaped", c, None)):
    J, T, amp, _ = predict(M, cc, u, gamma, C_fn=Cfn)
    print("%-7s aliasing J_inf=%.3e  predicted total=%.3e  1/C(edge)=%.2f" % (name, J, T, amp))
ks = np.arange(-N // 2, N // 2)
xg, wg = roots_legendre(M.deg + 12); tg, wgt = (xg + 1) / 2, wg / 2
Ck = np.real(np.einsum("g,xg->x", wgt, M.D_of(c, ks / n, tg)))
path = "coef_shape_s%g_m%d_h%d.txt" % (sigma, m, h)
write_coef(path, M, M.to_monomial(c))
with open(path, "a") as f:
    f.write("C %d\n" % N); f.write("\n".join("%.17e" % v for v in Ck) + "\n")
P, _ = M.values(c, np.array([0.5]))
print("   shaped centre values per cell: %s" % " ".join("%.3e" % v for v in P[0]))
