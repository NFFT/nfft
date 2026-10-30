"""Stock (Chebyshev) coefficients plus the polynomial's own deconvolution
factors psi_hat(k/n), computed in double by Gauss-Legendre: isolates the
accuracy of the PHI_HUT evaluation. Usage: kbownc.py sigma m h N"""
import sys
import numpy as np
from scipy.special import roots_legendre
from kbopt5 import SymModel, write_coef
sigma, m, h, N = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
n = int(round(sigma * N)); n += n % 2
M = SymModel(sigma, m, h); c0 = M.cheb_kb()
ks = np.arange(-N // 2, N // 2)
xg, wg = roots_legendre(M.deg + 12); tg, wgt = (xg + 1) / 2, wg / 2
Ck = np.real(np.einsum("g,xg->x", wgt, M.D_of(c0, ks / n, tg)))
path = "coef_ownc_s%g_m%d_h%d.txt" % (sigma, m, h)
write_coef(path, M, M.to_monomial(c0))
with open(path, "a") as f:
    f.write("C %d\n" % N); f.write("\n".join("%.17e" % v for v in Ck) + "\n")
print("wrote", path, " C(edge)=%.4f  I0-formula C(edge)=%.4f" % (Ck[0], M.C(np.array([ks[0] / n]))[0]))
