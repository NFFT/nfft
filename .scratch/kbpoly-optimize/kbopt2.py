"""Follow-ups on kbopt.py for one (sigma, m, degree):
 1. Lawson IRLS: J_inf of the best iterate after k iterations.
 2. True discrete minimax by LP (HiGHS) with a 32-gon modulus, as a check.
 3. Window shape: values at the stencil edges and the cell-boundary jumps,
    chebKB vs LS vs MM.
"""
import sys
import numpy as np
from scipy.optimize import linprog
from kbopt import Model, kb_window

sigma, m = float(sys.argv[1]), int(sys.argv[2])
d = Model(sigma, m, 1).lut_degree(53) + int(sys.argv[3])
M = Model(sigma, m, d)
C = M.C_lib(M.xi)
Phi = M.design(M.xi, M.t)
A = Phi / np.repeat(C, len(M.t))[:, None]
rhs = np.ones(A.shape[0])

print("=== sigma=%g m=%d deg=%d (LUT double=%d)" % (sigma, m, d, M.lut_degree(53)))

# 1. Lawson convergence
w = np.ones(len(rhs)) / len(rhs)
best = np.inf
marks = {1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 100, 250}
for it in range(1, 251):
    x = M._solve(A, rhs, np.sqrt(w))
    r = np.abs(A @ x - rhs)
    if r.max() < best:
        best, bx = r.max(), x
    w = w * r
    w /= w.sum()
    if it in marks:
        Jinf, J2 = M.metrics(bx.reshape(M.ncell, d + 1), M.C_lib)
        print("lawson it=%3d  design-grid max=%.3e   fine-grid J_inf=%.3e J_2=%.3e"
              % (it, best, Jinf, J2))
c_mm = bx.reshape(M.ncell, d + 1)

# 2. LP minimax: min tau s.t. Re(e^{-i th} (A x - 1)) <= tau for 32 angles
nx = A.shape[1]
ths = np.arange(32) * 2 * np.pi / 32
rows, rhs_lp = [], []
for th in ths:
    e = np.exp(-1j * th)
    G = np.real(e * A)
    rows.append(np.hstack([G, -np.ones((G.shape[0], 1))]))
    rhs_lp.append(np.real(e * rhs))
G = np.vstack(rows)
h = np.concatenate(rhs_lp)
cost = np.zeros(nx + 1)
cost[-1] = 1
res = linprog(cost, A_ub=G, b_ub=h, bounds=[(None, None)] * (nx + 1), method="highs")
if res.success:
    x_lp = res.x[:-1]
    c_lp = x_lp.reshape(M.ncell, d + 1)
    print("LP minimax (32-gon): design-grid tau=%.3e  fine-grid J_inf=%.3e J_2=%.3e"
          % (res.x[-1], *M.metrics(c_lp, M.C_lib)))
else:
    print("LP failed:", res.message)

# 3. Window shapes
c_kb = M.cheb_kb()
c_ls = M.fit_ls(C)
peak = kb_window(np.array([0.0]), M.b, M.R)[0]
def shape(c, name):
    left = M.values(c, np.array([0.0]))[0]          # p_a(0)
    right = M.values(c, np.array([1.0]))[0]         # p_a(1)
    # jumps between cell a and a+1: p_{a+1}(0) - p_a(1)
    jumps = left[1:] - right[:-1]
    print("%-7s edge p(-R)=%+.2e  p(R^-)=%+.2e  max interior jump=%.2e  (rel. to peak)"
          % (name, left[0] / peak, right[-1] / peak, np.abs(jumps).max() / peak))
    return left, right
shape(c_kb, "chebKB")
shape(c_ls, "LS")
shape(c_mm, "MM")
# true KB at the edge
print("KB      edge phi(R)/peak = %.2e, jumps are the Chebyshev truncation error only"
      % (kb_window(np.array([M.R - 1e-12]), M.b, M.R)[0] / peak))
# the difference window MM - KB on a fine grid: where does the optimizer move mass?
tt = np.linspace(0, 1, 11)[:-1]
dP = (M.values(c_mm, tt) - kb_window(tt[:, None] + M.cells[None, :], M.b, M.R)) / peak
print("max |MM - KB|/peak per cell:", " ".join("%.1e" % v for v in np.abs(dP).max(axis=0)))
