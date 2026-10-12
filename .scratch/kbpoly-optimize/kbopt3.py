"""Minimax on a finer design grid, to see what the discretised Lawson/LP of
kbopt2 leaves on the table. Args: sigma m deg nxi nt"""
import sys
import numpy as np
from scipy.optimize import linprog
from kbopt import Model

sigma, m, d, nxi, nt = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
M = Model(sigma, m, d, nxi=nxi, nt=nt, nfine=240)
C = M.C_lib(M.xi)
A = M.design(M.xi, M.t) / np.repeat(C, len(M.t))[:, None]
rhs = np.ones(A.shape[0])
print("=== sigma=%g m=%d deg=%d  design grid %d x %d, unknowns %d" % (sigma, m, d, nxi, nt, A.shape[1]))
c = M.fit_ls(C)
print("LS            fine J_inf=%.3e J_2=%.3e" % M.metrics(c, M.C_lib))
w = np.ones(len(rhs)) / len(rhs); best = np.inf
for it in range(1, 61):
    x = M._solve(A, rhs, np.sqrt(w)); r = np.abs(A @ x - rhs)
    if r.max() < best: best, bx = r.max(), x
    w *= r; w /= w.sum()
    if it in (5, 10, 20, 40, 60):
        print("Lawson it=%2d  design max=%.3e  fine J_inf=%.3e J_2=%.3e" % ((it, best) + M.metrics(bx.reshape(M.ncell, d + 1), M.C_lib)))
nx = A.shape[1]; rows = []; h = []
for th in np.arange(32) * 2 * np.pi / 32:
    e = np.exp(-1j * th); G = np.real(e * A)
    rows.append(np.hstack([G, -np.ones((G.shape[0], 1))])); h.append(np.real(e * rhs))
cost = np.zeros(nx + 1); cost[-1] = 1
res = linprog(cost, A_ub=np.vstack(rows), b_ub=np.concatenate(h), bounds=[(None, None)] * (nx + 1), method="highs")
if res.success:
    print("LP (32-gon)   design tau=%.3e  fine J_inf=%.3e J_2=%.3e  grow=%.1f" % ((res.x[-1],) + M.metrics(res.x[:-1].reshape(M.ncell, d + 1), M.C_lib) + (M.monomial_growth(res.x[:-1].reshape(M.ncell, d + 1)),)))
else:
    print("LP failed", res.message)
