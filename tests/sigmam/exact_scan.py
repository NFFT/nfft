"""Worst-case exact-arithmetic error of the width-(m+1) KB NFFT, max over theta and u.

Quad precision via numpy longdouble (binary128 here); phihat via mpmath.
"""
import numpy as np
from mpmath import mp, besseli, sqrt as msqrt, mpf, pi as mpi
mp.dps = 40
L = np.longdouble
PI = L("3.14159265358979323846264338327950288")

def win(b, w, x):
    a = w*w - x*x
    out = np.empty_like(x)
    i, o = a > 0, a < 0
    r = np.sqrt(a[i]); out[i] = np.sinh(b*r)/(PI*r)
    r = np.sqrt(-a[o]); out[o] = np.sin(b*r)/(PI*r)
    out[a == 0] = b/PI
    return out

def worst(s, m, NT=200, U=48):
    s = L(s); b = PI*(2 - 1/s); w = L(m + 1)
    u = np.arange(U + 1, dtype=L)/U
    x = u[:, None] + m - np.arange(2*m + 2, dtype=L)[None, :]
    phi = win(b, w, x)
    best, at = L(0), 0.0
    for j in range(NT + 1):
        th = PI/s*j/NT
        ph = L(mp.nstr(besseli(0, (m + 1)*msqrt(mpf(str(b))**2 - mpf(str(th))**2)), 36))
        re = (phi*np.cos(th*x)).sum(1); im = (phi*np.sin(th*x)).sum(1)
        e = np.sqrt((1 - re/ph)**2 + (im/ph)**2).max()
        if e > best: best, at = e, j/NT
    return best, at

if __name__ == "__main__":
    print("sigma,m,trunc,theta_frac,kappa")
    for s in [1.25, 1.3, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0, 2.125, 2.25, 2.375, 2.5]:
        b = mpi*(2 - 1/mpf(s)); D = 2*mpi*msqrt(1 - 1/mpf(s))
        for m in range(1, 25):
            e, at = worst(s, m)
            kap = besseli(0, (m + 1)*b)/besseli(0, (m + 1)*D)
            print("%g,%d,%.6e,%.3f,%s" % (s, m, float(e), at, mp.nstr(kap, 8)), flush=True)
