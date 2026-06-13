import itertools
import cmath

import mpmath
import pytest

from tests.refgen import transforms as T


def test_freqs_match_c_ordering():
    # NFFT even/odd, NFCT, NFST per-dim ranges
    assert T.freqs("nfft", [10]) == [list(range(-5, 5))]
    assert T.freqs("nfft", [5]) == [list(range(-2, 3))]
    assert T.freqs("nfft", [1]) == [[0]]
    assert T.freqs("nfct", [4]) == [[0, 1, 2, 3]]
    assert T.freqs("nfst", [4]) == [[1, 2, 3]]


def test_nn_counts():
    assert T.nn("nfft", [10, 20]) == 200
    assert T.nn("nfct", [10, 25]) == 250
    assert T.nn("nfst", [10, 25]) == 9 * 24


def test_ndft_equispaced_matches_dft():
    # On equispaced nodes x_j = j/N - 1/2 the NDFT reduces to a DFT we can
    # cross-check with cmath in double precision.
    mpmath.mp.dps = 30
    N, M = 8, 8
    x = [(mpmath.mpf(j) / N - mpmath.mpf(1) / 2,) for j in range(M)]
    f_hat = [mpmath.mpc(j + 1, -j) for j in range(N)]
    f = T.trafo("nfft", [N], M, x, f_hat)
    ks = list(range(-(N // 2), N - (N // 2)))
    for j in range(M):
        xj = float(x[j][0])
        ref = sum((j_c.real + 1j * j_c.imag) * cmath.exp(-2j * cmath.pi * k * xj)
                  for k, j_c in zip(ks, [complex(float(c.real), float(c.imag)) for c in f_hat]))
        got = complex(float(f[j].real), float(f[j].imag))
        assert abs(got - ref) < 1e-10


def test_adjoint_is_conjugate_transpose_nfft():
    mpmath.mp.dps = 30
    N, M = 5, 4
    x = [(mpmath.mpf(2 * j + 1) / (4 * M) - mpmath.mpf(1) / 2,) for j in range(M)]
    f = [mpmath.mpc(j + 1, j) for j in range(M)]
    f_hat = T.adjoint("nfft", [N], M, x, f)
    ks = list(range(-(N // 2), N - (N // 2)))
    for idx, k in enumerate(ks):
        ref = sum(f[j] * mpmath.expj(2 * mpmath.pi * k * x[j][0]) for j in range(M))
        assert abs(f_hat[idx] - ref) < mpmath.mpf(10) ** -25


def test_ndct_single_term_handcomputed():
    mpmath.mp.dps = 40
    # one coefficient f_hat[k=2]=3, one node x=1/8 -> f = 3*cos(2*pi*2*1/8)=3*cos(pi/2)=0
    x = [(mpmath.mpf(1) / 8,)]
    f_hat = [mpmath.mpf(0), mpmath.mpf(0), mpmath.mpf(3), mpmath.mpf(0)]
    f = T.trafo("nfct", [4], 1, x, f_hat)
    assert abs(f[0]) < mpmath.mpf(10) ** -30


def test_ndst_single_term_handcomputed():
    mpmath.mp.dps = 40
    # f_hat indexes k=1..3; coeff for k=2 is 5, node x=1/8 -> 5*sin(2*pi*2/8)=5*sin(pi/2)=5
    x = [(mpmath.mpf(1) / 8,)]
    f_hat = [mpmath.mpf(0), mpmath.mpf(5), mpmath.mpf(0)]  # k=1,2,3
    f = T.trafo("nfst", [4], 1, x, f_hat)
    assert abs(f[0] - 5) < mpmath.mpf(10) ** -30
