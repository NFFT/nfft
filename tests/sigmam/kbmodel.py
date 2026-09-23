#!/usr/bin/env python3
"""Worst-case error model of the 1-D Kaiser-Bessel NFFT and the sigma/m choice.

Measure: max_j |f_j - s_j| / ||f_hat||_1, sup over f_hat. Window half-width
w = m + 1 (stencil reach). For oversampling sigma and cut-off m

    b = pi (2 - 1/sigma),  D = 2 pi sqrt(1 - 1/sigma)

    E(m) = CT R / I0(w D)  +  eps (C1 I0(w b) / I0(w D) + C0)

    R    = max(|1 - beta exp(i pi w / sigma)|, G beta),   beta = b/pi = 2 - 1/sigma

The first term is truncation at the band edge, theta = pi/sigma, node on a grid
point. Two parts of size O(1) are left there after dividing by phihat =
I0(w D): the alias image, which lands exactly on the spectrum's edge where
phihat(b) = I0(0) = 1, and the dropped sample on the support edge, phi(w) =
b/pi, with phase exp(i pi w/sigma) against it. Their interference is the
ripple, period 2 sigma in m. Just inside the band edge the alias image is
gone, and the sample plus the rest of the sin tail, G beta, is the worst case. The second term is
round-off: the part the deconvolution amplifies by 1/phihat(band edge) =
I0(w b)/I0(w D), and the part it does not. E falls with
m, reaches a minimum, and rises again. Both building blocks use this one
model, so a sigma the first returns always has an m for the second.

Constants live in kbmodel.json; kbfit.py and kbcostfit.py refit them.
`python3 kbmodel.py --c` prints them as the #define block of
kernel/nfft/sigma_m.c.
"""

import json
import math
import os

EPS = {
    "float": 2.0**-23,
    "double": 2.0**-52,
    "x87": 2.0**-63,  # long double on x86
    "quad": 2.0**-112,  # long double on aarch64
}

M_MAX = 64

_here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(_here, "kbmodel.json")) as _f:
    K = json.load(_f)


def log_i0(x):
    """log I0(x), x >= 0, relative error below 1e-12."""
    if x < 30.0:
        t, s, k = 1.0, 1.0, 0
        q = 0.25 * x * x
        while t > 1e-17 * s:
            k += 1
            t *= q / (k * k)
            s += t
        return math.log(s)
    # asymptotic series, terms ((2k-1)!!)^2 / (k! 8^k x^k)
    t, s = 1.0, 1.0
    for k in range(1, 12):
        t *= (2 * k - 1) ** 2 / (8.0 * k * x)
        s += t
    return x - 0.5 * math.log(2.0 * math.pi * x) + math.log(s)


def ripple(sigma, w, k=None):
    """R: alias image against the support-edge sample, at the band edge."""
    beta = 2.0 - 1.0 / sigma
    c = 1.0 + beta * beta - 2.0 * beta * math.cos(math.pi * w / sigma)
    return max(math.sqrt(c), (k or K)["G"] * beta)


def error(m, sigma, eps, k=None):
    """Model worst-case error at integer cut-off m."""
    k = k or K
    b = math.pi * (2.0 - 1.0 / sigma)
    D = 2.0 * math.pi * math.sqrt(1.0 - 1.0 / sigma)
    w = m + 1
    ld = log_i0(w * D)
    return k["CT"] * ripple(sigma, w, k) * math.exp(-ld) + eps * (k["C1"] * math.exp(log_i0(w * b) - ld) + k["C0"])


def best_m(sigma, eps, k=None):
    """(m, E) at the model minimum over integer m."""
    prev = (1, error(1, sigma, eps, k))
    for m in range(2, M_MAX + 1):
        e = error(m, sigma, eps, k)
        if e >= prev[1]:
            break
        prev = (m, e)
    return prev


def floor_error(sigma, eps, k=None):
    """Best worst-case error reachable at sigma for any integer m."""
    return best_m(sigma, eps, k)[1]


def smallest_sigma(goal, eps, sigma_min=1.25, sigma_max=2.5, k=None):
    """Building block 1: smallest sigma in [sigma_min, sigma_max] whose floor
    meets goal. A goal below the floor at sigma_max is raised to it, so the
    result is always valid, and is sigma_max when nothing smaller suffices."""
    goal = max(goal, floor_error(sigma_max, eps, k))
    if floor_error(sigma_min, eps, k) <= goal:
        return sigma_min
    lo, hi = sigma_min, sigma_max  # floor(lo) > goal >= floor(hi)
    for _ in range(64):  # as the C port: fixed, any precision
        mid = 0.5 * (lo + hi)
        if floor_error(mid, eps, k) <= goal:
            hi = mid
        else:
            lo = mid
    return hi


def smallest_m(sigma, goal, eps, k=None):
    """Building block 2: smallest integer m >= 1 whose model error at sigma
    meets goal. The caller guarantees sigma can reach goal; if it cannot, the
    m at the model minimum is returned."""
    mb, _ = best_m(sigma, eps, k)
    for m in range(1, mb + 1):
        if error(m, sigma, eps, k) <= goal:
            return m
    return mb


def cost(N, M, n, m, k=None):
    """Model seconds per trafo: FFT, per node, per tap. Constants: kbcostfit.py."""
    c = (k or K)["cost"]
    return c["fft"] * n * math.log2(n) + c["M"] * M + c["tap"] * M * (2 * m + 2)


FFT_AMAX = 6  # 5 for complex long double: 2^a * sizeof(complex) <= 1 KB


def fft_sizes(lo, hi, amax=FFT_AMAX):
    """Sorted n in [lo, hi] of the form 2^a 3^b 5^c 7^d with 1 <= a <= amax.

    Only FFTW's specialised radices, and a power-of-two part of at most 1 KB
    in bytes: a quarter of the smallest common L1 set period (4 KB), so the
    FFT's strides cannot line up in a few cache sets on any usual host."""
    out = []
    p2 = 2
    for _ in range(amax):
        p3 = p2
        while p3 <= hi:
            p5 = p3
            while p5 <= hi:
                p7 = p5
                while p7 <= hi:
                    if p7 >= lo:
                        out.append(p7)
                    p7 *= 7
                p5 *= 5
            p3 *= 3
        p2 *= 2
    return sorted(out)


def fft_amax(eps):
    """a limit for the type: 5 for complex long double (32 bytes), else 6."""
    return 5 if eps < 1e-18 else FFT_AMAX


def fft_size(n, amax=FFT_AMAX):
    """n*: the smallest efficient FFT size >= n. At most 5 % above n for
    n >= 1024, median 0.8 %."""
    return fft_sizes(n, 2 * max(n, 2), amax)[0]


def candidates(N, M, goal, eps, sigma_min=1.25, sigma_max=2.5, k=None):
    """(n, m, cost) for every efficient FFT size n from the smallest usable
    sigma to sigma_max, each with the smallest m meeting the goal. Always at
    least one: fft_size of the lower end, even where it passes sigma_max."""
    goal = max(goal, floor_error(sigma_max, eps, k))
    lo = smallest_sigma(goal, eps, sigma_min, sigma_max, k)
    a = fft_amax(eps)
    ns = fft_sizes(math.ceil(lo * N), math.floor(sigma_max * N), a) or [fft_size(math.ceil(lo * N), a)]
    out = []
    for n in ns:
        m = smallest_m(n / N, goal, eps, k)
        out.append((n, m, cost(N, M, n, m, k)))
    return out


def regime(N, M, k=None):
    """"fft" (few nodes), "balanced", or "conv" (many nodes), by M/N.

    The bounds are where the cost model stops having one clear winner: below
    the band, near-ties all sit at the smallest sigma; above it there is one
    pair within the tie tolerance; inside, tied pairs lie far apart in sigma."""
    lo, hi = (k or K)["choice"]["balanced"]
    return "fft" if M < lo * N else ("conv" if M > hi * N else "balanced")


def sigma_for_m(m, goal, eps, lo, hi, k=None):
    """Smallest sigma in [lo, hi] at which cut-off m meets goal, or None.
    At fixed m the model error falls monotonically in sigma, so bisection."""
    if error(m, hi, eps, k) > goal:
        return None
    if error(m, lo, eps, k) <= goal:
        return lo
    for _ in range(64):  # as the C port: fixed, any precision
        mid = 0.5 * (lo + hi)
        if error(m, mid, eps, k) <= goal:
            hi = mid
        else:
            lo = mid
    return hi


def choose(N, M, goal, eps, sigma_min=1.25, sigma_max=2.5, k=None):
    """Building block 3: (sigma, n, m) for an N-coefficient, M-node 1-D problem
    that meets goal, floored at what sigma_max reaches. No search over n.

    Cost is a sawtooth in sigma: rising with n on each m plateau, falling
    where m steps down. The choice is therefore always the start of a plateau,
    and there is one per m.
      fft       the FFT dominates: the smallest usable sigma.
      conv      the taps dominate: the smallest m sigma_max allows, at the
                smallest sigma where it suffices.
      balanced  one candidate per m between those two; the smallest sigma
                within the tie tolerance of the cheapest."""
    goal = max(goal, floor_error(sigma_max, eps, k))
    lo = smallest_sigma(goal, eps, sigma_min, sigma_max, k)
    a = fft_amax(eps)
    n_lo = fft_size(math.ceil(lo * N), a)
    top = fft_sizes(n_lo, math.floor(sigma_max * N), a)
    n_hi = top[-1] if top else n_lo

    def at(m):
        """The plateau of m: its first efficient size, and the m there."""
        s = sigma_for_m(m, goal, eps, lo, n_hi / N, k)
        if s is None:
            return None
        n = min(fft_size(math.ceil(s * N), a), n_hi)
        mm = smallest_m(n / N, goal, eps, k)
        return n, mm, cost(N, M, n, mm, k)

    r = regime(N, M, k)
    if r == "fft":
        m = smallest_m(n_lo / N, goal, eps, k)
        return n_lo / N, n_lo, m
    m_min = smallest_m(n_hi / N, goal, eps, k)
    if r == "conv":
        n, m, _ = at(m_min)
        return n / N, n, m
    m_max = smallest_m(n_lo / N, goal, eps, k)
    c = [x for x in (at(m) for m in range(m_max, m_min - 1, -1)) if x]
    tie = (k or K)["choice"]["tie"]
    low = min(x[2] for x in c)
    n, m, _ = min((x for x in c if x[2] <= (1 + tie) * low), key=lambda x: x[0])
    return n / N, n, m


def choose_scan(N, M, goal, eps, sigma_min=1.25, sigma_max=2.5, k=None):
    """choose() by exhaustive scan over every efficient size; reference only."""
    c = candidates(N, M, goal, eps, sigma_min, sigma_max, k)
    low = min(x[2] for x in c)
    if regime(N, M, k) == "balanced":
        tie = (k or K)["choice"]["tie"]
        n, m, _ = min((x for x in c if x[2] <= (1 + tie) * low), key=lambda x: x[0])
    else:
        n, m, _ = min(c, key=lambda x: x[2])
    return n / N, n, m


def _check():
    # fft_size against brute force
    def ok(n, amax):
        a = (n & -n).bit_length() - 1
        r = n >> a
        for p in (3, 5, 7):
            while r % p == 0:
                r //= p
        return r == 1 and 1 <= a <= amax
    for amax in (5, 6):
        for n in range(1, 5000):
            f = fft_size(n, amax)
            assert f >= n and ok(f, amax) and not any(ok(q, amax) for q in range(n, f)), (n, f)
    gaps = [(fft_size(n + 1) - n) / n for n in fft_sizes(1024, 1 << 21)]
    assert max(gaps) <= 0.05 + 1e-12, max(gaps)
    # log_i0 against the power series at the switch-over and far out
    for x in (29.9, 30.0, 45.0):
        n = 400
        t, s = 1.0, 1.0
        for j in range(1, n):
            t *= 0.25 * x * x / (j * j)
            s += t
        assert abs(log_i0(x) - math.log(s)) < 1e-11, x
    for name, eps in EPS.items():
        grid = [1.25 + 1.25 * i / 200 for i in range(201)]
        fl = [floor_error(s, eps) for s in grid]
        assert all(a >= c for a, c in zip(fl, fl[1:])), name
        assert smallest_sigma(1.0, eps) == 1.25
        assert abs(smallest_sigma(0.0, eps) - 2.5) < 1e-12
        for g in (fl[0] * 2, fl[60], fl[150], fl[-1] * 1.0001):
            s = smallest_sigma(g, eps)
            assert floor_error(s, eps) <= g * (1 + 1e-12), (name, g, s)
            m = smallest_m(s, g, eps)
            assert error(m, s, eps) <= g * (1 + 1e-12), (name, g, s, m)
            assert m == 1 or error(m - 1, s, eps) > g, (name, g, s, m)
        for s in (1.25, 2.0, 2.5):
            mb, e = best_m(s, eps)
            assert smallest_m(s, 0.0, eps) == mb
            assert smallest_m(s, 1.0, eps) == 1
            assert error(mb + 1, s, eps) >= e and (mb == 1 or error(mb - 1, s, eps) > e)
        for N, M, g in ((1024, 1024, 1e-9), (4096, 40960, 1e-12), (65536, 1000, 1e-6), (3, 5, 1e-3),
                        (64, 8, 1e-18), (262144, 262144, 1e-12), (1000, 64000, 1e-7)):
            s, n, m = choose(N, M, g, eps)
            gg = max(g, floor_error(2.5, eps))
            assert error(m, s, eps) <= gg * (1 + 1e-12), (name, N, M, g)
            assert n == fft_size(n, fft_amax(eps)) and s >= smallest_sigma(g, eps), (name, N, M, g)
            # greedy against the exhaustive scan: the same pair where the
            # tie rule or a single winner decides; in the fft regime never a
            # larger sigma, and within the tie tolerance of the cheapest
            low = min(c[2] for c in candidates(N, M, g, eps))
            ref = choose_scan(N, M, g, eps)
            assert cost(N, M, n, m) <= (1 + K["choice"]["tie"]) * low * (1 + 1e-12)
            if regime(N, M) == "fft":
                assert n <= ref[1], (name, N, M, g)
            else:
                assert (n, m) == ref[1:], (name, N, M, g)


def c_defines(k=None):
    """The constant block of kernel/nfft/sigma_m.c."""
    k = k or K
    c, ch = k["cost"], k["choice"]
    return "\n".join([
        "#define CT K(%g)" % k["CT"],
        "#define C1 K(%g)" % k["C1"],
        "#define C0 K(%g)" % k["C0"],
        "#define G K(%g)" % k["G"],
        "",
        "#define COST_FFT K(%g)" % c["fft"],
        "#define COST_NODE K(%g)" % c["M"],
        "#define COST_TAP K(%g)" % c["tap"],
        "",
        "/* Regimes by M/N, and the tie tolerance of the balanced one. */",
        "#define BALANCED_LO K(%g)" % ch["balanced"][0],
        "#define BALANCED_HI K(%.1f)" % ch["balanced"][1],
        "#define TIE K(%g)" % ch["tie"],
    ])


if __name__ == "__main__":
    import sys
    if "--c" in sys.argv:
        print(c_defines())
        sys.exit(0)
    _check()
    goals = [1e-3, 1e-5, 1e-7, 1e-9, 1e-11, 1e-13, 1e-15, 1e-18, 1e-25, 1e-30]
    print("constants", K)
    print("\n%-7s %9s %9s   sigma, m for goal" % ("type", "fl@1.25", "fl@2.5"))
    print("%27s " % "" + " ".join("%11.0e" % g for g in goals))
    for name, eps in EPS.items():
        cells = []
        for g in goals:
            s = smallest_sigma(g, eps)
            cells.append("%6.3f,%-3d" % (s, smallest_m(s, max(g, floor_error(2.5, eps)), eps)))
        print("%-7s %9.2e %9.2e   " % (name, floor_error(1.25, eps), floor_error(2.5, eps))
              + " ".join("%11s" % c for c in cells))
    print("\nblock 3, double: sigma, m chosen")
    for N in (1024, 65536):
        for M in (N // 8, N, 8 * N, 16 * N):
            print("  N %6d M %7d %-8s  " % (N, M, regime(N, M)) + "  ".join(
                "%.0e: %.3f,%d" % ((g,) + choose(N, M, g, EPS["double"])[::2]) for g in (1e-5, 1e-9, 1e-13)))
