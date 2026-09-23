#!/usr/bin/env python3
"""Reference values for the attainable-accuracy probe.

The forward measure

    err = max_j |f_j - s_j| / ||f_hat||_1

has, over all f_hat, the supremum max_{j,k} |E_jk| where E is the error
matrix of the transform. It is attained at f_hat = e_k, so the probe runs
one transform per k and the reference is the closed form

    s_j = exp(-2*pi*i*k*x_j).

That needs more precision than the widest hardware type, so it is produced
here with mpmath and read back by floor.c. It depends on (N, k, x) only,
never on sigma or m, so one file serves the whole plane.

Nodes are a golden-ratio sequence quantised to a multiple of 2^-p. Pass the
working type's own mantissa width (24, 53, 113) for a node the caller could
actually hold; pass 24 for every type to suppress the rounding of n*x inside
the library, which displaces the node and is otherwise the larger error.

usage: floorgen.py <out-file> <mantissa-bits> [M]
"""

import os
import sys

from mpmath import mp, mpf, cos, sin, pi

DPS = 60
FULL_K_MAX = int(os.environ.get("FULL_K_MAX", 256))  # below this every k is probed, above it a subset

NS = [int(v) for v in os.environ.get("NS", "32,128,256,512,1024,4096").split(",")]


def nodes(M, p):
    """Golden-ratio nodes in [-1/2, 1/2), each an exact multiple of 2^-p."""
    q = 1 << p
    phi = (mp.sqrt(5) - 1) / 2
    out = []
    for j in range(M):
        t = mp.frac(mpf(j) * phi)
        num = int(mp.floor((t - mpf("0.5")) * q + mpf("0.5")))
        out.append(max(-q // 2, min(q // 2 - 1, num)))
    return out


def exact_decimal(num, p):
    """num / 2^p written out exactly: num * 5^p, scaled by 10^-p."""
    d = str(abs(num) * 5**p).rjust(p, "0")
    return "%s%s.%s" % ("-" if num < 0 else "", d[:-p] or "0", d[-p:])


def probe_k(N):
    h = N // 2
    if N <= FULL_K_MAX:
        return list(range(-h, h))
    ks = set()
    for base in (-h, -h // 2, 0, h // 2):
        for d in range(-2, 3):
            ks.add(base + d)
    for d in range(1, 5):
        ks.add(h - d)
    return sorted(k for k in ks if -h <= k < h)


def main():
    out_path = sys.argv[1]
    p = int(sys.argv[2])
    M = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    mp.dps = max(DPS, p // 3 + 20)

    q = 1 << p
    xs = nodes(M, p)
    two_pi = 2 * pi

    with open(out_path, "w") as f:
        f.write("M %d\n" % M)
        # A dyadic rational has an exact finite decimal expansion, so strtold
        # recovers the node bit for bit.
        f.write("X " + " ".join(exact_decimal(v, p) for v in xs) + "\n")
        for N in NS:
            ks = probe_k(N)
            f.write("N %d %d\n" % (N, len(ks)))
            f.write("K " + " ".join(str(k) for k in ks) + "\n")
            for k in ks:
                for num in xs:
                    # k*num is an exact integer; reduce mod q before the
                    # transcendental so no significance goes to the phase.
                    a = two_pi * mpf((k * num) % q) / q
                    f.write("%s %s\n" % (mp.nstr(cos(a), 40), mp.nstr(-sin(a), 40)))
            print("N=%d  k=%d" % (N, len(ks)), file=sys.stderr)


if __name__ == "__main__":
    main()
