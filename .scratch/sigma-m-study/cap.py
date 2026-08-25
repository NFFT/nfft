#!/usr/bin/env python3
"""The attainable-accuracy cap of the Kaiser-Bessel NFFT, per sigma.

At the minimum of the error curve the truncation term and the amplified
roundoff term balance, which puts the floor at eps**(D/b).
"""

import math

EPS = {"float": 2.0**-23, "double": 2.0**-52, "long double": 2.0**-112}


def main():
    hdr = " ".join(f"{k:>12}" for k in EPS)
    print(f"{'sigma':>6} {'D':>7} {'b':>7} {'A=b-D':>7} {'gamma=D/b':>10} | {hdr}")
    for s in [1.25, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 8.0, 1e6]:
        d = 2 * math.pi * math.sqrt(1 - 1 / s)
        b = 2 * math.pi * (1 - 1 / (2 * s))
        g = d / b
        fl = [math.exp(g * math.log(e)) for e in EPS.values()]
        tag = f"{s:6.2f}" if s < 1e5 else "   1e6"
        print(
            f"{tag} {d:7.3f} {b:7.3f} {b-d:7.3f} {g:10.4f} | "
            + " ".join(f"{v:12.2e}" for v in fl)
        )
    print()
    print("b^2 - (pi/sigma)^2 = D^2, so b > D and gamma < 1 at every finite sigma.")


if __name__ == "__main__":
    main()
