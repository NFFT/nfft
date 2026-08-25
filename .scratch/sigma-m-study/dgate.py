"""Head-to-head for the dyadic tuner, simulated on the gsweep measurements.

The dyadic tuner offers n = 2^j * next_power_of_2(N). With
t = next_power_of_2(N)/N in (1, 2], those grids have sigma = t, 2t, 4t. The
gsweep grid holds several sigma pairs in ratio 2, so for each such pair the
choice between j=0 and j=1 can be replayed against real measured error --
no new data, and the legacy grid is the j=1 member of the pair.

For each (precision, t, N, M, goal, direction):

  LEGACY  sigma = 2t, with m the smallest cut-off whose *measured* error
          meets the goal -- the oracle the legacy API does not provide.
  DYADIC  the model picks a cut-off at sigma = t and at sigma = 2t, drops a
          candidate whose predicted error cannot meet the goal, and takes
          the cheaper pair under n*log2(n) + (4/5)*M*(2m+2).

The dyadic choice is then checked against the measured error at the pair it
returned, so a miss is a real miss.

All three rungs are replayed: the sweep carries sigma = t, 2t and 4t for
every t below, so a triple exists at each geometry.

Usage: python3 dgate.py <gsweep-*.csv>...
"""
import math
import statistics as st
import sys
from collections import defaultdict

import gfit
import dfit

TS = (1.25, 1.5, 1.75)
RUNGS = (0, 1, 2)
NODE_WEIGHT = 0.8
M_MAX = 32
GOALS = [10.0 ** -e for e in range(2, 30)]


def cost(n, m, M):
    return n * math.log2(n) + NODE_WEIGHT * M * (2 * m + 2)


def model_m(f, sigma, eps, N, M, goal, m_cap):
    """Smallest cut-off the model rates sufficient, or None."""
    for m in range(1, m_cap + 1):
        if gfit.model(f, sigma, m, eps, N, M) <= goal:
            return m
    return None


def measured_m(series, goal):
    """Smallest cut-off whose measured error meets the goal, or None."""
    return next((m for m, e in series if e <= goal), None)


def shape_of(N, M):
    r = M / N
    if r <= 0.5:
        return "N-dominated"
    if r <= 1.5:
        return "balanced"
    return "M-dominated"


def run(pts, fits, label):
    by_key = defaultdict(list)
    for r in pts:
        by_key[(r["prec"], r["sigma"], r["N"], r["M"])].append(r)

    rows = []
    misses = 0
    for t in TS:
        sig = {j: t * (2.0 ** j) for j in RUNGS}
        for (prec, sigma, N, M), v in by_key.items():
            if abs(sigma - sig[1]) > 1e-9:
                continue
            series = {}
            for j in RUNGS:
                sv = by_key.get((prec, sig[j], N, M))
                if sv:
                    series[j] = sv
            if 1 not in series:
                continue
            eps = gfit.EPS[prec]
            for key, direction in (("fwd", "forward"), ("adj", "adjoint")):
                at = {}
                for j, sv in series.items():
                    d = sorted((r["m"], r[key]) for r in sv if r[key] > 0.0)
                    if d:
                        at[j] = d
                if 1 not in at:
                    continue
                for g in GOALS:
                    m_leg = measured_m(at[1], g)
                    if m_leg is None:
                        continue
                    cands = []
                    for j in RUNGS:
                        if j not in at:
                            continue
                        s = sig[j]
                        seen = dict(at[j])
                        cap = min(M_MAX, int(s * N) // 2 - 1, max(seen))
                        if cap < 1:
                            continue
                        m = model_m(fits[j][key], s, eps, N, M, g, cap)
                        if m is None:
                            continue
                        cands.append((cost(int(s * N), m, M), j, s, m, seen))
                    if not cands:
                        continue
                    cands.sort()
                    _, j, s, m, seen = cands[0]
                    if seen.get(m, float("inf")) > g:
                        misses += 1
                    rows.append(dict(prec=prec, N=N, M=M, t=t, goal=g,
                                     direction=direction, j=j, m=m,
                                     shape=shape_of(N, M),
                                     speedup=cost(int(sig[1] * N), m_leg, M)
                                     / cost(int(s * N), m, M)))

    print("--- %s : %d cases, %d accuracy misses" % (label, len(rows), misses))
    for shape in ("N-dominated", "balanced", "M-dominated"):
        r = [x for x in rows if x["shape"] == shape]
        if not r:
            continue
        sp = sorted(x["speedup"] for x in r)
        picks = defaultdict(int)
        for x in r:
            picks["j=%d" % x["j"]] += 1
        picks = {k: picks[k] for k in sorted(picks)}
        print("  %-12s n=%-5d median %.3f  p10 %.3f  min %.3f  "
              "ties-or-better %3.0f %%  picks %s"
              % (shape, len(sp), st.median(sp), sp[int(0.1 * len(sp))],
                 sp[0], 100.0 * sum(1 for x in sp if x >= 0.999) / len(sp),
                 dict(picks)))
    sp = sorted(x["speedup"] for x in rows)
    print("  %-12s n=%-5d median %.3f  p10 %.3f  min %.3f"
          % ("ALL", len(sp), st.median(sp), sp[int(0.1 * len(sp))], sp[0]))
    return rows


def main(paths):
    pts = gfit.load(paths)
    bands = {j: [r for r in pts if lo <= r["sigma"] <= hi]
             for j, (_, lo, hi) in enumerate(dfit.BANDS)}
    print("points: %d  (%s)"
          % (len(pts), ", ".join("band%d %d" % (j, len(b))
                                 for j, b in sorted(bands.items()))))
    print()

    shipped = {j: dfit.SHIPPED for j in RUNGS}
    perband = {j: {k: gfit.fit_direction(bands[j], k) for k in ("fwd", "adj")}
               for j in RUNGS}

    run(pts, shipped, "shipped constants")
    print()
    run(pts, perband, "per-band refit")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
