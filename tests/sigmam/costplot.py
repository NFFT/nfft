#!/usr/bin/env python3
"""Cost landscape of block 3: model cost against sigma, m stepping down.

usage: uv run --with matplotlib python costplot.py out.png [N] [type]
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import kbmodel as km

out = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 65536
typ = sys.argv[3] if len(sys.argv) > 3 else "double"
eps = km.EPS[typ]
goals = {"float": (1e-3, 1e-4, 1e-5), "double": (1e-5, 1e-9, 1e-13)}.get(typ, (1e-5, 1e-12, 1e-20))
Ms = (N // 8, N, 8 * N)
INK, INK2, GRID, S1, S2 = "#0b0b0b", "#52514e", "#ebebe6", "#2a78d6", "#eb6834"

fig, axes = plt.subplots(len(goals), len(Ms), figsize=(13, 9.5), sharex=True)
fig.patch.set_facecolor("#fcfcfb")
for i, g in enumerate(goals):
    for j, M in enumerate(Ms):
        ax = axes[i][j]
        ax.set_facecolor("#fcfcfb")
        c = km.candidates(N, M, g, eps)
        s = [n / N for n, _, _ in c]
        t = [1e3 * x for _, _, x in c]
        ax.plot(s, t, color=S1, lw=2, drawstyle="steps-post")
        low = min(x for _, _, x in c)
        near = {m for _, m, x in c if x <= 1.05 * low}
        n_c = km.choose(N, M, g, eps)[1]
        best = next(q for q, x in enumerate(c) if x[0] == n_c)
        # label each m plateau once, at its cheapest point
        seen, span = {}, {}
        for q, (n, m, x) in enumerate(c):
            if m not in seen or x < c[seen[m]][2]:
                seen[m] = q
            span[m] = span.get(m, 0) + (s[min(q + 1, len(s) - 1)] - s[q])
        for m, q in seen.items():
            if span[m] < 0.04 or q == best:
                continue
            ax.annotate("m=%d" % m, (s[q], t[q]), textcoords="offset points", xytext=(0, -13),
                        ha="left", fontsize=8, color=INK2)
        ax.scatter([s[best]], [t[best]], s=64, color=S2, zorder=3, edgecolors="#fcfcfb", linewidths=2)
        ax.annotate("best %.3f, m=%d" % (s[best], c[best][1]), (s[best], t[best]),
                    textcoords="offset points", xytext=(6, 8), fontsize=9, color=INK)
        ax.set_ylim(0, max(t) * 1.15)
        ax.set_xlim(1.2, 2.55)
        ax.grid(color=GRID, lw=1)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.set_title("goal %.0e, M = %s, %s   (within 5%%: m %s)" % (g, "N/8" if M < N else ("N" if M == N else "8N"),
                     km.regime(N, M), ",".join(map(str, sorted(near)))), fontsize=9.5, color=INK, loc="left")
        if j == 0:
            ax.set_ylabel("model ms per trafo", color=INK2, fontsize=9)
        if i == len(goals) - 1:
            ax.set_xlabel("sigma = n/N", color=INK2, fontsize=9)
fig.suptitle("Model cost of one 1-D trafo, N = %d, %s, n = fft_size (2^a 3^b 5^c 7^d, a <= 6). Orange: chosen (sigma, m)." % (N, typ),
             x=0.01, ha="left", fontsize=12, color=INK)
fig.tight_layout()
fig.savefig(out, dpi=110, facecolor=fig.get_facecolor())
