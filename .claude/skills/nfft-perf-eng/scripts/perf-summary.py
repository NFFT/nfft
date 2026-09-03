#!/usr/bin/env python3
"""Summary-assembly support — charts + completeness check. Deterministic, stdlib-only.

The Phase-F human write-up (`summary.html`) must present the WHOLE run to a reviewer and
link every asset. Two deterministic jobs that should not be hand-done:

  charts --taskdir DIR
        Generate self-contained inline-SVG charts from the artifacts (no JS, no network):
          * artifacts/chart-speedup-{d,f,l}.svg — per-case % faster (base→final), one chart
            per precision; green=faster / red=regressed / grey=within-noise, threshold drawn.
          * artifacts/chart-trend-{d,f,l}.svg   — log-log error-vs-N, baseline vs optimized,
            one per precision (only where artifacts/trend-{baseline,optimized}-{p}.dat exist).
        These are the REQUIRED visualizations (see ../details/deliverables.md#required-visualizations).

  check --taskdir DIR
        Verify summary.html links every deliverable + artifact in the task dir (nothing
        orphaned) and that the required charts exist and are referenced. Exit 1 on any gap.

Run with:  uv run python perf-summary.py <charts|check> --taskdir DIR
"""
import argparse
import glob
import json
import math
import os
import re
import sys

REGRESS_PCT = 0.02
REGRESS_SIGMA = 3.0


# ---------------------------------------------------------------- data loading
def load_cases(path):
    with open(path) as fh:
        data = json.load(fh)
    if isinstance(data, dict):
        data = [data]
    out = {}
    for el in data:
        items = el["benchmarks"] if isinstance(el, dict) and "benchmarks" in el else [el]
        for b in items:
            st = b.get("stats", b)
            out[b["name"]] = (float(st["median_ns"]), float(st.get("stdev_ns", 0.0)))
    return out


def compare_rows(taskdir):
    """[(prec, case, base, final, pct_faster, verdict)] over every precision present."""
    art = os.path.join(taskdir, "artifacts")
    rows = []
    for p in ("d", "f", "l"):
        bp = os.path.join(art, f"baseline-bench-{p}.json")
        fp = os.path.join(art, f"final-bench-{p}.json")
        if not (os.path.exists(bp) and os.path.exists(fp)):
            continue
        base, final = load_cases(bp), load_cases(fp)
        for name in sorted(base):
            if name not in final:
                continue
            bm, (fm, _) = base[name][0], (final[name][0], 0)
            bstd = base[name][1]
            delta = fm - bm
            thr = max(REGRESS_SIGMA * bstd, REGRESS_PCT * bm)
            pct_faster = (bm - fm) / bm * 100.0 if bm else 0.0
            if delta > thr:
                verdict = "regressed"
            elif -delta > thr:
                verdict = "faster"
            else:
                verdict = "noise"
            rows.append((p, name, bm, fm, pct_faster, verdict))
    return rows


def read_points(path):
    pts = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.replace(",", " ").split()
            if len(parts) >= 2:
                n, e = float(parts[0]), float(parts[1])
                if n > 0 and e > 0:
                    pts.append((n, e))
    return pts


def loglog_fit(pts):
    xs = [math.log10(n) for n, _ in pts]
    ys = [math.log10(e) for _, e in pts]
    k = len(xs)
    mx, my = sum(xs) / k, sum(ys) / k
    sxx = sum((x - mx) ** 2 for x in xs) or 1e-30
    p = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    b = my - p * mx
    return p, b  # slope, intercept (in log10 space)


# ---------------------------------------------------------------- tiny SVG kit
def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


COL = {"faster": "#5f8a59", "regressed": "#b4412f", "noise": "#9a8f7a"}


def speedup_svg(rows, title):
    """One precision's cases (rows = [(prec, case, base, final, pct, verdict), …])."""
    pad_l, pad_r, pad_t, pad_b = 250, 56, 46, 30
    rh = 22
    h = pad_t + pad_b + rh * max(1, len(rows))
    w = 760
    plot_w = w - pad_l - pad_r
    maxabs = max((abs(r[4]) for r in rows), default=1.0) or 1.0
    maxabs = max(maxabs, 1.0)
    zx = pad_l + plot_w * 0.42  # zero axis (room for faster bars to the right)
    sx = (w - pad_r - zx) / maxabs  # px per % on the faster side
    sxl = (zx - pad_l) / maxabs     # px per % on the slower side
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
           f'font-family="system-ui,Arial,sans-serif" font-size="12">']
    out.append(f'<text x="6" y="22" font-size="14" font-weight="600">{esc(title)}</text>')
    out.append(f'<line x1="{zx:.1f}" y1="{pad_t-6}" x2="{zx:.1f}" y2="{h-pad_b}" '
               f'stroke="#bbb" stroke-width="1"/>')
    out.append(f'<text x="{zx:.1f}" y="{pad_t-10}" text-anchor="middle" fill="#888">0%</text>')
    for i, (p, name, bm, fm, pct, verdict) in enumerate(rows):
        y = pad_t + i * rh
        out.append(f'<text x="6" y="{y+15}" fill="#333">{esc(name)}</text>')
        if pct >= 0:
            bw = pct * sx
            out.append(f'<rect x="{zx:.1f}" y="{y+4}" width="{bw:.1f}" height="{rh-8}" '
                       f'fill="{COL[verdict]}" opacity="0.9"/>')
            out.append(f'<text x="{zx+bw+5:.1f}" y="{y+15}" fill="{COL[verdict]}">{pct:+.1f}%</text>')
        else:
            bw = -pct * sxl
            out.append(f'<rect x="{zx-bw:.1f}" y="{y+4}" width="{bw:.1f}" height="{rh-8}" '
                       f'fill="{COL[verdict]}" opacity="0.9"/>')
            out.append(f'<text x="{zx-bw-5:.1f}" y="{y+15}" text-anchor="end" '
                       f'fill="{COL[verdict]}">{pct:+.1f}%</text>')
    out.append(f'<text x="6" y="{h-pad_b+18}" fill="#888">'
               f'bars → faster · green faster · red regressed · grey noise (threshold max(3σ, 2%))</text>')
    out.append('</svg>')
    return "\n".join(out)


def trend_svg(base, opt, title):
    w, h = 560, 420
    pad = 64
    allpts = base + opt
    xs = [math.log10(n) for n, _ in allpts]
    ys = [math.log10(e) for _, e in allpts]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    if x1 == x0:
        x1 += 1
    if y1 == y0:
        y1 += 1
    pb, pi = loglog_fit(base)
    po, oi = loglog_fit(opt)

    def X(lx):
        return pad + (lx - x0) / (x1 - x0) * (w - 2 * pad)

    def Y(ly):
        return h - pad - (ly - y0) / (y1 - y0) * (h - 2 * pad)

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
           f'font-family="system-ui,Arial,sans-serif" font-size="12">']
    out.append(f'<text x="{pad}" y="24" font-size="14" font-weight="600">{esc(title)}</text>')
    # axes
    out.append(f'<line x1="{pad}" y1="{h-pad}" x2="{w-pad}" y2="{h-pad}" stroke="#999"/>')
    out.append(f'<line x1="{pad}" y1="{pad}" x2="{pad}" y2="{h-pad}" stroke="#999"/>')
    for lx in range(math.floor(x0), math.ceil(x1) + 1):
        out.append(f'<text x="{X(lx):.1f}" y="{h-pad+18}" text-anchor="middle" fill="#888">'
                   f'10^{lx}</text>')
    for ly in range(math.floor(y0), math.ceil(y1) + 1):
        out.append(f'<text x="{pad-8}" y="{Y(ly)+4:.1f}" text-anchor="end" fill="#888">'
                   f'10^{ly}</text>')
    out.append(f'<text x="{w/2:.0f}" y="{h-16}" text-anchor="middle" fill="#555">N</text>')
    out.append(f'<text x="16" y="{h/2:.0f}" text-anchor="middle" fill="#555" '
               f'transform="rotate(-90 16 {h/2:.0f})">error</text>')
    for pts, col, slope, inter, lbl in (
        (base, "#26506e", pb, pi, f"baseline  p={pb:.2f}"),
        (opt, "#b4412f", po, oi, f"optimized p={po:.2f}"),
    ):
        # fit line across the x-range
        out.append(f'<line x1="{X(x0):.1f}" y1="{Y(slope*x0+inter):.1f}" '
                   f'x2="{X(x1):.1f}" y2="{Y(slope*x1+inter):.1f}" stroke="{col}" '
                   f'stroke-width="2" opacity="0.7"/>')
        for n, e in pts:
            out.append(f'<circle cx="{X(math.log10(n)):.1f}" cy="{Y(math.log10(e)):.1f}" '
                       f'r="3.5" fill="{col}"/>')
    # legend
    ly = pad + 6
    for col, lbl in (("#26506e", f"baseline  p={pb:.2f}"), ("#b4412f", f"optimized p={po:.2f}")):
        out.append(f'<rect x="{w-pad-150}" y="{ly-10}" width="12" height="12" fill="{col}"/>')
        out.append(f'<text x="{w-pad-134}" y="{ly}" fill="#333">{lbl}</text>')
        ly += 18
    out.append('</svg>')
    return "\n".join(out)


# ---------------------------------------------------------------- subcommands
def cmd_charts(args):
    art = os.path.join(args.taskdir, "artifacts")
    os.makedirs(art, exist_ok=True)
    written = []
    rows = compare_rows(args.taskdir)
    prec_name = {"d": "double", "f": "float", "l": "long double"}
    by_prec = {}
    for r in rows:
        by_prec.setdefault(r[0], []).append(r)
    # one readable chart per precision (kept separate so the full case list doesn't crowd)
    for p in ("d", "f", "l"):
        if by_prec.get(p):
            path = os.path.join(art, f"chart-speedup-{p}.svg")
            title = f"Performance · {prec_name[p]} — % faster (base→final)"
            with open(path, "w") as fh:
                fh.write(speedup_svg(by_prec[p], title))
            written.append(path)
    if not rows:
        sys.stderr.write("NOTE: no baseline/final benchmark pairs — skipping speedup charts "
                         "(expected for a Phase B/C hard-gate block)\n")
    any_trend = False
    for p in ("d", "f", "l"):
        tb = os.path.join(art, f"trend-baseline-{p}.dat")
        to = os.path.join(art, f"trend-optimized-{p}.dat")
        if os.path.exists(tb) and os.path.exists(to):
            base, opt = read_points(tb), read_points(to)
            if len(base) >= 2 and len(opt) >= 2:
                path = os.path.join(art, f"chart-trend-{p}.svg")
                with open(path, "w") as fh:
                    fh.write(trend_svg(base, opt, f"Accuracy · {prec_name[p]} — error vs N (log-log)"))
                written.append(path)
                any_trend = True
    if not any_trend:
        sys.stderr.write("NOTE: no trend-{baseline,optimized}-{d,f,l}.dat — skipping trend charts "
                         "(run the differential trend study per precision to produce them)\n")
    for p in written:
        print(f"wrote {p}")
    return 0 if written else 1


def cmd_check(args):
    td = args.taskdir
    summ = os.path.join(td, "summary.html")
    if not os.path.exists(summ):
        sys.stderr.write(f"error: {summ} not found\n")
        return 1
    html = open(summ, encoding="utf-8").read()
    refs = set(re.findall(r'(?:href|src)\s*=\s*["\']([^"\']+)["\']', html))
    refs = {r.split("#")[0].rstrip("/") for r in refs}

    # everything in the task dir that should be linked (deliverables + artifacts), except summary itself
    expected = []
    for f in sorted(os.listdir(td)):
        full = os.path.join(td, f)
        if os.path.isfile(full) and f != "summary.html" and (f.endswith(".md") or f.endswith(".html")):
            expected.append(f)
    for f in sorted(glob.glob(os.path.join(td, "artifacts", "*"))):
        if os.path.isfile(f):
            expected.append("artifacts/" + os.path.basename(f))

    ref_bases = {os.path.basename(r) for r in refs}
    missing = [e for e in expected if e not in refs and os.path.basename(e) not in ref_bases]
    charts = glob.glob(os.path.join(td, "artifacts", "chart-*.svg"))

    # broken links: local relative refs that don't resolve (e.g. template still points at a
    # chart/artifact that wasn't produced — prune the template to what exists).
    broken = []
    for r in sorted(refs):
        if not r or r.startswith(("http://", "https://", "mailto:")):
            continue
        if not os.path.exists(os.path.join(td, r)):
            broken.append(r)

    ok = True
    if missing:
        ok = False
        print("UNLINKED (present in task dir but not referenced by summary.html):")
        for m in missing:
            print(f"  - {m}")
    if broken:
        ok = False
        print("BROKEN LINKS (referenced by summary.html but absent — prune to what exists):")
        for b in broken:
            print(f"  - {b}")
    if not charts:
        print("NOTE: no charts (artifacts/chart-*.svg) — run `perf-summary.py charts` "
              "(required for a completed run; N/A for a hard-gate block).")
    else:
        for c in charts:
            rel = "artifacts/" + os.path.basename(c)
            if rel not in refs and os.path.basename(c) not in {os.path.basename(r) for r in refs}:
                ok = False
                print(f"UNLINKED CHART: {rel}")
    if ok:
        print(f"OK: summary.html links all {len(expected)} deliverable/artifact files"
              + (f" and {len(charts)} chart(s)" if charts else "") + "; no broken links")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="summary charts + completeness check")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("charts", cmd_charts), ("check", cmd_check)):
        p = sub.add_parser(name)
        p.add_argument("--taskdir", required=True)
        p.set_defaults(func=fn)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
