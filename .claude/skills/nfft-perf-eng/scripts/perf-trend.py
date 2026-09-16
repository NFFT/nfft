#!/usr/bin/env python3
"""Differential trend analysis — log-log least-squares error growth fit. Stdlib-only.

The strongest settle for an accuracy-for-speed risk (and the empirical half of the
Phase-D rounding-error analysis): fit the suite's error metric against a swept parameter
as  error ≈ C · N^p,  and compare the exponent p of the optimized code against the
baseline. A Δslope near zero ⇒ the order of growth is unchanged (risk retired); a clearly
steeper optimized slope ⇒ an order-of-growth regression (proven), even if every single
size is still under the flat bound. See ../details/extending-tests.md (differential trend
analysis) and ../details/phase-d-error-analysis.md.

Subcommands
-----------
  fit <data>
        Fit one (N, error) table; print exponent p, constant C, and R².

  compare <baseline.dat> <optimized.dat> [--tol 0.15]
        Fit both; print each exponent and Δp = p_opt - p_base with a verdict:
        |Δp| <= tol  -> retired (order of growth preserved)
        Δp  >  tol   -> PROVEN order-of-growth regression (hard no)
        Exit code 1 on a proven regression, else 0.

Each <data> file is whitespace- or comma-separated with two columns per line: `N error`.
Blank lines and `#` comments are ignored. Gather the data by running the suite at a
geometric sweep of N against a per-point oracle, measuring both the baseline and the
optimized code (`git stash` between runs) — this script only does the deterministic fit.

The oracle MUST be strictly more precise than the build under test: use refgen
(`tests.refgen.transforms`) at high mpmath `dps`, NOT a same-precision C reference. A C
`long double` oracle is fine for the double/float builds but cannot validate the
long-double build itself (quad cannot out-resolve quad). The verdict from `compare` is
`Δp = p_opt - p_base` (order-of-growth change), NOT the absolute `p`: the absolute exponent
often sits a little above the textbook √N (≈0.6-0.7) because of a sub-dominant
working-precision term (e.g. phase formation ~N·u), and that appears in the baseline too.

Run with:  uv run python perf-trend.py ...   (no third-party deps)
"""
import argparse
import math
import sys


def read_points(path):
    pts = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.replace(",", " ").split()
            if len(parts) < 2:
                continue
            n, err = float(parts[0]), float(parts[1])
            if n <= 0 or err <= 0:
                sys.stderr.write(f"WARN: dropping non-positive point ({n}, {err}) in {path}\n")
                continue
            pts.append((n, err))
    if len(pts) < 2:
        raise SystemExit(f"error: need >=2 positive (N, error) points in {path}")
    return pts


def loglog_fit(pts):
    """Least-squares fit ln(err) = ln(C) + p*ln(N). Returns (p, C, r2)."""
    xs = [math.log(n) for n, _ in pts]
    ys = [math.log(e) for _, e in pts]
    k = len(xs)
    mx = sum(xs) / k
    my = sum(ys) / k
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx == 0:
        raise SystemExit("error: all N values are identical — cannot fit a slope")
    p = sxy / sxx
    ln_c = my - p * mx
    # R^2
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (ln_c + p * x)) ** 2 for x, y in zip(xs, ys))
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot else 1.0
    return p, math.exp(ln_c), r2


def cmd_fit(args):
    p, c, r2 = loglog_fit(read_points(args.data))
    print(f"error ~ C * N^p   (log-log least squares, {args.data})")
    print(f"  exponent p = {p:.4f}")
    print(f"  constant C = {c:.4e}")
    print(f"  R^2        = {r2:.4f}")
    return 0


def cmd_compare(args):
    pb, cb, rb = loglog_fit(read_points(args.baseline))
    po, co, ro = loglog_fit(read_points(args.optimized))
    dp = po - pb
    print(f"baseline   : p = {pb:.4f}  (C = {cb:.4e}, R^2 = {rb:.4f})  [{args.baseline}]")
    print(f"optimized  : p = {po:.4f}  (C = {co:.4e}, R^2 = {ro:.4f})  [{args.optimized}]")
    print(f"Δslope     : {dp:+.4f}   (tolerance ±{args.tol})")
    if dp > args.tol:
        print("VERDICT: PROVEN order-of-growth regression — optimized error grows faster. "
              "Hard no: fix or revert (see ../details/risk-assessment.md).")
        return 1
    print("VERDICT: retired — order of growth preserved within tolerance.")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="nfft-perf-eng error-trend fit")
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fit", help="fit one (N, error) table")
    f.add_argument("data")
    f.set_defaults(func=cmd_fit)

    c = sub.add_parser("compare", help="compare baseline vs optimized exponents")
    c.add_argument("baseline")
    c.add_argument("optimized")
    c.add_argument("--tol", type=float, default=0.15, help="|Δslope| tolerance (default 0.15)")
    c.set_defaults(func=cmd_compare)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
