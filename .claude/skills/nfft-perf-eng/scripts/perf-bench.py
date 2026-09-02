#!/usr/bin/env python3
"""Benchmark-JSON helpers for the nfft-perf-eng skill — deterministic, stdlib-only.

Turns the codspeed walltime JSON into the skill's canonical markdown tables, and applies
the noise rule for the exit gate so the pass/fail call is reproducible, not eyeballed.

Subcommands
-----------
  snapshot <bench.json> [--prec d]
        Emit a "Benchmark snapshot" table (prec | case | median_ns | stdev_ns | rounds),
        one row per case. Used in Phase A / C / E deliverables.

  compare (--base B.json --final F.json [--prec d] | --taskdir DIR)
        Emit a "Comparison table" (prec | case | base median_ns | final median_ns | Δ% |
        threshold | verdict) and a one-line summary. With --taskdir it loops every precision,
        diffing artifacts/baseline-bench-<p>.json against artifacts/final-bench-<p>.json.
        Exit code is 1 if ANY case regressed beyond the noise rule (so a script can gate on it).

The noise rule (Phase E): a case counts as regressed only if
    final_median - base_median  >  max(3 * base_stdev, 0.02 * base_median).
A faster-or-equal case is "faster" when it clears the same band, else "within noise".

Accepts both collation shapes: a flat array of benchmark objects, or an array of process
objects each carrying a `benchmarks` list.

Run with:  uv run python perf-bench.py <subcommand> ...   (no third-party deps)
See ../details/measurement-modes.md (noise rule) and ../details/deliverables.md (formats).
"""
import argparse
import json
import os
import sys

REGRESS_PCT = 0.02   # 2% of the base median
REGRESS_SIGMA = 3.0  # 3 standard deviations
THRESHOLD_LABEL = "2%/3σ"


def iter_benchmarks(data):
    """Yield benchmark objects from either collation shape."""
    if isinstance(data, dict):
        data = [data]
    for el in data:
        if isinstance(el, dict) and "benchmarks" in el:
            for b in el["benchmarks"]:
                yield b
        else:
            yield el


def load_cases(path):
    """path -> {case_name: {'median','stdev','rounds'}} (last wins on duplicate names)."""
    with open(path) as fh:
        data = json.load(fh)
    out = {}
    for b in iter_benchmarks(data):
        name = b["name"]
        st = b.get("stats", b)
        out[name] = {
            "median": float(st["median_ns"]),
            "stdev": float(st.get("stdev_ns", 0.0)),
            "rounds": int(st.get("rounds", 0)),
        }
    return out


def fmt_int(x):
    return f"{int(round(x))}"


def cmd_snapshot(args):
    cases = load_cases(args.bench)
    prec = args.prec or "?"
    print("| prec | case | median_ns | stdev_ns | rounds |")
    print("|------|------|-----------|----------|--------|")
    for name in sorted(cases):
        c = cases[name]
        print(f"| {prec} | {name} | {fmt_int(c['median'])} | {fmt_int(c['stdev'])} | {c['rounds']} |")
    return 0


def compare_pair(base, final, prec, rows, regressed_out=None):
    """Append comparison rows for one precision; return (n_regressed, n_missing).

    If regressed_out is a list, append (prec, name) for each regressed case (used by the
    confirm re-run to re-measure only the flagged cases)."""
    n_reg = 0
    n_missing = 0
    for name in sorted(base):
        if name not in final:
            n_missing += 1
            sys.stderr.write(f"WARN: [{prec}] case missing from final: {name}\n")
            continue
        b, f = base[name], final[name]
        bm, fm = b["median"], f["median"]
        delta = fm - bm
        pct = (delta / bm * 100.0) if bm else 0.0
        threshold = max(REGRESS_SIGMA * b["stdev"], REGRESS_PCT * bm)
        if delta > threshold:
            verdict = "❌ regressed"
            n_reg += 1
            if regressed_out is not None:
                regressed_out.append((prec, name))
        elif -delta > threshold:
            verdict = "✅ faster"
        else:
            verdict = "≈ noise"
        rows.append(
            f"| {prec} | {name} | {fmt_int(bm)} | {fmt_int(fm)} | {pct:+.1f}% | {THRESHOLD_LABEL} | {verdict} |"
        )
    # cases only in final (new) are informational
    for name in sorted(final):
        if name not in base:
            sys.stderr.write(f"WARN: [{prec}] case only in final (new): {name}\n")
    return n_reg, n_missing


def cmd_compare(args):
    pairs = []  # (prec, base_path, final_path)
    if args.taskdir:
        art = os.path.join(args.taskdir, "artifacts")
        for p in ("d", "f", "l"):
            bp = os.path.join(art, f"baseline-bench-{p}.json")
            fp = os.path.join(art, f"final-bench-{p}.json")
            if os.path.exists(bp) and os.path.exists(fp):
                pairs.append((p, bp, fp))
            else:
                sys.stderr.write(f"WARN: [{p}] missing {bp} or {fp} — skipping precision\n")
        if not pairs:
            sys.stderr.write("error: no baseline/final benchmark pairs found under "
                             f"{art}\n")
            return 2
    else:
        if not (args.base and args.final):
            sys.stderr.write("error: pass --taskdir, or both --base and --final\n")
            return 2
        pairs.append((args.prec or "?", args.base, args.final))

    rows = []
    regressed = []
    total_reg = 0
    for prec, bp, fp in pairs:
        n_reg, _ = compare_pair(load_cases(bp), load_cases(fp), prec, rows, regressed)
        total_reg += n_reg

    if args.emit_regressed:
        # machine-readable: one "prec<TAB>case" line per regressed case (for perf-confirm.sh)
        for prec, name in regressed:
            print(f"{prec}\t{name}")
        return 1 if total_reg else 0

    print("| prec | case | base median_ns | final median_ns | Δ% | threshold | verdict |")
    print("|------|------|----------------|-----------------|------|-----------|---------|")
    for r in rows:
        print(r)
    print()
    if total_reg:
        print(f"**VERDICT: {total_reg} case(s) regressed beyond the noise rule "
              f"({THRESHOLD_LABEL}). Re-run them before believing it (noise rarely survives a "
              f"second run; a real regression does).**")
    else:
        print(f"**VERDICT: no case regressed beyond the noise rule ({THRESHOLD_LABEL}).**")
    return 1 if total_reg else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="nfft-perf-eng benchmark JSON helper")
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("snapshot", help="emit a Benchmark snapshot table")
    s.add_argument("bench", help="collated benchmark JSON")
    s.add_argument("--prec", help="precision label d|f|l for the prec column")
    s.set_defaults(func=cmd_snapshot)

    c = sub.add_parser("compare", help="emit a Comparison table + noise-rule verdict")
    c.add_argument("--base", help="baseline collated benchmark JSON")
    c.add_argument("--final", help="final collated benchmark JSON")
    c.add_argument("--prec", help="precision label for the single-pair form")
    c.add_argument("--taskdir", help="task dir; loops d/f/l over artifacts/{baseline,final}-bench-<p>.json")
    c.add_argument("--emit-regressed", action="store_true",
                   help="print only 'prec<TAB>case' for regressed cases (for perf-confirm.sh); no table")
    c.set_defaults(func=cmd_compare)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
