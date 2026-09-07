#!/usr/bin/env python3
"""Parse `checkall` output into the correctness net — deterministic, stdlib-only.

The mechanical half of Phases B and E: the agent injects the (target-specific) fault or
makes the optimization; this script turns the resulting `checkall` stdout into the net —
which cases flipped to `-> FAIL`/`-> ERROR`, in the skill's canonical table shape, and the
failing-case name set for diffing one run against another. No eyeballing of grep output.

The harness prints, per case (tests/nfft.c etc.):
    <filename>  d = .., N = [..], M = .., <init>, m = .., <trafo>  -> FAIL <err> (<bound>)
and simpler lines for util (`log2i(..) = .. -> FAIL`). This parser keys on the trailing
`-> (OK|FAIL|ERROR)` and the optional `<err> (<bound>)`.

Subcommands
-----------
  table <log|->            Correctness-net markdown rows for every FAIL/ERROR line
  names <log|->            one failing case id per line (sorted, unique — for diffing nets)
  check <log|-> [--expect FILE]
                           exit 0 iff GREEN (no FAIL/ERROR); with --expect, exit 0 iff the
                           failing-name set equals FILE's lines (the net is unchanged)

`<log>` is a file path or `-` for stdin. Run with:  uv run python perf-net.py ...
See ../details/phase-b-correctness-net.md and ../details/phase-e-inner-loop.md.
"""
import argparse
import re
import sys

# trailing "-> STATUS [err (bound)]"; err/bound optional (util lines carry none)
LINE_RE = re.compile(
    r"->\s+(?P<status>OK|FAIL|ERROR)\b"
    r"(?:\s+(?P<err>[-+0-9.eEnafiNA]+)\s*\(\s*(?P<bound>[-+0-9.eEnafiNA]+)\s*\))?"
)
# the other harness format prints the numbers BEFORE the arrow: "err_rel = <err> < <bound>"
ERRREL_RE = re.compile(
    r"err_rel\s*=\s*(?P<err>[-+0-9.eEnafiNA]+)\s*<\s*(?P<bound>[-+0-9.eEnafiNA]+)"
)


def case_id(label):
    """Compact, stable id: '<filename> … <trafo delegate>' (matches the docs' net column)."""
    label = " ".join(label.split())
    first = label.split()[0] if label.split() else label
    delegate = ""
    if "," in label:
        delegate = label.rsplit(",", 1)[-1].strip()
    return f"{first} … {delegate}" if delegate and delegate != first else first


def suite_guess(case):
    """Coarse suite from the filename prefix (nfft_… → nfft); '?' if unknown."""
    head = case.split()[0]
    for s in ("nfft", "nfct", "nfst", "nnfft", "nfsft", "nfsoft", "nsfft", "fpt", "solver", "util"):
        if head.startswith(s):
            return s
    return "?"


def parse(stream):
    """Yield (status, case_id, err, bound) for every classified line."""
    for raw in stream:
        m = LINE_RE.search(raw)
        if not m:
            continue
        label = raw[: m.start()].rstrip()
        err, bound = m.group("err"), m.group("bound")
        if err is None:  # numbers may precede the arrow (err_rel = <err> < <bound>)
            em = ERRREL_RE.search(label)
            if em:
                err, bound = em.group("err"), em.group("bound")
        yield m.group("status"), case_id(label), err or "—", bound or "—"


def failing(path):
    if path == "-":
        return [(c, e, b) for st, c, e, b in parse(sys.stdin) if st in ("FAIL", "ERROR")]
    with open(path) as f:
        return [(c, e, b) for st, c, e, b in parse(f) if st in ("FAIL", "ERROR")]


def cmd_table(args):
    fails = failing(args.log)
    print("| suite | case | error | bound |")
    print("|-------|------|-------|-------|")
    seen = set()
    for c, e, b in fails:
        if c in seen:
            continue
        seen.add(c)
        print(f"| {suite_guess(c)} | {c} | {e} | {b} |")
    print(f"\n**Net size:** {len(seen)} cases", file=sys.stderr)
    return 0


def cmd_names(args):
    names = sorted({c for c, _, _ in failing(args.log)})
    for n in names:
        print(n)
    return 0


def cmd_check(args):
    names = {c for c, _, _ in failing(args.log)}
    if args.expect:
        expected = {ln.strip() for ln in open(args.expect) if ln.strip()}
        if names == expected:
            print(f"OK: failing set matches expected ({len(names)} cases)")
            return 0
        only_now = names - expected
        only_exp = expected - names
        if only_now:
            print(f"UNEXPECTED failures ({len(only_now)}): " + ", ".join(sorted(only_now)))
        if only_exp:
            print(f"MISSING expected failures ({len(only_exp)}): " + ", ".join(sorted(only_exp)))
        return 1
    if names:
        print(f"NOT GREEN: {len(names)} failing case(s): " + ", ".join(sorted(names)))
        return 1
    print("GREEN: no FAIL/ERROR")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="parse checkall output into the correctness net")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn, help_ in (
        ("table", cmd_table, "Correctness-net markdown rows"),
        ("names", cmd_names, "failing case ids, one per line"),
        ("check", cmd_check, "exit 0 iff green (or matches --expect)"),
    ):
        p = sub.add_parser(name, help=help_)
        p.add_argument("log", help="checkall log file, or - for stdin")
        if name == "check":
            p.add_argument("--expect", help="file of expected failing case ids (net unchanged check)")
        p.set_defaults(func=fn)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
