"""Convert per-case accuracy NDJSON (emitted by the CUnit harness via
NFFT_BENCH_OUT) into aggregated Bencher Metric Format (BMF) JSON.

Each input line is one raw case:
    {"module","oracle","openmp","dim","N","M","init","trafo","accuracy","bound","ok"}

Output is one BMF object per *accuracy metric* -- a combination of the
error-shaping parameters (module, runtime serial/omp, oracle file/online, speed
direct/fast, direction forward/adjoint, dimension, init variant) -- with the
bound-absorbed parameters (N, M) collapsed via max:

    {"<module>/<runtime>/<oracle>/<speed>/<direction>/<dim>d/<init-slug>": {
        "accuracy-digits": {"value": -log10(max(err))},  # primary (higher=better)
        "max-error":       {"value": max(err)}}, ...}     # secondary (raw worst err)

The primary measure is the worst-case *accurate digits* (-log10 of the worst
relative error). Raw errors span ~14 orders of magnitude (1e-17 .. 1e-3), which
no linear scale displays well; the log transform reads cleanly (~3 .. 18) and a
regression lowers it -> use a Bencher lower-boundary threshold. The exact worst
error is kept as the secondary `max-error` measure.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys


def slug(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")


def group_key(rec):
    trafo = rec["trafo"]
    speed = "direct" if "direct" in trafo else "fast"
    direction = "adjoint" if trafo.startswith("adjoint") else "forward"
    # `runtime` (serial vs OpenMP) is error-shaping: the parallel reduction order
    # perturbs the low bits, so the two builds get distinct metrics. Records
    # without an `openmp` field are treated as serial.
    runtime = "omp" if rec.get("openmp") else "serial"
    return (rec["module"], runtime, rec["oracle"], speed, direction,
            int(rec["dim"]), slug(rec["init"]))


def metric_name(key):
    module, runtime, oracle, speed, direction, dim, init = key
    return f"{module}/{runtime}/{oracle}/{speed}/{direction}/{dim}d/{init}"


# Floor for the -log10 transform so an exact-zero (perfect) error stays finite;
# any error <= this is reported as ~30 accurate digits.
_ERR_FLOOR = 1e-30


def _measures(max_err):
    # accuracy-digits = -log10(worst error): readable across the ~14 orders of
    # magnitude the raw error spans, and higher = better (a regression lowers it).
    digits = -math.log10(max(max_err, _ERR_FLOOR))
    return {
        "accuracy-digits": {"value": digits},  # primary
        "max-error": {"value": max_err},        # secondary (exact worst error)
    }


def convert(records):
    agg = {}  # key -> max_err (worst error over the bound-absorbed N/M)
    for rec in records:
        err = float(rec["accuracy"])
        if err < 0.0:
            raise ValueError(f"negative accuracy in record: {rec!r}")
        key = group_key(rec)
        if key not in agg or err > agg[key]:
            agg[key] = err
    return {metric_name(key): _measures(max_err) for key, max_err in agg.items()}


def read_ndjson(text):
    records = []
    for lineno, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid JSON on line {lineno}: {exc}") from exc
    return records


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", help="NDJSON file written via NFFT_BENCH_OUT")
    ap.add_argument("output", help="BMF JSON output path")
    args = ap.parse_args(argv)

    with open(args.input, encoding="utf-8") as f:
        records = read_ndjson(f.read())
    bmf = convert(records)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(bmf, f, indent=2, sort_keys=True)
    print(f"wrote {len(bmf)} accuracy metrics to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
