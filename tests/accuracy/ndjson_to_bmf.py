"""Convert per-case accuracy NDJSON into the aggregated accuracy JSON consumed 
by the HTML report renderer.

Each input line is one raw case:
    {"module","oracle","openmp","dim","N","M","init","trafo","accuracy","bound","ok"}

Output is one BMF object per *accuracy metric* -- the error-shaping parameters
(module, runtime serial/omp, oracle file/online, speed direct/fast, direction
forward/adjoint, dimension, init variant) -- with the bound-absorbed parameters
(N, M) collapsed via max:

    {"<module>/<runtime>/<oracle>/<speed>/<direction>/<dim>d/<init-slug>": {
        "accuracy-digits": {"value": -log10(max(err))},  # primary (higher=better)
        "max-error":       {"value": max(err)},           # secondary (raw worst err)
        "bound-digits":    {"value": -log10(bound)}}, ...} # worst case's bound

Errors are reported as digits (-log10) because the raw values span ~14 orders of
magnitude (1e-17 .. 1e-3), which no linear color scale renders usefully.
`bound-digits` carries the worst-error record's analytic bound so the heatmap can
color by margin (accuracy-digits - bound-digits = digits beyond the bound), which
normalizes away precision and per-case difficulty.
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
    return (
        rec["module"],
        runtime,
        rec["oracle"],
        speed,
        direction,
        int(rec["dim"]),
        slug(rec["init"]),
    )


def metric_name(key):
    module, runtime, oracle, speed, direction, dim, init = key
    return f"{module}/{runtime}/{oracle}/{speed}/{direction}/{dim}d/{init}"


# Floor for the -log10 transform so an exact-zero error stays finite. It must sit
# below any representable bound (bounds reach ~1e-32) so the margin stays >= 0 for
# a passing case: flooring is monotonic, so err <= bound => margin >= 0. err == 0
# then maps to ~300 digits (darkest green; heatmap shows the cell as infinity).
_ERR_FLOOR = 1e-300


def _measures(max_err, bound):
    digits = -math.log10(max(max_err, _ERR_FLOOR))
    bound_digits = -math.log10(max(bound, _ERR_FLOOR))
    return {
        "accuracy-digits": {"value": digits},  # primary
        "max-error": {"value": max_err},  # secondary (exact worst error)
        "bound-digits": {"value": bound_digits},  # the worst case's bound
    }


def convert(records):
    agg = {}  # key -> (max_err, bound) of the worst-error record in the group
    for rec in records:
        err = float(rec["accuracy"])
        if err < 0.0:
            raise ValueError(f"negative accuracy in record: {rec!r}")
        bound = float(rec["bound"])
        key = group_key(rec)
        if key not in agg or err > agg[key][0]:
            agg[key] = (err, bound)
    return {
        metric_name(key): _measures(max_err, bound)
        for key, (max_err, bound) in agg.items()
    }


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
