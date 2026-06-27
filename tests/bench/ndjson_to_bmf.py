"""Convert per-case accuracy NDJSON (emitted by the CUnit harness via
NFFT_BENCH_OUT) into aggregated Bencher Metric Format (BMF) JSON.

Each input line is one raw case:
    {"module","oracle","dim","N","M","init","trafo","accuracy","bound","ok"}

Output is one BMF object per *accuracy metric* -- a combination of the
error-shaping parameters (module, oracle file/online, speed direct/fast,
direction forward/adjoint, dimension, init variant) -- with the bound-absorbed
parameters (N, M) collapsed via max:

    {"<module>/<oracle>/<speed>/<direction>/<dim>d/<init-slug>": {
        "tightness-ratio": {"value": max(err/bound)},   # primary
        "max-error":       {"value": max(err)}}, ...}    # secondary
"""
from __future__ import annotations

import argparse
import json
import re
import sys


def slug(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")


def group_key(rec):
    trafo = rec["trafo"]
    speed = "direct" if "direct" in trafo else "fast"
    direction = "adjoint" if trafo.startswith("adjoint") else "forward"
    return (rec["module"], rec["oracle"], speed, direction,
            int(rec["dim"]), slug(rec["init"]))


def metric_name(key):
    module, oracle, speed, direction, dim, init = key
    return f"{module}/{oracle}/{speed}/{direction}/{dim}d/{init}"


def convert(records):
    agg = {}  # key -> [max_ratio, max_err]
    for rec in records:
        bound = float(rec["bound"])
        if bound <= 0.0:
            raise ValueError(f"non-positive bound in record: {rec!r}")
        err = float(rec["accuracy"])
        ratio = err / bound
        key = group_key(rec)
        if key not in agg:
            agg[key] = [ratio, err]
        else:
            agg[key][0] = max(agg[key][0], ratio)
            agg[key][1] = max(agg[key][1], err)
    return {
        metric_name(key): {
            "tightness-ratio": {"value": r},
            "max-error": {"value": e},
        }
        for key, (r, e) in agg.items()
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
