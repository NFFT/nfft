"""CLI: regenerate reference data + generated C headers + the data EXTRA_DIST list.

    python -m tests.refgen.generate --module all --precision 64
"""
import argparse
import os
import random
import struct

import mpmath

from tests.refgen import grids as G
from tests.refgen import io_format as IO
from tests.refgen import registration as REG
from tests.refgen import transforms as T

_MODULES = ("nfft", "nfct", "nfst")


def _is_complex(module):
    return module == "nfft"


def _f32(v):
    """Round a Python float (double) to IEEE-754 single precision, returned as the
    exact double that holds that float value."""
    return struct.unpack("f", struct.pack("f", v))[0]


def _draw_inputs(module, kind, d, N, M, rng, variant=None):
    """Reproducible inputs. Nodes and coefficients are drawn as single-precision
    (float) values: a float is exactly representable in float, double, long double,
    and __float128, so the C reader reproduces the input bit-for-bit at *every*
    precision — including the float build — leaving no input-rounding term (only
    summation error). Summation itself is done at the configured high precision."""
    lo, hi = (-0.5, 0.5) if module == "nfft" else (0.0, 0.5)
    x = [tuple(mpmath.mpf(_f32(rng.uniform(lo, hi))) for _ in range(d)) for _ in range(M)]
    count = T.nn(module, N, variant) if kind == "trafo" else M
    if _is_complex(module):
        coeff = [mpmath.mpc(_f32(rng.uniform(-1, 1)), _f32(rng.uniform(-1, 1))) for _ in range(count)]
    else:
        coeff = [mpmath.mpf(_f32(rng.uniform(-1, 1))) for _ in range(count)]
    return x, coeff


def _gen_module(module, precision, seed, data_dir):
    mpmath.mp.dps = precision
    is_c = _is_complex(module)
    for kind in G.KINDS:
        for (d, N, M, variant) in G.GRIDS[module]:
            # Deterministic per-file seed so a single file regenerates identically
            # regardless of iteration order. random.Random(str) uses a stable
            # internal hash (NOT PYTHONHASHSEED), so this is reproducible across
            # processes — unlike the builtin hash() of a tuple.
            rng = random.Random("%d:%s" % (seed, G.basename(module, kind, d, N, M, variant)))
            x, coeff = _draw_inputs(module, kind, d, N, M, rng, variant)
            if kind == "trafo":
                f_hat = coeff
                f = T.trafo(module, N, M, x, f_hat, variant)
            else:
                f = coeff
                f_hat = T.adjoint(module, N, M, x, f, variant)
            path = os.path.join(data_dir, G.basename(module, kind, d, N, M, variant) + ".txt")
            # trafo: f_hat is the (float-drawn) input, f the computed output;
            # adjoint: the roles swap. Inputs get trailing zeros stripped.
            IO.write_testcase(path, d, N, M, x, f_hat, f,
                              is_complex=is_c, ndig=precision,
                              input_is_f_hat=(kind == "trafo"))
            print("wrote", path)


def _gen_ndft_accuracy(precision, seed, data_dir):
    """The NFFT_NDFT_ACCURACY trafo case: nodes are the float rounding of
    0.4999 - 0.0001*(j%7), f_hat float-drawn from the seeded RNG.

    Working precision is precision + 24 digits. With |k x| <= 4096 the phase
    2 pi k x has absolute error below 2.6e4 * 10^-dps, so each of the N = 8192
    terms (|f_hat| <= sqrt(2)) is off by at most 3.7e4 * 10^-dps and the sum,
    rounding included, by less than 5e8 * 10^-dps. That leaves 15 guard digits
    above the printed precision for |f[j]| >= 1 (this case has min |f[j]| > 30)."""
    d, N, M = G.NFFT_NDFT_ACCURACY
    base = G.basename("nfft", "trafo", d, N, M)
    rng = random.Random("%d:%s" % (seed, base))
    x = [(mpmath.mpf(_f32(0.4999 - 0.0001 * (j % 7))),) for j in range(M)]
    f_hat = [mpmath.mpc(_f32(rng.uniform(-1, 1)), _f32(rng.uniform(-1, 1)))
             for _ in range(T.nn("nfft", N))]
    with mpmath.workdps(precision + 24):
        f = T.trafo("nfft", N, M, x, f_hat)
    path = os.path.join(data_dir, base + ".txt")
    IO.write_testcase(path, d, N, M, x, f_hat, f, is_complex=True,
                      ndig=precision, input_is_f_hat=True)
    print("wrote", path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--module", choices=("nfft", "nfct", "nfst", "all"), default="all")
    ap.add_argument("--precision", type=int, default=64,
                    help="working/printed significant decimal digits (>=34 for quad)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--data-dir", default=os.path.join("tests", "data"))
    ap.add_argument("--header-dir", default=os.path.join("tests", "data", "generated"))
    args = ap.parse_args(argv)

    modules = _MODULES if args.module == "all" else (args.module,)
    os.makedirs(args.header_dir, exist_ok=True)
    for m in modules:
        _gen_module(m, args.precision, args.seed, args.data_dir)
        hpath = os.path.join(args.header_dir, "%s_testcases.h" % m)
        with open(hpath, "w") as fp:
            fp.write(REG.render_header(m))
        print("wrote", hpath)
        if m == "nfft":
            _gen_ndft_accuracy(args.precision, args.seed, args.data_dir)
            npath = os.path.join(args.header_dir, "nfft_native_testcases.h")
            with open(npath, "w") as fp:
                fp.write(REG.render_native_header())
            print("wrote", npath)

    # EXTRA_DIST always reflects all modules (the file lists every module's data).
    mk = os.path.join(args.data_dir, "Makefile.am")
    with open(mk, "w") as fp:
        fp.write(REG.render_extra_dist())
    print("wrote", mk)


if __name__ == "__main__":
    main()
