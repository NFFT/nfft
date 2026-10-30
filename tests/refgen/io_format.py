"""Serialization of testcases into the exact CUnit reference-data format.

File layout (blank line between sections, one scalar per line):
    d / N[0..d-1] / M / x (M*d node-major) / f_hat (NN) / f (M)
NFFT scalars are complex ("re im"); NFCT/NFST are real (one value).

Inputs (the nodes, and f_hat for a trafo or f for an adjoint case) are floats and
are written with the fewest digits (at most 9) that read back as the same float.
The computed reference is written with ndig significant digits.
"""
import struct
from fractions import Fraction

import mpmath


def fmt_scalar(value, ndig, is_complex, is_input=False):
    """Format one scalar. Complex -> 're im'; real -> 're'. Fixed-point, no exponent.
    An input is a float and is written as the shortest decimal that reads back
    as that float; otherwise ndig significant digits are written."""
    fmt = _fmt_float if is_input else (lambda v: _fmt_real(v, ndig))
    if is_complex:
        return "%s %s" % (fmt(value.real), fmt(value.imag))
    return fmt(value)


def _fmt_real(v, ndig, strip_zeros=False):
    # mpmath.nstr with no exponent: values here are O(1)..O(M), safe as fixed-point.
    s = mpmath.nstr(mpmath.mpf(v), ndig, strip_zeros=strip_zeros, min_fixed=-mpmath.inf,
                    max_fixed=mpmath.inf)
    return s


def _f32_bits(q):
    return struct.unpack("<I", struct.pack("<f", q))[0]


def _f32_from_bits(b):
    return Fraction(struct.unpack("<f", struct.pack("<I", b))[0])


def _fmt_float(v):
    """Shortest decimal (at most 9 significant digits) that a correctly rounding
    decimal-to-float conversion, such as C's scanf("%f"), maps back to the float v.
    The C reader reads inputs as float and widens them, so every build precision
    gets v exactly."""
    man, exp = mpmath.mpf(v).man_exp
    a = abs(Fraction(man) * Fraction(2) ** exp)
    b = _f32_bits(float(a))
    if _f32_from_bits(b) != a:
        raise ValueError("input %s is not a single-precision value" % v)
    # Half-way points to the neighbouring floats. Requiring the decimal to lie
    # strictly inside them avoids relying on the tie-breaking rule.
    lo = (a + (_f32_from_bits(b - 1) if b > 0 else -_f32_from_bits(1))) / 2
    hi = (a + _f32_from_bits(b + 1)) / 2
    for ndig in range(1, 10):
        s = _fmt_real(v, ndig, strip_zeros=True)
        if lo < abs(Fraction(s)) < hi:
            return s
    raise AssertionError("no decimal of at most 9 digits reads back as %s" % v)


def write_testcase(path, d, N, M, x, f_hat, f, is_complex, ndig, input_is_f_hat=None):
    # The nodes x are always an input; of the coefficient vectors, the input is f_hat
    # for a trafo case and f for an adjoint case (input_is_f_hat True/False); None
    # writes both with ndig digits. The computed output is written with ndig digits.
    in_f_hat = input_is_f_hat is True
    in_f = input_is_f_hat is False
    lines = []
    lines.append(str(d))
    lines.append("")
    for n in N:
        lines.append(str(n))
    lines.append("")
    lines.append(str(M))
    lines.append("")
    for xj in x:  # node-major: node j's d coords contiguous
        for coord in xj:
            lines.append(_fmt_float(coord))
    lines.append("")
    for v in f_hat:
        lines.append(fmt_scalar(v, ndig, is_complex, is_input=in_f_hat))
    lines.append("")
    for v in f:
        lines.append(fmt_scalar(v, ndig, is_complex, is_input=in_f))
    lines.append("")
    with open(path, "w") as fp:
        fp.write("\n".join(lines) + "\n")


def parse_testcase(path, is_complex):
    """Inverse of write_testcase (used by tests; mirrors C setup_file token order)."""
    with open(path) as fp:
        toks = fp.read().split()
    it = iter(toks)
    d = int(next(it))
    N = [int(next(it)) for _ in range(d)]
    M = int(next(it))
    x = []
    for _ in range(M):
        x.append(tuple(mpmath.mpf(next(it)) for _ in range(d)))

    def read_vec(count):
        out = []
        for _ in range(count):
            if is_complex:
                re = mpmath.mpf(next(it))
                im = mpmath.mpf(next(it))
                out.append(mpmath.mpc(re, im))
            else:
                out.append(mpmath.mpf(next(it)))
        return out

    # NN is whatever remains minus M values; compute from remaining token count.
    per = 2 if is_complex else 1
    remaining = sum(1 for _ in iter(lambda: next(it, None), None))
    # re-tokenize because the generator above exhausted `it`; simpler: re-split.
    toks2 = open(path).read().split()
    base = 1 + d + 1 + M * d
    tail = toks2[base:]
    total_scalars = len(tail) // per
    nn = total_scalars - M
    it2 = iter(tail)

    def read_vec2(count):
        out = []
        for _ in range(count):
            if is_complex:
                out.append(mpmath.mpc(mpmath.mpf(next(it2)), mpmath.mpf(next(it2))))
            else:
                out.append(mpmath.mpf(next(it2)))
        return out

    f_hat = read_vec2(nn)
    f = read_vec2(M)
    return d, N, M, x, f_hat, f
