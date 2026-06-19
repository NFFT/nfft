"""Serialization of testcases into the exact CUnit reference-data format.

File layout (blank line between sections, one scalar per line):
    d / N[0..d-1] / M / x (M*d node-major) / f_hat (NN) / f (M)
NFFT scalars are complex ("re im"); NFCT/NFST are real (one value).
"""
import mpmath


def fmt_scalar(value, ndig, is_complex, strip_zeros=False):
    """Format one scalar. Complex -> 're im'; real -> 're'. Fixed-point, no exponent."""
    if is_complex:
        return "%s %s" % (_fmt_real(value.real, ndig, strip_zeros),
                          _fmt_real(value.imag, ndig, strip_zeros))
    return _fmt_real(value, ndig, strip_zeros)


def _fmt_real(v, ndig, strip_zeros=False):
    # mpmath.nstr with no exponent: values here are O(1)..O(M), safe as fixed-point.
    # strip_zeros=True drops trailing zeros: inputs are drawn as floats, so their
    # exact decimal is short and padding to ndig is pure waste (lossless to strip).
    s = mpmath.nstr(mpmath.mpf(v), ndig, strip_zeros=strip_zeros, min_fixed=-mpmath.inf,
                    max_fixed=mpmath.inf)
    return s


def write_testcase(path, d, N, M, x, f_hat, f, is_complex, ndig, input_is_f_hat=None):
    # Inputs are drawn as floats (exact in every build precision), so their decimals
    # are short — strip trailing zeros to keep the files compact. The high-precision
    # computed *output* keeps the full ndig. The nodes x are always an input; of the
    # coefficient vectors, the input is f_hat for a trafo case and f for an adjoint
    # case (input_is_f_hat True/False); None keeps full precision for both.
    strip_f_hat = input_is_f_hat is True
    strip_f = input_is_f_hat is False
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
            lines.append(_fmt_real(coord, ndig, strip_zeros=True))
    lines.append("")
    for v in f_hat:
        lines.append(fmt_scalar(v, ndig, is_complex, strip_zeros=strip_f_hat))
    lines.append("")
    for v in f:
        lines.append(fmt_scalar(v, ndig, is_complex, strip_zeros=strip_f))
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
