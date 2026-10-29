# Copyright (c) 2026 The NFFT Authors
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, see <https://www.gnu.org/licenses/>.

"""Faithful fixed-precision emulation of the C runtime evaluation.

mpmath rounds every operation to the working precision with round-to-nearest,
ties-to-even, so evaluating at prec = MANT_DIG reproduces the C arithmetic of
the matching IEEE format. This is what covers the formats a given host cannot
compile -- quadruple everywhere, and Intel double extended on a machine whose
`long double` is binary128.
"""
import mpmath as mp


def horner(c, x):
    """The Horner chain as written in kernel/util/sinc.c."""
    r = c[-1]
    for j in range(len(c) - 2, -1, -1):
        r = r * x + c[j]
    return r


def log_sinc(x, prec, Q, split):
    """Candidate runtime evaluation of log|sinc(x)| at `prec` mantissa bits."""
    with mp.workprec(prec):
        a = +abs(x)
        if a <= split:
            y = a * a
            return y * horner(Q, y)
        return mp.log(abs(mp.sin(a) / a))


def ulp_error(got, ref, prec):
    if ref == 0:
        return abs(got)
    e = mp.frexp(abs(ref))[1]
    return abs(got - ref) / mp.mpf(2) ** (int(e) - prec)
