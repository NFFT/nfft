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

# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///

"""Checks for update_copyright.py. Run: uv run support/test_update_copyright.py"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import update_copyright as uc


def test_format_years():
    assert uc.format_years({2026}) == "2026"
    assert (
        uc.format_years({2002, 2003, 2004, 2007, 2009, 2010})
        == "2002-2004, 2007, 2009-2010"
    )


def test_is_trivial():
    old = " * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts"
    new = " * Copyright (c) 2026 The NFFT Authors"
    assert uc.is_trivial([old], [new])
    assert uc.is_trivial(["  x = 1;"], ["x=1;", ""])
    assert not uc.is_trivial(["x = 1;"], ["x = 2;"])
    assert not uc.is_trivial([], ["  1.5, 2.5,"])


def test_update_c_and_matlab():
    c = (
        ["/*"]
        + uc.render(" * ", "2002, 2017")[:1]
        + [" *" if not l else " * " + l for l in uc.BODY]
        + [" */", "int x;"]
    )
    c[1] = " * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts"
    out, problems = uc.update(c, "2005-2026")
    assert not problems
    assert out[1] == " * Copyright (c) 2005-2026 The NFFT Authors"
    assert out[-2:] == [" */", "int x;"] and len(out) == len(c)
    m = [
        "function f()",
        "%F help",
        "%   Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts",
        "",
        "x = 1;",
    ]
    out, _ = uc.update(m, "2009")
    assert out[2] == "%   Copyright (c) 2009 The NFFT Authors" and out[3:] == m[3:]


def test_changed_license_is_kept():
    c = [
        " * Copyright (c) 2017 The NFFT Authors",
        " *",
        " * This program is free software under version 3.",
    ]
    out, problems = uc.update(c, "2026")
    assert out == c and problems


def test_add():
    py = ["#!/usr/bin/env python3", '"""Doc."""']
    out = uc.add("x/a.py", py, "2026")
    assert out[0] == py[0] and out[1] == "# Copyright (c) 2026 The NFFT Authors"
    assert out[-1] == py[1]
    m = ["function f()", "%F help", "x = 1;"]
    out = uc.add("x/f.m", m, "2026")
    assert out[:3] == ["function f()", "%F help", ""] and out[-1] == "x = 1;"
    assert uc.add("x/t.h.in", ["int x;"], "2026")[0] == "/*"


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
    print("ok")
