"""Tests for the optimized window generator.

    uv run --with python-flint==0.9.0 --with numpy --with mpmath==1.3.0 --with pytest \
        python -m pytest gen/kbpolygen/tests -q
"""

import re

import pytest
from flint import arb

from gen.kbpolygen import optima


@pytest.mark.parametrize(
    "sigma, ref",
    [(1.25, 1.058e-5), (1.5, 6.496e-7), (2.0, 2.183e-8)],
)
def test_reference_rows(sigma, ref):
    # Appendix A of the plan: paper m 4, degree 10 is library m 3, h 5.
    e = float(optima.solve_cell(sigma, 3, 5)["e"])
    assert 0.97 * ref <= e <= 1.03 * ref


def test_header_roundtrip():
    res = optima.solve_cell(2.0, 1, 3)
    text = "\n".join(optima.cell_text(res, 0))
    coef = re.search(r"coef_1_0\[\] = \{(.*?)\};", text, re.DOTALL).group(1)
    vals = re.findall(r"K\(([^)]*)\)", coef)
    assert vals == [optima.fmt(v) for v in res["coef"]]
    assert all(
        len(v.lstrip("-").replace(".", "").split("e")[0].lstrip("0")) <= 36
        for v in vals
    )
    for v, c in zip(vals, res["coef"]):
        assert abs(arb(v) - c) <= 1e-35 * abs(c)
    pts = re.search(r"pts_1_0\[\] = \{(.*?)\};", text, re.DOTALL).group(1)
    flat = [int(v) for v in re.findall(r"\d+", pts)]
    pairs = list(zip(flat[::2], flat[1::2]))
    assert len(flat) % 2 == 0 and pairs == sorted(set(pairs)) == res["pts"]
    assert all(0 <= v <= 65535 for v in flat)
