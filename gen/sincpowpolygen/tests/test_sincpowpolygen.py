"""Tests for the polynomial sinc-power degree generator.

uv run --with mpmath==1.3.0 --with pytest python -m pytest gen/sincpowpolygen/tests -q
"""

from gen.sincpowpolygen import generate


def setup_module(_module):
    generate._init()


def test_degrees_at_known_points():
    assert generate.degrees(2.0, 4) == {24: 7, 64: 18, 113: 28, 53: 15}
    assert generate.degrees(2.0, 32) == {24: 5, 64: 13, 113: 22, 53: 11}


def test_row_steps():
    _, (d0, steps) = generate.row_steps(2)
    assert d0 == {24: 8, 64: 18, 113: 29, 53: 16}
    assert steps[24] == [13579]
    assert steps[64] == [12010, 26479]
