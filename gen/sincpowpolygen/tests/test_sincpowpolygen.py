"""Tests for the polynomial sinc-power degree generator.

uv run --with mpmath==1.3.0 --with pytest python -m pytest gen/sincpowpolygen/tests -q
"""

from gen.sincpowpolygen import generate


def setup_module(_module):
    generate._init()


def test_degrees_at_known_points():
    assert generate.degrees(2.0, 4) == {24: 7, 64: 17, 113: 27, 53: 14}
    assert generate.degrees(2.0, 32) == {24: 5, 64: 13, 113: 21, 53: 11}


def test_row_takes_running_maximum():
    # The exact float degree of m = 2 dips from 8 to 7 near sigma 1.5; the
    # row stays at 8 there.
    _, (d0, steps) = generate.row_steps(2)
    assert d0[24] == 7
    assert steps[24] == [13891, 22871]
