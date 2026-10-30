"""Tests for the polynomial Kaiser-Bessel degree generator.

    uv run --with mpmath==1.3.0 --with numpy --with pytest python -m pytest gen/kbpolygen/tests -q
"""
from gen.kbpolygen import generate


def setup_module(_module):
    generate._init()


def test_halves_at_known_points():
    # Half degrees measured through the error functional; float is at its
    # rounding floor from m 4 on.
    assert generate.halves(2.0, 2) == {24: 3, 64: 3, 113: 3, 53: 3}
    assert generate.halves(2.0, 4) == {24: 4, 64: 5, 113: 5, 53: 5}
    assert generate.halves(4.0, 2) == {24: 3, 64: 4, 113: 4, 53: 4}


def test_guard_admits_measured_range():
    # Where the plan-time correction still pays: double to m 6 at sigma 2
    # and m 8 at sigma 1.25, float to m 2 at sigma 2 (measured in the
    # kbpoly-optimize analysis; above these the error is the rounding floor).
    import math
    b = lambda s: math.pi * (2 - 1 / s)
    assert int(generate.br_ls(53) / b(2.0)) - 1 == 6
    assert int(generate.br_ls(53) / b(1.25)) - 1 == 8
    assert int(generate.br_ls(24) / b(2.0)) - 1 == 2


def test_m_opt_matches_measured_floor():
    # Measured knees of the transform error: double 10 at sigma 1.1, 8 at 2.
    assert generate.m_opt(1.1, 53) in (10, 11)
    assert generate.m_opt(2.0, 53) in (8, 9)
