"""Tests for the polynomial Kaiser-Bessel degree generator.

uv run --with mpmath==1.3.0 --with pytest python -m pytest gen/kbpolygen/tests -q
"""

from gen.kbpolygen import generate


def setup_module(_module):
    generate._init()


def test_degrees_at_known_points():
    # sigma 2, m 4, reach m + 1: truncation-limited in every format but float.
    assert generate.degrees(2.0, 4) == {24: 7, 64: 11, 113: 11, 53: 11}


def test_m_opt_matches_measured_floor():
    # Measured knees of the transform error: double 10 at sigma 1.1, 8 at 2.
    assert generate.m_opt(1.1, 53) in (10, 11)
    assert generate.m_opt(2.0, 53) in (8, 9)
