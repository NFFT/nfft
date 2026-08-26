"""Tests for the I0 coefficient generator.

    uv run --with mpmath==1.3.0 --with pytest python -m pytest tests/besselgen/tests -q
"""
import mpmath as mp
import pytest

from tests.besselgen import emulate, generate, remez, scheme


def setup_module(_module):
    mp.mp.dps = 40


def test_clenshaw_matches_direct_chebyshev_sum():
    c = [mp.mpf(s) for s in ("0.5", "-1.25", "2", "0.125", "-3")]
    for t in ("-1", "-0.3", "0", "0.7", "1"):
        t = mp.mpf(t)
        direct = sum(cj * tj for cj, tj in zip(c, remez.cheb_basis(t, len(c) - 1)))
        assert abs(remez.clenshaw(c, t) - direct) < mp.mpf(10) ** -30


def test_cheb_to_monomial_preserves_the_polynomial():
    a, b = mp.mpf(0), mp.mpf("56.25")
    c = [mp.mpf(s) for s in ("1", "0.5", "-0.25", "0.125", "2", "-0.5", "0.75")]
    mono = remez.cheb_to_monomial(c, a, b)
    for i in range(11):
        x = a + (b - a) * mp.mpf(i) / 10
        t = (2 * x - a - b) / (b - a)
        want = remez.clenshaw(c, t)
        got = mono[-1]
        for j in range(len(mono) - 2, -1, -1):
            got = got * x + mono[j]
        assert abs(got - want) <= abs(want) * mp.mpf(10) ** -28


def test_minimax_beats_chebyshev_interpolation_and_equioscillates():
    f = mp.exp
    w = lambda x: 1 / mp.exp(x)          # noqa: E731  -> relative error
    a, b = mp.mpf(0), mp.mpf(2)
    c, err = remez.minimax_poly(f, w, a, b, 6)
    cheb = remez.cheb_series(f, a, b, 6)

    def relerr(coeffs):
        worst = mp.mpf(0)
        for i in range(201):
            x = a + (b - a) * mp.mpf(i) / 200
            t = (2 * x - a - b) / (b - a)
            worst = max(worst, abs(remez.clenshaw(coeffs, t) - f(x)) / f(x))
        return worst

    assert relerr(c) <= relerr(cheb)
    assert abs(relerr(c) - err) <= err * mp.mpf("0.05")


def test_branch1_coefficients_are_positive_and_start_at_the_taylor_series():
    # (I0(2 sqrt y) - 1)/y = 1 + y/4 + y^2/36 + ...
    P1, _ = scheme.fit_branch1(8, 10)
    # The minimax correction rides on top of the Taylor coefficients and grows
    # with the term index, so the match loosens along the series.
    assert all(v > 0 for v in P1[:6])
    assert abs(P1[0] - 1) < mp.mpf("1e-8")
    assert abs(P1[1] - mp.mpf(1) / 4) < mp.mpf("1e-8")
    assert abs(P1[2] - mp.mpf(1) / 36) < mp.mpf("1e-7")


def test_branch2_starts_at_one_over_sqrt_two_pi():
    # g(0) = 1/sqrt(2 pi); the Chebyshev coefficients sum to it at t = -1.
    P2, _ = scheme.fit_branch2(10, 12)
    got = remez.clenshaw(P2, mp.mpf(-1))
    assert abs(got - 1 / mp.sqrt(2 * mp.pi)) < mp.mpf("1e-12")


def test_reported_error_is_the_error_of_the_shipped_form():
    """Regression: the fit runs in the Chebyshev basis on a search interval
    stepped just inside the real one (the branch-1 weight vanishes at y=0), but
    the coefficients returned must describe the original interval. Converting
    over the stepped interval instead left an error floor at the step size,
    invisible at double and fatal at quad."""
    mp.mp.dps = 70
    P1, err = scheme.fit_branch1(10, 24)
    assert err < mp.mpf(10) ** -30, "shipped monomial form is worse than reported"
    assert abs(scheme.verify_branch1(P1, 10) - err) <= err * mp.mpf("0.01")
    mp.mp.dps = 40


@pytest.mark.parametrize("prec", [24, 53])
def test_spec_degrees_meet_the_format(prec):
    """Each pinned degree must leave the design error well under the format
    epsilon, so that what is measured at run time is evaluation rounding."""
    spec = generate.SPEC[prec]
    mp.mp.dps = 50
    _, e1 = scheme.fit_branch1(spec["split"], spec["n1"])
    _, e2 = scheme.fit_branch2(spec["split"], spec["n2"])
    eps = mp.mpf(2) ** -prec
    assert e1 < eps / 8
    assert e2 < eps / 8
    mp.mp.dps = 40


@pytest.mark.parametrize("prec", [24, 53])
def test_emulated_evaluation_is_accurate_in_the_target_format(prec):
    spec = generate.SPEC[prec]
    mp.mp.dps = 50
    P1, _ = scheme.fit_branch1(spec["split"], spec["n1"])
    P2, _ = scheme.fit_branch2(spec["split"], spec["n2"])
    with mp.workprec(prec):
        P1 = [+v for v in P1]
        P2 = [+v for v in P2]
    worst = mp.mpf(0)
    for i in range(1, 200):
        x = mp.mpf(spec["split"]) * 3 * i / 200
        with mp.workprec(prec):
            x = +x
        got = emulate.i0_new(x, prec, P1, P2, spec["split"])
        worst = max(worst, emulate.ulp_error(got, mp.besseli(0, x), prec))
    assert worst < 8, "max %.2f ulp" % float(worst)
    mp.mp.dps = 40


def test_emitted_literals_round_trip():
    v = mp.mpf("1.2345678901234567890123456789")
    s = generate.cnum(v, 20)
    assert "e" in s
    assert abs(mp.mpf(s) - v) <= abs(v) * mp.mpf(10) ** -19
