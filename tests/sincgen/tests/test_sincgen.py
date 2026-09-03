"""Tests for the log|sinc| coefficient generator.

    uv run --with mpmath==1.3.0 --with pytest python -m pytest tests/sincgen/tests -q
"""
import mpmath as mp
import pytest

from tests.sincgen import emulate, generate, scheme


def setup_module(_module):
    mp.mp.dps = 40


def _horner(c, x):
    r = c[-1]
    for j in range(len(c) - 2, -1, -1):
        r = r * x + c[j]
    return r


def test_fq_extends_continuously_to_zero():
    assert scheme.fq(mp.mpf(0)) == -mp.mpf(1) / 6
    y = mp.mpf(10) ** -12
    series = -mp.mpf(1) / 6 - y / 180
    assert abs(scheme.fq(y) - series) < mp.mpf(10) ** -27


def test_coefficients_start_at_the_taylor_series():
    # log sinc(x) = -x^2/6 - x^4/180 - x^6/2835 - x^8/37800 - ...
    Q, _ = scheme.fit_branch(2, 19)
    want = [-mp.mpf(1) / 6, -mp.mpf(1) / 180, -mp.mpf(1) / 2835,
            -mp.mpf(1) / 37800]
    # The minimax correction rides on top of the Taylor coefficients and grows
    # with the term index, so the match loosens along the series.
    for j, (got, ref) in enumerate(zip(Q, want)):
        assert abs(got - ref) <= abs(ref) * mp.mpf(10) ** -(18 - 3 * j)


def test_coefficients_are_negative_so_horner_cannot_cancel():
    Q, _ = scheme.fit_branch(2, 19)
    # A few trailing coefficients are positive but negligible; what matters is
    # that the evaluated sum never loses digits.
    assert scheme.horner_growth(Q, 2) < mp.mpf("1.001")


def test_the_split_is_where_horner_stops_cancelling():
    # Pinning SPLIT = 2 is a conditioning choice, not an accuracy one: at 2.5
    # the same degree still fits but the Horner chain starts to cancel.
    Q2, _ = scheme.fit_branch(2, 24)
    Q25, _ = scheme.fit_branch(mp.mpf("2.5"), 24)
    assert scheme.horner_growth(Q2, 2) < mp.mpf("1.001")
    assert scheme.horner_growth(Q25, mp.mpf("2.5")) > 10


def test_shipped_form_reaches_the_target_at_every_pinned_degree():
    for prec, spec in generate.SPEC.items():
        # Same working precision the generator uses for this format; 40 digits
        # cannot resolve the quadruple-precision target at all.
        with mp.workdps(60 if prec <= 64 else 90):
            target = mp.mpf(2) ** (-prec - 8)
            _, err = scheme.fit_branch(spec["split"], spec["n"])
            assert err <= target, "MANT_DIG=%d: %s > %s" % (prec, err, target)


def test_verify_branch_measures_the_monomial_form_not_the_fit():
    Q, err = scheme.fit_branch(2, 19)
    worst = mp.mpf(0)
    for i in range(1, 201):
        x = 2 * mp.mpf(i) / 200
        ref = scheme.log_sinc(x)
        worst = max(worst, abs(x * x * _horner(Q, x * x) - ref) / abs(ref))
    assert abs(worst - err) <= err


@pytest.mark.parametrize("prec", [24, 53, 64, 113])
def test_emulated_evaluation_is_accurate_in_the_target_format(prec):
    """Covers the formats a given host cannot compile: quadruple everywhere,
    and Intel double extended wherever `long double` is binary128."""
    spec = generate.SPEC[prec]
    mp.mp.dps = 60 if prec <= 64 else 90
    Q, _ = scheme.fit_branch(spec["split"], spec["n"])
    with mp.workprec(prec):
        Q = [+v for v in Q]
    worst = mp.mpf(0)
    for i in range(1, 300):
        x = mp.mpf(spec["split"]) * mp.mpf("1.5") * i / 300
        with mp.workprec(prec):
            x = +x
        got = emulate.log_sinc(x, prec, Q, spec["split"])
        worst = max(worst, emulate.ulp_error(got, scheme.log_sinc(x), prec))
    assert worst < 8, "max %.2f ulp" % float(worst)
    mp.mp.dps = 40


@pytest.mark.parametrize("prec", [24, 53, 64, 113])
def test_emulated_power_beats_the_pow_form(prec):
    """The point of the block: EXP(2m*log_sinc(t)) must beat POW(sinc(t), 2m)
    over the B-spline PHI_HUT band, |t| <= pi/4 at the default oversampling."""
    spec = generate.SPEC[prec]
    mp.mp.dps = 60 if prec <= 64 else 90
    Q, _ = scheme.fit_branch(spec["split"], spec["n"])
    with mp.workprec(prec):
        Q = [+v for v in Q]
    m = 11
    wlog = wpow = mp.mpf(0)
    for i in range(1, 200):
        t = mp.pi / 4 * i / 200
        with mp.workprec(prec):
            t = +t
        ref = (mp.sin(t) / t) ** (2 * m)
        with mp.workprec(prec):
            vlog = +mp.exp(2 * m * emulate.log_sinc(t, prec, Q, spec["split"]))
            vpow = +((+(mp.sin(t) / t)) ** (2 * m))
        wlog = max(wlog, emulate.ulp_error(vlog, ref, prec))
        wpow = max(wpow, emulate.ulp_error(vpow, ref, prec))
    assert wlog < wpow, "log %.1f ulp, pow %.1f ulp" % (float(wlog), float(wpow))
    mp.mp.dps = 40
