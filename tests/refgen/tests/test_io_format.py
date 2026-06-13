import mpmath

from tests.refgen import io_format as IO


def test_format_real_has_enough_digits():
    mpmath.mp.dps = 64
    s = IO.fmt_scalar(mpmath.mpf(1) / 3, ndig=50, is_complex=False)
    # at least 40 significant digits, no exponent for an O(1) value
    digits = s.replace("-", "").replace(".", "").lstrip("0")
    assert len(digits) >= 40
    assert "e" not in s and "E" not in s


def test_format_complex_two_columns():
    mpmath.mp.dps = 64
    s = IO.fmt_scalar(mpmath.mpc("0.5", "-0.25"), ndig=40, is_complex=True)
    parts = s.split()
    assert len(parts) == 2
    assert parts[0].startswith("0.5")
    assert parts[1].startswith("-0.25")


def test_write_then_parse_roundtrip_real(tmp_path):
    mpmath.mp.dps = 64
    d, N, M = 1, [4], 2
    x = [(mpmath.mpf("0.1"),), (mpmath.mpf("0.2"),)]
    f_hat = [mpmath.mpf("0.3"), mpmath.mpf("0.4"), mpmath.mpf("0.5")]  # nfst N-1=3
    f = [mpmath.mpf("0.6"), mpmath.mpf("0.7")]
    p = tmp_path / "nfst_1d_4_2.txt"
    IO.write_testcase(str(p), d, N, M, x, f_hat, f, is_complex=False, ndig=50)
    rd, rN, rM, rx, rfh, rf = IO.parse_testcase(str(p), is_complex=False)
    assert (rd, rN, rM) == (d, N, M)
    assert abs(rx[0][0] - x[0][0]) < mpmath.mpf(10) ** -40
    assert abs(rfh[2] - f_hat[2]) < mpmath.mpf(10) ** -40
    assert abs(rf[1] - f[1]) < mpmath.mpf(10) ** -40


def test_write_then_parse_roundtrip_complex(tmp_path):
    mpmath.mp.dps = 64
    d, N, M = 1, [1], 2
    x = [(mpmath.mpf("-0.1"),), (mpmath.mpf("0.3"),)]
    f_hat = [mpmath.mpc("0.2", "0.1")]
    f = [mpmath.mpc("0.4", "-0.5"), mpmath.mpc("-0.6", "0.7")]
    p = tmp_path / "nfft_1d_1_2.txt"
    IO.write_testcase(str(p), d, N, M, x, f_hat, f, is_complex=True, ndig=50)
    rd, rN, rM, rx, rfh, rf = IO.parse_testcase(str(p), is_complex=True)
    assert (rd, rN, rM) == (d, N, M)
    assert abs(rf[0] - f[0]) < mpmath.mpf(10) ** -40
