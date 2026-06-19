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


from tests.refgen import registration as REG


def test_header_contains_decls_and_arrays():
    h = REG.render_header("nfft")
    assert "data/nfft_1d_1_1.txt" in h
    assert "testcase_delegate_file_t nfft_1d_1_1 =" in h
    assert "*testcases_1d_file[]" in h
    assert "*testcases_adjoint_1d_file[]" in h
    assert "*testcases_3d_file[]" in h
    # header guard
    assert "#ifndef" in h and "#endif" in h


def test_extra_dist_lists_all_modules():
    txt = REG.render_extra_dist()
    assert txt.startswith("EXTRA_DIST =")
    assert "nfft_1d_1_1.txt" in txt
    assert "nfct_2d_10_25_50.txt" in txt
    assert "nfst_adjoint_3d_10_10_10_10.txt" in txt


def test_quad_precision_digits(tmp_path):
    import mpmath
    from tests.refgen import generate as GEN
    # Generate one nfft file at quad-target precision and check digit count.
    GEN.main(["--module", "nfft", "--precision", "40",
              "--data-dir", str(tmp_path), "--header-dir", str(tmp_path)])
    # Use a file whose output is a genuine multi-term sum (N=10 > 1).
    p = tmp_path / "nfft_1d_10_1.txt"
    toks = p.read_text().split()
    # The computed *output* (last scalar of the f section) must carry >=34
    # significant digits (quad-sufficient). Inputs (nodes/coefficients) are now
    # drawn as floats and written compact, so the digit budget lives in the output.
    out = toks[-1]
    sig = out.replace("-", "").replace(".", "").lstrip("0")
    assert len(sig) >= 34
