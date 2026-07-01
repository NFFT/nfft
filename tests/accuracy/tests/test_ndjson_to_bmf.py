import math

import pytest

from tests.accuracy.ndjson_to_bmf import (
    convert,
    read_ndjson,
    group_key,
    metric_name,
    slug,
)


def _rec(**kw):
    base = {
        "module": "nfft",
        "oracle": "file",
        "openmp": 0,
        "dim": 1,
        "N": [16],
        "M": 8,
        "init": "init_guru ()",
        "trafo": "trafo_direct",
        "accuracy": 1e-14,
        "bound": 1e-13,
        "ok": 1,
    }
    base.update(kw)
    return base


def test_slug_normalizes_init_name():
    assert slug("init_guru (PRE PSI)") == "init_guru_PRE_PSI"


def test_group_key_collapses_N_and_M():
    a = group_key(_rec(N=[16], M=8))
    b = group_key(_rec(N=[64], M=99))
    assert a == b  # N and M are bound-absorbed -> same group


def test_speed_and_direction_derived_from_trafo():
    assert (
        metric_name(group_key(_rec(trafo="trafo")))
        == "nfft/serial/file/fast/forward/1d/init_guru"
    )
    assert (
        metric_name(group_key(_rec(trafo="adjoint_direct")))
        == "nfft/serial/file/direct/adjoint/1d/init_guru"
    )


def test_serial_and_omp_are_separate_metrics():
    bmf = convert([_rec(openmp=0), _rec(openmp=1)])
    assert set(bmf) == {
        "nfft/serial/file/direct/forward/1d/init_guru",
        "nfft/omp/file/direct/forward/1d/init_guru",
    }


def test_missing_openmp_field_treated_as_serial():
    rec = _rec()
    del rec["openmp"]
    assert metric_name(group_key(rec)).split("/")[1] == "serial"


def test_convert_takes_max_error_and_digits():
    recs = [
        _rec(N=[16], accuracy=1e-14),  # smaller error
        _rec(N=[64], accuracy=5e-13),
    ]  # worst error -> drives both measures
    bmf = convert(recs)
    ((name, measures),) = bmf.items()
    assert measures["max-error"]["value"] == pytest.approx(5e-13)
    # accuracy-digits is -log10 of the WORST error (fewest digits)
    assert measures["accuracy-digits"]["value"] == pytest.approx(-math.log10(5e-13))


def test_zero_error_uses_low_floor_so_margin_clamps_high():
    # err == 0 floors well below any bound, so the margin stays large-positive,
    # never a false red.
    measures = convert([_rec(accuracy=0.0)]).popitem()[1]
    assert measures["max-error"]["value"] == 0.0
    assert measures["accuracy-digits"]["value"] == pytest.approx(300.0)  # low floor


def test_emits_bound_digits():
    measures = convert([_rec(accuracy=1e-14, bound=1e-13)]).popitem()[1]
    assert measures["bound-digits"]["value"] == pytest.approx(13.0)


def test_bound_digits_aligned_with_worst_error():
    # The colored margin is worst-case; the bound that matters is the one paired
    # with the worst (largest) error, not some other record's tighter bound.
    recs = [
        _rec(N=[16], accuracy=1e-14, bound=1e-16),  # better err, tight bound
        _rec(N=[64], accuracy=5e-13, bound=1e-12),
    ]  # worst err -> its bound
    measures = convert(recs).popitem()[1]
    assert measures["bound-digits"]["value"] == pytest.approx(12.0)


def test_sub_floor_error_keeps_margin_nonnegative():
    # err far below bound but smaller than the OLD 1e-30 floor: the lowered floor
    # must keep margin = accuracy-digits - bound-digits correct (>= 0 here).
    m = convert([_rec(accuracy=1e-35, bound=1e-32)]).popitem()[1]
    margin = m["accuracy-digits"]["value"] - m["bound-digits"]["value"]
    assert margin == pytest.approx(3.0)  # log10(1e-32 / 1e-35)


def test_file_and_online_are_separate_metrics():
    bmf = convert([_rec(oracle="file"), _rec(oracle="online")])
    assert len(bmf) == 2


def test_negative_accuracy_raises():
    with pytest.raises(ValueError, match="accuracy"):
        convert([_rec(accuracy=-1e-13)])


def test_read_ndjson_skips_blank_and_reports_bad_line():
    assert len(read_ndjson('{"module":"nfft"}\n\n')) == 1
    with pytest.raises(ValueError, match="line 1"):
        read_ndjson("not json\n")
