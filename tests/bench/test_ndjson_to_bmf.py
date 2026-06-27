import pytest

from ndjson_to_bmf import convert, read_ndjson, group_key, metric_name, slug


def _rec(**kw):
    base = {"module": "nfft", "oracle": "file", "dim": 1, "N": [16], "M": 8,
            "init": "init_guru ()", "trafo": "trafo_direct",
            "accuracy": 1e-14, "bound": 1e-13, "ok": 1}
    base.update(kw)
    return base


def test_slug_normalizes_init_name():
    assert slug("init_guru (PRE PSI)") == "init_guru_PRE_PSI"


def test_group_key_collapses_N_and_M():
    a = group_key(_rec(N=[16], M=8))
    b = group_key(_rec(N=[64], M=99))
    assert a == b  # N and M are bound-absorbed -> same group


def test_speed_and_direction_derived_from_trafo():
    assert metric_name(group_key(_rec(trafo="trafo"))) == \
        "nfft/file/fast/forward/1d/init_guru"
    assert metric_name(group_key(_rec(trafo="adjoint_direct"))) == \
        "nfft/file/direct/adjoint/1d/init_guru"


def test_convert_takes_max_ratio_and_max_error():
    recs = [_rec(N=[16], accuracy=1e-14, bound=1e-13),   # ratio 0.1
            _rec(N=[64], accuracy=5e-13, bound=1e-12)]   # ratio 0.5, err 5e-13
    bmf = convert(recs)
    (name, measures), = bmf.items()
    assert measures["tightness-ratio"]["value"] == pytest.approx(0.5)
    assert measures["max-error"]["value"] == pytest.approx(5e-13)


def test_file_and_online_are_separate_metrics():
    bmf = convert([_rec(oracle="file"), _rec(oracle="online")])
    assert len(bmf) == 2


def test_nonpositive_bound_raises():
    with pytest.raises(ValueError, match="bound"):
        convert([_rec(bound=0.0)])


def test_read_ndjson_skips_blank_and_reports_bad_line():
    assert len(read_ndjson('{"module":"nfft"}\n\n')) == 1
    with pytest.raises(ValueError, match="line 1"):
        read_ndjson("not json\n")
