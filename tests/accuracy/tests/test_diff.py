import json

import pytest

from tests.accuracy.diff import diff, digits, load_bmf_tree


def _bmf(**digit_by_name):
    return {
        n: {"accuracy-digits": {"value": v}, "max-error": {"value": 10**-v}}
        for n, v in digit_by_name.items()
    }


def test_digits_reads_accuracy_digits():
    assert digits(_bmf(a=13.5), "a") == 13.5


def test_improvement_and_regression_classified_and_pct():
    pr = {"t1": _bmf(a=14.3, b=12.0)}  # a improved +0.8, b regressed -1.0
    base = {"t1": _bmf(a=13.5, b=13.0)}
    r = diff(pr, base, gate=0.5)
    assert len(r.improvements) == 1 and len(r.regressions) == 1
    imp = r.improvements[0]
    assert (imp.testbed, imp.name) == ("t1", "a")
    assert imp.delta_digits == pytest.approx(0.8)
    assert imp.pct == pytest.approx(100 * 0.8 / 13.5)
    assert r.regressions[0].name == "b"
    assert r.unchanged_count == 0


def test_below_gate_is_unchanged():
    pr = {"t1": _bmf(a=13.9)}
    base = {"t1": _bmf(a=13.5)}  # +0.4 < 0.5 gate
    r = diff(pr, base, gate=0.5)
    assert r.improvements == [] and r.regressions == []
    assert r.unchanged_count == 1
    assert r.by_testbed["t1"] == {"improved": 0, "regressed": 0, "unchanged": 1}


def test_groups_sorted_by_abs_delta_desc():
    pr = {"t1": _bmf(a=15.0, b=16.0)}  # a +1.0, b +3.0
    base = {"t1": _bmf(a=14.0, b=13.0)}
    r = diff(pr, base, gate=0.5)
    assert [c.name for c in r.improvements] == ["b", "a"]


def test_added_and_removed_tracked_not_counted():
    pr = {"t1": _bmf(a=13.5, c=13.5)}
    base = {"t1": _bmf(a=13.5, d=13.5)}
    r = diff(pr, base, gate=0.5)
    assert r.added == ["t1::c"] and r.removed == ["t1::d"]
    assert r.improvements == [] and r.regressions == []


def test_load_bmf_tree_keys_by_testbed(tmp_path):
    (tmp_path / "ubuntu_gcc_kb_double.bmf.json").write_text(json.dumps(_bmf(a=13.5)))
    tree = load_bmf_tree(str(tmp_path))
    assert set(tree) == {"ubuntu_gcc_kb_double"}
    assert tree["ubuntu_gcc_kb_double"]["a"]["accuracy-digits"]["value"] == 13.5
