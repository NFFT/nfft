from tests.refgen import transforms as T


def test_type_i_freqs_unchanged():
    for N in ([4], [10], [10, 20], [15], [25]):
        assert T.freqs("nfft", N) == [list(range(-(n // 2), n - (n // 2))) for n in N]


def test_type_ii_axis_ascending_even():
    # ascending range -n/2+1 .. n/2; +n/2 at the last slot
    assert T.freqs("nfft", [4], variant=[T.TYPE_II])[0] == [-1, 0, 1, 2]


def test_type_ii_collapses_for_odd():
    assert T.freqs("nfft", [5], variant=[T.TYPE_II]) == T.freqs("nfft", [5])
