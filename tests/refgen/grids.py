"""The (d, N, M, variant) test-case grids per module. Single source of truth for
which reference files exist. Edit here to add/remove cases, then regenerate."""

from tests.refgen.transforms import TYPE_I, TYPE_II

# Each entry: (d, N-as-list, M, variant). variant is None (legacy type-I, all
# axes) except for the nfft type-II entries, where it is a list of per-axis
# TYPE_I/TYPE_II tags aligned with N.

# The original nfft grid: all type-I, even N.
# This exact list is what the legacy header/registration (tests/nfft.c) must
# keep seeing — kept as a separate name so registration.py can filter on it
# precisely, rather than re-deriving "legacy" from N's parity (N=1 is odd and
# IS part of the legacy grid, so parity alone is not a valid discriminator).
NFFT_LEGACY_GRID = (
    [(1, [n], m, None) for n in (1, 2, 4, 10, 20, 50) for m in (1, 10, 20, 50)]
    + [(2, list(N), m, None) for N in ((10, 10), (10, 20), (20, 10), (20, 20)) for m in (20, 50)]
    + [(3, [10, 10, 10], 10, None)]
)

# Native-only additions: odd-N type-I, and even-N type-II. Not part
# of the legacy (tests/nfft.c) roster; consumed via the native header instead.
_NFFT_NATIVE_ONLY = (
    # odd-N type-I (native-only; legacy rejects odd N)
    [(1, [n], m, None) for n in (5, 15, 25) for m in (10, 50)]
    + [(2, list(N), 50, None) for N in ((5, 10), (10, 5), (5, 5))]
    # type-II even-N (native-only)
    + [(1, [n], m, [TYPE_II]) for n in (4, 10, 20) for m in (10, 50)]
    + [(2, list(N), 50, [TYPE_II, TYPE_I]) for N in ((10, 20), (20, 10))]
)

# Native-only: complete the type/parity coverage in 2D-4D, so every dimension
# has a case that is type-I even throughout (2D and 3D get theirs from the
# legacy grid), one type-II throughout, one odd throughout, and -- from d = 3,
# where there are enough axes -- one mixing all three. Sizes stay on the
# per-dimension scale of the existing grid so mpmath stays cheap.
_NFFT_NATIVE_TYPES = [
    (2, [10, 20], 50, [TYPE_II, TYPE_II]),
    (2, [5, 10], 50, [TYPE_I, TYPE_II]),  # odd + even-II
    (3, [4, 10, 20], 10, [TYPE_II, TYPE_II, TYPE_II]),
    (3, [5, 9, 15], 10, None),  # odd everywhere
    (3, [10, 5, 20], 10, [TYPE_I, TYPE_I, TYPE_II]),  # even-I + odd + even-II
    (4, [4, 6, 8, 10], 10, None),  # even, type-I
    (4, [4, 6, 8, 10], 10, [TYPE_II] * 4),  # even, type-II
    (4, [3, 5, 7, 9], 10, None),  # odd
    (4, [4, 5, 6, 7], 10, [TYPE_I, TYPE_I, TYPE_II, TYPE_I]),  # three-way mix
]

# Native-only: unit axes, which the new API elides from the problem geometry.
# One unit axis, several, and the all-unit case that compresses to rank 0, per
# dimension; 1D is already covered by the legacy N=1 files.
_NFFT_NATIVE_UNIT = [
    (2, [1, 10], 50, None),
    (2, [10, 1], 50, None),
    (2, [1, 1], 20, None),  # rank 0
    (3, [1, 10, 10], 10, None),
    (3, [10, 1, 10], 10, None),
    (3, [10, 10, 1], 10, None),
    (3, [1, 1, 10], 10, None),
    (3, [1, 1, 1], 10, None),  # rank 0
    (4, [1, 4, 6, 8], 10, None),
    (4, [4, 1, 6, 1], 10, None),
    (4, [4, 6, 8, 1], 10, None),
    # live axes keep their own variants across elision
    (4, [1, 10, 1, 20], 10, [TYPE_I, TYPE_II, TYPE_I, TYPE_II]),
    (4, [1, 1, 1, 1], 10, None),  # rank 0
]

GRIDS = {
    "nfft": NFFT_LEGACY_GRID + _NFFT_NATIVE_ONLY + _NFFT_NATIVE_TYPES + _NFFT_NATIVE_UNIT,
    "nfct": (
        [(1, [n], m, None) for n in (1, 2, 4, 10, 25, 50) for m in (1, 10, 25, 50)]
        + [(2, list(N), m, None) for N in ((10, 10), (10, 25), (25, 10), (25, 25)) for m in (25, 50)]
        + [(3, [10, 10, 10], 10, None)]
    ),
    "nfst": (
        [(1, [n], m, None) for n in (2, 4, 10, 25, 50) for m in (1, 10, 25, 50)]
        + [(2, list(N), m, None) for N in ((10, 10), (10, 25), (25, 10), (25, 25)) for m in (25, 50)]
        + [(3, [10, 10, 10], 10, None)]
    ),
}

KINDS = ("trafo", "adjoint")  # both forward and adjoint files for every grid entry


def basename(module, kind, d, N, M, variant=None):
    adj = "_adjoint" if kind == "adjoint" else ""
    tag = "_t2" + "".join(str(v) for v in variant) if (variant and any(v == TYPE_II for v in variant)) else ""
    return "%s%s_%dd_%s_%d%s" % (module, adj, d, "_".join(str(n) for n in N), M, tag)
