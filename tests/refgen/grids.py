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

GRIDS = {
    "nfft": NFFT_LEGACY_GRID + _NFFT_NATIVE_ONLY,
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
