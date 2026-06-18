"""The (d, N, M) test-case grids per module. Single source of truth for which
reference files exist. Edit here to add/remove cases, then regenerate."""

# Each entry: (d, N-as-list, M)
GRIDS = {
    "nfft": (
        [(1, [n], m) for n in (1, 2, 4, 10, 20, 50) for m in (1, 10, 20, 50)]
        + [(2, list(N), m) for N in ((10, 10), (10, 20), (20, 10), (20, 20)) for m in (20, 50)]
        + [(3, [10, 10, 10], 10)]
    ),
    "nfct": (
        [(1, [n], m) for n in (1, 2, 4, 10, 25, 50) for m in (1, 10, 25, 50)]
        + [(2, list(N), m) for N in ((10, 10), (10, 25), (25, 10), (25, 25)) for m in (25, 50)]
        + [(3, [10, 10, 10], 10)]
    ),
    "nfst": (
        [(1, [n], m) for n in (2, 4, 10, 25, 50) for m in (1, 10, 25, 50)]
        + [(2, list(N), m) for N in ((10, 10), (10, 25), (25, 10), (25, 25)) for m in (25, 50)]
        + [(3, [10, 10, 10], 10)]
    ),
}

KINDS = ("trafo", "adjoint")  # both forward and adjoint files for every grid entry


def basename(module, kind, d, N, M):
    adj = "_adjoint" if kind == "adjoint" else ""
    return "%s%s_%dd_%s_%d" % (module, adj, d, "_".join(str(n) for n in N), M)
