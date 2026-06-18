# Legacy Mathematica reference-data generators

These notebooks (`check_nfft`/`check_nfct`/`check_nfst`/`check_nfsft`, `.nb` + `.m`)
are the original Mathematica generators for the `tests/data/*.txt` reference data.
They are **superseded** by the Python generator in [`tests/refgen/`](../refgen/) (see
[ADR-0002](../../docs/adr/0002-python-reference-data-generator.md) and
[`docs/agents/test-methodology.md`](../../docs/agents/test-methodology.md)).

They are kept here only as provenance for the currently shipped data and are **not**
part of the build or `make check`. They are slated for removal in a later change once
the reference data is regenerated from `tests/refgen/`.
