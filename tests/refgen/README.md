# Reference-data generator

Computes the direct NDFT/NDCT/NDST at arbitrary precision (mpmath) and emits the
CUnit reference data + generated C registration headers + the data `Makefile.am`
list. Single source of truth — see `docs/agents/test-methodology.md`.

## Install (dev only; not a build dependency)

Managed by [uv](https://docs.astral.sh/uv/). From the repository root:

```bash
uv sync   # creates .venv from the locked deps (mpmath, pytest); byte-stable output
```

## Regenerate everything

```bash
# from the repository root
uv run python -m tests.refgen.generate --module all --precision 64
```

## Run the generator's own tests

```bash
uv run python -m pytest tests/refgen/tests -q
```

## Options

- `--module {nfft,nfct,nfst,all}`
- `--precision N`   working/printed significant decimal digits (default 64; ≥34 for quad)
- `--seed N`        PRNG seed (default 1)
- `--data-dir P`    output dir for *.txt (default tests/data)
- `--header-dir P`  output dir for generated headers (default tests/data/generated)
