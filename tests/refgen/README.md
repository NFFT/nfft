# Reference-data generator

Computes the direct NDFT/NDCT/NDST at arbitrary precision (mpmath) and emits the
CUnit reference data + generated C registration headers + the data `Makefile.am`
list. Single source of truth — see `docs/agents/test-methodology.md`.

## Dependencies (dev only; not a build dependency)

Run via [uv](https://docs.astral.sh/uv/). The only dependency (`mpmath`, pinned for
byte-stable output) is declared **inline on the `uv run` command** — there is no
`pyproject.toml` or lockfile to maintain. uv resolves the deps into a cached,
ephemeral environment on each run.

## Regenerate everything

```bash
# from the repository root
uv run --with mpmath==1.3.0 python -m tests.refgen.generate --module all --precision 64
```

## Run the generator's own tests

```bash
uv run --with mpmath==1.3.0 --with pytest python -m pytest tests/refgen/tests -q
```

## Options

- `--module {nfft,nfct,nfst,all}`
- `--precision N`   working/printed significant decimal digits (default 64; ≥34 for quad)
- `--seed N`        PRNG seed (default 1)
- `--data-dir P`    output dir for *.txt (default tests/data)
- `--header-dir P`  output dir for generated headers (default tests/data/generated)
