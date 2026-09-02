# 06 — Revisit `m2K` and `WINDOW_HELP_ESTIMATE_m` for the Gaussian

Status: needs-triage

Do this last: both depend on `01-shape-parameter.md` having landed.

## `m2K`

`kernel/util/window.c` tabulates, for the Gaussian,

```c
static const INT m2K_[] = {0, 1, 3, 6, 7, 9, 11, 13, 15, 17, 19, 21, 22, 23, 24};
```

so `PRE_LIN_PSI` at the double default `m = 13` allocates
`K + 1 = 125829121` reals — 960 MB. Measured on `N = 256`, `M = 20000`:

| flag          | m=10   | m=13    |
|---------------|--------|---------|
| `PRE_PSI`     | 0.5 ms | 0.6 ms  |
| `PRE_LIN_PSI` | 9.6 ms | 22.4 ms |

34x slower than the exact table it is meant to approximate, and slightly less
accurate. Once the run is O(1) per point (`04-phi-run.md`), interpolating the
Gaussian buys nothing at all: the exact evaluation is already two `EXP`s per
run.

Either retabulate `m2K_` for the shifted `b`, or document that `PRE_LIN_PSI` is
not a useful flag for this window. Worth checking whether the Kaiser-Bessel
table has the same shape before concluding this is Gaussian-specific.

## `WINDOW_HELP_ESTIMATE_m`

Current: 13 (double), 5 (float), 17 (long double / quad).

After `01`:

- double `m = 13` reaches 5.5e-15 where it reached 3.3e-14; `m = 12` already
  beats today's `m = 13`, for ~8% fewer taps.
- float: `PHI_RUN` hits the roundoff floor (~4.6e-07) from `m = 7`, while the
  default `m = 5` gives 1.0e-06.

Re-derive all three from the measured curves rather than nudging them.
