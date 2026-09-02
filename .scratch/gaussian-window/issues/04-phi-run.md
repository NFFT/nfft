# 04 — Give the Gaussian a `PHI_RUN`

Status: ready-for-agent

## Problem

The Gaussian has no `PHI_RUN` override, so it falls through to the generic loop
in `infft.h`: `2m+2` calls to `EXP` and `2m+2` divisions by `n` per run. That is
the path used by the no-precompute trafo and adjoint in 1-D, 2-D, 3-D and d-D,
and by `precompute_psi` and `precompute_full_psi`.

## Change

The Gaussian's neighbour ratio is itself a geometric sequence:

```
phi(nx0 - (l+1)) / phi(nx0 - l) = exp((2(nx0 - l) - 1)/b)
```

so stepping `l` multiplies the ratio by the per-axis constant `exp(-2/b)`.
Starting from the grid point `c` nearest the node — which also bounds the
running factor, the point of `02-centred-fast-gaussian.md` — the whole run is
two `EXP`s, one division and two multiplies per point in two branch-free loops:

```c
static inline void Y(gaussian_phi_run)(R *dst, R b_inv, R norm, R step,
  INT mi, R nx0)
{
  const INT last = 2 * mi + 1;
  INT c = (INT)LRINT(nx0);
  R r, p, g, w;
  INT l;

  if (c < 0) c = 0; else if (c > last) c = last;

  r = nx0 - (R)c;
  p = EXP(-(r * r) * b_inv) * norm;
  dst[c] = p;

  g = EXP((K(2.0) * r - K(1.0)) * b_inv);
  for (w = p, l = c + 1; l <= last; l++) { w *= g; g *= step; dst[l] = w; }

  g = step / EXP((K(2.0) * r - K(1.0)) * b_inv);
  for (w = p, l = c - 1; l >= 0; l--) { w *= g; g *= step; dst[l] = w; }
}
```

`step = EXP(-2/b)` joins `1/b` and `1/sqrt(pi b)` in `ths->b`. `c` is clamped so
a caller whose run sits differently against the support still gets one
evaluation per point, the way `kb_phi_run` derives its interior range rather
than assuming it.

## Measurements

Against a long double reference, worst over the run, relative to the run peak,
and time for one run (`-O3 -ffast-math -march=native`):

| m  | run err | direct err | run ns | direct ns |
|----|---------|------------|--------|-----------|
| 4  | 6.6e-17 | 1.5e-16    | 29.7   | 33.1      |
| 7  | 1.6e-16 | 1.9e-16    | 33.3   | 40.6      |
| 10 | 2.1e-16 | 2.4e-16    | 42.2   | 64.2      |
| 13 | 1.9e-16 | 2.8e-16    | 42.6   | 73.4      |
| 17 | 2.4e-16 | 2.1e-16    | 46.9   | 88.2      |

About 1 ulp — indistinguishable from evaluating every point directly — and the
error does not grow with `m`.

End to end, folded together with 01 and 03 (`N = 256`, `M = 20000`, sigma=2, no
precomputation, min of 60 runs, libraries built side by side and measured
interleaved):

| m  | baseline | prototype | speedup |
|----|----------|-----------|---------|
| 4  |  905 us  |  950 us   | 0.95x   |
| 7  | 1246 us  | 1104 us   | 1.13x   |
| 10 | 1784 us  | 1366 us   | 1.31x   |
| 13 | 2198 us  | 1497 us   | 1.47x   |

The FFT and the gather are inside those timings, so the window evaluation itself
gains more than the ratio shows. At `m = 4` the two fixed `EXP`s are a wash
against ten direct ones.

At `m = 13` the resulting no-precompute path (1497 us) is faster than `FG_PSI`
(1811 us in the same build) and more accurate (5.5e-15 against 6.6e-15) — worth
noting when deciding how much `FG_PSI` is still worth carrying.

## Tests

`tests/window.c`: the run must reproduce per-point `PHI` at every sub-cell
offset and at cut-offs past the default, with no non-finite entry. Shared with
`02-centred-fast-gaussian.md`.
