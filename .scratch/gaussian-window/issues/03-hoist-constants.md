# 03 — Hoist the Gaussian's per-axis constants out of `PHI` and `PHI_HUT`

Status: ready-for-agent

## Problem

```c
#define PHI(n,x,d) ((R)EXP(-POW((x)*((R)n),K(2.0)) / ths->b[d])/SQRT(KPI*ths->b[d]))
#define PHI_HUT(n,k,d) ((R)EXP(-(POW(KPI*(k)/n,K(2.0))*ths->b[d])))
```

Per evaluation: a `POW(y, 2.0)`, a division by `b`, and a `SQRT` of a
loop-invariant. `ths->b[d]` is a memory operand the compiler cannot reliably
hoist across a store to the destination array. This is what #261 and #263
removed from the other three windows.

## Change

Follow the established layout: widen `ths->b` and precompute per axis in
`WINDOW_HELP_INIT`.

```c
  #define GA_B(ax)     (ths->b[(ax)])
  #define GA_B_INV(ax) (ths->b[(ths->d) + (ax)])
  #define GA_NORM(ax)  (ths->b[2 * (ths->d) + (ax)])
  #define PHI_HUT(n,k,ax) (Y(gaussian_phi_hut)(GA_B(ax), (R)(k) * KPI / (R)(n)))
  #define PHI(n,x,ax)     (Y(gaussian_phi)(GA_B_INV(ax), GA_NORM(ax), (R)(n) * (x)))
```

with

```c
static inline R Y(gaussian_phi_hut)(R b, R t) { return EXP(-(t * t) * b); }
static inline R Y(gaussian_phi)(R b_inv, R norm, R nx)
{ return EXP(-(nx * nx) * b_inv) * norm; }
```

`PHI` then costs one multiply, one `EXP` and one multiply. Taking `nx` rather
than `x` also matches what every caller already has to hand, and is what
`PHI_RUN` needs (see `04-phi-run.md`).

`ths->b[ax]` must stay the raw `b` at index `ax`: the `FG_PSI` blocks in
`nfft.c`, `nfct.c` and `nfst.c` read it directly, and `PHI_HUT` still wants it.

## Notes

Prototyped together with 01 and 04; `make check` green, 2355 assertions. See
`../prototype.patch`.
