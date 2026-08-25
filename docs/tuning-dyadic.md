# The dyadic ladder against the legacy choice

Three arms over the same nodes, the same data and the same `plan_ng`
path, so the only difference is the parameter pair.

- **legacy** -- `n = 2*next_power_of_2(N)`, cut-off found by searching
  upward for the smallest `m` whose *measured* error meets the goal.
  The legacy API has no cut-off estimator, so this is an oracle it
  does not possess. It is also rung 1 of the ladder.
- **dyadic** -- `nfft_tune_plan_dyadic`, the three rungs
  `n = 2^j * next_power_of_2(N)`.

Speedups are against legacy, so above 1.00x is a win. The tuner runs
unaided: no measurement, no refinement.

The N set spans `t = next_power_of_2(N)/N`, which is what decides
whether the ladder has a rung below the legacy grid at all -- rung 0
needs `t >= 5/4`. Below that the ladder can only return the legacy
grid or a wider one.

| N | 140 | 160 | 200 | 243 | 250 | 251 | 255 | 256 | 320 | 512 | 600 | 1100 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| t | 1.83 | 1.60 | 1.28 | 1.05 | 1.02 | 1.02 | 1.00 | 1.00 | 1.60 | 1.00 | 1.71 | 1.86 |

## Overall

| arm | cases | median | p10 | min | faster or tied | goal met |
|---|---|---|---|---|---|---|
| dyadic | 576 | 1.01x | 0.88x | 0.71x | 349/576 | 576/576 |

## By shape

### N-dominated (192 cases)

| arm | cases | median | p10 | min | faster or tied | goal met |
|---|---|---|---|---|---|---|
| dyadic | 192 | 1.18x | 0.96x | 0.80x | 127/192 | 192/192 |

### balanced (192 cases)

| arm | cases | median | p10 | min | faster or tied | goal met |
|---|---|---|---|---|---|---|
| dyadic | 192 | 1.04x | 0.89x | 0.74x | 122/192 | 192/192 |

### M-dominated (192 cases)

| arm | cases | median | p10 | min | faster or tied | goal met |
|---|---|---|---|---|---|---|
| dyadic | 192 | 1.00x | 0.85x | 0.71x | 100/192 | 192/192 |

## By oversampling headroom

`t < 5/4` is the control: rung 0 is illegal, so the ladder returns the
legacy grid or a wider one and can only tie or win on the cut-off.
A control row below 1.00x is a cut-off regression, not a ladder one.

### t < 5/4 (control) (288 cases)

| arm | cases | median | p10 | min | faster or tied | goal met |
|---|---|---|---|---|---|---|
| dyadic | 288 | 0.98x | 0.85x | 0.71x | 99/288 | 288/288 |

### t >= 5/4 (rung 0 legal) (288 cases)

| arm | cases | median | p10 | min | faster or tied | goal met |
|---|---|---|---|---|---|---|
| dyadic | 288 | 1.42x | 0.97x | 0.74x | 250/288 | 288/288 |

## By precision

### float (144 cases)

| arm | cases | median | p10 | min | faster or tied | goal met |
|---|---|---|---|---|---|---|
| dyadic | 144 | 1.00x | 0.87x | 0.71x | 81/144 | 144/144 |

### double (216 cases)

| arm | cases | median | p10 | min | faster or tied | goal met |
|---|---|---|---|---|---|---|
| dyadic | 216 | 1.01x | 0.89x | 0.74x | 139/216 | 216/216 |

### long double (216 cases)

| arm | cases | median | p10 | min | faster or tied | goal met |
|---|---|---|---|---|---|---|
| dyadic | 216 | 1.02x | 0.89x | 0.79x | 129/216 | 216/216 |

## Which rung the ladder took

Rung 1 is the legacy grid. Rung 0 halves it, rung 2 doubles it.

| shape | rung 0 | rung 1 (legacy) | rung 2 | median speedup |
|---|---|---|---|---|
| N-dominated | 94 | 98 | 0 | 1.18x |
| balanced | 94 | 98 | 0 | 1.04x |
| M-dominated | 81 | 111 | 0 | 1.00x |

## How much of this is timing noise

In 147 of 576 cases the ladder returns the legacy pair itself, so the
two arms run identical parameters and their measured ratio is pure
noise. It should read 1.00x exactly:

| median | mean | p05 | p95 | min | max | within 2 % | within 5 % |
|---|---|---|---|---|---|---|---|
| 1.000 | 1.002 | 0.933 | 1.059 | 0.800 | 1.177 | 102/147 | 126/147 |

The median is exact and the mean is unbiased, so **the medians in
this document carry signal**. The per-case spread does not: a
single row at 0.80x or 1.20x is within what identical parameters
produce, so the min and max columns above, and the counts of rows
below 1.00x, are noise-dominated and should not be read as losses.

## Where the ladder loses

227 of 576 cases below 1.00x.

- 200 are on the legacy grid itself, of which 143 carry a larger
  cut-off than the oracle's. That is the envelope's cost, not the
  ladder's: the same cut-off gap the model has always had.
- 27 chose a different rung.

| prec | N | t | M | shape | goal | dir | rung | n | m | n_leg | m_leg | speedup |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| float | 256 | 1.00 | 1024 | M-dominated | 1e-04 | adjoint | 1 | 512 | 3 | 512 | 2 | 0.71x |
| double | 200 | 1.28 | 800 | M-dominated | 1e-08 | forward | 0 | 256 | 7 | 512 | 4 | 0.74x |
| float | 243 | 1.05 | 243 | balanced | 1e-02 | adjoint | 1 | 512 | 2 | 512 | 1 | 0.74x |
| double | 255 | 1.00 | 1020 | M-dominated | 1e-04 | adjoint | 1 | 512 | 3 | 512 | 2 | 0.78x |
| float | 255 | 1.00 | 1020 | M-dominated | 1e-04 | adjoint | 1 | 512 | 3 | 512 | 2 | 0.78x |
| float | 251 | 1.02 | 251 | balanced | 1e-02 | adjoint | 1 | 512 | 2 | 512 | 1 | 0.79x |
| long double | 250 | 1.02 | 250 | balanced | 1e-20 | forward | 1 | 512 | 11 | 512 | 10 | 0.79x |
| long double | 512 | 1.00 | 512 | balanced | 1e-08 | adjoint | 1 | 1024 | 5 | 1024 | 4 | 0.80x |
| long double | 512 | 1.00 | 128 | N-dominated | 1e-14 | forward | 1 | 1024 | 8 | 1024 | 7 | 0.80x |
| long double | 512 | 1.00 | 128 | N-dominated | 1e-20 | forward | 1 | 1024 | 11 | 1024 | 11 | 0.80x |
| long double | 255 | 1.00 | 255 | balanced | 1e-20 | forward | 1 | 512 | 11 | 512 | 10 | 0.80x |
| long double | 251 | 1.02 | 251 | balanced | 1e-20 | forward | 1 | 512 | 11 | 512 | 10 | 0.81x |
| float | 243 | 1.05 | 972 | M-dominated | 1e-02 | forward | 1 | 512 | 2 | 512 | 1 | 0.81x |
| float | 243 | 1.05 | 972 | M-dominated | 1e-04 | adjoint | 1 | 512 | 3 | 512 | 2 | 0.81x |
| long double | 250 | 1.02 | 250 | balanced | 1e-08 | forward | 1 | 512 | 5 | 512 | 4 | 0.81x |
| double | 512 | 1.00 | 2048 | M-dominated | 1e-08 | adjoint | 1 | 1024 | 5 | 1024 | 4 | 0.82x |
| float | 256 | 1.00 | 1024 | M-dominated | 1e-02 | adjoint | 1 | 512 | 2 | 512 | 1 | 0.82x |
| float | 512 | 1.00 | 2048 | M-dominated | 1e-04 | adjoint | 1 | 1024 | 3 | 1024 | 2 | 0.82x |
| double | 250 | 1.02 | 1000 | M-dominated | 1e-04 | adjoint | 1 | 512 | 3 | 512 | 2 | 0.82x |
| float | 251 | 1.02 | 1004 | M-dominated | 1e-04 | adjoint | 1 | 512 | 3 | 512 | 2 | 0.83x |

## Accuracy

| arm | goal met |
|---|---|
| legacy (oracle) | 576/576 |
| dyadic | 576/576 |

The input is drawn with real and imaginary parts on `[0, 1)`, matching
`nfft_vrand_unit_complex` and the sweep the error model is fitted to.
Centred data measures a forward error up to 2.6x smaller; see issue 05
in `.scratch/tune-dyadic/`.

