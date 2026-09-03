# PRD: Gaussian window — accuracy and evaluation cost

Status: needs-triage

Companion to the Kaiser-Bessel work (#261, merged) and the B-spline / sinc-power
work (#263). Those two PRs restructured `PHI`, `PHI_HUT` and `PHI_RUN` so the
per-axis constants are hoisted into `ths->b` at `WINDOW_HELP_INIT` and the whole
run of `2m+2` window values is produced by one specialised routine. The Gaussian
never got that treatment: it still evaluates `POW(...,2.0)` and a `SQRT` per
point, still falls through to the generic `PHI_RUN`, and — unlike the other
windows — carries a second, older fast path (`FG_PSI` / `PRE_FG_PSI`) that
predates all of it.

This document records what the analysis found. There are five actionable items;
one of them is a latent correctness bug, and the largest accuracy win is a
one-line change that costs nothing.

All measurements below are from this repository at `000f1f6`, gcc 14,
`--with-window=gaussian`, double precision unless stated, `E_inf` measured
against `nfft_trafo_direct` and normalised by `sum |f_hat|`.

---

## 1. The shape parameter `b` is set for the wrong truncation width

`WINDOW_HELP_INIT` sets

```
b = 2 sigma / (2 sigma - 1) * m / pi
```

which is the classical balance (Potts–Steidl–Tasche) between the two error terms
of the truncated Gaussian:

- truncation, `~ exp(-m^2/b)`, amplified by the deconvolution factor
  `exp(pi^2 b / (4 sigma^2))` at the band edge;
- aliasing, `~ exp(-b pi^2 (1 - 1/sigma))`.

Balancing those exponents gives exactly the formula above — with `m` as the
half-width at which the window is truncated.

**But the window is not truncated at `m`.** `uo()` puts the run at
`u = floor(n x) - m` and the run covers `2m + 2` points, so the largest
`|n x - l|` a run actually reaches is between `m` and `m + 1`, never `m`. The
stencil is one cell wider than the parameter assumes, and `b` is correspondingly
too small: the window is being pinched harder than the truncation requires.

Substituting `m -> m + 1`:

```
b = 2 sigma / (2 sigma - 1) * (m + 1) / pi
```

### Evidence

An independent mpmath model (50 digits, exact single-frequency error
`E(k,x)` for the same `2m+2` stencil, worst over `k` and a sub-cell sweep of `x`)
scanned `b` for its true minimiser. The minimiser tracks
`2 sigma/(2 sigma - 1) * (m + c)/pi` with `c` between 0.9 and 1.2 across every
`(sigma, m)` tested, and `c = 1` is within a few percent of optimal everywhere:

| sigma | m  | E(b current) | E(b optimal) | ratio |
|-------|----|--------------|--------------|-------|
| 1.25  | 8  | 8.14e-05     | 3.75e-05     |  2.2  |
| 1.50  | 8  | 1.59e-06     | 2.45e-07     |  6.5  |
| 2.00  | 8  | 3.13e-08     | 2.69e-09     | 11.7  |
| 2.00  | 13 | 6.41e-13     | 9.76e-14     |  6.6  |
| 3.00  | 13 | 3.88e-15     | 2.78e-16     | 14.0  |

Confirmed in the library itself (`N = 256`, `PRE_PSI`, deterministic sub-cell
sweep, so the figure is a worst case rather than a random draw):

| case                | `m/pi` (today) | `(m+1)/pi` | gain |
|---------------------|----------------|------------|------|
| sigma=2,    m=4     | 7.83e-06       | 1.14e-06   | 6.9x |
| sigma=2,    m=7     | 1.10e-08       | 1.43e-09   | 7.7x |
| sigma=2,    m=10    | 1.70e-11       | 2.79e-12   | 6.1x |
| sigma=2,    m=13    | 2.69e-14       | 4.53e-15   | 5.9x |
| sigma=1.25, m=8     | 3.64e-06       | 1.57e-06   | 2.3x |
| sigma=1.5,  m=8     | 6.57e-08       | 1.38e-08   | 4.8x |
| sigma=3,    m=8     | 5.73e-11       | 6.76e-12   | 8.5x |

A scan over the shift constant (0.0, 0.8, 1.0, 1.2, 1.4) puts the library
optimum at 1.0–1.2 in every configuration, so `m + 1` is both the best fit and
the one with a structural reason behind it.

Equivalently, at `sigma = 2` a user who wants today's `m = 13` accuracy can run
`m = 12`: about 8% fewer taps in the B step.

**Cost: none.** One line in `WINDOW_HELP_INIT`.

Issue: [`issues/01-shape-parameter.md`](issues/01-shape-parameter.md)

---

## 2. `FG_PSI` / `PRE_FG_PSI` overflow — NaN output in single precision

Fast Gaussian gridding factorises the run as

```
phi(nx0 - l) = phi(nx0) * q^l * e(l),   q = exp(2 nx0 / b),  e(l) = exp(-l^2/b)
```

anchored at `l = 0`, i.e. at the **far left end** of the run, where `phi` is at
its smallest. The running factor `q^l` therefore has to climb the entire
dynamic range of the window and back down again:

| sigma | m  | b       | q         | q^(2m+1)  | min e(l)  |
|-------|----|---------|-----------|-----------|-----------|
| 2.0   | 4  | 1.6977  | 3.61e+02  | 1.04e+23  | 1.90e-21  |
| 2.0   | 8  | 3.3953  | 2.01e+02  | 1.37e+39  | 1.08e-37  |
| 2.0   | 13 | 5.5174  | 1.60e+02  | 3.19e+59  | 4.15e-58  |
| 2.0   | 24 | 10.1859 | 1.35e+02  | 2.86e+104 | 4.26e-103 |

Single precision tops out at 3.4e38. **From `m = 8` the running factor
overflows to `+inf` and the companion table underflows to zero**; the product is
NaN.

Verified end to end against `libnfft3f` (float build, Gaussian, sigma=2,
N=256, M=2000):

```
  PHI_RUN      m= 7  E_inf = 4.6166e-07  non-finite outputs: 0/2000
  FG_PSI       m= 7  E_inf = 7.4587e-07  non-finite outputs: 0/2000
  PHI_RUN      m= 9  E_inf = 4.7523e-07  non-finite outputs: 0/2000
  FG_PSI       m= 9                      non-finite outputs: 2000/2000  <--
  PRE_FG_PSI   m= 9                      non-finite outputs: 2000/2000  <--
  FG_PSI       m=13                      non-finite outputs: 2000/2000  <--
```

The default float cut-off is `WINDOW_HELP_ESTIMATE_m = 5`, which is why CI has
never caught this: the tests only ever instantiate the default `m`. Any float
user who raises `m` past 7 gets NaN, silently.

Even inside its range the anchoring costs accuracy. Window values against a long
double reference, worst over the run, relative to the run peak:

| precision | m  | FG as implemented | centred | direct |
|-----------|----|-------------------|---------|--------|
| float     | 7  | 1.2e-06 (~10 ulp) | 5.6e-08 | 5.6e-08 |
| double    | 8  | 3.3e-15 (~15 ulp) | 1.6e-16 | 1.6e-16 |
| double    | 13 | 1.9e-14 (~85 ulp) | 8.3e-17 | 1.2e-16 |

At `m = 13` that 1.9e-14 is the same order as the whole transform error, so
`FG_PSI` is already the accuracy-limiting path — visible in `make check` today,
where the `nfst_online` 2-D case reports `PRE FG PSI` at 5.04e-14 against
`PRE PSI` at 1.45e-14 on identical inputs. Once item 1 lands and the target
drops by 6x, it dominates outright.

**Fix:** anchor the recurrence at the grid point nearest the node instead of at
the end of the run. With `c` that index and `r = nx0 - c`, `|r| <= 1/2`, the
running factor never leaves `exp(+-(m+1)/b)` — about 12 to 19, bounded in every
precision at every `m` — and it is stepped at most `m+1` times from the peak of
the run rather than `2m+1` times from its smallest value.

Issue: [`issues/02-centred-fast-gaussian.md`](issues/02-centred-fast-gaussian.md)

---

## 3. `PHI` costs a `POW`, a `SQRT` and two divisions per point

```c
#define PHI(n,x,d) ((R)EXP(-POW((x)*((R)n),K(2.0)) / ths->b[d])/SQRT(KPI*ths->b[d]))
#define PHI_HUT(n,k,d) ((R)EXP(-(POW(KPI*(k)/n,K(2.0))*ths->b[d])))
```

`SQRT(KPI * ths->b[d])` is a loop-invariant recomputed at every point, `ths->b[d]`
is a memory operand the compiler cannot always hoist across a store to `dst`, and
`POW(y, 2.0)` is a library call unless `-ffast-math` happens to fold it. This is
exactly what #261 and #263 removed from the other three windows.

Store `1/b` and `1/sqrt(pi b)` per axis in `ths->b` and take `nx` rather than `x`,
so `PHI` becomes one multiply, one `EXP` and one multiply.

Issue: [`issues/03-hoist-constants.md`](issues/03-hoist-constants.md)

---

## 4. The Gaussian has no `PHI_RUN`

It falls through to the generic loop in `infft.h`: `2m+2` calls to `EXP` and
`2m+2` divisions by `n`. This is the path taken by the no-precompute
trafo/adjoint in 1-D, 2-D, 3-D and d-D, and by `precompute_psi` and
`precompute_full_psi`.

A centred run costs **two `EXP`s and one division for the whole run**,
independent of `m`, plus two multiplies per point in two branch-free loops. Two
formulations were measured against a long double reference:

| m  | ratio-recurrence err | direct err | ratio ns/run | direct ns/run |
|----|----------------------|------------|--------------|---------------|
| 4  | 6.6e-17              | 1.5e-16    | 29.7         | 33.1          |
| 7  | 1.6e-16              | 1.9e-16    | 33.3         | 40.6          |
| 10 | 2.1e-16              | 2.4e-16    | 42.2         | 64.2          |
| 13 | 1.9e-16              | 2.8e-16    | 42.6         | 73.4          |
| 17 | 2.4e-16              | 2.1e-16    | 46.9         | 88.2          |

About 1 ulp — indistinguishable from evaluating every point directly — at up to
1.9x the speed, with the gap widening in `m`.

End to end (`N = 256`, `M = 20000`, sigma=2, no precomputation, min of 60 runs,
baseline and prototype libraries built side by side and measured interleaved):

| m  | baseline    | prototype   | speedup | baseline E_inf | prototype E_inf |
|----|-------------|-------------|---------|----------------|-----------------|
| 4  |    905 us   |    950 us   | 0.95x   | 9.80e-06       | 1.14e-06        |
| 7  |   1246 us   |   1104 us   | 1.13x   | 1.37e-08       | 1.78e-09        |
| 10 |   1784 us   |   1366 us   | 1.31x   | 2.08e-11       | 3.56e-12        |
| 13 |   2198 us   |   1497 us   | 1.47x   | 3.33e-14       | 5.50e-15        |

(The prototype folds items 1, 3 and 4 together; the FFT and the gather are in
these timings, so the window evaluation itself speeds up by more than the ratio
shows. At `m = 4` the run is short enough that the two fixed `EXP`s are a wash.)

Notably, at `m = 13` the prototype's no-precompute path (1497 us) is **faster
than `FG_PSI`** (1811 us in the same build) and more accurate. Whatever
`FG_PSI` was for, a proper `PHI_RUN` supersedes it; `PRE_FG_PSI` keeps its own
reason to exist (two reals per node instead of `2m+2`).

Issue: [`issues/04-phi-run.md`](issues/04-phi-run.md)

---

## 5. Smaller items

**`fg_exp_l` is built by a nested recurrence.** `nfft_1d_init_fg_exp_l` and its
copies in `nfft.c`, `nfct.c` and `nfst.c` build `exp(-l^2/b)` as
`fg_exp_l[l] = fg_exp_l[l-1] * exp(-1/b)^(2l-1)`, accumulating roughly `2l` ulp
over the table. It is built once per transform, `2m+2` entries — a direct
`EXP(-l*l/b)` is free at that scale and exact to 1 ulp.

**`m2K` for the Gaussian.** `PRE_LIN_PSI` at the default `m = 13` allocates
`K + 1 = 125829121` reals — **960 MB** — and measures 22 ms/transform against
0.6 ms for `PRE_PSI` on the same problem, i.e. 34x slower and less accurate.
Once the run is O(1) per point, linear interpolation has nothing left to buy for
this window. Either retabulate `m2K_` for the Gaussian or document that
`PRE_LIN_PSI` is not the right flag here.

**`WINDOW_HELP_ESTIMATE_m`.** After item 1 the double default `m = 13` reaches
5.5e-15 where it used to reach 3.3e-14; `m = 12` already beats today's `m = 13`.
In float, `PHI_RUN` hits the roundoff floor (~4.6e-07) at `m = 7` while the
default `m = 5` gives 1.0e-06. Both deserve a re-look once the rest lands.

Issue: [`issues/05-fg-exp-table.md`](issues/05-fg-exp-table.md),
[`issues/06-m2k-and-default-m.md`](issues/06-m2k-and-default-m.md)

---

## Sequencing

1. `03-hoist-constants` + `04-phi-run` — pure win, no semantic change, matches
   the shape of #261/#263. Land first; they are one commit's worth of work.
2. `02-centred-fast-gaussian` — fixes the float NaN. Independent of the others.
3. `01-shape-parameter` — semantic. Changes every Gaussian result, so it wants
   its own commit, a note in the accuracy report, and a look at whether any
   reference data in `tests/data/` is Gaussian-specific.
   **Done:** branch `feature/gaussian-shape-parameter` (`d5b233b`), off
   `develop`, green in all three precisions with no test-data changes.
4. `05`, `06` — cleanup, after the above.

Items 1, 3 and 4 have been prototyped together and pass the full CUnit suite
(2355 assertions, `make check` green) with no test-data changes. The prototype
diff is in [`prototype.patch`](prototype.patch); it is a starting point, not a
finished change — it has no unit tests of its own and does not touch the
`FG_PSI` recurrence (item 2) or the `fg_exp_l` tables (item 5).

Per #263, new window code is expected to come with `tests/window.c` checks
against an mpmath reference under `tests/windowref/`. Each issue below says what
that means for it.
