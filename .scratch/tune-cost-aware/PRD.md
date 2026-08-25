# Cost-aware NFFT parameter tuning — handoff

Written 2026-08-25. This document is self-contained: it carries everything a
fresh session needs, with no reliance on the conversation that produced it.

## 1. What already exists

Branch **`feature/tune-m`**, pushed to `origin`, branched from `feature/planner`
at `2a9e9054`. Two commits:

| commit | what |
|---|---|
| `a0107d95` | the three tuning entry points, tests, build wiring, the study |
| `9ef638e5` | the head-to-head against the legacy parameter choice |

### The API

Declared in `include/nfft3.h` inside `NFFT_DEFINE_PLANNER_API`, so the names
mangle per precision (`nfft_tune`, `nfftf_tune`, `nfftl_tune`, …). Implemented
in `kernel/nfft/tune.c`. All three are **1-D and Kaiser-Bessel only**, analytic
(no measurement), and cost microseconds.

```c
/* Smallest cut-off m reaching `goal` at a given geometry.
 * 1 = met, 0 = not met (*m is the best available), -1 = bad args. */
int X(tune)(NFFT_INT N, NFFT_INT n, int adjoint, R goal, int *m, R *attained);

/* Smallest oversampled size n at which some m reaches `goal`.
 * 1 = reachable, 0 = not even at sigma=4, -1 = bad args. */
int X(tune_sigma)(NFFT_INT N, int adjoint, R goal, NFFT_INT *n, R *attained);

/* Both at once: cap the goal at what the window can deliver, find the
 * smallest sufficient oversampling, round n up to the next even 5-smooth
 * size, re-derive and re-check m there.
 * 1 = goal met, 0 = goal was below the reachable floor and the capped goal
 * was met instead, -1 = bad args. */
int X(tune_plan)(NFFT_INT N, int adjoint, R goal,
                 NFFT_INT *n, int *m, R *attained);
```

`adjoint` picks the error measure — 0 for `max_j |f_j - s_j| / ||f_hat||_1`,
non-zero for `max_k |fhat_k - s_k| / ||f||_1`. These are the measures
`tests/nfft.c` uses.

### The error model

In `kernel/nfft/tune.c`. With `u = 1 - 1/sigma`, `sigma = n/N`:

```
b = 2*pi*(1 - 1/(2*sigma))     Kaiser-Bessel shape parameter
D = 2*pi*sqrt(u)               window truncation decays as exp(-D*m)
A = b - D                      deconvolution amplifies roundoff as exp(A*m)

T(m,s)   = a * u^p * m^r * exp(-D*m)
U(m,s,e) = c * e * u^q * exp(alpha*A*m)          e = machine epsilon
E        = T + U
```

`D` and `A` are exact, not fitted — they follow from the identity
`b^2 - (pi/sigma)^2 = D^2`. Only the prefactors are fitted:

| | a | p | r | c | alpha | q |
|---|---|---|---|---|---|---|
| forward | 3.0431 | 0.902205 | -0.0183106 | 68.9787 | 0.967705 | 1.67263 |
| adjoint | 2.03698 | 0.234342 | 0.401585 | 33.5633 | 0.994013 | 1.15639 |

Facts worth knowing before touching this:

- **`A > 0` at every finite sigma**, so the error always falls with `m`, bottoms
  out, and climbs again. There is no `sigma` at which arbitrary accuracy is
  reachable in floating point.
- The floor is near `eps^(D/b)` where `A` is large; by `sigma = 2` the
  amplification is weak enough that it is back to a small multiple of `eps` in
  every precision. Measured floors, worst over N, forward: at `sigma = 3/2`,
  32·eps float / 21·eps double / **6394·eps long double** — the widest mantissa
  degrades first as sigma drops. At `sigma = 2`: 16 / 0.9 / 4.2 eps.
- **`sigma` is floored at 5/4**, checked in exact integer arithmetic
  (`4*n >= 5*N`); below that the fit is out of range and the attainable
  accuracy collapses.
- `TUNE_M_MAX` is 40, the top of the fitted range.
- The fitted constants are an **upper envelope**: they were raised until the
  formula dominates every one of 4674 measured points. So the model
  systematically over-estimates error, which is why it hands out roughly one
  more `m` than a measured search would. See §4.

### The study

`.scratch/sigma-m-study/` holds the drivers, data and fits. `README.md` there
is the entry point. Data is gzipped CSV under `data/`. Regenerate with:

```sh
sh .scratch/sigma-m-study/run_sweeps.sh "$PWD" /tmp/sweeps 40
cd .scratch/sigma-m-study && python3 constants.py /tmp/sweeps/sweep-*.csv
```

`constants.py` prints the constants above plus the validation numbers.
`selectform.py` picks the model form by validation rather than in-sample fit —
use it if you change the functional form, because a form that fits better
in-sample can pick worse `m` (an `m^w` term in the roundoff branch fit well and
produced `+8` to `+15` outliers in the chosen `m`).

## 2. How to build and test

Three CMake trees, one per precision. `NFFT_ENABLE_FLOAT` and
`NFFT_ENABLE_LONG_DOUBLE` are mutually exclusive; omitting both gives double.

```sh
cmake -S . -B bx-d                             -DNFFT_ENABLE_TESTS=ON \
      -DNFFT_ENABLE_OPENMP=OFF -DNFFT_ENABLE_EXAMPLES=OFF \
      -DNFFT_ENABLE_APPLICATIONS=OFF -DCMAKE_BUILD_TYPE=Release
cmake -S . -B bx-f -DNFFT_ENABLE_FLOAT=ON       ...same...
cmake -S . -B bx-l -DNFFT_ENABLE_LONG_DOUBLE=ON ...same...

sh .scratch/sigma-m-study/run_tests.sh "$PWD" bx-d   # and bx-f, bx-l
```

`run_tests.sh` builds `checkall_ng`, runs it, and reports the CUnit failure
count. Expect `exit=0 cunit_failures=0` in all three.

The tuning tests are the `tune` suite in `tests/tune_ng.c`, registered in
`tests/check_ng.c`: `meets_goal`, `unreachable`, `geometries`, `bad_args`,
`sigma_agrees`, `sigma_limits`, `plan`, `plan_capped`. `meets_goal` and `plan`
run real transforms and compare against the planner's own direct NDFT.

Autotools must keep working too — but read issue 01 before running it in the
same tree, it will silently corrupt the CMake builds.

## 3. Gotchas that cost time in the previous session

1. **`include/config.h` shadowing.** An Autotools `configure` writes
   `include/config.h` into the source tree, and it then overrides every CMake
   build's own `config.h`, silently switching their precision. This is issue 01.
   Symptom: `libnfft3f.so` exporting `nfft_tune` instead of `nfftf_tune`.
2. **Long-double is software binary128 on aarch64** (`LDBL_MANT_DIG == 113`).
   `checkall_ng` in `bx-l` takes minutes; the long-double sweeps take ~40. Any
   test that runs an `O(N*M)` direct NDFT in long double is expensive — that is
   why `check_tune_plan` only runs its real-transform check for `N <= 128`.
3. **`-ffast-math` is on** (the CMake maxopt flags mirror `AX_CC_MAXOPT`).
   Evaluating the same expression at two call sites can differ in the last bit,
   so do not assert exact equality between a returned `attained` and a
   recomputed one. `check_tune_unreachable` nudges the goal by `1e-6` relative
   for exactly this reason.
4. **`select.py` is a stdlib module name.** The model-form selector is
   `selectform.py` for that reason.
5. Work in a git worktree; `.claude/worktrees/<name>` off `feature/tune-m`.

## 4. The measured baseline

`docs/tuning-vs-legacy.md` is the full report; raw data in
`.scratch/sigma-m-study/data/compare-*.csv.gz`, driver `compare.c`, summariser
`compare_report.py`.

It compares `tune_plan` against the legacy geometry `n = 2*next_power_of_2(N)`
where legacy's `m` is found by **searching upward for the smallest cut-off whose
measured error meets the goal** — an oracle the legacy API does not actually
provide. Both sides run through the same `plan_ng` path, so only the parameter
pair differs. 288 configurations: N in {243, 250, 251, 255, 256, 512} × three
shapes (M = N/4, N, 4N) × 2–3 goals × both directions × three precisions.

Result:

- **Accuracy: 288/288 for the tuner, unaided.** Never misses.
- **Speed: median 0.97x — not a win.** By shape:

| shape | M | median speedup | tuner faster |
|---|---|---|---|
| N-dominated | N/4 | 1.31x | 95/96 |
| balanced | N | 0.96x | 42/96 |
| M-dominated | 4N | 0.77x | 3/96 |

The cause: `tune_plan` minimises `n` first, landing at about **0.62×** the
legacy size and paying with a cut-off **2 to 3 larger**. The FFT costs
`O(n log n)`, the node convolution `O(M m)`. Trading `n` down for `m` up is
right only while the FFT dominates. That is issue 02.

Median headroom below the goal is **27x for the tuner against 3x for the oracle
search**. That gap is structural, not a margin that can be tuned away: the
envelope is a worst case over all measured geometries, the oracle uses the
actual error of the one geometry in front of it. `exp(-D*m)` with `D ≈ 4.4` at
`sigma = 2` means `ln(27)/4.4 ≈ 0.75`, i.e. about one extra `m`. Shrinking the
1.25 safety margin does **not** help — `ln(1.25)/4.4 ≈ 0.05` of an `m`.

## 5. The work

| issue | title | status |
|---|---|---|
| [01](issues/01-config-h-shadowing-guard.md) | Guard against `include/config.h` shadowing CMake builds | ready-for-agent |
| [02](issues/02-cost-aware-tune-plan.md) | Take `M` and choose `(n, m)` by cost | ready-for-agent |

Do 01 first — without it, any three-precision verification of 02 is untrustworthy.
