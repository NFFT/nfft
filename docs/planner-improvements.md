# Planner API NFFT: improvements over the legacy plans

What the new planner API (`nfft_plan_ng_*`, `kernel/nfft/`, `kernel/planner/`)
does better than the legacy `nfft_plan` API. The legacy API stays as it is; this
list records what a caller gains by moving to the new one. Keep adding to it.

## Geometry

### Unit axes are elided.

The guru removes every axis with `N[t] == 1` from the problem before planning,
and gathers the node coordinates in step. Such an axis contributes one frequency
`k[t] = 0` and a phase of 1, so it changes no result. The legacy odometer
carries on `k[t] == N[t]/2 - 1`, which for `N[t] == 1` is `-1` and never fires,
so the legacy transforms give wrong results for these geometries. The new API
accepts any mix of unit and non-unit axes.

### Rank-0 problems are supported.

When every axis has unit length, elision leaves rank 0: one coefficient and no
phase. A dedicated solver serves it — the forward transform broadcasts
`f_hat[0]` to all `M` nodes, the adjoint sums the `f_j`. The legacy API has no
such case.

### Odd axis lengths are supported.

The frequency range of an odd axis is the symmetric sum `k = -(N-1)/2 ..
(N-1)/2`. The odometer uses the threshold `(N[t]-1)/2` instead of the legacy
`N[t]/2 - 1`. The two agree for even `N`, so the even path is unchanged, but
only the first is right for odd `N`.

### Even axis lengths have two NDFT variants.

Per axis, `NFFT_NDFT_TYPE_I` sums `k = -N/2 .. N/2-1`, which is what the legacy
implementation computes. `NFFT_NDFT_TYPE_II` sums the shifted, ascending range
`k = -N/2+1 .. N/2`. Type-II is a uniform `+1` shift of the type-I frequencies,
so the odometer is shared; only the effective frequency changes. Odd axes have
one variant, defined as type-I. The variant is chosen per axis, so mixed
problems are legal.

## Accuracy and range

### The Kaiser-Bessel window is peak-normalized.

`phi_hut` self-scales by `exp(-log I0(m*b))`, so `phi_hut(0) == 1`. The
deconvolution divides the coefficients by `phi_hut` once per axis, so with the
raw window the grid tracks the `I0(m*b)` peak to the power `d`. In float that
leaves the representable range at `d >= 2` or `3`. The scale is uniform over `k`
and `x`, so it cancels through the transform and changes the result by nothing.

### The Kaiser-Bessel evaluation avoids cancellation.

The exponent `m*(sqrt(b^2 - t^2) - b)` subtracts two nearly equal numbers near
the middle of the band. It is evaluated in the equal form
`-m*t^2 / (sqrt(b^2 - t^2) + b)`, which adds instead of subtracts. The
`sinh`-like part of `phi` is formed from its two exponentials with the peak
already factored out, for the same reason.

### The direct NDFT uses a blocked phase recurrence.

One reduced phase serves a block of 32 terms; inside a block the terms advance
by a complex multiply with the fixed step `exp(-+i 2 pi x)`. This replaces two
transcendental calls per term with two per block, and it keeps the argument of
`COS`/`SIN` small, so the phase error does not grow with `N`. The legacy direct
transform calls `COS`/`SIN` once per term on an argument up to `pi*N`, where
libm argument reduction dominates its error. Blocking is the only direct NDFT in
the new API, at every rank. Measured in this dev container (double, `-O3
-ffast-math`): 1D `N = 8192`, `M = 64`, nodes near `+1/2` gives a forward error
of `7.7e-15` against `5.1e-13` for the legacy direct, and is about 8% faster; 2D
`N = 96^2`, `M = 3000` runs 4.2 times faster than the same kernel unblocked.

## Interface

### The window is a runtime argument.

`nfft_plan_ng_guru` takes an `NFFT_WINDOW_*` ordinal, and all four windows
(Kaiser-Bessel, Gaussian, B-spline, sinc-power) work in the fast path through a
window vtable. The legacy library fixes one window at configure time
(`--with-window=`), so a program that needs two windows needs two builds.
`nfft_get_window_id()` still returns the compile-time default for callers that
want it.

### The planner measures and picks the algorithm.

Planning races the applicable solvers on the caller's own nodes and keeps the
fastest, the FFTW model. `NFFT_ESTIMATE` picks by an analytic cost model
instead. The legacy API asks the caller to choose the algorithm and the
precompute strategy through `PRE_*` flags, with no measurement.

### Decisions are cached as wisdom.

A measured decision is keyed by a 128-bit hash of the structural size class of
the problem — never by node coordinates — and reused for every later problem in
that class. Wisdom exports to and imports from a string or a file. An import is
rejected unless the library's solver registry and release reproduce the
signature in the file, so a stale cache cannot select a solver that no longer
means the same thing.

### The chosen plan is printable.

`nfft_fprint_plan` prints the plan as an S-expression that names the winning
solver and its children, so the algorithm actually in use is visible. The legacy
plan has no such report.

### Bad arguments return NULL.

The guru validates its arguments and returns `NULL`, in release builds too. This
covers insufficient oversampling: unless `NFFT_NO_FAST_NATIVE` is set, every
axis with `N[t] > 1` must satisfy `n[t] > N[t]`. `sigma <= 1` is refused rather
than served silently by a direct transform, which would cost `O(N*M)`. The
legacy API guards such cases with assertions that release builds compile out.

### The plan copies the nodes.

`x` is copied into a plan-owned buffer at construction, so the caller keeps
ownership and may free or reuse the array as soon as the guru returns. The plan
never writes to or frees the caller's `x`.

### Transforms can run on caller-supplied arrays.

`nfft_execute_on` and `nfft_execute_adjoint_on` take `f_hat` and `f` explicitly,
so one plan serves many array pairs without a new plan and without touching the
arrays bound at construction.

### The fast NFFT is composed from planned parts.

The fast solver states the deconvolution, the FFT, and the convolution as their
own sub-problems and plans each one through the same planner. A new
deconvolution or convolution algorithm is registered once and becomes available
to every rank, without touching the fast solver.
