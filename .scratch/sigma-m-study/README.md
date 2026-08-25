# Kaiser-Bessel sigma/m study

Where the constants in `kernel/nfft/tune.c` come from, and what they are worth.

## The question

For a one-dimensional Kaiser-Bessel NFFT of bandwidth `N` at `M` nodes, pick
the oversampled size `n` and the window cut-off `m` that reach a requested
accuracy most cheaply, in the measure the test suite uses:

    forward   max_j |f_j - s_j| / ||f_hat||_1
    adjoint   max_k |fhat_k - s_k| / ||f||_1

## The dyadic ladder

`next_power_of_2(2N)` is `2*next_power_of_2(N)` for every `N`, so the legacy
default grid is one rung of

    n = 2^j * next_power_of_2(N),   j in {0, 1, 2}

Writing `t = next_power_of_2(N)/N`, which lies in `(1, 2]`, the rungs
oversample by `t`, `2t` and `4t` -- three disjoint bands. Rung 1 is the legacy
choice. Rung 0 halves the FFT work and pays with a larger cut-off, and is legal
only when `t >= 5/4`, the oversampling floor below which the attainable
accuracy collapses. Rung 2 doubles the FFT work to buy a smaller cut-off.

Which rung wins turns on the node count: the FFT costs `O(n log n)` and the
node convolution `O(M*(2m+2))`, so many nodes argue for a wider grid and a
smaller cut-off. `X(tune_plan_dyadic)` costs all three under
`n*log2(n) + (4/5)*M*(2m+2)` and returns the cheapest that reaches the goal.

Every rung is a power of two. That matters beyond tidiness: FFT
implementations optimise radix-2 hardest, and `nsweep`-style measurements put
other smooth sizes 5 to 40 % above a power of two, so `n = 480` loses to
`n = 512` though it is the smaller grid. A search over general smooth sizes has
to carry a tie-break against that effect. A dyadic ladder cannot meet it.

The weight `4/5` is an operation count, not a measurement: about `5*n*log2(n)`
real operations for the FFT, four per window sample. `costfit.c` measures
timings to check the ranking it produces, but a fitted weight would bind the
library to one machine's cache behaviour.

## The model

Two effects oppose each other, both set by `sigma = n/N` alone. With

    b = 2*pi*(1 - 1/(2*sigma))          Kaiser-Bessel shape parameter
    D = 2*pi*sqrt(1 - 1/sigma)          truncation decay rate
    A = b - D                           roundoff amplification rate

the window is truncated to `2m+2` samples and what is cut off decays like
`exp(-D*m)`. The deconvolution divides by the window's Fourier transform,
smallest at the band edge, which amplifies roundoff like `exp(A*m)`. Both
rates follow from `b^2 - (pi/sigma)^2 = D^2`.

The fitted formula keeps those rates exactly and fits only prefactors:

    T = a * u^p * m^r * exp(-D*m)      * N^tn * M^tm    u = 1 - 1/sigma
    U = c * e * u^q * exp(alpha*A*m)   * N^un * M^um    e = machine epsilon
    E = T + U

The `N` and `M` powers are the measure's own normalisation, not curve fitting.
The forward error is a max over `M` nodes of a sum of `N` terms whose phases
cancel, divided by an `l1` norm over `N` terms, so it falls roughly like
`N^-1/2` and rises slowly with `M`; the adjoint swaps the two.

`E` falls, reaches a minimum, and rises again. Since `b > D` at every finite
sigma, `A > 0` and the minimum always exists, so there is an accuracy floor no
cut-off can beat. Balancing the two terms puts it near `eps^(D/b)`; run
`cap.py` for the table.

**One coefficient row per rung.** Each row has to dominate only its own band of
sigma, where a single set would have to dominate all of `[5/4, 8]` at once.
The row is selected by `j` rather than by testing sigma, because rung 0 and
rung 1 meet at `sigma = 2` when `t = 1`.

On the widest band `alpha` is pinned at 1 rather than fitted. It corrects the
exactly derived rate `A`, and `A` falls from 0.96 at `sigma = 5/4` to 0.056 at
4 and 0.012 at 8, so over `m <= 32` the branch is flat and the least squares
reads noise. Fitting it there produced cut-offs up to 29 above the measured
need; pinning it brings the worst case back to +1.

## The input distribution

`gsweep.c` and `compare.c` draw both input vectors with real and imaginary
parts on `[0, 1)`, matching `Y(vrand_unit_complex)` -- what the library's own
accuracy tests present and what any caller's data may look like. Drawing them
centred on zero instead is **not** equivalent:

| geometry | forward, `[0,1)` | forward, centred | ratio |
|---|---|---|---|
| N=128, sigma=2, m=1 | 9.81e-03 | 3.77e-03 | 2.60 |
| N=128, sigma=2, m=3 | 1.14e-06 | 6.19e-07 | 1.83 |
| N=256, sigma=2, m=5 | 9.74e-11 | 5.52e-11 | 1.77 |

Uncentred data has a coherent component, so the forward transform peaks where
the phases line up and the error peaks with it, while the `l1` norm it is
divided by does not. The adjoint measure is within 10 % either way: its input
sits at the nodes and its norm runs over the same index, so a constant offset
scales both halves of the ratio alike.

An envelope fitted to the centred draw does not hold for the uncentred one --
measured on this sweep, such a fit misses the goal 193 times. That is why the
draw here matches the library's.

**Neither draw is a bound.** The measure `max_j |f_j - s_j| / ||fhat||_1` has a
finite supremum over all inputs -- it is the largest absolute entry of the
error matrix, which does not depend on the input at all -- and random data
approaches it from below whatever the distribution. The constants here are
therefore calibrated for random input of the shape the library presents, not
for every input. See `.scratch/tune-dyadic/issues/05-input-distribution-in-the-fit.md`.

## The scripts

| file | what |
|---|---|
| `gsweep.c`, `run_gsweep.sh` | measure the error over the `(sigma, N, M, m)` grid against a long-double direct NDFT, one binary per precision. The reference depends on the nodes and data but not on sigma or m, so it is computed once per `(N, M, trial)` and reused across the whole plane. The sigma list is the seventh argument. |
| `gfit.py` | the fit itself: split each series at its minimum, regress both branches in log space, then raise the prefactors to the tightest pair dominating every measurement, plus a 1.25 margin. `alpha_fixed` pins the roundoff exponent where `A*m` is too small to identify it. |
| `dfit.py` | fit once per rung's band, validate each against the measured smallest sufficient cut-off, and emit the two three-row C tables. |
| `dgate.py` | replay the rung choice against measured error, with legacy given an oracle cut-off. Answers whether the ladder pays before any of it is built. |
| `dvalidate.c` | take the pair the tuner returns and measure it on draws the fit never saw, counting misses. |
| `compare.c`, `run_compare.sh`, `dcompare_report.py` | the head-to-head against the legacy grid, rendered to `docs/tuning-dyadic.md`. |
| `cap.py` | the accuracy floor as a function of sigma and precision. |
| `costfit.c`, `costfit.py`, `run_costfit.sh` | check the `4/5` cost weight's ranking against measured timings. |

## Reproducing

Three CMake trees, one per precision:

```sh
cmake -S . -B bx-d                              -DNFFT_ENABLE_TESTS=ON \
      -DNFFT_ENABLE_OPENMP=OFF -DNFFT_ENABLE_EXAMPLES=OFF \
      -DNFFT_ENABLE_APPLICATIONS=OFF -DCMAKE_BUILD_TYPE=Release
cmake -S . -B bx-f -DNFFT_ENABLE_FLOAT=ON        ...same...
cmake -S . -B bx-l -DNFFT_ENABLE_LONG_DOUBLE=ON  ...same...
```

Then the sweep, the fit, and the head-to-head:

```sh
sh run_gsweep.sh "$PWD" /tmp/gs 32 5 "" "" \
   "1.25,1.3,1.35,1.4,1.5,1.6,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0"
python3 dfit.py /tmp/gs/gsweep-*.csv        # prints the C tables and validates
sh run_compare.sh "$PWD" /tmp/cmp 40
python3 dcompare_report.py ../../docs/tuning-dyadic.md /tmp/cmp/compare-*.csv
```

Long double is software binary128 on aarch64, so its sweep leg takes tens of
minutes; the other two are a couple of minutes each.

## Data

| file | what |
|---|---|
| `data/gsweep-unit-{d,f,l}.csv.gz` | the sweep the constants are fitted to, 17 sigma values covering all three bands |
| `data/compare-dyadic-{d,f,l}.csv.gz` | the head-to-head behind `docs/tuning-dyadic.md` |
| `data/dvalidate-{d,f,l}.csv.gz` | fresh-draw accuracy of the returned pairs |
| `data/costfit-{d,f,l}.csv.gz` | the timings behind the cost weight check |
