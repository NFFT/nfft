# hostcheck: is the sigma/m heuristic a win on this host?

`nfft_sigma_m_1d` (`kernel/nfft/sigma_m.c`) picks the grid size `n` and cut-off
`m` of a 1-D NFFT for an accuracy goal. Its cost model and its FFT size rule
were tuned on one machine. This check shows whether its grid choice also beats
the legacy grid on another host or architecture.

## Run

From the repository root:

```sh
tests/sigmam/hostcheck/run.sh quick    # 48 cases, a few minutes
tests/sigmam/hostcheck/run.sh full     # 420 cases, 30-60 minutes
```

The script builds its own copy of the library (CMake, double, serial) in
`build/hostcheck/`, compiles the driver, runs the cases and prints a summary.
Close other heavy programs while it runs, because the results are timings.

### Requirements

- a C99 compiler, CMake, and Python 3 (standard library only)
- FFTW3 in double precision, with headers (`libfftw3-dev` on Debian and
  Ubuntu, `fftw` in Homebrew)

If CMake or the compiler do not find FFTW, pass its install prefix:

```sh
FFTW_DIR=/opt/homebrew tests/sigmam/hostcheck/run.sh quick
```

`CC`, `CFLAGS` and `LDFLAGS` are honoured too. A second argument sets the
output directory in place of `build/hostcheck`.

## What it measures

Each case is a triple `(N, M, goal)`: `N` Fourier coefficients, `M` random
nodes, and a worst-case accuracy goal. There are two arms, run on the same
random input:

| arm | n | m |
|---|---|---|
| helper | from `nfft_sigma_m_1d(N, M, goal, 1.25, 2.5, ...)` | from the helper |
| legacy | `2 * next_power_of_2(N)`, the grid `nfft_init_1d` uses | the helper's `m` |

Both arms use `nfft_init_1d`'s flags (`PRE_PHI_HUT | PRE_PSI`,
`FFTW_ESTIMATE`). The time is trafo plus adjoint: the best of 5 means, each
mean over at least 50 ms. Precomputation is not timed.

The error is `max |f_j - s_j| / ||f_hat||_1` over 100 of the nodes, against an
exact reference. Nodes are multiples of `2^-30`, so the reference is exact in
double on every platform.

The cases cover:

- `N` from 64 to 262144
- `M = N/8`, `N` and `8N`
- goals `1e-4`, `1e-8`, `1e-12` and `1e-14`

## Reading the summary

- **helper error above the goal** must be 0. Any other count means the error
  model is wrong on this platform. Report those cases first.
- **win / tie / loss** compares times: the legacy arm's time divided by the
  helper's. Within 2 %, or when both arms chose the same `n`, it is a tie.
- **by helper grid** separates the two FFT size rules:
  - small grids (`n <= 4096`, 64 KB): `n` rounds up to a power of two when that
    is at most 1.25 times the smooth size.
  - large grids: `2^a 3^b 5^c 7^d` with a small power-of-two part.
- **losses** lists every case the helper lost, slowest first.

Losses concentrated in small grids point at the 1.25 rounding factor. Losses
in large grids point at the power-of-two limit or the cost model.

`reference-full.txt` is the full summary from the machine the heuristic was
tuned on (Linux VM on Apple Silicon, FFTW 3.3.10 without NEON). Compare yours
with it.

## Sharing results

Send `build/hostcheck/summary-full.txt` and `build/hostcheck/results-full.csv`.
Their header lines record the host, CPU, compiler, FFTW version and NFFT
commit.
