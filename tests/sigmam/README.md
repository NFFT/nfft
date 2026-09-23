# sigmam: the model behind `nfft_sigma_m_1d`

`kernel/nfft/sigma_m.c` picks the grid size `n` (`sigma = n/N`) and cut-off `m`
of a 1-D Kaiser-Bessel NFFT for a worst-case accuracy goal at least cost. Its
constants are fitted here. `kbmodel.py` is the reference implementation, and
the C port returns the same `(n, m)`.

## Model

Error measure `max_j |f_j - s_j| / ||f_hat||_1`, sup over `f_hat`, which is
attained at `f_hat = e_k`. Window half-width `w = m + 1`. With
`b = pi (2 - 1/sigma)`, `D = 2 pi sqrt(1 - 1/sigma)`, `beta = 2 - 1/sigma`:

    E(m) = CT R / I0(w D)  +  eps (C1 I0(w b) / I0(w D) + C0)
    R    = max(|1 - beta exp(i pi w / sigma)|, G beta)

- **Truncation, first term.** At the band edge the alias image lands where
  `phihat(b) = I0(0) = 1`. It interferes with the dropped sample on the
  support edge, `phi(w) = b/pi`. The result is a ripple with period `2 sigma`
  in `m`. Just inside the edge the alias image is gone, and the rest of the
  sin tail adds to the sample: `G beta`.
- **Round-off, second term.** One part is amplified by the deconvolution,
  `1/phihat(band edge)`. The other part is not.

`E` falls with `m`, has one minimum, and at fixed `m` falls monotonically in
`sigma`.

Cost per trafo, in seconds:

    t = fft n log2 n + node M + tap M (2m + 2)

`n` is always an efficient FFT size (`fft_size`): `2^a 3^b 5^c 7^d` with a
power-of-two part of at most 1 KB in bytes (`a <= 6` in double), a quarter of
the smallest common L1 set period. Larger power-of-two parts make FFTW's strides
fall into a few cache sets, which is 3-5x slower; the next efficient size is at
most 5 % above `n`. A grid of at most 64 KB fits in cache on any usual host,
where radix 2 is FFTW's fastest: there `n` rounds up to the next power of two
`P` when `P <= 1.25 n`. Within that factor the power of two measured faster,
beyond it the smaller smooth size.

## Choice

1. Raise the goal to the floor at `sigma_max`, then find the smallest usable
   `sigma` by bisection.
2. Choose by the regime of `M/N`:
   - **few nodes** (`M < N/4`): that smallest sigma.
   - **many nodes** (`M > 8N`): the smallest `m` reachable at `sigma_max`, at
     the smallest sigma where it suffices.
   - **balanced** (in between): one candidate per `m` at the start of its
     plateau, and the smallest sigma within 5 % of the cheapest.

## Files

| file | what |
|---|---|
| `kbmodel.py` | model, the three blocks, `fft_size`, self-check; `--c` prints the C constants |
| `kbmodel.json` | the constants |
| `kbfit.py` | fits `CT, C1, C0, G`; validates both blocks against the sweeps |
| `kbcostfit.py` | fits the cost constants from `costprobe.c` timings |
| `kbties.py` | the tie map behind the regime bounds (`choice`) |
| `run_worst.sh`, `floor.c`, `floorgen.py` | worst-case sweeps: one transform per `k`, mpmath references |
| `exact_scan.py` | exact-arithmetic truncation, max over theta and node offset |
| `costprobe.c`, `data/cost-configs.txt` | transform timings |
| `costplot.py` | cost landscape plot |
| `data/` | the sweeps and timings the constants were fitted to |

## Refit

After a change to the window, the transform, or the floating-point paths:

```sh
cd tests/sigmam
# 1. worst-case sweeps and exact truncation (tens of minutes; quad is slow)
./run_worst.sh <src> <tree-double> <tree-float> <tree-longdouble> /tmp/rw
# 2. error constants: prints the fit and validation, --write updates kbmodel.json
uv run --with scipy python kbfit.py --write
# 3. cost constants, from timings on the host of interest
gcc -O2 -std=gnu99 -I<src>/include costprobe.c -L<lib> -lnfft3 -lfftw3 -lm -o costprobe
./costprobe < data/cost-configs.txt | gzip > data/cost-trafo-d.csv.gz
uv run --with numpy python kbcostfit.py --write
# 4. regime bounds: set "choice" in kbmodel.json to enclose the * band
python3 kbties.py
# 5. self-check, then paste into kernel/nfft/sigma_m.c
python3 kbmodel.py
python3 kbmodel.py --c
```

Then rebuild and run `make check`. The test `nfft_sigma_m_1d` runs the chosen
`(n, m)` against every `k` and fails if the worst case misses the goal.

`kbfit.py` fits by linear program: the model, less a 5 % margin, must stay
above every measured point from `m = 1` to the model's best `m`. `G` enters
non-linearly, so it is scanned. `kbcostfit.py` leaves out FFT sizes that are
irregular on the timing host.
