# ADR-0008: CONV node sorting

## Status
Accepted (2026-09-14).

## Context
The CONV solvers processed nodes in caller order. On random nodes each node
fetches its tap block cold, so cost per tap grows with the grid; sorted, it
stays flat. Measured before this work (`.scratch/fused-pruned-fft/FINDINGS.md`
sections 8 to 10): 5.2x on CONV at 256^3, 2.5x on the whole 3D transform at
128^3, equal at 4 MiB. Legacy `nfft_plan` sorts behind `NFFT_SORT_NODES`,
sorts again on every trafo and adjoint call, and reaches `x`, `psi` and `f`
through the index.

## Decision
- **Key.** The row-major index of the node's grid cell, axis 0 most
  significant, each axis shifted by `n_t/2` (`Y(conv_node_key)`). This is the
  order of `g` in memory. The shift makes a caller order sorted by cell over
  `[-0.5, 0.5)` monotone. Plain lexicographic only; a tiled key or a different
  grid axis order (`.scratch/grid-axis-order/PRD.md`) waits for measurement.
- **A raced solver, not a flag.** Every CONV registrar, serial and threaded,
  registers caller order (`reg_id` 0) and sorted (`reg_id` 1). MEASURE and
  PATIENT race them; wisdom records the winner.
- **Index once.** The sorted plan builds `perm` with a stable LSD radix sort on
  every wake (ADR-0009), then builds `u` and `psi` in key order. `apply` reads
  `psi`, `u` and `g` in order and uses `perm` only for `f`. Measured before this
  work at 3D, M = 2.1e6: 26.9 ms for the node data this way against 252.6 ms
  through the index. The radix sort is a serial copy local to
  `kernel/conv/conv.c`, not the legacy `kernel/util/sort.c` helper: it opens no
  OpenMP region and carries no dependency on legacy code. A threaded sort waits
  until the `u`/`psi` fills themselves are threaded.
- **Sortedness is part of the CONV key.** `mkproblem_nfft` tests key order on
  its copy of `x` and stores `x_sorted`; the CONV child hashes it and the sorted
  solver declines when it is set. Without the bit, a race on presorted nodes
  blesses caller order and wisdom replays it on random nodes of the same size
  class. It is the one property of `x` the key records, as FFTW records array
  alignment. Paths that never read `x` store 1.
- **ESTIMATE thresholds, by rank.** A sorted plan adds `2M` to `other`. A
  caller-order plan on unsorted nodes adds one operation per tap when the grid
  exceeds `Y(conv_sort_grid_bytes)(d)`. A rank's value is the largest measured
  grid where that rank gains less than 1.10x in either direction, rounded up to
  a power of two; `HUGE_VAL` if it never gains. Rank 4 and up use rank 3.

Times are ns per tap; gain is caller over sorted.

| d | n/axis | grid MiB | M | fwd caller | fwd sorted | adj caller | adj sorted | fwd gain | adj gain |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 16384 | 0.25 | 599186 | 0.354 | 0.503 | 0.450 | 0.568 | 0.70x | 0.79x |
| 1 | 32768 | 0.50 | 599186 | 0.366 | 0.463 | 0.510 | 0.520 | 0.79x | 0.98x |
| 1 | 65536 | 1.00 | 599186 | 0.399 | 0.481 | 0.541 | 0.536 | 0.83x | 1.01x |
| 1 | 131072 | 2.00 | 599186 | 0.416 | 0.481 | 0.544 | 0.540 | 0.87x | 1.01x |
| 1 | 262144 | 4.00 | 599186 | 0.409 | 0.494 | 0.561 | 0.575 | 0.83x | 0.98x |
| 1 | 524288 | 8.00 | 599186 | 0.437 | 0.624 | 0.565 | 0.574 | 0.70x | 0.98x |
| 1 | 1048576 | 16.00 | 599186 | 0.674 | 0.723 | 0.846 | 0.604 | 0.93x | 1.40x |
| 1 | 2097152 | 32.00 | 599186 | 1.875 | 0.526 | 2.130 | 0.581 | 3.56x | 3.67x |
| 1 | 4194304 | 64.00 | 599186 | 2.514 | 0.625 | 2.770 | 0.654 | 4.02x | 4.24x |
| 1 | 8388608 | 128.00 | 599186 | 2.941 | 0.674 | 3.060 | 0.720 | 4.36x | 4.25x |
| 1 | 16777216 | 256.00 | 599186 | 3.202 | 0.887 | 3.535 | 1.132 | 3.61x | 3.12x |
| 2 | 128 | 0.25 | 42799 | 0.314 | 0.246 | 0.445 | 0.383 | 1.27x | 1.16x |
| 2 | 182 | 0.51 | 42799 | 0.319 | 0.242 | 0.440 | 0.366 | 1.32x | 1.20x |
| 2 | 256 | 1.00 | 42799 | 0.381 | 0.251 | 0.503 | 0.419 | 1.52x | 1.20x |
| 2 | 362 | 2.00 | 42799 | 0.357 | 0.243 | 0.459 | 0.367 | 1.47x | 1.25x |
| 2 | 512 | 4.00 | 42799 | 0.427 | 0.338 | 0.651 | 0.521 | 1.26x | 1.25x |
| 2 | 724 | 8.00 | 42799 | 0.366 | 0.265 | 0.471 | 0.396 | 1.38x | 1.19x |
| 2 | 1024 | 16.00 | 42799 | 0.887 | 0.672 | 1.845 | 1.606 | 1.32x | 1.15x |
| 2 | 1448 | 31.99 | 42799 | 1.503 | 0.387 | 1.671 | 0.491 | 3.89x | 3.40x |
| 2 | 2048 | 64.00 | 42799 | 2.031 | 0.772 | 2.432 | 1.790 | 2.63x | 1.36x |
| 2 | 2896 | 127.97 | 42799 | 2.725 | 0.985 | 2.682 | 1.227 | 2.77x | 2.19x |
| 2 | 4096 | 256.00 | 42799 | 2.569 | 1.253 | 2.868 | 2.108 | 2.05x | 1.36x |
| 3 | 26 | 0.27 | 4096 | 0.313 | 0.262 | 0.470 | 0.404 | 1.19x | 1.16x |
| 3 | 32 | 0.50 | 4096 | 0.322 | 0.324 | 0.489 | 0.483 | 0.99x | 1.01x |
| 3 | 40 | 0.98 | 4096 | 0.348 | 0.273 | 0.484 | 0.413 | 1.27x | 1.17x |
| 3 | 50 | 1.91 | 4096 | 0.303 | 0.267 | 0.477 | 0.422 | 1.13x | 1.13x |
| 3 | 64 | 4.00 | 4096 | 0.352 | 0.363 | 0.505 | 0.499 | 0.97x | 1.01x |
| 3 | 80 | 7.81 | 4096 | 0.343 | 0.331 | 0.485 | 0.455 | 1.04x | 1.07x |
| 3 | 102 | 16.19 | 4096 | 0.537 | 0.333 | 0.681 | 0.498 | 1.62x | 1.37x |
| 3 | 128 | 32.00 | 4096 | 1.198 | 0.449 | 1.497 | 0.696 | 2.67x | 2.15x |
| 3 | 162 | 64.87 | 4096 | 2.191 | 0.896 | 2.522 | 1.016 | 2.45x | 2.48x |
| 3 | 204 | 129.54 | 4096 | 2.599 | 1.327 | 2.771 | 1.590 | 1.96x | 1.74x |
| 3 | 256 | 256.00 | 4096 | 2.368 | 1.605 | 2.601 | 2.111 | 1.48x | 1.23x |

| d | threshold bytes |
|---|---|
| 1 | 16777216 |
| 2 | 131072 (the program printed 0: gains stayed above 1.10x at every measured grid, so the rule's largest-losing-grid search found no row) |
| 3 | 8388608 |

Host: aarch64, 10 cores, dev container; ratios only. ESTIMATE is a
heuristic on other hosts; MEASURE is not affected.

Table above is run 1 (`.scratch/conv-node-sort/calib.md`); a second sweep
(`.scratch/conv-node-sort/calib-rerun.md`) is kept in `.scratch`, not
committed, and agrees with run 1 on ranks 2 and 3.

**How the committed values were chosen** (`kernel/conv/conv.c`, rank 1 =
16777216.0, rank 2 = 131072.0, rank 3 = 8388608.0):
- **Rank 1: 16 MiB from run 1.** Run 2 gave 8 MiB instead, because the 16 MiB
  row measured 0.93x in run 1 and 1.29x in run 2. The larger value was kept so
  ESTIMATE does not sort where one run measured a loss.
- **Rank 2: never below 1.10x in either direction**, at any measured grid
  (2^18 to 2^28 bytes), in both runs. The value is set to 131072 (half the
  smallest measured grid), so ESTIMATE sorts every rank-2 grid on random
  nodes. `HUGE_VAL` does not apply here: that spelling is for a rank that
  never gains, and rank 2 always gains.
- **Rank 3: 8 MiB in both runs.** The decisive 7.81 MiB row was within 0.03 of
  1.10 in run 1 (1.07x adjoint), so the sweep was re-run; run 2 (1.06x
  adjoint) confirmed it.

## Consequences
- **Adjoint bits depend on the family.** Sorted and caller order accumulate
  each cell in a different node order; forward results agree to rounding.
  Under MEASURE the winner, and with it the adjoint's last bits, can change
  from run to run; FFTW has the same property. Within one family serial and
  threaded stay bit-exact (`compare_conv`, `tests/ng_omp.c`).
- **Wisdom for sorted and presorted nodes is separate.** A wisdom-only plan for
  presorted nodes misses on wisdom made from random ones, and the reverse.
  The vocabulary tag in `config_signature` moves to `"plnr-flags-v2"`, because
  the CONV hash gained `x_sorted`.
- **Measured planning costs one more CONV wake and apply per family.** Peak
  table memory stays one candidate's (ADR-0009).
- **No benchmark or accuracy report moves:** `plan_ng` is used by no benchmark
  and no binding, and the accuracy report comes from the legacy suites.
- **Legacy is unchanged.**

## Not decided here
- Tiled or Morton keys, and grid axis reordering.
- A public flag to force or forbid sorting.
- A threshold in M or tap count. The threshold depends only on grid bytes, each
  rank was calibrated at one M, and the rank-3 sweep is sparse; node density
  also moves the gain. This waits for measurement.
