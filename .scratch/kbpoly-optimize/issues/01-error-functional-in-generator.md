# 01: Error functional of a piecewise polynomial window in the generator

Status: done (feature/kbpoly-optimize)

Part of: `.scratch/kbpoly-optimize/PRD.md` (R1). Background: `analysis.md`
section 3; prototype `kbopt.py` (class `Model`).

## Task

Add a module `tests/kbpolygen/functional.py` that evaluates, for a window
given as one polynomial per unit cell on the stencil cells `-R .. R-1`
(`R = m + 1`), in grid units with offset `t` and frequency `xi = k / n`:

- `D(xi, t) = sum_a p_a(t) e^(-2 pi i xi (a + t))`,
- `C(xi) = I0(R sqrt(b^2 - (2 pi xi)^2)) / I0(R b)`, the library's `PHI_HUT`,
  with `b = pi (2 - 1/sigma)`,
- `J_inf = max |D/C - 1|` and `J_2 = sqrt(mean |D/C - 1|^2)` over
  `|xi| <= 1/(2 sigma)` and `t in [0, 1)`, on a grid of at least `160 x 160`,
- the polynomial part `max |D_poly - D_KB| / |C|`, where `D_KB` is the same
  sum with the exact truncated Kaiser-Bessel values, and the ratio of the
  transform's `J_inf` with the polynomial window to the one with the exact
  window.

Inputs are Chebyshev coefficients per cell (so the generator's existing
`degrees()` machinery can hand over its truncated series) or a callable for
the window values. Use `numpy` in double; the generator already runs under
`uv run --with mpmath==1.3.0`, add `--with numpy` to the documented command.
`mpmath` stays for the Chebyshev coefficients of the exact window.

## Acceptance

- For the exact truncated window, `J_inf / (1 / I0(2 pi R sqrt(1 - 1/sigma)))`
  lies in `[1.10, 2.81]` over `sigma` in `{1.1, 1.5, 2, 3, 4}` and `m` in
  `{2, 3, 4}` (the constant `c` the generator's README documents).
- Reproduces, to three significant digits, the `exact truncated KB` and
  `chebKB` columns of `results/res_kbopt_py_2_2.txt`, `..._2_3.txt`,
  `..._4_2.txt` and `..._1_5_4.txt` at every degree listed.
- A `pytest` test in `tests/kbpolygen/tests/` covers both points.
- No change to `kernel/util/kbpoly_data.h` in this issue.
