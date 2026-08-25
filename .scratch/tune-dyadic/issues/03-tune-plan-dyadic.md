# 03 — `X(tune_plan_dyadic)`, API, build, tests

Status: done

## Problem

Nothing in the library offers the dyadic ladder. The earlier full search
walked every even 5-smooth `n` in `[5/4 N, 4N]` and had to carry a radix-2
tie-break to undo the FFT-size penalty that walk creates.

## What to do

Add, in `kernel/nfft/tune.c` beside the existing entry points:

```c
int X(tune_plan_dyadic)(NFFT_INT N, NFFT_INT M, int adjoint, R goal,
                        NFFT_INT *n, int *m, R *attained);
```

Algorithm, with `2^k = next_power_of_2(N)`:

1. cap the goal at the accuracy floor of rung 2, the widest grid on offer;
2. for `j` in 0, 1, 2 with `n = 2^j * 2^k`: skip when `n <= N`, when
   `4*n < 5*N`, or when `m_max(n) < 1`; take the smallest cut-off whose
   band-`j` model error meets the capped goal, skipping the rung when none
   does;
3. return the rung with the smallest `n*log2(n) + (4/5)*M*(2m+2)`.

No tie-break — every candidate is a power of two.

Return codes: 1 goal met, 0 capped goal met instead, -1 bad arguments.

Also:

- declare it in `include/nfft3.h` inside `NFFT_DEFINE_PLANNER_API`;
- add the six per-band coefficient sets from issue 02;

## Tests

In `tests/tune_ng.c`, registered in `tests/check_ng.c`:

- `dyadic_is_dyadic` — the returned `n` is always `2^j * next_power_of_2(N)`
  for some `j` in {0, 1, 2}, over a spread of N covering `t` in (1, 2].
- `dyadic_meets_goal` — run the real transform against the planner's own
  direct NDFT at the returned pair and confirm the goal is met. Keep it to
  `N <= 128` in long double, as `check_tune_plan` does.
- `dyadic_cost` — the returned pair costs no more than any other legal rung,
  under the same cost expression.
- `dyadic_never_worse_than_legacy` — the returned cost is no greater than the
  cost of rung 1 with the cut-off the model gives there. This is the property
  that makes the acceptance bar structural rather than empirical.
- `dyadic_bad_args` — null pointers, `N < 1`, `M < 1`, non-positive goal.
- `dyadic_capped` — a goal below the floor at rung 2 returns 0 with a pair
  meeting the capped goal.

## Watch out

- `-ffast-math` is on: do not assert exact equality between a returned
  `attained` and one recomputed at another call site.
- `N` a power of two gives `t = 1`, so rung 0 has `n = N` and must be skipped
  by the `n <= N` test. Cover it.
- `N = 1` leaves rung 1 at `n = 2`, where `m_max = 0`. Cover it.
- Build wiring: `tune.c` is already in `kernel/nfft/Makefile.am` and
  `kernel/CMakeLists.txt`, so no new file means no new wiring.
