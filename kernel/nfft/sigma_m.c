/*
 * Copyright (c) 2002, 2026 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/* Oversampling and cut-off for a 1-D Kaiser-Bessel NFFT of half-width m + 1.
 * Worst-case error model, sigma = n/N, b = pi (2 - 1/sigma),
 * D = 2 pi sqrt(1 - 1/sigma), beta = 2 - 1/sigma, w = m + 1:
 *   E(m) = CT R / I0(w D) + eps (C1 I0(w b) / I0(w D) + C0)
 *   R    = max(|1 - beta exp(i pi w / sigma)|, G beta)
 * truncation, the alias image against the support-edge sample, plus
 * round-off. E has one minimum in m and falls monotonically in sigma.
 * Cost per trafo: COST_FFT n log2 n + COST_NODE M + COST_TAP M (2m + 2). */

#include "config.h"

#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3.h"
#include "infft.h"

#undef X
#define X(name) NFFT(name)

/* The window macros in infft.h fall through to Kaiser-Bessel. */
#if !defined(DIRAC_DELTA) && !defined(GAUSSIAN) && !defined(B_SPLINE) \
    && !defined(SINC_POWER)
#define WINDOW_IS_KAISER_BESSEL 1
#endif

#define CT K(1.138)
#define C1 K(2.445)
#define C0 K(1.352)
#define G K(1.2)

#define COST_FFT K(2.83e-10)
#define COST_NODE K(3.16e-10)
#define COST_TAP K(4.42e-10)

/* Regimes by M/N, and the tie tolerance of the balanced one. */
#define BALANCED_LO K(0.25)
#define BALANCED_HI K(8.0)
#define TIE K(0.05)

#define M_MAX 64
#define BISECT 64

/* The power-of-two part of an FFT size stays at or below 1 KB in bytes, a
 * quarter of the smallest common L1 set period, so FFTW's strides cannot line
 * up in a few cache sets. */
#if defined(NFFT_LDOUBLE)
#define FFT_AMAX 5
#else
#define FFT_AMAX 6
#endif

#if defined(WINDOW_IS_KAISER_BESSEL)

static R log_i0(const R x)
{
  return x + Y(bessel_i0_logtail)(x);
}

static R model_error(const int m, const R sigma)
{
  const R b = KPI * (K(2.0) - K(1.0) / sigma);
  const R D = K2PI * SQRT(K(1.0) - K(1.0) / sigma);
  const R beta = K(2.0) - K(1.0) / sigma;
  const R w = (R)(m + 1);
  const R ld = log_i0(w * D);
  const R c = K(1.0) + beta * beta - K(2.0) * beta * COS(KPI * w / sigma);
  const R ripple = FMAX(SQRT(FMAX(c, K(0.0))), G * beta);

  return CT * ripple * EXP(-ld)
    + EPSILON * (C1 * EXP(log_i0(w * b) - ld) + C0);
}

/* The m at the model minimum, and the error there. */
static int best_m(const R sigma, R *err)
{
  int m;
  R prev = model_error(1, sigma);

  for (m = 2; m <= M_MAX; m++)
  {
    const R e = model_error(m, sigma);
    if (e >= prev)
      break;
    prev = e;
  }
  *err = prev;
  return m - 1;
}

static R floor_error(const R sigma)
{
  R e;
  best_m(sigma, &e);
  return e;
}

/* Smallest sigma in [sigma_min, sigma_max] whose floor meets goal; the caller
 * floors goal at sigma_max's floor. */
static R smallest_sigma(const R goal, const R sigma_min, const R sigma_max)
{
  R lo = sigma_min, hi = sigma_max;
  int i;

  if (floor_error(sigma_min) <= goal)
    return sigma_min;
  for (i = 0; i < BISECT; i++)
  {
    const R mid = K(0.5) * (lo + hi);
    if (floor_error(mid) <= goal)
      hi = mid;
    else
      lo = mid;
  }
  return hi;
}

/* Smallest m >= 1 meeting goal at sigma; the m at the minimum if none does. */
static int smallest_m(const R sigma, const R goal)
{
  R e;
  const int mb = best_m(sigma, &e);
  int m;

  for (m = 1; m < mb; m++)
    if (model_error(m, sigma) <= goal)
      return m;
  return mb;
}

/* Smallest sigma in [lo, hi] at which cut-off m meets goal, or -1. */
static R sigma_for_m(const int m, const R goal, R lo, R hi)
{
  int i;

  if (model_error(m, hi) > goal)
    return K(-1.0);
  if (model_error(m, lo) <= goal)
    return lo;
  for (i = 0; i < BISECT; i++)
  {
    const R mid = K(0.5) * (lo + hi);
    if (model_error(m, mid) <= goal)
      hi = mid;
    else
      lo = mid;
  }
  return hi;
}

/* n = 2^a 3^b 5^c 7^d with 1 <= a <= FFT_AMAX: FFTW's specialised radices. */
static int is_fft_size(INT n)
{
  int a = 0;

  if (n < 2)
    return 0;
  while (n % 2 == 0)
  {
    n /= 2;
    a++;
  }
  if (a < 1 || a > FFT_AMAX)
    return 0;
  while (n % 3 == 0)
    n /= 3;
  while (n % 5 == 0)
    n /= 5;
  while (n % 7 == 0)
    n /= 7;
  return n == 1;
}

/* The smallest efficient FFT size >= n; at most 5 % above n for n >= 1024. */
static INT fft_size(INT n)
{
  while (!is_fft_size(n))
    n++;
  return n;
}

static R cost(const INT M, const INT n, const int m)
{
  return COST_FFT * (R)n * LOG2((R)n) + COST_NODE * (R)M
    + COST_TAP * (R)M * (R)(2 * m + 2);
}

#endif /* WINDOW_IS_KAISER_BESSEL */

void X(sigma_m_1d)(int N, int M, R goal, R sigma_min, R sigma_max, int *n,
  int *m)
{
#if defined(WINDOW_IS_KAISER_BESSEL)
  R lo, n_top_R;
  INT n_lo, n_hi, n_top;
  int m_min, m_max, k;
  INT cn[M_MAX];
  int cm[M_MAX];
  R cc[M_MAX], low;
  int nc = 0, best;

  goal = FMAX(goal, floor_error(sigma_max));
  lo = smallest_sigma(goal, sigma_min, sigma_max);
  n_lo = fft_size((INT)CEIL(lo * (R)N));

  /* The largest efficient size not above sigma_max N, if any is left. */
  n_top_R = FLOOR(sigma_max * (R)N);
  for (n_top = (INT)n_top_R; n_top >= n_lo && !is_fft_size(n_top); n_top--)
    ;
  n_hi = (n_top >= n_lo) ? n_top : n_lo;

  /* Few nodes: the FFT dominates, take the smallest usable sigma. */
  if ((R)M < BALANCED_LO * (R)N)
  {
    *n = (int)n_lo;
    *m = smallest_m((R)n_lo / (R)N, goal);
    return;
  }

  /* One candidate per m, each at the start of its plateau: the first
   * efficient size at which it meets the goal. */
  m_min = smallest_m((R)n_hi / (R)N, goal);
  m_max = ((R)M > BALANCED_HI * (R)N) ? m_min
    : smallest_m((R)n_lo / (R)N, goal);
  for (k = m_max; k >= m_min; k--)
  {
    const R s = sigma_for_m(k, goal, lo, (R)n_hi / (R)N);
    INT c;
    if (s < K(0.0))
      continue;
    c = fft_size((INT)CEIL(s * (R)N));
    if (c > n_hi)
      c = n_hi;
    cn[nc] = c;
    cm[nc] = smallest_m((R)c / (R)N, goal);
    cc[nc] = cost(M, c, cm[nc]);
    nc++;
  }

  /* Many nodes: one candidate, the smallest m. Balanced: the smallest sigma
   * within the tie tolerance of the cheapest. */
  /* nc >= 1: n_hi >= lo N, so m_min meets goal at n_hi / N. */
  low = cc[0];
  for (k = 1; k < nc; k++)
    low = FMIN(low, cc[k]);
  best = -1;
  for (k = 0; k < nc; k++)
    if (cc[k] <= (K(1.0) + TIE) * low && (best < 0 || cn[k] < cn[best]))
      best = k;
  *n = (int)cn[best];
  *m = cm[best];
#else
  /* The model is calibrated for the Kaiser-Bessel window only; other windows
   * keep the library defaults. */
  (void)M;
  (void)goal;
  (void)sigma_min;
  (void)sigma_max;
  *n = 2 * (int)Y(next_power_of_2)(N);
  *m = WINDOW_HELP_ESTIMATE_m;
#endif
}
