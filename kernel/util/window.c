/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "api.h"
#include "bessel_i0_data.h" /* NFFT_I0_ASYMP_SPLIT */
#include "iplanner.h"

#if defined(DIRAC_DELTA)
static const INT m2K_[] = {0};
#elif defined(GAUSSIAN)
static const INT m2K_[] = {0, 1, 3, 6, 7, 9, 11, 13, 15, 17, 19, 21, 22, 23, 24};
#elif defined(B_SPLINE)
static const INT m2K_[] = {0, 0, 4, 7, 10, 13, 15, 17, 19, 22, 24};
#elif defined(SINC_POWER)
static const INT m2K_[] = {0, 0, 2, 5, 8, 11, 12, 14, 16, 18, 21, 23, 24, 24};
#else /* Kaiser-Bessel */
static const INT m2K_[] = {1, 3, 7, 9, 14, 17, 20, 23, 24};
#endif

/** the PRE_LIN_PSI table size K for cut-off m */
INT Y(m2K)(const INT m)
{
  int j = MIN(((int)(m)), ((int)((sizeof(m2K_) / sizeof(m2K_[0])) - 1)));
  return (INT)((1U << m2K_[j]) * (m + 2));
}

/** the default window cut-off m */
NFFT_INT Y(get_default_window_cut_off)()
{
  return (NFFT_INT)(WINDOW_HELP_ESTIMATE_m);
}

const char *Y(get_window_name)()
{
  return STRINGIZE(WINDOW_NAME);
}

/* Kaiser-Bessel window, scaled so that phi_hut(0) = 1. The per-axis constants
 * are lg_tail = log I0(m b) - m b and i0e_peak_inv = 1/(exp(-m b) I0(m b)). */

/* 1 / I0(m b). The out-of-support branch has no cancellation to preserve, so it
 * divides by the peak, formed on demand: the convolution never reaches it. */
static inline R kb_inv_peak(R b, R lg_tail, R m)
{
  const R xpk = m * b;

  if (xpk <= NFFT_LOG_MAX_SAFE)
    return K(1.0) / Y(bessel_i0)(xpk);

  return EXP(-xpk - lg_tail);
}

/* i0e_peak_inv is cached per axis, so this carries no division. exp(-x)*I0(x)
 * lies in (0, 1] and its reciprocal grows like sqrt(2 pi x), so the scaling
 * cannot overflow at any m b. */
R Y(kb_phi_hut)(R b, R i0e_peak_inv, R m, R n, R k)
{
  const R t = K(2.0) * KPI * k / n;
  const R s = (b - t) * (b + t); /* not b*b - t*t; >= 0 in band, so SQRT is safe */
  const R ra = SQRT(s);
  const R a = m * ra;

  /* I0(a)/I0(m b) = exp(a - m b) * i0e(a) / i0e(m b), with
   * a - m b = -m*t^2/(ra + b) formed without the cancellation. */
  return EXP(-m * t * t / (ra + b)) * Y(bessel_i0_exp_scaled)(a) * i0e_peak_inv;
}

R Y(kb_phi)(R b, R lg_tail, R m, R n, R x)
{
  const R nx = n * x;
  const R a = (m - nx) * (m + nx);

  if (a > K(0.0))
  {
    const R ra = SQRT(a);
    const R u = b * ra;
    /* 0.5*(e^u - e^-u)/I0(m b) = e^(u - m b - lg_tail) * (-expm1(-2u))/2, with
     * u - m b = -b*nx^2/(ra + m). */
    const R umb = -b * (nx * nx) / (ra + m);
    return EXP(umb - lg_tail) * (-EXPM1(K(-2.0) * u)) / (K(2.0) * KPI * ra);
  }

  if (a < K(0.0))
  {
    const R rma = SQRT(-a);
    return (SIN(b * rma) / (KPI * rma)) * kb_inv_peak(b, lg_tail, m);
  }

  return (b / KPI) * kb_inv_peak(b, lg_tail, m);
}

int Y(get_window_id)(void)
{
#if defined(DIRAC_DELTA)
  return NFFT_WINDOW_DIRAC_DELTA;
#elif defined(GAUSSIAN)
  return NFFT_WINDOW_GAUSSIAN;
#elif defined(B_SPLINE)
  return NFFT_WINDOW_B_SPLINE;
#elif defined(SINC_POWER)
  return NFFT_WINDOW_SINC_POWER;
#else /* Kaiser-Bessel */
  return NFFT_WINDOW_KAISER_BESSEL;
#endif
}

/* Runtime window evaluation, dispatched on the NFFT_WINDOW_* ordinal:
* b = pi(2 - N/n) is the per-axis Kaiser-Bessel shape parameter, computed here from (n,N).
 *
 * The KB case is peak-normalized: it self-scales by exp(-log I0(m*b)) so that
 * phi_hut(0) = 1, which keeps the deconvolution grid O(1) instead of tracking
 * the raw I0(m*b) peak which can overflow.
 * The scale is uniform over k and x, so Kaiser-Bessel differs from PHI_HUT/PHI in
 * absolute magnitude only; the other windows are unnormalized, as in the
 * macros. Dirac and unrecognized ordinals return 0; callers must decline. */

/* Kaiser-Bessel evaluation against the per-axis constants of kb_consts: xpk = m*b,
 * lg_peak = log I0(xpk), lt_xpk = lg_peak - xpk == bessel_i0_logtail(xpk). */
static inline R kb_phi_hut_eval(R b, R xpk, R lg_peak, R lt_xpk, INT n, int m, INT k)
{
  R t = K(2.0) * KPI * (R)k / (R)n;
  R s = b * b - t * t; /* >= 0 in band, so SQRT is safe */
  R a = (R)m * SQRT(s);
  if (a > NFFT_I0_ASYMP_SPLIT && xpk > NFFT_I0_ASYMP_SPLIT)
  {
    R lin = -(R)m * t * t / (SQRT(s) + b);
    return EXP(lin + Y(bessel_i0_logtail)(a) - lt_xpk);
  }
  return Y(bessel_i0_scaled)(a, lg_peak);
}

static inline R kb_phi_eval(R b, R xpk, R lg_peak, R lt_xpk, INT n, int m, R arg)
{
  R nx = (R)n * arg;
  R a = (R)m * (R)m - nx * nx;
  if (a > K(0.0))
  {
    R ra = SQRT(a);
    if (xpk > NFFT_I0_ASYMP_SPLIT)
    {
      R umb = -b * (nx * nx) / (ra + (R)m);
      R u = b * ra;
      return K(0.5) * (EXP(umb - lt_xpk) - EXP(-u - xpk - lt_xpk)) / (KPI * ra);
    }
    {
      R u = b * ra;
      return K(0.5) * (EXP(u - lg_peak) - EXP(-u - lg_peak)) / (KPI * ra);
    }
  }
  else if (a < K(0.0))
  {
    R rma = SQRT(-a);
    return SIN(b * rma) * EXP(-lg_peak) / (KPI * rma);
  }
  else
    return (b / KPI) * EXP(-lg_peak);
}

static inline void kb_consts(int window, INT n, INT N, int m,
                             R *b, R *xpk, R *lg_peak, R *lt_xpk)
{
  *b = KPI * (K(2.0) - (R)N / (R)n);
  *xpk = (R)m * *b;
  *lg_peak = Y(bessel_i0_log)(*xpk);
  *lt_xpk = *lg_peak - *xpk;
}

void Y(window_phi_hut_apply)(int window, INT n, INT N, int m, INT k0,
                             R *out, INT count)
{
  if (window == NFFT_WINDOW_KAISER_BESSEL)
  {
    R b, xpk, lg_peak, lt_xpk;
    INT i;
    kb_consts(window, n, N, m, &b, &xpk, &lg_peak, &lt_xpk);
    for (i = 0; i < count; i++)
      out[i] = kb_phi_hut_eval(b, xpk, lg_peak, lt_xpk, n, m, k0 + i);
  }
  else
  {
    INT i;
    for (i = 0; i < count; i++)
      out[i] = Y(window_phi_hut)(window, n, N, m, k0 + i);
  }
}

void Y(window_phi_precompute)(int window, INT n, INT N, int m,
                              const R *x, INT x_stride, INT num_nodes,
                              R *out, INT out_stride)
{
  if (window == NFFT_WINDOW_KAISER_BESSEL)
  {
    R b, xpk, lg_peak, lt_xpk;
    INT j;
    kb_consts(window, n, N, m, &b, &xpk, &lg_peak, &lt_xpk);
    for (j = 0; j < num_nodes; j++)
    {
      R xj = x[(size_t)j * (size_t)x_stride];
      INT c = LRINT(FLOOR((R)n * xj));
      int lj;
      for (lj = 0; lj <= 2 * m + 1; lj++)
      {
        INT idx = c - m + (INT)lj;
        out[(size_t)j * (size_t)out_stride + (size_t)lj] =
            kb_phi_eval(b, xpk, lg_peak, lt_xpk, n, m, xj - (R)idx / (R)n);
      }
    }
  }
  else
  {
    INT j;
    for (j = 0; j < num_nodes; j++)
    {
      R xj = x[(size_t)j * (size_t)x_stride];
      INT c = LRINT(FLOOR((R)n * xj));
      int lj;
      for (lj = 0; lj <= 2 * m + 1; lj++)
      {
        INT idx = c - m + (INT)lj;
        out[(size_t)j * (size_t)out_stride + (size_t)lj] =
            Y(window_phi)(window, n, N, m, xj - (R)idx / (R)n);
      }
    }
  }
}

R Y(window_phi_hut)(int window, INT n, INT N, int m, INT k)
{
  switch (window)
  {
  case NFFT_WINDOW_KAISER_BESSEL:
  {
    R b, xpk, lg_peak, lt_xpk;
    kb_consts(window, n, N, m, &b, &xpk, &lg_peak, &lt_xpk);
    return kb_phi_hut_eval(b, xpk, lg_peak, lt_xpk, n, m, k);
  }
  case NFFT_WINDOW_GAUSSIAN:
  {
    R sigma = (R)n / (R)N;
    R b = (K(2.0) * sigma) / (K(2.0) * sigma - K(1.0)) * ((R)m / KPI);
    R t = KPI * (R)k / (R)n;
    return EXP(-(t * t) * b);
  }
  case NFFT_WINDOW_B_SPLINE:
  {
    if (k == 0)
      return K(1.0) / (R)n;
    {
      R a = (R)k * KPI / (R)n;
      return POW(SIN(a) / a, K(2.0) * (R)m) / (R)n;
    }
  }
  case NFFT_WINDOW_SINC_POWER:
  {
    R sigma = (R)n / (R)N;
    R arg = (K(2.0) * (R)m * (R)k) / ((K(2.0) * sigma - K(1.0)) * (R)n / sigma) + (R)m;
    return Y(bsplines)((INT)(2 * m), arg);
  }
  default:
    return K(0.0);
  }
}

R Y(window_phi)(int window, INT n, INT N, int m, R x)
{
  switch (window)
  {
  case NFFT_WINDOW_KAISER_BESSEL:
  {
    R b, xpk, lg_peak, lt_xpk;
    kb_consts(window, n, N, m, &b, &xpk, &lg_peak, &lt_xpk);
    return kb_phi_eval(b, xpk, lg_peak, lt_xpk, n, m, x);
  }
  case NFFT_WINDOW_GAUSSIAN:
  {
    R sigma = (R)n / (R)N;
    R b = (K(2.0) * sigma) / (K(2.0) * sigma - K(1.0)) * ((R)m / KPI);
    return EXP(-POW(x * (R)n, K(2.0)) / b) / SQRT(KPI * b);
  }
  case NFFT_WINDOW_B_SPLINE:
    return Y(bsplines)((INT)(2 * m), x * (R)n + (R)m) / (R)n;
  case NFFT_WINDOW_SINC_POWER:
  {
    R sigma = (R)n / (R)N;
    return ((R)n / sigma) * (K(2.0) * sigma - K(1.0)) / (K(2.0) * (R)m) * POW(Y(sinc)(KPI * (R)n / sigma * x * (K(2.0) * sigma - K(1.0)) / (K(2.0) * (R)m)), K(2.0) * (R)m) / (R)n;
  }
  default:
    return K(0.0);
  }
}
