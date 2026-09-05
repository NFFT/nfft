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
#include "iplanner.h"

#if defined(DIRAC_DELTA)
  static const INT m2K_[] = {0};
#elif defined(GAUSSIAN)
  static const INT m2K_[] = {0, 1, 3, 6, 7, 9, 11, 13, 15, 17, 19, 21, 22, 23, 24};
#elif defined(B_SPLINE)
  static const INT m2K_[] = {0, 0, 4, 7, 10, 13, 15, 17, 19, 22, 24};
#elif defined(SINC_POWER)
  static const INT m2K_[] = {0, 0, 2, 5, 8, 11, 12, 14, 16, 18, 21, 23, 24, 24};
#else /* Kaiser-Bessel is the default. */
  static const INT m2K_[] = {1, 3, 7, 9, 14, 17, 20, 23, 24};
#endif

/**
 * Returns an appropriate value of the parameter K used with the PRE_LIN_PSI
 * flag for a given value of the cut-off parameter m.
 */
INT Y(m2K)(const INT m)
{
  int j = MIN(((int)(m)), ((int)((sizeof(m2K_) / sizeof(m2K_[0])) - 1)));
  return (INT)((1U << m2K_[j]) * (m + 2));
}

/**
 * Returns the default window cut off m for the selected window
 */
NFFT_INT Y(get_default_window_cut_off)()
{
  return (NFFT_INT)(WINDOW_HELP_ESTIMATE_m);
}

const char *Y(get_window_name)()
{
  return STRINGIZE(WINDOW_NAME);
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
#else /* Kaiser-Bessel is the default. */
  return NFFT_WINDOW_KAISER_BESSEL;
#endif
}

/* Runtime window evaluation, dispatched on the NFFT_WINDOW_* ordinal: the
 * plan-less counterpart of the PHI/PHI_HUT macros, with the per-axis constants
 * they keep in ths->b formed here from (n, N, m). Kaiser-Bessel carries the
 * macros' 1/I0(m b) scaling, uniform over k and x, so deconvolution and
 * convolution cancel it. Dirac and unrecognized ordinals return 0; callers must
 * decline. */

/* b = pi(2 - N/n), lg_tail = log I0(m b) - m b,
 * i0e_peak_inv = 1/(exp(-m b) I0(m b)), peak_inv = 1/I0(m b). */
static inline void kb_consts(INT n, INT N, int m, R *b, R *lg_tail,
                             R *i0e_peak_inv, R *peak_inv)
{
  const R bb = KPI * (K(2.0) - (R)N / (R)n);
  const R xpk = (R)m * bb;

  *b = bb;
  *lg_tail = Y(bessel_i0_logtail)(xpk);
  *i0e_peak_inv = K(1.0) / Y(bessel_i0_exp_scaled)(xpk);
  *peak_inv = EXP(-xpk - *lg_tail);
}

/* The floor(n x) centring below is the one WINDOW_STENCIL_REACH assumes, so the
 * run reaches m + 1 grid spacings from the node. */
static inline R gaussian_b(INT n, INT N, int m)
{
  const R sigma = (R)n / (R)N;

  return (K(2.0) * sigma) / (K(2.0) * sigma - K(1.0))
    * (((R)m + K(1.0)) / KPI);
}

static inline R sincpow_w(INT n, INT N, int m)
{
  const R sigma = (R)n / (R)N;

  return (K(2.0) * sigma - K(1.0)) / (K(2.0) * (R)m * sigma);
}

R Y(window_phi_hut)(int window, INT n, INT N, int m, INT k)
{
  switch (window)
  {
  case NFFT_WINDOW_KAISER_BESSEL:
  {
    R b, lg_tail, i0e_peak_inv, peak_inv;
    kb_consts(n, N, m, &b, &lg_tail, &i0e_peak_inv, &peak_inv);
    return Y(kb_phi_hut)(b, i0e_peak_inv, (R)m, (R)n, (R)k);
  }
  case NFFT_WINDOW_GAUSSIAN:
  {
    const R b = gaussian_b(n, N, m);
    const R t = KPI * (R)k / (R)n;
    return EXP(-(t * t) * b);
  }
  case NFFT_WINDOW_B_SPLINE:
    return Y(bspline_phi_hut)((R)m, (R)n, (R)k);
  case NFFT_WINDOW_SINC_POWER:
  {
    /* de Boor rather than the macro's Chebyshev evaluator: deconvolution
     * divides by this. */
    const R w = sincpow_w(n, N, m);
    return Y(bsplines)((INT)(2 * m), (R)k / (w * (R)n) + (R)m);
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
    R b, lg_tail, i0e_peak_inv, peak_inv;
    kb_consts(n, N, m, &b, &lg_tail, &i0e_peak_inv, &peak_inv);
    return Y(kb_phi)(b, lg_tail, peak_inv, (R)m, (R)n * x);
  }
  case NFFT_WINDOW_GAUSSIAN:
  {
    const R b = gaussian_b(n, N, m);
    const R nx = (R)n * x;
    return EXP(-(nx * nx) / b) / SQRT(KPI * b);
  }
  case NFFT_WINDOW_B_SPLINE:
    return Y(bsplines)((INT)(2 * m), x * (R)n + (R)m) / (R)n;
  case NFFT_WINDOW_SINC_POWER:
  {
    const R w = sincpow_w(n, N, m);
    return Y(sincpow_phi)(w, (R)m, KPI * (R)n * w * x);
  }
  default:
    return K(0.0);
  }
}

void Y(window_phi_hut_apply)(int window, INT n, INT N, int m, INT k0,
                             R *out, INT count)
{
  INT i;

  if (window == NFFT_WINDOW_KAISER_BESSEL)
  {
    R b, lg_tail, i0e_peak_inv, peak_inv;
    kb_consts(n, N, m, &b, &lg_tail, &i0e_peak_inv, &peak_inv);

    for (i = 0; i < count; i++)
      out[i] = Y(kb_phi_hut)(b, i0e_peak_inv, (R)m, (R)n, (R)(k0 + i));

    return;
  }

  for (i = 0; i < count; i++)
    out[i] = Y(window_phi_hut)(window, n, N, m, k0 + i);
}

void Y(window_phi_precompute)(int window, INT n, INT N, int m,
                              const R *x, INT x_stride, INT num_nodes,
                              R *out, INT out_stride)
{
  INT j;

  if (window == NFFT_WINDOW_KAISER_BESSEL)
  {
    R b, lg_tail, i0e_peak_inv, peak_inv;
    kb_consts(n, N, m, &b, &lg_tail, &i0e_peak_inv, &peak_inv);

    for (j = 0; j < num_nodes; j++)
    {
      const R xj = x[(size_t)j * (size_t)x_stride];
      const INT u = LRINT(FLOOR((R)n * xj)) - (INT)m;

      Y(kb_phi_run)(out + (size_t)j * (size_t)out_stride, b, lg_tail, peak_inv,
          (R)m, (INT)m, (R)n * xj - (R)u);
    }

    return;
  }

  for (j = 0; j < num_nodes; j++)
  {
    const R xj = x[(size_t)j * (size_t)x_stride];
    const INT u = LRINT(FLOOR((R)n * xj)) - (INT)m;
    int l;

    for (l = 0; l <= 2 * m + 1; l++)
      out[(size_t)j * (size_t)out_stride + (size_t)l] =
          Y(window_phi)(window, n, N, m, xj - (R)(u + (INT)l) / (R)n);
  }
}
