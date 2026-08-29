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
