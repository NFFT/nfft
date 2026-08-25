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

/* Kaiser-Bessel window. lg_peak = log I0(m b), lg_tail = lg_peak - m b. */

R Y(kb_phi_hut)(R b, R lg_peak, R lg_tail, R m, R n, R k)
{
  const R xpk = m * b;
  const R t = K(2.0) * KPI * k / n;
  const R rs = SQRT(b * b - t * t);
  const R a = m * rs;

  if (a > NFFT_I0_ASYMP_SPLIT && xpk > NFFT_I0_ASYMP_SPLIT)
  {
    /* I0(a)/I0(xpk) with a - xpk rationalized, so small t does not cancel */
    const R da = -m * t * t / (rs + b);
    return EXP(da + Y(bessel_i0_logtail)(a) - lg_tail);
  }

  return Y(bessel_i0_scaled)(a, lg_peak);
}

R Y(kb_phi)(R b, R lg_peak, R lg_tail, R m, R n, R x)
{
  const R nx = n * x;
  const R a = m * m - nx * nx;

  if (a > K(0.0))
  {
    const R ra = SQRT(a);
    const R u = b * ra;
    /* sinh(u) exp(-lg_peak) as exp(u - m b - lg_tail) (1 - exp(-2u)) / 2, with
     * u - m b rationalized: no overflow for large u, none of the cancellation
     * a difference of exponentials has for small u */
    const R du = -b * nx * nx / (ra + m);
    return K(-0.5) * EXP(du - lg_tail) * EXPM1(K(-2.0) * u) / (KPI * ra);
  }

  if (a < K(0.0))
  {
    const R rma = SQRT(-a);
    return SIN(b * rma) * EXP(-lg_peak) / (KPI * rma);
  }

  return (b / KPI) * EXP(-lg_peak);
}
