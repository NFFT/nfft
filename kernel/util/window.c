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

/* Half-width the Gaussian is fitted to: the argmin of a model of the transform
 * error over the offset c in h = m + c, relative to the stencil reach. Below
 * saturation the optimum sits about 0.3 past the reach;
 * once round-off dominates it moves below the reach, since a narrower window
 * lowers the amplified floor. The model, one axis, relative L2, is
 *
 *   sqrt(d) [erfc(H / sqrt b) lambda(b) + alias(b)]
 *     + 10 eps (lambda(b0) / 4.6)^(d-1) (lambda(b) / lambda(b0))^(1.3 (d-1))
 *
 * with H the reach, b0 the shape parameter at the reach, lambda the rms of
 * 1 / phi_hat over the band and alias the rms folded-copy sum, both in closed
 * form. The round-off constants are measured: a one-axis floor of 10 eps
 * (9.4 and 11.9 seen), a per-dimension level of lambda / 4.6 (within 1.7x for
 * d <= 4 at sigma 2 and 2.286) and a slope of 1.3 per extra dimension. The
 * argmin matches swept optima to 0.12 rms over float, double and binary128,
 * d = 1..5; the closed forms shift it by at most 0.14. The search is bounded
 * to the range that was measured. */
static R gaussian_model(const R c, const R m, const R sigma, const INT d,
    const R reach, const R lam0)
{
  const R h = m + c, H = m + reach;
  const R b = K(2.0) * sigma * h / ((K(2.0) * sigma - K(1.0)) * KPI);
  const R lam = EXP(b * KPI * KPI / (K(4.0) * sigma * sigma))
      * sigma / (KPI * SQRT(b));
  const R alias = K(2.0) * EXP(-b * KPI * KPI * (K(1.0) - K(1.0) / sigma))
      * SQRT(sigma / (K(8.0) * b * KPI * KPI));
  const R x = H / SQRT(b);
  const R tail = EXP(-x * x) / (x * SQRT(KPI));
  const R dm1 = (R)(d - 1);

  return SQRT((R)d) * (tail * lam + alias)
      + K(10.0) * EPSILON * POW(lam0 / K(4.6), dm1)
          * POW(lam / lam0, K(1.3) * dm1);
}

R Y(gaussian_half_width)(const R m, const R sigma, const INT d, const R reach,
    const R corr)
{
  const R b0 = K(2.0) * sigma * (m + reach) / ((K(2.0) * sigma - K(1.0)) * KPI);
  const R lam0 = EXP(b0 * KPI * KPI / (K(4.0) * sigma * sigma))
      * sigma / (KPI * SQRT(b0));
  R lo = reach - K(1.5), hi = reach + corr + K(0.3);
  int i;

  for (i = 0; i < 40; i++)
  {
    const R a = lo + (hi - lo) / K(3.0), b = hi - (hi - lo) / K(3.0);

    if (gaussian_model(a, m, sigma, d, reach, lam0)
        < gaussian_model(b, m, sigma, d, reach, lam0))
      hi = b;
    else
      lo = a;
  }
  /* Leave the reach only for a clear predicted gain, so that where the model
   * is flat or marginal, at the floor, the parent's width is kept exactly
   * rather than one ulp of noise traded for another. Widening needs 15%: the
   * model puts the analytic optimum in the right place but at about 0.6 of
   * its measured depth, so 15% predicted is the 1.3x seen. Narrowing needs
   * 20%, above the 1.1x it buys in two dimensions, where the slope term is
   * least certain. */
  {
    const R c = (lo + hi) / K(2.0);
    const R margin = c > reach ? K(0.85) : K(0.80);

    if (gaussian_model(c, m, sigma, d, reach, lam0)
        > margin * gaussian_model(reach, m, sigma, d, reach, lam0))
      return m + reach;
    return m + c;
  }
}
