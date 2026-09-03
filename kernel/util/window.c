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

/* Half-width the Gaussian is fitted to: the stencil reach plus the measured
 * correction, withdrawn before the transform saturates, since past that point
 * a wider window only raises the round-off floor. At s = 0 the width is the
 * bare reach exactly, so saturated cases keep their old error.
 *
 * Measured in binary128, one axis of the corrected window leaves about
 * 0.70 exp(-6.54 (sigma - 1) h / (2 sigma - 1)). The round-off floor is near
 * 10 epsilon on one axis and grows with the deconvolution amplification
 * lambda (rms of 1 / phi_hat over the band) as lambda^((d - 1) / 2): measured
 * 4, 25, 186 times the one-axis floor for d = 2, 3, 4 at sigma = 2, against
 * 5, 25, 125 predicted. Full correction at 20 floors up, none at 2; the 2
 * covers that spread. */
R Y(gaussian_half_width)(const R m, const R sigma, const INT d, const R reach,
    const R corr)
{
  const R h = m + reach + corr;
  const R b = K(2.0) * sigma * h / ((K(2.0) * sigma - K(1.0)) * KPI);
  const R err = K(0.70)
      * EXP(-K(6.54) * (sigma - K(1.0)) * h / (K(2.0) * sigma - K(1.0)));
  const R lambda = EXP(b * KPI * KPI / (K(4.0) * sigma * sigma))
      * sigma / (KPI * SQRT(b));
  const R floor_d = K(10.0) * EPSILON * POW(lambda, ((R)(d - 1)) / K(2.0));
  R s = LOG10(err / (K(2.0) * floor_d));

  s = s < K(0.0) ? K(0.0) : (s > K(1.0) ? K(1.0) : s);
  return m + reach + s * corr;
}
