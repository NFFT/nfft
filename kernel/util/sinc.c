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

#include "infft.h"
#include "sinc_data.h"

R Y(sinc)(const R x)
{
  /* Based on sinc function from Boost C++ library. */
  const R b = EPSILON;
  const R bs = SQRT(b);
  const R bs2 = SQRT(bs);

  if (FABS(x) >= bs2)
    return SIN(x)/x;
  else
  {
    R r = K(1.0);

    if (FABS(x) >= b)
    {
      const R x2 = x * x;
      r -= x2 / K(6.0);

      if (FABS(x) >= bs)
        r += (x2 * x2) / K(120.0);
    }

    return r;
  }
}

/* log|sinc(x)|, for callers that form EXP(2m * log_sinc(x)) instead of
 * POW(sinc(x), 2m). Absolute value, not log(sinc): 2m is even, so the callers
 * need a value past the first zero of sinc, which PHI for the sinc-power
 * window reaches at large sigma. Coefficients and their rationale in
 * tests/sincgen. */
R Y(log_sinc)(const R x)
{
  const R a = FABS(x);

  if (a <= NFFT_LOG_SINC_SPLIT)
  {
    const R y = a * a;
    R q = NFFT_LOG_SINC_Q[SIZE(NFFT_LOG_SINC_Q) - 1];
    INT j;

    for (j = (INT)SIZE(NFFT_LOG_SINC_Q) - 2; j >= 0; j--)
      q = q * y + NFFT_LOG_SINC_Q[j];

    return y * q;
  }

  return LOG(FABS(SIN(a) / a));
}
