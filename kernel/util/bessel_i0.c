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
#include "bessel_i0_data.h"

/* Modified Bessel function I0 and its logarithmic forms. Two ranges, split at
 * NFFT_I0_ASYMP_SPLIT:
 *
 *   i0_small_m1(x) = I0(x) - 1                 for x <= split
 *   i0_asymp(x)    = sqrt(x) * exp(-x) * I0(x) for x >  split
 *
 * i0_small_m1 returns I0 - 1, not I0, so bessel_i0_log can use LOG1P to calculate
 * log I0(x) ~ x^2/4 as x -> 0. */

static const INT N1 = sizeof(NFFT_I0_P1) / sizeof(NFFT_I0_P1[0]);
static const INT N2 = sizeof(NFFT_I0_P2) / sizeof(NFFT_I0_P2[0]);

/* Every coefficient is positive, so the Horner chain cannot cancel. */
static inline R i0_small_m1(const R x)
{
  const R h = x * K(0.5);
  const R y = h * h;
  R r = NFFT_I0_P1[N1 - 1];
  INT j;

  A(N1 >= 2);

  /* Horner */
  for (j = N1 - 2; j >= 0; j--)
    r = r * y + NFFT_I0_P1[j];

  return y * r;
}

/* t maps (split, inf) onto [-1, 1), keeping Clenshaw in its stable domain. */
static inline R i0_asymp(const R x)
{
  const R t = (K(2.0) * NFFT_I0_ASYMP_SPLIT - x) / x;
  const R t2 = K(2.0) * t;
  R b1 = K(0.0), b2 = K(0.0), b0;
  INT j;

  A(N2 >= 2);

  /* Clenshaw */
  for (j = N2 - 1; j >= 1; j--)
  {
    b0 = t2 * b1 - b2 + NFFT_I0_P2[j];
    b2 = b1;
    b1 = b0;
  }

  return t * b1 - b2 + NFFT_I0_P2[0];
}

/* exp(x)/sqrt(x) * p. The 1/sqrt(x) and p < 1 keep I0 finite past the point
 * where exp(x) alone overflows. */
static inline R i0_from_asymp(const R x, const R p)
{
  if (x <= NFFT_LOG_MAX_SAFE)
    return (EXP(x) / SQRT(x)) * p;

  /* SPlit exp(x) into exp(x/2) and exp(x/2) to avoid overflow. */
  {
    const R e = EXP(x * K(0.5));
    return (e / SQRT(x)) * (e * p);
  }
}

R Y(bessel_i0)(R x)
{
  x = FABS(x);

  if (x <= NFFT_I0_ASYMP_SPLIT)
    return K(1.0) + i0_small_m1(x);

  return i0_from_asymp(x, i0_asymp(x));
}

R Y(bessel_i0_log)(R x)
{
  x = FABS(x);

  if (x <= NFFT_I0_ASYMP_SPLIT)
    return LOG1P(i0_small_m1(x));

  return x - K(0.5) * LOG(x) + LOG(i0_asymp(x));
}

R Y(bessel_i0_scaled)(R x, R lg_peak)
{
  x = FABS(x);

  if (x <= NFFT_I0_ASYMP_SPLIT)
    return (K(1.0) + i0_small_m1(x)) * EXP(-lg_peak);

  /* x - lg_peak cancels for a large peak, and the exponential turns the loss
   * into relative error of up to lg_peak*eps. */
  return (EXP(x - lg_peak) / SQRT(x)) * i0_asymp(x);
}

R Y(bessel_i0_logtail)(R x)
{
  x = FABS(x);

  if (x <= NFFT_I0_ASYMP_SPLIT)
    /* Subtracts like-sized quantities above x ~ 1, losing a little over three
     * bits at the splits in use. */
    return LOG1P(i0_small_m1(x)) - x;

  return LOG(i0_asymp(x)) - K(0.5) * LOG(x);
}
