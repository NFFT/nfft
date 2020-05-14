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

R Y(float_property)(const float_property p)
{
  const R base = FLT_RADIX;
  static R eps = K(1.0);
  const R t = MANT_DIG;
  const R emin = MIN_EXP;
  const R emax = MAX_EXP;
  const R prec = eps * base;
  static R rmin = K(1.0);
  static R rmax = K(1.0);
  const R rnd = FLTROUND;
  static R sfmin = K(-1.0);
  static short first = TRUE;

  if (first)
  {
    /* Compute eps = 2^(1-MANT_DIG).
     * The usual definition of EPSILON is too small for double-double arithmetic on PowerPC. */
    for (INT i=0; i<MANT_DIG-1; i++)
      eps /= K(2.0);

    /* Compute rmin */
    {
      const INT n = 1 - MIN_EXP;
      INT i;
      for (i = 0; i < n; i++)
        rmin /= base;
    }

    /* Compute rmax */
    {
      INT i;
      rmax -= eps;
      for (i = 0; i < emax; i++)
        rmax *= base;
    }

    /* Compute sfmin */
    {
      R small = K(1.0) / rmax;
      sfmin = rmin;
      if (small >= sfmin)
        sfmin = small * (eps + K(1.0));
    }

    first = FALSE;
  }

  if (p == NFFT_EPSILON)
    return eps;
  else if (p == NFFT_SAFE__MIN)
    return sfmin;
  else if (p == NFFT_BASE)
    return base;
  else if (p == NFFT_PRECISION)
    return prec;
  else if (p == NFFT_MANT_DIG)
    return t;
  else if (p == NFFT_FLTROUND)
    return rnd;
  else if (p == NFFT_E_MIN)
    return  emin;
  else if (p == NFFT_R_MIN)
    return rmin;
  else if (p == NFFT_E_MAX)
    return emax;
  else if (p == NFFT_R_MAX)
    return rmax;
  else
    CK(0 /* cannot happen */);

  return K(-1.0);
} /* dlamch_ */

/** Computes double /f$\prod_{t=0}^{d-1} v_t/f$. */
R Y(prod_real)(R *vec, INT d)
{
  INT t;
  R prod;

  prod = K(1.0);
  for (t = 0; t < d; t++)
    prod *= vec[t];

  return prod;
}
