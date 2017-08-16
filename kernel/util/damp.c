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

/**
 * Compute damping factor for modified Fejer kernel:
 * /f$\frac{2}{N}\left(1-\frac{\left|2k+1\right|}{N}\right)/f$
 */
R Y(modified_fejer)(const INT N, const INT kk)
{
  return (K(2.0) / ((R) (N * N))
      * (K(1.0) - FABS(K(2.0) * ((R) kk) + K(1.0) ) / ((R) N)));
}

/** Compute damping factor for modified Jackson kernel. */
R Y(modified_jackson2)(const INT N, const INT kk)
{
  INT kj;
  const R n = (((R) N) / K(2.0) + K(1.0) ) / K(2.0);
  R result, k;

  for (result = K(0.0), kj = kk; kj <= kk + 1; kj++)
  {
    k = (R)(ABS(kj));

    if (k / n < K(1.0) )
      result += K(1.0)
          - (K(3.0) * k + K(6.0) * n * POW(k, K(2.0) )
              - K(3.0) * POW(k, K(3.0) ))
              / (K(2.0) * n * (K(2.0) * POW(n, K(2.0) ) + K(1.0) ));
    else
      result += (K(2.0) * n - k) * (POW(2 * n - k, K(2.0) ) - K(1.0) )
          / (K(2.0) * n * (K(2.0) * POW(n, K(2.0) ) + K(1.0) ));
  }

  return result;
}

/** Compute damping factor for modified generalised Jackson kernel. */
R Y(modified_jackson4)(const INT N, const INT kk)
{
  INT kj;
  const R n = (((R) N) / K(2.0) + K(3.0) ) / K(4.0);
  const R normalisation = (K(2416.0) * POW(n, K(7.0) )
      + K(1120.0) * POW(n, K(5.0) ) + K(784.0) * POW(n, K(3.0) ) + K(720.0) * n);
  R result, k;

  for (result = K(0.0), kj = kk; kj <= kk + 1; kj++)
  {
    k = (R)(ABS(kj));

    if (k / n < K(1.0) )
      result += K(1.0)
          - (K(1260.0) * k
              + (K(1680.0) * POW(n, K(5.0) ) + K(2240.0) * POW(n, K(3.0) )
                  + K(2940.0) * n) * POW(k, K(2.0) )
              - K(1715.0) * POW(k, K(3.0) )
              - (K(560.0) * POW(n, K(3.0) ) + K(1400.0) * n) * POW(k, K(4.0) )
              + K(490.0) * POW(k, K(5.0) ) + K(140.0) * n * POW(k, K(6.0) )
              - K(35.0) * POW(k, K(7.0) )) / normalisation;

    if ((K(1.0) <= k / n) && (k / n < K(2.0) ))
      result += ((K(2472.0) * POW(n, K(7.0) ) + K(336.0) * POW(n, K(5.0) )
          + K(3528.0) * POW(n, K(3.0) ) - K(1296.0) * n)
          - (K(392.0) * POW(n, K(6.0) ) - K(3920.0) * POW(n, K(4.0) )
              + K(8232.0) * POW(n, K(2.0) ) - K(756.0) )*k
          - (K(504.0)*POW(n, K(5.0)) + K(10080.0)*POW(n, K(3.0))
              - K(5292.0)*n)*POW(k, K(2.0)) - (K(1960.0)*POW(n, K(4.0))
              - K(7840.0)*POW(n, K(2.0)) + K(1029.0))*POW(k, K(3.0))
          + (K(2520.0)*POW(n, K(3.0)) - K(2520.0)*n) * POW(k, K(4.0))
          - (K(1176.0)*POW(n, K(2.0)) - K(294.0)) * POW(k, K(5.0))
          + K(252.0)*n*POW(k, K(6.0)) - K(21.0)*POW(k, K(7.0)))/normalisation;

    if ((K(2.0) <= k / n) && (k / n < K(3.0) ))
      result += (-(K(1112.0) * POW(n, K(7.0) ) - K(12880.0) * POW(n, K(5.0) )
          + K(7448.0) * POW(n, K(3.0) ) - K(720.0) * n)
          + (K(12152.0) * POW(n, K(6.0) ) - K(27440.0) * POW(n, K(4.0) )
              + K(8232.0) * POW(n, K(2.0) ) - K(252.0) )*k
          - (K(19320.0)*POW(n, K(5.0)) - K(21280.0)*POW(n, K(3.0))
              + K(2940.0)*n)*POW(k, K(2.0)) + (K(13720.0)*POW(n, K(4.0))
              - K(7840.0)*POW(n, K(2.0)) + K(343.0))*POW(k, K(3.0))
          - (K(5320.0)*POW(n, K(3.0)) - K(1400.0)*n)*POW(k, K(4.0))
          + (K(1176.0)*POW(n, K(2.0)) - K(98.0))*POW(k, K(5.0))
          - K(140.0)*n*POW(k, K(6.0)) + K(7.0) * POW(k, K(7.0)))/normalisation;

    if ((K(3.0) <= k / n) && (k / n < K(4.0) ))
      result += ((4 * n - k)
          * (POW(4 * n - k, K(2.0) ) - K(1.0) )*(POW(4*n-k, K(2.0))
              - K(4.0))*(POW(4*n-k, K(2.0)) - K(9.0)))/normalisation;
        }

  return result;
}

/** Compute damping factor for modified Sobolev kernel. */
R Y(modified_sobolev)(const R mu, const INT kk)
{
  R result;
  INT kj, k;

  for (result = K(0.0), kj = kk; kj <= kk + 1; kj++)
  {
    k = ABS(kj);
    if (k == 0)
      result += K(1.0);
    else
      result += POW((R) k, -K(2.0) * mu);
  }

  return result;
}

/** Comput damping factor for modified multiquadric kernel. */
R Y(modified_multiquadric)(const R mu, const R c, const INT kk)
{
  R result;
  INT kj, k;

  for (result = K(0.0), kj = kk; kj <= kk + 1; kj++)
  {
    k = ABS(kj);
    result += POW((R)(k * k) + c * c, -mu);
  }

  return result;
}
