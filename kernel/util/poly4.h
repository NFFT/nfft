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

#ifndef NFFT_POLY4_H
#define NFFT_POLY4_H

/* Evaluates the polynomial c[0] + c[1]*u + ... + c[n-1]*u^(n-1).
 *
 * Uses for disjoint chain to enable vectorization comapred to standard Horner.
 *
 * The terms are the same, only the summation order differs, so the regrouping
 * is safe only where the sum cannot cancel. The table length must be a multiple
 * of four. */
static inline R poly4(const R *c, const INT n, const R u)
{
  const R v = u * u;
  const R w = v * v;
  R a0, a1, a2, a3;
  INT j = n - 4;

  A(n >= 4 && n % 4 == 0);

  a0 = c[j];
  a1 = c[j + 1];
  a2 = c[j + 2];
  a3 = c[j + 3];

  for (j -= 4; j >= 0; j -= 4)
  {
    a0 = a0 * w + c[j];
    a1 = a1 * w + c[j + 1];
    a2 = a2 * w + c[j + 2];
    a3 = a3 * w + c[j + 3];
  }

  return (a0 + u * a1) + v * (a2 + u * a3);
}

#endif /* NFFT_POLY4_H */
