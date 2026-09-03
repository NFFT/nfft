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

static inline void bspline_help(const INT k, const R x, R *scratch, const INT j,
  const INT ug, const INT og, const INT r)
{
  INT i; /* Row index of the de Boor scheme */
  INT idx; /* Index in scratch */
  R a; /* Alpha of de Boor scheme */

  /* computation of one column */
  for (i = og + r - k + 1, idx = og; idx >= ug; i--, idx--)
  {
    a = (x - (R)i) / ((R)(k - j));
    scratch[idx] = (K(1.0) - a) * scratch[idx - 1] + a * scratch[idx];
  }
}

/* Evaluate the cardinal B-Spline B_{n-1} supported on [0,n]. */
R Y(bsplines)(const INT k, const R _x)
{
  const R kk = (R)k;
  R result_value;
  INT r;
  INT g1, g2; /* boundaries */
  INT j, idx, ug, og; /* indices */
  R a; /* Alpha of de Boor scheme*/
  R x = _x;
  R scratch[k];

  result_value = K(0.0);

  if (K(0.0) < x && x < kk)
  {
    /* Exploit symmetry around k/2, maybe. */
    if ( (kk - x) < x)
    {
      x = kk - x;
    }

    r = (INT)LRINT(CEIL(x) - K(1.0));

    /* Do not use the explicit formula x^k / k! for first interval! De Boor's
     * algorithm is more accurate. See https://github.com/NFFT/nfft/issues/16.
     */

    for (idx = 0; idx < k; idx++)
      scratch[idx] = K(0.0);

    scratch[k-r-1] = K(1.0);

    /* Bounds of the algorithm. */
    g1 = r;
    g2 = k - 1 - r;
    ug = g2;

    /* g1 <= g2 */

    for (j = 1, og = g2 + 1; j <= g1; j++, og++)
    {
      a = (x + (R)(k - r - og - 1)) / ((R)(k - j));
      scratch[og] = (K(1.0) - a) * scratch[og-1];
      bspline_help(k, x, scratch, j, ug + 1, og - 1, r);
      a = (x + (R)(k - r - ug - 1)) / ((R)(k - j));
      scratch[ug] = a * scratch[ug];
    }

    for (og-- ; j <= g2; j++)
    {
      bspline_help(k, x, scratch, j, ug + 1, og, r);
      a = (x + (R)(k - r - ug - 1)) / ((R)(k - j));
      scratch[ug] = a * scratch[ug];
    }

    for(; j < k; j++)
    {
      ug++;
      bspline_help(k, x, scratch, j, ug, og, r);
    }

    result_value = scratch[k-1];
  }

  return(result_value);
}

/* Coefficient l of s * P in the Chebyshev basis, from
 * s T_l = (T_{l+1} + T_{l-1})/2. P has degree below deg, and the rows of the
 * table are packed with no slack, so coefficient deg must be read as zero
 * rather than off the end of the row. */
static inline R cheb_shift(const R *c, const INT l, const INT deg)
{
  const R hi = (l + 1 < deg) ? c[l + 1] : K(0.0);

  if (l == 0)
    return K(0.5) * hi;

  if (l == 1)
    return c[0] + K(0.5) * hi;

  return K(0.5) * (c[l - 1] + hi);
}

/* Chebyshev coefficients of every piece of B_k: row i of tab is the series of
 * B_k on [i, i+1) in s = 2(x - i) - 1, k rows of k reals.
 *
 * Built from B_1, the box, by lifting the order recurrence
 *
 *   B_k(x) = x/(k-1) B_{k-1}(x) + (k-x)/(k-1) B_{k-1}(x-1)
 *
 * one order at a time. On [i, i+1) both terms carry the same local s, because
 * the shift by one interval and the shift by one knot cancel, so a lift is a
 * multiplication by two linear polynomials and stays inside the basis. Nothing
 * here evaluates a B-spline, so the coefficients carry no evaluator's error.
 *
 * The rows are lifted from the top down, which lets a lift read the order k-1
 * rows it needs out of the rows it has not reached yet. */
void Y(bspline_cheb_init)(R *tab, const INT k)
{
  INT kk, i, l;

  for (i = 0; i < k * k; i++)
    tab[i] = K(0.0);

  tab[0] = K(1.0);

  for (kk = 2; kk <= k; kk++)
  {
    const R inv = K(1.0) / (R)(kk - 1);
    R row[kk];

    for (i = kk - 1; i >= 0; i--)
    {
      const R a = ((R)i + K(0.5)) * inv, b = K(0.5) * inv;
      const R g = ((R)kk - (R)i - K(0.5)) * inv, d = -K(0.5) * inv;
      const R *p = (i < kk - 1) ? tab + i * k : 0;       /* B_{kk-1}(x)   */
      const R *q = (i > 0) ? tab + (i - 1) * k : 0;      /* B_{kk-1}(x-1) */

      for (l = 0; l < kk; l++)
      {
        R v = K(0.0);

        if (p)
          v += a * p[l] + b * cheb_shift(p, l, kk);
        if (q)
          v += g * q[l] + d * cheb_shift(q, l, kk);

        row[l] = v;
      }

      for (l = 0; l < kk; l++)
        tab[i * k + l] = row[l];
    }
  }
}

/* The first piece on which Clenshaw keeps relative accuracy: the pieces are
 * admissible from i0 through k-1-i0, by symmetry of B_k about k/2.
 *
 * Clenshaw's error is bounded by the sum of the coefficients rather than by the
 * value, so a piece is admissible when that sum stays within thresh of the
 * smallest value on it. B_k is unimodal, so that smallest value sits at a knot.
 * The outer pieces fail by a wide margin -- B_k reaches zero -- and callers that
 * only need accuracy against the peak have no use for this bound. */
INT Y(bspline_cheb_guard)(const R *tab, const INT k, const R thresh)
{
  INT i, l, i0 = 0;

  for (i = 0; i < k / 2; i++)
  {
    const R e0 = Y(bsplines)(k, (R)i), e1 = Y(bsplines)(k, (R)(i + 1));
    const R lo = MIN(e0, e1);
    R sum = K(0.0);

    for (l = 0; l < k; l++)
      sum += FABS(tab[i * k + l]);

    if (lo <= K(0.0) || sum > thresh * lo)
      i0 = i + 1;
  }

  return i0;
}
