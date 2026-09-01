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

/* Modified Bessel function I0, in the two scaled forms the Kaiser-Bessel window
 * needs. Two ranges, split at NFFT_I0_ASYMP_SPLIT:
 *
 *   i0_small_m1(x) = I0(x) - 1                 for x <= split
 *   i0_asymp(u)    = sqrt(x) * exp(-x) * I0(x) for x >  split, u = 1/x
 *
 * i0_small_m1 returns I0 - 1, not I0, so bessel_i0_logtail can use LOG1P to
 * calculate log I0(x) ~ x^2/4 as x -> 0.
 *
 * I0 itself is not offered. Above the split it would be exp(x) times the
 * asymptotic form, and any such form inherits the relative error of EXP, which
 * on some targets grows with the argument rather than staying near an ulp: an
 * x87 exp expanded inline reduces x * log2(e) in one rounding, which costs
 * 0.5 * eps * x. Both forms below avoid that, one by never exponentiating above
 * the split and the other by taking a logarithm instead. */

static const INT N1 = sizeof(NFFT_I0_P1) / sizeof(NFFT_I0_P1[0]);
static const INT N2 = sizeof(NFFT_I0_P2) / sizeof(NFFT_I0_P2[0]);

/* Both tables are summed as four independent chains over the coefficients whose
 * index shares a residue mod 4, so the dependent-operation count is a quarter of
 * a Horner pass. Table lengths are multiples of four, which is what lets the
 * loop run without a prologue; tests/besselgen pins the degrees to keep that.
 *
 * The regrouping is safe because neither sum cancels: every P1 coefficient is
 * positive, and the split is chosen per format so that P2 has a growth factor
 * of 1 (branch2_growth in tests/besselgen/scheme.py). */
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

static inline R i0_small_m1(const R x)
{
  const R h = x * K(0.5);
  const R y = h * h;

  return y * poly4(NFFT_I0_P1, N1, y);
}

/* Monomial in u = 1/x, which lies in (0, 1/split]. Callers pass u rather than x
 * because they need it for SQRT(u) below as well. */
static inline R i0_asymp(const R u)
{
  return poly4(NFFT_I0_P2, N2, u);
}

/* exp(-x) * I0(x), in (0, 1]. A ratio of two of these needs no logarithm and
 * cannot overflow: I0(a)/I0(c) is exp(a-c) times the ratio, and the caller
 * forms a-c however its own geometry allows. Above the split the exponential
 * cancels against the asymptotic form and disappears. */
R Y(bessel_i0_exp_scaled)(R x)
{
  x = FABS(x);

  /* The only EXP left in this file. Its argument is bounded by the split, so
   * the growth noted above is bounded with it. */
  if (x <= NFFT_I0_ASYMP_SPLIT)
    return (K(1.0) + i0_small_m1(x)) * EXP(-x);

  {
    const R u = K(1.0) / x;
    return SQRT(u) * i0_asymp(u);
  }
}

R Y(bessel_i0_logtail)(R x)
{
  x = FABS(x);

  if (x <= NFFT_I0_ASYMP_SPLIT)
    /* Subtracts like-sized quantities above x ~ 1, losing a little over three
     * bits at the splits in use. */
    return LOG1P(i0_small_m1(x)) - x;

  {
    const R u = K(1.0) / x;
    return LOG(SQRT(u) * i0_asymp(u));
  }
}
