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

/* Polynomial form of the Kaiser-Bessel window.
 *
 * One polynomial per unit cell of the window argument, in the offset t in
 * [0, 1) into the cell. A run of 2m + 2 taps spans 2m + 2 cells at one t, so
 * it is one Horner sweep over 2m + 2 columns, in place of 2m + 2 calls to
 * EXP, EXPM1 and SIN. Each column is analytic in t: both branches of the
 * window are even in the radical and ra^2 is a polynomial in t, so there is no
 * branch point and a Chebyshev interpolant converges geometrically.
 *
 * The interpolant is built from the window itself, not from a rewrite of it:
 * a plain sinh(b ra)/ra ratio carries 20 times the error of Y(kb_phi) and that
 * error would cap the fit.
 */

#include "nfft3.h"
#include "infft.h"
#include "kbpoly_data.h"

/* Chebyshev coefficients of the values in val, into dst, one column per cell.
 * dst[0] comes out already halved, so the series is sum_j dst[j] T_j. */
static void cheb_of(R *dst, const R *val, const R *cm, const INT w,
    const INT deg)
{
  const INT n = deg + 1;
  INT j, k, l;

  for (j = 0; j <= deg; j++)
    for (l = 0; l < w; l++)
    {
      R s = K(0.0);

      for (k = 0; k < n; k++)
        s += val[k * w + l] * cm[j * n + k];

      dst[j * w + l] = (K(2.0) / (R)n) * s * IF(j == 0, K(0.5), K(1.0));
    }
}

/* T_j(2t - 1) in powers of t, so the run time can use plain Horner. */
static void cheb_to_monomial(R *coef, const R *cheb, R *tp, const INT w,
    const INT deg)
{
  const INT n = deg + 1;
  INT i, j, l;

  for (i = 0; i < n * n; i++)
    tp[i] = K(0.0);

  tp[0] = K(1.0);

  if (deg >= 1)
  {
    tp[n + 0] = K(-1.0);
    tp[n + 1] = K(2.0);
  }

  for (j = 2; j <= deg; j++)
    for (i = 0; i <= j; i++)
    {
      R v = -tp[(j - 2) * n + i] - K(2.0) * tp[(j - 1) * n + i];

      if (i > 0)
        v += K(4.0) * tp[(j - 1) * n + i - 1];

      tp[j * n + i] = v;
    }

  for (i = 0; i <= deg; i++)
    for (l = 0; l < w; l++)
    {
      R s = K(0.0);

      for (j = i; j <= deg; j++)
        s += cheb[j * w + l] * tp[j * n + i];

      coef[i * w + l] = s;
    }
}

INT Y(kb_poly_degree)(const R b, const R reach, const INT m)
{
  const R sigma = KPI / (K2PI - b);
  INT mt = (INT)CEIL(KB_POLY_BR_OPT / b) - 1, deg, k;

  A(reach == (R)m + K(1.0));
  (void)reach;

  /* Past m_opt the degree does not rise. */
  mt = IF(m < mt, m, mt);
  mt = IF(mt < KB_POLY_TAB_MMAX, mt, KB_POLY_TAB_MMAX);
  mt = IF(mt > 1, mt, 1);
  deg = kb_poly_tab_deg0[mt];

  /* Thresholds are rounded down. A sigma off [1.1, 4], which X(check)
   * rejects, gets the degree at the nearer end. */
  for (k = kb_poly_tab_off[mt]; k < kb_poly_tab_off[mt + 1]; k++)
    deg += IF(sigma >= K(1E-4) * (R)kb_poly_tab_brk[k], 1, 0);

  return IF(deg > KB_POLY_DEG_MAX, KB_POLY_DEG_MAX, deg);
}

void Y(kb_poly_fit)(R *coef, const R b, const R lg_tail, const R peak_inv,
    const R reach, const INT m, const INT deg)
{
  const INT w = KB_POLY_COLS(m), n = deg + 1;
  /* val, cheb, res are n by w; cm and tp are n by n; node is n; out is w. */
  R *mem = (R*) Y(malloc)((size_t)(3 * n * w + 2 * n * n + n + w) * sizeof(R));
  R *val = mem, *cheb = val + n * w, *res = cheb + n * w;
  R *cm = res + n * w, *tp = cm + n * n, *node = tp + n * n, *out = node + n;
  INT j, k, l;

  for (k = 0; k < n; k++)
    node[k] = (COS(KPI * ((R)k + K(0.5)) / (R)n) + K(1.0)) / K(2.0);

  for (j = 0; j <= deg; j++)
    for (k = 0; k < n; k++)
      cm[j * n + k] = COS(KPI * (R)j * ((R)k + K(0.5)) / (R)n);

  for (k = 0; k < n; k++)
    for (l = 0; l < w; l++)
      val[k * w + l] = Y(kb_phi)(b, lg_tail, peak_inv, reach,
          node[k] + (R)(m + 1 - l));

  cheb_of(cheb, val, cm, w, deg);
  cheb_to_monomial(coef, cheb, tp, w, deg);

  /* One refinement step, against the evaluator the run time will use, so the
   * correction absorbs the conversion's rounding as well as the fit's. The
   * residual is a difference of two values that already agree to a few units
   * in the last place, so it is formed without cancellation error and its own
   * fit costs eps of something tiny. A second step changes nothing. */
  for (k = 0; k < n; k++)
  {
    for (l = 0; l < w; l++)
      out[l] = coef[deg * w + l];

    for (j = deg - 1; j >= 0; j--)
      for (l = 0; l < w; l++)
        out[l] = out[l] * node[k] + coef[j * w + l];

    for (l = 0; l < w; l++)
      res[k * w + l] = val[k * w + l] - out[l];
  }

  cheb_of(cheb, res, cm, w, deg);
  cheb_to_monomial(res, cheb, tp, w, deg);

  for (j = 0; j <= deg; j++)
    for (l = 0; l < w; l++)
      coef[j * w + l] += res[j * w + l];

  Y(free)(mem);
}

/* The first axis with the same b as axis t, t itself if none. The other
 * window constants follow from b and the reach, so equal b means an equal
 * table. */
static INT first_alike(const R *b, const INT t)
{
  INT s;

  for (s = 0; s < t && b[s] != b[t]; s++)
    ;

  return s;
}

R *Y(kb_poly_init)(const R *b, const INT d, const INT m, const R reach)
{
  INT deg = 0, stride, t;
  R *tab;

  /* One degree for all axes, the one the hardest axis needs. */
  for (t = 0; t < d; t++)
  {
    const INT dt = Y(kb_poly_degree)(b[t], reach, m);

    deg = IF(dt > deg, dt, deg);
  }

  stride = (deg + 1) * KB_POLY_COLS(m);
  tab = (R*) Y(malloc)((size_t)(1 + stride * d) * sizeof(R));
  tab[0] = (R)deg;

  for (t = 0; t < d; t++)
  {
    const INT s = first_alike(b, t);

    if (s == t)
      Y(kb_poly_fit)(tab + 1 + t * stride, b[t], b[d + t], b[3 * d + t],
          reach, m, deg);
    else
      memcpy(tab + 1 + t * stride, tab + 1 + s * stride,
          (size_t)stride * sizeof(R));
  }

  return tab;
}
