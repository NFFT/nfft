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

/* Polynomial form of the Kaiser-Bessel and sinc-power windows.
 *
 * One polynomial per unit cell of the window argument, in the offset t in
 * [0, 1) into the cell. A run of 2m + 2 taps spans 2m + 2 cells at one t, so
 * it is one Horner sweep over 2m + 2 columns, in place of 2m + 2 calls to the
 * window's closed form. Each column is analytic in t, so a Chebyshev
 * interpolant converges geometrically:
 *
 *   Kaiser-Bessel: both branches of the window are even in the radical and
 *     ra^2 is a polynomial in t, so there is no branch point.
 *   sinc-power: w sinc(pi w nx)^(2m) is entire. Its zeros are of order 2m
 *     and cost the interpolant nothing.
 *
 * The interpolant is built from the window itself, not from a rewrite of it:
 * a plain sinh(b ra)/ra ratio carries 20 times the error of Y(kb_phi) and that
 * error would cap the fit. The sinc power's closed form goes through its
 * logarithm for the same reason, and the fit inherits that accuracy.
 */

#include "nfft3.h"
#include "infft.h"

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

/* phi(nx) for one window, with the constants it needs in par. */
typedef R (*window_eval)(const R *par, const R nx);

/* The fit of every window: coef as Y(kb_poly_fit) documents it. */
static void poly_fit(R *coef, const window_eval phi, const R *par, const INT m,
    const INT deg)
{
  const INT w = WINDOW_POLY_COLS(m), n = deg + 1;
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
      val[k * w + l] = phi(par, node[k] + (R)(m + 1 - l));

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

/* par: b, lg_tail, peak_inv, reach. */
static R kb_eval(const R *par, const R nx)
{
  return Y(kb_phi)(par[0], par[1], par[2], par[3], nx);
}

void Y(kb_poly_fit)(R *coef, const R b, const R lg_tail, const R peak_inv,
    const R reach, const INT m, const INT deg)
{
  const R par[4] = {b, lg_tail, peak_inv, reach};

  poly_fit(coef, kb_eval, par, m, deg);
}

R *Y(kb_poly_init)(const R *b, const INT d, const INT m, const R reach)
{
  const INT deg = Y(window_poly_degree)(m);
  const INT stride = (deg + 1) * WINDOW_POLY_COLS(m);
  R *tab = (R*) Y(malloc)((size_t)(stride * d) * sizeof(R));
  INT t;

  for (t = 0; t < d; t++)
    Y(kb_poly_fit)(tab + t * stride, b[t], b[d + t], b[3 * d + t], reach, m,
        deg);

  return tab;
}

/* par: w, m. The argument is formed as the closed-form run forms it. */
static R sincpow_eval(const R *par, const R nx)
{
  return Y(sincpow_phi)(par[0], par[1], KPI * par[0] * nx);
}

void Y(sincpow_poly_fit)(R *coef, const R w, const INT m, const INT deg)
{
  const R par[2] = {w, (R)m};

  poly_fit(coef, sincpow_eval, par, m, deg);
}

R *Y(sincpow_poly_init)(const R *b, const INT d, const INT m)
{
  const INT deg = Y(window_poly_degree)(m);
  const INT stride = (deg + 1) * WINDOW_POLY_COLS(m);
  R *tab = (R*) Y(malloc)((size_t)(stride * d) * sizeof(R));
  INT t;

  for (t = 0; t < d; t++)
    Y(sincpow_poly_fit)(tab + t * stride, b[t], m, deg);

  return tab;
}
