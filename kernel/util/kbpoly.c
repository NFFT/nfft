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

/* Polynomial form of the Kaiser-Bessel window: one Chebyshev interpolant per
 * unit cell of |nx|, converted to monomials in the offset from the cell
 * centre. Each cell is analytic: both branches of the window are even in the
 * radical and ra^2 is a polynomial in the offset, so the fit converges
 * geometrically. The fit samples Y(kb_phi); a plain sinh(b ra)/ra has 20 times
 * its error, which would cap the fit. Y(kb_poly_run) and Y(kb_poly_phi) in
 * infft.h evaluate the tables, a run of 2m + 2 taps as one Horner sweep in s^2
 * over m + 1 mirrored columns. */

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

/* Series in T_j(2s) to powers of s. tp is (deg + 1)^2 scratch. */
static void cheb_to_monomial(R *coef, const R *cheb, R *tp, const INT w,
    const INT deg)
{
  const INT n = deg + 1;
  INT i, j, l;

  for (i = 0; i < n * n; i++)
    tp[i] = K(0.0);

  tp[0] = K(1.0);

  if (deg >= 1)
    tp[n + 1] = K(2.0);

  for (j = 2; j <= deg; j++)
    for (i = 0; i <= j; i++)
    {
      R v = -tp[(j - 2) * n + i];

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

/* All columns at offset s, evaluated as Y(kb_poly_run) does. */
static void eval_split(R *out, const R *coef, const INT w, const INT deg,
    const R s)
{
  const INT h = KB_POLY_HALF(deg);
  const R u = s * s;
  INT j, l;

  for (l = 0; l < w; l++)
  {
    R e = coef[(2 * h) * w + l], o = coef[(2 * h - 1) * w + l];

    for (j = h - 1; j >= 1; j--)
    {
      e = e * u + coef[(2 * j) * w + l];
      o = o * u + coef[(2 * j - 1) * w + l];
    }

    e = e * u + coef[l];
    out[l] = e + s * o;
  }
}

INT Y(kb_poly_degree)(const R b, const R reach, const INT m)
{
  const R sigma = KPI / (K2PI - b);
  /* Clamped to the table before the cast, which overflows for a tiny or
   * nonpositive b (oversampling up to 1/2); the clamps below give the same. */
  R mt_r = CEIL(KB_POLY_BR_OPT / b) - K(1.0);
  INT mt, deg, k;

  A(reach == (R)m + K(1.0));
  (void)reach;

  mt_r = mt_r < (R)KB_POLY_TAB_MMAX ? mt_r : (R)KB_POLY_TAB_MMAX;
  mt_r = mt_r > K(1.0) ? mt_r : K(1.0);
  mt = (INT)mt_r;

  /* Past m_opt the degree does not rise. */
  mt = IF(m < mt, m, mt);
  mt = IF(mt < KB_POLY_TAB_MMAX, mt, KB_POLY_TAB_MMAX);
  mt = IF(mt > 1, mt, 1);
  deg = kb_poly_tab_deg0[mt];

  /* Thresholds are rounded down. A sigma off [1.1, 4] gets the degree at the
   * nearer end. */
  for (k = kb_poly_tab_off[mt]; k < kb_poly_tab_off[mt + 1]; k++)
    deg += IF(sigma >= K(1E-4) * (R)kb_poly_tab_brk[k], 1, 0);

  return IF(deg > KB_POLY_DEG_MAX, KB_POLY_DEG_MAX, deg);
}

void Y(kb_poly_fit)(R *coef, const R b, const R lg_tail, const R peak_inv,
    const R reach, const INT m, const INT deg)
{
  const INT w = KB_POLY_COLS(m), n = KB_POLY_ROWS(deg), dg = n - 1;
  /* val, cheb, res are n by w; cm and tp are n by n; node is n; out is w. */
  R *mem = (R*) Y(malloc)((size_t)(3 * n * w + 2 * n * n + n + w) * sizeof(R));
  R *val = mem, *cheb = val + n * w, *res = cheb + n * w;
  R *cm = res + n * w, *tp = cm + n * n, *node = tp + n * n, *out = node + n;
  INT j, k, l;

  for (k = 0; k < n; k++)
    node[k] = COS(KPI * ((R)k + K(0.5)) / (R)n) / K(2.0);

  for (j = 0; j < n; j++)
    for (k = 0; k < n; k++)
      cm[j * n + k] = COS(KPI * (R)j * ((R)k + K(0.5)) / (R)n);

  for (k = 0; k < n; k++)
    for (l = 0; l < w; l++)
      val[k * w + l] = Y(kb_phi)(b, lg_tail, peak_inv, reach,
          ((R)(m + 1 - l) + K(0.5)) + node[k]);

  cheb_of(cheb, val, cm, w, dg);
  cheb_to_monomial(coef, cheb, tp, w, dg);

  /* One refinement step against the run-time evaluator absorbs the rounding of
   * the fit and of the monomial conversion. The residual is a few ulp, so its
   * own fit error is negligible; a second step changes nothing. */
  for (k = 0; k < n; k++)
  {
    eval_split(out, coef, w, dg, node[k]);

    for (l = 0; l < w; l++)
      res[k * w + l] = val[k * w + l] - out[l];
  }

  cheb_of(cheb, res, cm, w, dg);
  cheb_to_monomial(res, cheb, tp, w, dg);

  for (j = 0; j < n; j++)
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

  /* One degree for all axes, the largest any axis needs. */
  for (t = 0; t < d; t++)
  {
    const INT dt = Y(kb_poly_degree)(b[t], reach, m);

    deg = IF(dt > deg, dt, deg);
  }

  stride = KB_POLY_ROWS(deg) * KB_POLY_COLS(m);
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

void Y(kb_poly_run_block)(R *dst, const INT ld, const R *tab, const INT m,
    const INT deg, const R n, const R *x, const INT *u, const INT nb)
{
  const INT p = m + 1, cs = KB_POLY_COLS(m), h = KB_POLY_HALF(deg);
  const R *coef = tab + 1;
  /* the block in the innermost loop, so the Horner sweeps of nb nodes run
   * side by side in the vector lanes; the runs land in tmp column by column
   * and are transposed into dst at the end */
  R s[NFFT_NODE_BLOCK], q[NFFT_NODE_BLOCK], e[NFFT_NODE_BLOCK],
    o[NFFT_NODE_BLOCK];
  R tmp[NFFT_NODE_BLOCK * 2 * (KB_POLY_RUN_BLOCK_MAX_M + 1)];
  INT v, l, j;

  if (m > KB_POLY_RUN_BLOCK_MAX_M)
  {
    for (v = 0; v < nb; v++)
      Y(kb_poly_run)(dst + v * ld, tab, m, deg, NX_SUB(n, x[v], u[v]));
    return;
  }

  for (v = 0; v < nb; v++)
  {
    s[v] = NX_SUB(n, x[v], u[v]) - ((R)m + K(0.5));
    q[v] = s[v] * s[v];
  }

  for (l = 0; l < p; l++)
  {
    const R ce = coef[(2 * h) * cs + l], co = coef[(2 * h - 1) * cs + l];

    for (v = 0; v < nb; v++)
    {
      e[v] = ce;
      o[v] = co;
    }

    for (j = h - 1; j >= 1; j--)
    {
      const R a = coef[(2 * j) * cs + l], b = coef[(2 * j - 1) * cs + l];

      for (v = 0; v < nb; v++)
      {
        e[v] = e[v] * q[v] + a;
        o[v] = o[v] * q[v] + b;
      }
    }

    for (v = 0; v < nb; v++)
    {
      const R ev = e[v] * q[v] + coef[l];

      tmp[l * NFFT_NODE_BLOCK + v] = ev + s[v] * o[v];
      tmp[(2 * p - 1 - l) * NFFT_NODE_BLOCK + v] = ev - s[v] * o[v];
    }
  }

  for (v = 0; v < nb; v++)
    for (l = 0; l < 2 * p; l++)
      dst[v * ld + l] = tmp[l * NFFT_NODE_BLOCK + v];
}
