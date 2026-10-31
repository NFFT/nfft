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
 * infft.h evaluate the tables, a run of 2m + 2 taps by Horner in s^2 over
 * m + 1 mirrored columns. Y(bspline_poly_init) and Y(sincpow_poly_init)
 * build the B-spline and sinc-power tables the same way. */

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

  /* Thresholds are rounded down. A sigma off [1.1, 4] gets the degree at the
   * nearer end. */
  for (k = kb_poly_tab_off[mt]; k < kb_poly_tab_off[mt + 1]; k++)
    deg += IF(sigma >= K(1E-4) * (R)kb_poly_tab_brk[k], 1, 0);

  return IF(deg > KB_POLY_DEG_MAX, KB_POLY_DEG_MAX, deg);
}

typedef R (*sampler)(const void *ctx, R x);

/* Samples phi at deg + 1 Chebyshev nodes per cell into val and fits cheb,
 * both n by w; fit reuses val, node and cm for its refinement step. */
static void cheb_fit(R *cheb, R *node, R *cm, R *val, const INT m,
    const INT deg, sampler phi, const void *ctx)
{
  const INT w = KB_POLY_COLS(m), n = deg + 1;
  INT j, k, l;

  for (k = 0; k < n; k++)
    node[k] = COS(KPI * ((R)k + K(0.5)) / (R)n) / K(2.0);

  for (j = 0; j < n; j++)
    for (k = 0; k < n; k++)
      cm[j * n + k] = COS(KPI * (R)j * ((R)k + K(0.5)) / (R)n);

  for (k = 0; k < n; k++)
    for (l = 0; l < w; l++)
      val[k * w + l] = phi(ctx, ((R)(m + 1 - l) + K(0.5)) + node[k]);

  cheb_of(cheb, val, cm, w, deg);
}

/* The lowest even degree d with 2 sum_{j>d} max_cells |cheb_j| <= thr, at
 * most KB_POLY_DEG_MAX. Even, because the leading error is then odd in the
 * cell and has no mean. */
static INT cut_degree(const R *cheb, const INT w, const INT dg, const R thr)
{
  R sum = K(0.0);
  INT deg, l;

  for (deg = dg; deg >= 0; deg--)
  {
    R tail = K(0.0);

    for (l = 0; l < w; l++)
      tail = MAX(tail, FABS(cheb[deg * w + l]));

    if (K(2.0) * (sum + tail) > thr)
      break;

    sum += tail;
  }

  deg = IF(deg < 0, 0, deg);
  deg += deg & 1;

  return IF(deg > KB_POLY_DEG_MAX, KB_POLY_DEG_MAX, deg);
}

static void fit(R *coef, const INT m, const INT deg, sampler phi,
    const void *ctx)
{
  const INT w = KB_POLY_COLS(m), n = KB_POLY_ROWS(deg), dg = n - 1;
  /* val, cheb, res are n by w; cm and tp are n by n; node is n. */
  R *mem = (R*) Y(malloc)((size_t)(3 * n * w + 2 * n * n + n) * sizeof(R));
  R *val = mem, *cheb = val + n * w, *res = cheb + n * w;
  R *cm = res + n * w, *tp = cm + n * n, *node = tp + n * n;
  INT j, k, l;

  cheb_fit(cheb, node, cm, val, m, dg, phi, ctx);
  cheb_to_monomial(coef, cheb, tp, w, dg);

  /* One refinement step absorbs the rounding of the fit and of the monomial
   * conversion. The residual is a few ulp, so its own fit error is
   * negligible; a second step changes nothing. */
  for (k = 0; k < n; k++)
    for (l = 0; l < w; l++)
    {
      R v, mirror;

      KB_POLY_HORNER(1, &v, &mirror, coef + l, w, KB_POLY_HALF(dg), node[k]);
      res[k * w + l] = val[k * w + l] - v;
    }

  cheb_of(cheb, res, cm, w, dg);
  cheb_to_monomial(res, cheb, tp, w, dg);

  for (j = 0; j < n; j++)
    for (l = 0; l < w; l++)
      coef[j * w + l] += res[j * w + l];

  Y(free)(mem);
}

struct kb_ctx { R b, lg_tail, peak_inv, reach; };

static R kb_sample(const void *ctx, const R x)
{
  const struct kb_ctx *c = (const struct kb_ctx*) ctx;

  return Y(kb_phi)(c->b, c->lg_tail, c->peak_inv, c->reach, x);
}

void Y(kb_poly_fit)(R *coef, const R b, const R lg_tail, const R peak_inv,
    const R reach, const INT m, const INT deg)
{
  const struct kb_ctx c = {b, lg_tail, peak_inv, reach};

  fit(coef, m, deg, kb_sample, &c);
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

struct bs_ctx { INT k; };

static R bs_sample(const void *ctx, const R x)
{
  const INT k = ((const struct bs_ctx*) ctx)->k;

  return Y(bsplines)(k, x + (R)k / K(2.0));
}

/* The pieces have degree k-1, so k Chebyshev nodes fit them exactly. The
 * series is cut by the gen/kbpolygen rule against E = 4 (2 sigma - 1)^-k, the
 * window's own transform error at the tightest axis. */
R *Y(bspline_poly_init)(const R *sigma, const INT d, const INT m)
{
  const INT k = 2 * m + 2, w = KB_POLY_COLS(m), n = k, dg = k - 1;
  const struct bs_ctx c = {k};
  /* val is n by w; cheb and coef n+1 by w; cm n by n; tp n+1 by n+1. */
  R *mem = (R*) Y(malloc)((size_t)(3 * (n + 1) * w + n * n
      + (n + 1) * (n + 1) + n) * sizeof(R));
  R *val = mem, *cheb = val + n * w, *coef = cheb + (n + 1) * w;
  R *cm = coef + (n + 1) * w, *tp = cm + n * n, *node = tp + (n + 1) * (n + 1);
  const R peak = Y(bsplines)(k, (R)k / K(2.0));
  R thr = K(0.0), *tab;
  INT l, t, deg;

  cheb_fit(cheb, node, cm, val, m, dg, bs_sample, &c);

  for (l = 0; l < w; l++)
    cheb[n * w + l] = K(0.0);

  for (t = 0; t < d; t++)
  {
    const R e = K(4.0) * POW(K(2.0) * sigma[t] - K(1.0), -(R)k);
    const R th = MAX(e / K(40.0), EPSILON / K(2.0)) * peak;

    thr = IF(t == 0 || th < thr, th, thr);
  }

  deg = cut_degree(cheb, w, dg, thr);
  /* An odd KB_POLY_DEG_MAX fills its last row from the next coefficient. */
  cheb_to_monomial(coef, cheb, tp, w, KB_POLY_ROWS(deg) - 1);

  tab = (R*) Y(malloc)((size_t)(1 + KB_POLY_ROWS(deg) * w) * sizeof(R));
  tab[0] = (R)deg;
  memcpy(tab + 1, coef, (size_t)(KB_POLY_ROWS(deg) * w) * sizeof(R));

  Y(free)(mem);

  return tab;
}

struct sp_ctx { R w, m; };

static R sp_sample(const void *ctx, const R x)
{
  const struct sp_ctx *c = (const struct sp_ctx*) ctx;

  return Y(sincpow_phi)(c->w, c->m, KPI * c->w * x);
}

/* The window is entire, so the fit converges geometrically, but a fit in
 * working precision cannot certify a cut: its coefficients carry 2 to 4 eps
 * of sampling noise where the true series reaches eps at about 14 of
 * KB_POLY_DEG_MAX 16 in double, and alike in every format. So the degree is
 * KB_POLY_DEG_MAX; no usable transform-error bound exists to cut lower. */
R *Y(sincpow_poly_init)(const R *w, const INT d, const INT m)
{
  const INT deg = KB_POLY_DEG_MAX, stride = KB_POLY_ROWS(deg) * KB_POLY_COLS(m);
  R *tab = (R*) Y(malloc)((size_t)(1 + stride * d) * sizeof(R));
  INT t;

  tab[0] = (R)deg;

  for (t = 0; t < d; t++)
  {
    const INT s = first_alike(w, t);
    const struct sp_ctx c = {w[t], (R)m};

    if (s == t)
      fit(tab + 1 + t * stride, m, deg, sp_sample, &c);
    else
      memcpy(tab + 1 + t * stride, tab + 1 + s * stride,
          (size_t)stride * sizeof(R));
  }

  return tab;
}
