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
 * m + 1 mirrored columns. */

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
  INT mt = (INT)CEIL(KB_POLY_BR_OPT / b) - 1, h, k;

  A(reach == (R)m + K(1.0));
  (void)reach;

  /* Past m_opt the degree does not rise. */
  mt = IF(m < mt, m, mt);
  mt = IF(mt < KB_POLY_TAB_MMAX, mt, KB_POLY_TAB_MMAX);
  mt = IF(mt > 1, mt, 1);
  h = kb_poly_tab_h0[mt];

  /* Thresholds are rounded down. A sigma off [1.1, 4] gets the degree at the
   * nearer end. */
  for (k = kb_poly_tab_off[mt]; k < kb_poly_tab_off[mt + 1]; k++)
    h += IF(sigma >= K(1E-4) * (R)kb_poly_tab_brk[k], 1, 0);

  /* The cap is a plateau degree; an odd one rounds up, as the plan did. */
  return IF(2 * h > KB_POLY_DEG_MAX, 2 * ((KB_POLY_DEG_MAX + 1) / 2), 2 * h);
}

#if E_MANT_DIG == 24
  #define KB_POLY_BR_LS KB_POLY_BR_LS_24
#elif E_MANT_DIG == 64
  #define KB_POLY_BR_LS KB_POLY_BR_LS_64
#elif E_MANT_DIG == 113
  #define KB_POLY_BR_LS KB_POLY_BR_LS_113
#else
  #define KB_POLY_BR_LS KB_POLY_BR_LS_53
#endif

/* The truncation error is still 30 times the round-off floor u I0(bR)/I0(aR)
 * of the type the fit runs in; past that the fit would see noise. */
int Y(kb_poly_correct_ok)(const R b, const INT m)
{
  return b * ((R)m + K(1.0)) <= KB_POLY_BR_LS;
}

void Y(kb_poly_fit)(R *coef, const R b, const R lg_tail, const R peak_inv,
    const R reach, const INT m, const INT deg)
{
  const INT w = KB_POLY_COLS(m), n = KB_POLY_ROWS(deg), dg = n - 1;
  /* val, cheb, res are n by w; cm and tp are n by n; node is n. */
  R *mem = (R*) Y(malloc)((size_t)(3 * n * w + 2 * n * n + n) * sizeof(R));
  R *val = mem, *cheb = val + n * w, *res = cheb + n * w;
  R *cm = res + n * w, *tp = cm + n * n, *node = tp + n * n;
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

/* Least-squares solution of the overdetermined system A x = y, A rows by
 * cols, row-major, by Householder QR in place; y is overwritten and x is
 * returned in its first cols entries. */
static void lstsq(E *A, E *y, const INT rows, const INT cols)
{
  INT i, j, k;

  for (k = 0; k < cols; k++)
  {
    E norm = EK(0.0), alpha, beta, s;

    for (i = k; i < rows; i++)
      norm += A[i * cols + k] * A[i * cols + k];

    norm = E_SQRT(norm);

    if (norm == EK(0.0))
      continue;

    alpha = IF(A[k * cols + k] > EK(0.0), -norm, norm);
    /* v = a_k - alpha e_k, stored in column k; beta = 1 / (v . v / 2). */
    A[k * cols + k] -= alpha;
    beta = EK(0.0);

    for (i = k; i < rows; i++)
      beta += A[i * cols + k] * A[i * cols + k];

    beta = EK(2.0) / beta;

    for (j = k + 1; j < cols; j++)
    {
      s = EK(0.0);

      for (i = k; i < rows; i++)
        s += A[i * cols + k] * A[i * cols + j];

      s *= beta;

      for (i = k; i < rows; i++)
        A[i * cols + j] -= s * A[i * cols + k];
    }

    s = EK(0.0);

    for (i = k; i < rows; i++)
      s += A[i * cols + k] * y[i];

    s *= beta;

    for (i = k; i < rows; i++)
      y[i] -= s * A[i * cols + k];

    A[k * cols + k] = alpha;
  }

  for (k = cols - 1; k >= 0; k--)
  {
    E s = y[k];

    for (j = k + 1; j < cols; j++)
      s -= A[k * cols + j] * y[j];

    y[k] = IF(A[k * cols + k] != EK(0.0), s / A[k * cols + k], EK(0.0));
  }
}

/* The transform sees the table only through the 2m + 2 taps at a node's
 * offset t and the deconvolution C(xi), xi = k / n; its relative error at
 * frequency xi is D(xi, t) / C(xi) - 1 with D the sum of the taps times
 * exp(-2 pi i xi (tap position)). D is linear in the coefficients, so a
 * least-squares fit of a correction to the residual D / C - 1 on a grid of
 * (xi, t) lowers that error directly, the function the table approximates no
 * longer mattering. The residual is the difference of quantities near one,
 * so the fit resolves it only while the error is well above eps_E: that is
 * what Y(kb_poly_correct_ok) guards, and why E rather than R where E is
 * wider. The unknowns are Chebyshev coefficients T_k(2s) per stencil column,
 * converted to monomials at the end as the fit itself is; column 0, outside
 * the stencil, keeps its interpolant. */
#define KB_POLY_CORR_NXI 40
#define KB_POLY_CORR_NT 32

void Y(kb_poly_correct)(R *coef, const R b, const R i0e_peak_inv,
    const INT m, const INT deg)
{
  const INT w = KB_POLY_COLS(m), h = KB_POLY_HALF(deg), nk = 2 * h + 1;
  const INT nc = m + 1, nu = nc * nk;
  const INT nxi = KB_POLY_CORR_NXI, nt = KB_POLY_CORR_NT;
  const INT rows = 2 * nxi * nt + nu, npt = nxi * nt;
  const E ximax = (EK(2.0) * E_PI - (E)b) / (EK(2.0) * E_PI); /* 1/(2 sigma) */
  E *A = (E*) Y(malloc)((size_t)(2 * rows * nu + 3 * rows + nt * nk + 2 * nc
      + npt + nu) * sizeof(E));
  E *As = A + rows * nu, *y = As + rows * nu, *ys = y + rows, *r = ys + rows;
  E *T = r + rows, *lo = T + nt * nk, *hi = lo + nc, *wt = hi + nc;
  E *best = wt + npt;
  R *mem = (R*) Y(malloc)((size_t)(2 * nk * w + nk * nk) * sizeof(R));
  R *cheb = mem, *res = cheb + nk * w, *tp = res + nk * w;
  E lambda = EK(0.0), best_max = EK(0.0), rms_ls = EK(0.0);
  INT i, q, c, k, j, row = 0, it;

  /* T_k(2s) at the offsets s = t - 1/2 of the grid. */
  for (q = 0; q < nt; q++)
  {
    const E x = EK(2.0) * (((E)q + EK(0.5)) / (E)nt - EK(0.5));

    T[q * nk] = EK(1.0);
    T[q * nk + 1] = x;

    for (k = 2; k < nk; k++)
      T[q * nk + k] = EK(2.0) * x * T[q * nk + k - 1] - T[q * nk + k - 2];
  }

  for (i = 0; i < nxi; i++)
  {
    /* Clustered towards the band ends, where the error peaks. */
    const E xi = ximax * (EK(1.0) - E_COS(E_PI * (E)i / (E)(nxi - 1))) / EK(2.0);
    const E C = (E)Y(kb_phi_hut)(b, i0e_peak_inv, (R)m + K(1.0), K(1.0), (R)xi);

    for (q = 0; q < nt; q++)
    {
      const E t = ((E)q + EK(0.5)) / (E)nt;
      const E ph = -EK(2.0) * E_PI * xi;
      E dre = EK(0.0), dim = EK(0.0);
      R lor, hir;

      /* The taps the transform uses at this offset: column c is the cell
       * j = m + 1 - c at s, and its mirror at -s. */
      for (c = 1; c <= nc; c++)
      {
        KB_POLY_HORNER(1, &lor, &hir, coef + c, w, h, (R)(t - EK(0.5)));
        lo[c - 1] = (E)lor;
        hi[c - 1] = (E)hir;
      }

      for (c = 1; c <= nc; c++)
      {
        const E jj = (E)(m + 1 - c);
        const E a1 = ph * (jj + t), a2 = ph * (t - jj - EK(1.0));
        const E c1 = E_COS(a1), s1 = E_SIN(a1), c2 = E_COS(a2), s2 = E_SIN(a2);

        dre += lo[c - 1] * c1 + hi[c - 1] * c2;
        dim += lo[c - 1] * s1 + hi[c - 1] * s2;

        for (k = 0; k < nk; k++)
        {
          const E tk = T[q * nk + k], sg = IF(k & 1, -tk, tk);

          A[row * nu + (c - 1) * nk + k] = (tk * c1 + sg * c2) / C;
          A[(row + 1) * nu + (c - 1) * nk + k] = (tk * s1 + sg * s2) / C;
        }
      }

      y[row] = -(dre / C - EK(1.0));
      y[row + 1] = -dim / C;
      row += 2;
    }
  }

  /* A small ridge, so that an over-provisioned degree cannot turn the
   * unconstrained directions into large coefficients. */
  for (i = 0; i < row * nu; i++)
    lambda += A[i] * A[i];

  lambda = E_SQRT(EK(1E-8) * lambda / (E)nu);

  for (j = 0; j < nu; j++)
  {
    for (k = 0; k < nu; k++)
      A[(row + j) * nu + k] = IF(j == k, lambda, EK(0.0));

    y[row + j] = EK(0.0);
  }

  /* Plain least squares first, the RMS optimum, then a few floored Lawson
   * reweightings towards the worst case: the residual of an iterate is
   * linear in the correction, r + A x, so no re-evaluation is needed. A
   * reweighted iterate replaces the best one only if it lowers the maximum
   * and raises the RMS by at most a tenth, so both measures end up at or
   * below the least-squares ones in the worst case and the RMS. */
  /* Row weights: the clustered frequency samples carry their quadrature
   * weight, so the least-squares iterate is the RMS optimum of the uniform
   * band, the measure random data see, while the samples still sit densely
   * where the error peaks. */
  for (i = 0; i < nxi; i++)
  {
    const E cc = E_SIN(E_PI * (E)i / (E)(nxi - 1));
    const E wq = IF(i == 0 || i == nxi - 1, EK(0.5), EK(1.0)) * (cc + EK(0.02));

    for (q = 0; q < nt; q++)
      wt[i * nt + q] = wq;
  }

  for (it = 0; it < 4; it++)
  {
    E mx = EK(0.0), ss = EK(0.0);

    for (i = 0; i < rows * nu; i++)
      As[i] = A[i];

    for (i = 0; i < rows; i++)
      ys[i] = y[i];

    for (i = 0; i < npt; i++)
    {
      const E sw = E_SQRT(wt[i]);

      for (k = 0; k < nu; k++)
      {
        As[(2 * i) * nu + k] *= sw;
        As[(2 * i + 1) * nu + k] *= sw;
      }

      ys[2 * i] *= sw;
      ys[2 * i + 1] *= sw;
    }

    lstsq(As, ys, rows, nu);

    for (i = 0; i < npt; i++)
    {
      E re = -y[2 * i], im = -y[2 * i + 1];

      for (k = 0; k < nu; k++)
      {
        re += A[(2 * i) * nu + k] * ys[k];
        im += A[(2 * i + 1) * nu + k] * ys[k];
      }

      r[i] = E_SQRT(re * re + im * im);
      mx = IF(r[i] > mx, r[i], mx);
      ss += re * re + im * im;
    }

    if (it == 0)
      rms_ls = ss;

    if (it == 0 || (mx < best_max && ss <= EK(1.21) * rms_ls))
    {
      best_max = mx;

      for (k = 0; k < nu; k++)
        best[k] = ys[k];
    }

    for (i = 0; i < npt; i++)
      wt[i] *= IF(r[i] > EK(0.5) * mx, r[i], EK(0.5) * mx) / mx;
  }

  /* The correction, Chebyshev per column, to monomials and onto the table. */
  for (k = 0; k < nk; k++)
    for (c = 0; c < w; c++)
      cheb[k * w + c] = IF(c == 0, K(0.0), (R)best[(c - 1) * nk + k]);

  cheb_to_monomial(res, cheb, tp, w, deg);

  for (k = 0; k < nk; k++)
    for (c = 1; c < w; c++)
      coef[k * w + c] += res[k * w + c];

  Y(free)(mem);
  Y(free)(A);
}

R *Y(kb_poly_init)(const R *b, const INT d, const INT m, const R reach,
    const unsigned opt)
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
    {
      Y(kb_poly_fit)(tab + 1 + t * stride, b[t], b[d + t], b[3 * d + t],
          reach, m, deg);

      if (opt && Y(kb_poly_correct_ok)(b[t], m))
        Y(kb_poly_correct)(tab + 1 + t * stride, b[t], b[2 * d + t], m, deg);
    }
    else
      memcpy(tab + 1 + t * stride, tab + 1 + s * stride,
          (size_t)stride * sizeof(R));
  }

  return tab;
}
