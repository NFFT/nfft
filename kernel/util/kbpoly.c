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
 * m + 1 mirrored columns. Y(kb_poly_correct) optimizes a table and its
 * deconvolution against the transform's error for OPT_WINDOW and
 * OPT_WINDOW_RMS. */

#include "nfft3.h"
#include "infft.h"
#include "kbpoly_data.h"
#include "kbpoly_opt_data.h"

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

#if MANT_DIG == 24
  #define KB_POLY_BR_LS KB_POLY_BR_LS_24
#elif MANT_DIG == 64
  #define KB_POLY_BR_LS KB_POLY_BR_LS_64
#elif MANT_DIG == 113
  #define KB_POLY_BR_LS KB_POLY_BR_LS_113
#else
  #define KB_POLY_BR_LS KB_POLY_BR_LS_53
#endif

/* The transform error is still 100 times the round-off floor
 * u I0(bR) / I0(aR) of R. */
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

/* Reflections are stored LAPACK style: v below the diagonal with v_j = 1
 * implied, tau in tau. The design matrix of the correction has condition
 * number up to 1e10, which rules out the normal equations in any E. */
INT Y(kb_poly_qr)(E *A, const INT rows, const INT cols, const E trunc,
    INT *perm, E *tau)
{
  INT i, j, k, p, n = IF(rows < cols, rows, cols);
  E *nrm = (E *)Y(malloc)((size_t)(2 * cols) * sizeof(E)), *nrm0 = nrm + cols;

  for (j = 0; j < cols; j++)
  {
    E s = EK(0.0);

    for (i = 0; i < rows; i++)
      s += A[j * rows + i] * A[j * rows + i];

    nrm[j] = nrm0[j] = E_SQRT(s);
    perm[j] = j;
  }

  for (j = 0; j < n; j++)
  {
    E *v, s = EK(0.0), beta;

    p = j;
    for (k = j + 1; k < cols; k++)
      if (nrm[k] > nrm[p])
        p = k;

    if (p != j)
    {
      E t;
      INT ti = perm[j];

      perm[j] = perm[p];
      perm[p] = ti;
      t = nrm[j];
      nrm[j] = nrm[p];
      nrm[p] = t;
      t = nrm0[j];
      nrm0[j] = nrm0[p];
      nrm0[p] = t;
      for (i = 0; i < rows; i++)
      {
        t = A[j * rows + i];
        A[j * rows + i] = A[p * rows + i];
        A[p * rows + i] = t;
      }
    }

    v = A + j * rows;
    for (i = j + 1; i < rows; i++)
      s += v[i] * v[i];

    if (s == EK(0.0))
      tau[j] = EK(0.0);
    else
    {
      const E nv = E_SQRT(v[j] * v[j] + s);

      beta = IF(v[j] > EK(0.0), -nv, nv);
      tau[j] = (beta - v[j]) / beta;
      s = EK(1.0) / (v[j] - beta);
      for (i = j + 1; i < rows; i++)
        v[i] *= s;
      v[j] = beta;
    }

    for (k = j + 1; k < cols; k++)
    {
      E *u = A + k * rows, w = u[j];

      for (i = j + 1; i < rows; i++)
        w += v[i] * u[i];

      w *= tau[j];
      u[j] -= w;
      for (i = j + 1; i < rows; i++)
        u[i] -= w * v[i];

      /* Downdate the column norm; recompute once it lost 90 %. */
      s = nrm[k] * nrm[k] - u[j] * u[j];
      if (s > EK(0.0) && E_SQRT(s) >= EK(0.1) * nrm0[k])
        nrm[k] = E_SQRT(s);
      else
      {
        s = EK(0.0);
        for (i = j + 1; i < rows; i++)
          s += u[i] * u[i];
        nrm[k] = nrm0[k] = E_SQRT(s);
      }
    }
  }

  Y(free)(nrm);

  for (k = 0; k < n; k++)
    if (E_FABS(A[k * rows + k]) <= trunc * E_FABS(A[0]))
      break;

  return k;
}

void Y(kb_poly_qr_basis)(const E *A, const INT rows, const INT cols,
    const INT k, const INT *perm, E *T)
{
  INT i, l, c;

  for (i = 0; i < cols * k; i++)
    T[i] = EK(0.0);

  /* Column c of R_k^{-1}, written to the row perm[i] of T. */
  for (c = 0; c < k; c++)
  {
    E *x = T + c * cols;

    x[perm[c]] = EK(1.0) / A[c * rows + c];
    for (i = c - 1; i >= 0; i--)
    {
      E s = EK(0.0);

      for (l = i + 1; l <= c; l++)
        s += A[l * rows + i] * x[perm[l]];

      x[perm[i]] = -s / A[i * rows + i];
    }
  }
}

void Y(kb_poly_ls)(E *A, E *y, E *x, const INT rows, const INT cols)
{
  INT i, j, k;
  INT *perm = (INT *)Y(malloc)((size_t)cols * sizeof(INT));
  E *tau = (E *)Y(malloc)((size_t)cols * sizeof(E));

  k = Y(kb_poly_qr)(A, rows, cols, EK(0.0), perm, tau);

  /* y := Q^T y */
  for (j = 0; j < k; j++)
  {
    E w = y[j];

    for (i = j + 1; i < rows; i++)
      w += A[j * rows + i] * y[i];

    w *= tau[j];
    y[j] -= w;
    for (i = j + 1; i < rows; i++)
      y[i] -= w * A[j * rows + i];
  }

  for (j = 0; j < cols; j++)
    x[j] = EK(0.0);

  for (j = k - 1; j >= 0; j--)
  {
    E s = y[j];

    for (i = j + 1; i < k; i++)
      s -= A[i * rows + j] * x[perm[i]];

    x[perm[j]] = s / A[j * rows + j];
  }

  Y(free)(perm);
  Y(free)(tau);
}

/* Gauss-Legendre nodes and weights on [0, 1]. */
static void gauleg(E *x, E *w, const INT n)
{
  INT i, j, s;

  for (i = 0; i < (n + 1) / 2; i++)
  {
    E z = E_COS(E_PI * ((E)i + EK(0.75)) / ((E)n + EK(0.5))), pp = EK(1.0);

    for (s = 0; s < 8; s++)
    {
      E p1 = EK(1.0), p2 = EK(0.0);

      for (j = 0; j < n; j++)
      {
        const E p3 = p2;

        p2 = p1;
        p1 = ((EK(2.0) * (E)j + EK(1.0)) * z * p2 - (E)j * p3)
            / ((E)j + EK(1.0));
      }

      pp = (E)n * (z * p1 - p2) / (z * z - EK(1.0));
      z -= p1 / pp;
    }

    x[i] = (EK(1.0) - z) / EK(2.0);
    x[n - 1 - i] = (EK(1.0) + z) / EK(2.0);
    w[i] = w[n - 1 - i] = EK(1.0) / ((EK(1.0) - z * z) * pp * pp);
  }
}

/* T_k(2 (t - 1/2)), k < nk. */
static void cheb_row(E *T, const INT nk, const E t)
{
  const E x = EK(2.0) * (t - EK(0.5));
  INT k;

  T[0] = EK(1.0);
  T[1] = x;

  for (k = 2; k < nk; k++)
    T[k] = EK(2.0) * x * T[k - 1] - T[k - 2];
}

/* The 2m + 2 taps at offset t: lo[c] at cell j = m - c (column c + 1 at
 * s = t - 1/2), hi[c] at the mirror cell -j - 1 (the same column at -s). */
static void taps(E *lo, E *hi, const R *coef, const INT m, const INT h,
    const E t)
{
  const INT w = KB_POLY_COLS(m);
  INT c;

  for (c = 0; c <= m; c++)
  {
    R lor, hir;

    KB_POLY_HORNER(1, &lor, &hir, coef + c + 1, w, h, (R)(t - EK(0.5)));
    lo[c] = (E)lor;
    hi[c] = (E)hir;
  }
}

/* Phases e^(-2 pi i xi y) of a sample set: ej, four per frequency and cell,
 * at the anchors y = j and y = -j - 1; et at the offsets; eg at the
 * quadrature nodes. */
static void phases(E *ej, E *et, E *eg, const E *xi, const INT nxi,
    const E *t, const INT ntt, const E *tg, const INT ng, const INT m)
{
  const E ph = -EK(2.0) * E_PI;
  INT i, c, q;

  for (i = 0; i < nxi; i++)
  {
    for (c = 0; c <= m; c++)
    {
      const E a = ph * xi[i] * (E)(m - c), am = -ph * xi[i] * (E)(m - c + 1);
      E *e = ej + 4 * (i * (m + 1) + c);

      e[0] = E_COS(a);
      e[1] = E_SIN(a);
      e[2] = E_COS(am);
      e[3] = E_SIN(am);
    }

    for (q = 0; q < ntt; q++)
    {
      const E a = ph * xi[i] * t[q];

      et[2 * (i * ntt + q)] = E_COS(a);
      et[2 * (i * ntt + q) + 1] = E_SIN(a);
    }

    for (q = 0; q < ng; q++)
    {
      const E a = ph * xi[i] * tg[q];

      eg[2 * (i * ng + q)] = E_COS(a);
      eg[2 * (i * ng + q) + 1] = E_SIN(a);
    }
  }
}

/* psi_hat(xi) of the table by quadrature in t; real, the window is even. */
static E psi_hat(const R *coef, const INT m, const INT h, const E *ej,
    const E *eg, const E *tg, const E *wg, const INT ng, E *lg, E *hg)
{
  E sum = EK(0.0);
  INT q, c;

  for (q = 0; q < ng; q++)
  {
    const E cg = eg[2 * q], sg = eg[2 * q + 1];

    taps(lg, hg, coef, m, h, tg[q]);

    for (c = 0; c <= m; c++)
      sum += wg[q] * (lg[c] * (ej[4 * c] * cg - ej[4 * c + 1] * sg)
          + hg[c] * (ej[4 * c + 2] * cg - ej[4 * c + 3] * sg));
  }

  return sum;
}

/* Error |D - C| / C of the table over a sample set with the phases of
 * phases(), through the run-time Horner; C = psi_hat by quadrature or, given
 * hut, from that table as the transform evaluates it. Returns the maximum,
 * negative if C is not positive or the error not finite somewhere; *rms the
 * root mean square with the frequency weights xw; den gets C, rr the error
 * per sample. */
static E grid_max(const R *coef, const R *hut, const INT m, const INT h,
    const E *xi, const E *xw, const INT nxi, const INT ntt, const E *t,
    const E *tg, const E *wg, const INT ng, const E *ej, const E *et,
    const E *eg, E *den, E *rr, E *lo, E *hi, E *lg, E *hg, E *rms)
{
  const INT nc = m + 1;
  E mx = EK(0.0), ss = EK(0.0), sw = EK(0.0);
  INT i, q, c;

  for (q = 0; q < ntt; q++)
    taps(lo + q * nc, hi + q * nc, coef, m, h, t[q]);

  for (i = 0; i < nxi; i++)
  {
    den[i] = IF(hut == 0,
        psi_hat(coef, m, h, ej + 4 * i * nc, eg + 2 * i * ng, tg, wg, ng, lg,
            hg),
        (E)Y(kb_poly_phi_hut)(hut, K(0.0), K(0.0), K(0.0), K(1.0), (R)xi[i]));

    if (!(den[i] > EK(0.0)))
      return EK(-1.0);

    for (q = 0; q < ntt; q++)
    {
      const INT row = i * ntt + q;
      const E ct = et[2 * row], st = et[2 * row + 1];
      E re = -den[i], im = EK(0.0), e;

      for (c = 0; c < nc; c++)
      {
        const E *ph = ej + 4 * (i * nc + c);

        re += lo[q * nc + c] * (ph[0] * ct - ph[1] * st)
            + hi[q * nc + c] * (ph[2] * ct - ph[3] * st);
        im += lo[q * nc + c] * (ph[0] * st + ph[1] * ct)
            + hi[q * nc + c] * (ph[2] * st + ph[3] * ct);
      }

      rr[2 * row] = re / den[i];
      rr[2 * row + 1] = im / den[i];
      e = E_SQRT(re * re + im * im) / den[i];

      if (!(e <= EK(1E30)))
        return EK(-1.0);

      mx = IF(e > mx, e, mx);
      ss += xw[i] * e * e;
      sw += xw[i];
    }
  }

  *rms = E_SQRT(ss / sw);

  return mx;
}

/* psi_hat(xi) per unknown into Bf (nu), with the phases of phases() at one
 * frequency. The mirror tap carries the mirror anchor phase and T_k(-2s) =
 * (-1)^k T_k(2s). Tg is nk scratch. */
static void psi_hat_basis(E *Bf, const E *ej, const E *eg, const INT nc,
    const INT nk, const E *tg, const E *wg, const INT ng, E *Tg)
{
  INT q, c, k;

  for (k = 0; k < nc * nk; k++)
    Bf[k] = EK(0.0);

  for (q = 0; q < ng; q++)
  {
    const E cg = eg[2 * q], sg = eg[2 * q + 1];

    cheb_row(Tg, nk, tg[q]);

    for (c = 0; c < nc; c++)
    {
      const E *e = ej + 4 * c;
      const E a1 = e[0] * cg - e[1] * sg, a2 = e[2] * cg - e[3] * sg;

      for (k = 0; k < nk; k++)
        Bf[c * nk + k] += wg[q] * Tg[k] * (a1 + IF(k & 1, -a2, a2));
    }
  }
}

/* Design set: npt points (xi[p], t[p]) with capacity cap. Per point the
 * complex row G (D - C per unknown, nr unknowns after the elimination of
 * unknown piv) and the real row B (psi_hat per unknown, reduced to nr the
 * same way). Rows act on the correction d of the Chebyshev unknowns; the
 * current window coef enters as gc (G c) and bc (B c). */
typedef struct
{
  INT npt, cap, nr, piv;
  E *xi, *t, *w; /* points and weights of the first LS pass */
  E *G; /* 2 * cap * nr, interleaved re/im, row-major */
  E *B; /* cap * nr */
  E *gc, *bc; /* 2 * cap, cap: G c and B c of the current window */
  E *el; /* nr: unknown piv = sum_j el[j] d_j, j over the others */
  const R *coef; /* the current window */
} kb_design;

/* Appends one point: phases, psi_hat per unknown (Gauss-Legendre nodes tg,
 * weights wg, ng of them), row G and row B, gc and bc through the
 * transform's taps and psi_hat of coef. */
static void design_add(kb_design *d, const E xi, const E t, const E w,
    const INT m, const INT h, const E *tg, const E *wg, const INT ng)
{
  const INT nc = m + 1, nk = 2 * h + 1, nu = nc * nk, nr = d->nr, p = d->npt;
  E *ws = (E*) Y(malloc)((size_t)(8 * nc + 2 + 2 * ng + 2 * nk + 3 * nu)
      * sizeof(E));
  E *ej = ws, *et = ej + 4 * nc, *eg = et + 2, *T = eg + 2 * ng;
  E *Tg = T + nk, *Bf = Tg + nk, *Gf = Bf + nu, *lo = Gf + 2 * nu;
  E *hi = lo + nc, *lg = hi + nc, *hg = lg + nc;
  E *gp = d->G + 2 * p * nr, *bp = d->B + p * nr;
  E C, re, im;
  INT c, k, j;

  A(p < d->cap);
  d->xi[p] = xi;
  d->t[p] = t;
  d->w[p] = w;

  phases(ej, et, eg, &xi, 1, &t, 1, tg, ng, m);
  psi_hat_basis(Bf, ej, eg, nc, nk, tg, wg, ng, Tg);
  cheb_row(T, nk, t);

  for (c = 0; c < nc; c++)
  {
    const E *e = ej + 4 * c;
    const E c1 = e[0] * et[0] - e[1] * et[1], s1 = e[0] * et[1] + e[1] * et[0];
    const E c2 = e[2] * et[0] - e[3] * et[1], s2 = e[2] * et[1] + e[3] * et[0];

    for (k = 0; k < nk; k++)
    {
      const E tk = T[k], sg = IF(k & 1, -tk, tk);

      Gf[2 * (c * nk + k)] = tk * c1 + sg * c2 - Bf[c * nk + k];
      Gf[2 * (c * nk + k) + 1] = tk * s1 + sg * s2;
    }
  }

  for (j = 0; j < nr; j++)
  {
    const INT u = j + IF(j >= d->piv, 1, 0);

    gp[2 * j] = Gf[2 * u] + Gf[2 * d->piv] * d->el[j];
    gp[2 * j + 1] = Gf[2 * u + 1] + Gf[2 * d->piv + 1] * d->el[j];
    bp[j] = Bf[u] + Bf[d->piv] * d->el[j];
  }

  taps(lo, hi, d->coef, m, h, t);
  C = psi_hat(d->coef, m, h, ej, eg, tg, wg, ng, lg, hg);
  re = -C;
  im = EK(0.0);

  for (c = 0; c < nc; c++)
  {
    const E *e = ej + 4 * c;

    re += lo[c] * (e[0] * et[0] - e[1] * et[1])
        + hi[c] * (e[2] * et[0] - e[3] * et[1]);
    im += lo[c] * (e[0] * et[1] + e[1] * et[0])
        + hi[c] * (e[2] * et[1] + e[3] * et[0]);
  }

  d->gc[2 * p] = re;
  d->gc[2 * p + 1] = im;
  d->bc[p] = C;
  d->npt++;

  Y(free)(ws);
}

/* max_i |gc_i + G_i dtot| / (bc_i + B_i dtot); returns -1 if a denominator
 * is <= 0 or the maximum not finite; fills gcn (2 npt), bcn (npt). */
static E design_ratio(const kb_design *d, const E *dtot, E *gcn, E *bcn)
{
  const INT nr = d->nr;
  E mx = EK(0.0);
  INT p, j;

  for (p = 0; p < d->npt; p++)
  {
    const E *g = d->G + 2 * p * nr, *bb = d->B + p * nr;
    E re = d->gc[2 * p], im = d->gc[2 * p + 1], den = d->bc[p], r;

    for (j = 0; j < nr; j++)
    {
      re += g[2 * j] * dtot[j];
      im += g[2 * j + 1] * dtot[j];
      den += bb[j] * dtot[j];
    }

    gcn[2 * p] = re;
    gcn[2 * p + 1] = im;
    bcn[p] = den;

    if (!(den > EK(0.0)))
      return EK(-1.0);

    r = E_SQRT(re * re + im * im) / den;
    mx = IF(r > mx, r, mx);
  }

  return IF(Y(e_finite)(mx), mx, EK(-1.0));
}

/* The transform sees the table through the 2m + 2 taps at a node's offset
 * t and the deconvolution C(xi); its relative error at frequency xi is
 * E(xi, t) = (D(xi, t) - C(xi)) / C(xi), D the tap sum with the phases
 * e^(-2 pi i xi (tap position)). With C the window's own transform psi_hat,
 * D and C are linear in the coefficients c. On a design of points the
 * problem is min over c of max_i |g_i c| / (b_i c), g = D - C, b = psi_hat;
 * the ratio is scale free, so sum_i b_i c over the initial design is fixed.
 * The level set |g_i c| <= tau b_i c is an intersection of second-order
 * cones: the problem is quasiconvex, its optimum the smallest tau with a
 * nonempty level set (Keiner, Potts: Optimized piecewise polynomial window
 * functions for the NFFT, 2026).
 *
 * Start: the interpolant plus one least-squares pass with the band's
 * quadrature weights, the RMS optimum. RMS mode repeats the pass with the
 * new denominators and keeps the best one. A table cell that is not an exact
 * hit replaces the start in the stencil columns and adds its exchange
 * points; the optimum of the cell above at the same degree, over 1.05, is
 * the first lower bound. An exact hit copies the cell and skips the design.
 *
 * Driver: bracket [lo, hi], hi the true ratio of the best window, lo the
 * largest level proven infeasible or left undecided. The model at the
 * current window scales the rows by 1 / (b_i c) and keeps the columns of a
 * pivoted QR with |R_kk| > TRUNC |R_00|; Y(socp3_solve) minimizes the slack
 * z of one cone per point at level tau. A Dinkelbach step solves at tau = hi and
 * accepts the window if its true ratio is below tau. After a step that
 * gains less than GAIN, bisection at sqrt(lo hi) (hi / 4 while lo = 0)
 * decides only the sign of z. Stop at hi <= (1 + DELTA) lo or after SOLVES
 * solves.
 *
 * Exchange: after the driver the ADD largest local maxima of |E| on a fine
 * grid join the design, lo carried over, until the fine maximum is within
 * DELTA_X of hi, after ROUNDS rounds, or after two rounds below FLOOR eps_E,
 * where the fine grid sees rounding noise. The round with the smallest fine
 * maximum wins.
 *
 * Unknowns are the coefficients of T_k(2s), k <= 2h, of the m + 1 stencil
 * columns, so an odd deg gets row 2h too; column 0, outside the stencil,
 * keeps its interpolant. Offsets are the limit t -> 0, where a window with
 * jumps peaks, and the midpoints of [0, 1/2]: E(xi, 1 - t) is the conjugate
 * of E(xi, t). The residual is a difference of quantities near one, so the
 * fit resolves it only well above eps_E; Y(kb_poly_correct_ok) guards that
 * coarsely, the hold-out exactly. The design capacity holds the initial
 * design, the points of a table cell and the exchange rounds. */
#define KB_POLY_CORR_TRUNC (EK(500.0) * E_EPSILON)
#define KB_POLY_CORR_NXI 24
#define KB_POLY_CORR_NT 16
#define KB_POLY_CORR_CELL 400
#define KB_POLY_CORR_ADD 40
#define KB_POLY_CORR_ROUNDS 10
#define KB_POLY_CORR_NIT 24
#define KB_POLY_CORR_STALL 4
#define KB_POLY_CORR_DELTA EK(0.01)
#define KB_POLY_CORR_GAIN EK(1.3)
#define KB_POLY_CORR_SOLVES 30
#define KB_POLY_CORR_MAXIT 60
#define KB_POLY_CORR_TOL (KB_POLY_CORR_DELTA / EK(100.0))
#define KB_POLY_CORR_FXI 401
#define KB_POLY_CORR_FT 64
#define KB_POLY_CORR_DELTA_X EK(0.02)
#define KB_POLY_CORR_FLOOR EK(500.0)
#define KB_POLY_CORR_ACCEPT EK(0.9)

/* E in the workspace of worst_points. */
static INT worst_work(const INT m, const INT ng)
{
  const INT nf = KB_POLY_CORR_FXI, nt = KB_POLY_CORR_FT + 1, nc = m + 1;

  return nf * (3 + 4 * nc + 4 * nt + 2 * ng) + nt * (1 + 2 * nc) + 2 * nc;
}

/* |E| of coef, C = psi_hat, on FXI frequencies uniform on [0, ximax] and the
 * offsets 0 (the limit) and (q + 0.5) / FT. The ADD largest local maxima
 * (>= 4 neighbours, t mirrored at the ends), each refined on an 11 x 11
 * patch of +-1 grid step (xi clipped, t mod 1), go to xi, t by decreasing
 * |E|, *nloc of them. Returns the fine-grid maximum before refinement, -1
 * if a denominator is <= 0. ws: worst_work(m, ng). */
static E worst_points(const R *coef, const INT m, const INT h, const E ximax,
    const E *tg, const E *wg, const INT ng, E *ws, E *xi, E *t, INT *nloc)
{
  const INT nf = KB_POLY_CORR_FXI, nt = KB_POLY_CORR_FT + 1, nc = m + 1;
  const E dx = ximax / (E)(nf - 1), dt = EK(1.0) / (E)KB_POLY_CORR_FT;
  E *fx = ws, *one = fx + nf, *den = one + nf, *ej = den + nf;
  E *et = ej + 4 * nf * nc, *rr = et + 2 * nf * nt, *eg = rr + 2 * nf * nt;
  E *ft = eg + 2 * nf * ng, *lo = ft + nt, *hi = lo + nt * nc;
  E *lg = hi + nt * nc, *hg = lg + nc;
  E ev[KB_POLY_CORR_ADD], rms, mx;
  INT i, q, a, k, n = 0;

  for (i = 0; i < nf; i++)
  {
    fx[i] = dx * (E)i;
    one[i] = EK(1.0);
  }

  for (q = 0; q < nt; q++)
    ft[q] = IF(q == 0, EK(0.0), ((E)q - EK(0.5)) * dt);

  phases(ej, et, eg, fx, nf, ft, nt, tg, ng, m);
  mx = grid_max(coef, 0, m, h, fx, one, nf, nt, ft, tg, wg, ng, ej, et, eg,
      den, rr, lo, hi, lg, hg, &rms);

  if (mx < EK(0.0))
    return EK(-1.0);

  for (i = 0; i < nf * nt; i++)
    rr[i] = E_SQRT(rr[2 * i] * rr[2 * i] + rr[2 * i + 1] * rr[2 * i + 1]);

  for (i = 0; i < nf; i++)
    for (q = 0; q < nt; q++)
    {
      const E v = rr[i * nt + q];
      const INT qd = IF(q == 0, 1, q - 1), qu = IF(q == nt - 1, q - 1, q + 1);

      if ((i > 0 && rr[(i - 1) * nt + q] > v)
          || (i < nf - 1 && rr[(i + 1) * nt + q] > v)
          || rr[i * nt + qd] > v || rr[i * nt + qu] > v
          || (n == KB_POLY_CORR_ADD && v <= ev[n - 1]))
        continue;

      for (a = IF(n < KB_POLY_CORR_ADD, n, n - 1); a > 0 && ev[a - 1] < v;
          a--)
      {
        ev[a] = ev[a - 1];
        xi[a] = xi[a - 1];
        t[a] = t[a - 1];
      }

      ev[a] = v;
      xi[a] = fx[i];
      t[a] = ft[q];
      n += IF(n < KB_POLY_CORR_ADD, 1, 0);
    }

  for (a = 0; a < n; a++)
  {
    for (k = 0; k < 11; k++)
    {
      const E x = xi[a] + (E)(k - 5) * dx / EK(5.0);
      const E s = t[a] + (E)(k - 5) * dt / EK(5.0);

      fx[k] = IF(x < EK(0.0), EK(0.0), IF(x > ximax, ximax, x));
      ft[k] = IF(s < EK(0.0), s + EK(1.0), IF(s >= EK(1.0), s - EK(1.0), s));
    }

    phases(ej, et, eg, fx, 11, ft, 11, tg, ng, m);

    if (grid_max(coef, 0, m, h, fx, one, 11, 11, ft, tg, wg, ng, ej, et, eg,
        den, rr, lo, hi, lg, hg, &rms) > ev[a])
      for (k = 0; k < 11 * 11; k++)
      {
        const E v = E_SQRT(rr[2 * k] * rr[2 * k]
            + rr[2 * k + 1] * rr[2 * k + 1]);

        if (v > ev[a])
        {
          ev[a] = v;
          xi[a] = fx[k / 11];
          t[a] = ft[k % 11];
        }
      }
  }

  for (a = 1; a < n; a++)
  {
    const E v = ev[a], x = xi[a], s = t[a];

    for (k = a; k > 0 && ev[k - 1] < v; k--)
    {
      ev[k] = ev[k - 1];
      xi[k] = xi[k - 1];
      t[k] = t[k - 1];
    }

    ev[k] = v;
    xi[k] = x;
    t[k] = s;
  }

  *nloc = n;

  return mx;
}

/* Preconditioner at the window with G c = gcn, B c = bcn: rows over bcn_i,
 * stacked re/im (As: 2 npt x nr, column-major), pivoted QR, T = P R_k^{-1}
 * (nr x k). Returns k; kappa = max |gcn_i| / bcn_i, E0n = gcn / (bcn kappa),
 * Gam = G T / bcn (re/im interleaved, npt x k), beta = B T kappa / bcn. */
static INT model(const kb_design *d, const E *gcn, const E *bcn, E *As,
    E *T, INT *perm, E *tau_h, E *kappa, E *E0n, E *Gam, E *beta)
{
  const INT npt = d->npt, nr = d->nr;
  E kp = EK(0.0);
  INT p, j, l, k;

  for (p = 0; p < npt; p++)
  {
    const E r = E_SQRT(gcn[2 * p] * gcn[2 * p]
        + gcn[2 * p + 1] * gcn[2 * p + 1]) / bcn[p];

    kp = IF(r > kp, r, kp);

    for (j = 0; j < nr; j++)
    {
      As[j * 2 * npt + 2 * p] = d->G[2 * (p * nr + j)] / bcn[p];
      As[j * 2 * npt + 2 * p + 1] = d->G[2 * (p * nr + j) + 1] / bcn[p];
    }
  }

  k = Y(kb_poly_qr)(As, 2 * npt, nr, KB_POLY_CORR_TRUNC, perm, tau_h);
  Y(kb_poly_qr_basis)(As, 2 * npt, nr, k, perm, T);

  for (p = 0; p < npt; p++)
  {
    const E *g = d->G + 2 * p * nr, *bb = d->B + p * nr;

    E0n[2 * p] = gcn[2 * p] / (bcn[p] * kp);
    E0n[2 * p + 1] = gcn[2 * p + 1] / (bcn[p] * kp);

    for (l = 0; l < k; l++)
    {
      const E *tl = T + l * nr;
      E re = EK(0.0), im = EK(0.0), be = EK(0.0);

      for (j = 0; j < nr; j++)
      {
        re += g[2 * j] * tl[j];
        im += g[2 * j + 1] * tl[j];
        be += bb[j] * tl[j];
      }

      Gam[2 * (p * k + l)] = re / bcn[p];
      Gam[2 * (p * k + l) + 1] = im / bcn[p];
      beta[p * k + l] = be * kp / bcn[p];
    }
  }

  *kappa = kp;

  return k;
}

/* Subproblem at level tau_n = tau / kappa, x = (w, z): min z subject to
 * |E0n_i + Gam_i w| <= tau_n (1 + beta_i w) + z. Cones h_i = (tau_n, Re
 * E0n_i, Im E0n_i), G_i = [[-tau_n beta_i, -1], [-Re Gam_i, 0], [-Im Gam_i,
 * 0]] (3 x (k + 1)), q = e_k. Layout as Y(socp3_solve) expects. */
static void subproblem(const INT npt, const INT k, const E tau_n,
    const E *E0n, const E *Gam, const E *beta, E *h, E *G, E *q)
{
  const INT n = k + 1;
  INT i, l;

  for (i = 0; i < npt; i++)
  {
    E *g0 = G + 3 * i * n, *g1 = g0 + n, *g2 = g1 + n;

    h[3 * i] = tau_n;
    h[3 * i + 1] = E0n[2 * i];
    h[3 * i + 2] = E0n[2 * i + 1];

    for (l = 0; l < k; l++)
    {
      g0[l] = -tau_n * beta[i * k + l];
      g1[l] = -Gam[2 * (i * k + l)];
      g2[l] = -Gam[2 * (i * k + l) + 1];
    }

    g0[k] = EK(-1.0);
    g1[k] = g2[k] = EK(0.0);
  }

  for (l = 0; l < k; l++)
    q[l] = EK(0.0);

  q[k] = EK(1.0);
}

/* dnew = dtot + kappa T w (w: k). */
static void apply_step(const INT nr, const INT k, const E *T, const E kappa,
    const E *w, const E *dtot, E *dnew)
{
  INT j, l;

  for (j = 0; j < nr; j++)
  {
    E s = EK(0.0);

    for (l = 0; l < k; l++)
      s += T[l * nr + j] * w[l];

    dnew[j] = dtot[j] + kappa * s;
  }
}

/* E in the driver's workspace for a design of capacity cap, nr unknowns. */
static INT driver_work(const INT cap, const INT nr)
{
  return cap * (11 + 8 * nr) + nr * (nr + 4) + 2 + Y(socp3_work)(cap, nr + 1);
}

/* Safeguarded Dinkelbach on the design d from dtot (in/out): hi the true
 * ratio of dtot, lo = *lo_io. Dinkelbach steps at tau = hi while they gain
 * KB_POLY_CORR_GAIN, else bisection at sqrt(lo hi), hi / 4 while lo = 0. A
 * window is accepted only if its true ratio is below tau; a level proven
 * infeasible, or undecided below hi, raises lo. y (3 npt) warm starts each
 * solve. ws: driver_work(cap, nr). Returns hi, -1 if dtot has no valid
 * ratio. */
static E driver(kb_design *d, E *dtot, E *lo_io, E *y, E *ws,
    const INT maxit)
{
  const INT cap = d->cap, nr = d->nr, npt = d->npt;
  E *As = ws, *Gam = As + 2 * cap * nr, *beta = Gam + 2 * cap * nr;
  E *G = beta + cap * nr, *h = G + 3 * cap * (nr + 1), *E0n = h + 3 * cap;
  E *gcn = E0n + 2 * cap, *bcn = gcn + 2 * cap, *T = bcn + cap;
  E *tau_h = T + nr * nr, *dnew = tau_h + nr, *q = dnew + nr;
  E *x = q + nr + 1, *work = x + nr + 1;
  INT *perm = (INT*) Y(malloc)((size_t)nr * sizeof(INT));
  E hi = design_ratio(d, dtot, gcn, bcn), lo = *lo_io, kappa = EK(0.0);
  E tau, r, pobj, dobj, pres, dres;
  INT k = 0, ns, j, it;
  int dink = 1, fresh = 0, st;

  for (ns = 0; hi >= EK(0.0) && ns < KB_POLY_CORR_SOLVES
      && hi > (EK(1.0) + KB_POLY_CORR_DELTA) * lo; ns++)
  {
    tau = IF(dink, hi, IF(lo > EK(0.0), E_SQRT(lo * hi), hi / EK(4.0)));

    if (!fresh)
    {
      k = model(d, gcn, bcn, As, T, perm, tau_h, &kappa, E0n, Gam, beta);
      fresh = 1;
    }

    subproblem(npt, k, tau / kappa, E0n, Gam, beta, h, G, q);

    for (j = 0; j <= k; j++)
      x[j] = EK(0.0);

    st = Y(socp3_solve)(npt, k + 1, h, G, q, x, y, 1, !dink,
        KB_POLY_CORR_TOL, E_EPSILON, maxit, work, &it, &pobj, &dobj, &pres,
        &dres);
    apply_step(nr, k, T, kappa, x, dtot, dnew);
    r = design_ratio(d, dnew, gcn, bcn);

    if (r >= EK(0.0) && r < tau)
    {
      dink = r * KB_POLY_CORR_GAIN < hi;
      hi = r;
      memcpy(dtot, dnew, (size_t)nr * sizeof(E));
      fresh = 0;
    }
    else
    {
      if (st == SOCP3_INFEASIBLE || tau < hi || (dres < KB_POLY_CORR_TOL
          && dobj > EK(10.0) * KB_POLY_CORR_TOL))
        lo = IF(tau > lo, tau, lo);

      dink = 0;
    }
  }

  Y(free)(perm);
  *lo_io = lo;

  return hi;
}

/* table = orig + the correction dv of the stencil columns; cheb, res
 * (nk w) and tp (nk nk) are scratch. */
static void corrected(R *table, const R *orig, const E *dv,
    const kb_design *d, const INT nk, const INT w, R *cheb, R *res, R *tp)
{
  INT k, c, j;

  for (k = 0; k < nk; k++)
    for (c = 0; c < w; c++)
      cheb[k * w + c] = K(0.0);

  for (j = 0; j < d->nr; j++)
  {
    const INT u = j + IF(j >= d->piv, 1, 0);

    cheb[(u % nk) * w + 1 + u / nk] = (R)dv[j];
    cheb[(d->piv % nk) * w + 1 + d->piv / nk] += (R)(d->el[j] * dv[j]);
  }

  cheb_to_monomial(res, cheb, tp, w, nk - 1);

  for (k = 0; k < nk; k++)
    for (c = 0; c < w; c++)
      table[k * w + c] = orig[k * w + c] + IF(c > 0, res[k * w + c], K(0.0));
}

/* sigma_cell == sigma, to the rounding of b in R. */
static int same_sigma(const R a, const R c)
{
  return FABS(a - c) <= (K(1E-9) + K(32.0) * EPSILON) * c;
}

INT Y(kb_poly_table)(const R b, const INT m, const INT deg,
    const kb_poly_opt_cell **exact, const kb_poly_opt_cell **below,
    const kb_poly_opt_cell **above)
{
  const R sigma = KPI / (K2PI - b);
  const INT h = KB_POLY_HALF(deg);
  const kb_poly_opt_cell *lo = NULL, *hi = NULL;
  INT r, i;

  for (r = 0; r < (INT)KB_POLY_OPT_NROWS; r++)
  {
    if (kb_poly_opt_rows[r].m != m)
      continue;

    for (i = 0; i < kb_poly_opt_rows[r].ncells; i++)
    {
      const kb_poly_opt_cell *c = kb_poly_opt_rows[r].cells + i;

      if (c->h > h)
        continue;

      if (c->sigma <= sigma || same_sigma(sigma, c->sigma))
        lo = IF(lo == NULL || c->sigma > lo->sigma
            || (c->sigma == lo->sigma && c->h > lo->h), c, lo);
      else
        hi = IF(hi == NULL || c->sigma < hi->sigma
            || (c->sigma == hi->sigma && c->h > hi->h), c, hi);
    }
  }

  *below = lo;
  *above = hi;
  *exact = IF(lo != NULL && lo->h == h && same_sigma(sigma, lo->sigma), lo,
      NULL);

  return IF(*exact != NULL, 1, 0) + IF(*below != NULL, 1, 0)
      + IF(*above != NULL, 1, 0);
}

INT Y(kb_poly_table_rows)(const kb_poly_opt_row **rows)
{
  *rows = kb_poly_opt_rows;

  return (INT)KB_POLY_OPT_NROWS;
}

/* Frequencies and offsets per column and half degree as gen/kbpolygen/optima.py
 * samples them; Gauss-Legendre with h + 24 nodes resolves psi_hat to the
 * 113-bit window errors. */
R Y(kb_poly_correct)(R *coef, R *hut, const R b, const INT m,
    const INT deg, const int minimax, const kb_poly_opt_cell *cell,
    const kb_poly_opt_cell *above)
{
  const INT w = KB_POLY_COLS(m), h = KB_POLY_HALF(deg), nk = 2 * h + 1;
  const INT nc = m + 1, nu = nc * nk, nr = nu - 1, n = nk * w;
  const INT nxi = IF(m <= 9, KB_POLY_CORR_NXI, 3 * (m + 1));
  const INT noff = IF(2 * h > KB_POLY_CORR_NT, 2 * h, KB_POLY_CORR_NT);
  const INT ntt = noff / 2 + 1;
  const INT ng = h + 24, npt = nxi * ntt, nh = KB_POLY_HUT_MAX;
  const INT cap = npt + KB_POLY_CORR_CELL
      + KB_POLY_CORR_ADD * KB_POLY_CORR_ROUNDS;
  const kb_poly_opt_cell *start = IF(minimax, cell, NULL);
  const int exact = start != NULL && start->h == h
      && same_sigma(KPI / (K2PI - b), start->sigma);
  const E ximax = (EK(2.0) * E_PI - (E)b) / (EK(2.0) * E_PI); /* 1/(2 sigma) */
  E *A = (E*) Y(malloc)((size_t)(cap * (6 + 3 * nr) + 2 * npt * nr + 5 * nr
      + 7 * npt + 2 * nu + 3 * nxi + 2 * ng + 4 * nxi * nc + 2 * nxi * ntt
      + 2 * nxi * ng + 2 * ntt * nc + 2 * ng * nc + 2 * nh + ntt + nk)
      * sizeof(E));
  E *Aw = A + cap * (6 + 3 * nr), *g = Aw + 2 * npt * nr, *scale = g + nr;
  E *dtot = scale + nr, *dbest = dtot + nr, *el = dbest + nr;
  E *gcn = el + nr, *bcn = gcn + 2 * npt, *yw = bcn + npt;
  E *rr = yw + 2 * npt, *Bsum = rr + 2 * nxi * ntt, *Bf = Bsum + nu;
  E *den = Bf + nu;
  E *xi = den + nxi, *xw = xi + nxi, *tg = xw + nxi, *wg = tg + ng;
  E *ej = wg + ng, *et = ej + 4 * nxi * nc, *eg = et + 2 * nxi * ntt;
  E *lo = eg + 2 * nxi * ng, *hi = lo + ntt * nc, *lg = hi + ntt * nc;
  E *hg = lg + ng * nc, *fh = hg + ng * nc, *fc = fh + nh, *ts = fc + nh;
  E *Tg = ts + ntt;
  R *mem = (R*) Y(malloc)((size_t)(5 * n + nk * nk) * sizeof(R));
  R *cheb = mem, *res = cheb + n, *tp = res + n, *best = tp + nk * nk;
  R *orig = best + n, *base = orig + n;
  kb_design d;
  E best_obj = EK(0.0), obj, plain, mx, rms;
  INT i, q, k, j, p, it, nt, stall = 0;

  d.npt = 0;
  d.cap = cap;
  d.nr = nr;
  d.xi = A;
  d.t = d.xi + cap;
  d.w = d.t + cap;
  d.G = d.w + cap;
  d.B = d.G + 2 * cap * nr;
  d.gc = d.B + cap * nr;
  d.bc = d.gc + 2 * cap;
  d.el = el;

  memcpy(orig, coef, (size_t)n * sizeof(R));

  /* The cell's rows 0 .. 2 h' over zero rows up to 2 h; a cell start keeps
   * the interpolant in column 0. */
  if (start != NULL)
    for (i = 0; i < n; i++)
      base[i] = IF(i % w == 0 && !exact, orig[i],
          IF(i < (2 * start->h + 1) * w, start->coef[i], K(0.0)));
  else
    memcpy(base, coef, (size_t)n * sizeof(R));

  d.coef = base;

  gauleg(tg, wg, ng);

  for (q = 0; q < ntt; q++)
    ts[q] = IF(q == 0, EK(0.0), ((E)q - EK(0.5)) / (E)noff);

  /* Frequencies clustered towards the band ends, where the error peaks,
   * with their quadrature weights. */
  for (i = 0; i < nxi; i++)
  {
    const E th = E_PI * (E)i / (E)(nxi - 1);

    xi[i] = ximax * (EK(1.0) - E_COS(th)) / EK(2.0);
    xw[i] = IF(i == 0 || i == nxi - 1, EK(0.5), EK(1.0))
        * (E_SIN(th) + EK(0.02));
  }

  if (exact)
    memcpy(best, base, (size_t)n * sizeof(R));
  else
  {
    /* The scale constraint over the initial design eliminates the unknown with
     * the largest |Bsum|. The band sum of T_0 of the outermost column changes
     * sign with sigma and m. */
    phases(ej, et, eg, xi, nxi, ts, 0, tg, ng, m);

    for (j = 0; j < nu; j++)
      Bsum[j] = EK(0.0);

    for (i = 0; i < nxi; i++)
    {
      psi_hat_basis(Bf, ej + 4 * i * nc, eg + 2 * i * ng, nc, nk, tg, wg, ng,
          Tg);

      for (j = 0; j < nu; j++)
        Bsum[j] += (E)ntt * Bf[j];
    }

    d.piv = 0;

    for (j = 1; j < nu; j++)
      d.piv = IF(E_FABS(Bsum[j]) > E_FABS(Bsum[d.piv]), j, d.piv);

    for (j = 0; j < nr; j++)
      el[j] = -Bsum[j + IF(j >= d.piv, 1, 0)] / Bsum[d.piv];

    for (i = 0; i < nxi; i++)
      for (q = 0; q < ntt; q++)
        design_add(&d, xi[i], ts[q], xw[i], m, h, tg, wg, ng);

    /* Columns scaled to unit norm for the LS. */
    for (j = 0; j < nr; j++)
    {
      E ss = EK(0.0);

      for (p = 0; p < npt; p++)
        ss += d.G[2 * (p * nr + j)] * d.G[2 * (p * nr + j)]
            + d.G[2 * (p * nr + j) + 1] * d.G[2 * (p * nr + j) + 1];

      scale[j] = IF(ss > EK(0.0), EK(1.0) / E_SQRT(ss), EK(1.0));
      dtot[j] = EK(0.0);
      dbest[j] = EK(0.0);
    }

    /* Weighted LS passes with the current denominators. The first is the RMS
     * optimum and the start of the minimax driver; RMS keeps the best pass. */
    for (it = 0; it <= KB_POLY_CORR_NIT; it++)
    {
      mx = design_ratio(&d, dtot, gcn, bcn);

      if (mx < EK(0.0))
        break;

      if (it > 0)
      {
        E ss = EK(0.0), sw = EK(0.0);

        if (minimax)
          break;

        for (p = 0; p < npt; p++)
        {
          ss += d.w[p] * (gcn[2 * p] * gcn[2 * p]
              + gcn[2 * p + 1] * gcn[2 * p + 1]) / (bcn[p] * bcn[p]);
          sw += d.w[p];
        }

        obj = E_SQRT(ss / sw);

        if (it == 1 || obj < best_obj)
        {
          best_obj = obj;
          stall = 0;
          memcpy(dbest, dtot, (size_t)nr * sizeof(E));
        }
        else
          stall++;

        if (it == KB_POLY_CORR_NIT || stall == KB_POLY_CORR_STALL)
          break;
      }

      /* Rows scaled by sqrt(w) / den, so the residual is the relative error. */
      for (p = 0; p < npt; p++)
      {
        const E wa = E_SQRT(d.w[p]) / bcn[p];

        for (j = 0; j < nr; j++)
        {
          Aw[j * 2 * npt + 2 * p] = wa * scale[j] * d.G[2 * (p * nr + j)];
          Aw[j * 2 * npt + 2 * p + 1] = wa * scale[j]
              * d.G[2 * (p * nr + j) + 1];
        }

        yw[2 * p] = -wa * gcn[2 * p];
        yw[2 * p + 1] = -wa * gcn[2 * p + 1];
      }

      Y(kb_poly_ls)(Aw, yw, g, 2 * npt, nr);

      for (j = 0; j < nr; j++)
        dtot[j] += g[j] * scale[j];
    }

    /* Exchange rounds: the fine-grid maxima join the design until the
     * fine-grid maximum is within DELTA_X of the design's; below the floor the
     * fine grid sees rounding noise, so two rounds at most. */
    if (minimax)
    {
      const INT nw = driver_work(cap, nr), nf = worst_work(m, ng);
      E *ws = (E*) Y(malloc)((size_t)(nw + 3 * cap + nf
          + 2 * KB_POLY_CORR_ADD) * sizeof(E));
      E *y = ws + nw, *fw = y + 3 * cap, *wx = fw + nf;
      E *wt = wx + KB_POLY_CORR_ADD, lo = EK(0.0), fbest = EK(-1.0), fine;
      INT r, nloc = 0;

      /* A cell's exchange points keep their place in the band; the optimum
       * at the larger sigma above, at the same degree, is nearly a lower
       * bound. The driver never lowers lo again. */
      if (start != NULL)
      {
        for (p = 0; p < start->npts && p < KB_POLY_CORR_CELL; p++)
          design_add(&d, ximax * (E)start->pts[2 * p] / EK(65535.0),
              (E)start->pts[2 * p + 1] / EK(65535.0), EK(1.0), m, h, tg, wg,
              ng);

        lo = IF(above != NULL && above->h == h, (E)above->e / EK(1.05),
            EK(0.0));
      }

      for (p = 0; p < cap; p++)
      {
        y[3 * p] = EK(1.0);
        y[3 * p + 1] = y[3 * p + 2] = EK(0.0);
      }

      for (r = 0; r < KB_POLY_CORR_ROUNDS; r++)
      {
        mx = driver(&d, dtot, &lo, y, ws, KB_POLY_CORR_MAXIT);
        corrected(best, base, dtot, &d, nk, w, cheb, res, tp);
        fine = worst_points(best, m, h, ximax, tg, wg, ng, fw, wx, wt, &nloc);

        if (r == 0 || (fine >= EK(0.0) && (fbest < EK(0.0) || fine < fbest)))
        {
          fbest = fine;
          memcpy(dbest, dtot, (size_t)nr * sizeof(E));
        }

        if (mx < EK(0.0) || fine < EK(0.0)
            || fine <= (EK(1.0) + KB_POLY_CORR_DELTA_X) * mx
            || (r >= 1 && mx < KB_POLY_CORR_FLOOR * E_EPSILON))
          break;

        for (p = 0; p < nloc; p++)
          design_add(&d, wx[p], wt[p], EK(1.0), m, h, tg, wg, ng);
      }

      Y(free)(ws);
    }

    corrected(best, base, dbest, &d, nk, w, cheb, res, tp);
  }

  /* The deconvolution table of the best pass (layout at KB_POLY_HUT_MAX),
   * trailing terms below the rounding of R dropped. */
  for (j = 0; j < nh; j++)
  {
    const E u = E_COS(E_PI * ((E)j + EK(0.5)) / (E)nh);
    const E x = ximax * E_SQRT((u + EK(1.0)) / EK(2.0));

    phases(ej, et, eg, &x, 1, ts, 0, tg, ng, m);
    fh[j] = psi_hat(best, m, h, ej, eg, tg, wg, ng, lg, hg);
  }

  for (k = 0; k < nh; k++)
  {
    E sum = EK(0.0);

    for (j = 0; j < nh; j++)
      sum += fh[j] * E_COS(E_PI * (E)k * ((E)j + EK(0.5)) / (E)nh);

    fc[k] = (EK(2.0) / (E)nh) * sum * IF(k == 0, EK(0.5), EK(1.0));
  }

  for (nt = nh; nt > 1 && FABS((R)fc[nt - 1])
      <= Y(float_property)(NFFT_EPSILON) * FABS((R)fc[0]) / K(8.0); nt--)
    ;

  hut[0] = (R)nt;
  hut[1] = (R)(EK(1.0) / (ximax * ximax));

  for (k = 0; k < nt; k++)
    hut[2 + k] = (R)fc[k];

  /* Hold-out: frequencies between the design nodes plus the band edge,
   * where the deconvolution amplifies most, offsets a quarter cell off the
   * design ones, C from the table as the transform evaluates it. A fit of
   * rounding noise does not carry over, so the optimized table replaces
   * the plain one only at KB_POLY_CORR_ACCEPT of its objective or less. */
  for (i = 0; i < nxi; i++)
    xi[i] = ximax * (EK(1.0) - E_COS(E_PI * ((E)i + EK(1.0)) / (E)nxi))
        / EK(2.0);

  for (q = 0; q < ntt; q++)
    ts[q] = ((E)q + EK(0.25)) / (E)noff;

  phases(ej, et, eg, xi, nxi, ts, ntt, tg, ng, m);
  mx = grid_max(best, hut, m, h, xi, xw, nxi, ntt, ts, tg, wg, ng, ej, et,
      eg, den, rr, lo, hi, lg, hg, &rms);
  obj = IF(minimax, mx, rms);
  plain = grid_max(orig, 0, m, h, xi, xw, nxi, ntt, ts, tg, wg, ng, ej, et,
      eg, den, rr, lo, hi, lg, hg, &rms);
  plain = IF(minimax, plain, rms);

  if (mx >= EK(0.0) && obj <= KB_POLY_CORR_ACCEPT * plain)
    memcpy(coef, best, (size_t)n * sizeof(R));
  else
  {
    memcpy(coef, orig, (size_t)n * sizeof(R));
    hut[0] = K(0.0);
    mx = EK(-1.0);
  }

  Y(free)(mem);
  Y(free)(A);

  return (R)mx;
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

  stride = KB_POLY_STRIDE(deg, m);
  tab = (R*) Y(malloc)((size_t)(1 + stride * d) * sizeof(R));
  tab[0] = (R)deg;

  for (t = 0; t < d; t++)
  {
    const INT s = first_alike(b, t);

    if (s == t)
    {
      R *coef = tab + 1 + t * stride;
      R *hut = coef + KB_POLY_ROWS(deg) * KB_POLY_COLS(m);

      Y(kb_poly_fit)(coef, b[t], b[d + t], b[3 * d + t], reach, m, deg);
      hut[0] = K(0.0);

      if ((opt & OPT_WINDOW) && Y(kb_poly_correct_ok)(b[t], m))
      {
        const kb_poly_opt_cell *exact, *below, *above;

        /* Off-grid from scratch: the cell start is slower at equal quality. */
        Y(kb_poly_table)(b[t], m, deg, &exact, &below, &above);
        Y(kb_poly_correct)(coef, hut, b[t], m, deg, 1, exact, NULL);
      }
      else if ((opt & OPT_WINDOW_RMS) && Y(kb_poly_correct_ok)(b[t], m))
        Y(kb_poly_correct)(coef, hut, b[t], m, deg, 0, NULL, NULL);
    }
    else
      memcpy(tab + 1 + t * stride, tab + 1 + s * stride,
          (size_t)stride * sizeof(R));
  }

  return tab;
}
