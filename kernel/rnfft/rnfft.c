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

/* Nonequispaced FFT of real-valued functions.
 *
 * The coefficients f_hat_k are given for the half spectrum H: 0 <= k_{d-1} <=
 * N_{d-1}/2 and -N_t/2 <= k_t < N_t/2 for t < d-1. With c_k = 1 for k_{d-1} in
 * {0, N_{d-1}/2} and c_k = 2 otherwise,
 *
 *   trafo:   f_j     = sum_{k in H} c_k Re(f_hat_k exp(-2 pi i k x_j)),
 *   adjoint: f_hat_k = c_k sum_j f_j exp(+2 pi i k x_j).
 *
 * The adjoint is the adjoint of the trafo for the inner products
 * sum_j f_j g_j and Re sum_k f_hat_k conj(g_hat_k). For a Hermitian spectrum,
 * f_hat_{-k} = conj(f_hat_k), the trafo is the complex NFFT.
 *
 * The fast trafo puts conj(c_k f_hat_k / phi_hat_k) / 2 on the half lattice of
 * the oversampled grid, Hermitian in the plane k_{d-1} = 0, so that FFTW's c2r
 * yields the real g_l = sum_{k in H} c_k Re(f_hat_k / phi_hat_k
 * exp(-2 pi i k l / n)), and sums f_j = sum_l g_l phi(x_j - l / n). The fast
 * adjoint is its exact transpose: real spreading, r2c, and f_hat_k =
 * c_k conj(g_hat_k) / phi_hat_k. */

#include "config.h"

#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3.h"
#include "infft.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#undef X
#define X(name) RNFFT(name)

/* Precompute flags whose window values depend on the nodes. */
#define PRE_NODE_PSI (PRE_PSI | PRE_FG_PSI | PRE_FULL_PSI)

/* Block size for the phase recurrence in the direct transforms. */
#define DIRECT_BLOCK 32

static inline INT intprod(const INT *vec, const INT d)
{
  INT t, p = 1;

  for (t = 0; t < d; t++)
    p *= vec[t];

  return p;
}

/* Number of rows of the oversampled grid, the dimensions before the last. */
static inline INT grid_rows(const X(plan) *ths)
{
  return ths->n_total / ths->n[ths->d - 1];
}

/* Row length of g_hat. */
static inline INT hat_len(const X(plan) *ths)
{
  return ths->n[ths->d - 1] / 2 + 1;
}

/* Row stride of g. In place, the rows hold the 2 (n_{d-1}/2 + 1) reals of the
 * row of g_hat they share memory with. */
static inline INT row_stride(const X(plan) *ths)
{
  return ths->flags & FFT_OUT_OF_PLACE ? ths->n[ths->d - 1] : 2 * hat_len(ths);
}

/* Grid index of node coordinate x shifted down by m. FLOOR alone, not
 * FLOOR(n x - m), so that the sort keys agree with the window runs. */
static inline INT grid_start(const R x, const INT n, const INT m)
{
  return (INT)FLOOR((R)n * x) - m;
}

/* frac(k x) in [-1/2, 1/2) */
static inline R reduced_frac(const R k, const R x)
{
  return FFMA(k, x, -RINT(k * x));
}

/* 2 pi sum_t k_t x_t, reduced mod 2 pi into [-pi, pi) */
static inline R reduced_omega(const INT d, const INT *k, const R *x)
{
  R s = K(0.0);
  INT t;

  for (t = 0; t < d; t++)
    s += reduced_frac((R)k[t], x[t]);

  return K2PI * (s - RINT(s));
}

/* Leading indices k_0, ..., k_{d-2} of row r of f_hat. */
static inline void row_index(const X(plan) *ths, INT r, INT *k)
{
  INT t;

  for (t = ths->d - 2; t >= 0; t--)
  {
    k[t] = r % ths->N[t] - ths->N[t] / 2;
    r /= ths->N[t];
  }
}

void X(trafo_direct)(const X(plan) *ths)
{
  const INT d = ths->d, Nh = ths->N[d - 1] / 2, H = Nh + 1,
    rows = ths->N_total / H;
  const C *f_hat = ths->f_hat;
  INT j;

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(j)
#endif
  for (j = 0; j < ths->M_total; j++)
  {
    const R *xj = ths->x + j * d;
    const R dphi = K2PI * xj[d - 1];
    const C dw = COS(dphi) - II * SIN(dphi);
    C v = K(0.0), e = K(0.0);
    INT k[d], r;

    row_index(ths, 0, k);
    for (r = 0; r < rows; r++)
    {
      const INT rb = r * H;
      INT i = 0, t;

      /* each block is seeded with an accurate phase at its centre */
      while (i < H)
      {
        const INT len = MIN(DIRECT_BLOCK, H - i);
        R omega;
        C w0;

        k[d - 1] = i + len / 2;
        omega = reduced_omega(d, k, xj);
        w0 = COS(omega) - II * SIN(omega);
#define TRAFO_STEP(ix, w) \
        { \
          const C s_ = f_hat[ix] * (w); \
          v += s_; \
          if ((ix) == rb || (ix) == rb + Nh) \
            e += s_; \
        }
        NFFT_RECUR_BLOCK(rb + i, len, w0, dw, TRAFO_STEP);
#undef TRAFO_STEP
        i += len;
      }

      for (t = d - 2; t >= 1 && k[t] == ths->N[t] / 2 - 1; t--)
        k[t] = -ths->N[t] / 2;
      if (t >= 0)
        k[t]++;
    }

    ths->f[j] = K(2.0) * CREAL(v) - CREAL(e);
  }
}

void X(adjoint_direct)(const X(plan) *ths)
{
  const INT d = ths->d, Nh = ths->N[d - 1] / 2, H = Nh + 1,
    rows = ths->N_total / H;
  C *f_hat = ths->f_hat;
  INT r;

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(r)
#endif
  for (r = 0; r < rows; r++)
  {
    C *fh = f_hat + r * H;
    INT k[d], j, l;

    row_index(ths, r, k);
    memset(fh, 0, (size_t)H * sizeof(C));

    for (j = 0; j < ths->M_total; j++)
    {
      const R *xj = ths->x + j * d, fj = ths->f[j];
      const R dphi = K2PI * xj[d - 1];
      const C dw = COS(dphi) + II * SIN(dphi);
      INT i = 0;

      while (i < H)
      {
        const INT len = MIN(DIRECT_BLOCK, H - i);
        R omega;
        C w0;

        k[d - 1] = i + len / 2;
        omega = reduced_omega(d, k, xj);
        w0 = COS(omega) + II * SIN(omega);
#define ADJOINT_STEP(ix, w) (fh[ix] += fj * (w))
        NFFT_RECUR_BLOCK(i, len, w0, dw, ADJOINT_STEP);
#undef ADJOINT_STEP
        i += len;
      }
    }

    for (l = 1; l < Nh; l++)
      fh[l] *= K(2.0);
  }
}

/* 1 / phi_hat(k) along leading dimension t, -N_t/2 <= k <= N_t/2 */
static inline R phi_inv(const X(plan) *ths, const INT t, const INT k)
{
  if ((ths->flags & PRE_PHI_HUT) && k < ths->N[t] / 2)
    return ths->c_phi_inv[t][k + ths->N[t] / 2];
  return K(1.0) / PHI_HUT(ths->n[t], k, t);
}

/* 1 / phi_hat(l) along the last dimension, 0 <= l <= N_{d-1}/2 */
static inline R phi_inv_last(const X(plan) *ths, const INT l)
{
  const INT t = ths->d - 1;

  if (ths->flags & PRE_PHI_HUT)
    return ths->c_phi_inv[t][l];
  return K(1.0) / PHI_HUT(ths->n[t], l, t);
}

/* Index of grid index q in the leading dimension t of f_hat, or -1. */
static inline INT coef_index(const X(plan) *ths, const INT t, const INT q)
{
  const INT n = ths->n[t], Nh = ths->N[t] / 2;

  if (q < Nh)
    return q + Nh;
  if (q >= n - Nh)
    return q - n + Nh;
  return -1;
}

/* g_hat from f_hat, row by row of the grid, each row written once. A grid
 * row q takes the coefficient row k' with q = k' mod n, if there is one, and
 * its plane l = 0 the Hermitian part (conj(v_q) + v_{-q}) / 2. */
static void D_trafo(X(plan) *ths)
{
  const INT d = ths->d, Nh = ths->N[d - 1] / 2, H = Nh + 1, h = hat_len(ths),
    rows = grid_rows(ths);
  INT q;

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(q)
#endif
  for (q = 0; q < rows; q++)
  {
    C *gq = ths->g_hat + q * h;
    INT r = 0, rm = 0, qq = q, rs = 1, t, l;
    R c = K(1.0);

    for (t = d - 2; t >= 0; t--)
    {
      const INT n = ths->n[t], qt = qq % n, i = coef_index(ths, t, qt),
        im = coef_index(ths, t, qt == 0 ? 0 : n - qt);

      qq /= n;
      r = (r < 0 || i < 0) ? -1 : r + i * rs;
      rm = (rm < 0 || im < 0) ? -1 : rm + im * rs;
      rs *= ths->N[t];
      if (i >= 0)
        c *= phi_inv(ths, t, i - ths->N[t] / 2);
      else if (im >= 0)
        c *= phi_inv(ths, t, im - ths->N[t] / 2);
    }

    if (r < 0 && rm < 0)
    {
      memset(gq, 0, (size_t)h * sizeof(C));
      continue;
    }

    memset(gq + H, 0, (size_t)(h - H) * sizeof(C));
    if (r >= 0)
    {
      const C *fh = ths->f_hat + r * H;

      for (l = 1; l < Nh; l++)
        gq[l] = CONJ(fh[l]) * (c * phi_inv_last(ths, l));
      gq[Nh] = CONJ(fh[Nh]) * (K(0.5) * c * phi_inv_last(ths, Nh));
    }
    else
      memset(gq + 1, 0, (size_t)Nh * sizeof(C));

    gq[0] = ((r >= 0 ? CONJ(ths->f_hat[r * H]) : K(0.0))
      + (rm >= 0 ? ths->f_hat[rm * H] : K(0.0)))
      * (K(0.5) * c * phi_inv_last(ths, 0));
  }
}

/* f_hat_k = c_k conj(g_hat_k) / phi_hat_k */
static void D_adjoint(X(plan) *ths)
{
  const INT d = ths->d, Nh = ths->N[d - 1] / 2, H = Nh + 1, h = hat_len(ths),
    rows = ths->N_total / H;
  INT r;

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(r)
#endif
  for (r = 0; r < rows; r++)
  {
    C *fh = ths->f_hat + r * H;
    const C *gq;
    INT k[d], q = 0, t, l;
    R c = K(1.0);

    row_index(ths, r, k);
    for (t = 0; t < d - 1; t++)
    {
      q = q * ths->n[t] + (k[t] < 0 ? k[t] + ths->n[t] : k[t]);
      c *= phi_inv(ths, t, k[t]);
    }
    gq = ths->g_hat + q * h;

    fh[0] = CONJ(gq[0]) * (c * phi_inv_last(ths, 0));
    for (l = 1; l < Nh; l++)
      fh[l] = CONJ(gq[l]) * (K(2.0) * c * phi_inv_last(ths, l));
    fh[Nh] = CONJ(gq[Nh]) * (c * phi_inv_last(ths, Nh));
  }
}

/* Window values psi[l] = phi(x_t - (u + l) / n_t), l < 2m+2, of node j along
 * dimension t, with u the grid start reduced to [0, n_t) in *u. The
 * precomputed tables hold node j in slot k, its place in the sorted order.
 * fg holds exp(-1/b) and its square per dimension, buf room for 2m+3
 * values. */
static inline const R *window(const X(plan) *ths, const INT j, const INT k,
  const INT t, const R *fg, R *buf, INT *u)
{
  const INT d = ths->d, n = ths->n[t], m = ths->m, w = 2 * m + 2;
  const R x = ths->x[j * d + t];
  const INT s = grid_start(x, n, m);
  INT l;

  *u = s < 0 ? s + n : s;

  if (ths->flags & PRE_PSI)
    return ths->psi + (k * d + t) * w;

  if (ths->flags & PRE_FG_PSI)
  {
    const R seed = ths->psi[2 * (k * d + t)];

    FG_RUN(buf, m, FABS(seed), ths->psi[2 * (k * d + t) + 1], fg[2 * t],
      fg[2 * t + 1]);
    return FG_RUN_START(buf, seed);
  }

  if (ths->flags & FG_PSI)
  {
    const R dx = NX_SUB(n, x, FLOOR((R)n * x + K(0.5)));

    FG_RUN(buf, m, PHI(n, dx / (R)n, t), EXP(K(2.0) * dx / ths->b[t]),
      fg[2 * t], fg[2 * t + 1]);
    return FG_RUN_START(buf, dx);
  }

  if (ths->flags & PRE_LIN_PSI)
  {
    const INT ip_s = ths->K / (m + 2);
    const R *tab = ths->psi + (ths->K + 1) * t;
    const R ip_y = FABS(NX_SUB(n, x, s)) * (R)ip_s;
    const INT ip_u = (INT)FLOOR(ip_y);
    const R ip_w = ip_y - (R)ip_u;

    for (l = 0; l < w; l++)
      buf[l] = tab[ABS(ip_u - l * ip_s)] * (K(1.0) - ip_w)
        + tab[ABS(ip_u - l * ip_s + 1)] * ip_w;
    return buf;
  }

  PHI_RUN(buf, n, x, s, t);
  return buf;
}

/* sum_l psi[l] row[(u + l) mod n] */
static inline R gather_dot(const R *restrict row, const R *restrict psi,
  const INT u, const INT n, const INT w)
{
  R s = K(0.0);
  INT l;

  if (u + w <= n)
  {
    for (l = 0; l < w; l++)
      s += psi[l] * row[u + l];
  }
  else
  {
    const INT k = n - u;

    for (l = 0; l < k; l++)
      s += psi[l] * row[u + l];
    for (; l < w; l++)
      s += psi[l] * row[l - k];
  }
  return s;
}

/* row[i] += a v[(i - u) mod n] for the indices i = (u + l) mod n, l < w, in
 * [lo, hi] */
static inline void scatter(R *restrict row, const R *restrict v, const R a,
  const INT u, const INT n, const INT w, const INT lo, const INT hi,
  const int atomic)
{
  const INT e0 = MIN(MIN(u + w, n) - 1, hi), e1 = MIN(u + w - n - 1, hi);
  INT i;

#ifdef _OPENMP
  if (atomic)
  {
    for (i = MAX(u, lo); i <= e0; i++)
    {
      #pragma omp atomic
      row[i] += a * v[i - u];
    }
    for (i = lo; i <= e1; i++)
    {
      #pragma omp atomic
      row[i] += a * v[i + n - u];
    }
    return;
  }
#else
  UNUSED(atomic);
#endif

  for (i = MAX(u, lo); i <= e0; i++)
    row[i] += a * v[i - u];
  for (i = lo; i <= e1; i++)
    row[i] += a * v[i + n - u];
}

/* f_j = sum_l g_l phi(x_j - l / n) */
static inline R node_trafo(const X(plan) *ths, const R *g, const INT gs,
  const INT j, const INT k, const R *fg)
{
  const INT d = ths->d, w = 2 * ths->m + 2, nl = ths->n[d - 1];
  R buf[d * (w + 1)], s = K(0.0);
  const R *psi[d];
  INT u[d], t, a, b;

  for (t = 0; t < d; t++)
    psi[t] = window(ths, j, k, t, fg, buf + t * (w + 1), &u[t]);

  if (d == 1)
    return gather_dot(g, psi[0], u[0], nl, w);

  if (d == 2)
  {
    const INT n0 = ths->n[0];

    for (a = 0; a < w; a++)
    {
      const INT i0 = u[0] + a < n0 ? u[0] + a : u[0] + a - n0;

      s += psi[0][a] * gather_dot(g + i0 * gs, psi[1], u[1], nl, w);
    }
  }
  else if (d == 3)
  {
    const INT n0 = ths->n[0], n1 = ths->n[1];

    for (a = 0; a < w; a++)
    {
      const INT i0 = u[0] + a < n0 ? u[0] + a : u[0] + a - n0;
      const R *g0 = g + i0 * n1 * gs;
      R s1 = K(0.0);

      for (b = 0; b < w; b++)
      {
        const INT i1 = u[1] + b < n1 ? u[1] + b : u[1] + b - n1;

        s1 += psi[1][b] * gather_dot(g0 + i1 * gs, psi[2], u[2], nl, w);
      }
      s += psi[0][a] * s1;
    }
  }
  else
  {
    const INT nb = ths->n[d - 2];
    INT c[d];

    for (t = 0; t < d - 2; t++)
      c[t] = 0;
    do
    {
      INT off = 0;
      R p = K(1.0), s1 = K(0.0);

      for (t = 0; t < d - 2; t++)
      {
        const INT i = u[t] + c[t] < ths->n[t] ? u[t] + c[t]
          : u[t] + c[t] - ths->n[t];

        off = off * ths->n[t] + i;
        p *= psi[t][c[t]];
      }
      off *= nb;

      for (b = 0; b < w; b++)
      {
        const INT i1 = u[d - 2] + b < nb ? u[d - 2] + b : u[d - 2] + b - nb;

        s1 += psi[d - 2][b] * gather_dot(g + (off + i1) * gs, psi[d - 1],
          u[d - 1], nl, w);
      }
      s += p * s1;

      for (t = d - 3; t >= 0 && ++c[t] == w; t--)
        c[t] = 0;
    } while (t >= 0);
  }

  return s;
}

/* g_l += f_j phi(x_j - l / n) for the grid points with first index in
 * [lo, hi] */
static inline void node_adjoint(const X(plan) *ths, R *g, const INT gs,
  const INT j, const INT k, const R *fg, const INT lo, const INT hi,
  const int atomic)
{
  const INT d = ths->d, w = 2 * ths->m + 2, nl = ths->n[d - 1];
  R buf[d * (w + 1)], v[w];
  const R *psi[d];
  INT u[d], t, l, a, b;

  for (t = 0; t < d; t++)
    psi[t] = window(ths, j, k, t, fg, buf + t * (w + 1), &u[t]);

  for (l = 0; l < w; l++)
    v[l] = ths->f[j] * psi[d - 1][l];

  if (d == 1)
    scatter(g, v, K(1.0), u[0], nl, w, lo, hi, atomic);
  else if (d == 2)
  {
    const INT n0 = ths->n[0];

    for (a = 0; a < w; a++)
    {
      const INT i0 = u[0] + a < n0 ? u[0] + a : u[0] + a - n0;

      if (i0 >= lo && i0 <= hi)
        scatter(g + i0 * gs, v, psi[0][a], u[1], nl, w, 0, nl - 1, atomic);
    }
  }
  else if (d == 3)
  {
    const INT n0 = ths->n[0], n1 = ths->n[1];

    for (a = 0; a < w; a++)
    {
      const INT i0 = u[0] + a < n0 ? u[0] + a : u[0] + a - n0;
      R *g0 = g + i0 * n1 * gs;

      if (i0 < lo || i0 > hi)
        continue;
      for (b = 0; b < w; b++)
      {
        const INT i1 = u[1] + b < n1 ? u[1] + b : u[1] + b - n1;

        scatter(g0 + i1 * gs, v, psi[0][a] * psi[1][b], u[2], nl, w, 0,
          nl - 1, atomic);
      }
    }
  }
  else
  {
    const INT nb = ths->n[d - 2];
    INT c[d];

    for (t = 0; t < d - 2; t++)
      c[t] = 0;
    do
    {
      const INT i0 = u[0] + c[0] < ths->n[0] ? u[0] + c[0]
        : u[0] + c[0] - ths->n[0];

      if (i0 >= lo && i0 <= hi)
      {
        INT off = 0;
        R p = K(1.0);

        for (t = 0; t < d - 2; t++)
        {
          const INT i = u[t] + c[t] < ths->n[t] ? u[t] + c[t]
            : u[t] + c[t] - ths->n[t];

          off = off * ths->n[t] + i;
          p *= psi[t][c[t]];
        }
        off *= nb;

        for (b = 0; b < w; b++)
        {
          const INT i1 = u[d - 2] + b < nb ? u[d - 2] + b
            : u[d - 2] + b - nb;

          scatter(g + (off + i1) * gs, v, p * psi[d - 2][b], u[d - 1], nl, w,
            0, nl - 1, atomic);
        }
      }

      for (t = d - 3; t >= 0 && ++c[t] == w; t--)
        c[t] = 0;
    } while (t >= 0);
  }
}

/* Sort the nodes by the linear index of their grid start. The keys go to
 * index_x[2j], the node indices to index_x[2j+1]; scratch is room for 2 M
 * more. */
static void sort_nodes(const X(plan) *ths, INT *scratch)
{
  const INT d = ths->d, M = ths->M_total;
  INT j, rhigh;

  for (j = 0; j < M; j++)
  {
    INT key = 0, t;

    for (t = 0; t < d; t++)
    {
      const INT s = grid_start(ths->x[j * d + t], ths->n[t], ths->m);

      key = key * ths->n[t] + (s < 0 ? s + ths->n[t] : s);
    }
    ths->index_x[2 * j] = key;
    ths->index_x[2 * j + 1] = j;
  }

  rhigh = (INT)LRINT(CEIL(LOG2((R)ths->n_total))) - 1;
  Y(sort_node_indices_radix_lsdf)(M, ths->index_x, scratch, MAX(rhigh, 0));
}

/* Sort at precompute time, with scratch of its own. */
static void sort_precompute(const X(plan) *ths)
{
  INT *scratch;

  if (!(ths->flags & NFFT_SORT_NODES))
    return;

  scratch = (INT*)Y(malloc)(2 * (size_t)ths->M_total * sizeof(INT));
  sort_nodes(ths, scratch);
  Y(free)(scratch);
}

/* Sort at transform time, if no precomputation did. The plan allocated
 * index_x with room for the scratch then. */
static void sort_transform(const X(plan) *ths)
{
  if ((ths->flags & NFFT_SORT_NODES) && !(ths->flags & PRE_NODE_PSI))
    sort_nodes(ths, ths->index_x + 2 * ths->M_total);
}

static inline INT node(const X(plan) *ths, const INT k)
{
  return (ths->flags & NFFT_SORT_NODES) ? ths->index_x[2 * k + 1] : k;
}

static void fg_init(const X(plan) *ths, R *fg)
{
  INT t;

  for (t = 0; t < ths->d; t++)
  {
    if (ths->flags & (FG_PSI | PRE_FG_PSI))
    {
      fg[2 * t] = EXP(K(-1.0) / ths->b[t]);
      fg[2 * t + 1] = fg[2 * t] * fg[2 * t];
    }
    else
      fg[2 * t] = fg[2 * t + 1] = K(0.0);
  }
}

static INT full_psi_len(const X(plan) *ths)
{
  INT t, p = 1;

  for (t = 0; t < ths->d; t++)
    p *= 2 * ths->m + 2;
  return p;
}

static void B_trafo(X(plan) *ths)
{
  const INT M = ths->M_total, gs = row_stride(ths);
  const R *g = ths->g;
  R fg[2 * ths->d];
  INT k;

  if (ths->flags & PRE_FULL_PSI)
  {
    const INT lprod = full_psi_len(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      const INT j = node(ths, k);
      const R *psi = ths->psi + k * lprod;
      const INT *ix = ths->psi_index_g + k * lprod;
      R s = K(0.0);
      INT l;

      for (l = 0; l < lprod; l++)
        s += psi[l] * g[ix[l]];
      ths->f[j] = s;
    }
    return;
  }

  sort_transform(ths);
  fg_init(ths, fg);

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    const INT j = node(ths, k);

    ths->f[j] = node_trafo(ths, g, gs, j, k, fg);
  }
}

#ifdef _OPENMP
static inline INT key_search(const INT *index_x, const INT M, const INT key)
{
  INT lo = 0, hi = M;

  while (lo < hi)
  {
    const INT mid = lo + (hi - lo) / 2;

    if (index_x[2 * mid] < key)
      lo = mid + 1;
    else
      hi = mid;
  }
  return lo;
}

/* The rows [*lo, *hi] of dimension 0 the calling thread writes, and the key
 * ranges [key[0], key[1]] and [key[2], key[3]] of the nodes whose window
 * reaches them. All -1 for a thread with no rows. */
static void block_init(const X(plan) *ths, INT *lo, INT *hi, INT *key)
{
  const INT n0 = ths->n[0], w = 2 * ths->m + 2, P = ths->n_total / n0,
    nt = MIN((INT)omp_get_num_threads(), n0), id = omp_get_thread_num();
  INT s;

  *lo = *hi = key[0] = key[1] = key[2] = key[3] = -1;
  if (id >= nt)
    return;

  *lo = id * (n0 / nt) + MIN(id, n0 % nt);
  *hi = *lo + n0 / nt + (id < n0 % nt ? 1 : 0) - 1;

  s = *lo - w + 1;
  key[1] = (*hi + 1) * P - 1;
  if (s >= 0)
    key[0] = s * P;
  else if (s + n0 <= *hi + 1)
  {
    key[0] = 0;
    key[1] = n0 * P - 1;
  }
  else
  {
    key[0] = 0;
    key[2] = (s + n0) * P;
    key[3] = n0 * P - 1;
  }
}
#endif

/* g_l += f_j psi_{j,l} over the precomputed products and grid indices */
static inline void full_psi_adjoint(X(plan) *ths, const INT j, const INT k,
  const INT lprod, const INT s0, const INT lo, const INT hi, const int atomic)
{
  const R *psi = ths->psi + k * lprod;
  const INT *ix = ths->psi_index_g + k * lprod;
  const R fj = ths->f[j];
  R *g = ths->g;
  INT l;

#ifdef _OPENMP
  if (atomic)
  {
    for (l = 0; l < lprod; l++)
    {
      #pragma omp atomic
      g[ix[l]] += psi[l] * fj;
    }
    return;
  }
#else
  UNUSED(atomic);
#endif

  for (l = 0; l < lprod; l++)
  {
    const INT i0 = ix[l] / s0;

    if (i0 >= lo && i0 <= hi)
      g[ix[l]] += psi[l] * fj;
  }
}

static void B_adjoint(X(plan) *ths)
{
  const INT M = ths->M_total, gs = row_stride(ths), n0 = ths->n[0],
    full = ths->flags & PRE_FULL_PSI,
    lprod = full ? full_psi_len(ths) : 0,
    s0 = ths->d == 1 ? 1 : ths->n_total / ths->n[0] / ths->n[ths->d - 1] * gs;
  R *g = ths->g;
  R fg[2 * ths->d];

  memset(g, 0, (size_t)(grid_rows(ths) * gs) * sizeof(R));

  if (!full)
  {
    sort_transform(ths);
    fg_init(ths, fg);
  }

#ifdef _OPENMP
  if (ths->flags & NFFT_OMP_BLOCKWISE_ADJOINT)
  {
    #pragma omp parallel default(shared)
    {
      INT lo, hi, key[4], part, k;

      block_init(ths, &lo, &hi, key);
      for (part = 0; part < 4; part += 2)
      {
        if (key[part] < 0)
          continue;
        for (k = key_search(ths->index_x, M, key[part]);
          k < M && ths->index_x[2 * k] <= key[part + 1]; k++)
        {
          const INT j = ths->index_x[2 * k + 1];

          if (full)
            full_psi_adjoint(ths, j, k, lprod, s0, lo, hi, 0);
          else
            node_adjoint(ths, g, gs, j, k, fg, lo, hi, 0);
        }
      }
    }
    return;
  }

  {
    INT k;

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < M; k++)
    {
      const INT j = node(ths, k);

      if (full)
        full_psi_adjoint(ths, j, k, lprod, s0, 0, n0 - 1, 1);
      else
        node_adjoint(ths, g, gs, j, k, fg, 0, n0 - 1, 1);
    }
  }
#else
  {
    INT k;

    for (k = 0; k < M; k++)
    {
      const INT j = node(ths, k);

      if (full)
        full_psi_adjoint(ths, j, k, lprod, s0, 0, n0 - 1, 0);
      else
        node_adjoint(ths, g, gs, j, k, fg, 0, n0 - 1, 0);
    }
  }
#endif
}

/* The fast transforms need a window that fits the grid. */
static int use_direct(const X(plan) *ths)
{
  INT t;

  for (t = 0; t < ths->d; t++)
    if (ths->N[t] <= ths->m || ths->n[t] <= 2 * ths->m + 2)
      return 1;
  return 0;
}

void X(trafo)(X(plan) *ths)
{
  if (use_direct(ths))
  {
    X(trafo_direct)(ths);
    return;
  }

  TIC(0)
  D_trafo(ths);
  TOC(0)

  TIC_FFTW(1)
  FFTW(execute)(ths->my_fftw_plan1);
  TOC_FFTW(1)

  TIC(2)
  B_trafo(ths);
  TOC(2)
}

void X(adjoint)(X(plan) *ths)
{
  if (use_direct(ths))
  {
    X(adjoint_direct)(ths);
    return;
  }

  TIC(2)
  B_adjoint(ths);
  TOC(2)

  TIC_FFTW(1)
  FFTW(execute)(ths->my_fftw_plan2);
  TOC_FFTW(1)

  TIC(0)
  D_adjoint(ths);
  TOC(0)
}

static void precompute_phi_hut(X(plan) *ths)
{
  const INT d = ths->d;
  INT t, k;

  ths->c_phi_inv = (R**)Y(malloc)((size_t)d * sizeof(R*));

  for (t = 0; t < d - 1; t++)
  {
    ths->c_phi_inv[t] = (R*)Y(malloc)((size_t)ths->N[t] * sizeof(R));
    for (k = 0; k < ths->N[t]; k++)
      ths->c_phi_inv[t][k] = K(1.0) / PHI_HUT(ths->n[t], k - ths->N[t] / 2, t);
  }

  ths->c_phi_inv[d - 1] = (R*)Y(malloc)((size_t)(ths->N[d - 1] / 2 + 1)
    * sizeof(R));
  for (k = 0; k <= ths->N[d - 1] / 2; k++)
    ths->c_phi_inv[d - 1][k] = K(1.0) / PHI_HUT(ths->n[d - 1], k, d - 1);
}

void X(precompute_lin_psi)(X(plan) *ths)
{
  INT t, j;

  for (t = 0; t < ths->d; t++)
  {
    const R step = (R)(ths->m + 2) / (R)(ths->K * ths->n[t]);

    for (j = 0; j <= ths->K; j++)
      ths->psi[(ths->K + 1) * t + j] = PHI(ths->n[t], (R)j * step, t);
  }
}

void X(precompute_fg_psi)(X(plan) *ths)
{
  const INT d = ths->d;
  INT t, k;

  sort_precompute(ths);

  for (t = 0; t < d; t++)
  {
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < ths->M_total; k++)
    {
      const R x = ths->x[node(ths, k) * d + t];
      const R dx = NX_SUB(ths->n[t], x, FLOOR((R)ths->n[t] * x + K(0.5)));
      const R seed = PHI(ths->n[t], dx / (R)ths->n[t], t);

      ths->psi[2 * (k * d + t)] = dx < K(0.0) ? -seed : seed;
      ths->psi[2 * (k * d + t) + 1] = EXP(K(2.0) * dx / ths->b[t]);
    }
  }
}

void X(precompute_psi)(X(plan) *ths)
{
  const INT d = ths->d, w = 2 * ths->m + 2;
  INT t, k;

  sort_precompute(ths);

  for (t = 0; t < d; t++)
  {
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < ths->M_total; k++)
    {
      const R x = ths->x[node(ths, k) * d + t];

      PHI_RUN(ths->psi + (k * d + t) * w, ths->n[t], x,
        grid_start(x, ths->n[t], ths->m), t);
    }
  }
}

/* For each node, in slot order, the (2m+2)^d products of window values and their offsets in
 * g, in the order of the tensor loops of node_trafo. */
void X(precompute_full_psi)(X(plan) *ths)
{
  const INT d = ths->d, w = 2 * ths->m + 2, lprod = full_psi_len(ths),
    gs = row_stride(ths);
  INT k;

  sort_precompute(ths);

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < ths->M_total; k++)
  {
    const INT j = node(ths, k);
    R buf[d * w];
    INT u[d], c[d], t, l;

    for (t = 0; t < d; t++)
    {
      const R x = ths->x[j * d + t];
      const INT s = grid_start(x, ths->n[t], ths->m);

      PHI_RUN(buf + t * w, ths->n[t], x, s, t);
      u[t] = s < 0 ? s + ths->n[t] : s;
      c[t] = 0;
    }

    for (l = 0; l < lprod; l++)
    {
      INT off = 0;
      R p = K(1.0);

      for (t = 0; t < d; t++)
      {
        const INT i = u[t] + c[t] < ths->n[t] ? u[t] + c[t]
          : u[t] + c[t] - ths->n[t];

        off = t < d - 1 ? off * ths->n[t] + i : off * gs + i;
        p *= buf[t * w + c[t]];
      }
      ths->psi[k * lprod + l] = p;
      ths->psi_index_g[k * lprod + l] = off;

      for (t = d - 1; t >= 0 && ++c[t] == w; t--)
        c[t] = 0;
    }

    ths->psi_index_f[k] = lprod;
  }
}

void X(precompute_one_psi)(X(plan) *ths)
{
  if (ths->flags & PRE_LIN_PSI)
    X(precompute_lin_psi)(ths);
  if (ths->flags & PRE_FG_PSI)
    X(precompute_fg_psi)(ths);
  if (ths->flags & PRE_PSI)
    X(precompute_psi)(ths);
  if (ths->flags & PRE_FULL_PSI)
    X(precompute_full_psi)(ths);
}

static void init_help(X(plan) *ths)
{
  const INT d = ths->d;
  INT t;

  if (ths->flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    ths->flags |= NFFT_SORT_NODES;

  ths->N_total = intprod(ths->N, d - 1) * (ths->N[d - 1] / 2 + 1);
  ths->n_total = intprod(ths->n, d);

  ths->sigma = (R*)Y(malloc)((size_t)d * sizeof(R));
  for (t = 0; t < d; t++)
    ths->sigma[t] = (R)ths->n[t] / (R)ths->N[t];

  WINDOW_HELP_INIT;

  if (ths->flags & MALLOC_X)
    ths->x = (R*)Y(malloc)((size_t)(d * ths->M_total) * sizeof(R));

  if (ths->flags & MALLOC_F_HAT)
    ths->f_hat = (C*)Y(malloc)((size_t)ths->N_total * sizeof(C));

  if (ths->flags & MALLOC_F)
    ths->f = (R*)Y(malloc)((size_t)ths->M_total * sizeof(R));

  if (ths->flags & PRE_PHI_HUT)
    precompute_phi_hut(ths);

  if (ths->flags & PRE_LIN_PSI)
  {
    if (ths->K == 0)
      ths->K = Y(m2K)(ths->m);
    ths->psi = (R*)Y(malloc)((size_t)((ths->K + 1) * d) * sizeof(R));
  }

  if (ths->flags & PRE_FG_PSI)
    ths->psi = (R*)Y(malloc)((size_t)(ths->M_total * d * 2) * sizeof(R));

  if (ths->flags & PRE_PSI)
    ths->psi = (R*)Y(malloc)((size_t)(ths->M_total * d * (2 * ths->m + 2))
      * sizeof(R));

  if (ths->flags & PRE_FULL_PSI)
  {
    const INT lprod = full_psi_len(ths);

    ths->psi = (R*)Y(malloc)((size_t)(ths->M_total * lprod) * sizeof(R));
    ths->psi_index_f = (INT*)Y(malloc)((size_t)ths->M_total * sizeof(INT));
    ths->psi_index_g = (INT*)Y(malloc)((size_t)(ths->M_total * lprod)
      * sizeof(INT));
  }

  if (ths->flags & FFTW_INIT)
  {
#ifdef _OPENMP
    INT nthreads = Y(get_num_threads)();
#endif

    ths->g_hat = (C*)Y(malloc)((size_t)(grid_rows(ths) * hat_len(ths))
      * sizeof(C));

    if (ths->flags & FFT_OUT_OF_PLACE)
      ths->g = (R*)Y(malloc)((size_t)ths->n_total * sizeof(R));
    else
      ths->g = (R*)ths->g_hat;

#if defined(_OPENMP) && defined(HAVE_FFTW_THREADS)
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
    FFTW(plan_with_nthreads)(nthreads);
#endif
    {
      int *_n = Y(malloc)((size_t)d * sizeof(int));

      for (t = 0; t < d; t++)
        _n[t] = (int)ths->n[t];

      ths->my_fftw_plan1 = FFTW(plan_dft_c2r)((int)d, _n, ths->g_hat, ths->g,
        ths->fftw_flags);
      ths->my_fftw_plan2 = FFTW(plan_dft_r2c)((int)d, _n, ths->g, ths->g_hat,
        ths->fftw_flags);
      Y(free)(_n);
    }
#if defined(_OPENMP) && defined(HAVE_FFTW_THREADS)
}
#endif
  }

  /* room for the sort scratch when the transforms sort */
  if (ths->flags & NFFT_SORT_NODES)
    ths->index_x = (INT*)Y(malloc)(sizeof(INT) * ((ths->flags & PRE_NODE_PSI)
      ? 2U : 4U) * (size_t)ths->M_total);
  else
    ths->index_x = NULL;

  ths->mv_trafo = (void (*)(void*))X(trafo);
  ths->mv_adjoint = (void (*)(void*))X(adjoint);
}

void X(init)(X(plan) *ths, int d, int *N, int M_total)
{
  INT t;

  ths->d = (INT)d;
  ths->M_total = (INT)M_total;

  ths->N = (INT*)Y(malloc)((size_t)d * sizeof(INT));
  ths->n = (INT*)Y(malloc)((size_t)d * sizeof(INT));
  for (t = 0; t < d; t++)
  {
    ths->N[t] = (INT)N[t];
    ths->n[t] = 2 * Y(next_power_of_2)(ths->N[t]);
  }

  ths->m = WINDOW_HELP_ESTIMATE_m;

  ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F
    | FFTW_INIT | NFFT_SORT_NODES;
#ifdef _OPENMP
  ths->flags |= NFFT_OMP_BLOCKWISE_ADJOINT;
#endif

  ths->fftw_flags = FFTW_ESTIMATE | FFTW_DESTROY_INPUT;

  ths->K = 0;
  init_help(ths);
}

void X(init_lin)(X(plan) *ths, int d, int *N, int M_total, int *n, int m,
  int K, unsigned flags, unsigned fftw_flags)
{
  INT t;

  ths->d = (INT)d;
  ths->M_total = (INT)M_total;

  ths->N = (INT*)Y(malloc)((size_t)d * sizeof(INT));
  ths->n = (INT*)Y(malloc)((size_t)d * sizeof(INT));
  for (t = 0; t < d; t++)
  {
    ths->N[t] = (INT)N[t];
    ths->n[t] = (INT)n[t];
  }

  ths->m = (INT)m;
  ths->flags = flags;
  ths->fftw_flags = fftw_flags;

  ths->K = K;
  init_help(ths);
}

void X(init_guru)(X(plan) *ths, int d, int *N, int M_total, int *n, int m,
  unsigned flags, unsigned fftw_flags)
{
  X(init_lin)(ths, d, N, M_total, n, m, 0, flags, fftw_flags);
}

void X(init_1d)(X(plan) *ths, int N1, int M_total)
{
  int N[1];

  N[0] = N1;
  X(init)(ths, 1, N, M_total);
}

void X(init_2d)(X(plan) *ths, int N1, int N2, int M_total)
{
  int N[2];

  N[0] = N1;
  N[1] = N2;
  X(init)(ths, 2, N, M_total);
}

void X(init_3d)(X(plan) *ths, int N1, int N2, int N3, int M_total)
{
  int N[3];

  N[0] = N1;
  N[1] = N2;
  N[2] = N3;
  X(init)(ths, 3, N, M_total);
}

const char* X(check)(X(plan) *ths)
{
  INT j;

  if (!ths->f)
    return "Member f not initialized.";

  if (!ths->x)
    return "Member x not initialized.";

  if (!ths->f_hat)
    return "Member f_hat not initialized.";

  if ((ths->flags & PRE_LIN_PSI) && ths->K < ths->M_total)
    return "Number of nodes too small to use PRE_LIN_PSI.";

  for (j = 0; j < ths->M_total * ths->d; j++)
    if (ths->x[j] < -K(0.5) || ths->x[j] >= K(0.5))
      return "ths->x out of range [-0.5,0.5)";

  for (j = 0; j < ths->d; j++)
  {
    if (ths->sigma[j] <= 1)
      return "Oversampling factor too small";

#ifdef WINDOW_IS_KAISER_BESSEL
    if (!KB_POLY_SIGMA_OK(ths->sigma[j]))
      return "Oversampling factor outside [1.1, 4]";
#endif

    if (ths->N[j] % 2 == 1 || ths->N[j] < 2)
      return "polynomial degree N has to be even and positive";
  }
  return 0;
}

void X(finalize)(X(plan) *ths)
{
  INT t;

  if (ths->flags & NFFT_SORT_NODES)
    Y(free)(ths->index_x);

  if (ths->flags & FFTW_INIT)
  {
#ifdef _OPENMP
    #pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
    FFTW(destroy_plan)(ths->my_fftw_plan2);
#ifdef _OPENMP
    #pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
    FFTW(destroy_plan)(ths->my_fftw_plan1);

    if (ths->flags & FFT_OUT_OF_PLACE)
      Y(free)(ths->g);

    Y(free)(ths->g_hat);
  }

  if (ths->flags & PRE_FULL_PSI)
  {
    Y(free)(ths->psi_index_g);
    Y(free)(ths->psi_index_f);
  }

  if (ths->flags & (PRE_FULL_PSI | PRE_PSI | PRE_FG_PSI | PRE_LIN_PSI))
    Y(free)(ths->psi);

  if (ths->flags & PRE_PHI_HUT)
  {
    for (t = 0; t < ths->d; t++)
      Y(free)(ths->c_phi_inv[t]);
    Y(free)(ths->c_phi_inv);
  }

  if (ths->flags & MALLOC_F)
    Y(free)(ths->f);

  if (ths->flags & MALLOC_F_HAT)
    Y(free)(ths->f_hat);

  if (ths->flags & MALLOC_X)
    Y(free)(ths->x);

  WINDOW_HELP_FINALIZE;

  Y(free)(ths->sigma);
  Y(free)(ths->n);
  Y(free)(ths->N);
}
