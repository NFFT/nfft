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

/* Nonequispaced FFT */

/* Authors: D. Potts, S. Kunis 2002-2009, Jens Keiner 2009, Toni Volkmer 2012 */

/* configure header */
#include "config.h"

/* complex datatype (maybe) */
#ifdef HAVE_COMPLEX_H
#include<complex.h>
#endif

/* NFFT headers */
#include "nfft3.h"
#include "infft.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef OMP_ASSERT
#include <assert.h>
#endif

#undef X
#define X(name) NFFT(name)

/** Compute aggregated product of integer array. */
static inline INT intprod(const INT *vec, const INT a, const INT d)
{
  INT t, p;

  p = 1;
  for (t = 0; t < d; t++)
    p *= vec[t] - a;

  return p;
}

/* handy shortcuts */
#define BASE(x) CEXP(x)

/** Index u wrapped into [0, n): one correction for u in [-n, 2n), the
 *  remainder else. */
static inline INT wrap_index(INT u, const INT n)
{
  if (u < 0)
    u += n;
  else if (u >= n)
    u -= n;

  return (u >= 0 && u < n) ? u : ((u % n) + n) % n;
}

/** Lower corner of the window of node x, wrapped into [0, n): the window
 *  covers the points uo_wrap, ..., uo_wrap+2m+1 of the padded grid of B. */
static inline INT uo_wrap(const R x, const INT n, const INT m)
{
  return wrap_index(LRINT(FLOOR(x * (R)(n))) - m, n);
}

/** Tiles of the grid of B (see Y(tile_layout)): the cells are the lower
 *  window corners uo_wrap in [0, n_t) on axis t. */
static INT tile_layout(const X(plan) *ths, INT *edge, INT *ntiles)
{
  return Y(tile_layout)(ths->d, 2 * ths->m + 2, ths->n, edge, ntiles);
}

/** Sort key of node i of ths->x: tile, then cell in the tile. */
static INT node_key(const X(plan) *ths, const INT i, const INT *edge,
    const INT *ntiles)
{
  const INT d = ths->d;
  INT c[d], t;

  for (t = 0; t < d; t++)
    c[t] = uo_wrap(ths->x[d * i + t], ths->n[t], ths->m);

  return Y(tile_key)(d, c, edge, ntiles);
}

/**
 * Sorts the nodes by the tile of the grid of B and, within a tile, by the
 * cell of the lower corner of their window, for better cache utilisation in
 * the multiplication with the matrix B and for the threaded B^T (see
 * Y(tiles_by_colour)). The result is written to ths->index_x as (key, node
 * index) pairs and the node range of each tile to ths->tile_ptr; the nodes
 * remain unchanged.
 */
static void sort0(X(plan) *ths)
{
  const INT d = ths->d, M = ths->M_total;
  INT edge[d], ntiles[d], nt, tilevol = 1, t, i, rhigh;
  INT *ar_x = ths->index_x, *ar_x_temp;

  nt = tile_layout(ths, edge, ntiles);
  for (t = 0; t < d; t++)
    tilevol *= edge[t];

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(i)
#endif
  for (i = 0; i < M; i++)
  {
    ar_x[2 * i] = node_key(ths, i, edge, ntiles);
    ar_x[2 * i + 1] = i;
  }

  rhigh = (INT) LRINT(CEIL(LOG2((R)(nt * tilevol)))) - 1;

  ar_x_temp = (INT*) Y(malloc)(2 * (size_t)(M) * sizeof(INT));
  Y(sort_node_indices_radix_lsdf)(M, ar_x, ar_x_temp, rhigh);
#ifdef OMP_ASSERT
  for (i = 1; i < M; i++)
    assert(ar_x[2 * (i - 1)] <= ar_x[2 * i]);
#endif
  Y(free)(ar_x_temp);

  /* nt tiles and the (empty) bucket of nodes outside the padded grid that
   * NFCT and NFST have */
  if (ths->tile_ptr == NULL)
    ths->tile_ptr = (INT*) Y(malloc)((size_t)(nt + 2) * sizeof(INT));
  Y(tile_ptr_from_keys)(M, ar_x, tilevol, nt + 1, ths->tile_ptr);
}

/**
 * Sorts the nodes for one transform, if flag NFFT_SORT_NODES is set and no
 * precomputation has taken a snapshot of them (see sort_nodes).
 *
 * \arg ths nfft_plan
 */
static inline void sort(X(plan) *ths)
{
  if ((ths->flags & NFFT_SORT_NODES) && ths->x_sorted == NULL)
    sort0(ths);
}

/**
 * Takes the snapshot of the nodes the transforms work on when flag
 * NFFT_SORT_NODES is set: sorts them into ths->index_x and copies them in that
 * order into ths->x_sorted. The psi tables are then built in the same order, so
 * that the matrix B is applied with sequential access to the node data and the
 * permutation is only used for ths->f. Without NFFT_SORT_NODES it does nothing.
 *
 * \arg ths nfft_plan
 */
static void sort_nodes(X(plan) *ths)
{
  const INT d = ths->d, M = ths->M_total;
  INT k;

  if (!(ths->flags & NFFT_SORT_NODES))
    return;

  sort0(ths);

  if (ths->x_sorted == NULL)
  {
    ths->x_sorted = (R*) Y(malloc)((size_t)(d * M) * sizeof(R));
    ths->f_sorted = (C*) Y(malloc)((size_t)(M) * sizeof(C));
  }

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    const R *xj = ths->x + ths->index_x[2 * k + 1] * d;
    INT t;
    for (t = 0; t < d; t++)
      ths->x_sorted[k * d + t] = xj[t];
  }
}

/** Nodes as the matrix B reads them: the sorted snapshot once it is taken,
 *  else the caller's array ths->x. */
static inline const R *node_x(const X(plan) *ths)
{
  return ths->x_sorted != NULL ? ths->x_sorted : ths->x;
}

/** Samples as the matrix B reads and writes them: the plan's scratch array in
 *  sorted order once the snapshot is taken, else the caller's array ths->f. */
static inline C *node_f(const X(plan) *ths)
{
  return ths->x_sorted != NULL ? ths->f_sorted : ths->f;
}

/** Row of the k-th node in processing order in the node data (node_x(ths),
 *  node_f(ths) and the psi tables). */
#define NODE_J(ths, k) (((ths)->x_sorted == NULL \
  && ((ths)->flags & NFFT_SORT_NODES)) ? (ths)->index_x[2 * (k) + 1] : (k))

/** Before the adjoint B: copy ths->f into the sorted scratch array. */
static void gather_f(X(plan) *ths)
{
  INT k;

  if (ths->x_sorted == NULL)
    return;

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < ths->M_total; k++)
    ths->f_sorted[k] = ths->f[ths->index_x[2 * k + 1]];
}

/** After the trafo B: copy the sorted scratch array back into ths->f. */
static void scatter_f(X(plan) *ths)
{
  INT k;

  if (ths->x_sorted == NULL)
    return;

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < ths->M_total; k++)
    ths->f[ths->index_x[2 * k + 1]] = ths->f_sorted[k];
}

/**
 * Halo of the grid of the matrix B: the window of a node covers the 2m+2 grid
 * points u, ..., u+2m+1 with the lower corner u wrapped into [0, n). The grid
 * is padded with this many points after the end of each dimension, so that
 * these indices never wrap around. The halo holds periodic copies of the
 * first points: filled from them before B (trafo), added back to them after
 * B^T (adjoint).
 */
#define HALO(m) (2 * (m) + 2)

/** Number of elements of the grid of ths padded with h points per dimension. */
static inline INT padded_total(const X(plan) *ths, const INT h)
{
  INT t, p = 1;

  for (t = 0; t < ths->d; t++)
    p *= ths->n[t] + h;

  return p;
}

/** Number of elements of the padded grid ths->g of the matrix B. */
static inline INT grid_B_total(const X(plan) *ths)
{
  return padded_total(ths, HALO(ths->m));
}

/** Number of elements of the FFTW grids ths->g1 and ths->g2. */
static inline INT grid_fftw_total(const X(plan) *ths)
{
  return padded_total(ths, ths->halo_hi);
}

/**
 * Copies the n_0 x ... x n_{d-1} grid from src, padded with hs points per
 * dimension, to dst, padded with hd points per dimension.
 */
static void copy_grid(const X(plan) *ths, C *dst, const INT hd, const C *src,
    const INT hs)
{
  const INT d = ths->d, nl = ths->n[d - 1];
  INT rows = 1, t, q;

  for (t = 0; t < d - 1; t++)
    rows *= ths->n[t];

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(q, t)
#endif
  for (q = 0; q < rows; q++)
  {
    INT r = q, od = 0, os = 0, sd = nl + hd, ss = nl + hs;

    for (t = d - 2; t >= 0; t--)
    {
      const INT i = r % ths->n[t];
      r /= ths->n[t];
      od += i * sd;
      os += i * ss;
      sd *= ths->n[t] + hd;
      ss *= ths->n[t] + hs;
    }

    memcpy(dst + od, src + os, (size_t)(nl) * sizeof(C));
  }
}

/** Length of the pieces of a row of the padded grid that fill_halo and
 *  fold_halo hand to one thread: they work on (outer row, piece) pairs, so
 *  that the first dimension, which has one outer row, runs in parallel too. */
#define HALO_BLOCK 1024

/**
 * Fills the halo of the padded grid gp periodically from its first points,
 * dimension by dimension: rows n_t, ..., n_t+h-1 of dimension t are copies of
 * rows 0, ..., h-1, over the whole padded range of the other dimensions, so
 * that corners come out right.
 */
static void fill_halo(const X(plan) *ths, C *gp)
{
  const INT d = ths->d, h = HALO(ths->m);
  INT t;

  for (t = 0; t < d; t++)
  {
    const INT nt = ths->n[t], pt = nt + h;
    INT outer = 1, inner = 1, nb, s, q;

    for (s = 0; s < t; s++)
      outer *= ths->n[s] + h;
    for (s = t + 1; s < d; s++)
      inner *= ths->n[s] + h;
    nb = (inner + HALO_BLOCK - 1) / HALO_BLOCK;

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(q)
#endif
    for (q = 0; q < outer * nb; q++)
    {
      const INT o = q / nb, i0 = (q - o * nb) * HALO_BLOCK;
      const INT i1 = MIN(i0 + HALO_BLOCK, inner);
      C *row = gp + o * pt * inner;
      INT r, i;

      for (r = 0; r < h; r++)
      {
        const C *src = row + (r < nt ? r : r % nt) * inner;
        C *dst = row + (nt + r) * inner;

        for (i = i0; i < i1; i++)
          dst[i] = src[i];
      }
    }
  }
}

/**
 * Adds the halo of the padded grid gp back to the points it is a periodic copy
 * of, dimension by dimension in reverse order (see fill_halo). Afterwards the
 * n_0 x ... x n_{d-1} part of gp holds the result.
 */
static void fold_halo(const X(plan) *ths, C *gp)
{
  const INT d = ths->d, h = HALO(ths->m);
  INT t;

  for (t = d - 1; t >= 0; t--)
  {
    const INT nt = ths->n[t], pt = nt + h;
    INT outer = 1, inner = 1, nb, s, q;

    for (s = 0; s < t; s++)
      outer *= ths->n[s] + h;
    for (s = t + 1; s < d; s++)
      inner *= ths->n[s] + h;
    nb = (inner + HALO_BLOCK - 1) / HALO_BLOCK;

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(q)
#endif
    for (q = 0; q < outer * nb; q++)
    {
      const INT o = q / nb, i0 = (q - o * nb) * HALO_BLOCK;
      const INT i1 = MIN(i0 + HALO_BLOCK, inner);
      C *row = gp + o * pt * inner;
      INT r, i;

      for (r = 0; r < h; r++)
      {
        C *dst = row + (r < nt ? r : r % nt) * inner;
        const C *src = row + (nt + r) * inner;

        for (i = i0; i < i1; i++)
          dst[i] += src[i];
      }
    }
  }
}

/** Points ths->g at the padded grid of the matrix B: ths->g2 when it is
 *  padded, else the plan's scratch grid ths->g_pad (allocated on first use and
 *  grown if n has grown). */
static void select_grid_B(X(plan) *ths)
{
  if (ths->halo_hi == HALO(ths->m))
  {
    ths->g = ths->g2;
    return;
  }

  {
    const INT size = grid_B_total(ths);
    if (ths->g_pad == NULL || ths->g_pad_size < size)
    {
      Y(free)(ths->g_pad);
      ths->g_pad = (C*) Y(malloc)((size_t)(size) * sizeof(C));
      ths->g_pad_size = size;
    }
  }
  ths->g = ths->g_pad;
}

/** After the FFT of the trafo: makes ths->g the padded grid of B, holding
 *  ths->g2 and its filled halo. */
static void grid_B_trafo(X(plan) *ths)
{
  select_grid_B(ths);
  if (ths->g != ths->g2)
    copy_grid(ths, ths->g, HALO(ths->m), ths->g2, ths->halo_hi);
  fill_halo(ths, ths->g);
}

/** After B^T of the adjoint: folds the halo of ths->g and leaves the result
 *  in ths->g2, the input of the FFT. */
static void grid_B_adjoint_done(X(plan) *ths)
{
  fold_halo(ths, ths->g);
  if (ths->g != ths->g2)
    copy_grid(ths, ths->g2, ths->halo_hi, ths->g, HALO(ths->m));
}

/** Sets the padded grid ths->g of B to zero, in parallel with OpenMP. */
static void zero_grid_B(const X(plan) *ths)
{
  Y(zero_reals_parallel)((R*)ths->g, 2 * grid_B_total(ths));
}

/** direct computation of non equispaced fourier transforms
 *  nfft_trafo_direct, ndft_conjugated, nfft_adjoint_direct, ndft_transposed
 *  require O(M_total N^d) arithemtical operations
 *
 * direct computation of the nfft_trafo_direct and ndft_conjugated, formula (1.1)
 * nfft_trafo_direct:
 * for j=0,...,M_total-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 * ndft_conjugated:
 * for j=0,...,M_total-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 *
 * direct computation of the nfft_adjoint_direct and ndft_transposed, formula (1.2)
 * nfft_adjoint_direct:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M_total-1} f[j] * exp(+2(pi) k x[j])
 * ndft_transposed:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M_total-1} f[j] * exp(-2(pi) k x[j])
 */
/* Block size for the phase recurrence in the direct transforms */
#define NFFT_DIRECT_RECURRENCE_BLOCK 32
#define TRAFO_STEP(i, w) (v += f_hat[i] * (w))
#define ADJOINT_STEP(i, w) (f_hat[i] += f[j] * (w))

/* Minimum innermost-dimension length for the multivariate recurrence to pay. */
#define NFFT_DIRECT_RECURRENCE_MIN_INNER 2

/* Accurate phase for exp(+-i 2pi k x): reduce k*x modulo 1 into ~[-1/2,1/2) so COS/SIN see a
 * small argument, error does not grow with N. Requires FMA single-rounding semantics. */
static inline R X(reduced_omega)(const R k, const R x)
{
  const R n = RINT(k * x);  // Nearest integer to k * x.
  return K2PI * FFMA(k, x, -n);  // Calculate k * x - n witha single rounding, then multiply 2 * pi.
}

/* One dimension's contribution to a multivariate phase: frac(k*x) in [-1/2,1/2). */
static inline R X(reduced_frac)(const R k, const R x)
{
  return FFMA(k, x, -RINT(k * x));
}

/* Accurate phase for exp(+-i 2pi sum_t k[t] x[t]), using frac(sum) == frac(sum of fracs) (mod 1).
 * Each x[t] is the raw node coordinate in [-1/2,1/2); accuracy ~d*u, independent of grid size.
 * For d==1 identical to reduced_omega. */
static inline R X(reduced_omega_multi)(const INT d, const INT *k, const R *x)
{
  R s = K(0.0);
  INT t;
  for (t = 0; t < d; t++)
    s += X(reduced_frac)((R)k[t], x[t]);
  return K2PI * (s - RINT(s));              /* reduce the sum mod 1, then scale by 2pi */
}

void X(trafo_direct)(const X(plan) *ths)
{
  C *f_hat = (C*)ths->f_hat, *f = (C*)ths->f;

  if (ths->d == 1)
  {
    /* specialize for univariate case, rationale: faster */
    const INT B = NFFT_DIRECT_RECURRENCE_BLOCK;
    INT j;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j)
#endif
    for (j = 0; j < ths->M_total; j++)
    {
      C v = K(0.0);
      const R x = ths->x[j];
      const R dphi = K2PI * x;                 /* |dphi| <= pi: accurate without reduction */
      const C dw = COS(dphi) - II * SIN(dphi); /* per-step phase factor exp(-i 2pi x)      */
      INT k_L = 0;
      while (k_L < ths->N_total)
      {
        /* Accurate seed exp(-i 2pi k x) at the block centre, then recur within the block. */
        const INT len = MIN(B, ths->N_total - k_L);
        const R omega = X(reduced_omega)((R)(k_L + len/2 - ths->N_total/2), x);
        const C w0 = COS(omega) - II * SIN(omega);
        NFFT_RECUR_BLOCK(k_L, len, w0, dw, TRAFO_STEP);
        k_L += len;
      }

      f[j] = v;
    }
  }
  else
  {
    /* multivariate case: recur the kernel phase along the contiguous innermost dimension
     * (length Nlast), re-seeded with an accurate phase at every row start and every B steps,
     * so the per-k_L COS/SIN of the 1d case is replaced by one complex multiply. The j loop
     * carries no dependency, so the same body serves serial and OpenMP (forward). */
    const INT B = NFFT_DIRECT_RECURRENCE_BLOCK;
    const INT Nlast = ths->N[ths->d - 1];
    const int use_rec = (Nlast >= NFFT_DIRECT_RECURRENCE_MIN_INNER);
    INT j;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j)
#endif
    for (j = 0; j < ths->M_total; j++)
    {
      C v = K(0.0);
      const R *xj = ths->x + (size_t)j * ths->d;          /* raw coordinates in [-1/2,1/2) */
      if (use_rec)
      {
        const R dphi = K2PI * xj[ths->d - 1];             /* innermost step, |dphi| <= pi  */
        const C dw = COS(dphi) - II * SIN(dphi);
        INT k[ths->d], t, k_L = 0;
        for (t = 0; t < ths->d; t++) k[t] = -ths->N[t]/2; /* multi-index, odometer (no div) */
        while (k_L < ths->N_total)                        /* one outer iteration per row    */
        {
          INT i = 0;
          while (i < Nlast)                               /* re-seed per row and per B steps */
          {
            const INT len = MIN(B, Nlast - i);
            k[ths->d - 1] = -Nlast/2 + i + len/2;
            const R omega = X(reduced_omega_multi)(ths->d, k, xj);
            const C w0 = COS(omega) - II * SIN(omega);
            NFFT_RECUR_BLOCK(k_L, len, w0, dw, TRAFO_STEP);
            i += len; k_L += len;
          }
          /* advance the outer multi-index (dims 0..d-2); innermost reset at next row's seed */
          for (t = ths->d - 2; (t >= 1) && (k[t] == ths->N[t]/2 - 1); t--) k[t] -= ths->N[t] - 1;
          k[t]++;
        }
      }
      else
      {
        /* short innermost dimension: per-k_L evaluation from the reduced-fraction prefix sums */
        R S[ths->d + 1];
        INT t, t2, k_L, k[ths->d];
        S[0] = K(0.0);
        for (t = 0; t < ths->d; t++)
        {
          k[t] = -ths->N[t]/2;
          S[t+1] = S[t] + X(reduced_frac)((R)k[t], xj[t]);
        }
        for (k_L = 0; k_L < ths->N_total; k_L++)
        {
          const R omega = K2PI * (S[ths->d] - RINT(S[ths->d]));
          v += f_hat[k_L] * (COS(omega) - II * SIN(omega));
          for (t = ths->d - 1; (t >= 1) && (k[t] == ths->N[t]/2 - 1); t--) k[t] -= ths->N[t] - 1;
          k[t]++;
          for (t2 = t; t2 < ths->d; t2++) S[t2+1] = S[t2] + X(reduced_frac)((R)k[t2], xj[t2]);
        }
      }

      f[j] = v;
    }
  }
}

/* Adjoint blocked recurrence over the row range [rlo,rhi); a row is the contiguous block of
 * N[d-1] frequencies the recurrence runs along. Callers must pass disjoint ranges -- that is
 * what keeps the f_hat accumulation race-free. */
static void X(adjoint_direct_rows)(const X(plan) *ths, const INT rlo, const INT rhi)
{
  C *f_hat = (C*)ths->f_hat, *f = (C*)ths->f;
  const INT B = NFFT_DIRECT_RECURRENCE_BLOCK;
  const INT Nlast = ths->N[ths->d - 1];
  INT k0[ths->d], t, j, r = rlo;

  if (rlo >= rhi)
    return;

  /* Decode the outer multi-index once; the odometer carries it, keeping division out of the
   * inner work. */
  k0[ths->d - 1] = -Nlast/2;
  for (t = ths->d - 2; t >= 0; t--)
  {
    k0[t] = r % ths->N[t] - ths->N[t]/2;
    r /= ths->N[t];
  }

  for (j = 0; j < ths->M_total; j++)
  {
    const R *xj = ths->x + (size_t)j * ths->d;  /* raw coordinates in [-1/2,1/2) */
    const R dphi = K2PI * xj[ths->d - 1];       /* innermost step, |dphi| <= pi */
    const C dw = COS(dphi) + II * SIN(dphi);
    INT k[ths->d], row, k_L = rlo * Nlast;
    for (t = 0; t < ths->d; t++) k[t] = k0[t];
    for (row = rlo; row < rhi; row++)
    {
      INT i = 0;
      while (i < Nlast)                         /* re-seed per row and per B steps */
      {
        const INT len = MIN(B, Nlast - i);
        k[ths->d - 1] = -Nlast/2 + i + len/2;
        const R omega = X(reduced_omega_multi)(ths->d, k, xj);
        const C w0 = COS(omega) + II * SIN(omega);
        NFFT_RECUR_BLOCK(k_L, len, w0, dw, ADJOINT_STEP);
        i += len; k_L += len;
      }
      /* advance the outer multi-index (dims 0..d-2); innermost reset at next row's seed */
      for (t = ths->d - 2; (t >= 1) && (k[t] == ths->N[t]/2 - 1); t--) k[t] -= ths->N[t] - 1;
      k[t]++;
    }
  }
}

void X(adjoint_direct)(const X(plan) *ths)
{
  C *f_hat = (C*)ths->f_hat, *f = (C*)ths->f;

  memset(f_hat, 0, (size_t)(ths->N_total) * sizeof(C));

  if (ths->d == 1)
  {
    /* specialize for univariate case, rationale: faster */
    const INT B = NFFT_DIRECT_RECURRENCE_BLOCK;
#ifdef _OPENMP
    if (ths->N_total > B)
    {
      /* Give each thread a disjoint, contiguous range of frequencies [klo,khi) (so the
       * f_hat[k] writes are race-free) and run the phase recurrence within it, re-seeded
       * every B steps. */
      #pragma omp parallel default(shared)
      {
        const int nt = omp_get_num_threads();
        const int tid = omp_get_thread_num();
        const INT klo = (INT)(((long long)ths->N_total * tid) / nt);
        const INT khi = (INT)(((long long)ths->N_total * (tid + 1)) / nt);
        INT j;
        for (j = 0; j < ths->M_total; j++)
        {
          const R x = ths->x[j];
          const R dphi = K2PI * x;
          const C dw = COS(dphi) + II * SIN(dphi);
          INT k_L = klo;
          while (k_L < khi)
          {
            const INT len = MIN(B, khi - k_L);
            const R omega = X(reduced_omega)((R)(k_L + len/2 - ths->N_total/2), x);
            const C w0 = COS(omega) + II * SIN(omega);
            NFFT_RECUR_BLOCK(k_L, len, w0, dw, ADJOINT_STEP);
            k_L += len;
          }
        }
      }
    }
    else
    {
      /* N <= B: the recurrence spans at most one block per thread-range, so its per-block
       * seed/setup costs more than it saves once threaded. Use the plain per-k
       * parallelisation. At these tiny N the per-entry phase error is in check. */
      INT k_L;
      #pragma omp parallel for default(shared) private(k_L)
      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        INT j;
        for (j = 0; j < ths->M_total; j++)
        {
          R omega = K2PI * ((R)(k_L - (ths->N_total/2))) * ths->x[j];
          f_hat[k_L] += f[j] * (COS(omega) + II * SIN(omega));
        }
      }
    }
#else
      INT j;
      for (j = 0; j < ths->M_total; j++)
      {
        const R x = ths->x[j];
        const R dphi = K2PI * x;
        const C dw = COS(dphi) + II * SIN(dphi);
        INT k_L = 0;
        while (k_L < ths->N_total)
        {
          const INT len = MIN(B, ths->N_total - k_L);
          const R omega = X(reduced_omega)((R)(k_L + len/2 - ths->N_total/2), x);
          const C w0 = COS(omega) + II * SIN(omega);
          NFFT_RECUR_BLOCK(k_L, len, w0, dw, ADJOINT_STEP);
          k_L += len;
        }
      }
#endif
  }
  else
  {
    /* multivariate case */
    const INT Nlast = ths->N[ths->d - 1];
    /* The recurrence spans a whole row, so the parallel unit is the row. Row count is not part
     * of the choice: the recurrence beats the per-k_L arm so heavily per unit of work that even
     * two rows across many threads stay ahead of an N_total-way naive split. */
    const INT nrows = ths->N_total / Nlast;

    if (Nlast >= NFFT_DIRECT_RECURRENCE_MIN_INNER)
    {
#ifdef _OPENMP
      #pragma omp parallel default(shared)
      {
        const int nt = omp_get_num_threads();
        const int tid = omp_get_thread_num();
        X(adjoint_direct_rows)(ths, (INT)(((long long)nrows * tid) / nt),
                                    (INT)(((long long)nrows * (tid + 1)) / nt));
      }
#else
      X(adjoint_direct_rows)(ths, 0, nrows);
#endif
    }
    else
    {
      /* Short row: still reduced per dimension, so both arms carry the same error. */
      INT j, k_L;
#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(j, k_L)
      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        INT k[ths->d], k_temp, t;

        k_temp = k_L;

        for (t = ths->d - 1; t >= 0; t--)
        {
          k[t] = k_temp % ths->N[t] - ths->N[t]/2;
          k_temp /= ths->N[t];
        }

        for (j = 0; j < ths->M_total; j++)
        {
          const R *xj = ths->x + (size_t)j * ths->d;
          const R omega = X(reduced_omega_multi)(ths->d, k, xj);
          f_hat[k_L] += f[j] * (COS(omega) + II * SIN(omega));
        }
      }
#else
      for (j = 0; j < ths->M_total; j++)
      {
        const R *xj = ths->x + (size_t)j * ths->d;
        R S[ths->d + 1];    /* prefix sums of the per-dimension reduced fractions */
        INT t, t2, k[ths->d];
        S[0] = K(0.0);
        for (t = 0; t < ths->d; t++)
        {
          k[t] = -ths->N[t]/2;
          S[t+1] = S[t] + X(reduced_frac)((R)k[t], xj[t]);
        }
        for (k_L = 0; k_L < ths->N_total; k_L++)
        {
          const R omega = K2PI * (S[ths->d] - RINT(S[ths->d]));
          f_hat[k_L] += f[j] * (COS(omega) + II * SIN(omega));
          for (t = ths->d - 1; (t >= 1) && (k[t] == ths->N[t]/2 - 1); t--) k[t] -= ths->N[t] - 1;
          k[t]++;
          for (t2 = t; t2 < ths->d; t2++) S[t2+1] = S[t2] + X(reduced_frac)((R)k[t2], xj[t2]);
        }
      }
#endif
    }
  }
}

/** fast computation of non-equispaced fourier transforms
 *  require O(N^d log(N) + M_total) arithmetical operations
 *
 * fast computation of the nfft_trafo and nfft_conjugated, formula (1.1)
 * nfft_trafo:
 * for j=0,...,M_total-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 * nfft_conjugated:
 * for j=0,...,M_total-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 *
 * direct computation of the nfft_adjoint and nfft_transposed, formula (1.2)
 * nfft_adjoint:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M_total-1} f[j] * exp(+2(pi) k x[j])
 * nfft_transposed:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M_total-1} f[j] * exp(-2(pi) k x[j])
 */

/** macros and small sub routines for the fast transforms
 */

/** computes 2m+2 indices for the matrix B, for node j of the node array xs
 */
static inline void uo(const X(plan) *ths, const R *xs, const INT j, INT *up,
  INT *op, const INT act_dim)
{
  const R xj = xs[j * ths->d + act_dim];
  INT c = LRINT(FLOOR(xj * (R)(ths->n[act_dim])));

  (*up) = c - (ths->m);
  (*op) = c + 1 + (ths->m);
}

static inline void uo2(INT *u, INT *o, const R x, const INT n, const INT m)
{
  INT c = LRINT(FLOOR(x * (R)(n)));

  *u = (c - m + n) % n;
  *o = (c + 1 + m + n) % n;
}

/* Helpers of the tensor product kernels Y(B_tensor_trafo) and
 * Y(B_tensor_adjoint) (infft.h) on the complex grid, 2 reals per point. */

/** wl[2l] = wl[2l+1] = p[l], l < L: weights in the row layout of the grid. */
static inline void dup_weights(R *restrict wl, const R *restrict p,
    const INT L)
{
  INT l;

  for (l = 0; l < L; l++)
    wl[2 * l] = wl[2 * l + 1] = p[l];
}

/** wf[2l] = p[l] Re(f), wf[2l+1] = p[l] Im(f), l < L: weights times the
 *  sample f in the row layout of the grid. */
static inline void scale_weights(R *restrict wf, const R *restrict p,
    const C f, const INT L)
{
  const R fr = CREAL(f), fi = CIMAG(f);
  INT l;

  for (l = 0; l < L; l++)
  {
    wf[2 * l] = p[l] * fr;
    wf[2 * l + 1] = p[l] * fi;
  }
}

static inline void zero_reals(R *a, const INT len)
{
  INT i;

  for (i = 0; i < len; i++)
    a[i] = K(0.0);
}

/** Sum of the L (re, im) pairs of acc. */
static inline C sum_pairs(const R *acc, const INT L)
{
  R re = K(0.0), im = K(0.0);
  INT l;

  for (l = 0; l < L; l++)
  {
    re += acc[2 * l];
    im += acc[2 * l + 1];
  }

  return re + II * im;
}

/** f_j = sum of psi[l] g[idx[l]], l < len (even): one node of B with
 *  PRE_FULL_PSI, into local accumulators. */
static inline C full_psi_trafo(const R *restrict psi, const INT *restrict idx,
    const C *g, const INT len)
{
  C a0 = K(0.0), a1 = K(0.0);
  INT l;

  for (l = 0; l < len; l += 2)
  {
    a0 += psi[l] * g[idx[l]];
    a1 += psi[l+1] * g[idx[l+1]];
  }

  return a0 + a1;
}

#define MACRO_D_compute_A \
{ \
  g_hat[k_plain[ths->d]] = f_hat[ks_plain[ths->d]] * c_phi_inv_k[ths->d]; \
}

#define MACRO_D_compute_T \
{ \
  f_hat[ks_plain[ths->d]] = g_hat[k_plain[ths->d]] * c_phi_inv_k[ths->d]; \
}

#define MACRO_D_init_result_A memset(g_hat, 0, (size_t)(grid_fftw_total(ths)) * sizeof(C));

#define MACRO_D_init_result_T memset(f_hat, 0, (size_t)(ths->N_total) * sizeof(C));

#define MACRO_with_PRE_PHI_HUT * ths->c_phi_inv[t2][ks[t2]];

#define MACRO_without_PRE_PHI_HUT / (PHI_HUT(ths->n[t2],ks[t2]-(ths->N[t2]/2),t2));

#define MACRO_init_k_ks \
{ \
  for (t = ths->d-1; 0 <= t; t--) \
  { \
    kp[t] = k[t] = 0; \
    ks[t] = ths->N[t]/2; \
  } \
  t++; \
}

#define MACRO_update_c_phi_inv_k(which_one) \
{ \
  for (t2 = t; t2 < ths->d; t2++) \
  { \
    c_phi_inv_k[t2+1] = c_phi_inv_k[t2] MACRO_ ##which_one; \
    ks_plain[t2+1] = ks_plain[t2]*ths->N[t2] + ks[t2]; \
    k_plain[t2+1] = k_plain[t2]*(ths->n[t2] + ths->halo_hi) + k[t2]; \
  } \
}

#define MACRO_count_k_ks \
{ \
  for (t = ths->d-1; (t > 0) && (kp[t] == ths->N[t]-1); t--) \
  { \
    kp[t] = k[t] = 0; \
    ks[t]= ths->N[t]/2; \
  } \
\
  kp[t]++; k[t]++; ks[t]++; \
  if(kp[t] == ths->N[t]/2) \
  { \
    k[t] = ths->n[t] - ths->N[t]/2; \
    ks[t] = 0; \
  } \
} \

/* sub routines for the fast transforms  matrix vector multiplication with D, D^T */
#define MACRO_D(which_one) \
static inline void D_serial_ ## which_one (X(plan) *ths) \
{ \
  C *f_hat, *g_hat; /* local copy */ \
  R c_phi_inv_k[ths->d+1]; /* postfix product of PHI_HUT */ \
  INT t, t2; /* index dimensions */ \
  INT k_L; /* plain index */ \
  INT kp[ths->d]; /* multi index (simple) */ \
  INT k[ths->d]; /* multi index in g_hat */ \
  INT ks[ths->d]; /* multi index in f_hat, c_phi_inv*/ \
  INT k_plain[ths->d+1]; /* postfix plain index */ \
  INT ks_plain[ths->d+1]; /* postfix plain index */ \
 \
  f_hat = (C*)ths->f_hat; g_hat = (C*)ths->g_hat; \
  MACRO_D_init_result_ ## which_one; \
\
  c_phi_inv_k[0] = K(1.0); \
  k_plain[0] = 0; \
  ks_plain[0] = 0; \
\
  MACRO_init_k_ks; \
\
  if (ths->flags & PRE_PHI_HUT) \
  { \
    for (k_L = 0; k_L < ths->N_total; k_L++) \
    { \
      MACRO_update_c_phi_inv_k(with_PRE_PHI_HUT); \
      MACRO_D_compute_ ## which_one; \
      MACRO_count_k_ks; \
    } \
  } \
  else \
  { \
    for (k_L = 0; k_L < ths->N_total; k_L++) \
    { \
      MACRO_update_c_phi_inv_k(without_PRE_PHI_HUT); \
      MACRO_D_compute_ ## which_one; \
      MACRO_count_k_ks; \
    } \
  } \
}

#ifdef _OPENMP
static inline void D_openmp_A(X(plan) *ths)
{
  C *f_hat, *g_hat;                     /**< local copy                     */
  INT k_L;                              /**< plain index                    */

  f_hat = (C*)ths->f_hat; g_hat = (C*)ths->g_hat;
  memset(g_hat, 0, (size_t)(grid_fftw_total(ths)) * sizeof(C));

  if (ths->flags & PRE_PHI_HUT)
  {
    #pragma omp parallel for default(shared) private(k_L)
    for (k_L = 0; k_L < ths->N_total; k_L++)
    {
      INT kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      INT k[ths->d];                        /**< multi index in g_hat           */
      INT ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      R c_phi_inv_k_val = K(1.0);
      INT k_plain_val = 0;
      INT ks_plain_val = 0;
      INT t;
      INT k_temp = k_L;

      for (t = ths->d-1; t >= 0; t--)
      {
        kp[t] = k_temp % ths->N[t];
        if (kp[t] >= ths->N[t]/2)
          k[t] = ths->n[t] - ths->N[t] + kp[t];
        else
          k[t] = kp[t];
        ks[t] = (kp[t] + ths->N[t]/2) % ths->N[t];
        k_temp /= ths->N[t];
      }

      for (t = 0; t < ths->d; t++)
      {
        c_phi_inv_k_val *= ths->c_phi_inv[t][ks[t]];
        ks_plain_val = ks_plain_val*ths->N[t] + ks[t];
        k_plain_val = k_plain_val*(ths->n[t] + ths->halo_hi) + k[t];
      }

      g_hat[k_plain_val] = f_hat[ks_plain_val] * c_phi_inv_k_val;
    } /* for(k_L) */
  } /* if(PRE_PHI_HUT) */
  else
  {
    #pragma omp parallel for default(shared) private(k_L)
    for (k_L = 0; k_L < ths->N_total; k_L++)
    {
      INT kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      INT k[ths->d];                        /**< multi index in g_hat           */
      INT ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      R c_phi_inv_k_val = K(1.0);
      INT k_plain_val = 0;
      INT ks_plain_val = 0;
      INT t;
      INT k_temp = k_L;

      for (t = ths->d-1; t >= 0; t--)
      {
        kp[t] = k_temp % ths->N[t];
        if (kp[t] >= ths->N[t]/2)
          k[t] = ths->n[t] - ths->N[t] + kp[t];
        else
          k[t] = kp[t];
        ks[t] = (kp[t] + ths->N[t]/2) % ths->N[t];
        k_temp /= ths->N[t];
      }

      for (t = 0; t < ths->d; t++)
      {
        c_phi_inv_k_val /= (PHI_HUT(ths->n[t],ks[t]-(ths->N[t]/2),t));
        ks_plain_val = ks_plain_val*ths->N[t] + ks[t];
        k_plain_val = k_plain_val*(ths->n[t] + ths->halo_hi) + k[t];
      }

      g_hat[k_plain_val] = f_hat[ks_plain_val] * c_phi_inv_k_val;
    } /* for(k_L) */
  } /* else(PRE_PHI_HUT) */
}
#endif

#ifndef _OPENMP
MACRO_D(A)
#endif

static inline void D_A(X(plan) *ths)
{
#ifdef _OPENMP
  D_openmp_A(ths);
#else
  D_serial_A(ths);
#endif
}

#ifdef _OPENMP
static void D_openmp_T(X(plan) *ths)
{
  C *f_hat, *g_hat;                     /**< local copy                     */
  INT k_L;                              /**< plain index                    */

  f_hat = (C*)ths->f_hat; g_hat = (C*)ths->g_hat;
  memset(f_hat, 0, ths->N_total * sizeof(C));

  if (ths->flags & PRE_PHI_HUT)
  {
    #pragma omp parallel for default(shared) private(k_L)
    for (k_L = 0; k_L < ths->N_total; k_L++)
    {
      INT kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      INT k[ths->d];                        /**< multi index in g_hat           */
      INT ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      R c_phi_inv_k_val = K(1.0);
      INT k_plain_val = 0;
      INT ks_plain_val = 0;
      INT t;
      INT k_temp = k_L;

      for (t = ths->d - 1; t >= 0; t--)
      {
        kp[t] = k_temp % ths->N[t];
        if (kp[t] >= ths->N[t]/2)
          k[t] = ths->n[t] - ths->N[t] + kp[t];
        else
          k[t] = kp[t];
        ks[t] = (kp[t] + ths->N[t]/2) % ths->N[t];
        k_temp /= ths->N[t];
      }

      for (t = 0; t < ths->d; t++)
      {
        c_phi_inv_k_val *= ths->c_phi_inv[t][ks[t]];
        ks_plain_val = ks_plain_val*ths->N[t] + ks[t];
        k_plain_val = k_plain_val*(ths->n[t] + ths->halo_hi) + k[t];
      }

      f_hat[ks_plain_val] = g_hat[k_plain_val] * c_phi_inv_k_val;
    } /* for(k_L) */
  } /* if(PRE_PHI_HUT) */
  else
  {
    #pragma omp parallel for default(shared) private(k_L)
    for (k_L = 0; k_L < ths->N_total; k_L++)
    {
      INT kp[ths->d];                       /**< multi index (simple)           */ //0..N-1
      INT k[ths->d];                        /**< multi index in g_hat           */
      INT ks[ths->d];                       /**< multi index in f_hat, c_phi_inv*/
      R c_phi_inv_k_val = K(1.0);
      INT k_plain_val = 0;
      INT ks_plain_val = 0;
      INT t;
      INT k_temp = k_L;

      for (t = ths->d-1; t >= 0; t--)
      {
        kp[t] = k_temp % ths->N[t];
        if (kp[t] >= ths->N[t]/2)
          k[t] = ths->n[t] - ths->N[t] + kp[t];
        else
          k[t] = kp[t];
        ks[t] = (kp[t] + ths->N[t]/2) % ths->N[t];
        k_temp /= ths->N[t];
      }

      for (t = 0; t < ths->d; t++)
      {
        c_phi_inv_k_val /= (PHI_HUT(ths->n[t],ks[t]-(ths->N[t]/2),t));
        ks_plain_val = ks_plain_val*ths->N[t] + ks[t];
        k_plain_val = k_plain_val*(ths->n[t] + ths->halo_hi) + k[t];
      }

      f_hat[ks_plain_val] = g_hat[k_plain_val] * c_phi_inv_k_val;
    } /* for(k_L) */
  } /* else(PRE_PHI_HUT) */
}
#endif

#ifndef _OPENMP
MACRO_D(T)
#endif

static void D_T(X(plan) *ths)
{
#ifdef _OPENMP
  D_openmp_T(ths);
#else
  D_serial_T(ths);
#endif
}

/* sub routines for the fast transforms matrix vector multiplication with B, B^T */

#define MACRO_with_FG_PSI fg_p[t2][lj[t2]]

#define MACRO_with_LIN_PSI fg_psi[t2][lj[t2]]

#define MACRO_with_PRE_PSI ths->psi[(j*ths->d+t2) * (2*ths->m+2)+lj[t2]]

#define MACRO_without_PRE_PSI_improved psij_const[t2 * (2*ths->m+2) + lj[t2]]

#define MACRO_without_PRE_PSI  PHI(ths->n[t2], \
  NX_SUB(ths->n[t2], xs[j*ths->d+t2], lj[t2]+u[t2])/((R)ths->n[t2]), t2)

#define MACRO_init_uo_l_lj_t \
INT l_all[ths->d*(2*ths->m+2)]; \
{ \
  for (t = ths->d-1; t >= 0; t--) \
  { \
    uo(ths, xs, j,&u[t],&o[t],t); \
    INT lj_t; \
    const INT uw_t = ((u[t] % ths->n[t]) + ths->n[t]) % ths->n[t]; \
    for (lj_t = 0; lj_t < 2*ths->m+2; lj_t++) \
      l_all[t*(2*ths->m+2) + lj_t] = uw_t + lj_t; \
    lj[t] = 0; \
  } \
  t++; \
}

#define MACRO_update_phi_prod_ll_plain(which_one) { \
  for (t2 = t; t2 < ths->d; t2++) \
    { \
      phi_prod[t2+1] = phi_prod[t2] * MACRO_ ## which_one; \
      ll_plain[t2+1] = ll_plain[t2] * (ths->n[t2] + HALO(ths->m)) + l_all[t2*(2*ths->m+2) + lj[t2]]; \
    } \
}

#define MACRO_count_uo_l_lj_t \
{ \
  for (t = ths->d-1; (t > 0) && (lj[t] == o[t]-u[t]); t--) \
  { \
    lj[t] = 0; \
  } \
 \
  lj[t]++; \
}

#define MACRO_COMPUTE_with_PRE_PSI MACRO_with_PRE_PSI
#define MACRO_COMPUTE_with_PRE_FG_PSI MACRO_with_FG_PSI
#define MACRO_COMPUTE_with_FG_PSI MACRO_with_FG_PSI
#define MACRO_COMPUTE_with_PRE_LIN_PSI MACRO_with_LIN_PSI
#define MACRO_COMPUTE_with_LIN_PSI MACRO_with_LIN_PSI
#define MACRO_COMPUTE_without_PRE_PSI MACRO_without_PRE_PSI_improved
#define MACRO_COMPUTE_without_PRE_PSI_improved MACRO_without_PRE_PSI_improved

/**
 * Window of node j for the psi mode of the plan, any d: the weights w[t][l],
 * l < 2m+2, on axis t, and the first window index lo[t] of axis t in the
 * padded grid of B. wbuf holds d rows of 2m+4 reals for the modes that form
 * the weights on the fly.
 */
static void node_window(const X(plan) *ths, const R *xs, const INT j,
    const R *fg_e, const R *fg_q, R *wbuf, const R **w, INT *lo)
{
  const INT d = ths->d, m = ths->m, L = 2 * m + 2;
  INT t;

  for (t = 0; t < d; t++)
  {
    R *b = wbuf + t * (L + 2);
    INT u, o;

    uo(ths, xs, j, &u, &o, t);
    lo[t] = wrap_index(u, ths->n[t]);

    if (ths->flags & PRE_PSI)
      w[t] = ths->psi + (j * d + t) * L;
    else if (ths->flags & PRE_FG_PSI)
    {
      FG_RUN(b, m, FABS(ths->psi[2 * (j * d + t)]),
          ths->psi[2 * (j * d + t) + 1], fg_e[t], fg_q[t]);
      w[t] = FG_RUN_START(b, ths->psi[2 * (j * d + t)]);
    }
    else if (ths->flags & FG_PSI)
    {
      const R fg_nx = (R)(ths->n[t]) * xs[j * d + t];
      const R fg_d = NX_SUB(ths->n[t], xs[j * d + t], FLOOR(fg_nx + K(0.5)));

      FG_RUN(b, m, (PHI(ths->n[t], fg_d / ((R)(ths->n[t])), t)),
          EXP(K(2.0) * fg_d / ths->b[t]), fg_e[t], fg_q[t]);
      w[t] = FG_RUN_START(b, fg_d);
    }
    else if (ths->flags & PRE_LIN_PSI)
    {
      /* the steps of the taps and of the node in the table alike, as in the
       * kernels for d = 1, 2, 3 */
      const INT ip_s = ths->K / (m + 2);
      const R y = NX_SUB(ths->n[t], xs[j * d + t], u) * ((R)ip_s);
      const INT ip_u = LRINT(FLOOR(y));
      const R ip_w = y - ip_u;
      INT l;

      for (l = 0; l < L; l++)
        b[l] = ths->psi[(ths->K + 1) * t + ABS(ip_u - l * ip_s)] * (1 - ip_w)
            + ths->psi[(ths->K + 1) * t + ABS(ip_u - l * ip_s + 1)] * (ip_w);
      w[t] = b;
    }
    else
    {
      PHI_RUN(b, ths->n[t], xs[j * d + t], u, t);
      w[t] = b;
    }
  }
}

/* B and B^T for node k of the node data, any d, with the tensor product
 * kernels of infft.h. Each in a function of its own so that the OpenMP region
 * of B_A holds no run-time sized arrays (see nfft_trafo_1d_B). */
static void B_node_A(const X(plan) *ths, const R *xs, C *fs, const INT *np,
    const R *fg_e, const R *fg_q, const INT k)
{
  const INT d = ths->d, L = 2 * ths->m + 2;
  const INT j = NODE_J(ths, k);
  R wbuf[d * (L + 2)], wl[2 * L], acc[2 * L];
  const R *w[d];
  INT lo[d];

  node_window(ths, xs, j, fg_e, fg_q, wbuf, w, lo);
  dup_weights(wl, w[d - 1], L);
  zero_reals(acc, 2 * L);
  Y(B_tensor_trafo)(d, L, 2, np, lo, w, wl, (const R*)ths->g, acc);
  fs[j] = sum_pairs(acc, L);
}

#ifndef _OPENMP
static void B_node_T(const X(plan) *ths, const R *xs, const C *fs,
    const INT *np, const R *fg_e, const R *fg_q, const INT k)
{
  const INT d = ths->d, L = 2 * ths->m + 2;
  const INT j = NODE_J(ths, k);
  R wbuf[d * (L + 2)], wf[2 * L];
  const R *w[d];
  INT lo[d];

  node_window(ths, xs, j, fg_e, fg_q, wbuf, w, lo);
  scale_weights(wf, w[d - 1], fs[j], L);
  Y(B_tensor_adjoint)(d, L, 2, np, lo, w, wf, (R*)ths->g);
}
#endif

/** Sizes of the padded grid of B and the factors of the FG modes per axis. */
static void B_init(const X(plan) *ths, INT *np, R *fg_e, R *fg_q)
{
  INT t;

  for (t = 0; t < ths->d; t++)
  {
    np[t] = ths->n[t] + HALO(ths->m);
    if (ths->flags & (PRE_FG_PSI | FG_PSI))
    {
      fg_e[t] = EXP(K(-1.0) / ths->b[t]);
      fg_q[t] = fg_e[t] * fg_e[t];
    }
  }
}

/** B for d >= 4: f_j = sum over the window of node j of psi * g. */
static void B_A(X(plan) *ths)
{
  const R *xs = node_x(ths);
  C *fs = node_f(ths);
  INT np[ths->d];
  R fg_e[ths->d], fg_q[ths->d];
  INT k;

  if (ths->flags & PRE_FULL_PSI)
  {
    INT lprod = 1;

    for (k = 0; k < ths->d; k++)
      lprod *= (2 * ths->m + 2);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < ths->M_total; k++)
    {
      const INT j = NODE_J(ths, k);
      fs[j] = full_psi_trafo(ths->psi + j * lprod, ths->psi_index_g + j * lprod,
          ths->g, lprod);
    }
    return;
  }

  sort(ths);
  B_init(ths, np, fg_e, fg_q);

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < ths->M_total; k++)
    B_node_A(ths, xs, fs, np, fg_e, fg_q, k);
}

#ifndef _OPENMP
/** B^T for d >= 4, serial: g = sum over the nodes j of psi * f_j on the
 *  window of node j. */
static void B_serial_T(X(plan) *ths)
{
  const R *xs = node_x(ths);
  const C *fs = node_f(ths);
  INT np[ths->d];
  R fg_e[ths->d], fg_q[ths->d];
  INT k;

  zero_grid_B(ths);

  if (ths->flags & PRE_FULL_PSI)
  {
    INT lprod = 1, l;

    for (k = 0; k < ths->d; k++)
      lprod *= (2 * ths->m + 2);

    for (k = 0; k < ths->M_total; k++)
    {
      const R *psi = ths->psi + k * lprod;
      const INT *idx = ths->psi_index_g + k * lprod;

      for (l = 0; l < lprod; l++)
        ths->g[idx[l]] += psi[l] * fs[k];
    }
    return;
  }

  sort(ths);
  B_init(ths, np, fg_e, fg_q);

  for (k = 0; k < ths->M_total; k++)
    B_node_T(ths, xs, fs, np, fg_e, fg_q, k);
}
#endif

/**
 * Calculates adjoint NFFT for flag PRE_FULL_PSI without sorted nodes.
 * Parallel calculation (OpenMP) with atomic operations.
 *
 * \arg lprod stride (2*m+2)^d
 *
 * \author Toni Volkmer
 */
static void nfft_adjoint_B_compute_full_psi(C *g, const INT *psi_index_g,
    const R *psi, const C *f, const INT M, const INT d, const INT *n,
    const INT m, const unsigned flags, const INT *index_x, const int sorted)
{
  INT k;
  INT lprod;

  UNUSED(n);
  {
    INT t;
    for(t = 0, lprod = 1; t < d; t++)
        lprod *= 2 * m + 2;
  }

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    INT l;
    INT j = (!sorted && (flags & NFFT_SORT_NODES)) ? index_x[2*k+1] : k;

    for (l = 0; l < lprod; l++)
    {
#ifdef _OPENMP
      C val = psi[j * lprod + l] * f[j];
      C *gref = g + psi_index_g[j * lprod + l];
      R *gref_real = (R*) gref;

      #pragma omp atomic
      gref_real[0] += CREAL(val);

      #pragma omp atomic
      gref_real[1] += CIMAG(val);
#else
      g[psi_index_g[j * lprod + l]] += psi[j * lprod + l] * f[j];
#endif
    }
  }
}



#ifdef _OPENMP
/* What the nodes of one tile of the threaded B^T need (adjoint_B_tile). */
typedef struct
{
  X(plan) *ths;
  const R *xs; /* node data, as node_x(ths) */
  const C *fs; /* samples, as node_f(ths) */
  const INT *np; /* sizes of the padded grid of B */
  const R *fg_e, *fg_q; /* factors of the FG modes per axis */
  INT lprod; /* (2m+2)^d, the PRE_FULL_PSI entries per node */
} adjoint_tile_ctx;

/** B^T for the sorted nodes k0, ..., k1-1 of one tile, straight into the
 *  padded grid of B: no other thread writes to their windows meanwhile (see
 *  Y(tiles_by_colour)). */
static void adjoint_B_tile(void *p, const INT k0, const INT k1)
{
  const adjoint_tile_ctx *c = (const adjoint_tile_ctx*)p;
  const X(plan) *ths = c->ths;
  const INT d = ths->d, L = 2 * ths->m + 2;
  C *g = ths->g;
  INT k;

  if (ths->flags & PRE_FULL_PSI)
  {
    for (k = k0; k < k1; k++)
    {
      const INT j = NODE_J(ths, k);
      const R *psi = ths->psi + j * c->lprod;
      const INT *idx = ths->psi_index_g + j * c->lprod;
      const C f = c->fs[j];
      INT l;

      for (l = 0; l < c->lprod; l++)
        g[idx[l]] += psi[l] * f;
    }
    return;
  }

  {
    R wbuf[d * (L + 2)], wf[2 * L];
    const R *w[d];
    INT lo[d];

    for (k = k0; k < k1; k++)
    {
      const INT j = NODE_J(ths, k);

      node_window(ths, c->xs, j, c->fg_e, c->fg_q, wbuf, w, lo);
      scale_weights(wf, w[d - 1], c->fs[j], L);
      Y(B_tensor_adjoint)(d, L, 2, c->np, lo, w, wf, (R*)g);
    }
  }
}

/**
 * Threaded B^T for sorted nodes (flag NFFT_SORT_NODES, also set by
 * NFFT_OMP_BLOCKWISE_ADJOINT), any d and psi mode: the tiles of the grid are
 * spread in colour phases, the threads of a phase taking whole tiles, without
 * atomics (see Y(tiles_by_colour)). The result does not depend on the number
 * of threads. ths->g must be zero. Returns 0 and does nothing without sorted
 * nodes.
 */
static int adjoint_B_tiled(X(plan) *ths)
{
  const INT d = ths->d;
  INT np[d], edge[d], ntiles[d], t;
  R fg_e[d], fg_q[d];
  adjoint_tile_ctx c;

  if (!(ths->flags & NFFT_SORT_NODES))
    return 0;

  sort(ths);
  B_init(ths, np, fg_e, fg_q);
  tile_layout(ths, edge, ntiles);

  c.ths = ths;
  c.xs = node_x(ths);
  c.fs = node_f(ths);
  c.np = np;
  c.fg_e = fg_e;
  c.fg_q = fg_q;
  for (t = 0, c.lprod = 1; t < d; t++)
    c.lprod *= 2 * ths->m + 2;

  Y(tiles_by_colour)(d, ntiles, ths->tile_ptr, adjoint_B_tile, &c);
  return 1;
}

#define MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_with_PRE_PSI
#define MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_with_PRE_PSI \
  MACRO_update_phi_prod_ll_plain(with_PRE_PSI);

#define MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_with_PRE_FG_PSI \
      R fg_psi[ths->d][2*ths->m+3]; \
      const R *fg_p[ths->d]; \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        FG_RUN(fg_psi[t2], ths->m, FABS(ths->psi[2*(j*ths->d+t2)]), \
            ths->psi[2*(j*ths->d+t2)+1], fg_e[t2], fg_q[t2]); \
        fg_p[t2] = FG_RUN_START(fg_psi[t2], ths->psi[2*(j*ths->d+t2)]); \
      }
#define MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_with_PRE_FG_PSI \
  MACRO_update_phi_prod_ll_plain(with_FG_PSI);

#define MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_with_FG_PSI \
      R fg_psi[ths->d][2*ths->m+3]; \
      const R *fg_p[ths->d]; \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        const R fg_nx = (R)(ths->n[t2]) * xs[j*ths->d+t2]; \
        const R fg_d = \
            NX_SUB(ths->n[t2], xs[j*ths->d+t2], FLOOR(fg_nx + K(0.5))); \
        FG_RUN(fg_psi[t2], ths->m, \
            (PHI(ths->n[t2], fg_d / ((R)(ths->n[t2])), t2)), \
            EXP(K(2.0) * fg_d / ths->b[t2]), \
            fg_e[t2], fg_q[t2]); \
        fg_p[t2] = FG_RUN_START(fg_psi[t2], fg_d); \
      }
#define MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_with_FG_PSI \
  MACRO_update_phi_prod_ll_plain(with_FG_PSI);

#define MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_with_PRE_LIN_PSI \
      R y[ths->d]; \
      R fg_psi[ths->d][2*ths->m+2]; \
      INT l_fg,lj_fg; \
      R ip_w; \
      INT ip_u; \
      INT ip_s = ths->K/(ths->m+2); \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        y[t2] = ((NX_SUB(ths->n[t2], xs[j*ths->d+t2], u[t2])) \
          * ((R)ths->K))/((R)ths->m+2); \
        ip_u  = LRINT(FLOOR(y[t2])); \
        ip_w  = y[t2]-ip_u; \
        for (l_fg = u[t2], lj_fg = 0; l_fg <= o[t2]; l_fg++, lj_fg++) \
        { \
          fg_psi[t2][lj_fg] = ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s)] \
            * (1-ip_w) + ths->psi[(ths->K+1)*t2 + ABS(ip_u-lj_fg*ip_s+1)] \
            * (ip_w); \
        } \
      }
#define MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_with_PRE_LIN_PSI \
  MACRO_update_phi_prod_ll_plain(with_LIN_PSI);

#define MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_without_PRE_PSI \
      R psij_const[ths->d * (2*ths->m+2)]; \
      for (t2 = 0; t2 < ths->d; t2++) \
      { \
        PHI_RUN(psij_const + t2 * (2*ths->m+2), ths->n[t2], \
            xs[j*ths->d+t2], u[t2], t2); \
      }
#define MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_without_PRE_PSI \
  MACRO_update_phi_prod_ll_plain(without_PRE_PSI_improved);

#define MACRO_adjoint_nd_B_OMP_COMPUTE(whichone) \
{ \
      INT u[ths->d], o[ths->d]; /* multi band with respect to x_j */ \
      INT l_L; /* index one row of B */ \
      INT lj[ths->d]; /* multi index 0<=lj<u+o+1 */ \
      INT ll_plain[ths->d+1]; /* postfix plain index in g */ \
      R phi_prod[ths->d+1]; /* postfix product of PHI */ \
 \
      phi_prod[0] = K(1.0); \
      ll_plain[0] = 0; \
 \
      MACRO_init_uo_l_lj_t; \
 \
      MACRO_adjoint_nd_B_OMP_COMPUTE_BEFORE_LOOP_ ## whichone \
 \
      if (ths->d == 4) \
      { \
        INT l0, l1, l2, l3; \
        for (l0 = 0; l0 < 2*ths->m+2; l0++) \
        { \
          lj[0] = l0; \
          t2 = 0; \
          phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
          ll_plain[t2+1] = ll_plain[t2] * (ths->n[t2] + HALO(ths->m)) + l_all[t2*(2*ths->m+2) + lj[t2]]; \
          for (l1 = 0; l1 < 2*ths->m+2; l1++) \
          { \
            lj[1] = l1; \
            t2 = 1; \
            phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
            ll_plain[t2+1] = ll_plain[t2] * (ths->n[t2] + HALO(ths->m)) + l_all[t2*(2*ths->m+2) + lj[t2]]; \
            for (l2 = 0; l2 < 2*ths->m+2; l2++) \
            { \
              lj[2] = l2; \
              t2 = 2; \
              phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
              ll_plain[t2+1] = ll_plain[t2] * (ths->n[t2] + HALO(ths->m)) + l_all[t2*(2*ths->m+2) + lj[t2]]; \
              for (l3 = 0; l3 < 2*ths->m+2; l3++) \
              { \
                lj[3] = l3; \
                t2 = 3; \
                phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                ll_plain[t2+1] = ll_plain[t2] * (ths->n[t2] + HALO(ths->m)) + l_all[t2*(2*ths->m+2) + lj[t2]]; \
 \
                C *lhs = ths->g + ll_plain[ths->d]; \
                R *lhs_real = (R*)lhs; \
                C val = phi_prod[ths->d] * fs[j]; \
 \
                _Pragma("omp atomic") \
                lhs_real[0] += CREAL(val); \
 \
                _Pragma("omp atomic") \
                lhs_real[1] += CIMAG(val); \
              } \
            } \
          } \
        } \
      } /* if(d==4) */ \
      else if (ths->d == 5) \
      { \
        INT l0, l1, l2, l3, l4; \
        for (l0 = 0; l0 < 2*ths->m+2; l0++) \
        { \
          lj[0] = l0; \
          t2 = 0; \
          phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
          ll_plain[t2+1] = ll_plain[t2] * (ths->n[t2] + HALO(ths->m)) + l_all[t2*(2*ths->m+2) + lj[t2]]; \
          for (l1 = 0; l1 < 2*ths->m+2; l1++) \
          { \
            lj[1] = l1; \
            t2 = 1; \
            phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
            ll_plain[t2+1] = ll_plain[t2] * (ths->n[t2] + HALO(ths->m)) + l_all[t2*(2*ths->m+2) + lj[t2]]; \
            for (l2 = 0; l2 < 2*ths->m+2; l2++) \
            { \
              lj[2] = l2; \
              t2 = 2; \
              phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
              ll_plain[t2+1] = ll_plain[t2] * (ths->n[t2] + HALO(ths->m)) + l_all[t2*(2*ths->m+2) + lj[t2]]; \
              for (l3 = 0; l3 < 2*ths->m+2; l3++) \
              { \
                lj[3] = l3; \
                t2 = 3; \
                phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                ll_plain[t2+1] = ll_plain[t2] * (ths->n[t2] + HALO(ths->m)) + l_all[t2*(2*ths->m+2) + lj[t2]]; \
                for (l4 = 0; l4 < 2*ths->m+2; l4++) \
                { \
                  lj[4] = l4; \
                  t2 = 4; \
                  phi_prod[t2+1] = phi_prod[t2] * MACRO_COMPUTE_ ## whichone; \
                  ll_plain[t2+1] = ll_plain[t2] * (ths->n[t2] + HALO(ths->m)) + l_all[t2*(2*ths->m+2) + lj[t2]]; \
 \
                  C *lhs = ths->g + ll_plain[ths->d]; \
                  R *lhs_real = (R*)lhs; \
                  C val = phi_prod[ths->d] * fs[j]; \
 \
                  _Pragma("omp atomic") \
                  lhs_real[0] += CREAL(val); \
 \
                  _Pragma("omp atomic") \
                  lhs_real[1] += CIMAG(val); \
                } \
              } \
            } \
          } \
        } \
      } /* if(d==5) */ \
      else { \
        for (l_L = 0; l_L < lprod; l_L++) \
        { \
          C *lhs; \
          R *lhs_real; \
          C val; \
 \
          MACRO_adjoint_nd_B_OMP_COMPUTE_UPDATE_ ## whichone \
 \
          lhs = ths->g + ll_plain[ths->d]; \
          lhs_real = (R*)lhs; \
          val = phi_prod[ths->d] * fs[j]; \
 \
          _Pragma("omp atomic") \
          lhs_real[0] += CREAL(val); \
 \
          _Pragma("omp atomic") \
          lhs_real[1] += CIMAG(val); \
 \
          MACRO_count_uo_l_lj_t; \
        } /* for(l_L) */ \
      } \
}

static inline void B_openmp_T(X(plan) *ths)
{
  const R *xs = node_x(ths);
  C *fs = node_f(ths);
  INT lprod; /* 'regular bandwidth' of matrix B  */
  INT k;

  zero_grid_B(ths);

  if (adjoint_B_tiled(ths))
    return;

  for (k = 0, lprod = 1; k < ths->d; k++)
    lprod *= (2*ths->m+2);

  if (ths->flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(ths->g, ths->psi_index_g, ths->psi, fs,
        ths->M_total, ths->d, ths->n, ths->m, ths->flags, ths->index_x,
        ths->x_sorted != NULL);
    return;
  }

  if (ths->flags & PRE_PSI)
  {

    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k < ths->M_total; k++)
    {
      INT t, t2; /* index dimensions */ \
      INT j = NODE_J(ths, k);
      MACRO_adjoint_nd_B_OMP_COMPUTE(with_PRE_PSI);
    } /* for(j) */
    return;
  } /* if(PRE_PSI) */

  if (ths->flags & PRE_FG_PSI)
  {
    INT t, t2; /* index dimensions */
    R fg_e[ths->d], fg_q[ths->d];
    for(t2 = 0; t2 < ths->d; t2++)
    {
      fg_e[t2] = EXP(K(-1.0)/ths->b[t2]);
      fg_q[t2] = fg_e[t2] * fg_e[t2];
    }


    #pragma omp parallel for default(shared) private(k,t,t2)
    for (k = 0; k < ths->M_total; k++)
    {
      INT j = NODE_J(ths, k);
      MACRO_adjoint_nd_B_OMP_COMPUTE(with_PRE_FG_PSI);
    } /* for(j) */
    return;
  } /* if(PRE_FG_PSI) */

  if (ths->flags & FG_PSI)
  {
    INT t, t2; /* index dimensions */
    R fg_e[ths->d], fg_q[ths->d];

    sort(ths);

    for (t2 = 0; t2 < ths->d; t2++)
    {
      fg_e[t2] = EXP(K(-1.0)/ths->b[t2]);
      fg_q[t2] = fg_e[t2] * fg_e[t2];
    }


    #pragma omp parallel for default(shared) private(k,t,t2)
    for (k = 0; k < ths->M_total; k++)
    {
      INT j = NODE_J(ths, k);
      MACRO_adjoint_nd_B_OMP_COMPUTE(with_FG_PSI);
    } /* for(j) */
    return;
  } /* if(FG_PSI) */

  if (ths->flags & PRE_LIN_PSI)
  {
    sort(ths);


    #pragma omp parallel for default(shared) private(k)
    for (k = 0; k<ths->M_total; k++)
    {
      INT t, t2; /* index dimensions */
      INT j = NODE_J(ths, k);
      MACRO_adjoint_nd_B_OMP_COMPUTE(with_PRE_LIN_PSI);
    } /* for(j) */
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  sort(ths);


  #pragma omp parallel for default(shared) private(k)
  for (k = 0; k < ths->M_total; k++)
  {
    INT t, t2; /* index dimensions */
    INT j = NODE_J(ths, k);
    MACRO_adjoint_nd_B_OMP_COMPUTE(without_PRE_PSI);
  } /* for(j) */
}
#endif

static void B_T(X(plan) *ths)
{
#ifdef _OPENMP
  B_openmp_T(ths);
#else
  B_serial_T(ths);
#endif
}

/* ## specialized version for d=1  ########################################### */

/* e = exp(-1/b) and q = e * e, the outward step factors FG_RUN needs. */
static void nfft_init_fg(R *e, R *q, const R b)
{
  *e = EXP(K(-1.0)/b);
  *q = (*e) * (*e);
}


static void nfft_trafo_1d_compute(C *fj, const C *g,const R *psij_const,
  const R *xj, const INT n, const INT m)
{
  const INT u = uo_wrap(*xj, n, m);
  const C *gj = g + u;
  C a0 = K(0.0), a1 = K(0.0);
  INT l;

  /* local accumulators (fj may alias g), two of them for the 2m+2 taps */
  for (l = 0; l < 2*m+2; l += 2)
  {
    a0 += psij_const[l] * gj[l];
    a1 += psij_const[l+1] * gj[l+1];
  }

  *fj = a0 + a1;
}

#ifndef _OPENMP
static void nfft_adjoint_1d_compute_serial(const C *restrict fj,
    C *restrict g, const R *restrict psij_const, const R *restrict xj,
    const INT n, const INT m)
{
  const INT u = uo_wrap(*xj, n, m);
  C *gj = g + u;
  const R *psij = psij_const;
  INT l;

  for (l = 0; l <= 2*m+1; l++)
    (*gj++) += (*psij++) * (*fj);
}
#endif

#ifdef _OPENMP
/* adjoint NFFT one-dimensional case with OpenMP atomic operations */
static void nfft_adjoint_1d_compute_omp_atomic(const C f, C *restrict g,
    const R *restrict psij_const, const R *restrict xj, const INT n,
    const INT m)
{
  const INT u = uo_wrap(*xj, n, m);
  INT l;

  for (l = 0; l <= 2*m+1; l++)
  {
    C *lhs = g+u+l;
    R *lhs_real = (R*)lhs;
    C val = psij_const[l] * f;
    #pragma omp atomic
    lhs_real[0] += CREAL(val);

    #pragma omp atomic
    lhs_real[1] += CIMAG(val);
  }
}
#endif

/* Window buffers come from the heap, once per thread: GCC 13 and newer abort
 * expanding this function's OpenMP regions on aarch64-apple-darwin when they
 * hold a run-time sized stack array. */
static void nfft_trafo_1d_B(X(plan) *ths)
{
  const R *xs = node_x(ths);
  C *fs = node_f(ths);
  const INT n = ths->n[0], M = ths->M_total, m = ths->m, m2p2 = 2*m+2;
  const C *g = (C*)ths->g;

  if (ths->flags & PRE_FULL_PSI)
  {
    INT k;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
      fs[j] = full_psi_trafo(ths->psi + j*m2p2, ths->psi_index_g + j*m2p2,
          g, m2p2);
    }
    return;
  } /* if(PRE_FULL_PSI) */

  if (ths->flags & PRE_PSI)
  {
    INT k;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
      nfft_trafo_1d_compute(&fs[j], g, ths->psi + j * (2 * m + 2),
        &xs[j], n, m);
    }
    return;
  } /* if(PRE_PSI) */

  if (ths->flags & PRE_FG_PSI)
  {
    R fg_e[1], fg_q[1];

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);

#ifdef _OPENMP
    #pragma omp parallel default(shared)
#endif
    {
    INT k;
    R *psij_const = (R*)Y(malloc)((size_t)(2 * m + 3) * sizeof(R));
#ifdef _OPENMP
    #pragma omp for
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
      const R *psij0;

      {
        FG_RUN(psij_const, m, FABS(ths->psi[2 * j]), ths->psi[2 * j + 1],
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, ths->psi[2 * j]);
      }

      nfft_trafo_1d_compute(&fs[j], g, psij0, &xs[j], n, m);
    }
    Y(free)(psij_const);
    }

    return;
  } /* if(PRE_FG_PSI) */

  if (ths->flags & FG_PSI)
  {
    R fg_e[1], fg_q[1];

    sort(ths);

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);

#ifdef _OPENMP
    #pragma omp parallel default(shared)
#endif
    {
    INT k;
    R *psij_const = (R*)Y(malloc)((size_t)(2 * m + 3) * sizeof(R));
#ifdef _OPENMP
    #pragma omp for
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
      const R *psij0;

      {
        const R fg_nx = (R)(n) * xs[j];
        const R fg_d =
            NX_SUB(n, xs[j], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const, m,
            (PHI(ths->n[0], fg_d / (R)(n), 0)),
            EXP(K(2.0) * fg_d / ths->b[0]),
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, fg_d);
      }

      nfft_trafo_1d_compute(&fs[j], g, psij0, &xs[j], n, m);
    }
    Y(free)(psij_const);
    }
    return;
  } /* if(FG_PSI) */

  if (ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K, ip_s = K / (m + 2);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel default(shared)
#endif
    {
    INT k;
    R *psij_const = (R*)Y(malloc)((size_t)(m2p2) * sizeof(R));
#ifdef _OPENMP
    #pragma omp for
#endif
    for (k = 0; k < M; k++)
    {
      INT u, o, l;
      R ip_y, ip_w;
      INT ip_u;
      INT j = NODE_J(ths, k);

      uo(ths, xs, (INT)j, &u, &o, (INT)0);

      ip_y = FABS(NX_SUB(n, xs[j], u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);

      for (l = 0; l < m2p2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)] * (K(1.0) - ip_w)
          + ths->psi[ABS(ip_u-l*ip_s+1)] * (ip_w);

      nfft_trafo_1d_compute(&fs[j], g, psij_const, &xs[j], n, m);
    }
    Y(free)(psij_const);
    }
    return;
  } /* if(PRE_LIN_PSI) */
  else
  {
    /* no precomputed psi at all */

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel default(shared)
#endif
    {
    INT k;
    R *psij_const = (R*)Y(malloc)((size_t)(m2p2) * sizeof(R));
#ifdef _OPENMP
    #pragma omp for
#endif
    for (k = 0; k < M; k++)
    {
      INT u, o;
      INT j = NODE_J(ths, k);

      uo(ths, xs, (INT)j, &u, &o, (INT)0);

      PHI_RUN(psij_const, ths->n[0], xs[j], u, 0);

      nfft_trafo_1d_compute(&fs[j], g, psij_const, &xs[j], n, m);
    }
    Y(free)(psij_const);
    }
  }
}



static void nfft_adjoint_1d_B(X(plan) *ths)
{
  const R *xs = node_x(ths);
  C *fs = node_f(ths);
  const INT n = ths->n[0], M = ths->M_total, m = ths->m;
  INT k;
  C *g = (C*)ths->g;

  zero_grid_B(ths);

#ifdef _OPENMP
  if (adjoint_B_tiled(ths))
    return;
#endif

  if (ths->flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(g, ths->psi_index_g, ths->psi, fs, M,
        (INT)1, ths->n, m, ths->flags, ths->index_x,
        ths->x_sorted != NULL);
    return;
  } /* if(PRE_FULL_PSI) */

  if (ths->flags & PRE_PSI)
  {
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
#ifdef _OPENMP
      nfft_adjoint_1d_compute_omp_atomic(fs[j], g, ths->psi + j * (2 * m + 2), xs + j, n, m);
#else
      nfft_adjoint_1d_compute_serial(fs + j, g, ths->psi + j * (2 * m + 2), xs + j, n, m);
#endif
    }

    return;
  } /* if(PRE_PSI) */

  if (ths->flags & PRE_FG_PSI)
  {
    R fg_e[1], fg_q[1];

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);


#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[2 * m + 3];
      const R *psij0;
      INT j = NODE_J(ths, k);

      {
        FG_RUN(psij_const, m, FABS(ths->psi[2 * j]), ths->psi[2 * j + 1],
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, ths->psi[2 * j]);
      }

#ifdef _OPENMP
      nfft_adjoint_1d_compute_omp_atomic(fs[j], g, psij0, xs + j, n, m);
#else
      nfft_adjoint_1d_compute_serial(fs + j, g, psij0, xs + j, n, m);
#endif
    }

    return;
  } /* if(PRE_FG_PSI) */

  if (ths->flags & FG_PSI)
  {
    R fg_e[1], fg_q[1];

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[2 * m + 3];
      const R *psij0;
      INT j = NODE_J(ths, k);

      {
        const R fg_nx = (R)(n) * xs[j];
        const R fg_d =
            NX_SUB(n, xs[j], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const, m,
            (PHI(ths->n[0], fg_d / (R)(n), 0)),
            EXP(K(2.0) * fg_d / ths->b[0]),
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, fg_d);
      }

#ifdef _OPENMP
      nfft_adjoint_1d_compute_omp_atomic(fs[j], g, psij0, xs + j, n, m);
#else
      nfft_adjoint_1d_compute_serial(fs + j, g, psij0, xs + j, n, m);
#endif
    }

    return;
  } /* if(FG_PSI) */

  if (ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K;
    const INT ip_s = K / (m + 2);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u,o,l;
      INT ip_u;
      R ip_y, ip_w;
      INT j = NODE_J(ths, k);
      R psij_const[2 * m + 2];

      uo(ths, xs, j, &u, &o, (INT)0);

      ip_y = FABS(NX_SUB(n, xs[j], u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for (l = 0; l < 2 * m + 2; l++)
        psij_const[l]
            = ths->psi[ABS(ip_u-l*ip_s)] * (K(1.0) - ip_w)
                + ths->psi[ABS(ip_u-l*ip_s+1)] * (ip_w);

#ifdef _OPENMP
      nfft_adjoint_1d_compute_omp_atomic(fs[j], g, psij_const, xs + j, n, m);
#else
      nfft_adjoint_1d_compute_serial(fs + j, g, psij_const, xs + j, n, m);
#endif
    }
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  sort(ths);

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    INT u,o;
    R psij_const[2 * m + 2];
    INT j = NODE_J(ths, k);

    uo(ths, xs, j, &u, &o, (INT)0);

    PHI_RUN(psij_const, ths->n[0], xs[j], u, 0);

#ifdef _OPENMP
    nfft_adjoint_1d_compute_omp_atomic(fs[j], g, psij_const, xs + j, n, m);
#else
    nfft_adjoint_1d_compute_serial(fs + j, g, psij_const, xs + j, n, m);
#endif
  }
}

void X(trafo_1d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->n[0] <= 2*ths->m+2))
  {
    X(trafo_direct)(ths);
    return;
  }
  
  const INT N = ths->N[0], N2 = N/2, n = ths->n[0];
  C *f_hat1 = (C*)ths->f_hat, *f_hat2 = (C*)&ths->f_hat[N2];

  ths->g_hat = ths->g1;
  ths->g = ths->g2;

  {
    C *g_hat1 = (C*)&ths->g_hat[n-N/2], *g_hat2 = (C*)ths->g_hat;
    R *c_phi_inv1, *c_phi_inv2;

    TIC(0)
#ifdef _OPENMP
    {
      INT k;
      #pragma omp parallel for default(shared) private(k)
      for (k = 0; k < grid_fftw_total(ths); k++)
        ths->g_hat[k] = 0.0;
    }
#else
    memset(ths->g_hat, 0, (size_t)(grid_fftw_total(ths)) * sizeof(C));
#endif
    if(ths->flags & PRE_PHI_HUT)
    {
      INT k;
      c_phi_inv1 = ths->c_phi_inv[0];
      c_phi_inv2 = &ths->c_phi_inv[0][N2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k)
#endif
      for (k = 0; k < N2; k++)
      {
        g_hat1[k] = f_hat1[k] * c_phi_inv1[k];
        g_hat2[k] = f_hat2[k] * c_phi_inv2[k];
      }
    }
    else
    {
      INT k;
#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k)
#endif
      for (k = 0; k < N2; k++)
      {
        g_hat1[k] = f_hat1[k] / (PHI_HUT(ths->n[0],k-N2,0));
        g_hat2[k] = f_hat2[k] / (PHI_HUT(ths->n[0],k,0));
      }
    }
    TOC(0)

    TIC_FFTW(1)
    FFTW(execute)(ths->my_fftw_plan1);
    TOC_FFTW(1);
    grid_B_trafo(ths);

    TIC(2);
    nfft_trafo_1d_B(ths);
    scatter_f(ths);
    TOC(2);
  }
}

void X(adjoint_1d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->n[0] <= 2*ths->m+2))
  {
    X(adjoint_direct)(ths);
    return;
  }
  
  INT n,N;
  C *g_hat1,*g_hat2,*f_hat1,*f_hat2;
  R *c_phi_inv1, *c_phi_inv2;

  N=ths->N[0];
  n=ths->n[0];

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  f_hat1=(C*)ths->f_hat;
  f_hat2=(C*)&ths->f_hat[N/2];
  g_hat1=(C*)&ths->g_hat[n-N/2];
  g_hat2=(C*)ths->g_hat;

  TIC(2)
  select_grid_B(ths);
  gather_f(ths);
  nfft_adjoint_1d_B(ths);
  grid_B_adjoint_done(ths);
  TOC(2)

  TIC_FFTW(1)
  FFTW(execute)(ths->my_fftw_plan2);
  TOC_FFTW(1);

  TIC(0)
  if(ths->flags & PRE_PHI_HUT)
  {
    INT k;
    c_phi_inv1=ths->c_phi_inv[0];
    c_phi_inv2=&ths->c_phi_inv[0][N/2];

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < N/2; k++)
    {
      f_hat1[k] = g_hat1[k] * c_phi_inv1[k];
      f_hat2[k] = g_hat2[k] * c_phi_inv2[k];
    }
  }
  else
  {
    INT k;

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < N/2; k++)
    {
      f_hat1[k] = g_hat1[k] / (PHI_HUT(ths->n[0],k-N/2,0));
      f_hat2[k] = g_hat2[k] / (PHI_HUT(ths->n[0],k,0));
    }
  }
  TOC(0)
}


/* ################################################ SPECIFIC VERSIONS FOR d=2 */


static void nfft_trafo_2d_compute(C *fj, const C *g, const R *psij_const0,
    const R *psij_const1, const R *xj0, const R *xj1, const INT n0,
    const INT n1, const INT m)
{
  const INT L = 2 * m + 2;
  const INT np[2] = {n0 + HALO(m), n1 + HALO(m)};
  const INT lo[2] = {uo_wrap(*xj0, n0, m), uo_wrap(*xj1, n1, m)};
  const R *w[1] = {psij_const0};
  R wl[2 * L], acc[2 * L];

  dup_weights(wl, psij_const1, L);
  zero_reals(acc, 2 * L);
  Y(B_tensor_trafo)(2, L, 2, np, lo, w, wl, (const R*)g, acc);
  *fj = sum_pairs(acc, L);
}

#ifdef _OPENMP
/* adjoint NFFT two-dimensional case with OpenMP atomic operations */
static void nfft_adjoint_2d_compute_omp_atomic(const C f, C *restrict g,
            const R *restrict psij_const0, const R *restrict psij_const1,
            const R *restrict xj0, const R *restrict xj1, const INT n0,
            const INT n1, const INT m)
{
  const INT p1 = n1 + HALO(m), L = 2 * m + 2;
  const INT u0 = uo_wrap(*xj0, n0, m), u1 = uo_wrap(*xj1, n1, m);
  INT l0, l1;
  R wf[2 * L];

  scale_weights(wf, psij_const1, f, L);

  for(l0=0; l0<L; l0++)
  {
    const R w0 = psij_const0[l0];
    R *row = (R*)(g + (u0+l0) * p1 + u1);

    for(l1=0; l1<L; l1++)
    {
      #pragma omp atomic
      row[2*l1] += w0 * wf[2*l1];

      #pragma omp atomic
      row[2*l1+1] += w0 * wf[2*l1+1];
    }
  }
}
#endif

#ifndef _OPENMP
static void nfft_adjoint_2d_compute_serial(const C *restrict fj,
            C *restrict g, const R *restrict psij_const0,
            const R *restrict psij_const1, const R *restrict xj0,
            const R *restrict xj1, const INT n0, const INT n1, const INT m)
{
  const INT L = 2 * m + 2;
  const INT np[2] = {n0 + HALO(m), n1 + HALO(m)};
  const INT lo[2] = {uo_wrap(*xj0, n0, m), uo_wrap(*xj1, n1, m)};
  const R *w[1] = {psij_const0};
  R wf[2 * L];

  scale_weights(wf, psij_const1, *fj, L);
  Y(B_tensor_adjoint)(2, L, 2, np, lo, w, wf, (R*)g);
}
#endif

static void nfft_trafo_2d_B(X(plan) *ths)
{
  const R *xs = node_x(ths);
  C *fs = node_f(ths);
  const C *g = (C*)ths->g;
  const INT n0 = ths->n[0];
  const INT n1 = ths->n[1];
  const INT M = ths->M_total;
  const INT m = ths->m;

  INT k;

  if(ths->flags & PRE_FULL_PSI)
  {
    const INT lprod = (2*m+2) * (2*m+2);
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
      fs[j] = full_psi_trafo(ths->psi + j*lprod, ths->psi_index_g + j*lprod,
          g, lprod);
    }
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->flags & PRE_PSI)
  {
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
      nfft_trafo_2d_compute(fs + j, g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), xs+2*j, xs+2*j+1, n0, n1, m);
    }

      return;
  } /* if(PRE_PSI) */

  if(ths->flags & PRE_FG_PSI)
  {
    R fg_e[2], fg_q[2];

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);
    nfft_init_fg(&fg_e[1], &fg_q[1], ths->b[1]);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[2*(2*m+3)];
      const R *psij0, *psij1;
      INT j = NODE_J(ths, k);

      {
        FG_RUN(psij_const, m, FABS(ths->psi[2*j*2]), ths->psi[2*j*2+1],
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, ths->psi[2*j*2]);
      }

      {
        FG_RUN(psij_const + 2*m+3, m, FABS(ths->psi[2*(j*2+1)]), ths->psi[2*(j*2+1)+1],
            fg_e[1], fg_q[1]);
        psij1 = FG_RUN_START(psij_const + 2*m+3, ths->psi[2*(j*2+1)]);
      }

      nfft_trafo_2d_compute(fs + j, g, psij0, psij1, xs+2*j, xs+2*j+1, n0, n1, m);
    }

    return;
  } /* if(PRE_FG_PSI) */

  if(ths->flags & FG_PSI)
  {
    R fg_e[2], fg_q[2];

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);
    nfft_init_fg(&fg_e[1], &fg_q[1], ths->b[1]);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[2*(2*m+3)];
      const R *psij0, *psij1;
      INT j = NODE_J(ths, k);

      {
        const R fg_nx = (R)(n0) * xs[2*j];
        const R fg_d =
            NX_SUB(n0, xs[2*j], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const, m,
            (PHI(ths->n[0], fg_d / (R)(n0), 0)),
            EXP(K(2.0) * fg_d / ths->b[0]),
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, fg_d);
      }

      {
        const R fg_nx = (R)(n1) * xs[2*j+1];
        const R fg_d =
            NX_SUB(n1, xs[2*j+1], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const + 2*m+3, m,
            (PHI(ths->n[1], fg_d / (R)(n1), 1)),
            EXP(K(2.0) * fg_d / ths->b[1]),
            fg_e[1], fg_q[1]);
        psij1 = FG_RUN_START(psij_const + 2*m+3, fg_d);
      }

      nfft_trafo_2d_compute(fs + j, g, psij0, psij1, xs+2*j, xs+2*j+1, n0, n1, m);
    }

    return;
  } /* if(FG_PSI) */

  if(ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K, ip_s = K / (m + 2);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u, o, l;
      R ip_y, ip_w;
      INT ip_u;
      R psij_const[2*(2*m+2)];
      INT j = NODE_J(ths, k);

      uo(ths, xs, j,&u,&o,(INT)0);
      ip_y = FABS(NX_SUB(n0, xs[2*j], u)) * ((R)ip_s);
      ip_u = (INT)LRINT(FLOOR(ip_y));
      ip_w = ip_y - (R)(ip_u);
      for (l = 0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) + ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths, xs, j,&u,&o,(INT)1);
      ip_y = FABS(NX_SUB(n1, xs[2*j+1], u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for (l = 0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) + ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_trafo_2d_compute(fs + j, g, psij_const, psij_const+2*m+2, xs+2*j, xs+2*j+1, n0, n1, m);
    }
      return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */

  sort(ths);

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    R psij_const[2*(2*m+2)];
    INT u, o;
    INT j = NODE_J(ths, k);

    uo(ths, xs, j,&u,&o,(INT)0);
    PHI_RUN(psij_const, ths->n[0], xs[2*j], u, 0);

    uo(ths, xs, j,&u,&o,(INT)1);
    PHI_RUN(psij_const+2*m+2, ths->n[1], xs[2*j+1], u, 1);

    nfft_trafo_2d_compute(fs + j, g, psij_const, psij_const+2*m+2, xs+2*j, xs+2*j+1, n0, n1, m);
  }
}



static void nfft_adjoint_2d_B(X(plan) *ths)
{
  const R *xs = node_x(ths);
  C *fs = node_f(ths);
  const INT n0 = ths->n[0];
  const INT n1 = ths->n[1];
  const INT M = ths->M_total;
  const INT m = ths->m;
  C* g = (C*) ths->g;
  INT k;

  zero_grid_B(ths);

#ifdef _OPENMP
  if (adjoint_B_tiled(ths))
    return;
#endif

  if(ths->flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(g, ths->psi_index_g, ths->psi, fs, M,
        (INT)2, ths->n, m, ths->flags, ths->index_x,
        ths->x_sorted != NULL);
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->flags & PRE_PSI)
  {
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
#ifdef _OPENMP
      nfft_adjoint_2d_compute_omp_atomic(fs[j], g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), xs+2*j, xs+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute_serial(fs + j, g, ths->psi+j*2*(2*m+2), ths->psi+(j*2+1)*(2*m+2), xs+2*j, xs+2*j+1, n0, n1, m);
#endif
    }
    return;
  } /* if(PRE_PSI) */

  if(ths->flags & PRE_FG_PSI)
  {
    R fg_e[2], fg_q[2];

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);
    nfft_init_fg(&fg_e[1], &fg_q[1], ths->b[1]);


#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[2*(2*m+3)];
      const R *psij0, *psij1;
      INT j = NODE_J(ths, k);

      {
        FG_RUN(psij_const, m, FABS(ths->psi[2*j*2]), ths->psi[2*j*2+1],
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, ths->psi[2*j*2]);
      }

      {
        FG_RUN(psij_const + 2*m+3, m, FABS(ths->psi[2*(j*2+1)]), ths->psi[2*(j*2+1)+1],
            fg_e[1], fg_q[1]);
        psij1 = FG_RUN_START(psij_const + 2*m+3, ths->psi[2*(j*2+1)]);
      }

#ifdef _OPENMP
      nfft_adjoint_2d_compute_omp_atomic(fs[j], g, psij0, psij1, xs+2*j, xs+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute_serial(fs + j, g, psij0, psij1, xs+2*j, xs+2*j+1, n0, n1, m);
#endif
    }

    return;
  } /* if(PRE_FG_PSI) */

  if(ths->flags & FG_PSI)
  {
    R fg_e[2], fg_q[2];

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);
    nfft_init_fg(&fg_e[1], &fg_q[1], ths->b[1]);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[2*(2*m+3)];
      const R *psij0, *psij1;
      INT j = NODE_J(ths, k);

      {
        const R fg_nx = (R)(n0) * xs[2*j];
        const R fg_d =
            NX_SUB(n0, xs[2*j], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const, m,
            (PHI(ths->n[0], fg_d / (R)(n0), 0)),
            EXP(K(2.0) * fg_d / ths->b[0]),
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, fg_d);
      }

      {
        const R fg_nx = (R)(n1) * xs[2*j+1];
        const R fg_d =
            NX_SUB(n1, xs[2*j+1], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const + 2*m+3, m,
            (PHI(ths->n[1], fg_d / (R)(n1), 1)),
            EXP(K(2.0) * fg_d / ths->b[1]),
            fg_e[1], fg_q[1]);
        psij1 = FG_RUN_START(psij_const + 2*m+3, fg_d);
      }

#ifdef _OPENMP
      nfft_adjoint_2d_compute_omp_atomic(fs[j], g, psij0, psij1, xs+2*j, xs+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute_serial(fs + j, g, psij0, psij1, xs+2*j, xs+2*j+1, n0, n1, m);
#endif
    }

    return;
  } /* if(FG_PSI) */

  if(ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K;
    const INT ip_s = K / (m + 2);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u,o,l;
      INT ip_u;
      R ip_y, ip_w;
      INT j = NODE_J(ths, k);
      R psij_const[2*(2*m+2)];

      uo(ths, xs, j,&u,&o,(INT)0);
      ip_y = FABS(NX_SUB(n0, (xs[2*j]), u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths, xs, j,&u,&o,(INT)1);
      ip_y = FABS(NX_SUB(n1, (xs[2*j+1]), u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

#ifdef _OPENMP
      nfft_adjoint_2d_compute_omp_atomic(fs[j], g, psij_const, psij_const+2*m+2, xs+2*j, xs+2*j+1, n0, n1, m);
#else
      nfft_adjoint_2d_compute_serial(fs + j, g, psij_const, psij_const+2*m+2, xs+2*j, xs+2*j+1, n0, n1, m);
#endif
  }
      return;
    } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  sort(ths);

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    INT u,o;
    R psij_const[2*(2*m+2)];
    INT j = NODE_J(ths, k);

    uo(ths, xs, j,&u,&o,(INT)0);
    PHI_RUN(psij_const, ths->n[0], xs[2*j], u, 0);

    uo(ths, xs, j,&u,&o,(INT)1);
    PHI_RUN(psij_const+2*m+2, ths->n[1], xs[2*j+1], u, 1);

#ifdef _OPENMP
    nfft_adjoint_2d_compute_omp_atomic(fs[j], g, psij_const, psij_const+2*m+2, xs+2*j, xs+2*j+1, n0, n1, m);
#else
    nfft_adjoint_2d_compute_serial(fs + j, g, psij_const, psij_const+2*m+2, xs+2*j, xs+2*j+1, n0, n1, m);
#endif
  }
}


void X(trafo_2d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->N[1] <= ths->m) || (ths->n[0] <= 2*ths->m+2) || (ths->n[1] <= 2*ths->m+2))
  {
    X(trafo_direct)(ths);
    return;
  }
  
  INT k0,k1,n0,n1,N0,N1,p1;
  C *g_hat,*f_hat;
  R *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12;
  R ck01, ck02, ck11, ck12;
  C *g_hat11,*f_hat11,*g_hat21,*f_hat21,*g_hat12,*f_hat12,*g_hat22,*f_hat22;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  n0=ths->n[0];
  n1=ths->n[1];
  p1=n1+ths->halo_hi; /* row length of g_hat */

  f_hat=(C*)ths->f_hat;
  g_hat=(C*)ths->g_hat;

  TIC(0)
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k0)
  for (k0 = 0; k0 < grid_fftw_total(ths); k0++)
    ths->g_hat[k0] = 0.0;
#else
  memset(ths->g_hat, 0, (size_t)(grid_fftw_total(ths)) * sizeof(C));
#endif
  if(ths->flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,c_phi_inv11,c_phi_inv12,g_hat11,f_hat11,g_hat21,f_hat21,g_hat12,f_hat12,g_hat22,f_hat22,ck11,ck12)
#endif
      for(k0=0;k0<N0/2;k0++)
      {
        ck01=c_phi_inv01[k0];
        ck02=c_phi_inv02[k0];

        c_phi_inv11=ths->c_phi_inv[1];
        c_phi_inv12=&ths->c_phi_inv[1][N1/2];

        g_hat11=g_hat + (n0-(N0/2)+k0)*p1+n1-(N1/2);
        f_hat11=f_hat + k0*N1;
        g_hat21=g_hat + k0*p1+n1-(N1/2);
        f_hat21=f_hat + ((N0/2)+k0)*N1;
        g_hat12=g_hat + (n0-(N0/2)+k0)*p1;
        f_hat12=f_hat + k0*N1+(N1/2);
        g_hat22=g_hat + k0*p1;
        f_hat22=f_hat + ((N0/2)+k0)*N1+(N1/2);

        for(k1=0;k1<N1/2;k1++)
        {
          ck11=c_phi_inv11[k1];
          ck12=c_phi_inv12[k1];

          g_hat11[k1] = f_hat11[k1] * ck01 * ck11;
          g_hat21[k1] = f_hat21[k1] * ck02 * ck11;
          g_hat12[k1] = f_hat12[k1] * ck01 * ck12;
          g_hat22[k1] = f_hat22[k1] * ck02 * ck12;
        }
      }
    }
  else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,ck11,ck12)
#endif
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(ths->n[0],k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(ths->n[0],k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(ths->n[1],k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(ths->n[1],k1,1));
      g_hat[(n0-N0/2+k0)*p1+n1-N1/2+k1] = f_hat[k0*N1+k1]             * ck01 * ck11;
      g_hat[k0*p1+n1-N1/2+k1]           = f_hat[(N0/2+k0)*N1+k1]      * ck02 * ck11;
      g_hat[(n0-N0/2+k0)*p1+k1]         = f_hat[k0*N1+N1/2+k1]        * ck01 * ck12;
      g_hat[k0*p1+k1]                   = f_hat[(N0/2+k0)*N1+N1/2+k1] * ck02 * ck12;
    }
      }

  TOC(0)

  TIC_FFTW(1)
  FFTW(execute)(ths->my_fftw_plan1);
  TOC_FFTW(1);
  grid_B_trafo(ths);

  TIC(2);
  nfft_trafo_2d_B(ths);
  scatter_f(ths);
  TOC(2);
}

void X(adjoint_2d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->N[1] <= ths->m) || (ths->n[0] <= 2*ths->m+2) || (ths->n[1] <= 2*ths->m+2))
  {
    X(adjoint_direct)(ths);
    return;
  }
  
  INT k0,k1,n0,n1,N0,N1,p1;
  C *g_hat,*f_hat;
  R *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12;
  R ck01, ck02, ck11, ck12;
  C *g_hat11,*f_hat11,*g_hat21,*f_hat21,*g_hat12,*f_hat12,*g_hat22,*f_hat22;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  n0=ths->n[0];
  n1=ths->n[1];
  p1=n1+ths->halo_hi; /* row length of g_hat */

  f_hat=(C*)ths->f_hat;
  g_hat=(C*)ths->g_hat;

  TIC(2);
  select_grid_B(ths);
  gather_f(ths);
  nfft_adjoint_2d_B(ths);
  grid_B_adjoint_done(ths);
  TOC(2);

  TIC_FFTW(1)
  FFTW(execute)(ths->my_fftw_plan2);
  TOC_FFTW(1);

  TIC(0)
  if(ths->flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,c_phi_inv11,c_phi_inv12,g_hat11,f_hat11,g_hat21,f_hat21,g_hat12,f_hat12,g_hat22,f_hat22,ck11,ck12)
#endif
      for(k0=0;k0<N0/2;k0++)
      {
        ck01=c_phi_inv01[k0];
        ck02=c_phi_inv02[k0];

        c_phi_inv11=ths->c_phi_inv[1];
        c_phi_inv12=&ths->c_phi_inv[1][N1/2];

        g_hat11=g_hat + (n0-(N0/2)+k0)*p1+n1-(N1/2);
        f_hat11=f_hat + k0*N1;
        g_hat21=g_hat + k0*p1+n1-(N1/2);
        f_hat21=f_hat + ((N0/2)+k0)*N1;
        g_hat12=g_hat + (n0-(N0/2)+k0)*p1;
        f_hat12=f_hat + k0*N1+(N1/2);
        g_hat22=g_hat + k0*p1;
        f_hat22=f_hat + ((N0/2)+k0)*N1+(N1/2);

        for(k1=0;k1<N1/2;k1++)
        {
          ck11=c_phi_inv11[k1];
          ck12=c_phi_inv12[k1];

          f_hat11[k1] = g_hat11[k1] * ck01 * ck11;
          f_hat21[k1] = g_hat21[k1] * ck02 * ck11;
          f_hat12[k1] = g_hat12[k1] * ck01 * ck12;
          f_hat22[k1] = g_hat22[k1] * ck02 * ck12;
        }
      }
    }
  else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k0,k1,ck01,ck02,ck11,ck12)
#endif
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(ths->n[0],k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(ths->n[0],k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(ths->n[1],k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(ths->n[1],k1,1));
      f_hat[k0*N1+k1]             = g_hat[(n0-N0/2+k0)*p1+n1-N1/2+k1] * ck01 * ck11;
      f_hat[(N0/2+k0)*N1+k1]      = g_hat[k0*p1+n1-N1/2+k1]           * ck02 * ck11;
      f_hat[k0*N1+N1/2+k1]        = g_hat[(n0-N0/2+k0)*p1+k1]         * ck01 * ck12;
      f_hat[(N0/2+k0)*N1+N1/2+k1] = g_hat[k0*p1+k1]                   * ck02 * ck12;
    }
      }
  TOC(0)
}

/* ################################################ SPECIFIC VERSIONS FOR d=3 */


static void nfft_trafo_3d_compute(C *fj, const C *g, const R *psij_const0,
    const R *psij_const1, const R *psij_const2, const R *xj0, const R *xj1,
    const R *xj2, const INT n0, const INT n1, const INT n2, const INT m)
{
  const INT L = 2 * m + 2;
  const INT np[3] = {n0 + HALO(m), n1 + HALO(m), n2 + HALO(m)};
  const INT lo[3] = {uo_wrap(*xj0, n0, m), uo_wrap(*xj1, n1, m),
    uo_wrap(*xj2, n2, m)};
  const R *w[2] = {psij_const0, psij_const1};
  R wl[2 * L], acc[2 * L];

  dup_weights(wl, psij_const2, L);
  zero_reals(acc, 2 * L);
  Y(B_tensor_trafo)(3, L, 2, np, lo, w, wl, (const R*)g, acc);
  *fj = sum_pairs(acc, L);
}

#ifdef _OPENMP
/* adjoint NFFT three-dimensional case with OpenMP atomic operations */
static void nfft_adjoint_3d_compute_omp_atomic(const C f, C *restrict g,
    const R *restrict psij_const0, const R *restrict psij_const1,
    const R *restrict psij_const2, const R *restrict xj0,
    const R *restrict xj1, const R *restrict xj2, const INT n0, const INT n1,
    const INT n2, const INT m)
{
  const INT p1 = n1 + HALO(m), p2 = n2 + HALO(m), L = 2 * m + 2;
  const INT u0 = uo_wrap(*xj0, n0, m), u1 = uo_wrap(*xj1, n1, m),
    u2 = uo_wrap(*xj2, n2, m);
  INT l0, l1, l2;
  R wf[2 * L];

  scale_weights(wf, psij_const2, f, L);

  for(l0=0; l0<L; l0++)
  {
    for(l1=0; l1<L; l1++)
    {
      const R w01 = psij_const0[l0] * psij_const1[l1];
      R *row = (R*)(g + ((u0+l0) * p1 + u1+l1) * p2 + u2);

      for(l2=0; l2<L; l2++)
      {
#pragma omp atomic
        row[2*l2] += w01 * wf[2*l2];

#pragma omp atomic
        row[2*l2+1] += w01 * wf[2*l2+1];
      }
    }
  }
}
#endif

#ifndef _OPENMP
static void nfft_adjoint_3d_compute_serial(const C *restrict fj,
    C *restrict g, const R *restrict psij_const0,
    const R *restrict psij_const1, const R *restrict psij_const2,
    const R *restrict xj0, const R *restrict xj1, const R *restrict xj2,
    const INT n0, const INT n1, const INT n2, const INT m)
{
  const INT L = 2 * m + 2;
  const INT np[3] = {n0 + HALO(m), n1 + HALO(m), n2 + HALO(m)};
  const INT lo[3] = {uo_wrap(*xj0, n0, m), uo_wrap(*xj1, n1, m),
    uo_wrap(*xj2, n2, m)};
  const R *w[2] = {psij_const0, psij_const1};
  R wf[2 * L];

  scale_weights(wf, psij_const2, *fj, L);
  Y(B_tensor_adjoint)(3, L, 2, np, lo, w, wf, (R*)g);
}
#endif

static void nfft_trafo_3d_B(X(plan) *ths)
{
  const R *xs = node_x(ths);
  C *fs = node_f(ths);
  const INT n0 = ths->n[0];
  const INT n1 = ths->n[1];
  const INT n2 = ths->n[2];
  const INT M = ths->M_total;
  const INT m = ths->m;

  const C* g = (C*) ths->g;

  INT k;

  if(ths->flags & PRE_FULL_PSI)
  {
    const INT lprod = (2*m+2) * (2*m+2) * (2*m+2);
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
      fs[j] = full_psi_trafo(ths->psi + j*lprod, ths->psi_index_g + j*lprod,
          g, lprod);
    }
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->flags & PRE_PSI)
  {
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
      nfft_trafo_3d_compute(fs + j, g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
    }
    return;
  } /* if(PRE_PSI) */

  if(ths->flags & PRE_FG_PSI)
  {
    R fg_e[3], fg_q[3];

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);
    nfft_init_fg(&fg_e[1], &fg_q[1], ths->b[1]);
    nfft_init_fg(&fg_e[2], &fg_q[2], ths->b[2]);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
      R psij_const[3*(2*m+3)];
      const R *psij0, *psij1, *psij2;

      {
        FG_RUN(psij_const, m, FABS(ths->psi[2*j*3]), ths->psi[2*j*3+1],
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, ths->psi[2*j*3]);
      }

      {
        FG_RUN(psij_const + 2*m+3, m, FABS(ths->psi[2*(j*3+1)]), ths->psi[2*(j*3+1)+1],
            fg_e[1], fg_q[1]);
        psij1 = FG_RUN_START(psij_const + 2*m+3, ths->psi[2*(j*3+1)]);
      }

      {
        FG_RUN(psij_const + 2*(2*m+3), m, FABS(ths->psi[2*(j*3+2)]), ths->psi[2*(j*3+2)+1],
            fg_e[2], fg_q[2]);
        psij2 = FG_RUN_START(psij_const + 2*(2*m+3), ths->psi[2*(j*3+2)]);
      }

      nfft_trafo_3d_compute(fs + j, g, psij0, psij1, psij2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
    }

    return;
  } /* if(PRE_FG_PSI) */

  if(ths->flags & FG_PSI)
  {
    R fg_e[3], fg_q[3];

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);
    nfft_init_fg(&fg_e[1], &fg_q[1], ths->b[1]);
    nfft_init_fg(&fg_e[2], &fg_q[2], ths->b[2]);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
      R psij_const[3*(2*m+3)];
      const R *psij0, *psij1, *psij2;

      {
        const R fg_nx = (R)(n0) * xs[3*j];
        const R fg_d =
            NX_SUB(n0, xs[3*j], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const, m,
            (PHI(ths->n[0], fg_d / (R)(n0), 0)),
            EXP(K(2.0) * fg_d / ths->b[0]),
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, fg_d);
      }

      {
        const R fg_nx = (R)(n1) * xs[3*j+1];
        const R fg_d =
            NX_SUB(n1, xs[3*j+1], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const + 2*m+3, m,
            (PHI(ths->n[1], fg_d / (R)(n1), 1)),
            EXP(K(2.0) * fg_d / ths->b[1]),
            fg_e[1], fg_q[1]);
        psij1 = FG_RUN_START(psij_const + 2*m+3, fg_d);
      }

      {
        const R fg_nx = (R)(n2) * xs[3*j+2];
        const R fg_d =
            NX_SUB(n2, xs[3*j+2], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const + 2*(2*m+3), m,
            (PHI(ths->n[2], fg_d / (R)(n2), 2)),
            EXP(K(2.0) * fg_d / ths->b[2]),
            fg_e[2], fg_q[2]);
        psij2 = FG_RUN_START(psij_const + 2*(2*m+3), fg_d);
      }

      nfft_trafo_3d_compute(fs + j, g, psij0, psij1, psij2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
    }

    return;
  } /* if(FG_PSI) */

  if(ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K, ip_s = K / (m + 2);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u, o, l;
      R ip_y, ip_w;
      INT ip_u;
      R psij_const[3*(2*m+2)];
      INT j = NODE_J(ths, k);

      uo(ths, xs, j,&u,&o,(INT)0);
      ip_y = FABS(NX_SUB(n0, xs[3*j+0], u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths, xs, j,&u,&o,(INT)1);
      ip_y = FABS(NX_SUB(n1, xs[3*j+1], u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths, xs, j,&u,&o,(INT)2);
      ip_y = FABS(NX_SUB(n2, xs[3*j+2], u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[2*(2*m+2)+l] = ths->psi[2*(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[2*(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      nfft_trafo_3d_compute(fs + j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
    }
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */

  sort(ths);

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    R psij_const[3*(2*m+2)];
    INT u, o;
    INT j = NODE_J(ths, k);

    uo(ths, xs, j,&u,&o,(INT)0);
    PHI_RUN(psij_const, ths->n[0], xs[3*j], u, 0);

    uo(ths, xs, j,&u,&o,(INT)1);
    PHI_RUN(psij_const+2*m+2, ths->n[1], xs[3*j+1], u, 1);

    uo(ths, xs, j,&u,&o,(INT)2);
    PHI_RUN(psij_const+2*(2*m+2), ths->n[2], xs[3*j+2], u, 2);

    nfft_trafo_3d_compute(fs + j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
  }
}


static void nfft_adjoint_3d_B(X(plan) *ths)
{
  const R *xs = node_x(ths);
  C *fs = node_f(ths);
  INT k;
  const INT n0 = ths->n[0];
  const INT n1 = ths->n[1];
  const INT n2 = ths->n[2];
  const INT M = ths->M_total;
  const INT m = ths->m;

  C* g = (C*) ths->g;

  zero_grid_B(ths);

#ifdef _OPENMP
  if (adjoint_B_tiled(ths))
    return;
#endif

  if(ths->flags & PRE_FULL_PSI)
  {
    nfft_adjoint_B_compute_full_psi(g, ths->psi_index_g, ths->psi, fs, M,
        (INT)3, ths->n, m, ths->flags, ths->index_x,
        ths->x_sorted != NULL);
    return;
  } /* if(PRE_FULL_PSI) */

  if(ths->flags & PRE_PSI)
  {
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
#ifdef _OPENMP
      nfft_adjoint_3d_compute_omp_atomic(fs[j], g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute_serial(fs + j, g, ths->psi+j*3*(2*m+2), ths->psi+(j*3+1)*(2*m+2), ths->psi+(j*3+2)*(2*m+2), xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
#endif
    }
    return;
  } /* if(PRE_PSI) */

  if(ths->flags & PRE_FG_PSI)
  {
    R fg_e[3], fg_q[3];

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);
    nfft_init_fg(&fg_e[1], &fg_q[1], ths->b[1]);
    nfft_init_fg(&fg_e[2], &fg_q[2], ths->b[2]);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      R psij_const[3*(2*m+3)];
      const R *psij0, *psij1, *psij2;
      INT j = NODE_J(ths, k);

      {
        FG_RUN(psij_const, m, FABS(ths->psi[2*j*3]), ths->psi[2*j*3+1],
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, ths->psi[2*j*3]);
      }

      {
        FG_RUN(psij_const + 2*m+3, m, FABS(ths->psi[2*(j*3+1)]), ths->psi[2*(j*3+1)+1],
            fg_e[1], fg_q[1]);
        psij1 = FG_RUN_START(psij_const + 2*m+3, ths->psi[2*(j*3+1)]);
      }

      {
        FG_RUN(psij_const + 2*(2*m+3), m, FABS(ths->psi[2*(j*3+2)]), ths->psi[2*(j*3+2)+1],
            fg_e[2], fg_q[2]);
        psij2 = FG_RUN_START(psij_const + 2*(2*m+3), ths->psi[2*(j*3+2)]);
      }

#ifdef _OPENMP
      nfft_adjoint_3d_compute_omp_atomic(fs[j], g, psij0, psij1, psij2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute_serial(fs + j, g, psij0, psij1, psij2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
#endif
    }

    return;
  } /* if(PRE_FG_PSI) */

  if(ths->flags & FG_PSI)
  {
    R fg_e[3], fg_q[3];

    nfft_init_fg(&fg_e[0], &fg_q[0], ths->b[0]);
    nfft_init_fg(&fg_e[1], &fg_q[1], ths->b[1]);
    nfft_init_fg(&fg_e[2], &fg_q[2], ths->b[2]);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT j = NODE_J(ths, k);
      R psij_const[3*(2*m+3)];
      const R *psij0, *psij1, *psij2;

      {
        const R fg_nx = (R)(n0) * xs[3*j];
        const R fg_d =
            NX_SUB(n0, xs[3*j], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const, m,
            (PHI(ths->n[0], fg_d / (R)(n0), 0)),
            EXP(K(2.0) * fg_d / ths->b[0]),
            fg_e[0], fg_q[0]);
        psij0 = FG_RUN_START(psij_const, fg_d);
      }

      {
        const R fg_nx = (R)(n1) * xs[3*j+1];
        const R fg_d =
            NX_SUB(n1, xs[3*j+1], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const + 2*m+3, m,
            (PHI(ths->n[1], fg_d / (R)(n1), 1)),
            EXP(K(2.0) * fg_d / ths->b[1]),
            fg_e[1], fg_q[1]);
        psij1 = FG_RUN_START(psij_const + 2*m+3, fg_d);
      }

      {
        const R fg_nx = (R)(n2) * xs[3*j+2];
        const R fg_d =
            NX_SUB(n2, xs[3*j+2], FLOOR(fg_nx + K(0.5)));
        FG_RUN(psij_const + 2*(2*m+3), m,
            (PHI(ths->n[2], fg_d / (R)(n2), 2)),
            EXP(K(2.0) * fg_d / ths->b[2]),
            fg_e[2], fg_q[2]);
        psij2 = FG_RUN_START(psij_const + 2*(2*m+3), fg_d);
      }

#ifdef _OPENMP
      nfft_adjoint_3d_compute_omp_atomic(fs[j], g, psij0, psij1, psij2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute_serial(fs + j, g, psij0, psij1, psij2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
#endif
    }

    return;
  } /* if(FG_PSI) */

  if(ths->flags & PRE_LIN_PSI)
  {
    const INT K = ths->K;
    const INT ip_s = K / (m + 2);

    sort(ths);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k)
#endif
    for (k = 0; k < M; k++)
    {
      INT u,o,l;
      INT ip_u;
      R ip_y, ip_w;
      INT j = NODE_J(ths, k);
      R psij_const[3*(2*m+2)];

      uo(ths, xs, j,&u,&o,(INT)0);
      ip_y = FABS(NX_SUB(n0, xs[3*j+0], u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[l] = ths->psi[ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths, xs, j,&u,&o,(INT)1);
      ip_y = FABS(NX_SUB(n1, xs[3*j+1], u)) * ((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[2*m+2+l] = ths->psi[(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

      uo(ths, xs, j,&u,&o,(INT)2);
      ip_y = FABS(NX_SUB(n2, xs[3*j+2], u))*((R)ip_s);
      ip_u = (INT)(LRINT(FLOOR(ip_y)));
      ip_w = ip_y - (R)(ip_u);
      for(l=0; l < 2*m+2; l++)
        psij_const[2*(2*m+2)+l] = ths->psi[2*(K+1)+ABS(ip_u-l*ip_s)]*(K(1.0)-ip_w) +
          ths->psi[2*(K+1)+ABS(ip_u-l*ip_s+1)]*(ip_w);

#ifdef _OPENMP
      nfft_adjoint_3d_compute_omp_atomic(fs[j], g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
#else
      nfft_adjoint_3d_compute_serial(fs + j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
#endif
    }
    return;
  } /* if(PRE_LIN_PSI) */

  /* no precomputed psi at all */
  sort(ths);

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k)
#endif
  for (k = 0; k < M; k++)
  {
    INT u,o;
    R psij_const[3*(2*m+2)];
    INT j = NODE_J(ths, k);

    uo(ths, xs, j,&u,&o,(INT)0);
    PHI_RUN(psij_const, ths->n[0], xs[3*j], u, 0);

    uo(ths, xs, j,&u,&o,(INT)1);
    PHI_RUN(psij_const+2*m+2, ths->n[1], xs[3*j+1], u, 1);

    uo(ths, xs, j,&u,&o,(INT)2);
    PHI_RUN(psij_const+2*(2*m+2), ths->n[2], xs[3*j+2], u, 2);

#ifdef _OPENMP
    nfft_adjoint_3d_compute_omp_atomic(fs[j], g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
#else
    nfft_adjoint_3d_compute_serial(fs + j, g, psij_const, psij_const+2*m+2, psij_const+(2*m+2)*2, xs+3*j, xs+3*j+1, xs+3*j+2, n0, n1, n2, m);
#endif
  }
}


void X(trafo_3d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->N[1] <= ths->m) || (ths->N[2] <= ths->m) || (ths->n[0] <= 2*ths->m+2) || (ths->n[1] <= 2*ths->m+2) || (ths->n[2] <= 2*ths->m+2))
  {
    X(trafo_direct)(ths);
    return;
  }
  
  INT k0,k1,k2,n0,n1,n2,N0,N1,N2,p1,p2;
  C *g_hat,*f_hat;
  R *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12, *c_phi_inv21, *c_phi_inv22;
  R ck01, ck02, ck11, ck12, ck21, ck22;
  C *g_hat111,*f_hat111,*g_hat211,*f_hat211,*g_hat121,*f_hat121,*g_hat221,*f_hat221;
  C *g_hat112,*f_hat112,*g_hat212,*f_hat212,*g_hat122,*f_hat122,*g_hat222,*f_hat222;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  N2=ths->N[2];
  n0=ths->n[0];
  n1=ths->n[1];
  n2=ths->n[2];
  p1=n1+ths->halo_hi; /* row lengths of g_hat */
  p2=n2+ths->halo_hi;

  f_hat=(C*)ths->f_hat;
  g_hat=(C*)ths->g_hat;

  TIC(0)
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(k0)
  for (k0 = 0; k0 < grid_fftw_total(ths); k0++)
    ths->g_hat[k0] = 0.0;
#else
  memset(ths->g_hat, 0, (size_t)(grid_fftw_total(ths)) * sizeof(C));
#endif

  if(ths->flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,c_phi_inv11,c_phi_inv12,ck11,ck12,c_phi_inv21,c_phi_inv22,g_hat111,f_hat111,g_hat211,f_hat211,g_hat121,f_hat121,g_hat221,f_hat221,g_hat112,f_hat112,g_hat212,f_hat212,g_hat122,f_hat122,g_hat222,f_hat222,ck21,ck22)
#endif
      for(k0=0;k0<N0/2;k0++)
  {
    ck01=c_phi_inv01[k0];
    ck02=c_phi_inv02[k0];
    c_phi_inv11=ths->c_phi_inv[1];
    c_phi_inv12=&ths->c_phi_inv[1][N1/2];

    for(k1=0;k1<N1/2;k1++)
      {
        ck11=c_phi_inv11[k1];
        ck12=c_phi_inv12[k1];
        c_phi_inv21=ths->c_phi_inv[2];
        c_phi_inv22=&ths->c_phi_inv[2][N2/2];

        g_hat111=g_hat + ((n0-(N0/2)+k0)*p1+n1-(N1/2)+k1)*p2+n2-(N2/2);
        f_hat111=f_hat + (k0*N1+k1)*N2;
        g_hat211=g_hat + (k0*p1+n1-(N1/2)+k1)*p2+n2-(N2/2);
        f_hat211=f_hat + (((N0/2)+k0)*N1+k1)*N2;
        g_hat121=g_hat + ((n0-(N0/2)+k0)*p1+k1)*p2+n2-(N2/2);
        f_hat121=f_hat + (k0*N1+(N1/2)+k1)*N2;
        g_hat221=g_hat + (k0*p1+k1)*p2+n2-(N2/2);
        f_hat221=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2;

        g_hat112=g_hat + ((n0-(N0/2)+k0)*p1+n1-(N1/2)+k1)*p2;
        f_hat112=f_hat + (k0*N1+k1)*N2+(N2/2);
        g_hat212=g_hat + (k0*p1+n1-(N1/2)+k1)*p2;
        f_hat212=f_hat + (((N0/2)+k0)*N1+k1)*N2+(N2/2);
        g_hat122=g_hat + ((n0-(N0/2)+k0)*p1+k1)*p2;
        f_hat122=f_hat + (k0*N1+N1/2+k1)*N2+(N2/2);
        g_hat222=g_hat + (k0*p1+k1)*p2;
        f_hat222=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2+(N2/2);

        for(k2=0;k2<N2/2;k2++)
    {
      ck21=c_phi_inv21[k2];
      ck22=c_phi_inv22[k2];

      g_hat111[k2] = f_hat111[k2] * ck01 * ck11 * ck21;
      g_hat211[k2] = f_hat211[k2] * ck02 * ck11 * ck21;
      g_hat121[k2] = f_hat121[k2] * ck01 * ck12 * ck21;
      g_hat221[k2] = f_hat221[k2] * ck02 * ck12 * ck21;

      g_hat112[k2] = f_hat112[k2] * ck01 * ck11 * ck22;
      g_hat212[k2] = f_hat212[k2] * ck02 * ck11 * ck22;
      g_hat122[k2] = f_hat122[k2] * ck01 * ck12 * ck22;
      g_hat222[k2] = f_hat222[k2] * ck02 * ck12 * ck22;
    }
      }
  }
    }
  else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,ck11,ck12,ck21,ck22)
#endif
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(ths->n[0],k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(ths->n[0],k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(ths->n[1],k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(ths->n[1],k1,1));

      for(k2=0;k2<N2/2;k2++)
        {
    ck21=K(1.0)/(PHI_HUT(ths->n[2],k2-N2/2,2));
    ck22=K(1.0)/(PHI_HUT(ths->n[2],k2,2));

    g_hat[((n0-N0/2+k0)*p1+n1-N1/2+k1)*p2+n2-N2/2+k2] = f_hat[(k0*N1+k1)*N2+k2]                  * ck01 * ck11 * ck21;
    g_hat[(k0*p1+n1-N1/2+k1)*p2+n2-N2/2+k2]           = f_hat[((N0/2+k0)*N1+k1)*N2+k2]           * ck02 * ck11 * ck21;
    g_hat[((n0-N0/2+k0)*p1+k1)*p2+n2-N2/2+k2]         = f_hat[(k0*N1+N1/2+k1)*N2+k2]             * ck01 * ck12 * ck21;
    g_hat[(k0*p1+k1)*p2+n2-N2/2+k2]                   = f_hat[((N0/2+k0)*N1+N1/2+k1)*N2+k2]      * ck02 * ck12 * ck21;

    g_hat[((n0-N0/2+k0)*p1+n1-N1/2+k1)*p2+k2]         = f_hat[(k0*N1+k1)*N2+N2/2+k2]             * ck01 * ck11 * ck22;
    g_hat[(k0*p1+n1-N1/2+k1)*p2+k2]                   = f_hat[((N0/2+k0)*N1+k1)*N2+N2/2+k2]      * ck02 * ck11 * ck22;
    g_hat[((n0-N0/2+k0)*p1+k1)*p2+k2]                 = f_hat[(k0*N1+N1/2+k1)*N2+N2/2+k2]        * ck01 * ck12 * ck22;
    g_hat[(k0*p1+k1)*p2+k2]                           = f_hat[((N0/2+k0)*N1+N1/2+k1)*N2+N2/2+k2] * ck02 * ck12 * ck22;
        }
    }
      }

  TOC(0)

  TIC_FFTW(1)
  FFTW(execute)(ths->my_fftw_plan1);
  TOC_FFTW(1);
  grid_B_trafo(ths);

  TIC(2);
  nfft_trafo_3d_B(ths);
  scatter_f(ths);
  TOC(2);
}

void X(adjoint_3d)(X(plan) *ths)
{
  if((ths->N[0] <= ths->m) || (ths->N[1] <= ths->m) || (ths->N[2] <= ths->m) || (ths->n[0] <= 2*ths->m+2) || (ths->n[1] <= 2*ths->m+2) || (ths->n[2] <= 2*ths->m+2))
  {
    X(adjoint_direct)(ths);
    return;
  }
  
  INT k0,k1,k2,n0,n1,n2,N0,N1,N2,p1,p2;
  C *g_hat,*f_hat;
  R *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12, *c_phi_inv21, *c_phi_inv22;
  R ck01, ck02, ck11, ck12, ck21, ck22;
  C *g_hat111,*f_hat111,*g_hat211,*f_hat211,*g_hat121,*f_hat121,*g_hat221,*f_hat221;
  C *g_hat112,*f_hat112,*g_hat212,*f_hat212,*g_hat122,*f_hat122,*g_hat222,*f_hat222;

  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  N0=ths->N[0];
  N1=ths->N[1];
  N2=ths->N[2];
  n0=ths->n[0];
  n1=ths->n[1];
  n2=ths->n[2];
  p1=n1+ths->halo_hi; /* row lengths of g_hat */
  p2=n2+ths->halo_hi;

  f_hat=(C*)ths->f_hat;
  g_hat=(C*)ths->g_hat;

  TIC(2);
  select_grid_B(ths);
  gather_f(ths);
  nfft_adjoint_3d_B(ths);
  grid_B_adjoint_done(ths);
  TOC(2);

  TIC_FFTW(1)
  FFTW(execute)(ths->my_fftw_plan2);
  TOC_FFTW(1);

  TIC(0)
  if(ths->flags & PRE_PHI_HUT)
    {
      c_phi_inv01=ths->c_phi_inv[0];
      c_phi_inv02=&ths->c_phi_inv[0][N0/2];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,c_phi_inv11,c_phi_inv12,ck11,ck12,c_phi_inv21,c_phi_inv22,g_hat111,f_hat111,g_hat211,f_hat211,g_hat121,f_hat121,g_hat221,f_hat221,g_hat112,f_hat112,g_hat212,f_hat212,g_hat122,f_hat122,g_hat222,f_hat222,ck21,ck22)
#endif
      for(k0=0;k0<N0/2;k0++)
  {
    ck01=c_phi_inv01[k0];
    ck02=c_phi_inv02[k0];
    c_phi_inv11=ths->c_phi_inv[1];
    c_phi_inv12=&ths->c_phi_inv[1][N1/2];

    for(k1=0;k1<N1/2;k1++)
      {
        ck11=c_phi_inv11[k1];
        ck12=c_phi_inv12[k1];
        c_phi_inv21=ths->c_phi_inv[2];
        c_phi_inv22=&ths->c_phi_inv[2][N2/2];

        g_hat111=g_hat + ((n0-(N0/2)+k0)*p1+n1-(N1/2)+k1)*p2+n2-(N2/2);
        f_hat111=f_hat + (k0*N1+k1)*N2;
        g_hat211=g_hat + (k0*p1+n1-(N1/2)+k1)*p2+n2-(N2/2);
        f_hat211=f_hat + (((N0/2)+k0)*N1+k1)*N2;
        g_hat121=g_hat + ((n0-(N0/2)+k0)*p1+k1)*p2+n2-(N2/2);
        f_hat121=f_hat + (k0*N1+(N1/2)+k1)*N2;
        g_hat221=g_hat + (k0*p1+k1)*p2+n2-(N2/2);
        f_hat221=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2;

        g_hat112=g_hat + ((n0-(N0/2)+k0)*p1+n1-(N1/2)+k1)*p2;
        f_hat112=f_hat + (k0*N1+k1)*N2+(N2/2);
        g_hat212=g_hat + (k0*p1+n1-(N1/2)+k1)*p2;
        f_hat212=f_hat + (((N0/2)+k0)*N1+k1)*N2+(N2/2);
        g_hat122=g_hat + ((n0-(N0/2)+k0)*p1+k1)*p2;
        f_hat122=f_hat + (k0*N1+(N1/2)+k1)*N2+(N2/2);
        g_hat222=g_hat + (k0*p1+k1)*p2;
        f_hat222=f_hat + (((N0/2)+k0)*N1+(N1/2)+k1)*N2+(N2/2);

        for(k2=0;k2<N2/2;k2++)
    {
      ck21=c_phi_inv21[k2];
      ck22=c_phi_inv22[k2];

      f_hat111[k2] = g_hat111[k2] * ck01 * ck11 * ck21;
      f_hat211[k2] = g_hat211[k2] * ck02 * ck11 * ck21;
      f_hat121[k2] = g_hat121[k2] * ck01 * ck12 * ck21;
      f_hat221[k2] = g_hat221[k2] * ck02 * ck12 * ck21;

      f_hat112[k2] = g_hat112[k2] * ck01 * ck11 * ck22;
      f_hat212[k2] = g_hat212[k2] * ck02 * ck11 * ck22;
      f_hat122[k2] = g_hat122[k2] * ck01 * ck12 * ck22;
      f_hat222[k2] = g_hat222[k2] * ck02 * ck12 * ck22;
    }
      }
  }
    }
  else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(k0,k1,k2,ck01,ck02,ck11,ck12,ck21,ck22)
#endif
    for(k0=0;k0<N0/2;k0++)
      {
  ck01=K(1.0)/(PHI_HUT(ths->n[0],k0-N0/2,0));
  ck02=K(1.0)/(PHI_HUT(ths->n[0],k0,0));
  for(k1=0;k1<N1/2;k1++)
    {
      ck11=K(1.0)/(PHI_HUT(ths->n[1],k1-N1/2,1));
      ck12=K(1.0)/(PHI_HUT(ths->n[1],k1,1));

      for(k2=0;k2<N2/2;k2++)
        {
    ck21=K(1.0)/(PHI_HUT(ths->n[2],k2-N2/2,2));
    ck22=K(1.0)/(PHI_HUT(ths->n[2],k2,2));

    f_hat[(k0*N1+k1)*N2+k2]                  = g_hat[((n0-N0/2+k0)*p1+n1-N1/2+k1)*p2+n2-N2/2+k2] * ck01 * ck11 * ck21;
    f_hat[((N0/2+k0)*N1+k1)*N2+k2]           = g_hat[(k0*p1+n1-N1/2+k1)*p2+n2-N2/2+k2]           * ck02 * ck11 * ck21;
    f_hat[(k0*N1+N1/2+k1)*N2+k2]             = g_hat[((n0-N0/2+k0)*p1+k1)*p2+n2-N2/2+k2]         * ck01 * ck12 * ck21;
    f_hat[((N0/2+k0)*N1+N1/2+k1)*N2+k2]      = g_hat[(k0*p1+k1)*p2+n2-N2/2+k2]                   * ck02 * ck12 * ck21;

    f_hat[(k0*N1+k1)*N2+N2/2+k2]             = g_hat[((n0-N0/2+k0)*p1+n1-N1/2+k1)*p2+k2]         * ck01 * ck11 * ck22;
    f_hat[((N0/2+k0)*N1+k1)*N2+N2/2+k2]      = g_hat[(k0*p1+n1-N1/2+k1)*p2+k2]                   * ck02 * ck11 * ck22;
    f_hat[(k0*N1+N1/2+k1)*N2+N2/2+k2]        = g_hat[((n0-N0/2+k0)*p1+k1)*p2+k2]                 * ck01 * ck12 * ck22;
    f_hat[((N0/2+k0)*N1+N1/2+k1)*N2+N2/2+k2] = g_hat[(k0*p1+k1)*p2+k2]                           * ck02 * ck12 * ck22;
        }
    }
      }

  TOC(0)
}

/** user routines
 */
void X(trafo)(X(plan) *ths)
{
  /* use direct transform if degree N is too low */
  for (int j = 0; j < ths->d; j++)
  {
    if((ths->N[j] <= ths->m) || (ths->n[j] <= 2*ths->m+2))
    {
      X(trafo_direct)(ths);
      return;
    }
  }
  
  switch(ths->d)
  {
    case 1: X(trafo_1d)(ths); break;
    case 2: X(trafo_2d)(ths); break;
    case 3: X(trafo_3d)(ths); break;
    default:
    {
      /* use ths->my_fftw_plan1 */
      ths->g_hat = ths->g1;
      ths->g = ths->g2;

      /** form \f$ \hat g_k = \frac{\hat f_k}{c_k\left(\phi\right)} \text{ for }
       *  k \in I_N \f$
       */
      TIC(0)
      D_A(ths);
      TOC(0)

      /** compute by d-variate discrete Fourier transform
       *  \f$ g_l = \sum_{k \in I_N} \hat g_k {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
       *  \text{ for } l \in I_n \f$
       */
      TIC_FFTW(1)
      FFTW(execute)(ths->my_fftw_plan1);
      TOC_FFTW(1)
      grid_B_trafo(ths);

      /** set \f$ f_j =\sum_{l \in I_n,m(x_j)} g_l \psi\left(x_j-\frac{l}{n}\right)
       *  \text{ for } j=0,\dots,M_total-1 \f$
       */
      TIC(2)
      B_A(ths);
      scatter_f(ths);
      TOC(2)
    }
  }
} /* nfft_trafo */

void X(adjoint)(X(plan) *ths)
{
  /* use direct transform if degree N is too low */
  for (int j = 0; j < ths->d; j++)
  {
    if((ths->N[j] <= ths->m) || (ths->n[j] <= 2*ths->m+2))
    {
      X(adjoint_direct)(ths);
      return;
    }
  }
  
  switch(ths->d)
  {
    case 1: X(adjoint_1d)(ths); break;
    case 2: X(adjoint_2d)(ths); break;
    case 3: X(adjoint_3d)(ths); break;
    default:
    {
      /* use ths->my_fftw_plan2 */
      ths->g_hat=ths->g1;
      ths->g=ths->g2;
      
      /** set \f$ g_l = \sum_{j=0}^{M_total-1} f_j \psi\left(x_j-\frac{l}{n}\right)
       *  \text{ for } l \in I_n,m(x_j) \f$
       */
      TIC(2)
      select_grid_B(ths);
      gather_f(ths);
      B_T(ths);
      grid_B_adjoint_done(ths);
      TOC(2)

      /** compute by d-variate discrete Fourier transform
       *  \f$ \hat g_k = \sum_{l \in I_n} g_l {\rm e}^{+2\pi {\rm i} \frac{kl}{n}}
       *  \text{ for }  k \in I_N\f$
       */
      TIC_FFTW(1)
      FFTW(execute)(ths->my_fftw_plan2);
      TOC_FFTW(1)

      /** form \f$ \hat f_k = \frac{\hat g_k}{c_k\left(\phi\right)} \text{ for }
       *  k \in I_N \f$
       */
      TIC(0)
      D_T(ths);
      TOC(0)
    }
  }
} /* nfft_adjoint */


/** initialisation of direct transform
 */
static void precompute_phi_hut(X(plan) *ths)
{
  INT ks[ths->d]; /* index over all frequencies */
  INT t; /* index over all dimensions */

  ths->c_phi_inv = (R**) Y(malloc)((size_t)(ths->d) * sizeof(R*));

  for (t = 0; t < ths->d; t++)
  {
    ths->c_phi_inv[t] = (R*)Y(malloc)((size_t)(ths->N[t]) * sizeof(R));

    for (ks[t] = 0; ks[t] < ths->N[t]; ks[t]++)
    {
      ths->c_phi_inv[t][ks[t]]= K(1.0) / (PHI_HUT(ths->n[t], ks[t] - ths->N[t] / 2,t));
    }
  }
} /* nfft_phi_hut */

/** create a lookup table, but NOT for each node
 *  good idea K=2^xx
 *  assumes an EVEN window function
 */
void X(precompute_lin_psi)(X(plan) *ths)
{
  INT t;                                /**< index over all dimensions       */
  INT j;                                /**< index over all nodes            */
  R step;                          /**< step size in [0,(m+2)/n]        */

  for (t=0; t<ths->d; t++)
    {
      step = ((R)(ths->m+2)) / ((R)(ths->K * ths->n[t]));
      for(j = 0;j <= ths->K; j++)
  {
    ths->psi[(ths->K+1)*t + j] = PHI(ths->n[t], (R)(j) * step,t);
  } /* for(j) */
    } /* for(t) */
}

static void precompute_fg_psi_sorted(X(plan) *ths)
{
  INT t;                                /**< index over all dimensions       */
  const R *xs;

  xs = node_x(ths);

  for (t=0; t<ths->d; t++)
  {
    INT j;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j)
#endif
    for (j = 0; j < ths->M_total; j++)
      {
        /* anchored at the grid point nearest the node, offset d from it. The
         * seed carries the sign of d, which is all FG_RUN_START needs to
         * place the run; a node on a grid point has d = 0 and a positive
         * seed. FLOOR rather than LRINT keeps this loop vectorisable. */
        const R nx = (R)(ths->n[t]) * xs[j*ths->d+t];
        const R d = NX_SUB(ths->n[t], xs[j*ths->d+t], FLOOR(nx + K(0.5)));
        const R seed = (PHI(ths->n[t], d / (R)(ths->n[t]), t));

        ths->psi[2*(j*ths->d+t)] = (d < K(0.0)) ? -seed : seed;
        ths->psi[2*(j*ths->d+t)+1] = EXP(K(2.0) * d / ths->b[t]);
      } /* for(j) */
  }
  /* for(t) */
} /* nfft_precompute_fg_psi */

void X(precompute_fg_psi)(X(plan) *ths)
{
  sort_nodes(ths);
  precompute_fg_psi_sorted(ths);
}

static void precompute_psi_sorted(X(plan) *ths)
{
  INT t; /* index over all dimensions */
  INT u, o; /* depends on x_j */
  const R *xs;

  xs = node_x(ths);

  for (t=0; t<ths->d; t++)
  {
    INT j;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j,u,o)
#endif
    for (j = 0; j < ths->M_total; j++)
    {
      uo(ths, xs, j, &u,&o,t);

      PHI_RUN(ths->psi + (j * ths->d + t) * (2 * ths->m + 2), ths->n[t],
          xs[j*ths->d+t], u, t);
    } /* for(j) */
  }
  /* for(t) */
} /* nfft_precompute_psi */

void X(precompute_psi)(X(plan) *ths)
{
  sort_nodes(ths);
  precompute_psi_sorted(ths);
}

#ifdef _OPENMP
static void nfft_precompute_full_psi_omp(X(plan) *ths)
{
  INT j;                                /**< index over all nodes            */
  INT lprod;                            /**< 'bandwidth' of matrix B         */
  const R *xs = node_x(ths);

  {
    INT t;
    for(t=0,lprod = 1; t<ths->d; t++)
        lprod *= 2*ths->m+2;
  }

  #pragma omp parallel for default(shared) private(j)
  for(j=0; j<ths->M_total; j++)
    {
      INT t,t2;                             /**< index over all dimensions       */
      INT l_L;                              /**< plain index 0<=l_L<lprod        */
      INT lj[ths->d];                       /**< multi index 0<=lj<u+o+1         */
      INT ll_plain[ths->d+1];               /**< postfix plain index             */

      INT u[ths->d], o[ths->d];             /**< depends on x_j                  */

      R phi_prod[ths->d+1];
      INT ix = j*lprod;

      phi_prod[0]=1;
      ll_plain[0]=0;

      R psij_const[ths->d * (2*ths->m+2)];

      MACRO_init_uo_l_lj_t;

      for (t2 = 0; t2 < ths->d; t2++)
        PHI_RUN(psij_const + t2 * (2*ths->m+2), ths->n[t2],
            xs[j*ths->d+t2], u[t2], t2);

      for(l_L=0; l_L<lprod; l_L++, ix++)
      {
        MACRO_update_phi_prod_ll_plain(without_PRE_PSI_improved);

        ths->psi_index_g[ix]=ll_plain[ths->d];
        ths->psi[ix]=phi_prod[ths->d];

        MACRO_count_uo_l_lj_t;
      } /* for(l_L) */

      ths->psi_index_f[j]=lprod;
    } /* for(j) */
}
#endif

static void precompute_full_psi_sorted(X(plan) *ths)
{
#ifdef _OPENMP
  nfft_precompute_full_psi_omp(ths);
#else
  INT t, t2; /* index over all dimensions */
  INT j; /* index over all nodes */
  INT l_L; /* plain index 0 <= l_L < lprod */
  INT lj[ths->d]; /* multi index 0<=lj<u+o+1 */
  INT ll_plain[ths->d+1]; /* postfix plain index */
  INT lprod; /* 'bandwidth' of matrix B */
  INT u[ths->d], o[ths->d]; /* depends on x_j */

  R phi_prod[ths->d+1];

  INT ix, ix_old;
  const R *xs;

  xs = node_x(ths);

  phi_prod[0] = K(1.0);
  ll_plain[0] = 0;

  for (t = 0, lprod = 1; t < ths->d; t++)
    lprod *= 2 * ths->m + 2;

  for (j = 0, ix = 0, ix_old = 0; j < ths->M_total; j++)
  {
    R psij_const[ths->d * (2 * ths->m + 2)];

    MACRO_init_uo_l_lj_t;

    for (t2 = 0; t2 < ths->d; t2++)
      PHI_RUN(psij_const + t2 * (2 * ths->m + 2), ths->n[t2],
          xs[j * ths->d + t2], u[t2], t2);

    for (l_L = 0; l_L < lprod; l_L++, ix++)
    {
      MACRO_update_phi_prod_ll_plain(without_PRE_PSI_improved);

      ths->psi_index_g[ix] = ll_plain[ths->d];
      ths->psi[ix] = phi_prod[ths->d];

      MACRO_count_uo_l_lj_t;
    } /* for(l_L) */

    ths->psi_index_f[j] = ix - ix_old;
    ix_old = ix;
  } /* for(j) */
#endif
}

void X(precompute_full_psi)(X(plan) *ths)
{
  sort_nodes(ths);
  precompute_full_psi_sorted(ths);
}

void X(precompute_one_psi)(X(plan) *ths)
{
  /* one snapshot of the nodes for all psi tables (NFFT_SORT_NODES only) */
  sort_nodes(ths);

  if(ths->flags & PRE_LIN_PSI)
    X(precompute_lin_psi)(ths);
  if(ths->flags & PRE_FG_PSI)
    precompute_fg_psi_sorted(ths);
  if(ths->flags & PRE_PSI)
    precompute_psi_sorted(ths);
  if(ths->flags & PRE_FULL_PSI)
    precompute_full_psi_sorted(ths);
}

static void init_help(X(plan) *ths)
{
  INT t; /* index over all dimensions */
  INT lprod; /* 'bandwidth' of matrix B */

  if (ths->flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    ths->flags |= NFFT_SORT_NODES;

  ths->N_total = intprod(ths->N, 0, ths->d);
  ths->n_total = intprod(ths->n, 0, ths->d);

  ths->sigma = (R*) Y(malloc)((size_t)(ths->d) * sizeof(R));

  for(t = 0;t < ths->d; t++)
    ths->sigma[t] = ((R)ths->n[t]) / (R)(ths->N[t]);

  WINDOW_HELP_INIT;

  if(ths->flags & MALLOC_X)
    ths->x = (R*)Y(malloc)((size_t)(ths->d * ths->M_total) * sizeof(R));

  if(ths->flags & MALLOC_F_HAT)
    ths->f_hat = (C*)Y(malloc)((size_t)(ths->N_total) * sizeof(C));

  if(ths->flags & MALLOC_F)
    ths->f = (C*)Y(malloc)((size_t)(ths->M_total) * sizeof(C));

  if(ths->flags & PRE_PHI_HUT)
    precompute_phi_hut(ths);

  if (ths->flags & PRE_LIN_PSI)
  {
      if (ths->K == 0)
      {
        ths->K = Y(m2K)(ths->m);
      }
      ths->psi = (R*) Y(malloc)((size_t)((ths->K+1) * ths->d) * sizeof(R));
  }

  if(ths->flags & PRE_FG_PSI)
    ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * ths->d * 2) * sizeof(R));

  if(ths->flags & PRE_PSI)
    ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * ths->d * (2 * ths->m + 2)) * sizeof(R));

  if(ths->flags & PRE_FULL_PSI)
  {
      for (t = 0, lprod = 1; t < ths->d; t++)
        lprod *= 2 * ths->m + 2;

      ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * lprod) * sizeof(R));

      ths->psi_index_f = (INT*) Y(malloc)((size_t)(ths->M_total) * sizeof(INT));
      ths->psi_index_g = (INT*) Y(malloc)((size_t)(ths->M_total * lprod) * sizeof(INT));
  }

  ths->halo_hi = 0;

  if(ths->flags & FFTW_INIT)
  {
#ifdef _OPENMP
    INT nthreads = Y(get_num_threads)();
#endif

    /* g1 and g2 are the padded grid of B, the FFT works on their
     * n_0 x ... x n_{d-1} part */
    ths->halo_hi = HALO(ths->m);

    ths->g1 = (C*)Y(malloc)((size_t)(grid_fftw_total(ths)) * sizeof(C));

    if(ths->flags & FFT_OUT_OF_PLACE)
      ths->g2 = (C*) Y(malloc)((size_t)(grid_fftw_total(ths)) * sizeof(C));
    else
      ths->g2 = ths->g1;

#if defined(_OPENMP) && defined(HAVE_FFTW_THREADS)
#pragma omp critical (nfft_omp_critical_fftw_plan)
{
    FFTW(plan_with_nthreads)(nthreads);
#endif
    {
      FFTW(iodim64) *dims = Y(malloc)((size_t)(ths->d) * sizeof(FFTW(iodim64)));
      ptrdiff_t stride = 1;

      for (t = ths->d - 1; t >= 0; t--)
      {
        dims[t].n = (ptrdiff_t)(ths->n[t]);
        dims[t].is = stride;
        dims[t].os = stride;
        stride *= (ptrdiff_t)(ths->n[t] + ths->halo_hi);
      }

      ths->my_fftw_plan1 = FFTW(plan_guru64_dft)((int)ths->d, dims, 0, NULL,
          ths->g1, ths->g2, FFTW_FORWARD, ths->fftw_flags);
      ths->my_fftw_plan2 = FFTW(plan_guru64_dft)((int)ths->d, dims, 0, NULL,
          ths->g2, ths->g1, FFTW_BACKWARD, ths->fftw_flags);
      Y(free)(dims);
    }
#if defined(_OPENMP) && defined(HAVE_FFTW_THREADS)
}
#endif
  }

  if(ths->flags & NFFT_SORT_NODES)
    ths->index_x = (INT*) Y(malloc)(sizeof(INT) * 2U * (size_t)(ths->M_total));
  else
    ths->index_x = NULL;
  ths->x_sorted = NULL;
  ths->tile_ptr = NULL;
  ths->f_sorted = NULL;
  ths->g_pad = NULL;
  ths->g_pad_size = 0;

  ths->mv_trafo = (void (*) (void* ))X(trafo);
  ths->mv_adjoint = (void (*) (void* ))X(adjoint);
}

void X(init)(X(plan) *ths, int d, int *N, int M_total)
{
  INT t; /* index over all dimensions */

  ths->d = (INT)d;

  ths->N = (INT*) Y(malloc)((size_t)(d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->N[t] = (INT)N[t];

  ths->M_total = (INT)M_total;

  ths->n = (INT*) Y(malloc)((size_t)(d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->n[t] = 2 * (Y(next_power_of_2)(ths->N[t]));

  ths->m = WINDOW_HELP_ESTIMATE_m;

  if (d > 1)
  {
#ifdef _OPENMP
    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | NFFT_SORT_NODES |
                 NFFT_OMP_BLOCKWISE_ADJOINT;
#else
    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | NFFT_SORT_NODES;
#endif
  }
  else
  {
    /* In 1D the unsorted OpenMP adjoint has to use atomic operations, the
     * sorted one splits the grid among the threads. Serially, sorting costs
     * more in the precomputation than the 1D transforms gain from it, unless
     * the plan is used for about ten transforms or more. */
#ifdef _OPENMP
    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | FFT_OUT_OF_PLACE | NFFT_SORT_NODES |
                 NFFT_OMP_BLOCKWISE_ADJOINT;
#else
    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | FFT_OUT_OF_PLACE;
#endif
  }

  ths->fftw_flags= FFTW_ESTIMATE| FFTW_DESTROY_INPUT;

  ths->K = 0;
  init_help(ths);
}

void X(init_guru)(X(plan) *ths, int d, int *N, int M_total, int *n, int m,
  unsigned flags, unsigned fftw_flags)
{
  INT t; /* index over all dimensions */

  ths->d = (INT)d;
  ths->M_total = (INT)M_total;
  ths->N = (INT*)Y(malloc)((size_t)(ths->d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->N[t] = (INT)N[t];

  ths->n = (INT*)Y(malloc)((size_t)(ths->d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->n[t] = (INT)n[t];

  ths->m = (INT)m;

  ths->flags = flags;
  ths->fftw_flags = fftw_flags;

  ths->K = 0;
  init_help(ths);
}

void X(init_lin)(X(plan) *ths, int d, int *N, int M_total, int *n, int m, int K,
  unsigned flags, unsigned fftw_flags)
{
  INT t; /* index over all dimensions */

  ths->d = (INT)d;
  ths->M_total = (INT)M_total;
  ths->N = (INT*)Y(malloc)((size_t)(ths->d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->N[t] = (INT)N[t];

  ths->n = (INT*)Y(malloc)((size_t)(ths->d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->n[t] = (INT)n[t];

  ths->m = (INT)m;

  ths->flags = flags;
  ths->fftw_flags = fftw_flags;

  ths->K = K;
  init_help(ths);
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
  {
    if ((ths->x[j]<-K(0.5)) || (ths->x[j]>= K(0.5)))
    {
      return "ths->x out of range [-0.5,0.5)";
    }
  }

  for (j = 0; j < ths->d; j++)
  {
    if (ths->sigma[j] <= 1)
      return "Oversampling factor too small";

#ifdef WINDOW_IS_KAISER_BESSEL
    if (!KB_POLY_SIGMA_OK(ths->sigma[j]))
      return "Oversampling factor outside [1.1, 4]";
#endif
    
    /* Automatically calls trafo_direct if 
    if(ths->N[j] <= ths->m)
      return "Polynomial degree N is <= cut-off m";
    */

    if(ths->N[j]%2 == 1)
      return "polynomial degree N has to be even";
  }
  return 0;
}

void X(finalize)(X(plan) *ths)
{
  INT t; /* index over dimensions */

  if(ths->flags & NFFT_SORT_NODES)
    Y(free)(ths->index_x);

  if (ths->x_sorted != NULL)
  {
    Y(free)(ths->x_sorted);
    Y(free)(ths->f_sorted);
  }

  if (ths->tile_ptr != NULL)
    Y(free)(ths->tile_ptr);

  if (ths->g_pad != NULL)
    Y(free)(ths->g_pad);

  if(ths->flags & FFTW_INIT)
  {
#ifdef _OPENMP
    #pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
    FFTW(destroy_plan)(ths->my_fftw_plan2);
#ifdef _OPENMP
    #pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
    FFTW(destroy_plan)(ths->my_fftw_plan1);

    if(ths->flags & FFT_OUT_OF_PLACE)
      Y(free)(ths->g2);

    Y(free)(ths->g1);
  }

  if(ths->flags & PRE_FULL_PSI)
  {
    Y(free)(ths->psi_index_g);
    Y(free)(ths->psi_index_f);
    Y(free)(ths->psi);
  }

  if(ths->flags & PRE_PSI)
    Y(free)(ths->psi);

  if(ths->flags & PRE_FG_PSI)
    Y(free)(ths->psi);

  if(ths->flags & PRE_LIN_PSI)
    Y(free)(ths->psi);

  if(ths->flags & PRE_PHI_HUT)
  {
    for (t = 0; t < ths->d; t++)
        Y(free)(ths->c_phi_inv[t]);
    Y(free)(ths->c_phi_inv);
  }

  if(ths->flags & MALLOC_F)
    Y(free)(ths->f);

  if(ths->flags & MALLOC_F_HAT)
    Y(free)(ths->f_hat);

  if(ths->flags & MALLOC_X)
    Y(free)(ths->x);

  WINDOW_HELP_FINALIZE;

  Y(free)(ths->sigma);
  Y(free)(ths->n);
  Y(free)(ths->N);
}
