/*
 * Copyright (c) 2002, 2020 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* Nonequispaced fast cosine transform */

/* Author: Steffen Klatt 2004-2006, Jens Keiner 2010 */

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
#define X(name) NFST(name)

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
#define BASE(x) SIN(x)
#define NN(x) (x + 1)
#define OFFSET 1
#define FOURIER_TRAFO FFTW_RODFT00
#define FFTW_DEFAULT_FLAGS FFTW_ESTIMATE | FFTW_DESTROY_INPUT

#define NODE(p,r) (xs[(p) * ths->d + (r)])

/* Block size for the phase recurrence in the direct transforms */
#define NFFT_DIRECT_RECURRENCE_BLOCK 32
#define TRAFO_STEP(i, w) (v += f_hat[i] * BASEPART(w))
#define ADJOINT_STEP(i, w) (f_hat[i] += f[j] * BASEPART(w))
#define ROW_TRAFO_STEP(i, w) (s += fr[i] * BASEPART(w))
#define ROW_ADJOINT_STEP(i, w) (fr[i] += a * BASEPART(w))

/* Rows share their block seeds, since the outer axes enter through P; seeds of the first
 * NFFT_DIRECT_ROW_SEEDS blocks are computed once per node, later ones per row. */
#define NFFT_DIRECT_ROW_SEEDS 8

/* Minimum innermost-dimension length for the multivariate recurrence to pay. */
#define NFFT_DIRECT_RECURRENCE_MIN_INNER 2


/* Threaded adjoint: rows shorter than NFFT_DIRECT_SHORT_ROW, and fewer than
 * NFFT_DIRECT_ROWS_PER_THREAD per thread, are shared out by frequency instead of by row
 * to avoid cache contention. */
#define NFFT_DIRECT_SHORT_ROW 8
#define NFFT_DIRECT_ROWS_PER_THREAD 4

/* The (co)sine value is a part of the complex phase exp(+i 2pi (k+OFFSET) x): NDCT reads the
 * real part (cos), NDST reads the imaginary part (sin). */
#define BASEPART(w) CIMAG(w)

/* Accurate phase for exp(+i 2pi k x): reduce k*x modulo 1 into ~[-1/2,1/2) so COS/SIN see a
 * small argument, error does not grow with N. Requires FMA single-rounding semantics. */
static inline R X(reduced_omega)(const R k, const R x)
{
  const R n = RINT(k * x);      // Nearest integer to k * x.
  return K2PI * FFMA(k, x, -n); // Calculate k * x - n with a single rounding, then multiply 2 * pi.
}

/* Seed of block b of a row of length nl: the factor at the block centre. */
static inline C X(row_seed)(const INT b, const INT nl, const R xl)
{
  const INT i = b * NFFT_DIRECT_RECURRENCE_BLOCK;
  const INT len = MIN(NFFT_DIRECT_RECURRENCE_BLOCK, nl - i);
  const R omega = X(reduced_omega)((R)(i + len/2 + OFFSET), xl);
  return COS(omega) + II * SIN(omega);
}

/* BASE(2pi (k + OFFSET) x) from the reduced phase. */
static inline R X(base_reduced)(const INT k, const R x)
{
  return BASE(X(reduced_omega)((R)(k + OFFSET), x));
}

#define MACRO_without_PRE_PSI psij_const[t * (2 * ths->m + 2) + lj[t]]

/* Window runs of node j on all axes, for MACRO_without_PRE_PSI. */
#define MACRO_fill_psij_const \
{ \
  for (t = 0; t < ths->d; t++) \
    PHI_RUN(psij_const + t * (2 * ths->m + 2), (2 * NN(ths->n[t])), \
        xs[j * ths->d + t], u[t], t); \
}
#define MACRO_compute_PSI PHI((2 * NN(ths->n[t])), \
  NX_SUB((2 * NN(ths->n[t])), NODE(j,t), lj[t] + u[t]) \
  / ((R)(2 * NN(ths->n[t]))), t)

/**
 * Direct computation of non equispaced sine transforms
 *  nfst_trafo_direct,  nfst_adjoint_direct
 *  require O(M N^d) arithemtical operations
 *
 * direct computation of the nfst_trafo_direct, formula (1.1)
 * nfst_trafo_direct:
 * for j=0,...,M-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * sin(2 (pi) k x[j])
 *
 * direct computation of the nfft_adjoint_direct, formula (1.2)
 * nfst_adjoint_direct:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * sin(2 (pi) k x[j])
 */
void X(trafo_direct)(const X(plan) *ths)
{
  R *f_hat = (R*)ths->f_hat, *f = (R*)ths->f;

  memset(f, 0, (size_t)(ths->M_total) * sizeof(R));

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
      R v = K(0.0);
      const R x = ths->x[j];
      const R dphi = K2PI * x;                 /* |dphi| <= pi: accurate without reduction */
      const C dw = COS(dphi) + II * SIN(dphi); /* per-step phase factor exp(+i 2pi x)      */
      INT k_L = 0;
      while (k_L < ths->N_total)
      {
        /* Accurate seed exp(+i 2pi (k + OFFSET) x) at the block centre, then recur within it. */
        const INT len = MIN(B, ths->N_total - k_L);
        const R omega = X(reduced_omega)((R)(k_L + len/2 + OFFSET), x);
        const C w0 = COS(omega) + II * SIN(omega);
        NFFT_RECUR_BLOCK(k_L, len, w0, dw, TRAFO_STEP);
        k_L += len;
      }

      f[j] = v;
    }
  }
  else
  {
    /* Multivariate: Recur phase along innermost axis and re-seed at every row start and after B steps. */
    const INT B = NFFT_DIRECT_RECURRENCE_BLOCK;
    const INT nl = ths->N[ths->d - 1] - OFFSET;
    const INT nrows = nl > 0 ? ths->N_total / nl : 0;
    const int use_rec = (nl >= NFFT_DIRECT_RECURRENCE_MIN_INNER);
    INT j;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(j)
#endif
    for (j = 0; j < ths->M_total; j++)
    {
      const R *xj = ths->x + (size_t)j * ths->d;
      const R xl = xj[ths->d - 1];
      const R dphi = K2PI * xl;
      const C dw = use_rec ? COS(dphi) + II * SIN(dphi) : K(0.0);
      R P[ths->d], v = K(0.0);
      C seed[NFFT_DIRECT_ROW_SEEDS];
      INT k[ths->d], t, t2, r, i, b;

      for (b = 0; use_rec && b < NFFT_DIRECT_ROW_SEEDS && b * B < nl; b++)
        seed[b] = X(row_seed)(b, nl, xl);

      /* P[t + 1] is the product of the outer axes' factors up to axis t. */
      P[0] = K(1.0);
      for (t = 0; t < ths->d - 1; t++)
      {
        k[t] = 0;
        P[t + 1] = P[t] * X(base_reduced)(0, xj[t]);
      }

      for (r = 0; r < nrows; r++)
      {
        const R *fr = f_hat + r * nl;
        R s = K(0.0);

        if (use_rec)
        {
          for (b = 0, i = 0; i < nl; b++, i += B)
          {
            const C w0 = b < NFFT_DIRECT_ROW_SEEDS ? seed[b] : X(row_seed)(b, nl, xl);
            NFFT_RECUR_BLOCK(i, MIN(B, nl - i), w0, dw, ROW_TRAFO_STEP);
          }
        }
        else
        {
          for (i = 0; i < nl; i++)
            s += fr[i] * X(base_reduced)(i, xl);
        }

        v += P[ths->d - 1] * s;

        if (r + 1 < nrows)
        {
          for (t = ths->d - 2; (t >= 1) && (k[t] == ths->N[t] - OFFSET - 1); t--)
            k[t] = 0;

          k[t]++;

          for (t2 = t; t2 < ths->d - 1; t2++)
            P[t2 + 1] = P[t2] * X(base_reduced)(k[t2], xj[t2]);
        }
      }

      f[j] = v;
    }
  }
}

/* Adjoint over the rows [rlo,rhi) of f_hat, a row being the N[d-1] - OFFSET contiguous
 * frequencies of the innermost axis. */
static void X(adjoint_direct_rows)(const X(plan) *ths, const INT rlo, const INT rhi)
{
  R *f_hat = (R*)ths->f_hat, *f = (R*)ths->f;
  const INT B = NFFT_DIRECT_RECURRENCE_BLOCK;
  const INT nl = ths->N[ths->d - 1] - OFFSET;
  const int use_rec = (nl >= NFFT_DIRECT_RECURRENCE_MIN_INNER);
  INT k0[ths->d], t, j, r = rlo;

  if (rlo >= rhi)
    return;

  /* Decode the outer multi-index of row rlo once; the row loop advances it from there. */
  for (t = ths->d - 2; t >= 0; t--)
  {
    k0[t] = r % (ths->N[t] - OFFSET);
    r /= ths->N[t] - OFFSET;
  }

  for (j = 0; j < ths->M_total; j++)
  {
    const R *xj = ths->x + (size_t)j * ths->d;
    const R xl = xj[ths->d - 1];
    const R dphi = K2PI * xl;
    const C dw = use_rec ? COS(dphi) + II * SIN(dphi) : K(0.0);
    R P[ths->d];
    C seed[NFFT_DIRECT_ROW_SEEDS];
    INT k[ths->d], t2, i, b;

    for (b = 0; use_rec && b < NFFT_DIRECT_ROW_SEEDS && b * B < nl; b++)
      seed[b] = X(row_seed)(b, nl, xl);

    /* P[t + 1] is f[j] times the outer axes' factors up to axis t. */
    P[0] = f[j];
    for (t = 0; t < ths->d - 1; t++)
    {
      k[t] = k0[t];
      P[t + 1] = P[t] * X(base_reduced)(k[t], xj[t]);
    }

    for (r = rlo; r < rhi; r++)
    {
      R *fr = f_hat + r * nl;
      const R a = P[ths->d - 1];

      if (use_rec)
      {
        for (b = 0, i = 0; i < nl; b++, i += B)
        {
          const C w0 = b < NFFT_DIRECT_ROW_SEEDS ? seed[b] : X(row_seed)(b, nl, xl);
          NFFT_RECUR_BLOCK(i, MIN(B, nl - i), w0, dw, ROW_ADJOINT_STEP);
        }
      }
      else
      {
        for (i = 0; i < nl; i++)
          fr[i] += a * X(base_reduced)(i, xl);
      }

      if (r + 1 < rhi)
      {
        for (t = ths->d - 2; (t >= 1) && (k[t] == ths->N[t] - OFFSET - 1); t--)
          k[t] = 0;

        k[t]++;

        for (t2 = t; t2 < ths->d - 1; t2++)
          P[t2 + 1] = P[t2] * X(base_reduced)(k[t2], xj[t2]);
      }
    }
  }
}

void X(adjoint_direct)(const X(plan) *ths)
{
  R *f_hat = (R*)ths->f_hat, *f = (R*)ths->f;

  memset(f_hat, 0, (size_t)(ths->N_total) * sizeof(R));

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
            const R omega = X(reduced_omega)((R)(k_L + len/2 + OFFSET), x);
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
          R omega = K2PI * ((R)(k_L + OFFSET)) * ths->x[j];
          f_hat[k_L] += f[j] * BASE(omega);
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
          const R omega = X(reduced_omega)((R)(k_L + len/2 + OFFSET), x);
          const C w0 = COS(omega) + II * SIN(omega);
          NFFT_RECUR_BLOCK(k_L, len, w0, dw, ADJOINT_STEP);
          k_L += len;
        }
      }
#endif
  }
  else
  {
    /* Multivariate case: Recurrence runs along a row, parallel unit is the row */
    const INT nl = ths->N[ths->d - 1] - OFFSET;
    const INT nrows = nl > 0 ? ths->N_total / nl : 0;
#ifdef _OPENMP
    if (nl >= NFFT_DIRECT_RECURRENCE_MIN_INNER && (nl >= NFFT_DIRECT_SHORT_ROW
        || nrows >= NFFT_DIRECT_ROWS_PER_THREAD * omp_get_max_threads()))
    {
      #pragma omp parallel default(shared)
      {
        const int nt = omp_get_num_threads();
        const int tid = omp_get_thread_num();
        X(adjoint_direct_rows)(ths, (INT)(((long long)nrows * tid) / nt),
                                    (INT)(((long long)nrows * (tid + 1)) / nt));
      }
    }
    else
    {
      /* Few short rows: parallelize over frequencies. */
      INT k_L;
      #pragma omp parallel for default(shared) private(k_L)
      for (k_L = 0; k_L < ths->N_total; k_L++)
      {
        INT k[ths->d], k_temp = k_L, t, j;

        for (t = ths->d - 1; t >= 0; t--)
        {
          k[t] = k_temp % (ths->N[t] - OFFSET);
          k_temp /= ths->N[t] - OFFSET;
        }

        for (j = 0; j < ths->M_total; j++)
        {
          R v = f[j];

          for (t = 0; t < ths->d; t++)
            v *= X(base_reduced)(k[t], ths->x[j * ths->d + t]);

          f_hat[k_L] += v;
        }
      }
    }
#else
    X(adjoint_direct_rows)(ths, 0, nrows);
#endif
  }
}

/** fast computation of non equispaced sine transforms
 *  require O(N^d log(N) + M) arithemtical operations
 *
 * fast computation of the nfst_trafo, formula (1.1)
 * nfst_trafo:
 * for j=0,...,M-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * sin(2 (pi) k x[j])
 *
 * direct computation of the nfst_adjoint, formula (1.2)
 * nfst_adjoint:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * sin(2 (pi) k x[j])
 */

/** macros and small sub routines for the fast transforms
 */

/** computes 2m+2 indices for the matrix B, for node j of the node array xs
 */
static inline void uo(const X(plan) *ths, const R *xs, const INT j, INT *up,
  INT *op, const INT act_dim)
{
  const R xj = xs[j * ths->d + act_dim];
  INT c = (INT)FLOOR(xj * (2 * NN(ths->n[(act_dim)])));

  (*up) = c - (ths->m);
  (*op) = c + 1 + (ths->m);
}

/** Tiles of the grid of B (see Y(tile_layout)): the cells are the cells c in
 *  [0, NN(n_t)] of the nodes in [0, 0.5] on axis t. */
static INT tile_layout(const X(plan) *ths, INT *edge, INT *ntiles)
{
  const INT d = ths->d;
  INT ext[d], t;

  for (t = 0; t < d; t++)
    ext[t] = NN(ths->n[t]) + 1;

  return Y(tile_layout)(d, 2 * ths->m + 2, ext, edge, ntiles);
}

/** Sort key of node i of ths->x: tile, then cell in the tile; nodes outside
 *  [0, 0.5] on some axis, whose window may leave the padded grid, get the key
 *  nt * tilevol of the bucket after the nt tiles. */
static INT node_key(const X(plan) *ths, const INT i, const INT *edge,
    const INT *ntiles, const INT overflow)
{
  const INT d = ths->d;
  INT c[d], t;

  for (t = 0; t < d; t++)
  {
    /* as uo() */
    c[t] = (INT)FLOOR(ths->x[i * d + t] * (2 * NN(ths->n[t])));
    if (c[t] < 0 || c[t] > NN(ths->n[t]))
      return overflow;
  }

  return Y(tile_key)(d, c, edge, ntiles);
}

/**
 * Sorts the nodes by the tile of the grid of B and, within a tile, by the
 * cell they lie in, for better cache utilisation in the multiplication with
 * the matrix B and for the threaded B^T (see Y(tiles_by_colour)). The result
 * is written to ths->index_x as (key, node index) pairs and the node range of
 * each tile, and of the bucket of nodes outside [0, 0.5] after them, to
 * ths->tile_ptr; the nodes remain unchanged.
 */
static void sort0(X(plan) *ths)
{
  const INT d = ths->d, M = ths->M_total;
  INT edge[d], ntiles[d], nt, tilevol = 1, j, t, rhigh;
  INT *ar_x = ths->index_x, *ar_x_temp;

  nt = tile_layout(ths, edge, ntiles);
  for (t = 0; t < d; t++)
    tilevol *= edge[t];

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(j)
#endif
  for (j = 0; j < M; j++)
  {
    ar_x[2 * j] = node_key(ths, j, edge, ntiles, nt * tilevol);
    ar_x[2 * j + 1] = j;
  }

  rhigh = (INT) LRINT(CEIL(LOG2((R)(nt * tilevol + 1)))) - 1;

  ar_x_temp = (INT*) Y(malloc)(2 * (size_t)(M) * sizeof(INT));
  Y(sort_node_indices_radix_lsdf)(M, ar_x, ar_x_temp, rhigh);
  Y(free)(ar_x_temp);

  if (ths->tile_ptr == NULL)
    ths->tile_ptr = (INT*) Y(malloc)((size_t)(nt + 2) * sizeof(INT));
  Y(tile_ptr_from_keys)(M, ar_x, tilevol, nt + 1, ths->tile_ptr);
}

/** Sorts the nodes for one transform, if flag NFFT_SORT_NODES is set and no
 *  precomputation has taken a snapshot of them (see sort_nodes). */
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
    ths->f_sorted = (R*) Y(malloc)((size_t)(M) * sizeof(R));
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
static inline R *node_f(const X(plan) *ths)
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
 * Halo of the grid of the matrix B. The window of a node covers the extended
 * indices e = u - OFFSET, ..., u - OFFSET + 2m+1 (u from uo()), of the even
 * (NFCT) or odd (NFST) periodic continuation of the grid, which run up to m+1
 * (NFST: m+2) points past either end. The grid is padded with HALO_LO points
 * before and HALO_HI points after each dimension (one spare point on either
 * side), so that B never reflects an index. The halo holds the continuation:
 * filled from the grid before B (trafo), added back to it after B^T
 * (adjoint).
 */
#define HALO_LO(m) ((m) + 2)
#define HALO_HI(m) ((m) + 3)
#define HALO(m) (HALO_LO(m) + HALO_HI(m))

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
  return padded_total(ths, ths->halo_lo + ths->halo_hi);
}

/** Offset of grid point 0 in a grid of ths padded with lo points before and
 *  hi points after each dimension. */
static inline INT padded_offset(const X(plan) *ths, const INT lo, const INT hi)
{
  INT t, off = 0;

  for (t = 0; t < ths->d; t++)
    off = off * (ths->n[t] + lo + hi) + lo;

  return off;
}

/**
 * Source of the extended index e of the continuation of a grid of n points:
 * returns the grid index and sets *sign, or returns -1 for a zero (NFST).
 */
static inline INT halo_source(const INT e, const INT n, R *sign)
{
  const INT period = 2 * NN(n);
  const INT q = ((e + OFFSET) % period + period) % period;

#if OFFSET == 0
  *sign = K(1.0);
  return q <= NN(n) ? q : period - q;
#else
  if (q == 0 || q == NN(n))
  {
    *sign = K(0.0);
    return -1;
  }
  if (q < NN(n))
  {
    *sign = K(1.0);
    return q - 1;
  }
  *sign = K(-1.0);
  return period - q - 1;
#endif
}

/**
 * Copies the n_0 x ... x n_{d-1} grid from src, padded with (los, his) points
 * per dimension, to dst, padded with (lod, hid).
 */
static void copy_grid(const X(plan) *ths, R *dst, const INT lod, const INT hid,
    const R *src, const INT los, const INT his)
{
  const INT d = ths->d, nl = ths->n[d - 1];
  R *dst0 = dst + padded_offset(ths, lod, hid);
  const R *src0 = src + padded_offset(ths, los, his);
  INT rows = 1, t, q;

  for (t = 0; t < d - 1; t++)
    rows *= ths->n[t];

#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(q, t)
#endif
  for (q = 0; q < rows; q++)
  {
    INT r = q, od = 0, os = 0, sd = nl + lod + hid, ss = nl + los + his;

    for (t = d - 2; t >= 0; t--)
    {
      const INT i = r % ths->n[t];
      r /= ths->n[t];
      od += i * sd;
      os += i * ss;
      sd *= ths->n[t] + lod + hid;
      ss *= ths->n[t] + los + his;
    }

    memcpy(dst0 + od, src0 + os, (size_t)(nl) * sizeof(R));
  }
}

/** Length of the pieces of a row of the padded grid that fill_halo and
 *  fold_halo hand to one thread: they work on (outer row, piece) pairs, so
 *  that the first dimension, which has one outer row, runs in parallel too. */
#define HALO_BLOCK 1024

/**
 * The h halo rows of a dimension of nt points, padded with lo points before:
 * dst[q] is the padded index of halo row q, src[q] the grid row it continues
 * (-1 for a zero) and sign[q] its sign.
 */
static inline void halo_rows(const INT nt, const INT lo, const INT h,
    INT *src, INT *dst, R *sign)
{
  INT q;

  for (q = 0; q < h; q++)
  {
    const INT e = q < lo ? q - lo : nt + q - lo;
    dst[q] = e + lo;
    src[q] = halo_source(e, nt, &sign[q]);
  }
}

/**
 * Fills the halo of the padded grid gp with the continuation of the grid,
 * dimension by dimension: each halo row of dimension t is (plus or minus) a
 * row of the grid, or zero, over the whole padded range of the other
 * dimensions, so that corners come out right.
 */
static void fill_halo(const X(plan) *ths, R *gp)
{
  const INT d = ths->d, lo = HALO_LO(ths->m), hi = HALO_HI(ths->m);
  const INT h = lo + hi;
  INT src[h], dst[h];
  R sign[h];
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

    halo_rows(nt, lo, h, src, dst, sign);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(q)
#endif
    for (q = 0; q < outer * nb; q++)
    {
      const INT o = q / nb, i0 = (q - o * nb) * HALO_BLOCK;
      const INT i1 = MIN(i0 + HALO_BLOCK, inner);
      R *row = gp + o * pt * inner;
      INT r, i;

      for (r = 0; r < h; r++)
      {
        R *rd = row + dst[r] * inner;

        if (src[r] < 0)
          for (i = i0; i < i1; i++)
            rd[i] = K(0.0);
        else
        {
          const R *rs = row + (src[r] + lo) * inner;
          const R sg = sign[r];
          for (i = i0; i < i1; i++)
            rd[i] = sg * rs[i];
        }
      }
    }
  }
}

/**
 * Adds the halo of the padded grid gp back to the grid points it continues,
 * dimension by dimension in reverse order (see fill_halo). The NFCT
 * additionally doubles the first and last grid row of each dimension: its
 * matrix B^T weights the inner points by 1/2, which D^T applies to all.
 */
static void fold_halo(const X(plan) *ths, R *gp)
{
  const INT d = ths->d, lo = HALO_LO(ths->m), hi = HALO_HI(ths->m);
  const INT h = lo + hi;
  INT src[h], dst[h];
  R sign[h];
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

    halo_rows(nt, lo, h, src, dst, sign);

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(q)
#endif
    for (q = 0; q < outer * nb; q++)
    {
      const INT o = q / nb, i0 = (q - o * nb) * HALO_BLOCK;
      const INT i1 = MIN(i0 + HALO_BLOCK, inner);
      R *row = gp + o * pt * inner;
      INT r, i;

      for (r = 0; r < h; r++)
        if (src[r] >= 0)
        {
          R *rd = row + (src[r] + lo) * inner;
          const R *rs = row + dst[r] * inner;
          const R sg = sign[r];
          for (i = i0; i < i1; i++)
            rd[i] += sg * rs[i];
        }
#if OFFSET == 0
      {
        R *r0 = row + lo * inner, *r1 = row + (nt - 1 + lo) * inner;
        for (i = i0; i < i1; i++)
        {
          r0[i] *= K(2.0);
          r1[i] *= K(2.0);
        }
      }
#endif
    }
  }
}

/** Points ths->g at the padded grid of the matrix B: the FFTW grid gf (g1 or
 *  g2) when it is padded, else the plan's scratch grid ths->g_pad (allocated
 *  on first use and grown if n has grown). */
static void select_grid_B(X(plan) *ths, R *gf)
{
  if (ths->halo_lo == HALO_LO(ths->m) && ths->halo_hi == HALO_HI(ths->m))
  {
    ths->g = gf;
    return;
  }

  {
    const INT size = grid_B_total(ths);
    if (ths->g_pad == NULL || ths->g_pad_size < size)
    {
      Y(free)(ths->g_pad);
      ths->g_pad = (R*) Y(malloc)((size_t)(size) * sizeof(R));
      ths->g_pad_size = size;
    }
  }
  ths->g = ths->g_pad;
}

/** After the FFT of the trafo: makes ths->g the padded grid of B, holding the
 *  FFT output gf and its filled halo. */
static void grid_B_trafo(X(plan) *ths, R *gf)
{
  select_grid_B(ths, gf);
  if (ths->g != gf)
    copy_grid(ths, ths->g, HALO_LO(ths->m), HALO_HI(ths->m), gf,
        ths->halo_lo, ths->halo_hi);
  fill_halo(ths, ths->g);
}

/** After B^T of the adjoint: folds the halo of ths->g and leaves the result
 *  in the FFT input gf. */
static void grid_B_adjoint_done(X(plan) *ths, R *gf)
{
  fold_halo(ths, ths->g);
  if (ths->g != gf)
    copy_grid(ths, gf, ths->halo_lo, ths->halo_hi, ths->g, HALO_LO(ths->m),
        HALO_HI(ths->m));
}

#define MACRO_D_compute_A \
{ \
  g_hat[kg_plain[ths->d]] = f_hat[k_L] * c_phi_inv_k[ths->d]; \
}

#define MACRO_D_compute_T \
{ \
  f_hat[k_L] = g_hat[kg_plain[ths->d]] * c_phi_inv_k[ths->d]; \
}

#define MACRO_D_init_result_A memset(ths->g_hat, 0, (size_t)(grid_fftw_total(ths)) * sizeof(R));

#define MACRO_D_init_result_T memset(f_hat, 0, (size_t)(ths->N_total) * sizeof(R));

#define MACRO_with_PRE_PHI_HUT ths->c_phi_inv[t][kg[t]]

#define MACRO_compute_PHI_HUT_INV (K(1.0) / (PHI_HUT((2 * NN(ths->n[t])), kg[t] + OFFSET, t)))

#define MACRO_init_k_ks \
{ \
  for (t = 0; t < ths->d; t++) \
  { \
    kg[t] = 0; \
  } \
  i = 0; \
}

#define MACRO_update_c_phi_inv_k(what_kind, which_phi) \
{ \
  for (t = i; t < ths->d; t++) \
  { \
    MACRO_update_c_phi_inv_k_ ## what_kind(which_phi); \
    kg_plain[t+1] = kg_plain[t] * (ths->n[t] + ths->halo_lo + ths->halo_hi) \
      + kg[t]; \
  } \
}

#define MACRO_update_c_phi_inv_k_A(which_phi) \
{ \
  c_phi_inv_k[t+1] = K(0.5) * c_phi_inv_k[t] * MACRO_ ## which_phi; \
}

#define MACRO_update_c_phi_inv_k_T(which_phi) \
{ \
  c_phi_inv_k[t+1] = K(0.5) * c_phi_inv_k[t] * MACRO_ ## which_phi; \
}

#define MACRO_count_k_ks \
{ \
  kg[ths->d - 1]++; \
  i = ths->d - 1; \
\
  while ((kg[i] == ths->N[i] - 1) && (i > 0)) \
  { \
    kg[i - 1]++; \
    kg[i] = 0; \
    i--; \
  } \
}

/* sub routines for the fast transforms  matrix vector multiplication with D, D^T */
#define MACRO_D(which_one) \
static inline void D_ ## which_one (X(plan) *ths) \
{ \
  R *g_hat, *f_hat; /* local copy */ \
  R c_phi_inv_k[ths->d+1]; /* postfix product of PHI_HUT */ \
  INT t; /* index dimensions */ \
  INT i; \
  INT k_L; /* plain index */ \
  INT kg[ths->d]; /* multi index in g_hat */ \
  INT kg_plain[ths->d+1]; /* postfix plain index */ \
\
  f_hat = (R*)ths->f_hat; \
  g_hat = (R*)ths->g_hat + padded_offset(ths, ths->halo_lo, ths->halo_hi); \
  MACRO_D_init_result_ ## which_one; \
\
  c_phi_inv_k[0] = K(1.0); \
  kg_plain[0] = 0; \
\
  MACRO_init_k_ks; \
\
  if (ths->flags & PRE_PHI_HUT) \
  { \
    for (k_L = 0; k_L < ths->N_total; k_L++) \
    { \
      MACRO_update_c_phi_inv_k(which_one, with_PRE_PHI_HUT); \
      MACRO_D_compute_ ## which_one; \
      MACRO_count_k_ks; \
    } \
  } \
  else \
  { \
    for (k_L = 0; k_L < ths->N_total; k_L++) \
    { \
      MACRO_update_c_phi_inv_k(which_one,compute_PHI_HUT_INV); \
      MACRO_D_compute_ ## which_one; \
      MACRO_count_k_ks; \
    } \
  } \
}

MACRO_D(A)
MACRO_D(T)

/* sub routines for the fast transforms matrix vector multiplication with B, B^T */
#define MACRO_B_init_result_A memset(f, 0, (size_t)(ths->M_total) * sizeof(R));
#define MACRO_B_init_result_T Y(zero_reals_parallel)(g, grid_B_total(ths));

/* one node of B and B^T with PRE_FULL_PSI: node k in the node data */
#define MACRO_B_PRE_FULL_PSI_node_A \
{ \
  const R *psik = ths->psi + k * lprod; \
  const INT *ixk = ths->psi_index_g + k * lprod; \
  R a = K(0.0); \
  INT l_L; \
 \
  for (l_L = 0; l_L < lprod; l_L++) \
    a += psik[l_L] * g[ixk[l_L]]; \
  f[k] = a; \
}

#define MACRO_B_PRE_FULL_PSI_node_T \
{ \
  const R *psik = ths->psi + k * lprod; \
  const INT *ixk = ths->psi_index_g + k * lprod; \
  INT l_L; \
 \
  for (l_L = 0; l_L < lprod; l_L++) \
    g[ixk[l_L]] += psik[l_L] * f[k]; \
}

/* l[t]: index in the padded grid of B of the first point of the window of
 * node j on axis t; fits: whether the windows on all axes lie in the padded
 * grid, which holds for nodes in [0, 0.5] */
#define MACRO_init_uo_l_lj_t \
{ \
  fits = 1; \
  for (t2 = 0; t2 < ths->d; t2++) \
  { \
    uo(ths, xs, j, &u[t2], &o[t2], t2); \
    l[t2] = u[t2] - OFFSET + HALO_LO(ths->m); \
    fits = fits && l[t2] >= 0 \
      && l[t2] <= ths->n[t2] + HALO(ths->m) - (2 * ths->m + 2); \
    lj[t2] = 0; \
   } \
   t2 = 0; \
}

/* the halo of the grid of B holds the reflected (NFST: and negated) values,
 * so B is a plain tensor product, for B^T alike (which_one) */
#define MACRO_update_phi_prod_ll_plain(which_one,which_psi) \
{ \
  for (t = t2; t < ths->d; t++) \
  { \
    phi_prod[t+1] = phi_prod[t] * (MACRO_ ## which_psi); \
    ll_plain[t+1]  = ll_plain[t] * (ths->n[t] + HALO(ths->m)) + l[t] + lj[t]; \
  } \
}

#define MACRO_count_uo_l_lj_t \
{ \
  lj[ths->d - 1]++; \
  t2 = ths->d - 1; \
 \
  while ((lj[t2] == (2 * ths->m + 2)) && (t2 > 0)) \
  { \
    lj[t2 - 1]++; \
    lj[t2] = 0; \
    t2--; \
  } \
}

/* for a node whose window leaves the padded grid (a node outside [0, 0.5]):
 * each index is reflected into the grid as the halo would be, with its sign */
#define MACRO_update_phi_prod_ll_reflected(which_psi) \
{ \
  for (t = t2; t < ths->d; t++) \
  { \
    R sg; \
    const INT s = halo_source(u[t] - OFFSET + lj[t], ths->n[t], &sg); \
    phi_prod[t+1] = phi_prod[t] * (MACRO_ ## which_psi) * sg; \
    ll_plain[t+1] = ll_plain[t] * (ths->n[t] + HALO(ths->m)) \
      + HALO_LO(ths->m) + (s < 0 ? 0 : s); \
  } \
}

/* all 2m+2 to the power d points of the window of node j, for a node whose
 * window leaves the padded grid (a node outside [0, 0.5]): each index is
 * reflected into the grid as the halo would be, with its sign. B (trafo):
 * returns the sum; B^T (adjoint, fj != 0 on entry is added): adds into g. */
static R taps_reflected(const X(plan) *ths, const INT *u, const R *const *w,
    R *g, const int adjoint, const R fj)
{
  const INT d = ths->d, L = 2 * ths->m + 2;
  INT lj[d], t;
  R sum = K(0.0);

  for (t = 0; t < d; t++)
    lj[t] = 0;

  for (;;)
  {
    INT ix = 0;
    R prod = K(1.0);

    for (t = 0; t < d; t++)
    {
      R sg;
      const INT s = halo_source(u[t] - OFFSET + lj[t], ths->n[t], &sg);
      prod *= w[t][lj[t]] * sg;
      ix = ix * (ths->n[t] + HALO(ths->m)) + HALO_LO(ths->m)
        + (s < 0 ? 0 : s);
    }

    if (adjoint)
      g[ix] += prod * fj;
    else
      sum += prod * g[ix];

    for (t = d - 1; t >= 0 && ++lj[t] == L; t--)
      lj[t] = 0;
    if (t < 0)
      break;
  }

  return sum;
}

/* B (A) and B^T (T) for node j, whose window lies in the padded grid: the
 * tensor product kernels of infft.h on the rows of the last dimension */
#define MACRO_B_node_A \
{ \
  R acc[L]; \
  INT i; \
 \
  for (i = 0; i < L; i++) \
    acc[i] = K(0.0); \
  Y(B_tensor_trafo)(d, L, 1, np, l, w, w[d - 1], g, acc); \
  for (i = 0; i < L; i++) \
    f[j] += acc[i]; \
}

#define MACRO_B_node_T \
{ \
  R wf[L]; \
  INT i; \
 \
  for (i = 0; i < L; i++) \
    wf[i] = w[d - 1][i] * f[j]; \
  Y(B_tensor_adjoint)(d, L, 1, np, l, w, wf, g); \
}

#define MACRO_B_node_reflected_A f[j] = taps_reflected(ths, u, w, g, 0, K(0.0));
#define MACRO_B_node_reflected_T taps_reflected(ths, u, w, g, 1, f[j]);

/* B is parallel over the nodes, each node writes its own f[j]. B^T with
 * OpenMP and sorted nodes spreads tiles of the grid in colour phases
 * (adjoint_B_tiled), else it runs serially. */
#ifdef _OPENMP
#define MACRO_B_parallel_A _Pragma("omp parallel for default(shared) private(k)")
static int adjoint_B_tiled(X(plan) *ths);
#define MACRO_B_tiled_T \
  if (adjoint_B_tiled(ths)) \
    return;
#else
#define MACRO_B_parallel_A
#define MACRO_B_tiled_T
#endif
#define MACRO_B_parallel_T
#define MACRO_B_tiled_A

#define MACRO_B(which_one) \
/* B (B^T) for node k of the node data, in a function of its own so that the \
 * OpenMP region of B_ holds no run-time sized arrays (see nfft_trafo_1d_B) */ \
static void B_node_ ## which_one (const X(plan) *ths, const R *xs, R *f, \
    R *g, const INT *np, const R *fg_e, const R *fg_q, const INT k) \
{ \
  const INT d = ths->d, m = ths->m, L = 2 * m + 2; \
  INT t; \
  INT u[d], o[d]; /* multi band with respect to x_j */ \
  INT l[d]; /* first index of the window in the padded grid of B */ \
  INT lj[d]; /* set, but not needed, by MACRO_init_uo_l_lj_t */ \
  int fits; /* window of node j inside the padded grid */ \
  INT t2; \
  /* window weights per axis; 2m+4, not the 2m+3 FG_RUN fills: an even \
   * row stride */ \
  R wbuf[d][2 * m + 4]; \
  const R *w[d]; /* w[t][lj], lj < 2m+2: weights of node j on axis t */ \
  const INT j = NODE_J(ths, k); \
\
  MACRO_init_uo_l_lj_t; \
  UNUSED(lj); \
\
  for (t = 0; t < d; t++) \
  { \
    if (ths->flags & PRE_PSI) \
      w[t] = ths->psi + (j * d + t) * L; \
    else if (ths->flags & PRE_FG_PSI) \
    { \
      /* uo() anchors the run at the grid point below the node, run \
       * index m, which FG_RUN puts at buffer index m + 2 for a fill \
       * based at wbuf[t] + 1 */ \
      FG_RUN(wbuf[t] + 1, m, ths->psi[2 * (j * d + t)], \
          ths->psi[2 * (j * d + t) + 1], fg_e[t], fg_q[t]); \
      w[t] = wbuf[t] + 2; \
    } \
    else if (ths->flags & FG_PSI) \
    { \
      const INT fg_c = u[t] + m; \
      FG_RUN(wbuf[t] + 1, m, \
          (PHI((2 * NN(ths->n[t])), \
               NX_SUB((2 * NN(ths->n[t])), xs[j * d + t], fg_c) \
               / ((R)(2 * NN(ths->n[t]))), (t))), \
          EXP(K(2.0) \
              * NX_SUB((2 * NN(ths->n[t])), xs[j * d + t], fg_c) \
              / ths->b[t]), \
          fg_e[t], fg_q[t]); \
      w[t] = wbuf[t] + 2; \
    } \
    else if (ths->flags & PRE_LIN_PSI) \
    { \
      const INT ip_s = ths->K / (m + 2); \
      const R y = ((NX_SUB((2 * NN(ths->n[t])), xs[j * d + t], u[t])) \
              * ((R)ths->K)) / (m + 2); \
      const INT ip_u = LRINT(FLOOR(y)); \
      const R ip_w = y - ip_u; \
      INT lj_fg; \
\
      for (lj_fg = 0; lj_fg < L; lj_fg++) \
        wbuf[t][lj_fg] = ths->psi[(ths->K + 1) * t \
            + ABS(ip_u - lj_fg * ip_s)] * (1 - ip_w) \
          + ths->psi[(ths->K + 1) * t + ABS(ip_u - lj_fg * ip_s + 1)] \
            * (ip_w); \
      w[t] = wbuf[t]; \
    } \
    else \
    { \
      PHI_RUN(wbuf[t], (2 * NN(ths->n[t])), xs[j * d + t], u[t], t); \
      w[t] = wbuf[t]; \
    } \
  } \
\
  if (fits) \
    MACRO_B_node_ ## which_one \
  else \
    MACRO_B_node_reflected_ ## which_one \
} \
\
static inline void B_ ## which_one (X(plan) *ths) \
{ \
  const INT d = ths->d, m = ths->m, L = 2 * m + 2; \
  const R *xs = node_x(ths); \
  R *f = node_f(ths), *g = (R*)ths->g; \
  INT np[d]; /* sizes of the padded grid of B */ \
  R fg_e[d], fg_q[d]; \
  INT k, t; \
\
  MACRO_B_init_result_ ## which_one \
  MACRO_B_tiled_ ## which_one \
\
  if (ths->flags & PRE_FULL_PSI) \
  { \
    INT lprod = 1; \
\
    for (t = 0; t < d; t++) \
      lprod *= L; \
\
    MACRO_B_parallel_ ## which_one \
    for (k = 0; k < ths->M_total; k++) \
      MACRO_B_PRE_FULL_PSI_node_ ## which_one; \
    return; \
  } \
\
  sort(ths); \
\
  for (t = 0; t < d; t++) \
  { \
    np[t] = ths->n[t] + HALO(m); \
    if (ths->flags & (PRE_FG_PSI | FG_PSI)) \
    { \
      fg_e[t] = EXP(K(-1.0) / ths->b[t]); \
      fg_q[t] = fg_e[t] * fg_e[t]; \
    } \
  } \
\
  MACRO_B_parallel_ ## which_one \
  for (k = 0; k < ths->M_total; k++) \
    B_node_ ## which_one (ths, xs, f, g, np, fg_e, fg_q, k); \
} /* B */

MACRO_B(A)
MACRO_B(T)

#ifdef _OPENMP
/* What the nodes of one tile of the threaded B^T need (adjoint_B_tile). */
typedef struct
{
  X(plan) *ths;
  const R *xs; /* node data, as node_x(ths) */
  R *f; /* samples, as node_f(ths) */
  const INT *np; /* sizes of the padded grid of B */
  const R *fg_e, *fg_q; /* factors of the FG modes per axis */
  INT lprod; /* (2m+2)^d, the PRE_FULL_PSI entries per node */
} adjoint_tile_ctx;

/** B^T for the sorted nodes k0, ..., k1-1 of one tile (or of the bucket of
 *  nodes outside [0, 0.5]), straight into the padded grid of B. */
static void adjoint_B_tile(void *p, const INT k0, const INT k1)
{
  const adjoint_tile_ctx *c = (const adjoint_tile_ctx*)p;
  const X(plan) *ths = c->ths;
  R *g = (R*)ths->g, *f = c->f;
  INT k;

  if (ths->flags & PRE_FULL_PSI)
  {
    const INT lprod = c->lprod;

    for (k = k0; k < k1; k++)
      MACRO_B_PRE_FULL_PSI_node_T;
    return;
  }

  for (k = k0; k < k1; k++)
    B_node_T(ths, c->xs, c->f, g, c->np, c->fg_e, c->fg_q, k);
}

/**
 * Threaded B^T for sorted nodes (flag NFFT_SORT_NODES), any d and psi mode:
 * the tiles of the grid are spread in colour phases, the threads of a phase
 * taking whole tiles, without atomics (see Y(tiles_by_colour)); then the
 * nodes outside [0, 0.5], whose windows are reflected anywhere into the
 * grid, one by one. The result does not depend on the number of threads.
 * ths->g must be zero. Returns 0 and does nothing without sorted nodes.
 */
static int adjoint_B_tiled(X(plan) *ths)
{
  const INT d = ths->d, m = ths->m;
  INT np[d], edge[d], ntiles[d], nt, t;
  R fg_e[d], fg_q[d];
  adjoint_tile_ctx c;

  if (!(ths->flags & NFFT_SORT_NODES))
    return 0;

  sort(ths);
  nt = tile_layout(ths, edge, ntiles);

  for (t = 0; t < d; t++)
  {
    np[t] = ths->n[t] + HALO(m);
    if (ths->flags & (PRE_FG_PSI | FG_PSI))
    {
      fg_e[t] = EXP(K(-1.0) / ths->b[t]);
      fg_q[t] = fg_e[t] * fg_e[t];
    }
  }

  c.ths = ths;
  c.xs = node_x(ths);
  c.f = node_f(ths);
  c.np = np;
  c.fg_e = fg_e;
  c.fg_q = fg_q;
  for (t = 0, c.lprod = 1; t < d; t++)
    c.lprod *= 2 * m + 2;

  Y(tiles_by_colour)(d, ntiles, ths->tile_ptr, adjoint_B_tile, &c);
  adjoint_B_tile(&c, ths->tile_ptr[nt], ths->tile_ptr[nt + 1]);
  return 1;
}
#endif

/**
 * user routines
 */
void X(trafo)(X(plan) *ths)
{
  switch(ths->d)
  {
    default:
    {
      /* use ths->my_fftw_r2r_plan */
      ths->g_hat = ths->g1;
      ths->g = ths->g2;

      /* form \f$ \hat g_k = \frac{\hat f_k}{c_k\left(\phi\right)} \text{ for }
       * k \in I_N \f$ */
      TIC(0)
      D_A(ths);
      TOC(0)

      /* Compute by d-variate discrete Fourier transform
       * \f$ g_l = \sum_{k \in I_N} \hat g_k {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
       * \text{ for } l \in I_n \f$ */
      TIC_FFTW(1)
      FFTW(execute)(ths->my_fftw_r2r_plan);
      TOC_FFTW(1)
      grid_B_trafo(ths, ths->g2);

      /* Set \f$ f_j = \sum_{l \in I_n,m(x_j)} g_l \psi\left(x_j-\frac{l}{n}\right)
       * \text{ for } j=0,\dots,M-1 \f$ */
      TIC(2)
      B_A(ths);
      scatter_f(ths);
      TOC(2)
    }
  }
} /* trafo */

void X(adjoint)(X(plan) *ths)
{
  switch(ths->d)
  {
    default:
    {
      /* use ths->my_fftw_plan */
      ths->g_hat = ths->g2;
      ths->g = ths->g1;

      /* Set \f$ g_l = \sum_{j=0}^{M-1} f_j \psi\left(x_j-\frac{l}{n}\right)
       * \text{ for } l \in I_n,m(x_j) \f$ */
      TIC(2)
      select_grid_B(ths, ths->g1);
      gather_f(ths);
      B_T(ths);
      grid_B_adjoint_done(ths, ths->g1);
      TOC(2)

      /* Compute by d-variate discrete cosine transform
       * \f$ \hat g_k = \sum_{l \in I_n} g_l {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
       * \text{ for }  k \in I_N\f$ */
      TIC_FFTW(1)
      FFTW(execute)(ths->my_fftw_r2r_plan);
      TOC_FFTW(1)

      /* Form \f$ \hat f_k = \frac{\hat g_k}{c_k\left(\phi\right)} \text{ for }
       * k \in I_N \f$ */
      TIC(0)
      D_T(ths);
      TOC(0)
    }
  }
} /* adjoint */

/** initialisation of direct transform
 */
static inline void precompute_phi_hut(X(plan) *ths)
{
  INT ks[ths->d]; /* index over all frequencies */
  INT t; /* index over all dimensions */

  ths->c_phi_inv = (R**) Y(malloc)((size_t)(ths->d) * sizeof(R*));

  for (t = 0; t < ths->d; t++)
  {
    ths->c_phi_inv[t] = (R*)Y(malloc)((size_t)(ths->N[t] - OFFSET) * sizeof(R));

    for (ks[t] = 0; ks[t] < ths->N[t] - OFFSET; ks[t]++)
    {
      ths->c_phi_inv[t][ks[t]] = (K(1.0) / (PHI_HUT((2 * NN(ths->n[t])), ks[t] + OFFSET, t)));
    }
  }
} /* phi_hut */

/** create a lookup table, but NOT for each node
 *  good idea K=2^xx
 *  TODO: estimate K, call from init
 *  assumes an EVEN window function
 */
void X(precompute_lin_psi)(X(plan) *ths)
{
  INT t; /**< index over all dimensions */
  INT j; /**< index over all nodes */
  R step; /**< step size in [0,(m+2)/n] */

  for (t = 0; t < ths->d; t++)
  {
    step = ((R)(ths->m+2)) / (((R)ths->K) * (2 * NN(ths->n[t])));

    for (j = 0; j <= ths->K; j++)
    {
      ths->psi[(ths->K + 1) * t + j] = PHI((2 * NN(ths->n[t])), (j * step), t);
    } /* for(j) */
  } /* for(t) */
}

static void precompute_fg_psi_sorted(X(plan) *ths)
{
  INT t; /* index over all dimensions */
  INT u, o, c; /* depends on x_j */
  const R *xs;

  xs = node_x(ths);

  for (t = 0; t < ths->d; t++)
  {
    INT j;
//    #pragma omp parallel for default(shared) private(j,u,o)
    for (j = 0; j < ths->M_total; j++)
    {
      uo(ths, xs, j, &u, &o, t);

      /* anchored at the grid point below the node, which uo() puts at run
       * index m */
      c = u + ths->m;

      ths->psi[2 * (j*ths->d + t)] = (PHI((2 * NN(ths->n[t])),
          NX_SUB((2 * NN(ths->n[t])), xs[j * ths->d + t], c)
          / ((R)(2 * NN(ths->n[t]))), (t)));
      ths->psi[2 * (j*ths->d + t) + 1] = EXP(K(2.0)
          * NX_SUB((2 * NN(ths->n[t])), xs[j * ths->d + t], c) / ths->b[t]);
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

  for (t = 0; t < ths->d; t++)
  {
    INT j;

    for (j = 0; j < ths->M_total; j++)
    {
      uo(ths, xs, j, &u, &o, t);

      PHI_RUN(ths->psi + (j * ths->d + t) * (2 * ths->m + 2),
          (2 * NN(ths->n[t])), xs[(j) * ths->d + (t)], u, t);
    } /* for (j) */
  } /* for (t) */
} /* precompute_psi */

void X(precompute_psi)(X(plan) *ths)
{
  sort_nodes(ths);
  precompute_psi_sorted(ths);
}

static void precompute_full_psi_sorted(X(plan) *ths)
{
  INT t, t2; /* index over all dimensions */
  INT j; /* index over all nodes */
  INT l_L; /* plain index 0 <= l_L < lprod */
  INT l[ths->d]; /* first index of the window in the padded grid of B */
  INT lj[ths->d]; /* multi index 0<=lj<u+o+1 */
  int fits; /* window of node j inside the padded grid */
  INT ll_plain[ths->d+1]; /* postfix plain index */
  INT lprod; /* 'bandwidth' of matrix B */
  INT u[ths->d], o[ths->d]; /* depends on x_j */

  R phi_prod[ths->d+1];

  INT ix, ix_old;
  const R *xs;

  xs = node_x(ths);

  phi_prod[0] = K(1.0);
  ll_plain[0]  = 0;

  for (t = 0, lprod = 1; t < ths->d; t++)
    lprod *= 2 * ths->m + 2;

  for (j = 0, ix = 0, ix_old = 0; j < ths->M_total; j++)
  {
    R psij_const[ths->d * (2 * ths->m + 2)];

    MACRO_init_uo_l_lj_t;
    MACRO_fill_psij_const;

    for (l_L = 0; l_L < lprod; l_L++, ix++)
    {
      if (fits)
        MACRO_update_phi_prod_ll_plain(A, without_PRE_PSI)
      else
        MACRO_update_phi_prod_ll_reflected(without_PRE_PSI)

      ths->psi_index_g[ix] = ll_plain[ths->d];
      ths->psi[ix] = phi_prod[ths->d];

      MACRO_count_uo_l_lj_t;
    } /* for (l_L) */

    ths->psi_index_f[j] = ix - ix_old;
    ix_old = ix;
  } /* for(j) */
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

  if(ths->flags & PRE_PSI)
    precompute_psi_sorted(ths);
  if(ths->flags & PRE_FULL_PSI)
    precompute_full_psi_sorted(ths);
  if(ths->flags & PRE_FG_PSI)
    precompute_fg_psi_sorted(ths);
  if(ths->flags & PRE_LIN_PSI)
    X(precompute_lin_psi)(ths);
}

static inline void init_help(X(plan) *ths)
{
  INT t; /* index over all dimensions */
  INT lprod; /* 'bandwidth' of matrix B */

  if (ths->flags & NFFT_OMP_BLOCKWISE_ADJOINT)
    ths->flags |= NFFT_SORT_NODES;

  ths->N_total = intprod(ths->N, OFFSET, ths->d);
  ths->n_total = intprod(ths->n, 0, ths->d);

  ths->sigma = (R*)Y(malloc)((size_t)(ths->d) * sizeof(R));

  for (t = 0; t < ths->d; t++)
    ths->sigma[t] = ((R)NN(ths->n[t])) / ths->N[t];

  /* Assign r2r transform kinds for each dimension */
  ths->r2r_kind = (FFTW(r2r_kind)*)Y(malloc)((size_t)(ths->d) * sizeof (FFTW(r2r_kind)));
  for (t = 0; t < ths->d; t++)
    ths->r2r_kind[t] = FOURIER_TRAFO;

  WINDOW_HELP_INIT;

  if (ths->flags & MALLOC_X)
    ths->x = (R*)Y(malloc)((size_t)(ths->d * ths->M_total) * sizeof(R));

  if (ths->flags & MALLOC_F_HAT)
    ths->f_hat = (R*)Y(malloc)((size_t)(ths->N_total) * sizeof(R));

  if (ths->flags & MALLOC_F)
    ths->f = (R*)Y(malloc)((size_t)(ths->M_total) * sizeof(R));

  if (ths->flags & PRE_PHI_HUT)
    precompute_phi_hut(ths);

  if(ths->flags & PRE_LIN_PSI)
  {
      ths->K = (1U<< 10) * (ths->m+2);
      ths->psi = (R*) Y(malloc)((size_t)((ths->K + 1) * ths->d) * sizeof(R));
  }

  if(ths->flags & PRE_FG_PSI)
    ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * ths->d * 2) * sizeof(R));

  if (ths->flags & PRE_PSI)
    ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * ths->d * (2 * ths->m + 2 )) * sizeof(R));

  if(ths->flags & PRE_FULL_PSI)
  {
      for (t = 0, lprod = 1; t < ths->d; t++)
        lprod *= 2 * ths->m + 2;

      ths->psi = (R*) Y(malloc)((size_t)(ths->M_total * lprod) * sizeof(R));

      ths->psi_index_f = (INT*) Y(malloc)((size_t)(ths->M_total) * sizeof(INT));
      ths->psi_index_g = (INT*) Y(malloc)((size_t)(ths->M_total * lprod) * sizeof(INT));
  }

  ths->halo_lo = 0;
  ths->halo_hi = 0;

  if (ths->flags & FFTW_INIT)
  {
    /* g1 and g2 are the padded grid of B, the FFT works on their
     * n_0 x ... x n_{d-1} part */
    ths->halo_lo = HALO_LO(ths->m);
    ths->halo_hi = HALO_HI(ths->m);

    ths->g1 = (R*)Y(malloc)((size_t)(grid_fftw_total(ths)) * sizeof(R));

    if (ths->flags & FFT_OUT_OF_PLACE)
      ths->g2 = (R*) Y(malloc)((size_t)(grid_fftw_total(ths)) * sizeof(R));
    else
      ths->g2 = ths->g1;

    {
      const INT off = padded_offset(ths, ths->halo_lo, ths->halo_hi);
      FFTW(iodim64) *dims = Y(malloc)((size_t)(ths->d) * sizeof(FFTW(iodim64)));
      ptrdiff_t stride = 1;

      for (t = ths->d - 1; t >= 0; t--)
      {
        dims[t].n = (ptrdiff_t)(ths->n[t]);
        dims[t].is = stride;
        dims[t].os = stride;
        stride *= (ptrdiff_t)(ths->n[t] + ths->halo_lo + ths->halo_hi);
      }

      ths->my_fftw_r2r_plan = FFTW(plan_guru64_r2r)((int)ths->d, dims, 0,
          NULL, ths->g1 + off, ths->g2 + off, ths->r2r_kind, ths->fftw_flags);
      Y(free)(dims);
    }
  }

  if (ths->flags & NFFT_SORT_NODES)
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
  int t; /* index over all dimensions */

  ths->d = (INT)d;

  ths->N = (INT*) Y(malloc)((size_t)(d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->N[t] = (INT)N[t];

  ths->M_total = (INT)M_total;

  ths->n = (INT*) Y(malloc)((size_t)(d) * sizeof(INT));

  for (t = 0; t < d; t++)
    ths->n[t] = 2 * (Y(next_power_of_2)(ths->N[t]) - 1) + OFFSET;

  ths->m = WINDOW_HELP_ESTIMATE_m;

  if (d > 1)
  {
//#ifdef _OPENMP
//    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
//                      FFTW_INIT | NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT;
//#else
    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | NFFT_SORT_NODES;
//#endif
  }
  else
    /* sorted, as for d > 1: the transforms gain more than the sort costs */
    ths->flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X| MALLOC_F_HAT | MALLOC_F |
                      FFTW_INIT | FFT_OUT_OF_PLACE | NFFT_SORT_NODES;

  ths->fftw_flags = FFTW_ESTIMATE | FFTW_DESTROY_INPUT;

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

  for (j = 0; j < ths->M_total * ths->d; j++)
  {
    if ((ths->x[j] < K(0.0)) || (ths->x[j] >= K(0.5)))
    {
      return "ths->x out of range [0.0,0.5)";
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

    if(ths->N[j] - 1 <= ths->m)
      return "Polynomial degree N is smaller than cut-off m";
  }
  return 0;
}

void X(finalize)(X(plan) *ths)
{
  INT t; /* index over dimensions */

  if (ths->flags & NFFT_SORT_NODES)
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

  if (ths->flags & FFTW_INIT)
  {
#ifdef _OPENMP
    #pragma omp critical (nfft_omp_critical_fftw_plan)
#endif
    FFTW(destroy_plan)(ths->my_fftw_r2r_plan);

    if (ths->flags & FFT_OUT_OF_PLACE)
      Y(free)(ths->g2);

    Y(free)(ths->g1);
  }

  if(ths->flags & PRE_FULL_PSI)
  {
    Y(free)(ths->psi_index_g);
    Y(free)(ths->psi_index_f);
    Y(free)(ths->psi);
  }

  if (ths->flags & PRE_PSI)
    Y(free)(ths->psi);

  if(ths->flags & PRE_FG_PSI)
    Y(free)(ths->psi);

  if(ths->flags & PRE_LIN_PSI)
    Y(free)(ths->psi);

  if (ths->flags & PRE_PHI_HUT)
  {
    for (t = 0; t < ths->d; t++)
      Y(free)(ths->c_phi_inv[t]);
    Y(free)(ths->c_phi_inv);
  }

  if (ths->flags & MALLOC_F)
    Y(free)(ths->f);

  if(ths->flags & MALLOC_F_HAT)
    Y(free)(ths->f_hat);

  if (ths->flags & MALLOC_X)
    Y(free)(ths->x);

  WINDOW_HELP_FINALIZE;

  Y(free)(ths->N);
  Y(free)(ths->n);
  Y(free)(ths->sigma);

  Y(free)(ths->r2r_kind);
} /* finalize */
