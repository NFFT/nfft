/*
 * Copyright (c) 2026 Jens Keiter, Stefan Kunis, Daniel Potts
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

/* NFFT_PROBLEM_CONV -- Step C (node convolution, matrix B) of the NFFT
 * decomposition. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "conv.h"

#include <math.h> /* pow, HUGE_VAL */

static void hash(const problem *p, md5 *ctx)
{
  const problem_conv *ego = (const problem_conv *)p;
  int t;
  Y(md5_put_str)(ctx, "conv");
  Y(md5_put_int)(ctx, ego->sign);
  Y(tensor_md5)(ctx, ego->sz);
  Y(tensor_md5)(ctx, ego->vecsz);
  Y(md5_put_int)(ctx, (int)Y(log2i)(ego->M));
  Y(md5_put_int)(ctx, ego->m);
  for (t = 0; t < ego->sz->rnk; t++)
    Y(md5_put_INT)(ctx, ego->N[t]);
  Y(md5_put_int)(ctx, ego->window);
  /* x itself is not hashed: the key is a size class. Its one recorded
   * property is whether the rows are already in key order, since a sorted
   * winner blessed on random nodes must not replay on presorted ones, nor the
   * reverse. */
  Y(md5_put_int)(ctx, ego->x_sorted);
}

static void print(const problem *p, printer *pr)
{
  const problem_conv *ego = (const problem_conv *)p;
  pr->print(pr, "(conv sign=%d m=%d M=%D ", ego->sign, ego->m, ego->M);
  Y(tensor_print)(ego->sz, pr);
  pr->putchr(pr, ')');
}

static void destroy(problem *p)
{
  problem_conv *ego = (problem_conv *)p;
  Y(tensor_destroy)(ego->sz);
  Y(tensor_destroy)(ego->vecsz);
  Y(free)(ego->N);
  /* x is a borrowed caller array */
}

static const problem_adt conv_adt = {NFFT_PROBLEM_CONV, hash, print, destroy};

problem *Y(mkproblem_conv)(int d, const INT *n, const INT *N, INT M, int m,
                        int window, int sign, int x_sorted, R *x, C *g, C *f)
{
  problem_conv *ego;
  tensor *grid;
  int t;

  A(sign == 1 || sign == -1);
  A(d >= 1);
  A(M >= (INT)1);

  /* square grid geometry n_t -> n_t (row-major strides). */
  grid = Y(tensor_create)(d);
  {
    INT s = (INT)1;
    for (t = d - 1; t >= 0; --t) {
      A(n[t] >= (INT)1);
      grid->dims[t] = Y(mvdim_square)(n[t], s, s);
      s *= n[t];
    }
  }

  ego = (problem_conv *)Y(problem_create)(sizeof(problem_conv), &conv_adt);
  ego->sz = grid;
  ego->vecsz = Y(tensor_create)(0);
  ego->N = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  for (t = 0; t < d; t++) {
    A(N[t] >= (INT)1);
    ego->N[t] = N[t];
  }
  ego->M = M;
  ego->m = m;
  ego->window = window;
  ego->sign = sign;
  ego->x_sorted = x_sorted;
  ego->x = x; /* borrowed alias */
  ego->g = g; /* borrowed alias */
  ego->f = f; /* borrowed alias */
  return (problem *)ego;
}

INT Y(problem_conv_n)(const problem *p, int t)
{
  const problem_conv *ego = (const problem_conv *)p;
  A(t >= 0 && t < ego->sz->rnk);
  return ego->sz->dims[t].n_in; /* square: n_in == n_out == n_t */
}

INT Y(problem_conv_N)(const problem *p, int t)
{
  const problem_conv *ego = (const problem_conv *)p;
  A(t >= 0 && t < ego->sz->rnk);
  return ego->N[t];
}

INT Y(problem_conv_ntot)(const problem *p)
{
  const problem_conv *ego = (const problem_conv *)p;
  return Y(tensor_sz_in)(ego->sz);
}

/* One node's tap block, rank >= 4. Axis 0 runs outside the odometer so that it
 * can be clipped to [lo, hi) and its psi factor hoisted; axes 1..d-1 keep the
 * postfix odometer the nd solver has always used. */
C Y(conv_node_nd)(C *g, const R *psi, const INT *u, const INT *n, int d, INT len,
               INT lo, INT hi, C fj, int forward)
{
  INT tof0[2], gof0[2], rl0[2];
  INT lj[d]; /* tap multi index over axes 1..d-1 */
  INT ll[d + 1]; /* postfix plain index within a plane */
  R pp[d + 1]; /* postfix product of psi over axes 1..d-1 */
  INT l_all[d * len]; /* wrapped grid indices per axis/tap */
  INT stride0 = 1, lprod = 1, i, t, t2, l_L;
  C acc = K(0.0);
  int a;

  A((size_t)d * (size_t)len * sizeof(INT) <= (size_t)(64 * 1024));

  for (t = 1; t < d; t++) {
    stride0 *= n[t];
    lprod *= len;
    for (i = 0; i < len; i++)
      l_all[t * len + i] = (u[t] + i) % n[t];
  }

  Y(conv_runs)(u[0], n[0], len, tof0, gof0, rl0);
  for (a = 0; a < 2; a++) {
    Y(conv_clip)(lo, hi, &tof0[a], &gof0[a], &rl0[a]);
    for (i = 0; i < rl0[a]; i++) {
      const R p0 = psi[tof0[a] + i];
      C *gp = g + (gof0[a] + i) * stride0;
      const C fp = fj * p0;
      C sub = K(0.0);
      for (t = 1; t < d; t++)
        lj[t] = 0;
      t = 1;
      pp[1] = K(1.0);
      ll[1] = 0;
      for (l_L = 0; l_L < lprod; l_L++) {
        for (t2 = t; t2 < d; t2++) {
          pp[t2 + 1] = pp[t2] * psi[t2 * len + lj[t2]];
          ll[t2 + 1] = ll[t2] * n[t2] + l_all[t2 * len + lj[t2]];
        }
        if (forward)
          sub += pp[d] * gp[ll[d]];
        else
          gp[ll[d]] += pp[d] * fp;
        for (t = d - 1; (t > 1) && (lj[t] == len - 1); t--)
          lj[t] = 0;
        lj[t]++;
      }
      if (forward)
        acc += sub * p0;
    }
  }
  return acc;
}

/* The wrapped window start u would show a sorted caller order as a rotation;
 * the shifted cell does not. c matches fill_u's expression, so key and u
 * agree on every node. */
INT Y(conv_node_key)(const R *xj, const INT *n, int d)
{
  INT key = 0;
  int t;
  for (t = 0; t < d; t++) {
    const INT c = LRINT(FLOOR(xj[t] * (R)n[t]));
    key = key * n[t] + ((c + n[t] / 2) % n[t] + n[t]) % n[t];
  }
  return key;
}

int Y(conv_nodes_sorted)(const R *x, const INT *n, int d, INT M)
{
  INT j, prev = 0;
  for (j = 0; j < M; j++) {
    const INT key = Y(conv_node_key)(x + j * d, n, d);
    if (key < prev)
      return 0;
    prev = key;
  }
  return 1;
}

/* Stable LSD radix sort of (key, node) pairs. The bit count is an integer
 * count of ntot - 1; legacy sort0 takes it from a floating LOG2. */
void Y(conv_sort_nodes)(const R *x, const INT *n, int d, INT M, INT *perm)
{
  INT *keys = (INT *)Y(malloc)(4 * (size_t)M * sizeof(INT));
  INT ntot = 1, j, rhigh = -1;
  int t;
  for (t = 0; t < d; t++)
    ntot *= n[t];
  for (j = ntot - 1; j > 0; j >>= 1)
    rhigh++;
  for (j = 0; j < M; j++) {
    keys[2 * j] = Y(conv_node_key)(x + j * d, n, d);
    keys[2 * j + 1] = j;
  }
  Y(sort_node_indices_radix_lsdf)(M, keys, keys + 2 * M, rhigh);
  for (j = 0; j < M; j++)
    perm[j] = keys[2 * j + 1];
  Y(free)(keys);
}

void Y(conv_fill_u)(const R *x, const INT *perm, const INT *n, int d, int m,
                 INT M, INT *u)
{
  INT k;
  int t;
  for (t = 0; t < d; t++) {
    const INT nt = n[t];
    for (k = 0; k < M; k++) {
      const INT j = perm ? perm[k] : k;
      INT c = LRINT(FLOOR(x[j * d + t] * (R)nt));
      u[k * d + t] = (((c - m) % nt) + nt) % nt;
    }
  }
}

/* In key order the nodes are gathered a block at a time, so the window
 * evaluation still runs over a contiguous array and no copy of x is kept. */
#define CONV_PSI_BLOCK 256

void Y(conv_fill_psi)(int window, const R *x, const INT *perm, const INT *n,
                   const INT *N, int d, int m, INT M, R *psi)
{
  const INT len = 2 * (INT)m + 2;
  R xb[CONV_PSI_BLOCK];
  INT k0, k, cnt;
  int t;
  if (!perm) {
    for (t = 0; t < d; t++)
      Y(window_phi_precompute)(window, n[t], N[t], m, x + t, d, M, psi + t * len,
                            d * len);
    return;
  }
  for (k0 = 0; k0 < M; k0 += cnt) {
    cnt = (M - k0 < CONV_PSI_BLOCK) ? M - k0 : CONV_PSI_BLOCK;
    for (t = 0; t < d; t++) {
      for (k = 0; k < cnt; k++)
        xb[k] = x[perm[k0 + k] * d + t];
      Y(window_phi_precompute)(window, n[t], N[t], m, xb, 1, cnt,
                            psi + (k0 * d + t) * len, d * len);
    }
  }
}

double Y(conv_sort_grid_bytes)(int d)
{
  static const double bytes[3] = {16777216.0, 131072.0, 8388608.0};
  A(d >= 1);
  return bytes[(d > 3 ? 3 : d) - 1];
}

/* Sorted: one load of perm and one access of f per node. Caller order on
 * random nodes: every node fetches its tap block cold once the grid outgrows
 * the cache, counted as one extra memory operation per tap. That under-counts
 * the measured gain, but ESTIMATE only needs the two candidates ranked. */
void Y(conv_order_ops)(opcnt *ops, INT M, int d, int m, INT ntot, int sorted,
                    int x_sorted)
{
  if (sorted)
    ops->other += 2.0 * (double)M;
  else if (!x_sorted
           && (double)ntot * (double)sizeof(C) > Y(conv_sort_grid_bytes)(d))
    ops->other += (double)M * pow((double)(2 * m + 2), (double)d);
}
