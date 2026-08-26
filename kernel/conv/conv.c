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

static void cv_hash(const problem *p, md5 *ctx) {
  const problem_conv *ego = (const problem_conv *)p;
  int t;
  Y(md5_put_str)
  (ctx, "conv");
  Y(md5_put_int)
  (ctx, ego->sign);
  Y(tensor_md5)
  (ctx, ego->sz);
  Y(tensor_md5)
  (ctx, ego->vecsz);
  Y(md5_put_int)
  (ctx, (int)Y(log2i)(ego->M));
  Y(md5_put_int)
  (ctx, ego->m);
  for (t = 0; t < ego->sz->rnk; t++)
    Y(md5_put_int)
  (ctx, ego->N[t]);
  Y(md5_put_int)
  (ctx, ego->window);
  /* x is deliberately not hashed. */
}

static void cv_print(const problem *p, printer *pr) {
  const problem_conv *ego = (const problem_conv *)p;
  pr->print(pr, "(conv sign=%d m=%d M=%D ", ego->sign, ego->m, ego->M);
  Y(tensor_print)
  (ego->sz, pr);
  pr->putchr(pr, ')');
}

static void cv_destroy(problem *p) {
  problem_conv *ego = (problem_conv *)p;
  Y(tensor_destroy)
  (ego->sz);
  Y(tensor_destroy)
  (ego->vecsz);
  Y(free)
  (ego->N);
  /* x is a borrowed caller array */
}

static const problem_adt conv_adt = {
    NFFT_PROBLEM_CONV, cv_hash, cv_print, cv_destroy};

problem *Y(mkproblem_conv)(int d, const INT *n, const INT *N, INT M, int m,
                           int window, int sign, R *x, C *g, C *f) {
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
  ego->x = x; /* borrowed alias */
  ego->g = g; /* borrowed alias */
  ego->f = f; /* borrowed alias */
  return (problem *)ego;
}

INT Y(problem_conv_n)(const problem *p, int t) {
  const problem_conv *ego = (const problem_conv *)p;
  A(t >= 0 && t < ego->sz->rnk);
  return ego->sz->dims[t].n_in; /* square: n_in == n_out == n_t */
}

INT Y(problem_conv_N)(const problem *p, int t) {
  const problem_conv *ego = (const problem_conv *)p;
  A(t >= 0 && t < ego->sz->rnk);
  return ego->N[t];
}

INT Y(problem_conv_ntot)(const problem *p) {
  const problem_conv *ego = (const problem_conv *)p;
  return Y(tensor_sz_in)(ego->sz);
}
