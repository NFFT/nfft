/*
 * Copyright (c) 2026 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* NFFT_PROBLEM_DECONV -- Step A (deconvolution) of the NFFT decomposition. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

static void dv_hash(const problem *p, md5 *ctx) {
  const problem_deconv *ego = (const problem_deconv *)p;
  int t;
  Y(md5_put_str)
  (ctx, "deconv");
  Y(md5_put_int)
  (ctx, ego->sign);
  Y(tensor_md5)
  (ctx, ego->sz);
  Y(tensor_md5)
  (ctx, ego->vecsz);
  for (t = 0; t < ego->sz->rnk; t++)
    Y(md5_put_int)
  (ctx, ego->variant[t]);
  Y(md5_put_int)
  (ctx, ego->m);
  Y(md5_put_int)
  (ctx, ego->window);
  /* No M, no x, no fftw_flags: Step A is node- and FFT-independent. */
}

static void dv_print(const problem *p, printer *pr) {
  const problem_deconv *ego = (const problem_deconv *)p;
  int t;
  pr->print(pr, "(deconv sign=%d m=%d ", ego->sign, ego->m);
  Y(tensor_print)
  (ego->sz, pr);
  pr->print(pr, " variant=");
  for (t = 0; t < ego->sz->rnk; t++)
    pr->print(pr, "%d", ego->variant[t]);
  pr->putchr(pr, ')');
}

static void dv_destroy(problem *p) {
  problem_deconv *ego = (problem_deconv *)p;
  Y(tensor_destroy)
  (ego->sz);
  Y(tensor_destroy)
  (ego->vecsz);
  Y(free)
  (ego->variant);
}

static const problem_adt deconv_adt = {
    NFFT_PROBLEM_DECONV, dv_hash, dv_print, dv_destroy};

problem *Y(mkproblem_deconv)(int d, const INT *N, const int *variant,
                             const INT *n, int m, int window, int sign,
                             C *f_hat, C *g) {
  problem_deconv *ego;
  tensor *fwd;
  int t;

  A(sign == 1 || sign == -1);
  A(d >= 1);

  /* forward tensor: factor t has n_in = N[t], n_out = n[t] (row-major strides),
   * same shape as problem_nfft's forward tensor. */
  fwd = Y(tensor_create)(d);
  {
    INT is_acc = (INT)1, os_acc = (INT)1;
    for (t = d - 1; t >= 0; --t) {
      A(N[t] >= (INT)1);
      A(n[t] >= (INT)1);
      fwd->dims[t].n_in = N[t];
      fwd->dims[t].is = is_acc;
      fwd->dims[t].n_out = n[t];
      fwd->dims[t].os = os_acc;
      is_acc *= N[t];
      os_acc *= n[t];
    }
  }

  ego = (problem_deconv *)Y(problem_create)(sizeof(problem_deconv), &deconv_adt);

  if (sign == -1) {
    tensor *adj = Y(tensor_adjoint)(fwd);
    Y(tensor_destroy)
    (fwd);
    ego->sz = adj;
  } else
    ego->sz = fwd;

  ego->vecsz = Y(tensor_create)(0);
  ego->m = m;
  ego->window = window;
  ego->sign = sign;

  ego->variant = (int *)Y(malloc)((size_t)d * sizeof(int));
  for (t = 0; t < d; t++) {
    /* variant == NULL means "all type-I" (mirrors mkproblem_nfft). */
    int vt = variant ? variant[t] : NFFT_NDFT_TYPE_I;
    A(vt == NFFT_NDFT_TYPE_I || vt == NFFT_NDFT_TYPE_II);
    /* type-I and type-II coincide for odd N: normalize so they share a key. */
    ego->variant[t] = (N[t] % 2 == 1) ? NFFT_NDFT_TYPE_I : vt;
  }

  ego->f_hat = f_hat; /* borrowed alias */
  ego->g = g;         /* borrowed alias */

  return (problem *)ego;
}

/* direction-aware accessors: N_t = n_in when forward, n_out when adjoint. */
INT Y(problem_deconv_N)(const problem *p, int t) {
  const problem_deconv *ego = (const problem_deconv *)p;
  A(t >= 0 && t < ego->sz->rnk);
  return ego->sign == 1 ? ego->sz->dims[t].n_in : ego->sz->dims[t].n_out;
}

INT Y(problem_deconv_n)(const problem *p, int t) {
  const problem_deconv *ego = (const problem_deconv *)p;
  A(t >= 0 && t < ego->sz->rnk);
  return ego->sign == 1 ? ego->sz->dims[t].n_out : ego->sz->dims[t].n_in;
}

INT Y(problem_deconv_Ntot)(const problem *p) {
  const problem_deconv *ego = (const problem_deconv *)p;
  return ego->sign == 1 ? Y(tensor_sz_in)(ego->sz) : Y(tensor_sz_out)(ego->sz);
}

INT Y(problem_deconv_ntot)(const problem *p) {
  const problem_deconv *ego = (const problem_deconv *)p;
  return ego->sign == 1 ? Y(tensor_sz_out)(ego->sz) : Y(tensor_sz_in)(ego->sz);
}

int Y(problem_deconv_variant)(const problem *p, int t) {
  const problem_deconv *ego = (const problem_deconv *)p;
  A(t >= 0 && t < ego->sz->rnk);
  return ego->variant[t];
}
