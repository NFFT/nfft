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

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

/* NFFT problem kind. A problem_nfft describes one direction of an NFFT
 * computation. sz is oriented along dataflow: forward (sign +1) has
 * n_in = N_t, n_out = n_t; adjoint (sign -1) stores the tensor_adjoint of it.
 * The direction-aware accessors below are the only sanctioned view of sz. */

/* floor(log2(v)) for v >= 1, without libm. */
static int floor_log2_int(INT v) {
  int r = 0;
  A(v >= (INT)1);
  while (v > (INT)1) {
    v >>= 1;
    ++r;
  }
  return r;
}

/* Unit-axis elision: live[] receives the caller indices of the "live" axes
 * (N_t != 1) in caller order.  Drop-only, no canonical stride sort (unlike
 * FFTW's X(tensor_compress)): x columns pair positionally with axes. Returns
 * their count, which may be 0 (rank-0 problem). */
static int live_axes(int d, const INT *N, int *live) {
  int t, k = 0;
  for (t = 0; t < d; t++)
    if (N[t] != (INT)1)
      live[k++] = t;
  return k;
}

/* problem hash */
static void hash(const problem *p, md5 *ctx) {
  const problem_nfft *ego = (const problem_nfft *)p;

  Y(md5_put_str)
  (ctx, "nfft");
  Y(md5_put_int)
  (ctx, ego->sign);
  Y(tensor_md5)
  (ctx, ego->sz);
  Y(tensor_md5)
  (ctx, ego->vecsz);
  /* variant[] is canonicalized in mkproblem_nfft (odd-N axes forced to
   * type-I), and its run length equals sz->rnk, hashed above -- no marker. */
  {
    int t;
    for (t = 0; t < ego->sz->rnk; t++)
      Y(md5_put_int)
    (ctx, ego->variant[t]);
  }
  Y(md5_put_int)
  (ctx, floor_log2_int(ego->M));
  Y(md5_put_int)
  (ctx, ego->m);
  Y(md5_put_int)
  (ctx, ego->window);
  Y(md5_put_unsigned)
  (ctx, ego->fftw_flags);

  /* x/f_hat/f are not hashed: the wisdom key stays data-blind, so any
   * correctly-shaped arrays reuse the same cached decision.
   * TODO: hash the pointers' alignment class (FFTW parity, never the address);
   * new-array execute must then respect the blessed alignment. */
}

static void print(const problem *p, printer *pr) {
  const problem_nfft *ego = (const problem_nfft *)p;

  pr->print(pr, "(nfft sign=%d m=%d M=%D ", ego->sign, ego->m, ego->M);
  Y(tensor_print)
  (ego->sz, pr);
  {
    int t;
    pr->print(pr, " variant=");
    for (t = 0; t < ego->sz->rnk; t++)
      pr->print(pr, "%d", ego->variant[t]);
  }
  pr->putchr(pr, ')');
}

static void destroy(problem *p) {
  problem_nfft *ego = (problem_nfft *)p;

  Y(tensor_destroy)
  (ego->sz);
  Y(tensor_destroy)
  (ego->vecsz);
  Y(free)
  (ego->variant);
  /* Free x iff this problem owns its copy (top-level).  f_hat/f are borrowed
   * caller arrays. */
  if (ego->x_owned) {
    Y(free)
    (ego->x);
  }
  /* The base allocation is freed by problem_destroy. */
}

static const problem_adt nfft_problem_adt = {
    NFFT_PROBLEM_NFFT,
    hash,
    print,
    destroy};

/* constructor */
problem *Y(mkproblem_nfft)(int d, const INT *N, const int *variant,
                           const INT *n, INT M, int m, int window, int sign,
                           unsigned fftw_flags, R *x, int copy_x, C *f_hat, C *f) {
  problem_nfft *ego;
  tensor *fwd;
  int t;
  int *live;
  int dl; /* live (compressed) rank; may be 0 (all axes unit) */

  A(sign == 1 || sign == -1);
  A(d >= 1); /* caller rank >= 1; only the compressed rank may be 0 */
  A(M >= (INT)1);

  /* The copy_x path elides unit axes from N/n/variant and the x copy below, 
   * so sz is the compressed forward tensor in caller order.
   * The borrowed (copy_x==0) path cannot gather x into a fresh buffer, so it
   * keeps full rank. */
  live = (int *)Y(malloc)((size_t)d * sizeof(int));
  if (copy_x)
    dl = live_axes(d, N, live);
  else {
    for (t = 0; t < d; t++)
      live[t] = t;
    dl = d;
  }

  /* Forward tensor over the live axes. A dropped unit axis contributes
   * stride factor 1, so the surviving strides and the flat f_hat/f layout are
   * unchanged (zero-copy f_hat). */
  fwd = Y(tensor_create)(dl);
  {
    INT is_acc = (INT)1;
    INT os_acc = (INT)1;
    for (t = dl - 1; t >= 0; --t) {
      int c = live[t];
      A(N[c] >= (INT)1);
      A(n[c] >= (INT)1);
      fwd->dims[t].n_in = N[c];
      fwd->dims[t].is = is_acc;
      fwd->dims[t].n_out = n[c];
      fwd->dims[t].os = os_acc;
      is_acc *= N[c];
      os_acc *= n[c];
    }
  }

  ego = (problem_nfft *)Y(problem_create)(sizeof(problem_nfft), &nfft_problem_adt);

  if (sign == -1) {
    tensor *adj = Y(tensor_adjoint)(fwd);
    Y(tensor_destroy)
    (fwd);
    ego->sz = adj;
  } else
    ego->sz = fwd;

  ego->vecsz = Y(tensor_create)(0);
  ego->M = M;
  ego->m = m;
  ego->window = window;
  ego->sign = sign;
  ego->fftw_flags = fftw_flags;

  /* The copy gathers only the live columns.  dl==0 (all axes unit) has no
   * columns: leave x NULL/unowned (the rank-0 solver never reads x). */
  if (copy_x && dl > 0) {
    R *xc = (R *)Y(malloc)((size_t)dl * (size_t)M * sizeof(R));
    INT jj;
    for (jj = 0; jj < M; jj++)
      for (t = 0; t < dl; t++)
        xc[jj * dl + t] = x[jj * d + live[t]];
    ego->x = xc;
    ego->x_owned = 1;
  } else if (copy_x) /* dl == 0 */
  {
    ego->x = 0;
    ego->x_owned = 0;
  } else {
    ego->x = x;
    ego->x_owned = 0;
  }
  ego->f_hat = f_hat;
  ego->f = f;

  ego->variant = (dl > 0) ? (int *)Y(malloc)((size_t)dl * sizeof(int)) : 0;
  for (t = 0; t < dl; t++) {
    int c = live[t];
    int vt = variant ? variant[c] : NFFT_NDFT_TYPE_I;
    A(vt == NFFT_NDFT_TYPE_I || vt == NFFT_NDFT_TYPE_II);
    ego->variant[t] = (N[c] % 2 == 1) ? NFFT_NDFT_TYPE_I : vt;
  }

  Y(free)
  (live);
  return (problem *)ego;
}

/* Direction-aware geometry accessors. Forward: n_in = N_t, n_out = n_t;
 * adjoint (stored after tensor_adjoint): n_in = n_t, n_out = N_t. t indexes
 * the axis list in caller order, range [0, sz->rnk) -- the unit-elided live
 * axes on the copy_x=1 path, all d axes on the borrowed path. */

INT Y(problem_nfft_N)(const problem *p, int t) {
  const problem_nfft *ego = (const problem_nfft *)p;

  A(t >= 0 && t < ego->sz->rnk);
  if (ego->sign == 1)
    return ego->sz->dims[t].n_in;
  else
    return ego->sz->dims[t].n_out;
}

INT Y(problem_nfft_n)(const problem *p, int t) {
  const problem_nfft *ego = (const problem_nfft *)p;

  A(t >= 0 && t < ego->sz->rnk);
  if (ego->sign == 1)
    return ego->sz->dims[t].n_out;
  else
    return ego->sz->dims[t].n_in;
}

INT Y(problem_nfft_Ntot)(const problem *p) {
  const problem_nfft *ego = (const problem_nfft *)p;

  if (ego->sign == 1)
    return Y(tensor_sz_in)(ego->sz);
  else
    return Y(tensor_sz_out)(ego->sz);
}

INT Y(problem_nfft_ntot)(const problem *p) {
  const problem_nfft *ego = (const problem_nfft *)p;

  if (ego->sign == 1)
    return Y(tensor_sz_out)(ego->sz);
  else
    return Y(tensor_sz_in)(ego->sz);
}

int Y(problem_nfft_variant)(const problem *p, int t) {
  const problem_nfft *ego = (const problem_nfft *)p;
  A(t >= 0 && t < ego->sz->rnk);
  return ego->variant[t];
}

int Y(problem_nfft_has_unit_axis)(const problem *p) {
  const problem_nfft *ego = (const problem_nfft *)p;
  int t;
  for (t = 0; t < ego->sz->rnk; t++)
    if (Y(problem_nfft_N)(p, t) == (INT)1)
      return 1;
  return 0;
}
