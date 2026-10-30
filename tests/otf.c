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

/* Tests of the window evaluated a block of nodes at a time, by the
 * transforms without a psi table and by the fill of the PRE_PSI table.
 *
 * - PHI_RUN_BLOCK, the window runs of a block of nodes at once, against
 *   PHI_RUN node by node, for every block length and for cut-off parameters on
 *   both sides of the limit of the Kaiser-Bessel block evaluator.
 * - The plans of *_init, *_init_1d, *_init_2d and *_init_3d hold a PRE_PSI
 *   table and sort the nodes, except for the 1D NFFT of the serial library.
 *   They and the same plans without a psi flag give the same trafo and
 *   adjoint as the same plans with PRE_FULL_PSI, whose table is filled node
 *   by node with PHI_RUN, so neither side shares the block evaluator with the
 *   reference. In 1D the plans without a psi flag are also tried with the
 *   nodes sorted the other way round (unsorted where the default sorts). The
 *   node counts include partial blocks and fewer nodes than one block.
 * - For d > 1 the PRE_PSI plan, whose table is filled a block of nodes at a
 *   time, and the plan without a psi flag, which forms the window node by
 *   node with PHI_RUN, add the same taps in the same order; they must agree
 *   to a few hundred ulps.
 * - *_precompute_psi on a plan without PRE_PSI only sorts the nodes.
 * Tolerances (otf_tol): against PRE_FULL_PSI, which adds the (2m+2)^d taps
 * of a node in another order than the tensor-product kernel, the rounding of
 * these sums shows through their cancellation, up to about ten ulps per tap
 * with the wide windows (Gaussian, m = 13 and 17 in long double; sinc, m = 11
 * and 13) on the small grids here; hence 32 ulps per tap, and no less than
 * 1000 ulps. Plans that add in the same order: 500 ulps. Each test prints its
 * largest distance / tolerance. The PRE_FULL_PSI plans are validated against
 * the direct transforms by the other tests. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "otf.h"

#define OTF_SEED 20261029L

/* Agreement of two plans of dimension d and cut-off m, relative l2 distance
 * (otf_err_*): same_order whether they add the taps of a node in the same
 * order. */
static R otf_tol(const int d, const int m, const int same_order)
{
  if (same_order)
    return K(500.0) * EPSILON;
  return FMAX(K(1000.0), K(32.0) * POW((R)(2 * m + 2), (R)d)) * EPSILON;
}

/* The largest distance / tolerance of the current test, against PRE_FULL_PSI
 * [0] and in the same order [1], printed at its end. */
static R otf_worst[2];

/* Agreement of a window run filled block-wise and node by node, relative to
 * the largest value of the run. */
#define OTF_RUN_TOL (K(16.0) * EPSILON)

/* l2 distance of the outputs a and b of length n of a transform of dimension
 * d, relative to the l2 norm of b, or to that expected for the random input
 * in of length nin when it is larger: each output adds about nin random
 * terms, with a mean square of at least 2^-d (the cosines and sines of the
 * NFCT and NFST), so the norm of b is about |in| sqrt(n) 2^(-d/2) or more.
 * Without that floor one output that cancels (M = 1) would inflate the
 * distance relative to it. */
static R otf_err_c(const C *a, const C *b, const INT n, const C *in,
  const INT nin, const int d)
{
  R num = K(0.0), den = K(0.0), inn = K(0.0);
  INT j;

  for (j = 0; j < n; j++)
  {
    num += CABS(a[j] - b[j]) * CABS(a[j] - b[j]);
    den += CABS(b[j]) * CABS(b[j]);
  }
  for (j = 0; j < nin; j++)
    inn += CABS(in[j]) * CABS(in[j]);
  den = FMAX(den, inn * (R)n * POW(K(0.5), (R)d));

  return (den > K(0.0)) ? SQRT(num / den) : SQRT(num);
}

static R otf_err_r(const R *a, const R *b, const INT n, const R *in,
  const INT nin, const int d)
{
  R num = K(0.0), den = K(0.0), inn = K(0.0);
  INT j;

  for (j = 0; j < n; j++)
  {
    num += (a[j] - b[j]) * (a[j] - b[j]);
    den += b[j] * b[j];
  }
  for (j = 0; j < nin; j++)
    inn += in[j] * in[j];
  den = FMAX(den, inn * (R)n * POW(K(0.5), (R)d));

  return (den > K(0.0)) ? SQRT(num / den) : SQRT(num);
}

static int otf_report(const char *module, const int d, const int M,
  const char *stage, const char *which, const R err, const R tol,
  const int same_order)
{
  const int ok = (err <= tol);

  otf_worst[same_order] = FMAX(otf_worst[same_order], err / tol);
  if (!ok)
    printf("%s_otf d = %d, M = %4d, %-26s, %-7s -> FAIL %.3e > %.3e\n",
      module, d, M, stage, which, (double)err, (double)tol);

  return ok;
}

/* PHI_RUN_BLOCK against PHI_RUN for the 1D plan ths, all block lengths nb, a
 * row stride ld > 2m+2; the entries between the runs must stay untouched. */
static int otf_run_block_one(const NFFT(plan) *ths)
{
  const INT m = ths->m, L = 2 * m + 2, ld = L + 3, n = ths->n[0];
  R x[NFFT_NODE_BLOCK], *blk, *one;
  INT u[NFFT_NODE_BLOCK], nb, v, l;
  int ok = 1;

  blk = (R*)Y(malloc)((size_t)(NFFT_NODE_BLOCK * ld) * sizeof(R));
  one = (R*)Y(malloc)((size_t)L * sizeof(R));

  for (nb = 1; nb <= NFFT_NODE_BLOCK; nb++)
  {
    for (v = 0; v < nb; v++)
    {
      x[v] = Y(drand48)() - K(0.5);
      u[v] = (INT)FLOOR(x[v] * (R)n) - m;
    }
    for (l = 0; l < NFFT_NODE_BLOCK * ld; l++)
      blk[l] = K(-7.0);

    PHI_RUN_BLOCK(blk, ld, n, x, u, nb, 0);

    for (v = 0; v < NFFT_NODE_BLOCK; v++)
    {
      R err = K(0.0), peak = K(0.0);

      if (v < nb)
      {
        PHI_RUN(one, n, x[v], u[v], 0);
        for (l = 0; l < L; l++)
        {
          err = FMAX(err, FABS(blk[v * ld + l] - one[l]));
          peak = FMAX(peak, FABS(one[l]));
        }
        if (err > OTF_RUN_TOL * peak)
        {
          printf("otf_run_block m = %d, nb = %d, node %d -> FAIL %.3e > "
            "%.3e\n", (int)m, (int)nb, (int)v, (double)err,
            (double)(OTF_RUN_TOL * peak));
          ok = 0;
        }
      }
      for (l = (v < nb) ? L : 0; l < ld; l++)
        if (blk[v * ld + l] != K(-7.0))
        {
          printf("otf_run_block m = %d, nb = %d, node %d, entry %d "
            "written -> FAIL\n", (int)m, (int)nb, (int)v, (int)l);
          ok = 0;
        }
    }
  }

  Y(free)(one);
  Y(free)(blk);
  return ok;
}

void Y(check_window_run_block)(void)
{
  /* both sides of KB_POLY_RUN_BLOCK_MAX_M */
  static const int ms[] = {2, 4, 6, 8, 16, 17, 20};
  int i, ok = 1;

  Y(srand48)(OTF_SEED);

  for (i = 0; i < (int)(sizeof(ms) / sizeof(ms[0])); i++)
  {
    NFFT(plan) p;
    int N = 128, n = 256;

    NFFT(init_guru)(&p, 1, &N, 1, &n, ms[i], PRE_PHI_HUT | MALLOC_X
      | MALLOC_F_HAT | MALLOC_F, FFTW_ESTIMATE);
    {
      const NFFT(plan) *ths = &p;
      ok &= otf_run_block_one(ths);
    }
    NFFT(finalize)(&p);
  }

  CU_ASSERT(ok);
}

/* The check for one module. T is the sample type (C or R), ERR its distance,
 * XLO and XLEN the node interval [XLO, XLO + XLEN). */
#define OTF_CHECK(MODULE, NAME, T, ERR, XLO, XLEN, SORT1D) \
static void MODULE(otf_random)(T *v, const INT n) \
{ \
  INT j; \
  for (j = 0; j < n; j++) \
    v[j] = (T)(Y(drand48)() - K(0.5)); \
} \
\
/* trafo and adjoint of plans a and b, on the same input; same_order as for \
 * otf_tol */ \
static int MODULE(otf_compare)(MODULE(plan) *a, MODULE(plan) *b, \
  const T *f_hat, const T *f, const int d, const char *stage, \
  const int same_order) \
{ \
  int ok = 1; \
  const int M = (int)a->M_total; \
  const R tol = otf_tol(d, (int)a->m, same_order); \
  memcpy(a->f_hat, f_hat, (size_t)(a->N_total) * sizeof(T)); \
  memcpy(b->f_hat, f_hat, (size_t)(b->N_total) * sizeof(T)); \
  MODULE(trafo)(a); \
  MODULE(trafo)(b); \
  ok &= otf_report(NAME, d, M, stage, "trafo", \
    ERR((const T*)a->f, (const T*)b->f, a->M_total, f_hat, a->N_total, d), \
    tol, same_order); \
  memcpy(a->f, f, (size_t)(a->M_total) * sizeof(T)); \
  memcpy(b->f, f, (size_t)(b->M_total) * sizeof(T)); \
  MODULE(adjoint)(a); \
  MODULE(adjoint)(b); \
  ok &= otf_report(NAME, d, M, stage, "adjoint", \
    ERR((const T*)a->f_hat, (const T*)b->f_hat, a->N_total, f, a->M_total, \
    d), tol, same_order); \
  return ok; \
} \
\
static void MODULE(otf_init_default)(MODULE(plan) *p, const int d, int *N, \
  const int M, const int which) \
{ \
  switch (which) \
  { \
    case 1: \
      if (d == 1) { MODULE(init_1d)(p, N[0], M); return; } \
      if (d == 2) { MODULE(init_2d)(p, N[0], N[1], M); return; } \
      MODULE(init_3d)(p, N[0], N[1], N[2], M); return; \
    default: \
      MODULE(init)(p, d, N, M); \
  } \
} \
\
static int MODULE(otf_one)(const int d, int *N, const int M, \
  const int which, const int toggle_sort) \
{ \
  MODULE(plan) a, b, c; \
  int n[3], t, ok = 1; \
  unsigned flags0; \
  T *f_hat, *f; \
  R *x; \
  const char *init = which ? "init_<d>d" : "init"; \
\
  MODULE(otf_init_default)(&a, d, N, M, which); \
\
  /* the default: a PRE_PSI table, the nodes sorted for d > 1, and in 1D if \
   * SORT1D */ \
  if (!(a.flags & PRE_PSI) || !(a.flags & NFFT_SORT_NODES) != (d == 1 \
    && !(SORT1D))) \
  { \
    printf("%s_otf d = %d, %s: default flags %#x -> FAIL\n", NAME, d, \
      init, a.flags); \
    ok = 0; \
  } \
\
  /* the same plan, with the nodes sorted the other way round if \
   * toggle_sort, without a psi table, and with PRE_FULL_PSI \
   * (NFFT_OMP_BLOCKWISE_ADJOINT implies NFFT_SORT_NODES) */ \
  flags0 = a.flags & ~(unsigned)PRE_ONE_PSI; \
  if (toggle_sort) \
    flags0 = (flags0 & NFFT_SORT_NODES) ? flags0 & ~(unsigned)(NFFT_SORT_NODES \
      | NFFT_OMP_BLOCKWISE_ADJOINT) : flags0 | NFFT_SORT_NODES; \
  for (t = 0; t < d; t++) \
    n[t] = (int)a.n[t]; \
  MODULE(init_guru)(&b, d, N, M, n, (int)a.m, flags0, a.fftw_flags); \
  MODULE(init_guru)(&c, d, N, M, n, (int)a.m, flags0 | PRE_FULL_PSI, \
    a.fftw_flags); \
\
  x = (R*)Y(malloc)((size_t)(d * M) * sizeof(R)); \
  f_hat = (T*)Y(malloc)((size_t)(a.N_total) * sizeof(T)); \
  f = (T*)Y(malloc)((size_t)M * sizeof(T)); \
  for (t = 0; t < d * M; t++) \
    x[t] = XLO + XLEN * Y(drand48)(); \
  MODULE(otf_random)(f_hat, a.N_total); \
  MODULE(otf_random)(f, (INT)M); \
\
  memcpy(a.x, x, (size_t)(d * M) * sizeof(R)); \
  memcpy(b.x, x, (size_t)(d * M) * sizeof(R)); \
  memcpy(c.x, x, (size_t)(d * M) * sizeof(R)); \
  MODULE(precompute_one_psi)(&a); \
  MODULE(precompute_one_psi)(&b); \
  MODULE(precompute_one_psi)(&c); \
  ok &= MODULE(otf_compare)(&a, &c, f_hat, f, d, init, 0); \
  ok &= MODULE(otf_compare)(&b, &c, f_hat, f, d, "on the fly", 0); \
  if (d > 1 && !toggle_sort) \
    ok &= MODULE(otf_compare)(&a, &b, f_hat, f, d, "on the fly, same order", \
      1); \
\
  /* new nodes; precompute_psi fills the table of a and only sorts the \
   * nodes of b */ \
  for (t = 0; t < d * M; t++) \
    x[t] = XLO + XLEN * Y(drand48)(); \
  memcpy(a.x, x, (size_t)(d * M) * sizeof(R)); \
  memcpy(b.x, x, (size_t)(d * M) * sizeof(R)); \
  memcpy(c.x, x, (size_t)(d * M) * sizeof(R)); \
  MODULE(precompute_psi)(&a); \
  MODULE(precompute_psi)(&b); \
  MODULE(precompute_full_psi)(&c); \
  ok &= MODULE(otf_compare)(&a, &c, f_hat, f, d, "precompute_psi", 0); \
  ok &= MODULE(otf_compare)(&b, &c, f_hat, f, d, \
    "precompute_psi, no table", 0); \
  if (d > 1 && !toggle_sort) \
    ok &= MODULE(otf_compare)(&a, &b, f_hat, f, d, \
      "precompute_psi, same order", 1); \
\
  MODULE(finalize)(&a); \
  MODULE(finalize)(&b); \
  MODULE(finalize)(&c); \
  Y(free)(x); \
  Y(free)(f_hat); \
  Y(free)(f); \
  return ok; \
} \
\
void MODULE(check_otf)(void) \
{ \
  /* N per dimension, for d = 1, 2, 3; M below, at and above one block */ \
  static const int NN[] = {64, 16, 8}; \
  static const int MM[] = {1, 31, 32, 33, 100}; \
  int d, ok = 1; \
  size_t i; \
\
  Y(srand48)(OTF_SEED); \
  otf_worst[0] = otf_worst[1] = K(0.0); \
\
  for (d = 1; d <= 3; d++) \
  { \
    int N[3], t, which; \
    for (t = 0; t < d; t++) \
      N[t] = NN[d - 1]; \
    for (i = 0; i < sizeof(MM) / sizeof(MM[0]); i++) \
      for (which = 0; which <= 1; which++) \
      { \
        ok &= MODULE(otf_one)(d, N, MM[i], which, 0); \
        if (d == 1) \
          ok &= MODULE(otf_one)(d, N, MM[i], which, 1); \
      } \
  } \
\
  printf("%s_otf: largest distance / tolerance %.3f, in the same order " \
    "%.3f\n", NAME, (double)otf_worst[0], (double)otf_worst[1]); \
  CU_ASSERT(ok); \
}

#undef X
#define X(name) NFFT(name)
/* the 1D NFFT sorts by default only in the OpenMP library */
#ifdef _OPENMP
#define OTF_NFFT_SORT1D 1
#else
#define OTF_NFFT_SORT1D 0
#endif
OTF_CHECK(X, "nfft", C, otf_err_c, K(-0.5), K(1.0), OTF_NFFT_SORT1D)

#ifdef HAVE_NFCT
#undef X
#define X(name) NFCT(name)
OTF_CHECK(X, "nfct", R, otf_err_r, K(0.0), K(0.5), 1)
#endif

#ifdef HAVE_NFST
#undef X
#define X(name) NFST(name)
OTF_CHECK(X, "nfst", R, otf_err_r, K(0.0), K(0.5), 1)
#endif
