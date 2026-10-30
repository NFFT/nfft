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
 *   table, and sort the nodes for d > 1. They and the same plans without a
 *   psi flag give the same trafo and adjoint as the same plans with
 *   PRE_FULL_PSI, whose table is filled node by node with PHI_RUN, so neither
 *   side shares the block evaluator with the reference. In 1D the plans
 *   without a psi flag are also tried with NFFT_SORT_NODES. The node counts
 *   include partial blocks and fewer nodes than one block.
 * - *_precompute_psi on a plan without PRE_PSI only sorts the nodes.
 * The plans are compared with each other at a tolerance of a few hundred ulps,
 * as in the tests of NFFT_SORT_NODES; the PRE_FULL_PSI plans are validated
 * against the direct transforms by the other tests. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "otf.h"

#define OTF_SEED 20261029L

/* Agreement of two plans, relative l2 distance. */
#define OTF_TOL (K(500.0) * EPSILON)

/* Agreement of a window run filled block-wise and node by node, relative to
 * the largest value of the run. */
#define OTF_RUN_TOL (K(16.0) * EPSILON)

static R otf_err_c(const C *a, const C *b, const INT n)
{
  R num = K(0.0), den = K(0.0);
  INT j;

  for (j = 0; j < n; j++)
  {
    num += CABS(a[j] - b[j]) * CABS(a[j] - b[j]);
    den += CABS(b[j]) * CABS(b[j]);
  }

  return (den > K(0.0)) ? SQRT(num / den) : SQRT(num);
}

static R otf_err_r(const R *a, const R *b, const INT n)
{
  R num = K(0.0), den = K(0.0);
  INT j;

  for (j = 0; j < n; j++)
  {
    num += (a[j] - b[j]) * (a[j] - b[j]);
    den += b[j] * b[j];
  }

  return (den > K(0.0)) ? SQRT(num / den) : SQRT(num);
}

static int otf_report(const char *module, const int d, const int M,
  const char *stage, const char *which, const R err)
{
  const int ok = (err <= OTF_TOL);

  if (!ok)
    printf("%s_otf d = %d, M = %4d, %-24s, %-7s -> FAIL %.3e > %.3e\n",
      module, d, M, stage, which, (double)err, (double)OTF_TOL);

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
#define OTF_CHECK(MODULE, NAME, T, ERR, XLO, XLEN) \
static void MODULE(otf_random)(T *v, const INT n) \
{ \
  INT j; \
  for (j = 0; j < n; j++) \
    v[j] = (T)(Y(drand48)() - K(0.5)); \
} \
\
/* trafo and adjoint of plans a and b, on the same input */ \
static int MODULE(otf_compare)(MODULE(plan) *a, MODULE(plan) *b, \
  const T *f_hat, const T *f, const int d, const char *stage) \
{ \
  int ok = 1; \
  const int M = (int)a->M_total; \
  memcpy(a->f_hat, f_hat, (size_t)(a->N_total) * sizeof(T)); \
  memcpy(b->f_hat, f_hat, (size_t)(b->N_total) * sizeof(T)); \
  MODULE(trafo)(a); \
  MODULE(trafo)(b); \
  ok &= otf_report(NAME, d, M, stage, "trafo", \
    ERR((const T*)a->f, (const T*)b->f, a->M_total)); \
  memcpy(a->f, f, (size_t)(a->M_total) * sizeof(T)); \
  memcpy(b->f, f, (size_t)(b->M_total) * sizeof(T)); \
  MODULE(adjoint)(a); \
  MODULE(adjoint)(b); \
  ok &= otf_report(NAME, d, M, stage, "adjoint", \
    ERR((const T*)a->f_hat, (const T*)b->f_hat, a->N_total)); \
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
  const int which, const unsigned extra) \
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
  /* the default: a PRE_PSI table, the nodes sorted for d > 1 */ \
  if (!(a.flags & PRE_PSI) || !(a.flags & NFFT_SORT_NODES) != (d == 1)) \
  { \
    printf("%s_otf d = %d, %s: default flags %#x -> FAIL\n", NAME, d, \
      init, a.flags); \
    ok = 0; \
  } \
\
  /* the same plan, plus the flags extra, without a psi table, and with \
   * PRE_FULL_PSI */ \
  flags0 = (a.flags & ~(unsigned)PRE_ONE_PSI) | extra; \
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
  ok &= MODULE(otf_compare)(&a, &c, f_hat, f, d, init); \
  ok &= MODULE(otf_compare)(&b, &c, f_hat, f, d, "on the fly"); \
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
  ok &= MODULE(otf_compare)(&a, &c, f_hat, f, d, "precompute_psi"); \
  ok &= MODULE(otf_compare)(&b, &c, f_hat, f, d, \
    "precompute_psi, no table"); \
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
\
  for (d = 1; d <= 3; d++) \
  { \
    int N[3], t, which; \
    for (t = 0; t < d; t++) \
      N[t] = NN[d - 1]; \
    for (i = 0; i < sizeof(MM) / sizeof(MM[0]); i++) \
      for (which = 0; which <= 1; which++) \
      { \
        ok &= MODULE(otf_one)(d, N, MM[i], which, 0U); \
        if (d == 1) \
          ok &= MODULE(otf_one)(d, N, MM[i], which, NFFT_SORT_NODES); \
      } \
  } \
\
  CU_ASSERT(ok); \
}

#undef X
#define X(name) NFFT(name)
OTF_CHECK(X, "nfft", C, otf_err_c, K(-0.5), K(1.0))

#ifdef HAVE_NFCT
#undef X
#define X(name) NFCT(name)
OTF_CHECK(X, "nfct", R, otf_err_r, K(0.0), K(0.5))
#endif

#ifdef HAVE_NFST
#undef X
#define X(name) NFST(name)
OTF_CHECK(X, "nfst", R, otf_err_r, K(0.0), K(0.5))
#endif
