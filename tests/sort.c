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

/* Tests of flag NFFT_SORT_NODES in NFFT, NFCT and NFST.
 *
 * A plan with NFFT_SORT_NODES must give the same trafo and adjoint as the same
 * plan without it, for every psi mode:
 * - after the precomputation, which sorts the nodes and keeps a sorted copy,
 * - after the nodes have changed and the precomputation has run again,
 *   by *_precompute_one_psi and by the precomputation routine of the mode,
 * - without any precomputation of the nodes (psi modes that allow it), where
 *   the transforms sort the current nodes themselves, also after the nodes
 *   have changed.
 * The trafo of both plans agrees exactly per node and the adjoint up to the
 * order of summation, so the plans are compared with each other at a
 * tolerance of a few hundred ulps, which holds for every window and
 * precision. The unsorted plans are validated against the direct transforms
 * by the other tests. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "sort.h"

#define SORT_SEED 20261028L

/* Cut-off m and oversampling n = 2N, small enough for d = 4 and valid for
 * every window. */
#define SORT_M 4

/* Agreement of the sorted and the unsorted plan, relative l2 distance. */
#define SORT_TOL (K(500.0) * EPSILON)

/* Psi modes, as precompute flags. */
static const unsigned sort_modes[] = {0U, PRE_PSI, PRE_FULL_PSI, FG_PSI,
  FG_PSI | PRE_FG_PSI, PRE_LIN_PSI};
static const char *sort_mode_names[] = {"no psi", "PRE_PSI", "PRE_FULL_PSI",
  "FG_PSI", "PRE_FG_PSI", "PRE_LIN_PSI"};
#define SORT_NMODES (sizeof(sort_modes) / sizeof(sort_modes[0]))

/* Modes whose transforms read no table that depends on the nodes. */
#define SORT_MODE_NODE_FREE(mode) (((mode) & (PRE_PSI | PRE_FULL_PSI \
  | PRE_FG_PSI)) == 0U)

/* Sorted variants: the sort alone and, for NFFT, with the blockwise adjoint
 * (which implies the sort). */
static const unsigned sort_variants[] = {NFFT_SORT_NODES,
  NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT};
static const char *sort_variant_names[] = {"SORT", "SORT|BLOCKWISE"};

/* Runs CHECK(all_stages). The OpenMP library runs all stages on one thread
 * and then, on the default number of threads, all stages in 1D, which covers
 * sorting new nodes again on several threads, but only the first stages for
 * d > 1. The plans are tiny, so threaded
 * transforms and FFTs spend their time in the barriers of the OpenMP runtime,
 * which takes many minutes when other processes compete for the CPUs. FFTW
 * stays on one thread throughout: NFFT plans set the FFTW thread count
 * themselves, but NFCT and NFST plans use the count the last NFFT plan left,
 * which is restored at the end. FFTW 3.3.10 and later report that count;
 * older versions do not, and the count every NFFT, NFSFT and FPT plan sets,
 * the number of threads of the library, is restored instead. */
#ifdef _OPENMP
#if defined(HAVE_FFTW_THREADS)
#define SORT_FFTW_THREADS(n) FFTW(plan_with_nthreads)((int)(n))
#if defined(HAVE_FFTW_PLANNER_NTHREADS)
#define SORT_FFTW_NTHREADS() FFTW(planner_nthreads)()
#else
#define SORT_FFTW_NTHREADS() ((int)Y(get_num_threads)())
#endif
#else
#define SORT_FFTW_THREADS(n)
#define SORT_FFTW_NTHREADS() 1
#endif
#define SORT_THREADED_RUNS(CHECK) \
{ \
  const INT sort_nthreads = Y(get_num_threads)(); \
  const int sort_fftw_nthreads = SORT_FFTW_NTHREADS(); \
  Y(set_num_threads)(1); \
  SORT_FFTW_THREADS(1); \
  CHECK(1); \
  Y(set_num_threads)(sort_nthreads); \
  if (sort_nthreads > 1) \
    CHECK(0); \
  SORT_FFTW_THREADS(sort_fftw_nthreads); \
  (void)sort_fftw_nthreads; \
}
#else
#define SORT_THREADED_RUNS(CHECK) CHECK(1);
#endif

#define SORT_FLAGS (PRE_PHI_HUT | MALLOC_X | MALLOC_F_HAT | MALLOC_F \
  | FFTW_INIT | FFT_OUT_OF_PLACE)
#define SORT_FFTW_FLAGS (FFTW_ESTIMATE | FFTW_DESTROY_INPUT)

static R sort_err_c(const C *a, const C *b, const INT n)
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

static R sort_err_r(const R *a, const R *b, const INT n)
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

static int sort_report(const char *module, const int d, const char *variant,
  const char *mode, const char *stage, const char *which, const R err)
{
  const int ok = (err <= SORT_TOL);

  if (!ok)
    printf("%s_sort d = %d, %-14s, %-12s, %-22s, %-7s -> FAIL %.3e > %.3e\n",
      module, d, variant, mode, stage, which, (double)err, (double)SORT_TOL);

  return ok;
}

/* The check for one module. T is the sample type (C or R), ERR its distance,
 * XLO and XLEN the node interval [XLO, XLO + XLEN), NVARIANTS the number of
 * sorted variants to test. */
#define SORT_CHECK(MODULE, NAME, T, ERR, XLO, XLEN, NVARIANTS) \
static void MODULE(sort_fill)(MODULE(plan) *p, const R *x) \
{ \
  memcpy(p->x, x, (size_t)(p->d * p->M_total) * sizeof(R)); \
} \
\
static void MODULE(sort_random_x)(R *x, const INT n) \
{ \
  INT j; \
  for (j = 0; j < n; j++) \
    x[j] = XLO + XLEN * Y(drand48)(); \
} \
\
static void MODULE(sort_random)(T *v, const INT n) \
{ \
  INT j; \
  for (j = 0; j < n; j++) \
    v[j] = (T)(Y(drand48)() - K(0.5)); \
} \
\
/* the precomputation routine of a single psi mode */ \
static void MODULE(sort_precompute_mode)(MODULE(plan) *p, \
  const unsigned mode) \
{ \
  if (mode & PRE_PSI) \
    MODULE(precompute_psi)(p); \
  else if (mode & PRE_FULL_PSI) \
    MODULE(precompute_full_psi)(p); \
  else if (mode & PRE_FG_PSI) \
    MODULE(precompute_fg_psi)(p); \
} \
\
/* trafo and adjoint of the sorted plan s and the unsorted plan u, on the \
 * same input */ \
static int MODULE(sort_compare)(MODULE(plan) *s, MODULE(plan) *u, \
  const T *f_hat, const T *f, const int d, const char *variant, \
  const char *mode, const char *stage) \
{ \
  int ok = 1; \
  memcpy(s->f_hat, f_hat, (size_t)(s->N_total) * sizeof(T)); \
  memcpy(u->f_hat, f_hat, (size_t)(u->N_total) * sizeof(T)); \
  MODULE(trafo)(s); \
  MODULE(trafo)(u); \
  ok &= sort_report(NAME, d, variant, mode, stage, "trafo", \
    ERR((const T*)s->f, (const T*)u->f, s->M_total)); \
  memcpy(s->f, f, (size_t)(s->M_total) * sizeof(T)); \
  memcpy(u->f, f, (size_t)(u->M_total) * sizeof(T)); \
  MODULE(adjoint)(s); \
  MODULE(adjoint)(u); \
  ok &= sort_report(NAME, d, variant, mode, stage, "adjoint", \
    ERR((const T*)s->f_hat, (const T*)u->f_hat, s->N_total)); \
  /* the adjoint reads f in the caller's order and leaves it as it is */ \
  ok &= sort_report(NAME, d, variant, mode, stage, "f kept", \
    ERR((const T*)s->f, f, s->M_total)); \
  return ok; \
} \
\
static int MODULE(sort_one)(const int d, int *N, int *n, const int M, \
  const unsigned mode, const char *mode_name, const unsigned variant, \
  const char *variant_name, const int all_stages) \
{ \
  MODULE(plan) s, u; \
  INT N_total = 1; \
  int t, ok = 1; \
  R *x0, *x1; \
  T *f_hat, *f; \
\
  for (t = 0; t < d; t++) \
    N_total *= N[t]; \
\
  x0 = (R*)Y(malloc)((size_t)(d * M) * sizeof(R)); \
  x1 = (R*)Y(malloc)((size_t)(d * M) * sizeof(R)); \
  f_hat = (T*)Y(malloc)((size_t)N_total * sizeof(T)); \
  f = (T*)Y(malloc)((size_t)M * sizeof(T)); \
  MODULE(sort_random_x)(x0, (INT)(d * M)); \
  MODULE(sort_random_x)(x1, (INT)(d * M)); \
  MODULE(sort_random)(f_hat, N_total); \
  MODULE(sort_random)(f, (INT)M); \
\
  /* precomputed: the sorted copy of the nodes, then changed nodes and a new \
   * precomputation */ \
  MODULE(init_guru)(&s, d, N, M, n, SORT_M, SORT_FLAGS | mode | variant, \
    SORT_FFTW_FLAGS); \
  MODULE(init_guru)(&u, d, N, M, n, SORT_M, SORT_FLAGS | mode, \
    SORT_FFTW_FLAGS); \
  MODULE(sort_fill)(&s, x0); \
  MODULE(sort_fill)(&u, x0); \
  MODULE(precompute_one_psi)(&s); \
  MODULE(precompute_one_psi)(&u); \
  ok &= MODULE(sort_compare)(&s, &u, f_hat, f, d, variant_name, mode_name, \
    "precomputed"); \
  if (all_stages) \
  { \
    /* twice, as an iterative solver does */ \
    ok &= MODULE(sort_compare)(&s, &u, f_hat, f, d, variant_name, mode_name, \
      "precomputed, again"); \
    MODULE(sort_fill)(&s, x1); \
    MODULE(sort_fill)(&u, x1); \
    MODULE(precompute_one_psi)(&s); \
    MODULE(precompute_one_psi)(&u); \
    ok &= MODULE(sort_compare)(&s, &u, f_hat, f, d, variant_name, mode_name, \
      "new nodes, precomputed"); \
    /* the precomputation routine of the mode on its own, which sorts too */ \
    if (mode & (PRE_PSI | PRE_FULL_PSI | PRE_FG_PSI)) \
    { \
      MODULE(sort_fill)(&s, x0); \
      MODULE(sort_fill)(&u, x0); \
      MODULE(sort_precompute_mode)(&s, mode); \
      MODULE(sort_precompute_mode)(&u, mode); \
      ok &= MODULE(sort_compare)(&s, &u, f_hat, f, d, variant_name, \
        mode_name, "new nodes, own routine"); \
    } \
  } \
  MODULE(finalize)(&s); \
\
  /* never precomputed: the transforms sort the current nodes */ \
  if (SORT_MODE_NODE_FREE(mode)) \
  { \
    MODULE(init_guru)(&s, d, N, M, n, SORT_M, SORT_FLAGS | mode | variant, \
      SORT_FFTW_FLAGS); \
    if (mode & PRE_LIN_PSI) \
      MODULE(precompute_lin_psi)(&s); \
    MODULE(sort_fill)(&s, x0); \
    MODULE(sort_fill)(&u, x0); \
    ok &= MODULE(sort_compare)(&s, &u, f_hat, f, d, variant_name, \
      mode_name, "not precomputed"); \
    if (all_stages) \
    { \
      MODULE(sort_fill)(&s, x1); \
      MODULE(sort_fill)(&u, x1); \
      ok &= MODULE(sort_compare)(&s, &u, f_hat, f, d, variant_name, \
        mode_name, "not precomputed, new"); \
    } \
    MODULE(finalize)(&s); \
  } \
\
  MODULE(finalize)(&u); \
  Y(free)(x0); \
  Y(free)(x1); \
  Y(free)(f_hat); \
  Y(free)(f); \
  return ok; \
} \
\
/* every d, psi mode and sorted variant; all_stages also changes the nodes \
 * and runs the precomputation routines of the modes for d > 1, which 1D \
 * always does */ \
static int MODULE(sort_all)(const int all_stages) \
{ \
  /* N per dimension and M, for d = 1, ..., 4; n = 2N */ \
  static const int NN[] = {64, 16, 8, 8}; \
  static const int MM[] = {300, 300, 200, 150}; \
  int d, ok = 1; \
  size_t i, v; \
\
  for (d = 1; d <= 4; d++) \
  { \
    int N[4], n[4], t; \
    for (t = 0; t < d; t++) \
    { \
      N[t] = NN[d - 1]; \
      n[t] = 2 * NN[d - 1]; \
    } \
    for (i = 0; i < SORT_NMODES; i++) \
      for (v = 0; v < (size_t)(NVARIANTS); v++) \
        ok &= MODULE(sort_one)(d, N, n, MM[d - 1], sort_modes[i], \
          sort_mode_names[i], sort_variants[v], sort_variant_names[v], \
          all_stages || d == 1); \
  } \
\
  return ok; \
} \
\
void MODULE(check_sort)(void) \
{ \
  int ok = 1; \
\
  Y(srand48)(SORT_SEED); \
  SORT_THREADED_RUNS(ok &= MODULE(sort_all)) \
\
  CU_ASSERT(ok); \
}

#undef X
#define X(name) NFFT(name)
SORT_CHECK(X, "nfft", C, sort_err_c, K(-0.5), K(1.0), 2)

#ifdef HAVE_NFCT
#undef X
#define X(name) NFCT(name)
SORT_CHECK(X, "nfct", R, sort_err_r, K(0.0), K(0.5), 1)
#endif

#ifdef HAVE_NFST
#undef X
#define X(name) NFST(name)
SORT_CHECK(X, "nfst", R, sort_err_r, K(0.0), K(0.5), 1)
#endif

#ifdef SORT_NNFFT
/* NNFFT: its inner NFFT plan gets NFFT_SORT_NODES from
 * NFFT_OMP_BLOCKWISE_ADJOINT and is precomputed with the nodes divided by the
 * oversampling factors, which the transforms divide by again. Compares with
 * the same NNFFT plan whose inner plan has the sort switched off. */
static void sort_nnfft_x(double *x, const int n)
{
  int j;
  for (j = 0; j < n; j++)
    x[j] = -0.5 + Y(drand48)();
}

static int sort_nnfft_compare(nnfft_plan *s, nnfft_plan *u,
  const double _Complex *f_hat, const double _Complex *f, const int d,
  const char *mode, const char *stage)
{
  int ok = 1;
  memcpy(s->f_hat, f_hat, (size_t)(s->N_total) * sizeof(double _Complex));
  memcpy(u->f_hat, f_hat, (size_t)(u->N_total) * sizeof(double _Complex));
  nnfft_trafo(s);
  nnfft_trafo(u);
  ok &= sort_report("nnfft", d, "SORT", mode, stage, "trafo",
    sort_err_c((const C*)s->f, (const C*)u->f, s->M_total));
  memcpy(s->f, f, (size_t)(s->M_total) * sizeof(double _Complex));
  memcpy(u->f, f, (size_t)(u->M_total) * sizeof(double _Complex));
  nnfft_adjoint(s);
  nnfft_adjoint(u);
  ok &= sort_report("nnfft", d, "SORT", mode, stage, "adjoint",
    sort_err_c((const C*)s->f_hat, (const C*)u->f_hat, s->N_total));
  return ok;
}

void nnfft_check_sort(void)
{
  static const unsigned modes[] = {0U, PRE_PSI, PRE_FULL_PSI};
  static const char *mode_names[] = {"no psi", "PRE_PSI", "PRE_FULL_PSI"};
  const unsigned unsorted = NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT;
  const int M = 200, N_total = 60;
  int d, ok = 1;
  size_t i;

  Y(srand48)(SORT_SEED);

  for (d = 1; d <= 3; d++)
  {
    int N[3], N1[3], t, j;
    double *x0 = (double*)Y(malloc)((size_t)(d * M) * sizeof(double));
    double *x1 = (double*)Y(malloc)((size_t)(d * M) * sizeof(double));
    double *v = (double*)Y(malloc)((size_t)(d * N_total) * sizeof(double));
    double _Complex *f_hat = (double _Complex*)Y(malloc)((size_t)N_total
      * sizeof(double _Complex));
    double _Complex *f = (double _Complex*)Y(malloc)((size_t)M
      * sizeof(double _Complex));

    for (t = 0; t < d; t++)
    {
      N[t] = 12;
      N1[t] = 24;
    }
    sort_nnfft_x(x0, d * M);
    sort_nnfft_x(x1, d * M);
    sort_nnfft_x(v, d * N_total);
    for (j = 0; j < N_total; j++)
      f_hat[j] = (Y(drand48)() - 0.5) + (Y(drand48)() - 0.5) * _Complex_I;
    for (j = 0; j < M; j++)
      f[j] = (Y(drand48)() - 0.5) + (Y(drand48)() - 0.5) * _Complex_I;

    for (i = 0; i < sizeof(modes) / sizeof(modes[0]); i++)
    {
      nnfft_plan s, u;
      const unsigned flags = PRE_PHI_HUT | MALLOC_X | MALLOC_V | MALLOC_F
        | MALLOC_F_HAT | modes[i];
      unsigned u_flags;

      nnfft_init_guru(&s, d, N_total, M, N, N1, SORT_M, flags);
      nnfft_init_guru(&u, d, N_total, M, N, N1, SORT_M, flags);
      ok &= ((s.direct_plan->flags & NFFT_SORT_NODES) != 0U);
      u_flags = u.direct_plan->flags;
      u.direct_plan->flags &= ~unsorted;

      memcpy(s.v, v, (size_t)(d * N_total) * sizeof(double));
      memcpy(u.v, v, (size_t)(d * N_total) * sizeof(double));
      memcpy(s.x, x0, (size_t)(d * M) * sizeof(double));
      memcpy(u.x, x0, (size_t)(d * M) * sizeof(double));
      nnfft_precompute_one_psi(&s);
      nnfft_precompute_one_psi(&u);
      ok &= sort_nnfft_compare(&s, &u, f_hat, f, d, mode_names[i],
        "precomputed");
      memcpy(s.x, x1, (size_t)(d * M) * sizeof(double));
      memcpy(u.x, x1, (size_t)(d * M) * sizeof(double));
      nnfft_precompute_one_psi(&s);
      nnfft_precompute_one_psi(&u);
      ok &= sort_nnfft_compare(&s, &u, f_hat, f, d, mode_names[i],
        "new nodes, precomputed");

      u.direct_plan->flags = u_flags;
      nnfft_finalize(&s);
      nnfft_finalize(&u);
    }

    Y(free)(x0);
    Y(free)(x1);
    Y(free)(v);
    Y(free)(f_hat);
    Y(free)(f);
  }

  CU_ASSERT(ok);
}
#endif
