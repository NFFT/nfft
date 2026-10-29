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

/* Tests of the tensor product kernels of the matrices B and B^T in NFFT, NFCT
 * and NFST.
 *
 * The kernels run over the rows of the window along the last dimension with
 * loops whose length is the row length (2m+2 points), which the compiler
 * splits into vector code and a remainder. For a sweep of window lengths m
 * and d = 1, ..., 4, the fast transforms with the kernels (no psi table and
 * PRE_PSI, with and without sorting, for NFFT also the blockwise adjoint)
 * must agree to rounding with the PRE_FULL_PSI transform of the same plan
 * parameters, whose B is a plain gather over the precomputed indices and
 * products of the same window values. A tap dropped, counted twice or given
 * a wrong weight shows as an error of the order of the samples, whatever the
 * accuracy of the window for that m. Only failures are printed (all results
 * with -DKERNEL_VERBOSE). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "kernel.h"

#define KERNEL_SEED 20261030L
#define KERNEL_M 37 /* nodes */

/* agreement of two transforms that sum the same terms in a different order */
#define KERNEL_TOL (K(1000.0) * EPSILON)

/* window lengths: d <= 3, and d = 4 */
static const int kernel_m[] = {1, 2, 3, 4, 5, 6, 7, 8, 11};
static const int kernel_m_4d[] = {1, 2, 3};

/* psi modes and variants compared with PRE_FULL_PSI (blockwise: NFFT only,
 * the last one) */
static const unsigned kernel_modes[] = {0U, PRE_PSI};
static const char *kernel_mode_names[] = {"no psi", "PRE_PSI"};
static const unsigned kernel_variants[] = {0U, NFFT_SORT_NODES,
  NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT};
static const char *kernel_variant_names[] = {"", "SORT", "SORT|BLOCKWISE"};

#define KERNEL_FLAGS (PRE_PHI_HUT | MALLOC_X | MALLOC_F_HAT | MALLOC_F \
  | FFTW_INIT | FFT_OUT_OF_PLACE)
#define KERNEL_FFTW_FLAGS (FFTW_ESTIMATE | FFTW_DESTROY_INPUT)

static R kernel_err_c(const C *a, const C *b, const INT n)
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

static R kernel_err_r(const R *a, const R *b, const INT n)
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

static int kernel_report(const char *module, const int d, const int m,
  const char *mode, const char *variant, const char *which, const R err)
{
  const int ok = (err <= KERNEL_TOL);

#ifndef KERNEL_VERBOSE
  if (ok)
    return ok;
#endif

  printf("%s_kernel d = %d, m = %2d, %-8s %-15s %-7s -> %-4s %.3e (%.3e)\n",
    module, d, m, mode, variant, which, ok ? "OK" : "FAIL", (double)err,
    (double)KERNEL_TOL);

  return ok;
}

/* The check for one module. T is the sample type (C or R), ERR its distance,
 * NVARIANTS the number of variants, OFS = n - 2N per dimension; the nodes
 * lie in [XLO, XLO + XW). */
#define KERNEL_CHECK(MODULE, NAME, T, ERR, NVARIANTS, OFS, XLO, XW) \
static void MODULE(kernel_random)(T *v, const INT n) \
{ \
  INT j; \
  for (j = 0; j < n; j++) \
    v[j] = (T)(Y(drand48)() - K(0.5)); \
} \
\
/* trafo and adjoint of a plan with the given precompute flags; returns the \
 * number of coefficients (NFST: N_t - 1 per dimension) */ \
static INT MODULE(kernel_run)(const int d, int *N, int *n, const int m, \
  const unsigned flags, const R *x, const T *f_hat, const T *f, T *tr, \
  T *ad) \
{ \
  MODULE(plan) p; \
  INT n_hat; \
\
  MODULE(init_guru)(&p, d, N, KERNEL_M, n, m, KERNEL_FLAGS | flags, \
    KERNEL_FFTW_FLAGS); \
  memcpy(p.x, x, (size_t)(d * KERNEL_M) * sizeof(R)); \
  MODULE(precompute_one_psi)(&p); \
\
  memcpy(p.f_hat, f_hat, (size_t)(p.N_total) * sizeof(T)); \
  MODULE(trafo)(&p); \
  memcpy(tr, p.f, (size_t)(KERNEL_M) * sizeof(T)); \
\
  memcpy(p.f, f, (size_t)(KERNEL_M) * sizeof(T)); \
  MODULE(adjoint)(&p); \
  memcpy(ad, p.f_hat, (size_t)(p.N_total) * sizeof(T)); \
  n_hat = p.N_total; \
\
  MODULE(finalize)(&p); \
  return n_hat; \
} \
\
static int MODULE(kernel_one)(const int d, const int m) \
{ \
  int N[4], n[4], t, ok = 1; \
  INT N_total = 1, n_hat; \
  size_t i, v; \
  R *x; \
  T *f_hat, *f, *tr0, *ad0, *tr, *ad; \
\
  /* the smallest even N above m; the fast NFFT needs N > m and \
   * n > 2m+2, so n = 2N+2 there */ \
  for (t = 0; t < d; t++) \
  { \
    N[t] = 2 * (m / 2 + 1); \
    n[t] = 2 * N[t] + (OFS); \
    N_total *= N[t]; \
  } \
\
  x = (R*)Y(malloc)((size_t)(d * KERNEL_M) * sizeof(R)); \
  f_hat = (T*)Y(malloc)((size_t)N_total * sizeof(T)); \
  f = (T*)Y(malloc)((size_t)(KERNEL_M) * sizeof(T)); \
  tr0 = (T*)Y(malloc)((size_t)(KERNEL_M) * sizeof(T)); \
  tr = (T*)Y(malloc)((size_t)(KERNEL_M) * sizeof(T)); \
  ad0 = (T*)Y(malloc)((size_t)N_total * sizeof(T)); \
  ad = (T*)Y(malloc)((size_t)N_total * sizeof(T)); \
\
  for (t = 0; t < d * KERNEL_M; t++) \
    x[t] = (XLO) + (XW) * Y(drand48)(); \
  MODULE(kernel_random)(f_hat, N_total); \
  MODULE(kernel_random)(f, (INT)(KERNEL_M)); \
\
  n_hat = MODULE(kernel_run)(d, N, n, m, PRE_FULL_PSI, x, f_hat, f, tr0, \
    ad0); \
\
  for (i = 0; i < sizeof(kernel_modes) / sizeof(kernel_modes[0]); i++) \
    for (v = 0; v < (size_t)(NVARIANTS); v++) \
    { \
      MODULE(kernel_run)(d, N, n, m, kernel_modes[i] | kernel_variants[v], \
        x, f_hat, f, tr, ad); \
      ok &= kernel_report(NAME, d, m, kernel_mode_names[i], \
        kernel_variant_names[v], "trafo", ERR(tr, tr0, (INT)(KERNEL_M))); \
      ok &= kernel_report(NAME, d, m, kernel_mode_names[i], \
        kernel_variant_names[v], "adjoint", ERR(ad, ad0, n_hat)); \
    } \
\
  Y(free)(x); \
  Y(free)(f_hat); \
  Y(free)(f); \
  Y(free)(tr0); \
  Y(free)(tr); \
  Y(free)(ad0); \
  Y(free)(ad); \
  return ok; \
} \
\
void MODULE(check_kernel)(void) \
{ \
  int d, ok = 1; \
  size_t i; \
\
  Y(srand48)(KERNEL_SEED); \
\
  for (d = 1; d <= 3; d++) \
    for (i = 0; i < sizeof(kernel_m) / sizeof(kernel_m[0]); i++) \
      ok &= MODULE(kernel_one)(d, kernel_m[i]); \
\
  for (i = 0; i < sizeof(kernel_m_4d) / sizeof(kernel_m_4d[0]); i++) \
    ok &= MODULE(kernel_one)(4, kernel_m_4d[i]); \
\
  CU_ASSERT(ok); \
}

#undef X
#define X(name) NFFT(name)
KERNEL_CHECK(X, "nfft", C, kernel_err_c, 3, 2, K(-0.5), K(1.0))

#ifdef HAVE_NFCT
#undef X
#define X(name) NFCT(name)
KERNEL_CHECK(X, "nfct", R, kernel_err_r, 2, 1, K(0.0), K(0.5))
#endif

#ifdef HAVE_NFST
#undef X
#define X(name) NFST(name)
KERNEL_CHECK(X, "nfst", R, kernel_err_r, 2, -1, K(0.0), K(0.5))
#endif
