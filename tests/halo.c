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

/* Tests of the halo of the grid of the matrix B in NFFT, NFCT and NFST.
 *
 * B works on the oversampled grid padded with a halo, which holds the periodic
 * (NFFT), even (NFCT) or odd (NFST) continuation of the grid, so that the
 * window of a node never wraps around or is reflected. The tests put the
 * halo to work as hard as possible:
 * - grids barely larger than the window (NFFT: n = 2m+4, so that nearly every
 *   window wraps around; NFCT, NFST: N = 4, so that windows are reflected at
 *   both ends at once), and regular grids,
 * - nodes on and next to the ends of the node interval,
 * - every psi mode, in place and out of place, with and without sorting, and
 *   for NFFT with the blockwise adjoint,
 * - plans that allocate their grids (FFTW_INIT, B works in place on the padded
 *   FFTW grid) and plans whose caller provides plain grids and FFTW plans (no
 *   FFTW_INIT, B works on a scratch grid of the plan).
 * Each fast transform is compared with the direct transform. The tolerance
 * is loose: a wrong halo point gives an error of the order of the samples,
 * far above the window error. The two kinds of plans differ only in the
 * layout of their grids and must agree to rounding. Only failures are
 * printed (all results with -DHALO_VERBOSE). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "halo.h"

#define HALO_SEED 20261029L

/* Fast against direct transform, relative l2 error. */
#if defined(NFFT_SINGLE)
#define HALO_TOL K(1.0E-3)
#else
#define HALO_TOL K(1.0E-6)
#endif
/* PRE_LIN_PSI interpolates the window linearly between table points. */
#define HALO_TOL_LIN K(1.0E-3)
/* d = 4 runs with m = 4: window error 1e-5 to 1e-3, depending on the window */
#define HALO_TOL_4D K(1.0E-2)
#define HALO_M_4D 4

/* Plans with plain and with padded grids, relative l2 distance. */
#define HALO_TOL_LAYOUT (K(500.0) * EPSILON)

#define HALO_M 64

/* Psi modes, as precompute flags; the Gaussian window only for FG_PSI. */
static const unsigned halo_modes[] = {0U, PRE_PSI, PRE_FULL_PSI, PRE_LIN_PSI
#if defined(GAUSSIAN)
  , FG_PSI, FG_PSI | PRE_FG_PSI
#endif
};
static const char *halo_mode_names[] = {"no psi", "PRE_PSI", "PRE_FULL_PSI",
  "PRE_LIN_PSI"
#if defined(GAUSSIAN)
  , "FG_PSI", "PRE_FG_PSI"
#endif
};
#define HALO_NMODES (sizeof(halo_modes) / sizeof(halo_modes[0]))

/* Variants: sorting, FFT in place or out of place, the blockwise adjoint
 * (NFFT only, the last one). */
static const unsigned halo_variants[] = {FFT_OUT_OF_PLACE, 0U,
  NFFT_SORT_NODES | FFT_OUT_OF_PLACE, NFFT_SORT_NODES,
  NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT};
static const char *halo_variant_names[] = {"out of place", "in place",
  "SORT, out of place", "SORT, in place", "SORT|BLOCKWISE"};

#define HALO_FLAGS (PRE_PHI_HUT | MALLOC_X | MALLOC_F_HAT | MALLOC_F)
#define HALO_FFTW_FLAGS (FFTW_ESTIMATE | FFTW_DESTROY_INPUT)

static R halo_err_c(const C *a, const C *b, const INT n)
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

static R halo_err_r(const R *a, const R *b, const INT n)
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

static int halo_report(const char *module, const int d, const int *n,
  const char *variant, const char *mode, const char *grid, const char *which,
  const R err, const R tol)
{
  const int ok = (err <= tol);
  int t;

#ifndef HALO_VERBOSE
  if (ok)
    return ok;
#endif

  printf("%s_halo d = %d, n = [", module, d);
  for (t = 0; t < d; t++)
    printf("%s%3d", t > 0 ? ", " : "", n[t]);
  printf("], %-18s, %-12s, %-6s, %-13s -> %-4s %.3e (%.3e)\n", variant, mode,
    grid, which, ok ? "OK" : "FAIL", (double)err, (double)tol);

  return ok;
}

/* The check for one module. T is the sample type (C or R), ERR its distance,
 * XLO and XHI the node interval [XLO, XHI), NVARIANTS the number of variants
 * to test, PLAIN_GRIDS(p, nn, out_of_place) gives a plan without FFTW_INIT
 * its grids and FFTW plans, FREE_GRIDS(p, out_of_place) frees them. */
#define HALO_CHECK(MODULE, NAME, T, ERR, XLO, XHI, NVARIANTS, PLAIN_GRIDS, \
  FREE_GRIDS, INIT) \
static void MODULE(halo_nodes)(R *x, const int d, const int M) \
{ \
  /* nodes on and next to the ends of the interval, then random ones */ \
  const R edge[] = {XLO, XLO + EPSILON, XHI - K(4.0) * EPSILON, XHI \
    - K(0.25) / (R)(M), XLO + K(0.25) / (R)(M), (XLO + XHI) / K(2.0)}; \
  const int nedge = (int)(sizeof(edge) / sizeof(edge[0])); \
  int j, t; \
\
  for (j = 0; j < M; j++) \
    for (t = 0; t < d; t++) \
    { \
      if (j < nedge * nedge) \
        x[j * d + t] = edge[(t == 0 ? j : j / nedge + t) % nedge]; \
      else \
        x[j * d + t] = XLO + (XHI - XLO) * Y(drand48)(); \
    } \
} \
\
static void MODULE(halo_random)(T *v, const INT n) \
{ \
  INT j; \
  for (j = 0; j < n; j++) \
    v[j] = (T)(Y(drand48)() - K(0.5)); \
} \
\
/* trafo and adjoint of plan p, fast and direct, on the same input; the fast \
 * results are left in tr and ad */ \
static int MODULE(halo_run)(MODULE(plan) *p, const R *x, const T *f_hat, \
  const T *f, T *tr, T *ad, const int d, const int *n, const char *variant, \
  const char *mode, const char *grid, const R tol) \
{ \
  int ok = 1; \
  T *ref_tr = (T*)Y(malloc)((size_t)(p->M_total) * sizeof(T)); \
  T *ref_ad = (T*)Y(malloc)((size_t)(p->N_total) * sizeof(T)); \
\
  memcpy(p->x, x, (size_t)(p->d * p->M_total) * sizeof(R)); \
  MODULE(precompute_one_psi)(p); \
\
  memcpy(p->f_hat, f_hat, (size_t)(p->N_total) * sizeof(T)); \
  MODULE(trafo_direct)(p); \
  memcpy(ref_tr, p->f, (size_t)(p->M_total) * sizeof(T)); \
  MODULE(trafo)(p); \
  memcpy(tr, p->f, (size_t)(p->M_total) * sizeof(T)); \
  ok &= halo_report(NAME, d, n, variant, mode, grid, "trafo", \
    ERR(tr, ref_tr, p->M_total), tol); \
\
  memcpy(p->f, f, (size_t)(p->M_total) * sizeof(T)); \
  MODULE(adjoint_direct)(p); \
  memcpy(ref_ad, p->f_hat, (size_t)(p->N_total) * sizeof(T)); \
  MODULE(adjoint)(p); \
  memcpy(ad, p->f_hat, (size_t)(p->N_total) * sizeof(T)); \
  ok &= halo_report(NAME, d, n, variant, mode, grid, "adjoint", \
    ERR(ad, ref_ad, p->N_total), tol); \
\
  Y(free)(ref_tr); \
  Y(free)(ref_ad); \
  return ok; \
} \
\
static int MODULE(halo_one)(const int d, int *N, int *n, const int m, \
  const unsigned mode, const char *mode_name, const unsigned variant, \
  const char *variant_name, const R tol) \
{ \
  MODULE(plan) p; \
  INT N_total = 1; \
  int t, ok = 1; \
  const int M = HALO_M; \
  R *x; \
  T *f_hat, *f, *tr0, *ad0, *tr1, *ad1; \
\
  for (t = 0; t < d; t++) \
    N_total *= N[t]; \
\
  x = (R*)Y(malloc)((size_t)(d * M) * sizeof(R)); \
  f_hat = (T*)Y(malloc)((size_t)N_total * sizeof(T)); \
  f = (T*)Y(malloc)((size_t)M * sizeof(T)); \
  tr0 = (T*)Y(malloc)((size_t)M * sizeof(T)); \
  tr1 = (T*)Y(malloc)((size_t)M * sizeof(T)); \
  ad0 = (T*)Y(malloc)((size_t)N_total * sizeof(T)); \
  ad1 = (T*)Y(malloc)((size_t)N_total * sizeof(T)); \
  MODULE(halo_nodes)(x, d, M); \
  MODULE(halo_random)(f_hat, N_total); \
  MODULE(halo_random)(f, (INT)M); \
\
  /* the plan allocates its grids: B works on the padded FFTW grid */ \
  INIT(&p, d, N, M, n, m, HALO_FLAGS | FFTW_INIT | mode | variant, \
    HALO_FFTW_FLAGS); \
  ok &= MODULE(halo_run)(&p, x, f_hat, f, tr0, ad0, d, n, variant_name, \
    mode_name, "padded", tol); \
  MODULE(finalize)(&p); \
\
  /* the caller provides plain grids: B works on a scratch grid */ \
  INIT(&p, d, N, M, n, m, HALO_FLAGS | mode | variant, HALO_FFTW_FLAGS); \
  { \
    int nn[4]; \
    for (t = 0; t < d; t++) \
      nn[t] = n[t]; \
    PLAIN_GRIDS(p, nn, variant & FFT_OUT_OF_PLACE); \
  } \
  ok &= MODULE(halo_run)(&p, x, f_hat, f, tr1, ad1, d, n, variant_name, \
    mode_name, "plain", tol); \
  ok &= halo_report(NAME, d, n, variant_name, mode_name, "both", \
    "trafo agrees", ERR(tr1, tr0, (INT)M), HALO_TOL_LAYOUT); \
  ok &= halo_report(NAME, d, n, variant_name, mode_name, "both", \
    "adj. agrees", ERR(ad1, ad0, p.N_total), HALO_TOL_LAYOUT); \
  FREE_GRIDS(p, variant & FFT_OUT_OF_PLACE); \
  MODULE(finalize)(&p); \
\
  Y(free)(x); \
  Y(free)(f_hat); \
  Y(free)(f); \
  Y(free)(tr0); \
  Y(free)(tr1); \
  Y(free)(ad0); \
  Y(free)(ad1); \
  return ok; \
} \
\
/* all modes and variants for one grid */ \
static int MODULE(halo_grid)(const int d, int *N, int *n, const int m, \
  const R tol) \
{ \
  int ok = 1; \
  size_t i, v; \
\
  for (i = 0; i < HALO_NMODES; i++) \
    for (v = 0; v < (size_t)(NVARIANTS); v++) \
      ok &= MODULE(halo_one)(d, N, n, m, halo_modes[i], halo_mode_names[i], \
        halo_variants[v], halo_variant_names[v], \
        (halo_modes[i] & PRE_LIN_PSI) ? FMAX(tol, HALO_TOL_LIN) : tol); \
\
  return ok; \
}

/* NFFT: plain grids of n_total elements, FFTW plans as FFTW_INIT made them
 * before the halo. */
#define NFFT_PLAIN_GRIDS(p, nn, oop) \
{ \
  (p).g1 = Y(malloc)((size_t)((p).n_total) * sizeof(C)); \
  (p).g2 = (oop) ? Y(malloc)((size_t)((p).n_total) * sizeof(C)) \
    : (void*)(p).g1; \
  (p).my_fftw_plan1 = FFTW(plan_dft)((int)(p).d, nn, (p).g1, (p).g2, \
    FFTW_FORWARD, (p).fftw_flags); \
  (p).my_fftw_plan2 = FFTW(plan_dft)((int)(p).d, nn, (p).g2, (p).g1, \
    FFTW_BACKWARD, (p).fftw_flags); \
}
#define NFFT_FREE_GRIDS(p, oop) \
{ \
  FFTW(destroy_plan)((p).my_fftw_plan1); \
  FFTW(destroy_plan)((p).my_fftw_plan2); \
  if (oop) \
    Y(free)((p).g2); \
  Y(free)((p).g1); \
}

/* NFFT: PRE_LIN_PSI with a table of 2^12 (m+2) points per axis; the
 * default of init_guru grows to 2^24 (m+2) for large m. */
#define NFFT_INIT(p, d, N, M, n, m, flags, fftw_flags) \
{ \
  if ((flags) & PRE_LIN_PSI) \
    NFFT(init_lin)(p, d, N, M, n, m, (1 << 12) * ((m) + 2), flags, \
      fftw_flags); \
  else \
    NFFT(init_guru)(p, d, N, M, n, m, flags, fftw_flags); \
}

#undef X
#define X(name) NFFT(name)
HALO_CHECK(X, "nfft", C, halo_err_c, K(-0.5), K(0.5), 5, NFFT_PLAIN_GRIDS,
  NFFT_FREE_GRIDS, NFFT_INIT)

void X(check_halo)(void)
{
  /* smallest even N above m (the fast transforms need N > m), with the
   * tightest grid the fast transforms take, n = 2m+4, and with n = 4N */
  const int m = WINDOW_HELP_ESTIMATE_m;
  const int N0 = 2 * ((m + 2) / 2);
  int d, ok = 1;

  Y(srand48)(HALO_SEED);

  for (d = 1; d <= 3; d++)
  {
    int N[3], n[3], t;
    for (t = 0; t < d; t++)
    {
      N[t] = N0;
      n[t] = 2 * m + 4;
    }
    ok &= X(halo_grid)(d, N, n, m, HALO_TOL);
    for (t = 0; t < d; t++)
      n[t] = 4 * N0;
    ok &= X(halo_grid)(d, N, n, m, HALO_TOL);
  }

  /* generic d */
  {
    int N[4], n[4], t;
    for (t = 0; t < 4; t++)
    {
      N[t] = 6;
      n[t] = 2 * HALO_M_4D + 4;
    }
    ok &= X(halo_grid)(4, N, n, HALO_M_4D, HALO_TOL_4D);
  }

  CU_ASSERT(ok);
}

/* NFCT, NFST: plain grids of n_total elements, r2r FFTW plan. */
#define TRIG_PLAIN_GRIDS(p, nn, oop) \
{ \
  (p).g1 = (R*)Y(malloc)((size_t)((p).n_total) * sizeof(R)); \
  (p).g2 = (oop) ? (R*)Y(malloc)((size_t)((p).n_total) * sizeof(R)) \
    : (p).g1; \
  (p).my_fftw_r2r_plan = FFTW(plan_r2r)((int)(p).d, nn, (p).g1, (p).g2, \
    (p).r2r_kind, (p).fftw_flags); \
}
#define TRIG_FREE_GRIDS(p, oop) \
{ \
  FFTW(destroy_plan)((p).my_fftw_r2r_plan); \
  if (oop) \
    Y(free)((p).g2); \
  Y(free)((p).g1); \
}

/* NFCT and NFST: N = 4, where windows are reflected at both ends at once,
 * and N = 16; sigma = 2. OFS is n - 2N (+1 NFCT, -1 NFST). */
#define TRIG_CHECK_HALO(MODULE, OFS) \
void MODULE(check_halo)(void) \
{ \
  const int m = WINDOW_HELP_ESTIMATE_m; \
  static const int NN[] = {4, 16}; \
  int d, ok = 1; \
  size_t i; \
\
  Y(srand48)(HALO_SEED); \
\
  for (d = 1; d <= 3; d++) \
    for (i = 0; i < sizeof(NN) / sizeof(NN[0]); i++) \
    { \
      int N[3], n[3], t; \
      if (d == 3 && NN[i] > 8) \
        continue; \
      for (t = 0; t < d; t++) \
      { \
        N[t] = NN[i]; \
        n[t] = 2 * NN[i] + (OFS); \
      } \
      ok &= MODULE(halo_grid)(d, N, n, m, HALO_TOL); \
    } \
\
  CU_ASSERT(ok); \
}

#ifdef HAVE_NFCT
#undef X
#define X(name) NFCT(name)
HALO_CHECK(X, "nfct", R, halo_err_r, K(0.0), K(0.5), 4, TRIG_PLAIN_GRIDS,
  TRIG_FREE_GRIDS, X(init_guru))
TRIG_CHECK_HALO(X, 1)
#endif

#ifdef HAVE_NFST
#undef X
#define X(name) NFST(name)
HALO_CHECK(X, "nfst", R, halo_err_r, K(0.0), K(0.5), 4, TRIG_PLAIN_GRIDS,
  TRIG_FREE_GRIDS, X(init_guru))
TRIG_CHECK_HALO(X, -1)
#endif
