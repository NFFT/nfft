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

/* Tests of the tiles of the grid of B in NFFT, NFCT and NFST.
 *
 * With NFFT_SORT_NODES the nodes are sorted by tile, then by cell, and the
 * threaded B^T spreads the tiles in colour phases (Y(tiles_by_colour)). The
 * grids here are large enough for several tiles per axis (odd and even
 * counts, tiles clamped at the end of an axis, edges clamped to 2m+1), so
 * that tiles of different colours share grid points. For d = 1, ..., 4 and
 * the psi modes, the sorted transforms (NFFT: also with
 * NFFT_OMP_BLOCKWISE_ADJOINT) must agree to rounding with the unsorted ones;
 * the tile ranges must be ordered, cover all nodes and match the sort keys;
 * and the sorted adjoint must give the same bits with one thread as with
 * all of them. The NFCT and NFST nodes lie in [-0.25, 0.75), so that the
 * bucket of nodes outside [0, 0.5] after the tiles is covered too. Only
 * failures are printed (all results with -DTILE_VERBOSE). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "tile.h"

#define TILE_SEED 20261031L

/* agreement of two transforms that sum the same terms in a different order */
#define TILE_TOL (K(1000.0) * EPSILON)

/* d, N per axis, m, number of nodes: tile edges 1024 (d = 1), 64 (d = 2),
 * 16, 12 and 2m+1 = 17 (d = 3, columns over the last axis) and 2m+2 (d = 4) */
static const int tile_cases[][4] = {
  {1, 1100, 3, 4000}, {1, 600, 6, 3000},
  {2, 80, 2, 3000}, {2, 150, 4, 4000},
  {3, 20, 1, 3000}, {3, 16, 4, 3000}, {3, 20, 8, 1500},
  {4, 6, 1, 2000}, {4, 8, 2, 1000}
};

static const unsigned tile_modes[] = {0U, PRE_PSI, PRE_FULL_PSI, PRE_LIN_PSI
#if defined(GAUSSIAN)
  , FG_PSI, FG_PSI | PRE_FG_PSI
#endif
};
static const char *tile_mode_names[] = {"no psi", "PRE_PSI", "FULL_PSI",
  "LIN_PSI"
#if defined(GAUSSIAN)
  , "FG_PSI", "PRE_FG"
#endif
};

#define TILE_FLAGS (PRE_PHI_HUT | MALLOC_X | MALLOC_F_HAT | MALLOC_F \
  | FFTW_INIT | FFT_OUT_OF_PLACE)
#define TILE_FFTW_FLAGS (FFTW_ESTIMATE | FFTW_DESTROY_INPUT)

static R tile_err_c(const C *a, const C *b, const INT n)
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

static R tile_err_r(const R *a, const R *b, const INT n)
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

static int tile_report(const char *module, const int *cs, const char *mode,
  const char *what, const R err, const R tol)
{
  const int ok = (err <= tol);

#ifndef TILE_VERBOSE
  if (ok)
    return ok;
#endif

  printf("%s_tile d = %d, N = %4d, m = %d, %-8s %-22s -> %-4s %.3e (%.3e)\n",
    module, cs[0], cs[1], cs[2], mode, what, ok ? "OK" : "FAIL", (double)err,
    (double)tol);

  return ok;
}

/* Tile ranges of a sorted plan: ordered, from 0 to M over the nt tiles and
 * the bucket after them, and node k in range b exactly when its key lies in
 * [b, b+1) times the tile volume. */
static int tile_ranges_ok(const INT d, const INT L, const INT *ext,
  const INT M, const INT *index_x, const INT *tile_ptr)
{
  INT edge[4], ntiles[4], nt, tilevol = 1, t, b, k;

  nt = Y(tile_layout)(d, L, ext, edge, ntiles);
  for (t = 0; t < d; t++)
  {
    /* an edge of at least 2m+1 keeps the windows of a colour apart */
    if (edge[t] < L - 1 && edge[t] != ext[t])
      return 0;
    tilevol *= edge[t];
  }

  if (tile_ptr[0] != 0 || tile_ptr[nt + 1] != M)
    return 0;

  for (b = 0; b <= nt; b++)
  {
    if (tile_ptr[b] > tile_ptr[b + 1])
      return 0;
    for (k = tile_ptr[b]; k < tile_ptr[b + 1]; k++)
      if (index_x[2 * k] / tilevol != b)
        return 0;
  }

  return 1;
}

/* The check for one module. T is the sample type (C or R), ERR its distance,
 * NVARIANTS the sorted variants (NFFT: also blockwise), OFS = n - 2N per
 * dimension, EXT(n) the cells per axis; the nodes lie in [XLO, XLO + XW). */
#define TILE_CHECK(MODULE, NAME, T, ERR, NVARIANTS, OFS, EXT, XLO, XW) \
static void MODULE(tile_random)(T *v, const INT n) \
{ \
  INT j; \
  for (j = 0; j < n; j++) \
    v[j] = (T)(Y(drand48)() - K(0.5)); \
} \
\
/* trafo and adjoint of a plan with the given flags, the adjoint with \
 * nthreads threads (0: as set); checks the tile ranges of sorted plans; \
 * returns the number of coefficients (NFST: N_t - 1 per dimension) */ \
static INT MODULE(tile_run)(const int *cs, int *N, int *n, \
  const unsigned flags, const R *x, const T *f_hat, const T *f, T *tr, \
  T *ad, const int nthreads, int *ranges_ok) \
{ \
  const int d = cs[0], m = cs[2], M = cs[3]; \
  const INT nt0 = Y(get_num_threads)(); \
  MODULE(plan) p; \
  INT n_hat; \
\
  MODULE(tile_init)(&p, d, N, M, n, m, TILE_FLAGS | flags); \
  memcpy(p.x, x, (size_t)(d * M) * sizeof(R)); \
  MODULE(precompute_one_psi)(&p); \
\
  if (flags & NFFT_SORT_NODES) \
  { \
    INT ext[4], t; \
    for (t = 0; t < d; t++) \
      ext[t] = EXT(p.n[t]); \
    *ranges_ok &= tile_ranges_ok(d, 2 * m + 2, ext, M, p.index_x, \
      p.tile_ptr); \
  } \
\
  memcpy(p.f_hat, f_hat, (size_t)(p.N_total) * sizeof(T)); \
  MODULE(trafo)(&p); \
  memcpy(tr, p.f, (size_t)(M) * sizeof(T)); \
\
  if (nthreads > 0) \
    Y(set_num_threads)(nthreads); \
  memcpy(p.f, f, (size_t)(M) * sizeof(T)); \
  MODULE(adjoint)(&p); \
  memcpy(ad, p.f_hat, (size_t)(p.N_total) * sizeof(T)); \
  Y(set_num_threads)(nt0); \
  n_hat = p.N_total; \
\
  MODULE(finalize)(&p); \
  return n_hat; \
} \
\
static int MODULE(tile_one)(const int *cs) \
{ \
  const int d = cs[0], M = cs[3]; \
  const unsigned variants[] = {NFFT_SORT_NODES, \
    NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT}; \
  const char *variant_names[] = {"SORT", "SORT|BLOCKWISE"}; \
  int N[4], n[4], t, ok = 1, ranges_ok = 1; \
  INT N_total = 1, n_hat; \
  size_t i, v; \
  R *x; \
  T *f_hat, *f, *tr0, *ad0, *tr, *ad, *ad1; \
\
  for (t = 0; t < d; t++) \
  { \
    N[t] = cs[1]; \
    n[t] = 2 * N[t] + (OFS); \
    N_total *= N[t]; \
  } \
\
  x = (R*)Y(malloc)((size_t)(d * M) * sizeof(R)); \
  f_hat = (T*)Y(malloc)((size_t)N_total * sizeof(T)); \
  f = (T*)Y(malloc)((size_t)(M) * sizeof(T)); \
  tr0 = (T*)Y(malloc)((size_t)(M) * sizeof(T)); \
  tr = (T*)Y(malloc)((size_t)(M) * sizeof(T)); \
  ad0 = (T*)Y(malloc)((size_t)N_total * sizeof(T)); \
  ad = (T*)Y(malloc)((size_t)N_total * sizeof(T)); \
  ad1 = (T*)Y(malloc)((size_t)N_total * sizeof(T)); \
\
  for (t = 0; t < d * M; t++) \
    x[t] = (XLO) + (XW) * Y(drand48)(); \
  MODULE(tile_random)(f_hat, N_total); \
  MODULE(tile_random)(f, (INT)(M)); \
\
  for (i = 0; i < sizeof(tile_modes) / sizeof(tile_modes[0]); i++) \
  { \
    n_hat = MODULE(tile_run)(cs, N, n, tile_modes[i], x, f_hat, f, tr0, \
      ad0, 0, &ranges_ok); \
\
    for (v = 0; v < (size_t)(NVARIANTS); v++) \
    { \
      char what[64]; \
\
      MODULE(tile_run)(cs, N, n, tile_modes[i] | variants[v], x, f_hat, f, \
        tr, ad, 0, &ranges_ok); \
      sprintf(what, "%s trafo", variant_names[v]); \
      ok &= tile_report(NAME, cs, tile_mode_names[i], what, \
        ERR(tr, tr0, (INT)(M)), TILE_TOL); \
      sprintf(what, "%s adjoint", variant_names[v]); \
      ok &= tile_report(NAME, cs, tile_mode_names[i], what, \
        ERR(ad, ad0, n_hat), TILE_TOL); \
\
      /* the same bits with one thread */ \
      MODULE(tile_run)(cs, N, n, tile_modes[i] | variants[v], x, f_hat, f, \
        tr, ad1, 1, &ranges_ok); \
      sprintf(what, "%s adjoint 1 thread", variant_names[v]); \
      ok &= tile_report(NAME, cs, tile_mode_names[i], what, \
        memcmp(ad, ad1, (size_t)n_hat * sizeof(T)) == 0 ? K(0.0) \
        : K(1.0), K(0.0)); \
    } \
  } \
\
  ok &= tile_report(NAME, cs, "", "tile ranges", ranges_ok ? K(0.0) \
    : K(1.0), K(0.0)); \
\
  Y(free)(x); \
  Y(free)(f_hat); \
  Y(free)(f); \
  Y(free)(tr0); \
  Y(free)(tr); \
  Y(free)(ad0); \
  Y(free)(ad); \
  Y(free)(ad1); \
  return ok; \
} \
\
void MODULE(check_tile)(void) \
{ \
  int ok = 1; \
  size_t i; \
\
  Y(srand48)(TILE_SEED); \
\
  for (i = 0; i < sizeof(tile_cases) / sizeof(tile_cases[0]); i++) \
    ok &= MODULE(tile_one)(tile_cases[i]); \
\
  CU_ASSERT(ok); \
}

/* PRE_LIN_PSI with a table of K = 2^10 (m+2) + 1, not a multiple of m+2, so
 * that every kernel of B and B^T must step through the table alike */
static void NFFT(tile_init)(NFFT(plan) *p, const int d, int *N, const int M,
  int *n, const int m, const unsigned flags)
{
  if (flags & PRE_LIN_PSI)
    NFFT(init_lin)(p, d, N, M, n, m, (1 << 10) * (m + 2) + 1, flags,
      TILE_FFTW_FLAGS);
  else
    NFFT(init_guru)(p, d, N, M, n, m, flags, TILE_FFTW_FLAGS);
}

#define NFFT_EXT(n) (n)

#undef X
#define X(name) NFFT(name)
TILE_CHECK(X, "nfft", C, tile_err_c, 2, 2, NFFT_EXT, K(-0.5), K(1.0))

#ifdef HAVE_NFCT
static void NFCT(tile_init)(NFCT(plan) *p, const int d, int *N, const int M,
  int *n, const int m, const unsigned flags)
{
  NFCT(init_guru)(p, d, N, M, n, m, flags, TILE_FFTW_FLAGS);
}

/* cells 0, ..., n-1 of the nodes in [0, 0.5] (cell width 1/(2(n-1))) */
#define NFCT_EXT(n) (n)

#undef X
#define X(name) NFCT(name)
TILE_CHECK(X, "nfct", R, tile_err_r, 1, 1, NFCT_EXT, K(-0.25), K(1.0))
#endif

#ifdef HAVE_NFST
static void NFST(tile_init)(NFST(plan) *p, const int d, int *N, const int M,
  int *n, const int m, const unsigned flags)
{
  NFST(init_guru)(p, d, N, M, n, m, flags, TILE_FFTW_FLAGS);
}

/* cells 0, ..., n+1 of the nodes in [0, 0.5] (cell width 1/(2(n+1))) */
#define NFST_EXT(n) ((n) + 2)

#undef X
#define X(name) NFST(name)
TILE_CHECK(X, "nfst", R, tile_err_r, 1, -1, NFST_EXT, K(-0.25), K(1.0))
#endif
