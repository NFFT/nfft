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

#include <complex.h> /* before nfft3.h so fftw_complex is C-compatible */
#include <math.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "deconv.h"
#include "conv.h"
#include "ng_omp.h"
#include "alloc_count.h"

/* A private planner carrying both DECONV rosters, so a case can pick the
 * serial or the threaded solver by setting nthr and the impatience bound
 * rather than by reaching into the registry. nthr < 2 makes the threaded
 * solver decline; PLNR_NO_NONTHREADED with nthr > 1 makes the serial ones
 * decline. */
static planner *both_rosters(int nthr, int force_threaded)
{
  planner *pl = Y(planner_create)();
  pl->nthr = nthr;
  Y(deconv_solvers_register)(pl);
  Y(nfft_threads_conf_standard)(pl);
  if (force_threaded)
    pl->flags.l |= PLNR_NO_NONTHREADED;
  return pl;
}

/* Plan p on pl, wake it, run one direction, destroy it. Returns 0 when no
 * solver applied, so a case can assert applicability separately. */
static int run_once(planner *pl, problem *p, int forward)
{
  plan *pln = Y(planner_mkplan)(pl, p);
  if (!pln)
    return 0;
  Y(plan_awake)(pln, PLNR_AWAKE);
  if (forward)
    pln->adt->apply(pln, p);
  else
    pln->adt->apply_adjoint(pln, p);
  Y(plan_destroy)(pln);
  return 1;
}

/* The add-on library's contract now that the roster carries the threaded
 * DECONV solver: init_threads succeeds, the count can be raised, and
 * cleanup_threads puts it back. */
void Y(check_ng_omp_roster)(void)
{
  CU_ASSERT_EQUAL(Y(nfft_threads_roster_size)(), 2);
  CU_ASSERT_EQUAL(NFFT(init_threads)(), 1);

  NFFT(plan_with_nthreads)(4);
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 4);

  NFFT(cleanup_threads)();
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1);

  /* Re-initialising after a cleanup works; the second raise takes effect. */
  NFFT(plan_with_nthreads)(2);
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 2);
  NFFT(cleanup_threads)();
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 1);
}

/* Forward and adjoint agree with the serial solver of the same rank. For
 * d <= 3 the factor association matches the serial solvers' exactly, so the
 * tolerance is there only for the d >= 4 case against deconv-nd.c. */
static void compare_rank(int d, const INT *N, const INT *n)
{
  const int w = Y(get_window_id)();
  INT Ntot = 1, ntot = 1;
  C *fh, *g_ser, *g_thr, *fh_ser, *fh_thr, *g_in;
  problem *p;
  planner *pl;
  INT i;
  int t;
  R tol = K(1.0e-13);

  /* max(1e-13, 1e4*eps): rank >= 4 goes through a different factor
   * association (see check_ng_omp_deconv_4d), and the fixed floor alone is
   * below float's own epsilon -- same idiom as check_ng_omp_nfft_all_levels. */
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;

  for (t = 0; t < d; t++) {
    Ntot *= N[t];
    ntot *= n[t];
  }

  fh = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  fh_ser = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  fh_thr = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  g_in = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  g_ser = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  g_thr = (C *)Y(malloc)((size_t)ntot * sizeof(C));

  /* Sentinel: the adjoint solvers must write every slot themselves. Without
   * this, an incomplete write would pass by coincidence if both solvers left
   * the same uninitialised bytes behind. */
  for (i = 0; i < Ntot; i++)
    fh_ser[i] = fh_thr[i] = K(-999.0);

  for (i = 0; i < Ntot; i++)
    fh[i] = (R)(i + 1) / (R)Ntot + K(0.25) * II;
  for (i = 0; i < ntot; i++)
    g_in[i] = (R)(ntot - i) / (R)ntot - K(0.5) * II;

  /* forward: f_hat -> g, serial then threaded */
  p = Y(mkproblem_deconv)(d, N, 0, n, 6, w, 1, fh, g_ser);
  pl = both_rosters(1, 0);
  CU_ASSERT_EQUAL_FATAL(run_once(pl, p, 1), 1);
  Y(planner_destroy)(pl);
  Y(problem_destroy)(p);

  p = Y(mkproblem_deconv)(d, N, 0, n, 6, w, 1, fh, g_thr);
  pl = both_rosters(4, 1);
  CU_ASSERT_EQUAL_FATAL(run_once(pl, p, 1), 1);
  Y(planner_destroy)(pl);
  Y(problem_destroy)(p);

  for (i = 0; i < ntot; i++)
    CU_ASSERT_DOUBLE_EQUAL(CABS(g_thr[i] - g_ser[i]), K(0.0),
                           tol * (CABS(g_ser[i]) + K(1.0)));

  /* adjoint: g -> f_hat, serial then threaded */
  memcpy(g_ser, g_in, (size_t)ntot * sizeof(C));
  p = Y(mkproblem_deconv)(d, N, 0, n, 6, w, -1, fh_ser, g_ser);
  pl = both_rosters(1, 0);
  CU_ASSERT_EQUAL_FATAL(run_once(pl, p, 0), 1);
  Y(planner_destroy)(pl);
  Y(problem_destroy)(p);

  memcpy(g_thr, g_in, (size_t)ntot * sizeof(C));
  p = Y(mkproblem_deconv)(d, N, 0, n, 6, w, -1, fh_thr, g_thr);
  pl = both_rosters(4, 1);
  CU_ASSERT_EQUAL_FATAL(run_once(pl, p, 0), 1);
  Y(planner_destroy)(pl);
  Y(problem_destroy)(p);

  for (i = 0; i < Ntot; i++)
    CU_ASSERT_DOUBLE_EQUAL(CABS(fh_thr[i] - fh_ser[i]), K(0.0),
                           tol * (CABS(fh_ser[i]) + K(1.0)));

  Y(free)(g_thr);
  Y(free)(g_ser);
  Y(free)(g_in);
  Y(free)(fh_thr);
  Y(free)(fh_ser);
  Y(free)(fh);
}

void Y(check_ng_omp_deconv_1d)(void)
{
  const INT N[1] = {64}, n[1] = {128};
  compare_rank(1, N, n);
}

void Y(check_ng_omp_deconv_2d)(void)
{
  const INT N[2] = {16, 24}, n[2] = {32, 48};
  compare_rank(2, N, n);
}

void Y(check_ng_omp_deconv_3d)(void)
{
  const INT N[3] = {8, 10, 12}, n[3] = {16, 20, 24};
  compare_rank(3, N, n);
}

/* Rank 4 goes against deconv-nd.c, whose factor association differs, so this
 * is the case the 1e-13 tolerance exists for. */
void Y(check_ng_omp_deconv_4d)(void)
{
  const INT N[4] = {6, 6, 8, 8}, n[4] = {12, 12, 16, 16};
  compare_rank(4, N, n);
}

/* An odd N and a type-II axis both move Nneg, which is what decides where a
 * row's two runs fall. Covered here rather than in compare_rank because the
 * variant array has to be passed through. */
void Y(check_ng_omp_deconv_variants)(void)
{
  const INT N[2] = {15, 16}, n[2] = {32, 32};
  const int variant[2] = {NFFT_NDFT_TYPE_I, NFFT_NDFT_TYPE_II};
  const int w = Y(get_window_id)();
  const INT Ntot = 15 * 16, ntot = 32 * 32;
  C *fh = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *g_ser = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  C *g_thr = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  problem *p;
  planner *pl;
  INT i;
  R tol = K(1.0e-13);

  /* max(1e-13, 1e4*eps): same idiom as check_ng_omp_nfft_all_levels, kept as
   * an absolute bound since this comparison already is one. */
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;

  for (i = 0; i < Ntot; i++)
    fh[i] = (R)(i + 1) / (R)Ntot;

  p = Y(mkproblem_deconv)(2, N, variant, n, 6, w, 1, fh, g_ser);
  pl = both_rosters(1, 0);
  CU_ASSERT_EQUAL_FATAL(run_once(pl, p, 1), 1);
  Y(planner_destroy)(pl);
  Y(problem_destroy)(p);

  p = Y(mkproblem_deconv)(2, N, variant, n, 6, w, 1, fh, g_thr);
  pl = both_rosters(4, 1);
  CU_ASSERT_EQUAL_FATAL(run_once(pl, p, 1), 1);
  Y(planner_destroy)(pl);
  Y(problem_destroy)(p);

  for (i = 0; i < ntot; i++)
    CU_ASSERT_DOUBLE_EQUAL(CABS(g_thr[i] - g_ser[i]), K(0.0), tol);

  Y(free)(g_thr);
  Y(free)(g_ser);
  Y(free)(fh);
}

/* Every output cell is written by exactly one thread, so the results must not
 * move with the team size. A reduction or an atomic would show up here. */
void Y(check_ng_omp_thread_count_invariance)(void)
{
  const INT N[3] = {8, 10, 12}, n[3] = {16, 20, 24};
  const int w = Y(get_window_id)();
  const INT Ntot = 8 * 10 * 12, ntot = 16 * 20 * 24;
  const int counts[3] = {2, 4, 7};
  C *fh = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *ref = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  C *got = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  INT i;
  int c;

  for (i = 0; i < Ntot; i++)
    fh[i] = (R)(i + 1) / (R)Ntot - K(0.3) * II;

  for (c = 0; c < 3; c++) {
    problem *p =
         Y(mkproblem_deconv)(3, N, 0, n, 6, w, 1, fh, c == 0 ? ref : got);
    planner *pl = both_rosters(counts[c], 1);
    CU_ASSERT_EQUAL_FATAL(run_once(pl, p, 1), 1);
    Y(planner_destroy)(pl);
    Y(problem_destroy)(p);
    if (c > 0)
      for (i = 0; i < ntot; i++)
        CU_ASSERT_EQUAL(memcmp(&got[i], &ref[i], sizeof(C)), 0);
  }

  Y(free)(got);
  Y(free)(ref);
  Y(free)(fh);
}

/* A one-thread budget has nothing to distribute, so the threaded solver must
 * decline and leave the plan to the serial roster. */
void Y(check_ng_omp_declines_serial_budget)(void)
{
  const INT N = 64, n = 128;
  const int w = Y(get_window_id)();
  C *fh = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *g = (C *)Y(malloc)((size_t)n * sizeof(C));
  planner *pl = Y(planner_create)();
  problem *p;
  INT i;

  for (i = 0; i < N; i++)
    fh[i] = K(1.0);

  /* Only the threaded roster, and a one-thread budget: nothing applies. */
  pl->nthr = 1;
  Y(nfft_threads_conf_standard)(pl);
  p = Y(mkproblem_deconv)(1, &N, 0, &n, 6, w, 1, fh, g);
  CU_ASSERT_PTR_NULL(Y(planner_mkplan)(pl, p));
  Y(problem_destroy)(p);
  Y(planner_destroy)(pl);

  /* Two threads, same roster: it applies. */
  pl = Y(planner_create)();
  pl->nthr = 2;
  Y(nfft_threads_conf_standard)(pl);
  p = Y(mkproblem_deconv)(1, &N, 0, &n, 6, w, 1, fh, g);
  CU_ASSERT_EQUAL(run_once(pl, p, 1), 1);
  Y(problem_destroy)(p);
  Y(planner_destroy)(pl);

  Y(free)(g);
  Y(free)(fh);
}

/* Both CONV rosters on one planner, with the same nthr / force_threaded
 * convention as both_rosters. */
static planner *both_conv_rosters(int nthr, int force_threaded)
{
  planner *pl = Y(planner_create)();
  pl->nthr = nthr;
  Y(conv_solvers_register)(pl);
  Y(nfft_threads_conf_standard)(pl);
  if (force_threaded)
    pl->flags.l |= PLNR_NO_NONTHREADED;
  return pl;
}

/* Tables live only between a wake and the return to SLEEPY (ADR-0009), so a
 * measured race holds one candidate's tables at a time. Every CONV candidate
 * this planner has, serial and threaded, races here. */
void Y(check_ng_omp_race_memory)(void)
{
  const INT N[2] = {128, 128}, n[2] = {256, 256}, M = 65536;
  const int w = Y(get_window_id)(), d = 2;
  const size_t len = 14;
  /* One candidate's M-sized tables at most: psi, u, the threaded solver's
   * border, and a sorted solver's perm. */
  const size_t tables = (size_t)M
                        * ((size_t)d * len * sizeof(R)
                           + (size_t)d * sizeof(INT) + 2 * sizeof(INT));
  R *x = (R *)malloc((size_t)(M * d) * sizeof(R));
  C *g = (C *)malloc((size_t)(n[0] * n[1]) * sizeof(C));
  C *f = (C *)malloc((size_t)M * sizeof(C));
  size_t peak;
  unsigned nmeasured;
  int planned;
  INT i;

  for (i = 0; i < M * d; i++)
    x[i] = (R)((i * 7919) % 1000) / K(1000.0) - K(0.5);
  for (i = 0; i < n[0] * n[1]; i++)
    g[i] = K(1.0);
  for (i = 0; i < M; i++)
    f[i] = K(0.0);

  alloc_count_begin();
  {
    planner *pl = both_conv_rosters(4, 0);
    problem *p = Y(mkproblem_conv)(d, n, N, M, 6, w, 1, 0, x, g, f);
    plan *pln;
    pl->nmeasured = 0;
    pln = Y(planner_mkplan)(pl, p);
    planned = pln != 0;
    nmeasured = pl->nmeasured;
    if (pln)
      Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
    Y(planner_destroy)(pl);
  }
  peak = alloc_count_peak;
  alloc_count_end();

  CU_ASSERT_TRUE(planned);
  CU_ASSERT(nmeasured >= 2u); /* a real race */
  CU_ASSERT(peak < tables + tables / 2);

  free(f);
  free(g);
  free(x);
}

/* Plan and run one CONV family: registration id 0 is caller order, 1 is
 * sorted. A planner here offers only serial or only threaded CONV candidates,
 * so the id picks exactly one. Returns 0 when that family declines. A race
 * between families would make the memcmp below depend on which one won. */
static int run_family(planner *pl, problem *p, int forward, int sorted)
{
  plan *cands[8];
  unsigned slv[8];
  plan *pln = 0;
  int nc = Y(planner_candidates)(pl, p, cands, slv, 8), i;
  for (i = 0; i < nc; i++) {
    if (pln == 0 && pl->slvdescs[slv[i]].reg_id == sorted)
      pln = cands[i];
    else
      Y(plan_destroy)(cands[i]);
  }
  if (!pln)
    return 0;
  Y(plan_awake)(pln, PLNR_AWAKE);
  if (forward)
    pln->adt->apply(pln, p);
  else
    pln->adt->apply_adjoint(pln, p);
  Y(plan_destroy)(pln);
  return 1;
}

/* The threaded CONV solver reproduces the serial one to the last bit: the same
 * kernels, the same per-node order. Anything less means an operation moved.
 * counts covers a team larger than the band count and one that is not a
 * divisor of it. clustered picks the node fill: spread across the whole axis,
 * or packed into a 2% slice of it. */
static void compare_conv_body(int d, const INT *N, const INT *n, INT M,
                              int forward, int clustered)
{
  const int w = Y(get_window_id)();
  const int counts[3] = {2, 4, 7};
  INT ntot = 1, i;
  R *x;
  C *g_ref, *g_got, *f_ref, *f_got;
  int t, c;

  for (t = 0; t < d; t++)
    ntot *= n[t];

  x = (R *)Y(malloc)((size_t)M * (size_t)d * sizeof(R));
  g_ref = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  g_got = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  f_ref = (C *)Y(malloc)((size_t)M * sizeof(C));
  f_got = (C *)Y(malloc)((size_t)M * sizeof(C));

  if (clustered)
    for (i = 0; i < M * d; i++)
      x[i] = K(0.02) * ((R)((i * 7919) % 1000) / K(1000.0)) - K(0.49);
  else
    for (i = 0; i < M * d; i++)
      x[i] = (R)((i * 7919) % 1000) / K(1000.0) - K(0.5);
  for (i = 0; i < ntot; i++)
    g_ref[i] = (R)(i + 1) / (R)ntot - K(0.25) * II;
  for (i = 0; i < M; i++)
    f_ref[i] = (R)(M - i) / (R)M + K(0.75) * II;
  memcpy(g_got, g_ref, (size_t)ntot * sizeof(C));
  memcpy(f_got, f_ref, (size_t)M * sizeof(C));

  {
    problem *p = Y(mkproblem_conv)(d, n, N, M, 6, w, 1, 0, x, g_ref, f_ref);
    planner *pl = both_conv_rosters(1, 0);
    CU_ASSERT_EQUAL_FATAL(run_family(pl, p, forward, 0), 1);
    Y(planner_destroy)(pl);
    Y(problem_destroy)(p);
  }

  for (c = 0; c < 3; c++) {
    problem *p;
    planner *pl;
    memcpy(g_got, g_ref, (size_t)ntot * sizeof(C));
    if (!forward)
      for (i = 0; i < ntot; i++)
        g_got[i] = K(0.0);
    p = Y(mkproblem_conv)(d, n, N, M, 6, w, 1, 0, x, g_got, f_got);
    pl = both_conv_rosters(counts[c], 1);
    CU_ASSERT_EQUAL_FATAL(run_family(pl, p, forward, 0), 1);
    Y(planner_destroy)(pl);
    Y(problem_destroy)(p);

    if (forward) {
      CU_ASSERT_EQUAL(memcmp(f_got, f_ref, (size_t)M * sizeof(C)), 0);
    } else {
      CU_ASSERT_EQUAL(memcmp(g_got, g_ref, (size_t)ntot * sizeof(C)), 0);
    }
  }

  Y(free)(f_got);
  Y(free)(f_ref);
  Y(free)(g_got);
  Y(free)(g_ref);
  Y(free)(x);
}

static void compare_conv(int d, const INT *N, const INT *n, INT M, int forward)
{
  compare_conv_body(d, N, n, M, forward, 0);
}

static void compare_conv_clustered(int d, const INT *N, const INT *n, INT M,
                                   int forward)
{
  compare_conv_body(d, N, n, M, forward, 1);
}

void Y(check_ng_omp_conv_forward)(void)
{
  {
    const INT N[1] = {64}, n[1] = {128};
    compare_conv(1, N, n, 200, 1);
  }
  {
    const INT N[2] = {16, 24}, n[2] = {48, 64};
    compare_conv(2, N, n, 150, 1);
  }
  {
    const INT N[3] = {8, 10, 12}, n[3] = {32, 32, 32};
    compare_conv(3, N, n, 120, 1);
  }
  {
    const INT N[4] = {6, 6, 8, 8}, n[4] = {24, 24, 24, 24};
    compare_conv(4, N, n, 60, 1);
  }
}

void Y(check_ng_omp_conv_adjoint)(void)
{
  {
    const INT N[1] = {64}, n[1] = {128};
    compare_conv(1, N, n, 200, 0);
  }
  {
    const INT N[2] = {16, 24}, n[2] = {48, 64};
    compare_conv(2, N, n, 150, 0);
  }
  {
    const INT N[3] = {8, 10, 12}, n[3] = {32, 32, 32};
    compare_conv(3, N, n, 120, 0);
  }
  {
    /* n[0] = 48 against a support of 2*6+2 = 14 rows gives nb = 3 bands, so
     * the band clip in Y(conv_node_nd) and the adjoint's two-list merge are
     * actually exercised for rank >= 4 (nb was 1, a no-op, at n[0] = 24). */
    const INT N[4] = {6, 6, 8, 8}, n[4] = {48, 24, 24, 24};
    compare_conv(4, N, n, 60, 0);
  }
}

/* n[0] = 20 against a support of 2*6+2 = 14 rows leaves room for one band, so
 * the adjoint runs on one thread while the forward still uses the team. The
 * result must be the serial one either way. With nb == 1 the adjoint's own
 * `ia = ea` special case skips the cyclic-predecessor list entirely, so what
 * this actually exercises is Y(conv_runs)'s own axis wrap, not the
 * cyclic-predecessor merge -- check_ng_omp_conv_wrap covers that. */
void Y(check_ng_omp_conv_narrow_grid)(void)
{
  const INT N[2] = {8, 8}, n[2] = {20, 20};
  compare_conv(2, N, n, 80, 0);
  compare_conv(2, N, n, 80, 1);
}

/* Nodes clustered into a tenth of the grid: most bands hold nothing and one
 * holds nearly everything. This is where an off-by-one in the merge shows,
 * since a band's predecessor list is then either empty or the whole set. It is
 * a correctness case only -- the balance claim it motivates is a performance
 * claim and belongs in the ADR's measured table. */
void Y(check_ng_omp_conv_clustered)(void)
{
  const INT N[2] = {16, 16}, n[2] = {64, 64};
  compare_conv_clustered(2, N, n, 400, 0);
  compare_conv_clustered(2, N, n, 400, 1);
}

/* A band count that does not divide the grid, and nodes crowded at the wrap,
 * so that the last band is wider than the rest and band 0's predecessor is the
 * last one. */
void Y(check_ng_omp_conv_wrap)(void)
{
  const INT N[1] = {16}, n[1] = {45};
  compare_conv(1, N, n, 300, 0);
}

/* With both rosters live and two threads, a DECONV child problem has two
 * applicable solvers, so both are measured and the winner carries a measured
 * cost. This is the race the whole branch exists to enable. */
void Y(check_ng_omp_child_race)(void)
{
  const INT N = 4096, n = 8192;
  const int w = Y(get_window_id)();
  C *fh = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *g = (C *)Y(malloc)((size_t)n * sizeof(C));
  planner *pl = both_rosters(4, 0);
  problem *p;
  plan *pln;
  INT i;

  for (i = 0; i < N; i++)
    fh[i] = K(1.0);
  for (i = 0; i < n; i++)
    g[i] = K(0.0);

  pl->nmeasured = 0;
  p = Y(mkproblem_deconv)(1, &N, 0, &n, 6, w, 1, fh, g);
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  CU_ASSERT_EQUAL(pl->nmeasured, 2u);
  CU_ASSERT_TRUE(pln->pcost > 0.0);
  Y(plan_destroy)(pln);
  Y(problem_destroy)(p);
  Y(planner_destroy)(pl);

  Y(free)(g);
  Y(free)(fh);
}

/* Whether p's printed plan tree contains needle. */
static int plan_tree_has(Y(plan_ng) *p, const char *needle)
{
  FILE *tmp = tmpfile();
  long len;
  char *buf;
  int found;

  if (!tmp)
    return 0;
  NFFT(fprint_plan)(p, tmp);
  fseek(tmp, 0, SEEK_END);
  len = ftell(tmp);
  if (len < 0)
    len = 0;
  rewind(tmp);
  buf = (char *)Y(malloc)((size_t)len + 1);
  if (len > 0 && fread(buf, 1, (size_t)len, tmp) != (size_t)len)
    len = 0;
  buf[len] = '\0';
  found = strstr(buf, needle) != 0;
  Y(free)(buf);
  fclose(tmp);
  return found;
}

/* End to end through the public API with the count raised, at every patience
 * level. Below NFFT_PATIENT the mapping sets PLNR_NO_NONTHREADED, so both
 * children must be the threaded ones and fast_native must not decline. The
 * results must match the one-thread results. */
void Y(check_ng_omp_nfft_all_levels)(void)
{
  const INT N = 4096, n = 8192, M = 512;
  const unsigned levels[3] = {NFFT_ESTIMATE, NFFT_MEASURE, NFFT_PATIENT};
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f_ref = (C *)Y(malloc)((size_t)M * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  INT j;
  int i;
  R tol;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);

  /* Reference on one thread. */
  {
    Y(plan_ng) *p =
         NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT(get_window_id)(), x,
                              (FC *)f_hat, (FC *)f_ref, 0u, NFFT_ESTIMATE);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);
    for (j = 0; j < N; j++)
      f_hat[j] = (R)(j + 1) / (R)N;
    NFFT(precompute)(p);
    NFFT(execute)(p);
    NFFT(plan_ng_destroy)(p);
  }

  NFFT(plan_with_nthreads)(4);
  CU_ASSERT_EQUAL(NFFT(planner_nthreads)(), 4);

  /* max(1e-12, 1e4*eps): trafo vs. a genuinely different code path (threaded
   * composer vs. the nthr=1 reference), so real precision-scale rounding is
   * expected; the fixed floor alone is below float's own epsilon. */
  tol = K(1.0e-12);
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;

  for (i = 0; i < 3; i++) {
    Y(plan_ng) *p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6,
                                       NFFT(get_window_id)(), x,
                              (FC *)f_hat, (FC *)f, 0u, levels[i]);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);
    /* Below PATIENT, NO_NONTHREADED forces every serial DECONV/CONV solver to
     * decline once nthr > 1, so both threaded solvers are structurally
     * guaranteed to be in the tree -- this is checked below only for
     * NFFT_ESTIMATE and NFFT_MEASURE. At NFFT_PATIENT the serial children
     * compete again and Task 7's race picks a winner by real measured cost;
     * either may legitimately win, so only correctness is checked there, not
     * which solver answered. */
    if (levels[i] != NFFT_PATIENT) {
      CU_ASSERT_TRUE(plan_tree_has(p, "deconv_solver_threaded"));
      CU_ASSERT_TRUE(plan_tree_has(p, "conv_solver_threaded"));
    }
    for (j = 0; j < N; j++)
      f_hat[j] = (R)(j + 1) / (R)N;
    NFFT(precompute)(p);
    NFFT(execute)(p);
    for (j = 0; j < M; j++)
      CU_ASSERT_DOUBLE_EQUAL(CABS(f[j] - f_ref[j]), K(0.0),
                             tol * (CABS(f_ref[j]) + K(1.0)));
    NFFT(execute_adjoint)(p);
    NFFT(plan_ng_destroy)(p);
  }

  NFFT(cleanup_threads)();
  Y(free)(f);
  Y(free)(f_ref);
  Y(free)(f_hat);
  Y(free)(x);
}

/* NFFT_NO_FAST_NATIVE forbids the only composer, and the direct NDFT solvers
 * genuinely own their loops, so below NFFT_PATIENT a raised count has nothing
 * left to plan. Pinned so the refusal stays deliberate. */
void Y(check_ng_omp_no_fast_native_refuses)(void)
{
  const INT N = 64, n = 128, M = 100;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  INT j;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);

  NFFT(plan_with_nthreads)(4);
  CU_ASSERT_PTR_NULL(NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6,
                                        NFFT(get_window_id)(), x,
                                  (FC *)f_hat, (FC *)f, 0u,
                                  NFFT_MEASURE | NFFT_NO_FAST_NATIVE));
  /* At NFFT_PATIENT the serial NDFT competes again and it plans. */
  {
    Y(plan_ng) *p = NFFT(plan_ng_guru)(1, &N, 0, &n, M, 6,
                                       NFFT(get_window_id)(), x,
                              (FC *)f_hat, (FC *)f, 0u,
                              NFFT_PATIENT | NFFT_NO_FAST_NATIVE);
    CU_ASSERT_PTR_NOT_NULL(p);
    NFFT(plan_ng_destroy)(p);
  }
  NFFT(cleanup_threads)();
  Y(free)(f);
  Y(free)(f_hat);
  Y(free)(x);
}
