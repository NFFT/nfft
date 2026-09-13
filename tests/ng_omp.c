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
#include "ng_omp.h"

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
  CU_ASSERT_EQUAL(Y(nfft_threads_roster_size)(), 1);
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
                           K(1.0e-13) * (CABS(g_ser[i]) + K(1.0)));

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
                           K(1.0e-13) * (CABS(fh_ser[i]) + K(1.0)));

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
