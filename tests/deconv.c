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

/* Unit tests for the DECONV leaf solver (kernel/deconv/): Step A of the fast
 * NFFT decomposition, the division by phi_hut onto the zero-padded oversampled
 * grid.  Outputs are asserted as values from a clean input -- never as a
 * round-trip, which the inverse would also pass.  The per-rank tests sweep the
 * axis-type space (even type-I, even type-II, odd), because the frequency
 * layout, and therefore the index mapping under test, differs per type. */

#include <complex.h> /* before nfft3.h so fftw_complex is C-compatible */
#include <stdio.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "config.h" /* ABS_SRCDIR */
#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "deconv.h"

void Y(check_deconv_problem)(void) {
  planner *pl = Y(planner_create)();
  INT N1 = 16, n1 = 32;
  int w = Y(get_window_id)();
  problem *pf, *pa, *p2;
  md5sig s1, s2;
  char buf[160];
  printer *pr;
  C fh = K(1.0), g = K(0.0), fh2 = K(2.0);

  /* forward: sz oriented N -> n, no M/x in the type; variant NULL = all type-I */
  pf = Y(mkproblem_deconv)(1, &N1, 0, &n1, 6, w, 1, &fh, &g);
  CU_ASSERT_EQUAL(pf->adt->kind, NFFT_PROBLEM_DECONV);
  CU_ASSERT_EQUAL(Y(problem_deconv_N)(pf, 0), (INT)16);
  CU_ASSERT_EQUAL(Y(problem_deconv_n)(pf, 0), (INT)32);
  CU_ASSERT_EQUAL(Y(problem_deconv_Ntot)(pf), (INT)16);
  CU_ASSERT_EQUAL(Y(problem_deconv_ntot)(pf), (INT)32);
  CU_ASSERT_EQUAL(Y(problem_deconv_variant)(pf, 0), NFFT_NDFT_TYPE_I);

  pr = Y(printer_create_str)(buf);
  pr->print(pr, "%P", pf);
  Y(printer_destroy)
  (pr);
  CU_ASSERT_STRING_EQUAL(buf,
                         "(deconv sign=1 m=6 (tensor 1 (16 1 32 1)) variant=0)");

  /* adjoint: tensor_adjoint orientation, distinct key from forward */
  pa = Y(mkproblem_deconv)(1, &N1, 0, &n1, 6, w, -1, &fh, &g);
  CU_ASSERT_EQUAL(Y(problem_deconv_N)(pa, 0), (INT)16); /* direction-aware */
  CU_ASSERT_EQUAL(Y(problem_deconv_n)(pa, 0), (INT)32);

  Y(problem_md5)
  (pl, pf, s1);
  Y(problem_md5)
  (pl, pa, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);

  /* geometry and m both shift the key */
  p2 = Y(mkproblem_deconv)(1, &N1, 0, &n1, 8, w, 1, &fh, &g);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);
  Y(problem_destroy)
  (p2);

  /* key is data-blind: a different f_hat/g pointer yields the same key */
  p2 = Y(mkproblem_deconv)(1, &N1, 0, &n1, 6, w, 1, &fh2, &g);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2] && s1[3] == s2[3]);
  Y(problem_destroy)
  (p2);

  /* an even-N type-II axis shifts the key; odd N normalizes type-II -> type-I */
  {
    int v2 = NFFT_NDFT_TYPE_II;
    INT No = 15, no = 32;
    problem *pii = Y(mkproblem_deconv)(1, &N1, &v2, &n1, 6, w, 1, &fh, &g);
    problem *poii = Y(mkproblem_deconv)(1, &No, &v2, &no, 6, w, 1, &fh, &g);
    Y(problem_md5)
    (pl, pii, s2);
    CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);
    CU_ASSERT_EQUAL(Y(problem_deconv_variant)(pii, 0), NFFT_NDFT_TYPE_II);
    CU_ASSERT_EQUAL(Y(problem_deconv_variant)(poii, 0), NFFT_NDFT_TYPE_I);
    Y(problem_destroy)
    (pii);
    Y(problem_destroy)
    (poii);
  }

  Y(problem_destroy)
  (pf);
  Y(problem_destroy)
  (pa);
  Y(planner_destroy)
  (pl);
}

/* The DECONV solver, planned directly through planner_mkplan (there is no
 * deconv guru). Forward, adjoint and the type-II shift are each asserted as
 * values from a clean input, never as a round-trip: a round-trip also passes
 * for the inverse, which the adjoint of a real diagonal is not.
 * Uses the process-global planner, so it tears it down at the end. */
void Y(check_deconv_solver)(void) {
  INT N = 16, n = 32;
  int m = 6, w = Y(get_window_id)();
  INT ks;
  C *f_hat, *g;
  problem *pf, *pa;
  plan *pln, *pln_adj;
  R phi_hut0;

  f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  g = (C *)Y(malloc)((size_t)n * sizeof(C));

  Y(deconv_ensure_registered)
  ();

  /* forward: single spike at ks = N/2 (frequency k = 0) */
  for (ks = 0; ks < N; ks++)
    f_hat[ks] = K(0.0);
  f_hat[N / 2] = K(1.0);
  for (ks = 0; ks < n; ks++)
    g[ks] = K(0.0);

  pf = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, 1, f_hat, g);
  pln = Y(planner_mkplan)(Y(the_planner)(), pf);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  Y(plan_awake)
  (pln, PLNR_AWAKE);
  pln->adt->apply(pln, pf);

  phi_hut0 = Y(window_phi_hut)(w, n, N, m, (INT)0);
  CU_ASSERT(CABS(g[0] - K(1.0) / phi_hut0) < K(1e-12));
  for (ks = N; ks < n; ks++)
    CU_ASSERT(CABS(g[ks]) < K(1e-12));

  Y(plan_destroy)
  (pln);
  Y(problem_destroy)
  (pf);

  /* D is a real diagonal scale-and-pad, so D^H multiplies by the same
   * 1/phi_hut and gathers. */
  pa = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, -1, f_hat, g);
  pln_adj = Y(planner_mkplan)(Y(the_planner)(), pa);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln_adj);
  Y(plan_awake)
  (pln_adj, PLNR_AWAKE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln_adj->adt->apply_adjoint);

  /* frequency k=0 lives at pos 0: g[0]=1 -> f_hat[N/2] = 1/phi_hut(0). */
  for (ks = 0; ks < N; ks++)
    f_hat[ks] = K(0.0);
  for (ks = 0; ks < n; ks++)
    g[ks] = K(0.0);
  g[0] = K(1.0);
  pln_adj->adt->apply_adjoint(pln_adj, pa);
  CU_ASSERT(CABS(f_hat[N / 2] - K(1.0) / phi_hut0) < K(1e-12));

  /* frequency k=1 lives at pos 1: g[1]=1 -> f_hat[N/2 + 1] = 1/phi_hut(1).
   * The inverse would give phi_hut(1) here and fail. */
  for (ks = 0; ks < N; ks++)
    f_hat[ks] = K(0.0);
  for (ks = 0; ks < n; ks++)
    g[ks] = K(0.0);
  g[1] = K(1.0);
  pln_adj->adt->apply_adjoint(pln_adj, pa);
  {
    R phi_hut1 = Y(window_phi_hut)(w, n, N, m, (INT)1);
    CU_ASSERT(CABS(f_hat[N / 2 + 1] - K(1.0) / phi_hut1) < K(1e-12));
  }

  Y(plan_destroy)
  (pln_adj);
  Y(problem_destroy)
  (pa);

  /* Type-II: the +1 frequency shift must reach the precomputed window
   * envelope, not just pos. A spike at ks=N/2 is frequency 1, landing at pos
   * 1, so the output divides by phi_hut(1). An awake that drops the shift
   * would divide by phi_hut(0) and fail. */
  {
    int v2 = NFFT_NDFT_TYPE_II;
    problem *p2;
    plan *pln2;
    R phi_hut1 = Y(window_phi_hut)(w, n, N, m, (INT)1);
    for (ks = 0; ks < N; ks++)
      f_hat[ks] = K(0.0);
    f_hat[N / 2] = K(1.0);
    for (ks = 0; ks < n; ks++)
      g[ks] = K(0.0);
    p2 = Y(mkproblem_deconv)(1, &N, &v2, &n, m, w, 1, f_hat, g);
    pln2 = Y(planner_mkplan)(Y(the_planner)(), p2);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln2);
    Y(plan_awake)
    (pln2, PLNR_AWAKE);
    pln2->adt->apply(pln2, p2);
    /* Relative tolerance: 1/phi_hut(0) and 1/phi_hut(1) are both ~7e-12 and
     * differ by only ~2.4%, so an absolute 1e-12 test cannot tell them apart. */
    CU_ASSERT(CABS(g[1] - K(1.0) / phi_hut1) < K(1e-9) * (K(1.0) / phi_hut1));
    CU_ASSERT(CABS(g[0]) < K(1e-9) * (K(1.0) / phi_hut1)); /* pos 0 stays empty */
    Y(plan_destroy)
    (pln2);
    Y(problem_destroy)
    (p2);
  }

  Y(free)
  (f_hat);
  Y(free)
  (g);

  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}

/* Plan a DECONV problem through the process-global planner and wake it. */
static plan *deconv_awake_plan(problem *p)
{
  plan *pln = Y(planner_mkplan)(Y(the_planner)(), p);
  if (pln)
    Y(plan_awake)
  (pln, PLNR_AWAKE);
  return pln;
}

/* The value a unit spike at frequency k must produce: prod_t 1/phi_hut(k[t]). */
static R deconv_scale(int w, int d, const INT *N, const INT *n, int m,
                      const INT *k)
{
  R s = K(1.0);
  int t;
  for (t = 0; t < d; t++)
    s /= Y(window_phi_hut)(w, n[t], N[t], m, k[t]);
  return s;
}

/* The DECONV spike tests compare against a computed scale over a few chained
 * multiplies, so the tolerance scales with eps. Observed float relerr peaks
 * near 1 eps; 64x leaves six binary orders of headroom. */
static R deconv_general_tol(void)
{
  return K(64.0) * Y(float_property)(NFFT_EPSILON);
}

/* Cells whose magnitude exceeds eps. A spike must produce exactly one. */
static INT count_above(const C *v, INT len, R eps)
{
  INT j, c = 0;
  for (j = 0; j < len; j++)
    if (CABS(v[j]) > eps)
      c++;
  return c;
}

/* Rank-1 DECONV for odd N and for type-II: the slot -> grid map is
 * pos(ks) = ks < Nneg ? n - Nneg + ks : ks - Nneg, with Nneg = N/2 - type_ii.
 * The grid is pre-dirtied with the expected peak so a stale zero-pad cell shows
 * up as a second above-threshold cell, not as a silent leftover. */
void Y(check_deconv_1d)(void)
{
  const int m = 6, w = Y(get_window_id)();
  Y(deconv_ensure_registered)
  ();

  /* (a) odd N = 15, n = 32: Nneg = 7, Npos = 8. */
  {
    INT N = 15, n = 32, ks;
    INT spikes[3] = {0, 7, 14}; /* slots */
    INT freqs[3] = {-7, 0, 7}; /* k(ks) = ks - 7 */
    INT poss[3] = {25, 0, 7}; /* pos(ks) */
    int i;
    C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *g = (C *)Y(malloc)((size_t)n * sizeof(C));
    for (i = 0; i < 3; i++) {
      R sc = deconv_scale(w, 1, &N, &n, m, &freqs[i]);
      problem *p;
      plan *pln;
      for (ks = 0; ks < N; ks++)
        f_hat[ks] = K(0.0);
      f_hat[spikes[i]] = K(1.0);
      for (ks = 0; ks < n; ks++)
        g[ks] = sc; /* pre-dirty at the peak magnitude */
      p = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, 1, f_hat, g);
      pln = deconv_awake_plan(p);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
      pln->adt->apply(pln, p);
      CU_ASSERT(CABS(g[poss[i]] - sc) < deconv_general_tol() * sc);
      CU_ASSERT_EQUAL(count_above(g, n, deconv_general_tol() * sc), (INT)1);
      Y(plan_destroy)
      (pln);
      Y(problem_destroy)
      (p);
    }
    /* adjoint: g[25] = 1 gathers into slot 0 only. */
    {
      INT k = -7;
      R sc = deconv_scale(w, 1, &N, &n, m, &k);
      problem *p;
      plan *pln;
      for (ks = 0; ks < N; ks++)
        f_hat[ks] = K(0.0);
      for (ks = 0; ks < n; ks++)
        g[ks] = K(0.0);
      g[25] = K(1.0);
      p = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, -1, f_hat, g);
      pln = deconv_awake_plan(p);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
      pln->adt->apply_adjoint(pln, p);
      CU_ASSERT(CABS(f_hat[0] - sc) < deconv_general_tol() * sc);
      CU_ASSERT_EQUAL(count_above(f_hat, N, deconv_general_tol() * sc), (INT)1);
      Y(plan_destroy)
      (pln);
      Y(problem_destroy)
      (p);
    }
    Y(free)
    (f_hat);
    Y(free)
    (g);
  }

  /* (b) even N = 16 type-II, n = 32: Nneg = 7, Npos = 9. */
  {
    INT N = 16, n = 32, ks;
    int v = NFFT_NDFT_TYPE_II;
    INT spikes[2] = {0, 15};
    INT freqs[2] = {-7, 8};
    INT poss[2] = {25, 8};
    int i;
    C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *g = (C *)Y(malloc)((size_t)n * sizeof(C));
    for (i = 0; i < 2; i++) {
      R sc = deconv_scale(w, 1, &N, &n, m, &freqs[i]);
      problem *p;
      plan *pln;
      for (ks = 0; ks < N; ks++)
        f_hat[ks] = K(0.0);
      f_hat[spikes[i]] = K(1.0);
      for (ks = 0; ks < n; ks++)
        g[ks] = sc;
      p = Y(mkproblem_deconv)(1, &N, &v, &n, m, w, 1, f_hat, g);
      pln = deconv_awake_plan(p);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
      pln->adt->apply(pln, p);
      CU_ASSERT(CABS(g[poss[i]] - sc) < deconv_general_tol() * sc);
      CU_ASSERT_EQUAL(count_above(g, n, deconv_general_tol() * sc), (INT)1);
      Y(plan_destroy)
      (pln);
      Y(problem_destroy)
      (p);
    }
    Y(free)
    (f_hat);
    Y(free)
    (g);
  }

  /* (c) n < N is declined, not planned: the zero-pad length would wrap. */
  {
    INT N = 32, n = 16;
    C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *g = (C *)Y(malloc)((size_t)n * sizeof(C));
    problem *p = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, 1, f_hat, g);
    CU_ASSERT_PTR_NULL(Y(planner_mkplan)(Y(the_planner)(), p));
    Y(problem_destroy)
    (p);
    Y(free)
    (f_hat);
    Y(free)
    (g);
  }

  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}

/* Rank-2 DECONV with a type-II axis and an odd axis. Slots are checked as
 * single spikes, so a one-cell index slip shows up as a miss plus a stray. */
void Y(check_deconv_2d)(void) {
  const int m = 6, w = Y(get_window_id)();
  INT N[2] = {16, 15}, n[2] = {32, 32};
  int variant[2] = {NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I};
  /* axis 0: Nneg = 7, Npos = 9. axis 1: Nneg = 7, Npos = 8. */
  INT Ntot = N[0] * N[1], ntot = n[0] * n[1], j;
  C *f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *g = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  INT slots[3][2] = {{0, 0}, {7, 7}, {15, 14}};
  INT freqs[3][2] = {{-7, -7}, {0, 0}, {8, 7}};
  INT cells[3][2] = {{25, 25}, {0, 0}, {8, 7}};
  int i;

  Y(deconv_ensure_registered)
  ();

  for (i = 0; i < 3; i++) {
    R sc = deconv_scale(w, 2, N, n, m, freqs[i]);
    INT fs = slots[i][0] * N[1] + slots[i][1];
    INT gs = cells[i][0] * n[1] + cells[i][1];
    problem *p;
    plan *pln;
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    f_hat[fs] = K(1.0);
    for (j = 0; j < ntot; j++)
      g[j] = sc;
    p = Y(mkproblem_deconv)(2, N, variant, n, m, w, 1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply(pln, p);
    CU_ASSERT(CABS(g[gs] - sc) < deconv_general_tol() * sc);
    CU_ASSERT_EQUAL(count_above(g, ntot, deconv_general_tol() * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  /* adjoint: the corner cell gathers into slot (0,0) only. */
  {
    R sc = deconv_scale(w, 2, N, n, m, freqs[0]);
    problem *p;
    plan *pln;
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    for (j = 0; j < ntot; j++)
      g[j] = K(0.0);
    g[cells[0][0] * n[1] + cells[0][1]] = K(1.0);
    p = Y(mkproblem_deconv)(2, N, variant, n, m, w, -1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply_adjoint(pln, p);
    CU_ASSERT(CABS(f_hat[0] - sc) < deconv_general_tol() * sc);
    CU_ASSERT_EQUAL(count_above(f_hat, Ntot, deconv_general_tol() * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  Y(free)
  (f_hat);
  Y(free)
  (g);
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}

/* Rank-3 DECONV with two type-II axes and one odd axis. */
void Y(check_deconv_3d)(void) {
  const int m = 6, w = Y(get_window_id)();
  INT N[3] = {16, 15, 10}, n[3] = {32, 32, 20};
  int variant[3] = {NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I, NFFT_NDFT_TYPE_II};
  /* axis 0: Nneg 7 / Npos 9. axis 1: Nneg 7 / Npos 8. axis 2: Nneg 4 / Npos 6. */
  INT Ntot = N[0] * N[1] * N[2], ntot = n[0] * n[1] * n[2], j;
  C *f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *g = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  INT slots[3][3] = {{0, 0, 0}, {7, 7, 4}, {15, 14, 9}};
  INT freqs[3][3] = {{-7, -7, -4}, {0, 0, 0}, {8, 7, 5}};
  INT cells[3][3] = {{25, 25, 16}, {0, 0, 0}, {8, 7, 5}};
  int i;

  Y(deconv_ensure_registered)
  ();

  for (i = 0; i < 3; i++) {
    R sc = deconv_scale(w, 3, N, n, m, freqs[i]);
    INT fs = (slots[i][0] * N[1] + slots[i][1]) * N[2] + slots[i][2];
    INT gs = (cells[i][0] * n[1] + cells[i][1]) * n[2] + cells[i][2];
    problem *p;
    plan *pln;
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    f_hat[fs] = K(1.0);
    for (j = 0; j < ntot; j++)
      g[j] = sc;
    p = Y(mkproblem_deconv)(3, N, variant, n, m, w, 1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply(pln, p);
    CU_ASSERT(CABS(g[gs] - sc) < deconv_general_tol() * sc);
    CU_ASSERT_EQUAL(count_above(g, ntot, deconv_general_tol() * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  /* adjoint: the corner cell gathers into slot (0,0,0) only. */
  {
    R sc = deconv_scale(w, 3, N, n, m, freqs[0]);
    INT gs = (cells[0][0] * n[1] + cells[0][1]) * n[2] + cells[0][2];
    problem *p;
    plan *pln;
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    for (j = 0; j < ntot; j++)
      g[j] = K(0.0);
    g[gs] = K(1.0);
    p = Y(mkproblem_deconv)(3, N, variant, n, m, w, -1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply_adjoint(pln, p);
    CU_ASSERT(CABS(f_hat[0] - sc) < deconv_general_tol() * sc);
    CU_ASSERT_EQUAL(count_above(f_hat, Ntot, deconv_general_tol() * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  Y(free)
  (f_hat);
  Y(free)
  (g);
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}

/* Rank-4 DECONV with mixed type-II and odd axes, driving the carry odometer. */
void Y(check_deconv_nd)(void) {
  const int m = 3, w = Y(get_window_id)();
  INT N[4] = {8, 7, 6, 5}, n[4] = {16, 16, 16, 16};
  int variant[4] = {NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I, NFFT_NDFT_TYPE_II,
                    NFFT_NDFT_TYPE_I};
  /* Nneg / Npos per axis: 3/5, 3/4, 2/4, 2/3. */
  INT Ntot = N[0] * N[1] * N[2] * N[3];
  INT ntot = n[0] * n[1] * n[2] * n[3], j;
  C *f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *g = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  INT slots[4][4] = {{0, 0, 0, 0}, {3, 3, 2, 2}, {7, 6, 5, 4}, {0, 3, 0, 4}};
  INT freqs[4][4] = {
      {-3, -3, -2, -2}, {0, 0, 0, 0}, {4, 3, 3, 2}, {-3, 0, -2, 2}};
  INT cells[4][4] = {
      {13, 13, 14, 14}, {0, 0, 0, 0}, {4, 3, 3, 2}, {13, 0, 14, 2}};
  int i, t;

  Y(deconv_ensure_registered)
  ();

  for (i = 0; i < 4; i++) {
    R sc = deconv_scale(w, 4, N, n, m, freqs[i]);
    INT fs = 0, gs = 0;
    problem *p;
    plan *pln;
    for (t = 0; t < 4; t++) {
      fs = fs * N[t] + slots[i][t];
      gs = gs * n[t] + cells[i][t];
    }
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    f_hat[fs] = K(1.0);
    for (j = 0; j < ntot; j++)
      g[j] = sc;
    p = Y(mkproblem_deconv)(4, N, variant, n, m, w, 1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply(pln, p);
    CU_ASSERT(CABS(g[gs] - sc) < deconv_general_tol() * sc);
    CU_ASSERT_EQUAL(count_above(g, ntot, deconv_general_tol() * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  /* adjoint: the corner cell gathers into slot (0,0,0,0) only. */
  {
    R sc = deconv_scale(w, 4, N, n, m, freqs[0]);
    INT gs = 0;
    problem *p;
    plan *pln;
    for (t = 0; t < 4; t++)
      gs = gs * n[t] + cells[0][t];
    for (j = 0; j < Ntot; j++)
      f_hat[j] = K(0.0);
    for (j = 0; j < ntot; j++)
      g[j] = K(0.0);
    g[gs] = K(1.0);
    p = Y(mkproblem_deconv)(4, N, variant, n, m, w, -1, f_hat, g);
    pln = deconv_awake_plan(p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    pln->adt->apply_adjoint(pln, p);
    CU_ASSERT(CABS(f_hat[0] - sc) < deconv_general_tol() * sc);
    CU_ASSERT_EQUAL(count_above(f_hat, Ntot, deconv_general_tol() * sc), (INT)1);
    Y(plan_destroy)
    (pln);
    Y(problem_destroy)
    (p);
  }

  Y(free)
  (f_hat);
  Y(free)
  (g);
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}
