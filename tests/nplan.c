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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "nfft3.h" /* legacy API (the direct oracle) + Y(malloc)/Y(free) */
#include "infft.h"
#include "iplanner.h"
#include "plan_ng_test.h" /* test-only plan_ng_test_awake_state hook */
#include "nplan.h"

/* The test-only raceable permuting native solver in tests/nplan_perm.c. */
extern int Y(nfft_perm_break_restore);
plan *Y(nfft_perm_test_mkplan)(const problem *p);
void Y(nfft_solver_perm_test_register)(planner *pl);
extern INT Y(nfft_slow_test_applies);
void Y(nfft_solver_slow_test_register)(planner *pl);

void Y(check_nplan_problem)(void) {
  planner *pl = Y(planner_create)();
  INT N1 = 16, n1 = 32;
  problem *p1, *p2;
  md5sig s1, s2;
  char buf[128];
  printer *pr;

  p1 = Y(mkproblem_nfft)(1, &N1, 0, &n1, (INT)1000, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, 0, 0, 0, 0);
  CU_ASSERT_EQUAL(p1->adt->kind, NFFT_PROBLEM_NFFT);

  pr = Y(printer_create_str)(buf);
  pr->print(pr, "%P", p1);
  Y(printer_destroy)
  (pr);
  CU_ASSERT_STRING_EQUAL(buf, "(nfft sign=1 m=6 M=1000 (tensor 1 (16 1 32 1)) variant=0)");

  /* M-bucketing: 1000 and 1023 share a bucket; 1024 does not */
  Y(problem_md5)
  (pl, p1, s1);
  p2 = Y(mkproblem_nfft)(1, &N1, 0, &n1, (INT)1023, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, 0, 0, 0, 0);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2] && s1[3] == s2[3]);
  Y(problem_destroy)
  (p2);
  p2 = Y(mkproblem_nfft)(1, &N1, 0, &n1, (INT)1024, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, 0, 0, 0, 0);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);
  Y(problem_destroy)
  (p2);

  /* the adjoint problem carries the adjoint tensor and its own key */
  p2 = Y(mkproblem_nfft)(1, &N1, 0, &n1, (INT)1000, 6, NFFT_WINDOW_KAISER_BESSEL, -1, 0u, 0, 0, 0, 0);
  pr = Y(printer_create_str)(buf);
  pr->print(pr, "%P", p2);
  Y(printer_destroy)
  (pr);
  CU_ASSERT_STRING_EQUAL(buf, "(nfft sign=-1 m=6 M=1000 (tensor 1 (32 1 16 1)) variant=0)");
  CU_ASSERT_EQUAL(Y(problem_nfft_N)(p2, 0), (INT)16); /* direction-aware */
  CU_ASSERT_EQUAL(Y(problem_nfft_n)(p2, 0), (INT)32);
  CU_ASSERT_EQUAL(Y(problem_nfft_Ntot)(p2), (INT)16);
  CU_ASSERT_EQUAL(Y(problem_nfft_ntot)(p2), (INT)32);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);
  Y(problem_destroy)
  (p2);

  /* forward accessors read the forward orientation */
  CU_ASSERT_EQUAL(Y(problem_nfft_N)(p1, 0), (INT)16);
  CU_ASSERT_EQUAL(Y(problem_nfft_n)(p1, 0), (INT)32);

  /* geometry, m, and fftw_flags all shift the key */
  {
    INT N2 = 32;
    p2 = Y(mkproblem_nfft)(1, &N2, 0, &n1, (INT)1000, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, 0, 0, 0, 0);
    Y(problem_md5)
    (pl, p2, s2);
    CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);
    Y(problem_destroy)
    (p2);
  }
  p2 = Y(mkproblem_nfft)(1, &N1, 0, &n1, (INT)1000, 8, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, 0, 0, 0, 0);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);
  Y(problem_destroy)
  (p2);
  p2 = Y(mkproblem_nfft)(1, &N1, 0, &n1, (INT)1000, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 64u, 0, 0, 0, 0);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);
  Y(problem_destroy)
  (p2);

  /* symmetric geometry (N == n): the tensors coincide under adjoint and
   * only sign separates the directions' keys */
  {
    INT Ns = 16;
    problem *pf = Y(mkproblem_nfft)(1, &Ns, 0, &Ns, (INT)1000, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, 0, 0, 0, 0);
    problem *pa = Y(mkproblem_nfft)(1, &Ns, 0, &Ns, (INT)1000, 6, NFFT_WINDOW_KAISER_BESSEL, -1, 0u, 0, 0, 0, 0);
    md5sig sf, sa;
    Y(problem_md5)
    (pl, pf, sf);
    Y(problem_md5)
    (pl, pa, sa);
    CU_ASSERT(sf[0] != sa[0] || sf[1] != sa[1] || sf[2] != sa[2] || sf[3] != sa[3]);
    Y(problem_destroy)
    (pf);
    Y(problem_destroy)
    (pa);
  }

  Y(problem_destroy)
  (p1);
  Y(planner_destroy)
  (pl);
}

/* The compressed rank and surviving-axis order are pinned numerically by
 * check_nplan_unit_axis_correct, not here. */
void Y(check_nplan_elides_unit_axes_geometry)(void) {
  /* (a) A borrowed x cannot be gathered into a fresh buffer, so the copy_x=0
   * path keeps full rank and a unit-padded problem hashes distinctly from the
   * lower-rank one. On the copy_x=1 path the same geometry compresses to
   * exactly the lower-rank problem, so the two hashes coincide. */
  {
    INT N3[3] = {64, 64, 1}, n3[3] = {128, 128, 1};
    INT N2[2] = {64, 64}, n2[2] = {128, 128};
    R x3[3 * 100] = {K(0.0)};
    planner *pl = Y(planner_create)();
    problem *p, *p2, *p3;
    md5sig s1, s2, s3;

    p = Y(mkproblem_nfft)(3, N3, 0, n3, (INT)100, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, 0, 0, 0, 0);
    CU_ASSERT_EQUAL(((const problem_nfft *)p)->sz->rnk, 3); /* borrowed: full rank kept */
    CU_ASSERT_EQUAL(Y(problem_nfft_N)(p, 2), (INT)1);       /* unit axis visible */
    Y(problem_md5)
    (pl, p, s1);
    p2 = Y(mkproblem_nfft)(2, N2, 0, n2, (INT)100, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, 0, 0, 0, 0);
    Y(problem_md5)
    (pl, p2, s2);
    CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);

    p3 = Y(mkproblem_nfft)(3, N3, 0, n3, (INT)100, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, x3, 1, 0, 0);
    CU_ASSERT_EQUAL(((const problem_nfft *)p3)->sz->rnk, 2); /* top-level: compressed */
    Y(problem_md5)
    (pl, p3, s3);
    CU_ASSERT(s3[0] == s2[0] && s3[1] == s2[1] && s3[2] == s2[2] && s3[3] == s2[3]);

    Y(problem_destroy)
    (p);
    Y(problem_destroy)
    (p2);
    Y(problem_destroy)
    (p3);
    Y(planner_destroy)
    (pl);
  }

  /* (b) Elision is drop-only, with no canonical sort, so axis reordering
   * yields distinct problems and hashes: the x-column pairing differs. */
  {
    INT Na[2] = {32, 64}, na[2] = {64, 128};
    INT Nb[2] = {64, 32}, nb[2] = {128, 64};
    planner *pl = Y(planner_create)();
    problem *p, *p2;
    md5sig s1, s2;
    p = Y(mkproblem_nfft)(2, Na, 0, na, (INT)100, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, 0, 0, 0, 0);
    p2 = Y(mkproblem_nfft)(2, Nb, 0, nb, (INT)100, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, 0, 0, 0, 0);
    Y(problem_md5)
    (pl, p, s1);
    Y(problem_md5)
    (pl, p2, s2);
    CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);
    Y(problem_destroy)
    (p);
    Y(problem_destroy)
    (p2);
    Y(planner_destroy)
    (pl);
  }
}

void Y(check_nplan_variant_key)(void) {
  INT N1 = 16, n1 = 32, No = 15, no = 32;
  int v_i = NFFT_NDFT_TYPE_I, v_ii = NFFT_NDFT_TYPE_II;
  R x = K(0.1);
  C fh = K(1.0), f = K(0.0);
  md5sig k_i, k_ii, ko_i, ko_ii;
  /* x is a single scalar, not a d*M array: copy_x=0 so mkproblem_nfft never
   * reads past it. These calls probe the key, not x's contents. */
  problem *pi = Y(mkproblem_nfft)(1, &N1, &v_i, &n1, 1000, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, &x, 0, &fh, &f);
  problem *pii = Y(mkproblem_nfft)(1, &N1, &v_ii, &n1, 1000, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, &x, 0, &fh, &f);
  Y(problem_md5)
  (Y(the_planner)(), pi, k_i);
  Y(problem_md5)
  (Y(the_planner)(), pii, k_ii);
  CU_ASSERT(k_i[0] != k_ii[0] || k_i[1] != k_ii[1] ||
            k_i[2] != k_ii[2] || k_i[3] != k_ii[3]);
  {
    problem *pn = Y(mkproblem_nfft)(1, &N1, 0, &n1, 1000, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, &x, 0, &fh, &f);
    md5sig kn;
    Y(problem_md5)
    (Y(the_planner)(), pn, kn);
    CU_ASSERT(kn[0] == k_i[0] && kn[1] == k_i[1] &&
              kn[2] == k_i[2] && kn[3] == k_i[3]); /* NULL == explicit type-I */
    Y(problem_destroy)
    (pn);
  }
  CU_ASSERT_EQUAL(Y(problem_nfft_variant)(pii, 0), NFFT_NDFT_TYPE_II);
  problem *poi = Y(mkproblem_nfft)(1, &No, &v_i, &no, 1000, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, &x, 0, &fh, &f);
  problem *poii = Y(mkproblem_nfft)(1, &No, &v_ii, &no, 1000, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, &x, 0, &fh, &f);
  Y(problem_md5)
  (Y(the_planner)(), poi, ko_i);
  Y(problem_md5)
  (Y(the_planner)(), poii, ko_ii);
  CU_ASSERT(ko_i[0] == ko_ii[0] && ko_i[1] == ko_ii[1] &&
            ko_i[2] == ko_ii[2] && ko_i[3] == ko_ii[3]); /* odd collapses */
  CU_ASSERT_EQUAL(Y(problem_nfft_variant)(poii, 0), NFFT_NDFT_TYPE_I);
  Y(problem_destroy)
  (pi);
  Y(problem_destroy)
  (pii);
  Y(problem_destroy)
  (poi);
  Y(problem_destroy)
  (poii);
}

/* mkproblem_nfft borrows the caller's f_hat/f and never frees them. x on a
 * copy_x=1 problem is instead a private copy, equal in content and freed with
 * the problem. */
void Y(check_nplan_problem_owns_data)(void) {
  INT N = 16, n = 32, M = 100;
  static R x[100];
  static C f_hat[16], f[100];
  problem *p;
  const problem_nfft *pn;
  INT j;

  p = Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, x, /*copy_x=*/1, f_hat, f);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  pn = (const problem_nfft *)p;
  CU_ASSERT_PTR_NOT_EQUAL(pn->x, x); /* x is copied, not aliased */
  CU_ASSERT_TRUE(pn->x_owned);
  for (j = 0; j < M; j++)
    CU_ASSERT_EQUAL(pn->x[j], x[j]); /* same content */
  CU_ASSERT_PTR_EQUAL(pn->f_hat, f_hat);
  CU_ASSERT_PTR_EQUAL(pn->f, f);
  Y(problem_destroy)
  (p);
}

/* The wisdom key is data-blind: two problems of identical shape hash the same
 * whatever their data pointers, so the cache is keyed on geometry. */
void Y(check_nplan_key_is_data_blind)(void) {
  planner *pl = Y(planner_create)();
  INT N = 16, n = 32, M = 100;
  static R x[100];
  static C f_hat[16], f[100];
  problem *p0, *p1;
  md5sig s0, s1;

  p0 = Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, 0, 0, 0, 0);
  p1 = Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, 1, 0u, x, /*copy_x=*/1, f_hat, f);
  Y(problem_md5)
  (pl, p0, s0);
  Y(problem_md5)
  (pl, p1, s1);
  CU_ASSERT(s0[0] == s1[0] && s0[1] == s1[1] && s0[2] == s1[2] && s0[3] == s1[3]);

  Y(problem_destroy)
  (p0);
  Y(problem_destroy)
  (p1);
  Y(planner_destroy)
  (pl);
}

/* Both helpers request the compile-time window: the native fast solver is
 * window-gated (kernel/nfft/nfft-nd.c), so a hardcoded KB problem would make it
 * decline under a non-KB build and the expected winner would change. */
static problem *mk_problem(int d, INT Nv, INT nv, INT M, int m, int sign) {
  INT N[4], n[4];
  int t;
  for (t = 0; t < d; t++) {
    N[t] = Nv;
    n[t] = nv;
  }
  return Y(mkproblem_nfft)(d, N, 0, n, M, m, Y(get_window_id)(), sign, 0u, 0, 0, 0, 0);
}

static problem *mk_problem_arr(int d, const INT *N, const INT *n, INT M, int m,
                          int sign) {
  return Y(mkproblem_nfft)(d, N, 0, n, M, m, Y(get_window_id)(), sign, 0u, 0, 0, 0, 0);
}

/* Assert the winner is exactly solver `nam`. The print format is
 * "(<name> pcost=<int>)", so bracketing the name with '(' and ' ' pins the
 * token. */
static void expect_winner(planner *pl, problem *p, const char *nam) {
  plan *pln = Y(planner_mkplan)(pl, p);
  /* 512: printer_create_str is unbounded, and the bundle print nests both
   * children plus the FFTW plan description, so a smaller buffer overflows
   * the stack instead of truncating. */
  char buf[512];
  char tok[64];
  printer *pr;
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  pr = Y(printer_create_str)(buf);
  pr->print(pr, "%p", pln);
  Y(printer_destroy)
  (pr);
  tok[0] = '(';
  strncpy(tok + 1, nam, sizeof(tok) - 3);
  tok[sizeof(tok) - 3] = '\0';
  strcat(tok, " ");
  CU_ASSERT_PTR_NOT_NULL(strstr(buf, tok));
  Y(plan_destroy)
  (pln);
}

/* The fast solver declines sigma = n/N <= 1 (guards_ok, kernel/nfft/nfft-nd.c)
 * and a direct native takes those. The guru rejects sigma <= 1 up front, so the
 * decline is reachable only through mkproblem_nfft, as below. */
void Y(check_nplan_guard_declines)(void) {
  planner *pl;
  problem *p;
  Y(nfft_ensure_registered)
  ();
  pl = Y(the_planner)();
  Y(planner_forget)
  (pl, PLNR_FORGET_ALL);
  {
    R x = K(0.1);
    C fh = K(0.0), f = K(0.0);
    INT N = 63, n = 128; /* odd, big M */
    /* x is a scalar, not an M-element array: copy_x=0 so mkproblem_nfft never
     * reads past it. */
    p = Y(mkproblem_nfft)(1, &N, 0, &n, 200000, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, &x, 0, &fh, &f);
    expect_winner(pl, p, "nfft_solver_fast_native");
    Y(problem_destroy)
    (p);
  }
  {
    R x = K(0.1);
    C fh = K(0.0), f = K(0.0);
    INT N = 64, n = 128;
    int v = NFFT_NDFT_TYPE_II; /* type-II, big M */
    p = Y(mkproblem_nfft)(1, &N, &v, &n, 200000, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, &x, 0, &fh, &f);
    expect_winner(pl, p, "nfft_solver_fast_native");
    Y(problem_destroy)
    (p);
  }
  {
    R x = K(0.1);
    C fh = K(0.0), f = K(0.0);
    INT N = 64, n = 64; /* sigma == 1: declined even at big M */
    p = Y(mkproblem_nfft)(1, &N, 0, &n, 200000, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, &x, 0, &fh, &f);
    expect_winner(pl, p, "nfft_solver_ndft_1d_blocked");
    Y(problem_destroy)
    (p);
  }
  /* NO_DIRECT with no applicable fast plan -> NULL. N=4 fails the fast guard's
   * N > m at m=6, and NO_DIRECT removes the NDFT fallback. sigma <= 1 can no
   * longer reach the search; the guru rejects it first. */
  {
    INT N = 4, n = 16, M = 1000;
    R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
    C *fh = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
    INT j;
    for (j = 0; j < M; j++)
      x[j] = K(0.0);
    CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, fh, f, 0u,
                                       NFFT_ESTIMATE | NFFT_NO_DIRECT));
    Y(free)
    (x);
    Y(free)
    (fh);
    Y(free)
    (f);
  }
  /* NO_DIRECT + odd: the fast serves it. */
  {
    INT N = 63, n = 128, M = 1000;
    R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
    C *fh = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
    Y(plan_ng) * pg;
    INT j;
    for (j = 0; j < M; j++)
      x[j] = K(0.0);
    pg = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, fh, f, 0u,
                         NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL(pg);
    Y(plan_ng_destroy)
    (pg);
    Y(free)
    (x);
    Y(free)
    (fh);
    Y(free)
    (f);
  }
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}

void Y(check_nplan_solvers)(void) {
  planner *pl = Y(planner_create)();
  problem *p;

  /* mkplan_native_fast builds its DECONV/CONV children through the
   * process-global planner, not the local pl, so that instance must have them
   * registered or the fast solver declines every problem below. */
  Y(nfft_ensure_registered)
  ();

  Y(nfft_solvers_register)
  (pl);

  /* (a) big M, 1d: fast wins the estimate */
  p = mk_problem(1, 256, 512, 4096, 6, 1);
  expect_winner(pl, p, "nfft_solver_fast_native");
  Y(problem_destroy)
  (p);

  /* (b) degenerate size: guards decline fast -> direct-class (blocked NDFT) */
  p = mk_problem(1, 4, 16, 8, 6, 1);
  expect_winner(pl, p, "nfft_solver_ndft_1d_blocked");
  Y(problem_destroy)
  (p);

  /* (c) tiny M: estimate picks direct-class (blocked NDFT; 4x oversampling
   * keeps the margin ~4x) */
  p = mk_problem(1, 256, 1024, 1, 6, 1);
  expect_winner(pl, p, "nfft_solver_ndft_1d_blocked");
  Y(problem_destroy)
  (p);

  /* (d) NO_DIRECT: fast where legal; NULL where not */
  pl->flags.l = PLNR_NO_DIRECT;
  pl->flags.u = PLNR_NO_DIRECT;
  p = mk_problem(1, 256, 1024, 1, 6, 1);
  expect_winner(pl, p, "nfft_solver_fast_native");
  Y(problem_destroy)
  (p);
  p = mk_problem(1, 4, 16, 8, 6, 1);
  CU_ASSERT_PTR_NULL(Y(planner_mkplan)(pl, p));
  Y(problem_destroy)
  (p);
  pl->flags.l = 0;
  pl->flags.u = 0;

  /* (e) 2d analog of (a) */
  p = mk_problem(2, 64, 128, 4096, 6, 1);
  expect_winner(pl, p, "nfft_solver_fast_native");
  Y(problem_destroy)
  (p);

  /* adjoint direction selects independently (its own wisdom entry) */
  p = mk_problem(1, 256, 512, 4096, 6, -1);
  expect_winner(pl, p, "nfft_solver_fast_native");
  Y(problem_destroy)
  (p);

  /* d=3/d=4: the planner-native fast wins both */
  p = mk_problem(3, 16, 32, 4096, 3, 1);
  expect_winner(pl, p, "nfft_solver_fast_native");
  Y(problem_destroy)
  (p);
  p = mk_problem(4, 8, 16, 4096, 3, 1);
  expect_winner(pl, p, "nfft_solver_fast_native");
  Y(problem_destroy)
  (p);

  /* (f) unit axis among non-unit axes: every rank >= 1 solver declines, so no
   * plan exists. mk_problem_arr borrows x, so the problem keeps full rank --
   * the only way to build such a problem. Leading, trailing and adjoint all
   * decline. */
  {
    INT Nu3[3] = {64, 64, 1}, nu3[3] = {128, 128, 1};
    INT Nu2[2] = {1, 64}, nu2[2] = {1, 128};
    p = mk_problem_arr(3, Nu3, nu3, 4096, 6, 1);
    CU_ASSERT_EQUAL(Y(problem_nfft_has_unit_axis)(p), 1);
    CU_ASSERT_PTR_NULL(Y(planner_mkplan)(pl, p));
    Y(problem_destroy)
    (p);
    p = mk_problem_arr(2, Nu2, nu2, 4096, 6, 1);
    CU_ASSERT_PTR_NULL(Y(planner_mkplan)(pl, p));
    Y(problem_destroy)
    (p);
    p = mk_problem_arr(3, Nu3, nu3, 4096, 6, -1);
    CU_ASSERT_PTR_NULL(Y(planner_mkplan)(pl, p));
    Y(problem_destroy)
    (p);
  }

  Y(planner_destroy)
  (pl);

  /* Registration re-arms across global-planner generations; pointer identity
   * is not a valid generation check, since a recreated planner can reuse the
   * freed address. nslvdesc counts descriptors of every kind: the NFFT roster
   * in kernel/nfft/conf.c is 5 (composed native fast, two native 1D direct,
   * generic nD NDFT, rank-0 base case), plus the DECONV and CONV rosters of 4
   * each, which ensure_registered installs first because the composed fast
   * recurses into them. */
  Y(nfft_ensure_registered)
  ();
  CU_ASSERT_EQUAL(Y(the_planner)()->nslvdesc, 13u);
  Y(the_planner_destroy)
  ();
  Y(nfft_ensure_registered)
  ();
  CU_ASSERT_EQUAL(Y(the_planner)()->nslvdesc, 13u); /* fresh generation */
  Y(the_planner_destroy)
  ();
}

/* dispatch: for d == 1 the two planner-native NDFT solvers compete and each
 * per-variant flag prunes one. */
void Y(check_nplan_ndft_dispatch)(void) {
  planner *pl = Y(planner_create)();
  problem *p;
  /* See check_nplan_solvers: the fast solver's DECONV/CONV children recurse
   * via the process-global planner, so it must be registered too. */
  Y(nfft_ensure_registered)
  ();
  Y(nfft_solvers_register)
  (pl);

  /* Tiny M keeps a direct-class solver ahead of fast; blocked beats plain. */
  p = mk_problem(1, 256, 1024, 1, 6, 1);
  expect_winner(pl, p, "nfft_solver_ndft_1d_blocked");
  Y(problem_destroy)
  (p);

  /* NO_NDFT_BLOCKED leaves the plain solver. */
  pl->flags.l = pl->flags.u = PLNR_NO_NDFT_BLOCKED;
  p = mk_problem(1, 256, 1024, 1, 6, 1);
  expect_winner(pl, p, "nfft_solver_ndft_1d");
  Y(problem_destroy)
  (p);

  /* NO_NDFT_PLAIN leaves the blocked solver. */
  pl->flags.l = pl->flags.u = PLNR_NO_NDFT_PLAIN;
  p = mk_problem(1, 256, 1024, 1, 6, 1);
  expect_winner(pl, p, "nfft_solver_ndft_1d_blocked");
  Y(problem_destroy)
  (p);

  /* Both disabled + fast inapplicable (degenerate size) -> no plan at all. */
  pl->flags.l = pl->flags.u = PLNR_NO_NDFT_PLAIN | PLNR_NO_NDFT_BLOCKED;
  p = mk_problem(1, 4, 16, 8, 6, 1);
  CU_ASSERT_PTR_NULL(Y(planner_mkplan)(pl, p));
  Y(problem_destroy)
  (p);

  /* NO_DIRECT subsumes both native solvers; degenerate 1D -> no plan. */
  pl->flags.l = pl->flags.u = PLNR_NO_DIRECT;
  p = mk_problem(1, 4, 16, 8, 6, 1);
  CU_ASSERT_PTR_NULL(Y(planner_mkplan)(pl, p));
  Y(problem_destroy)
  (p);
  /* ...but where fast is legal it still wins under NO_DIRECT. */
  p = mk_problem(1, 256, 512, 4096, 6, 1);
  expect_winner(pl, p, "nfft_solver_fast_native");
  Y(problem_destroy)
  (p);

  /* Adjoint direction dispatches identically (its own problem/sign). */
  pl->flags.l = pl->flags.u = 0;
  p = mk_problem(1, 256, 1024, 1, 6, -1);
  expect_winner(pl, p, "nfft_solver_ndft_1d_blocked");
  Y(problem_destroy)
  (p);

  Y(planner_destroy)
  (pl);
}

/* dispatch: for d >= 2 the single generic native multivariate NDFT solver
 * wins whenever a direct plan is legal. NO_FAST_NATIVE must be set explicitly
 * to steer away from the fast solver. */
void Y(check_nplan_ndft_multivariate_dispatch)(void) {
  planner *pl = Y(planner_create)();
  problem *p;
  /* See check_nplan_solvers: the fast solver's DECONV/CONV children recurse
   * via the process-global planner, so it must be registered too. */
  Y(nfft_ensure_registered)
  ();
  Y(nfft_solvers_register)
  (pl);

  /* d=2/d=3: NO_FAST_NATIVE declines the native fast; the direct native NDFT
   * remains and wins. */
  pl->flags.l = pl->flags.u = PLNR_NO_FAST_NATIVE;
  p = mk_problem(2, 32, 64, 4096, 6, 1);
  expect_winner(pl, p, "nfft_solver_ndft_nd");
  Y(problem_destroy)
  (p);
  p = mk_problem(2, 32, 64, 4096, 6, -1); /* adjoint dispatches on its own sign */
  expect_winner(pl, p, "nfft_solver_ndft_nd");
  Y(problem_destroy)
  (p);
  p = mk_problem(3, 16, 32, 4096, 3, 1);
  expect_winner(pl, p, "nfft_solver_ndft_nd");
  Y(problem_destroy)
  (p);
  p = mk_problem(3, 16, 32, 4096, 3, -1);
  expect_winner(pl, p, "nfft_solver_ndft_nd");
  Y(problem_destroy)
  (p);
  pl->flags.l = pl->flags.u = 0;

  /* d=4: the native fast is not disabled here, so force the direct native by
   * cost -- tiny M makes the O(Ntot*M) direct estimate beat the fast estimate. */
  p = mk_problem(4, 8, 16, 1, 6, 1);
  expect_winner(pl, p, "nfft_solver_ndft_nd");
  Y(problem_destroy)
  (p);
  /* big-M d=4 still picks the native fast */
  p = mk_problem(4, 8, 16, 4096, 3, 1);
  expect_winner(pl, p, "nfft_solver_fast_native");
  Y(problem_destroy)
  (p);

  Y(planner_destroy)
  (pl);
}

/* deterministic LCG so references are reproducible across runs */
static unsigned lcg_state;
static R lcg_unit(void) {
  lcg_state = lcg_state * 1664525u + 1013904223u;
  return ((R)(lcg_state >> 8) / (R)16777216.0) - (R)0.5;
}

static void fill_nodes(R *x, int d, INT M, unsigned seed) {
  INT j;
  lcg_state = seed;
  for (j = 0; j < d * M; j++)
    x[j] = lcg_unit();
}

static void fill_fhat(C *f_hat, INT Ntot, unsigned seed) {
  INT j;
  lcg_state = seed;
  for (j = 0; j < Ntot; j++)
    f_hat[j] = lcg_unit() + _Complex_I * lcg_unit();
}

static R rel_max_err(const C *a, const C *b, INT len) {
  R num = (R)0, den = (R)0;
  INT j;
  for (j = 0; j < len; j++) {
    R e = CABS(a[j] - b[j]);
    if (e > num)
      num = e;
    if (CABS(b[j]) > den)
      den = CABS(b[j]);
  }
  return den > (R)0 ? num / den : num;
}

/* Per-window accuracy bound for the native fast pipeline against the exact
 * NDFT oracle. The a/b calibration follows err_trafo in tests/nfft.c, with the
 * fast pipeline's 56*eps round-off floor. A single KB-calibrated tolerance is
 * far too tight for gaussian/bspline/sinc, whose true error at m=6/sigma=2 is
 * ~1e-6..1e-7. */
static R err_bound(int window, R m, R s) {
  R eps = Y(float_property)(NFFT_EPSILON), a, b, err;
  switch (window) {
  case NFFT_WINDOW_GAUSSIAN:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    a = K(0.41);
    b = K(50.0);
#else
    a = K(0.95);
    b = K(50.0);
#endif
    err = EXP(-m * KPI * (K(1.0) - K(1.0) / (K(2.0) * K(2.0) - K(1.0))));
    break;
  case NFFT_WINDOW_B_SPLINE:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    a = K(1.0);
    b = K(2000.0);
#else
    a = K(0.3);
    b = K(50.0);
#endif
    err = K(3000.0) * K(4.0) * POW(K(1.0) / (K(2.0) * s - K(1.0)), K(2.0) * m);
    break;
  case NFFT_WINDOW_SINC_POWER:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    a = K(1.0);
    b = K(2000.0);
#else
    a = K(0.3);
    b = K(50.0);
#endif
    err = (K(1.0) / (m - K(1.0))) * ((K(2.0) / POW(s, K(2.0) * m)) + POW(s / (K(2.0) * s - K(1.0)), K(2.0) * m));
    break;
  case NFFT_WINDOW_KAISER_BESSEL:
  default:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    a = K(0.3);
    b = K(2100.0);
#else
    a = K(1.5);
    b = K(50.0);
#endif
    err = KPI * (SQRT(m) + m) * SQRT(SQRT(K(1.0) - K(1.0) / K(2.0))) * EXP(-K2PI * m * SQRT(K(1.0) - K(1.0) / K(2.0)));
    break;
  }
  return FMAX(FMAX(a * err, b * eps), K(56.0) * eps);
}

static void check_case_against_direct(int d, INT Nv, INT nv, INT M,
                                      unsigned seed, unsigned planning) {
  INT N[4], n[4], Ntot = 1;
  int t;
  INT j;
  R *xin;
  C *f_hat, *f;
  Y(plan_ng) * p;
  NFFT(plan)
  ref;
  R err;

  for (t = 0; t < d; t++) {
    N[t] = Nv;
    n[t] = nv;
    Ntot *= Nv;
  }
  xin = (R *)Y(malloc)((size_t)(d * M) * sizeof(R));
  fill_nodes(xin, d, M, seed);
  f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  f = (C *)Y(malloc)((size_t)M * sizeof(C));

  p = Y(plan_ng_guru)(d, N, 0, n, M, 6, NFFT_WINDOW_KAISER_BESSEL, xin, f_hat, f, 0u, planning);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);
  /* fill data after planning: a measured guru clobbers f_hat/f */
  fill_fhat(f_hat, Ntot, seed + 1000u);
  Y(execute)
  (p);

  /* legacy direct oracle on identical inputs */
  {
    int Ni[4], ni[4];
    for (t = 0; t < d; t++) {
      Ni[t] = (int)N[t];
      ni[t] = (int)n[t];
    }
    NFFT(init_guru)
    (&ref, d, Ni, (int)M, ni, 6, MALLOC_X | MALLOC_F_HAT | MALLOC_F, 0u);
  }
  for (j = 0; j < d * M; j++)
    ref.x[j] = xin[j];
  for (j = 0; j < Ntot; j++)
    ref.f_hat[j] = f_hat[j];
  NFFT(trafo_direct)
  (&ref);
  err = rel_max_err(f, ref.f, M);
  /* max(1e-10, 1e4*eps): window truncation at m=6 is precision-independent,
   * ~4e-12 in 1d and ~1e-11 in 2d */
  {
    R tol = (R)1.0e-10;
    if ((R)1.0e4 * EPSILON > tol)
      tol = (R)1.0e4 * EPSILON;
    CU_ASSERT(err < tol);
  }

  /* adjoint: same bundle, planner-chosen adjoint vs direct oracle; the
   * trafo result serves as the shared adjoint input (already in ref.f) */
  for (j = 0; j < M; j++)
    f[j] = ref.f[j];
  Y(execute_adjoint)
  (p);
  NFFT(adjoint_direct)
  (&ref);
  err = rel_max_err(f_hat, ref.f_hat, Ntot);
  {
    R tol = (R)1.0e-10;
    if ((R)1.0e4 * EPSILON > tol)
      tol = (R)1.0e4 * EPSILON;
    CU_ASSERT(err < tol);
  }

  NFFT(finalize)
  (&ref);
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (xin);
  Y(free)
  (f_hat);
  Y(free)
  (f);
}

/* Brute-force direct NDFT reference. k[] is decoded per k_L by modular
 * arithmetic (row-major, last axis fastest), which stays correct for unit axes
 * (k[unit] == 0). The legacy X(trafo_direct)/X(adjoint_direct) odometer has a
 * broken carry test for a non-outermost unit axis, so the unit-axis cases in
 * check_case_against_direct_arr must validate against this, not against the
 * legacy oracle. */
static void ndft_ref_trafo(int d, const INT *N, const int *variant, INT Ntot,
                           INT M, const R *x, const C *f_hat, C *f) {
  INT j, kl;
  for (j = 0; j < M; j++) {
    C v = K(0.0);
    for (kl = 0; kl < Ntot; kl++) {
      INT tmp = kl;
      R omega = K(0.0);
      int t;
      for (t = d - 1; t >= 0; t--) {
        INT slot = tmp % N[t];    /* array slot 0..N-1 on this axis */
        INT kt = slot - N[t] / 2; /* type-I; correct for odd (symmetric) */
        tmp /= N[t];
        /* type-II even axis: ascending range, uniform +1 shift (odd axes are
         * normalized to type-I upstream, so the N even guard suffices here). */
        if (variant && variant[t] == NFFT_NDFT_TYPE_II && N[t] % 2 == 0)
          kt += 1;
        omega += (R)kt * K2PI * x[j * d + t];
      }
      v += f_hat[kl] * (COS(omega) - II * SIN(omega));
    }
    f[j] = v;
  }
}

static void ndft_ref_adjoint(int d, const INT *N, const int *variant,
                             INT Ntot, INT M, const R *x, const C *f,
                             C *f_hat) {
  INT j, kl;
  for (kl = 0; kl < Ntot; kl++) {
    INT tmp = kl, kt[d];
    int t;
    C acc = K(0.0);
    for (t = d - 1; t >= 0; t--) {
      INT slot = tmp % N[t];   /* array slot 0..N-1 on this axis */
      kt[t] = slot - N[t] / 2; /* type-I; correct for odd (symmetric) */
      tmp /= N[t];
      /* type-II even axis: ascending range, uniform +1 shift (odd axes are
       * normalized to type-I upstream, so the N even guard suffices here). */
      if (variant && variant[t] == NFFT_NDFT_TYPE_II && N[t] % 2 == 0)
        kt[t] += 1;
    }
    for (j = 0; j < M; j++) {
      R omega = K(0.0);
      for (t = 0; t < d; t++)
        omega += (R)kt[t] * K2PI * x[j * d + t];
      acc += f[j] * (COS(omega) + II * SIN(omega));
    }
    f_hat[kl] = acc;
  }
}

/* Array-geometry variant of check_case_against_direct: per-axis N[]/n[], so a
 * unit axis can be embedded. Validates both directions against the
 * unit-axis-safe reference above. */
static void check_case_against_direct_arr(int d, const INT *N, const INT *n,
                                          INT M, unsigned seed,
                                          unsigned planning) {
  INT Ntot = 1;
  int t;
  INT j;
  R *xin;
  C *f_hat, *f, *ref_f, *ref_fhat;
  Y(plan_ng) * p;
  R err, tol;

  for (t = 0; t < d; t++)
    Ntot *= N[t];

  xin = (R *)Y(malloc)((size_t)(d * M) * sizeof(R));
  fill_nodes(xin, d, M, seed);
  f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  f = (C *)Y(malloc)((size_t)M * sizeof(C));
  ref_f = (C *)Y(malloc)((size_t)M * sizeof(C));
  ref_fhat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));

  p = Y(plan_ng_guru)(d, N, 0, n, M, 6, NFFT_WINDOW_KAISER_BESSEL, xin, f_hat, f, 0u, planning);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);
  fill_fhat(f_hat, Ntot, seed + 1000u); /* after planning: measured clobbers */
  Y(execute)
  (p);

  ndft_ref_trafo(d, N, 0, Ntot, M, xin, f_hat, ref_f);
  err = rel_max_err(f, ref_f, M);
  tol = (R)1.0e-10;
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;
  CU_ASSERT(err < tol);

  /* adjoint input = forward output (from the correct reference) */
  for (j = 0; j < M; j++)
    f[j] = ref_f[j];
  Y(execute_adjoint)
  (p);
  ndft_ref_adjoint(d, N, 0, Ntot, M, xin, ref_f, ref_fhat);
  err = rel_max_err(f_hat, ref_fhat, Ntot);
  CU_ASSERT(err < tol);

  Y(plan_ng_destroy)
  (p);
  Y(free)
  (xin);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(free)
  (ref_f);
  Y(free)
  (ref_fhat);
}

void Y(check_nplan_correct)(void) {
  check_case_against_direct(1, 256, 512, 4096, 42u, NFFT_ESTIMATE);
  check_case_against_direct(1, 4, 16, 8, 43u, NFFT_ESTIMATE);
  check_case_against_direct(2, 32, 64, 512, 44u, NFFT_ESTIMATE);
  /* the measured bundle stays correct, whatever won the race */
  check_case_against_direct(1, 64, 128, 2048, 45u, NFFT_MEASURE);
  /* lone-direct measured race (N<=m, so fast declines): the value-blind
   * zeroing must not read the uninitialised psi pointer of a core with no psi
   * allocated */
  check_case_against_direct(1, 4, 16, 8, 46u, NFFT_MEASURE);
  /* d=2 fast measured race: the union core carries NFFT_SORT_NODES, so the
   * value-blind race reads index_x as a gather permutation and must zero it,
   * or the SLEEPY race reads out of bounds */
  check_case_against_direct(2, 32, 64, 512, 48u, NFFT_MEASURE);

  /* NO_DIRECT + degenerate: both directions absent, so the guru returns NULL.
   * x/f_hat/f are still required arguments, though never read here. */
  {
    INT N = 4, n = 16, M = 8;
    R x[8];
    C f_hat[4], f[8];
    fill_nodes(x, 1, M, 900u);
    CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u,
                                       NFFT_ESTIMATE | NFFT_NO_DIRECT));
  }

  /* unit axis embedded among non-unit axes: the guru elides it at
   * construction and the compressed problem stays correct in both directions
   * and both modes, whichever solver wins. */
  {
    INT Nu3[3] = {64, 64, 1}, nu3[3] = {128, 128, 1};
    INT Nu2[2] = {1, 64}, nu2[2] = {1, 128};
    check_case_against_direct_arr(3, Nu3, nu3, 512, 61u, NFFT_ESTIMATE);
    check_case_against_direct_arr(3, Nu3, nu3, 512, 62u, NFFT_MEASURE);
    check_case_against_direct_arr(2, Nu2, nu2, 512, 63u, NFFT_ESTIMATE);
    check_case_against_direct_arr(2, Nu2, nu2, 512, 64u, NFFT_MEASURE);
  }

  /* Force each native 1D variant -- NO_FAST_NATIVE keeps the fast from
   * winning, one NO_NDFT_* flag pins the survivor -- and validate both
   * directions against the legacy direct oracle. */
  check_case_against_direct(1, 256, 512, 4096, 51u,
                            NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE | NFFT_NO_NDFT_BLOCKED); /* plain, estimate */
  check_case_against_direct(1, 256, 512, 4096, 52u,
                            NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE | NFFT_NO_NDFT_PLAIN); /* blocked, estimate */
  check_case_against_direct(1, 128, 256, 1024, 53u,
                            NFFT_MEASURE | NFFT_NO_FAST_NATIVE | NFFT_NO_NDFT_BLOCKED); /* plain, measured */
  check_case_against_direct(1, 128, 256, 1024, 54u,
                            NFFT_MEASURE | NFFT_NO_FAST_NATIVE | NFFT_NO_NDFT_PLAIN); /* blocked, measured */
  check_case_against_direct(1, 16, 64, 256, 55u,
                            NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE | NFFT_NO_NDFT_PLAIN); /* Ntot<=B blocked adjoint path */

  /* Force the generic multivariate NDFT: at d=2/3 NO_FAST_NATIVE declines the
   * fast, and at d=4 a tiny M makes the direct native win on cost. */
  check_case_against_direct(2, 32, 64, 512, 71u,
                            NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
  check_case_against_direct(2, 24, 48, 300, 72u,
                            NFFT_MEASURE | NFFT_NO_FAST_NATIVE);
  check_case_against_direct(3, 12, 24, 200, 73u,
                            NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
  check_case_against_direct(3, 10, 20, 128, 74u,
                            NFFT_MEASURE | NFFT_NO_FAST_NATIVE);
  check_case_against_direct(4, 8, 16, 4, 75u, NFFT_ESTIMATE);
  check_case_against_direct(4, 6, 12, 4, 76u, NFFT_MEASURE);

  Y(the_planner_destroy)
  ();
}

void Y(check_nplan_wisdom_memo)(void) {
  INT N = 256, n = 512, M = 4096;
  static R x[4096];
  static C f_hat[256], f[4096];
  /* 512: printer_create_str is unbounded, and the bundle print nests both
   * children plus the FFTW plan description, so a smaller buffer overflows
   * the stack instead of truncating. */
  char b1[512], b2[512];
  printer *pr;
  unsigned nelem;
  Y(plan_ng) * p1, *p2;

  fill_nodes(x, 1, M, 700u);

  p1 = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, Y(get_window_id)(), x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p1);
  nelem = Y(the_planner)()->htab_unblessed.nelem;
  /* estimate flows never bless: blessing is measured-mode only */
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_blessed.nelem, 0u);
  pr = Y(printer_create_str)(b1);
  Y(plan_ng_print)
  (p1, pr);
  Y(printer_destroy)
  (pr);
  CU_ASSERT_PTR_NOT_NULL(strstr(b1, "(fwd (nfft_solver_fast_native"));
  /* forward-only race -- the adjoint half is dormant (adj (null)). */
  CU_ASSERT_PTR_NOT_NULL(strstr(b1, "(adj (null))"));

  p2 = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, Y(get_window_id)(), x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p2);
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_unblessed.nelem, nelem);
  pr = Y(printer_create_str)(b2);
  Y(plan_ng_print)
  (p2, pr);
  Y(printer_destroy)
  (pr);
  CU_ASSERT_STRING_EQUAL(b1, b2);

  {
    /* build_core forces FFTW_PRESERVE_INPUT off whatever the caller passes,
     * so it must not fragment the wisdom key: keyable_fftw_flags strips it
     * before hashing, and a stray key would grow htab_unblessed.nelem. */
    Y(plan_ng) *p3 = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, Y(get_window_id)(), x, f_hat, f,
                                     FFTW_PRESERVE_INPUT, NFFT_ESTIMATE);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p3);
    CU_ASSERT_EQUAL(Y(the_planner)()->htab_unblessed.nelem, nelem);
    Y(plan_ng_destroy)
    (p3);
  }

  Y(plan_ng_destroy)
  (p1);
  Y(plan_ng_destroy)
  (p2);

  /* global teardown re-arms registration; the fresh generation really
   * re-searches (its wisdom starts empty and grows) */
  Y(the_planner_destroy)
  ();
  p1 = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, Y(get_window_id)(), x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p1);
  CU_ASSERT(Y(the_planner)()->htab_unblessed.nelem > 0u);
  Y(plan_ng_destroy)
  (p1);
  Y(the_planner_destroy)
  ();
}

/* Measured planning */

static char *export_to_string(planner *pl) {
  size_t cnt;
  char *s;
  printer *pr = Y(printer_create_cnt)(&cnt);
  Y(planner_export)
  (pl, pr);
  Y(printer_destroy)
  (pr);
  s = (char *)Y(malloc)(cnt + 1);
  pr = Y(printer_create_str)(s);
  Y(planner_export)
  (pl, pr);
  Y(printer_destroy)
  (pr);
  return s;
}

static int import_from_string(planner *pl, const char *s) {
  scanner *sc = Y(scanner_create_str)(s);
  int ret = Y(planner_import)(pl, sc);
  Y(scanner_destroy)
  (sc);
  return ret;
}

/* Extract the per-direction winning solver names from the bundle print
 * "(nfft-plan-ng (fwd (NAME pcost=...)) (adj (null)))". The measured winner is
 * machine- and thread-dependent, so measured tests must not hard-code one:
 * they capture the actual choice and assert only that fwd names a real solver
 * and that a repeat or imported plan reproduces it. Which solver wins at which
 * size is pinned by the estimate-mode check_nplan_solvers instead.
 * The race is forward-only, so adj is captured verbatim and reads "(null))". */
static void bundle_winners(Y(plan_ng) * p, char *fwd, char *adj) {
  /* 512: printer_create_str is unbounded, and the bundle print nests both
   * children plus the FFTW plan description, so a smaller buffer overflows
   * the stack instead of truncating. */
  char buf[512];
  const char *q;
  printer *pr = Y(printer_create_str)(buf);
  Y(plan_ng_print)
  (p, pr);
  Y(printer_destroy)
  (pr);
  q = strstr(buf, "(fwd (");
  CU_ASSERT_PTR_NOT_NULL_FATAL(q);
  sscanf(q + 6, "%95s", fwd);
  q = strstr(buf, "(adj (");
  CU_ASSERT_PTR_NOT_NULL_FATAL(q);
  sscanf(q + 6, "%95s", adj);
  /* only the forward direction selects a real solver; adj is dormant */
  CU_ASSERT_PTR_NOT_NULL(strstr(fwd, "nfft_solver_"));
}

/* The bundle print carries the forward direction's solver registrar name.
 * The race is forward-only, so (adj ...) always reads (null). */
void Y(check_nplan_print_includes_registrar_names)(void) {
  INT N = 64, n = 128, M = 8192;
  static R x[8192];
  static C f_hat[64], f[8192];
  Y(plan_ng) * p;
  /* 512: printer_create_str is unbounded, and the bundle print nests both
   * children plus the FFTW plan description, so a smaller buffer overflows
   * the stack instead of truncating. */
  char buf[512];
  printer *pr;

  Y(the_planner_destroy)
  ();
  fill_nodes(x, 1, M, 122u);

  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  pr = Y(printer_create_str)(buf);
  Y(plan_ng_print)
  (p, pr);
  Y(printer_destroy)
  (pr);

  CU_ASSERT_PTR_NOT_NULL(strstr(buf, "(fwd (nfft_solver_"));
  CU_ASSERT_PTR_NOT_NULL(strstr(buf, "(adj (null))"));

  Y(plan_ng_destroy)
  (p);
  Y(the_planner_destroy)
  ();
}

void Y(check_nplan_measured)(void) {
  INT N = 64, n = 128, M = 8192;
  static R x[8192];
  static C f_hat[64], f[8192];
  Y(plan_ng) * p, *p2;
  unsigned ub;
  char f0[96], a0[96], f1[96], a1[96];

  Y(the_planner_destroy)
  (); /* fresh store: absolute counts below */
  fill_nodes(x, 1, M, 51u);

  /* (a) A large-M measured race blesses exactly one size-class entry,
   * whatever wins, because only the forward problem is raced. The unblessed
   * table may still grow: the composed native fast memoises its DECONV/CONV
   * children on their first plan. Once those are memoised a repeat query is
   * idempotent -- neither table grows and the same winners come back. */
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  if (Y(the_planner)()->htab_blessed.nelem == 0u) {
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
    CU_FAIL("exhaustive mode requires a working cycle counter; measured path unverifiable");
#else
    CU_PASS("no usable cycle counter: measured degraded to estimate-grade (known-skipped)");
#endif
    /* no usable clock: measured degrades to estimate-grade selection, which
     * memoises unblessed only, and the bundle still works */
    CU_ASSERT(Y(the_planner)()->htab_unblessed.nelem > 0u);
    Y(precompute)
    (p);
    fill_fhat(f_hat, N, 52u);
    Y(execute)
    (p);
    Y(plan_ng_destroy)
    (p);
    Y(the_planner_destroy)
    ();
    return;
  }
  /* the forward direction raced and blessed (the adjoint half is dormant) */
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_blessed.nelem, 1u);
  ub = Y(the_planner)()->htab_unblessed.nelem; /* may be >0: DECONV/CONV children */
  bundle_winners(p, f0, a0);                   /* capture the actual per-direction winners */

  /* repeat measured guru: blessed hit, no race, no growth, same winners */
  p2 = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p2);
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_blessed.nelem, 1u);
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_unblessed.nelem, ub);
  bundle_winners(p2, f1, a1);
  CU_ASSERT_STRING_EQUAL(f1, f0);
  CU_ASSERT_STRING_EQUAL(a1, a0);
  Y(plan_ng_destroy)
  (p2);

  /* an estimate-mode guru sees the measured entries: the blessed hit answers
   * the estimate query, with no further growth and the same winners. A
   * blessed-hit lookup never reads x/f_hat/f, so the main arrays are reused. */
  p2 = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p2);
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_unblessed.nelem, ub);
  bundle_winners(p2, f1, a1);
  CU_ASSERT_STRING_EQUAL(f1, f0);
  CU_ASSERT_STRING_EQUAL(a1, a0);
  Y(plan_ng_destroy)
  (p2);
  Y(plan_ng_destroy)
  (p);

  /* (b) a tiny-M measured race on a fast-capable union core: both candidates
   * genuinely run and the planner blesses whichever measured fastest, so
   * assert only that a valid bundle came out */
  {
    INT Nd = 16, nd = 1024, Md = 1;
    R xd[1];
    C f_hatd[16], fd[1];
    fill_nodes(xd, 1, Md, 53u);
    p = Y(plan_ng_guru)(1, &Nd, 0, &nd, Md, 6, NFFT_WINDOW_KAISER_BESSEL, xd, f_hatd, fd, 0u, NFFT_MEASURE);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);
    bundle_winners(p, f1, a1); /* asserts both are real solver names */
    Y(plan_ng_destroy)
    (p);
  }

  /* a gated query misses the blessed {F=0} entries, so it searches and
   * memoises */
  ub = Y(the_planner)()->htab_unblessed.nelem;
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u,
                      NFFT_ESTIMATE | NFFT_NO_DIRECT);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  CU_ASSERT(Y(the_planner)()->htab_unblessed.nelem > ub); /* new memos */
  Y(plan_ng_destroy)
  (p);

  Y(the_planner_destroy)
  ();
}

/* Measured planning is destructive by default, as FFTW_MEASURE is: the
 * provisional race core aliases the caller's f_hat/f rather than allocating
 * scratch, and races on their real bytes un-zeroed. The race is forward-only,
 * so it reads f_hat and clobbers f; filling after planning is the documented
 * lifecycle. */
void Y(check_nplan_destructive_default)(void) {
  /* M small enough that the blocked NDFT survives the estimate gate
   * (PLNR_PRUNE_RATIO), so two candidates race. */
  INT N = 64, n = 128, M = 64;
  static R x[8192];
  static C f_hat[64], f[8192], saved[8192];
  Y(plan_ng) * p;
  NFFT(plan)
  ref;
  INT j;
  R err, tol;

  tol = (R)1.0e-10;
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;

  Y(the_planner_destroy)
  ();
  fill_nodes(x, 1, M, 91u);
  fill_fhat(f_hat, N, 92u); /* pre-fill before planning: read by the race */
  for (j = 0; j < M; j++)   /* pre-fill before planning: written by the race */
    f[j] = (C)((R)1 + (R)j);
  memcpy(saved, f, (size_t)M * sizeof(C));

  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);

  if (Y(the_planner)()->htab_blessed.nelem == 0u) {
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
    CU_FAIL("exhaustive mode requires a working cycle counter; measured path unverifiable");
#else
    CU_PASS("no usable cycle counter: measured degraded to estimate-grade (known-skipped)");
#endif
    /* No usable clock: no race ran and nothing was clobbered, so only the
     * fill-after lifecycle is exercised here. */
    Y(precompute)
    (p);
    fill_fhat(f_hat, N, 93u);
    Y(execute)
    (p);
    Y(plan_ng_destroy)
    (p);
    Y(the_planner_destroy)
    ();
    return;
  }

  /* The provisional core aliases the real f and the trafo candidates wrote
   * it, so the pre-fill is gone. This needs >= 2 raced candidates: a lone
   * candidate is adopted untimed, without ever running apply. */
  CU_ASSERT_NOT_EQUAL(memcmp(f, saved, (size_t)M * sizeof(C)), 0);

  /* Correct when f_hat is (re)filled after planning. */
  Y(precompute)
  (p);
  fill_fhat(f_hat, N, 93u);
  Y(execute)
  (p);
  {
    int Ni = (int)N, ni = (int)n;
    NFFT(init_guru)
    (&ref, 1, &Ni, (int)M, &ni, 6, MALLOC_X | MALLOC_F_HAT | MALLOC_F, 0u);
  }
  for (j = 0; j < M; j++)
    ref.x[j] = x[j];
  for (j = 0; j < N; j++)
    ref.f_hat[j] = f_hat[j];
  NFFT(trafo_direct)
  (&ref);
  err = rel_max_err(f, ref.f, M);
  CU_ASSERT(err < tol);
  NFFT(finalize)
  (&ref);

  Y(plan_ng_destroy)
  (p);
  Y(the_planner_destroy)
  ();
}

/* New-array execute: nfft_execute_on / _adjoint_on run the plan on
 * caller-supplied f_hat/f, with x fixed at plan time, without disturbing the
 * plan-time bindings. */
void Y(check_nplan_execute_on)(void) {
  INT N = 64, n = 128, M = 512;
  static R x[512];
  static C f_hat[64], f[512];     /* plan-time bound arrays */
  static C f_hat_b[64], f_b[512]; /* independent new arrays */
  static C f_b_snap[512];
  Y(plan_ng) * p;
  NFFT(plan)
  ref;
  INT j;
  R err, tol;
  int Ni, ni, w;
  R sigma;

  /* window-aware tolerance from err_bound; sigma = n/N = 2.0 */
  w = Y(get_window_id)();
  sigma = ((R)n) / ((R)N);
  tol = err_bound(w, (R)6, sigma);
  Ni = (int)N;
  ni = (int)n;

  Y(the_planner_destroy)
  ();
  fill_nodes(x, 1, M, 81u);
  /* Request the compile-time window: under a non-KB build a hardcoded KB
   * request runs the native with the wrong window, and the bound above no
   * longer applies. */
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, Y(get_window_id)(), x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);

  /* (a) forward new-array execute on f_hat_b -> f_b, not the bound arrays */
  fill_fhat(f_hat_b, N, 82u);
  Y(execute_on)
  (p, f_hat_b, f_b);
  NFFT(init_guru)
  (&ref, 1, &Ni, (int)M, &ni, 6, MALLOC_X | MALLOC_F_HAT | MALLOC_F, 0u);
  for (j = 0; j < M; j++)
    ref.x[j] = x[j];
  for (j = 0; j < N; j++)
    ref.f_hat[j] = f_hat_b[j];
  NFFT(trafo_direct)
  (&ref);
  err = rel_max_err(f_b, ref.f, M);
  CU_ASSERT(err < tol);
  memcpy(f_b_snap, f_b, (size_t)M * sizeof(C));

  /* (b) plan-time bound execute is independent: nfft_execute(p) writes the
   * bound f and leaves f_b untouched */
  fill_fhat(f_hat, N, 777u);
  Y(execute)
  (p);
  for (j = 0; j < N; j++)
    ref.f_hat[j] = f_hat[j];
  NFFT(trafo_direct)
  (&ref);
  err = rel_max_err(f, ref.f, M);
  CU_ASSERT(err < tol);
  /* f_b from the _on call must be untouched by nfft_execute(p) */
  CU_ASSERT_EQUAL(memcmp(f_b, f_b_snap, (size_t)M * sizeof(C)), 0);

  /* (c) adjoint new-array execute: reads f_b, writes f_hat_b */
  Y(execute_adjoint_on)
  (p, f_hat_b, f_b);
  for (j = 0; j < M; j++)
    ref.f[j] = f_b[j];
  NFFT(adjoint_direct)
  (&ref);
  err = rel_max_err(f_hat_b, ref.f_hat, N);
  CU_ASSERT(err < tol);

  NFFT(finalize)
  (&ref);
  Y(plan_ng_destroy)
  (p);
  Y(the_planner_destroy)
  ();
}

void Y(check_nplan_measured_wisdom)(void) {
  INT N = 64, n = 128, M = 8192;
  static R x[8192];
  static C f_hat[64], f[8192];
  Y(plan_ng) * p;
  char *wis;
  char f0[96], a0[96], f1[96], a1[96];

  Y(the_planner_destroy)
  ();
  fill_nodes(x, 1, M, 61u);

  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  if (Y(the_planner)()->htab_blessed.nelem == 0u) {
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
    CU_FAIL("exhaustive mode requires a working cycle counter; measured path unverifiable");
#else
    CU_PASS("no usable cycle counter: measured degraded to estimate-grade (known-skipped)");
#endif
    Y(plan_ng_destroy)
    (p);
    Y(the_planner_destroy)
    ();
    return;
  }
  /* remember the per-direction winners chosen on this machine */
  bundle_winners(p, f0, a0);
  wis = export_to_string(Y(the_planner)());
  Y(plan_ng_destroy)
  (p);

  /* fresh process simulation: new generation, import, re-plan measured */
  Y(the_planner_destroy)
  ();
  Y(nfft_ensure_registered)
  (); /* import checks the config signature */
  CU_ASSERT(import_from_string(Y(the_planner)(), wis));
  Y(free)
  (wis);
  /* the race is forward-only, so only the sign=+1 entry is exported */
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_blessed.nelem, 1u);

  /* The measured guru adopts the imported blessed winner: the top-level
   * problem is a blessed hit, so the same forward choice comes back. Its
   * DECONV/CONV children are never blessed or exported, so they re-race on
   * the imported re-plan and each leaves an unblessed memo. That is the
   * native winner rebuilding its children in a fresh generation, not a wisdom
   * leak: the count is stable across the cycle. */
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_blessed.nelem, 1u);
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_unblessed.nelem, 2u);
  bundle_winners(p, f1, a1);
  CU_ASSERT_STRING_EQUAL(f1, f0);
  CU_ASSERT_STRING_EQUAL(a1, a0);
  Y(plan_ng_destroy)
  (p);

  Y(the_planner_destroy)
  ();
}

/* A zero-second timelimit forces every measured race to bail out before a
 * measurement completes, so the guru takes the same estimate-grade restart
 * path as the no-usable-clock case. Estimate never blesses, so the blessed
 * store stays empty; the degraded bundle still produces unblessed memos. */
void Y(check_nplan_timelimit_tight_degrades_to_estimate)(void) {
  /* M small enough that the blocked NDFT survives the estimate gate
   * (PLNR_PRUNE_RATIO), so two candidates race. */
  INT N = 64, n = 128, M = 64;
  static R x[8192];
  static C f_hat[64], f[8192];
  Y(plan_ng) * p;

  Y(the_planner_destroy)
  (); /* fresh store: absolute counts below */
  fill_nodes(x, 1, M, 81u);

  Y(set_timelimit)
  (0.0); /* tight: every candidate race bails at once */
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  /* No candidate finished before expiry: degraded to estimate-grade
   * selection; estimate never blesses. */
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_blessed.nelem, 0u);
  /* Estimate-grade selection memoises unblessed entries. */
  CU_ASSERT(Y(the_planner)()->htab_unblessed.nelem > 0u);
  Y(plan_ng_destroy)
  (p);

  Y(set_timelimit)
  (-1.0); /* restore unlimited for the rest of the suite */
  Y(the_planner_destroy)
  ();
}

/* With the default unlimited timelimit (-1.0) the race completes and blesses
 * one new entry for the forward direction. */
void Y(check_nplan_timelimit_unset_measures_and_blesses)(void) {
  INT N = 64, n = 128, M = 8192;
  static R x[8192];
  static C f_hat[64], f[8192];
  Y(plan_ng) * p;

  Y(the_planner_destroy)
  (); /* fresh store: absolute counts below */
  fill_nodes(x, 1, M, 82u);

  CK(Y(planner_timelimit)(Y(the_planner)()) == -1.0);
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  if (Y(the_planner)()->htab_blessed.nelem == 0u) {
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
    CU_FAIL("exhaustive mode requires a working cycle counter; measured path unverifiable");
#else
    CU_PASS("no usable cycle counter: measured degraded to estimate-grade (known-skipped)");
#endif
    /* no usable clock on this platform: degradation, not timelimit. */
    Y(plan_ng_destroy)
    (p);
    Y(the_planner_destroy)
    ();
    return;
  }
  /* the forward direction's race finished and blessed. */
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_blessed.nelem, 1u);
  Y(plan_ng_destroy)
  (p);

  Y(the_planner_destroy)
  ();
}

/* The public nfft_set_timelimit symbol forwards to the internal global-planner
 * timelimit (one nfft_ symbol serves all three kinds). */
void Y(check_nplan_set_timelimit_roundtrip)(void) {
  Y(set_timelimit)
  (0.5);
  CK(Y(planner_timelimit)(Y(the_planner)()) == 0.5);
  Y(set_timelimit)
  (-1.0); /* restore unlimited for the rest of the suite */
}

/* Exercise the public surface (nfft3.h) through the Y() mangle. */
void Y(check_nplan_public_api)(void) {
  INT N = 64, n = 128, M = 512;
  static R x[512];
  static C f_hat[64], f[512];
  Y(plan_ng) * p;
  char *s;
  FILE *tmp;
  char buf[256];
  size_t got;
  unsigned blessed;

  Y(the_planner_destroy)
  ();
  fill_nodes(x, 1, M, 71u);
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);

  /* fprint: the plan tree lands in the stream and names a registrar */
  tmp = tmpfile();
  CU_ASSERT_PTR_NOT_NULL_FATAL(tmp);
  Y(fprint_plan)
  (p, tmp);
  rewind(tmp);
  got = fread(buf, 1, sizeof(buf) - 1, tmp);
  buf[got] = '\0';
  fclose(tmp);
  CU_ASSERT(got > 0);
  CU_ASSERT_PTR_NOT_NULL(strstr(buf, "nfft_solver_"));

  blessed = Y(the_planner)()->htab_blessed.nelem; /* 1 with a clock, else 0 */

  /* string wisdom roundtrip through the public functions */
  s = Y(export_wisdom_to_string)();
  CU_ASSERT_PTR_NOT_NULL_FATAL(s);
  Y(forget_wisdom)
  ();
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_blessed.nelem, 0u);
  CU_ASSERT(Y(import_wisdom_from_string)(s));
  Y(free)
  (s);
  CU_ASSERT_EQUAL(Y(the_planner)()->htab_blessed.nelem, blessed);

  /* filename wisdom roundtrip */
  {
    const char *path = "nplan_public_wisdom.tmp";
    CU_ASSERT(Y(export_wisdom_to_filename)(path));
    Y(forget_wisdom)
    ();
    CU_ASSERT_EQUAL(Y(the_planner)()->htab_blessed.nelem, 0u);
    CU_ASSERT(Y(import_wisdom_from_filename)(path));
    CU_ASSERT_EQUAL(Y(the_planner)()->htab_blessed.nelem, blessed);
    remove(path);
  }

  Y(plan_ng_destroy)
  (p);
  Y(the_planner_destroy)
  ();
}

/* Forward-only race: the bundle selects one plan (dir[FWD]) and the adjoint
 * reuses it via apply_adjoint, so plan_ng_print always shows a real forward
 * solver and a (null) adjoint, in both planning modes. */
void Y(check_nplan_forward_only_race)(void) {
  INT N = 64, n = 128, M = 2048;
  static R x[2048];
  static C f_hat[64], f[2048];
  Y(plan_ng) * p;
  /* 512: printer_create_str is unbounded, and the bundle print nests both
   * children plus the FFTW plan description, so a smaller buffer overflows
   * the stack instead of truncating. */
  char buf[512];
  printer *pr;

  Y(the_planner_destroy)
  ();
  fill_nodes(x, 1, M, 501u);

  /* estimate mode: one forward plan, adjoint dormant */
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  pr = Y(printer_create_str)(buf);
  Y(plan_ng_print)
  (p, pr);
  Y(printer_destroy)
  (pr);
  CU_ASSERT_PTR_NOT_NULL(strstr(buf, "(fwd (nfft_solver_"));
  CU_ASSERT_PTR_NOT_NULL(strstr(buf, "(adj (null))"));
  Y(plan_ng_destroy)
  (p);

  /* measured mode: same single-plan shape */
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  pr = Y(printer_create_str)(buf);
  Y(plan_ng_print)
  (p, pr);
  Y(printer_destroy)
  (pr);
  CU_ASSERT_PTR_NOT_NULL(strstr(buf, "(fwd (nfft_solver_"));
  CU_ASSERT_PTR_NOT_NULL(strstr(buf, "(adj (null))"));
  Y(plan_ng_destroy)
  (p);
  Y(the_planner_destroy)
  ();
}

/* execute_adjoint on the forward-winning plan produces the correct adjoint;
 * wrong wiring would run the forward transform and mismatch the reference.
 * Measured, so a fast solver wins and its apply_adjoint is exercised. */
void Y(check_nplan_apply_adjoint)(void) {
  INT N = 64, n = 128, M = 2048, Ntot = N;
  R *xin;
  C *f_hat, *f;
  Y(plan_ng) * p;
  NFFT(plan)
  ref;
  R err, tol;

  tol = (R)1.0e-10;
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;

  Y(the_planner_destroy)
  ();
  xin = (R *)Y(malloc)((size_t)M * sizeof(R));
  fill_nodes(xin, 1, M, 511u);
  f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  f = (C *)Y(malloc)((size_t)M * sizeof(C));

  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, xin, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);

  {
    int Ni = (int)N, ni = (int)n;
    NFFT(init_guru)
    (&ref, 1, &Ni, (int)M, &ni, 6,
     MALLOC_X | MALLOC_F_HAT | MALLOC_F, 0u);
  }
  memcpy(ref.x, xin, (size_t)M * sizeof(R));
  fill_fhat(ref.f, M, 512u);
  memcpy(f, ref.f, (size_t)M * sizeof(C));
  Y(execute_adjoint)
  (p);
  NFFT(adjoint_direct)
  (&ref);
  err = rel_max_err(f_hat, ref.f_hat, Ntot);
  CU_ASSERT(err < tol);

  NFFT(finalize)
  (&ref);
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (xin);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* plan_ng_guru copies the caller's x, so the plan runs on the value x had at
 * guru time and a mutation before precompute/execute stays invisible. This is
 * the assertion that separates the copy contract from an aliasing one. N <= m
 * makes the fast solver decline, which is incidental: x-copying is a
 * problem-level property. */
void Y(check_nplan_x_copied_not_aliased)(void) {
  INT N = 4, n = 16, M = 1;
  R x[1];
  C f_hat[4], f[1];
  Y(plan_ng) * p;
  NFFT(plan)
  ref;
  R err, tol;

  tol = (R)1.0e-10;
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;

  Y(the_planner_destroy)
  ();
  x[0] = (R)0.0; /* the value the plan copies */
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);

  /* Mutate x after guru, before precompute: the plan already has its copy,
   * so this write must stay invisible to the transform. */
  x[0] = (R)0.3125;

  Y(precompute)
  (p);
  fill_fhat(f_hat, N, 998u);
  Y(execute)
  (p);

  {
    int Ni = (int)N, ni = (int)n;
    NFFT(init_guru)
    (&ref, 1, &Ni, (int)M, &ni, 6,
     MALLOC_X | MALLOC_F_HAT | MALLOC_F, 0u);
  }
  ref.x[0] = (R)0.0; /* the guru-time value, not the post-guru clobber */
  memcpy(ref.f_hat, f_hat, (size_t)N * sizeof(C));
  NFFT(trafo_direct)
  (&ref);
  err = rel_max_err(f, ref.f, M);
  CU_ASSERT(err < tol);

  NFFT(finalize)
  (&ref);
  Y(plan_ng_destroy)
  (p);
  Y(the_planner_destroy)
  ();
}

/* A wrapper plan's core carries PRE_PSI and builds psi once per awake period;
 * a native winner stays coreless. */
void Y(check_nplan_per_plan_core)(void) {
  INT N = 256, n = 512, M = 4096;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  Y(plan_ng) * p;
  fill_nodes(x, 1, M, 700u);
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);
  Y(execute)
  (p);
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
}

/* No plan allocates its own data arrays: f_hat/f are caller-owned and aliased,
 * and x is aliased to the problem's own copy, so build_core
 * (kernel/nfft/nsolver.c) must never request MALLOC_X | MALLOC_F_HAT |
 * MALLOC_F. The winner here is coreless, so there are no core flags left to
 * inspect and this only covers planning such a problem.
 * TODO: assert the core flags on a plan that builds a core. */
void Y(check_nplan_core_owns_no_data_arrays)(void) {
  INT N = 64, n = 128, M = 512;
  static R x[512];
  static C f_hat[64], f[512];
  Y(plan_ng) * p;

  Y(the_planner_destroy)
  ();
  fill_nodes(x, 1, M, 731u);
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);

  Y(plan_ng_destroy)
  (p);
  Y(the_planner_destroy)
  ();
}

/* The blocked recurrence is strictly more accurate than plain per-term in the
 * ill-conditioned regime (large N, nodes near +-1/2), where plain's error is
 * dominated by libm reducing huge phase arguments. That loss is present in
 * every precision. The reference sums in long double with argument reduction,
 * independent of R. */
void Y(check_nplan_ndft_accuracy)(void) {
  enum {
    N = 8192,
    M = 64
  };
  const INT Nn = N, nn = 2 * N, Mm = M;
  static R x[M];
  static C fh[N], f_plain[M], f_blocked[M];
  long double _Complex ref[M];
  const long double TWO_PI = 6.28318530717958647692528676655900577L;
  Y(plan_ng) * pp;
  INT j, k;
  long double num_p = 0.0L, num_b = 0.0L, den = 0.0L;
  R err_plain, err_blocked;

  Y(the_planner_destroy)
  ();

  for (j = 0; j < Mm; j++)
    x[j] = (R)0.4999 - (R)0.0001 * (R)((int)j % 7); /* near +1/2: worst for plain */
  fill_fhat(fh, Nn, 4242u);

  for (j = 0; j < Mm; j++) {
    long double _Complex v = 0.0L;
    long double xj = (long double)x[j];
    for (k = 0; k < Nn; k++) {
      long double kk = (long double)(k - Nn / 2);
      long double t = kk * xj;
      long double r = t - rintl(t);
      long double om = TWO_PI * r;
      v += ((long double _Complex)fh[k]) * (cosl(om) - _Complex_I * sinl(om));
    }
    ref[j] = v;
  }

  /* Without NFFT_NO_FAST_NATIVE the composed fast would outrun the direct
   * natives here and return an approximate result. */
  pp = Y(plan_ng_guru)(1, &Nn, 0, &nn, Mm, 6, NFFT_WINDOW_KAISER_BESSEL, x, fh, f_plain, 0u,
                       NFFT_ESTIMATE | NFFT_NO_NDFT_BLOCKED |
                           NFFT_NO_FAST_NATIVE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pp);
  Y(precompute)
  (pp);
  Y(execute)
  (pp);
  Y(plan_ng_destroy)
  (pp);

  pp = Y(plan_ng_guru)(1, &Nn, 0, &nn, Mm, 6, NFFT_WINDOW_KAISER_BESSEL, x, fh, f_blocked, 0u,
                       NFFT_ESTIMATE | NFFT_NO_NDFT_PLAIN |
                           NFFT_NO_FAST_NATIVE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pp);
  Y(precompute)
  (pp);
  Y(execute)
  (pp);
  Y(plan_ng_destroy)
  (pp);

  for (j = 0; j < Mm; j++) {
    long double dp = cabsl((long double _Complex)f_plain[j] - ref[j]);
    long double db = cabsl((long double _Complex)f_blocked[j] - ref[j]);
    long double da = cabsl(ref[j]);
    if (dp > num_p)
      num_p = dp;
    if (db > num_b)
      num_b = db;
    if (da > den)
      den = da;
  }
  err_plain = (R)(den > 0.0L ? num_p / den : num_p);
  err_blocked = (R)(den > 0.0L ? num_b / den : num_b);

  /* Blocked strictly more accurate than plain, in every build precision. */
  CU_ASSERT(err_blocked < err_plain);

  Y(the_planner_destroy)
  ();
}

/* Core elision: an NDFT-only bundle needs no legacy core, so it builds none
 * -- no FFTW planning, no phi_hut, no g1/g2 -- and still executes. */
void Y(check_nplan_core_elision)(void) {
  Y(the_planner_destroy)
  ();

  /* (a) tiny M -> blocked NDFT beats fast; coreless; still executes. */
  {
    INT N = 256, n = 1024, M = 1;
    static R x[1];
    static C f_hat[256], f[1];
    Y(plan_ng) * p;
    NFFT(plan)
    ref;
    R err, tol;
    INT j;

    tol = (R)1.0e-10;
    if ((R)1.0e4 * EPSILON > tol)
      tol = (R)1.0e4 * EPSILON;

    fill_nodes(x, 1, M, 11u);
    p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_ESTIMATE);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    Y(precompute)
    (p);
    fill_fhat(f_hat, N, 111u);
    Y(execute)
    (p);
    {
      int Ni = (int)N, ni = (int)n;
      NFFT(init_guru)
      (&ref, 1, &Ni, (int)M, &ni, 6,
       MALLOC_X | MALLOC_F_HAT | MALLOC_F, 0u);
    }
    for (j = 0; j < M; j++)
      ref.x[j] = x[j];
    for (j = 0; j < N; j++)
      ref.f_hat[j] = f_hat[j];
    NFFT(trafo_direct)
    (&ref);
    err = rel_max_err(f, ref.f, M);
    CU_ASSERT(err < tol);
    NFFT(finalize)
    (&ref);
    Y(plan_ng_destroy)
    (p);
  }

  /* (c) a pure multivariate NDFT-only bundle is coreless too. */
  {
    INT N[2] = {32, 32}, n[2] = {64, 64};
    INT M = 4;
    static R x[8];
    static C f_hat[1024], f[4];
    Y(plan_ng) * p;
    NFFT(plan)
    ref;
    R err, tol;
    INT j;

    tol = (R)1.0e-10;
    if ((R)1.0e4 * EPSILON > tol)
      tol = (R)1.0e4 * EPSILON;

    fill_nodes(x, 2, M, 21u);
    p = Y(plan_ng_guru)(2, N, 0, n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u,
                        NFFT_ESTIMATE);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    Y(precompute)
    (p);
    fill_fhat(f_hat, 1024, 121u);
    Y(execute)
    (p);
    {
      int Ni[2] = {32, 32}, ni[2] = {64, 64};
      NFFT(init_guru)
      (&ref, 2, Ni, (int)M, ni, 6, MALLOC_X | MALLOC_F_HAT | MALLOC_F, 0u);
    }
    for (j = 0; j < 2 * M; j++)
      ref.x[j] = x[j];
    for (j = 0; j < 1024; j++)
      ref.f_hat[j] = f_hat[j];
    NFFT(trafo_direct)
    (&ref);
    err = rel_max_err(f, ref.f, M);
    CU_ASSERT(err < tol);
    NFFT(finalize)
    (&ref);
    Y(plan_ng_destroy)
    (p);
  }

  Y(the_planner_destroy)
  ();
}

/* The guru accepts a per-axis variant and builds an executable type-II
 * plan (numerics checked in check_nplan_type_ii_1d). */
void Y(check_nplan_variant_guru)(void) {
  INT N = 16, n = 32, M = 100, j;
  int variant = NFFT_NDFT_TYPE_II;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  Y(plan_ng) * p;
  for (j = 0; j < M; j++)
    x[j] = K(0.3) * ((R)j / (R)M) - K(0.15);
  p = Y(plan_ng_guru)(1, &N, &variant, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(plan_ng_destroy)
  (p);
  Y(the_planner_destroy)
  ();
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
}

/* The plan tree must name the given solver, so a test that means to exercise
 * the fast path fails loudly if a direct solver served it instead. */
static void assert_plan_names(Y(plan_ng) * p, const char *needle)
{
  FILE *tmp = tmpfile();
  char buf[1024];
  size_t got;
  CU_ASSERT_PTR_NOT_NULL_FATAL(tmp);
  Y(fprint_plan)
  (p, tmp);
  rewind(tmp);
  got = fread(buf, 1, sizeof(buf) - 1, tmp);
  buf[got] = '\0';
  fclose(tmp);
  CU_ASSERT_PTR_NOT_NULL(strstr(buf, needle));
}

/* One geometry through the composed fast solver, forced with NFFT_NO_DIRECT
 * and checked against the exact NDFT reference at the window-aware err_bound
 * (m = 6, sigma = min axis oversampling), widened 4x for -ffast-math. */
static void fast_case(int d, const INT *N, const int *variant, const INT *n,
                      INT M, unsigned seed)
{
  INT Ntot = 1, j;
  int t;
  R *xin;
  C *f_hat, *f, *ref_f, *ref_fhat;
  Y(plan_ng) * p;
  R sigma = (R)0, tol, err;
  for (t = 0; t < d; t++) {
    R st = (R)n[t] / (R)N[t];
    if (sigma == (R)0 || st < sigma)
      sigma = st;
  }
  tol = K(4.0) * err_bound(NFFT_WINDOW_KAISER_BESSEL, (R)6, sigma);
  for (t = 0; t < d; t++)
    Ntot *= N[t];
  xin = (R *)Y(malloc)((size_t)(d * M) * sizeof(R));
  f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  f = (C *)Y(malloc)((size_t)M * sizeof(C));
  ref_f = (C *)Y(malloc)((size_t)M * sizeof(C));
  ref_fhat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  fill_nodes(xin, d, M, seed);
  p = Y(plan_ng_guru)(d, N, variant, n, M, 6, NFFT_WINDOW_KAISER_BESSEL, xin,
                      f_hat, f, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  assert_plan_names(p, "nfft_solver_fast_native");
  Y(precompute)
  (p);
  fill_fhat(f_hat, Ntot, seed + 1000u);
  Y(execute)
  (p);
  ndft_ref_trafo(d, N, variant, Ntot, M, xin, f_hat, ref_f);
  err = rel_max_err(f, ref_f, M);
  CU_ASSERT(err <= tol);
  for (j = 0; j < M; j++)
    f[j] = ref_f[j];
  Y(execute_adjoint)
  (p);
  ndft_ref_adjoint(d, N, variant, Ntot, M, xin, ref_f, ref_fhat);
  err = rel_max_err(f_hat, ref_fhat, Ntot);
  CU_ASSERT(err <= tol);
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (xin);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(free)
  (ref_f);
  Y(free)
  (ref_fhat);
}

/* Odd per-axis N for the nD native, forced through the direct native with
 * NFFT_NO_FAST_NATIVE and validated against ndft_ref_trafo / ndft_ref_adjoint
 * above rather than the legacy oracle, which carries its own copy of the
 * odometer bug. */
static void odd_case(int d, const INT *N, const INT *n, INT M, unsigned seed) {
  INT Ntot = 1, j;
  int t;
  R *xin;
  C *f_hat, *f, *ref_f, *ref_fhat;
  Y(plan_ng) * p;
  R tol = (R)1.0e-10;
  R err;
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;
  for (t = 0; t < d; t++)
    Ntot *= N[t];
  xin = (R *)Y(malloc)((size_t)(d * M) * sizeof(R));
  f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  f = (C *)Y(malloc)((size_t)M * sizeof(C));
  ref_f = (C *)Y(malloc)((size_t)M * sizeof(C));
  ref_fhat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  fill_nodes(xin, d, M, seed);
  p = Y(plan_ng_guru)(d, N, 0, n, M, 6, NFFT_WINDOW_KAISER_BESSEL, xin, f_hat,
                      f, 0u, NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);
  fill_fhat(f_hat, Ntot, seed + 1000u);
  Y(execute)
  (p);
  ndft_ref_trafo(d, N, 0, Ntot, M, xin, f_hat, ref_f);
  err = rel_max_err(f, ref_f, M);
  CU_ASSERT(err < tol);
  for (j = 0; j < M; j++)
    f[j] = ref_f[j];
  Y(execute_adjoint)
  (p);
  ndft_ref_adjoint(d, N, 0, Ntot, M, xin, ref_f, ref_fhat);
  err = rel_max_err(f_hat, ref_fhat, Ntot);
  CU_ASSERT(err < tol);
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (xin);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(free)
  (ref_f);
  Y(free)
  (ref_fhat);
}

void Y(check_nplan_odd_n)(void) {
  {
    INT N[1] = {15}, n[1] = {32};
    odd_case(1, N, n, 500, 11u);
  }
  {
    INT N[1] = {25}, n[1] = {64};
    odd_case(1, N, n, 500, 12u);
  }
  {
    INT N[2] = {5, 10}, n[2] = {16, 20};
    odd_case(2, N, n, 400, 13u);
  }
  {
    INT N[2] = {10, 5}, n[2] = {20, 16};
    odd_case(2, N, n, 400, 14u);
  } /* inner odd */
  {
    INT N[2] = {5, 5}, n[2] = {16, 16};
    odd_case(2, N, n, 400, 15u);
  }
  {
    INT N[3] = {5, 1, 4}, n[3] = {16, 1, 16};
    odd_case(3, N, n, 200, 16u);
  } /* odd+unit */
  /* The same shapes through the composed fast solver. */
  {
    INT N[1] = {15}, n[1] = {32};
    fast_case(1, N, 0, n, 500, 41u);
  }
  {
    INT N[1] = {25}, n[1] = {64};
    fast_case(1, N, 0, n, 500, 42u);
  }
  {
    INT N[2] = {15, 10}, n[2] = {32, 24};
    fast_case(2, N, 0, n, 400, 43u);
  }
  {
    INT N[2] = {9, 9}, n[2] = {20, 20};
    fast_case(2, N, 0, n, 400, 44u);
  }
  {
    INT N[3] = {15, 10, 7}, n[3] = {32, 24, 16};
    fast_case(3, N, 0, n, 200, 45u);
  }
  {
    INT N[4] = {7, 7, 7, 7}, n[4] = {16, 16, 16, 16};
    fast_case(4, N, 0, n, 100, 46u); /* rank >= 4 odometer */
  }
  {
    INT N[3] = {15, 16, 9}, n[3] = {32, 32, 20};
    int variant[3] = {NFFT_NDFT_TYPE_I, NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I};
    fast_case(3, N, variant, n, 200, 47u); /* odd and type-II on one plan */
  }
  Y(the_planner_destroy)
  ();
}

/* type-II 1D (ascending, +N/2 at the last slot): matches the type-II
 * reference, satisfies the modulation identity f_II(x) = exp(-2pi i x) f_I(x),
 * and differs from type-I. Direct and fast share the 1e-10 tolerance. */
static void type_ii_1d_once(unsigned steer) {
  INT N = 32, n = 64, M = 400, Ntot = N, j;
  int v_ii = NFFT_NDFT_TYPE_II;
  R *xin = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  C *ref = (C *)Y(malloc)((size_t)M * sizeof(C));
  C *ref_i = (C *)Y(malloc)((size_t)M * sizeof(C));
  C *ref_fhat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  Y(plan_ng) * p;
  R tol = (R)1.0e-10;
  R err;
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;
  fill_nodes(xin, 1, M, 21u);
  p = Y(plan_ng_guru)(1, &N, &v_ii, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, xin, f_hat, f, 0u,
                      NFFT_ESTIMATE | steer);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);
  fill_fhat(f_hat, Ntot, 22u);
  Y(execute)
  (p);
  ndft_ref_trafo(1, &N, &v_ii, Ntot, M, xin, f_hat, ref);
  err = rel_max_err(f, ref, M);
  CU_ASSERT(err < tol);
  /* modulation identity and difference vs type-I */
  ndft_ref_trafo(1, &N, 0, Ntot, M, xin, f_hat, ref_i);
  {
    C *mod = (C *)Y(malloc)((size_t)M * sizeof(C));
    for (j = 0; j < M; j++)
      mod[j] = (COS(K2PI * xin[j]) - II * SIN(K2PI * xin[j])) * ref_i[j];
    CU_ASSERT(rel_max_err(ref, mod, M) < tol); /* f_II = e^{-2pi i x} f_I */
    Y(free)
    (mod);
  }
  CU_ASSERT(rel_max_err(ref, ref_i, M) > (R)1.0e-3); /* differs from type-I */
  for (j = 0; j < M; j++)
    f[j] = ref[j];
  Y(execute_adjoint)
  (p);
  ndft_ref_adjoint(1, &N, &v_ii, Ntot, M, xin, ref, ref_fhat);
  err = rel_max_err(f_hat, ref_fhat, Ntot);
  CU_ASSERT(err < tol);
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (xin);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(free)
  (ref);
  Y(free)
  (ref_i);
  Y(free)
  (ref_fhat);
}

void Y(check_nplan_type_ii_1d)(void) {
  type_ii_1d_once(NFFT_NO_FAST_NATIVE
                  | NFFT_NO_NDFT_BLOCKED); /* plain direct */
  type_ii_1d_once(NFFT_NO_FAST_NATIVE
                  | NFFT_NO_NDFT_PLAIN); /* blocked direct */
  type_ii_1d_once(NFFT_NO_DIRECT); /* composed fast */
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}

/* per-axis type-II in the nD native (mixed type-I/type-II axes), forced
 * through the given side with steer. The 1e-10 tolerance holds direct and
 * fast alike. */
static void type_ii_nd_once(unsigned steer)
{
  INT N[2] = {16, 16}, n[2] = {32, 32}, M = 400, Ntot = 256, j;
  int variant[2] = {NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I};
  R *xin = (R *)Y(malloc)((size_t)(2 * M) * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  C *ref = (C *)Y(malloc)((size_t)M * sizeof(C));
  C *ref_fhat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  Y(plan_ng) * p;
  R tol = (R)1.0e-10;
  R err;
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;
  fill_nodes(xin, 2, M, 31u);
  p = Y(plan_ng_guru)(2, N, variant, n, M, 6, NFFT_WINDOW_KAISER_BESSEL, xin,
                      f_hat, f, 0u, NFFT_ESTIMATE | steer);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);
  fill_fhat(f_hat, Ntot, 32u);
  Y(execute)
  (p);
  ndft_ref_trafo(2, N, variant, Ntot, M, xin, f_hat, ref);
  err = rel_max_err(f, ref, M);
  CU_ASSERT(err < tol);
  for (j = 0; j < M; j++)
    f[j] = ref[j];
  Y(execute_adjoint)
  (p);
  ndft_ref_adjoint(2, N, variant, Ntot, M, xin, ref, ref_fhat);
  err = rel_max_err(f_hat, ref_fhat, Ntot);
  CU_ASSERT(err < tol);
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (xin);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(free)
  (ref);
  Y(free)
  (ref_fhat);
}

void Y(check_nplan_type_ii_nd)(void)
{
  type_ii_nd_once(NFFT_NO_FAST_NATIVE); /* direct native */
  type_ii_nd_once(NFFT_NO_DIRECT); /* composed fast */
  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}

/* The 3-state wakefulness enum. SLEEPY < AWAKE_ZERO < AWAKE
 * so ">= PLNR_AWAKE_ZERO" reads "runnable" and "== PLNR_AWAKE" reads
 * "correct". */
void Y(check_nplan_awake_states)(void) {
  CU_ASSERT_EQUAL(PLNR_SLEEPY, 0);
  CU_ASSERT_EQUAL(PLNR_AWAKE_ZERO, 1);
  CU_ASSERT_EQUAL(PLNR_AWAKE, 2);
  CU_ASSERT(PLNR_SLEEPY < PLNR_AWAKE_ZERO && PLNR_AWAKE_ZERO < PLNR_AWAKE);
}

/* mkproblem_nfft copies x for the top-level problem, so the user's array is
 * independent of the plan as soon as guru returns: clobbering it afterwards
 * must not disturb precompute or execute. */
void Y(check_nplan_user_x_pristine)(void) {
  INT N = 32, n = 64, M = 128;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  R *snap = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  Y(plan_ng) * p;
  INT j;
  fill_nodes(x, 1, M, 4242u);
  for (j = 0; j < M; j++)
    snap[j] = x[j];
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  for (j = 0; j < M; j++)
    CU_ASSERT_EQUAL(x[j], snap[j]); /* never written */
  for (j = 0; j < M; j++)
    x[j] = (R)-1; /* clobber user array post-plan */
  Y(precompute)
  (p);
  Y(execute)
  (p); /* must still run on the copy, not the clobbered user x */
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (x);
  Y(free)
  (snap);
  Y(free)
  (f_hat);
  Y(free)
  (f); /* user frees own x */
}

/* Direct-drive the permuting solver through SLEEPY -> AWAKE_ZERO -> SLEEPY
 * with no planner and no race: x must be permuted at AWAKE_ZERO and restored
 * exactly on the way back. copy_x=1, so pr owns a private permutable copy and
 * the test frees its own array separately. */
void Y(check_nplan_awake_zero_restore)(void) {
  INT N = 16, n = 32, M = 8, j;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  R snap[8];
  problem *pr;
  plan *pl;
  fill_nodes(x, 1, M, 55u);
  for (j = 0; j < M; j++)
    snap[j] = x[j];
  pr = Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, x, /*copy_x=*/1, 0, 0);
  pl = Y(nfft_perm_test_mkplan)(pr); /* borrows pr->x (the copy) */
  Y(nfft_perm_break_restore) = 0;
  Y(plan_awake)
  (pl, PLNR_AWAKE_ZERO);
  {
    int any = 0;
    const R *xc = ((const problem_nfft *)pr)->x;
    for (j = 0; j < M; j++)
      if (xc[j] != snap[j])
        any = 1;
    CU_ASSERT(any);
  }
  Y(plan_awake)
  (pl, PLNR_SLEEPY);
  {
    const R *xc = ((const problem_nfft *)pr)->x;
    for (j = 0; j < M; j++)
      CU_ASSERT_EQUAL(xc[j], snap[j]); /* restored */
  }
  Y(plan_destroy)
  (pl);
  Y(problem_destroy)
  (pr); /* frees the copy */
  Y(free)
  (x); /* test still owns its own array */
}

/* The restore guard Y(nfft_x_verify), the same one wrapping the measured race
 * in plan_ng.c, detects a broken restore. Gated on NFFT_DEBUG because the
 * guard is an A() and exists only in debug builds. */
void Y(check_nplan_restore_guard_fires)(void) {
#ifdef NFFT_DEBUG
  INT N = 16, n = 32, M = 8;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  md5sig sig;
  problem *pr;
  plan *pl;
  fill_nodes(x, 1, M, 77u);
  pr = Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u, x, 1, 0, 0);
  Y(nfft_x_md5)
  (((const problem_nfft *)pr)->x, M, sig); /* snapshot the copy */
  pl = Y(nfft_perm_test_mkplan)(pr);
  Y(nfft_perm_break_restore) = 1;
  Y(plan_awake)
  (pl, PLNR_AWAKE_ZERO); /* permute */
  Y(plan_awake)
  (pl, PLNR_SLEEPY); /* restore skipped */
  CU_ASSERT_FALSE(Y(nfft_x_verify)(((const problem_nfft *)pr)->x, M, sig));
  Y(nfft_perm_break_restore) = 0; /* undo for a clean teardown */
  Y(plan_awake)
  (pl, PLNR_AWAKE_ZERO);
  Y(plan_awake)
  (pl, PLNR_SLEEPY);
  Y(plan_destroy)
  (pl);
  Y(problem_destroy)
  (pr);
  Y(free)
  (x);
#else
  CU_PASS("restore guard is A(...)-gated; build with ./configure --enable-debug");
#endif
}

/* A candidate whose pcost exceeds PLNR_PRUNE_RATIO times the cheapest estimate
 * is never timed. The slow test solver accepts every NFFT problem with pcost
 * 1e18, so it joins the candidate set, must be pruned, and must never win. */
void Y(check_nplan_measured_prunes_by_estimate)(void) {
  INT N = 64, n = 128, M = 64;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  char buf[512];
  printer *pr;
  Y(plan_ng) * p;

  Y(the_planner_destroy)();
  fill_nodes(x, 1, M, 90u);
  Y(nfft_ensure_registered)();
  Y(nfft_solver_slow_test_register)(Y(the_planner)());
  Y(nfft_slow_test_applies) = 0;

  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 4, NFFT_WINDOW_KAISER_BESSEL, x, f_hat,
                      f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  CU_ASSERT_EQUAL(Y(nfft_slow_test_applies), (INT)0);

  pr = Y(printer_create_str)(buf);
  Y(plan_ng_print)(p, pr);
  Y(printer_destroy)(pr);
  CU_ASSERT_PTR_NULL(strstr(buf, "slow_test"));

  Y(plan_ng_destroy)(p);
  Y(the_planner_destroy)();
  Y(free)(x);
  Y(free)(f_hat);
  Y(free)(f);
}

/* Register the permuting solver into the process-global planner alongside the
 * real roster: its mkplan accepts unconditionally, so every NFFT problem has
 * >= 2 candidates and the race loop, with its debug x-restore guard between
 * candidates, actually runs. No abort under --enable-debug proves the restore
 * is exact inside a real race, not just under direct drive. */
void Y(check_nplan_in_race_guard_passes)(void) {
  INT N = 16, n = 32, M = 8192;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  Y(plan_ng) * p;

  Y(the_planner_destroy)
  ();
  fill_nodes(x, 1, M, 88u);

  Y(nfft_ensure_registered)
  ();
  Y(nfft_solver_perm_test_register)
  (Y(the_planner)());

  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_MEASURE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);
  fill_fhat(f_hat, N, 89u);
  Y(execute)
  (p);
  Y(plan_ng_destroy)
  (p);

  /* Destroy the process-global planner so the test-only solver does not
   * leak into later tests (a fresh generation re-registers only the real
   * roster via nfft_ensure_registered). */
  Y(the_planner_destroy)
  ();

  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
}

/* AWAKE_ZERO is internal to planning: a fresh guru's winner is returned
 * SLEEPY, and precompute is what awakens it to AWAKE. */
void Y(check_nplan_awake_zero_internal)(void) {
  INT N = 32, n = 64, M = 128;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  Y(plan_ng) * p;
  fill_nodes(x, 1, M, 8u);
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, f, 0u, NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  CU_ASSERT_EQUAL(Y(plan_ng_test_awake_state)(p), PLNR_SLEEPY); /* winner SLEEPY */
  Y(precompute)
  (p);
  CU_ASSERT_EQUAL(Y(plan_ng_test_awake_state)(p), PLNR_AWAKE);
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
}

/* A unit axis is elided at construction, so the surviving axes go through the
 * fast algorithm. */
void Y(check_nplan_elides_unit_axes)(void) {
  const int d = 3;
  INT N[3] = {16, 1, 8}; /* middle axis unit */
  INT n[3] = {32, 1, 16};
  const INT M = 20;
  const int m = 2;
  INT j;
  int t;
  R *x = (R *)Y(malloc)((size_t)d * M * sizeof(R));
  C *fh = (C *)Y(malloc)((size_t)(N[0] * N[1] * N[2]) * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  Y(plan_ng) * p;
  for (j = 0; j < M; j++)
    for (t = 0; t < d; t++)
      x[j * d + t] = (R)(-0.4) + (R)((j * d + t) % 7) / (R)7.0;
  for (j = 0; j < N[0] * N[1] * N[2]; j++)
    fh[j] = (R)1.0 / (R)(j + 1);
  p = Y(plan_ng_guru)(d, N, 0, n, M, m, Y(get_window_id)(), x, fh, f, 0u,
                      NFFT_NO_DIRECT);
  /* A surviving unit axis fails guards_ok (N_t > m), so under NFFT_NO_DIRECT a
   * non-NULL plan proves the middle axis was elided. */
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (x);
  Y(free)
  (fh);
  Y(free)
  (f);
}

/* All axes unit -> a rank-0 problem, served by the ungated exact base case
 * nfft_solver_const_0d: forward broadcasts f_hat[0], adjoint reduces
 * sum_j f[j]. */
void Y(check_nplan_rank0_solver)(void) {
  const int d = 3;
  INT N[3] = {1, 1, 1}, n[3] = {1, 1, 1};
  const INT M = 12;
  const int m = 2;
  INT j;
  R *x = (R *)Y(malloc)((size_t)d * M * sizeof(R));
  C *fh = (C *)Y(malloc)(sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  Y(plan_ng) * p;
  C acc = K(0.0);
  for (j = 0; j < d * M; j++)
    x[j] = (R)0.1 * (R)(j % 5);
  fh[0] = (R)2.5 - II * (R)1.25;
  p = Y(plan_ng_guru)(d, N, 0, n, M, m, Y(get_window_id)(), x, fh, f, 0u,
                      NFFT_NO_DIRECT);
  /* only nfft_solver_const_0d can serve a rank-0 problem under NFFT_NO_DIRECT */
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);
  Y(execute)
  (p);
  for (j = 0; j < M; j++)
    CU_ASSERT_DOUBLE_EQUAL(CABS(f[j] - fh[0]), 0.0, 1e-12); /* broadcast */
  for (j = 0; j < M; j++)
    f[j] = (R)(j + 1) - II * (R)j;
  Y(execute_adjoint)
  (p);
  for (j = 0; j < M; j++)
    acc += (R)(j + 1) - II * (R)j;
  CU_ASSERT_DOUBLE_EQUAL(CABS(fh[0] - acc), 0.0, 1e-10); /* reduce */
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (x);
  Y(free)
  (fh);
  Y(free)
  (f);
}

/* Direct NDFT oracle over the full d-tensor (row-major, k_t in [-N_t/2,
 * N_t/2)). A unit axis contributes a single k=0 term of phase 1, so it drops
 * out. */
static void ndft_ref(int d, const INT *N, INT M, const R *x, const C *fh,
                     C *out) {
  INT Ntot = 1, j, l;
  int t;
  for (t = 0; t < d; t++)
    Ntot *= N[t];
  for (j = 0; j < M; j++) {
    C acc = K(0.0);
    for (l = 0; l < Ntot; l++) {
      INT rem = l;
      R ph = K(0.0);
      for (t = d - 1; t >= 0; t--) {
        INT kt = rem % N[t];
        rem /= N[t];
        ph += ((R)(kt - N[t] / 2)) * x[j * d + t];
      }
      {
        R w = K2PI * ph;
        acc += fh[l] * (COS(w) - II * SIN(w));
      }
    }
    out[j] = acc;
  }
}

/* Adjoint oracle: fh[l] = sum_j f[j] * exp(+2 pi i k.x_j). */
static void ndft_adj_ref(int d, const INT *N, INT M, const R *x, const C *f,
                         C *fh) {
  INT Ntot = 1, j, l;
  int t;
  for (t = 0; t < d; t++)
    Ntot *= N[t];
  for (l = 0; l < Ntot; l++) {
    C acc = K(0.0);
    for (j = 0; j < M; j++) {
      INT rem = l;
      R ph = K(0.0);
      for (t = d - 1; t >= 0; t--) {
        INT kt = rem % N[t];
        rem /= N[t];
        ph += ((R)(kt - N[t] / 2)) * x[j * d + t];
      }
      {
        R w = K2PI * ph;
        acc += f[j] * (COS(w) + II * SIN(w));
      }
    }
    fh[l] = acc;
  }
}

/* Both directions across every unit-axis position (leading, middle, trailing,
 * two-unit), against the direct NDFT oracle above under the window-aware
 * bound. With flags 0u a direct solver would serve these shapes, so the
 * accuracy check alone is not an elision gate; each shape is planned a second
 * time with NFFT_NO_DIRECT, where a non-NULL plan proves elision happened. */
void Y(check_nplan_unit_axis_correct)(void) {
  INT shapes[4][3] = {{1, 16, 8}, {16, 1, 8}, {16, 8, 1}, {16, 1, 1}};
  INT nshape[4][3] = {{1, 32, 16}, {32, 1, 16}, {32, 16, 1}, {32, 1, 1}};
  const int d = 3, m = 6;
  const INT M = 24;
  int s, t;
  INT j;
  int w = Y(get_window_id)();
  for (s = 0; s < 4; s++) {
    INT *N = shapes[s], *n = nshape[s], Ntot = N[0] * N[1] * N[2];
    R sigma = (R)2.0, tol = err_bound(w, (R)m, sigma), err;
    R *x = (R *)Y(malloc)((size_t)d * M * sizeof(R));
    C *fh = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
    C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
    C *ref = (C *)Y(malloc)((size_t)M * sizeof(C));
    C *fhref = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
    Y(plan_ng) * p;
    fill_nodes(x, d, M, (unsigned)(100u + s));
    fill_fhat(fh, Ntot, (unsigned)(200u + s));
    p = Y(plan_ng_guru)(d, N, 0, n, M, m, w, x, fh, f, 0u, 0u);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);
    Y(precompute)
    (p);
    Y(execute)
    (p);
    ndft_ref(d, N, M, x, fh, ref);
    err = rel_max_err(f, ref, M);
    /* -ffast-math reassociation pushes the intrinsic pipeline error to ~1.1x
     * the theoretical bound, so the assertions here widen it by 4x. */
    CU_ASSERT(err <= K(4.0) * tol); /* forward */
    /* The transpose-convolution round-off is larger than the forward's, so
     * the adjoint takes the same widened bound. */
    fill_fhat(f, M, (unsigned)(300u + s));
    Y(execute_adjoint)
    (p);
    ndft_adj_ref(d, N, M, x, f, fhref);
    err = rel_max_err(fh, fhref, Ntot);
    CU_ASSERT(err <= K(4.0) * tol); /* adjoint */
    /* the fast path engages only after elision */
    {
      Y(plan_ng) *pf = Y(plan_ng_guru)(d, N, 0, n, M, m, w, x, fh, f, 0u,
                                       NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL(pf);
      if (pf)
        Y(plan_ng_destroy)
      (pf);
    }
    Y(plan_ng_destroy)
    (p);
    Y(free)
    (x);
    Y(free)
    (fh);
    Y(free)
    (f);
    Y(free)
    (ref);
    Y(free)
    (fhref);
  }
}

/* New-array execute on a plan served by nfft_solver_fast_native, whose
 * DECONV/CONV children must forward the swapped problem pointers instead of
 * the ones cached at construction. check_nplan_execute_on plans a 1D shape
 * served by a different solver, so it does not reach this path. Compared
 * against ndft_ref / ndft_adj_ref under the window-aware bound. */
void Y(check_nplan_newarray_native_fast)(void) {
  /* This geometry passes guards_ok (kernel/nfft/nfft-nd.c), so the native fast
   * is admitted; NFFT_NO_DIRECT then forces it, since direct would otherwise
   * serve a shape this small. */
  INT N[2] = {16, 16}, n[2] = {32, 32}, M = 200, Ntot = N[0] * N[1], j;
  const int d = 2, m = 6;
  int w = Y(get_window_id)();
  R sigma = (R)2.0, tol = err_bound(w, (R)m, sigma), err;
  R *x = (R *)Y(malloc)((size_t)d * M * sizeof(R));
  C *fh = (C *)Y(malloc)((size_t)Ntot * sizeof(C)); /* plan-time bound arrays */
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  C *fh_b = (C *)Y(malloc)((size_t)Ntot * sizeof(C)); /* independent "new" */
  C *f_b = (C *)Y(malloc)((size_t)M * sizeof(C));
  C *ref_f = (C *)Y(malloc)((size_t)M * sizeof(C));
  C *ref_fh = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  Y(plan_ng) * p;
  fill_nodes(x, d, M, 611u);
  p = Y(plan_ng_guru)(d, N, 0, n, M, m, w, x, fh, f, 0u,
                      NFFT_NO_DIRECT);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);

  /* (a) forward new-array execute on fh_b -> f_b, not the bound arrays */
  fill_fhat(fh_b, Ntot, 612u);
  Y(execute_on)
  (p, fh_b, f_b);
  ndft_ref(d, N, M, x, fh_b, ref_f);
  err = rel_max_err(f_b, ref_f, M);
  CU_ASSERT(err <= tol); /* must match the oracle, not be all-zero */

  /* (b) adjoint new-array execute: reads f_b, writes fh_b. */
  Y(execute_adjoint_on)
  (p, fh_b, f_b);
  ndft_adj_ref(d, N, M, x, f_b, ref_fh);
  err = rel_max_err(fh_b, ref_fh, Ntot);
  CU_ASSERT(err <= K(4.0) * tol); /* adjoint: the same widened bound */

  Y(plan_ng_destroy)
  (p);
  Y(free)
  (x);
  Y(free)
  (fh);
  Y(free)
  (f);
  Y(free)
  (fh_b);
  Y(free)
  (f_b);
  Y(free)
  (ref_f);
  Y(free)
  (ref_fh);
  Y(the_planner_destroy)
  ();
}

/* New-array execute on a unit-axis plan: the unit axis is elided at
 * construction, so a replacement f_hat in the full d layout must read
 * correctly through the compressed strides. Compared against ndft_ref under
 * the window-aware bound. */
void Y(check_nplan_unit_axis_execute_on)(void) {
  const int d = 3;
  INT N[3] = {16, 1, 8}, n[3] = {32, 1, 16};
  const INT M = 20;
  const int m = 6;
  INT Ntot = N[0] * N[1] * N[2], j;
  int t;
  int w = Y(get_window_id)();
  R sigma = (R)2.0, tol = err_bound(w, (R)m, sigma), err;
  R *x = (R *)Y(malloc)((size_t)d * M * sizeof(R));
  C *fh0 = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *fh1 = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  C *ref = (C *)Y(malloc)((size_t)M * sizeof(C));
  Y(plan_ng) * p;
  for (j = 0; j < M; j++)
    for (t = 0; t < d; t++)
      x[j * d + t] = (R)(-0.4) + (R)((j * 3 + t) % 9) / (R)9.0;
  for (j = 0; j < Ntot; j++) {
    fh0[j] = K(0.0);
    fh1[j] = (R)1.0 / (R)(j + 3);
  }
  p = Y(plan_ng_guru)(d, N, 0, n, M, m, w, x, fh0, f, 0u, 0u);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  Y(precompute)
  (p);
  Y(execute_on)
  (p, fh1, f); /* run on the new f_hat */
  ndft_ref(d, N, M, x, fh1, ref);
  err = rel_max_err(f, ref, M);
  /* The residual here is the pipeline's own window-approximation error, not a
   * new-array defect: the decaying fh1 = 1/(j+3) lands it at ~0.8x the
   * theoretical KB bound. Random data is no safer, scattering around the
   * nominal tol, so this takes the same widened bound as the other cases. */
  CU_ASSERT(err <= K(4.0) * tol);
  Y(plan_ng_destroy)
  (p);
  Y(free)
  (x);
  Y(free)
  (fh0);
  Y(free)
  (fh1);
  Y(free)
  (f);
  Y(free)
  (ref);
}

/* The guru must return NULL rather than dereference on NULL or zero
 * arguments, including in release builds where A() is a no-op. */
void Y(check_nplan_guru_rejects_null_args)(void) {
  INT N = 4, n = 8, M = 10;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *fh = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  INT j;
  for (j = 0; j < M; j++)
    x[j] = K(0.0);

  CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(0, &N, 0, &n, M, 2, NFFT_WINDOW_KAISER_BESSEL,
                                     x, fh, f, 0u, NFFT_ESTIMATE)); /* d <= 0 */
  CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(1, 0, 0, &n, M, 2, NFFT_WINDOW_KAISER_BESSEL,
                                     x, fh, f, 0u, NFFT_ESTIMATE)); /* N == NULL */
  CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(1, &N, 0, 0, M, 2, NFFT_WINDOW_KAISER_BESSEL,
                                     x, fh, f, 0u, NFFT_ESTIMATE)); /* n == NULL */
  CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(1, &N, 0, &n, M, 2, NFFT_WINDOW_KAISER_BESSEL,
                                     0, fh, f, 0u, NFFT_ESTIMATE)); /* x == NULL */
  CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(1, &N, 0, &n, M, 2, NFFT_WINDOW_KAISER_BESSEL,
                                     x, 0, f, 0u, NFFT_ESTIMATE)); /* f_hat == NULL */
  CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(1, &N, 0, &n, M, 2, NFFT_WINDOW_KAISER_BESSEL,
                                     x, fh, 0, 0u, NFFT_ESTIMATE)); /* f == NULL */

  Y(free)
  (x);
  Y(free)
  (fh);
  Y(free)
  (f);
}

/* The guru must return NULL on non-positive per-axis geometry, in release
 * builds too. A unit axis (N[t] == 1) stays valid. */
void Y(check_nplan_guru_rejects_bad_geometry)(void) {
  INT Nzero[2] = {4, 0}, nok[2] = {8, 8}; /* N[1] == 0 */
  INT Nok[2] = {4, 4}, nzero[2] = {8, 0}; /* n[1] == 0 */
  INT M = 10;
  R *x = (R *)Y(malloc)((size_t)2 * (size_t)M * sizeof(R));
  C *fh = (C *)Y(malloc)((size_t)16 * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  INT j;
  for (j = 0; j < 2 * M; j++)
    x[j] = K(0.0);

  CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(2, Nzero, 0, nok, M, 2,
                                     NFFT_WINDOW_KAISER_BESSEL, x, fh, f, 0u, NFFT_ESTIMATE)); /* N[1] == 0 */
  CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(2, Nok, 0, nzero, M, 2,
                                     NFFT_WINDOW_KAISER_BESSEL, x, fh, f, 0u, NFFT_ESTIMATE)); /* n[1] == 0 */

  /* With the fast solver in play the guru refuses sigma = n/N <= 1 on any axis
   * surviving elision; NFFT_NO_FAST_NATIVE makes it legal again. */
  {
    INT Neq[2] = {4, 4}, neq[2] = {8, 4};   /* sigma == 1 on axis 1 */
    INT Nlt[2] = {4, 8}, nlt[2] = {8, 4};   /* sigma  < 1 on axis 1 */
    INT Nun[2] = {4, 1}, nun[2] = {8, 1};   /* unit axis: exempt, elided */
    Y(plan_ng) * pu;
    CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(2, Neq, 0, neq, M, 2,
                                       NFFT_WINDOW_KAISER_BESSEL, x, fh, f, 0u, NFFT_ESTIMATE));
    CU_ASSERT_PTR_NULL(Y(plan_ng_guru)(2, Nlt, 0, nlt, M, 2,
                                       NFFT_WINDOW_KAISER_BESSEL, x, fh, f, 0u, NFFT_ESTIMATE));
    /* Same geometry, fast path excluded by the caller. */
    pu = Y(plan_ng_guru)(2, Neq, 0, neq, M, 2, NFFT_WINDOW_KAISER_BESSEL, x, fh,
                         f, 0u, NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
    CU_ASSERT_PTR_NOT_NULL(pu);
    Y(plan_ng_destroy)
    (pu);
    /* Unit axes are elided, so n[t] == N[t] == 1 stays legal. */
    pu = Y(plan_ng_guru)(2, Nun, 0, nun, M, 2,
                         NFFT_WINDOW_KAISER_BESSEL, x, fh, f, 0u, NFFT_ESTIMATE);
    CU_ASSERT_PTR_NOT_NULL(pu);
    Y(plan_ng_destroy)
    (pu);
  }

  Y(free)
  (x);
  Y(free)
  (fh);
  Y(free)
  (f);
}
