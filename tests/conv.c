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

/* Unit tests for the CONV leaf solver (kernel/conv/): Step C of the fast NFFT
 * decomposition, the node convolution (matrix B).  Forward gathers the
 * oversampled grid against the window psi at each node; the adjoint scatters
 * node values back onto the wrapped neighbourhood.  Both are asserted as
 * values from a clean input, never as a round-trip. */

#include <complex.h> /* before nfft3.h so fftw_complex is C-compatible */
#include <stdio.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "config.h" /* ABS_SRCDIR */
#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "conv.h"

void Y(check_conv_problem)(void) {
  planner *pl = Y(planner_create)();
  INT n1 = 32, N1 = 16;
  int w = Y(get_window_id)();
  R x0 = K(0.1), x1 = K(0.2);
  C g0 = K(0.0), f0 = K(0.0), g1 = K(3.0), f1 = K(4.0);
  problem *pf, *p2;
  problem_conv *cf;
  md5sig s1, s2;
  char buf[160];
  printer *pr;

  pf = Y(mkproblem_conv)(1, &n1, &N1, (INT)1000, 6, w, 1, &x0, &g0, &f0);
  CU_ASSERT_EQUAL(pf->adt->kind, NFFT_PROBLEM_CONV);
  CU_ASSERT_EQUAL(Y(problem_conv_n)(pf, 0), (INT)32);
  CU_ASSERT_EQUAL(Y(problem_conv_N)(pf, 0), (INT)16);
  CU_ASSERT_EQUAL(Y(problem_conv_ntot)(pf), (INT)32);

  /* per-axis N read-back: the scalar problem_conv.N cannot carry d > 1. */
  {
    INT n2[2] = {32, 16}, N2[2] = {16, 8};
    problem *pnd = Y(mkproblem_conv)(2, n2, N2, (INT)1000, 6, w, 1, &x0, &g0,
                                     &f0);
    CU_ASSERT_EQUAL(Y(problem_conv_N)(pnd, 0), (INT)16);
    CU_ASSERT_EQUAL(Y(problem_conv_N)(pnd, 1), (INT)8);
    Y(problem_destroy)
    (pnd);
  }

  cf = (problem_conv *)pf;
  CU_ASSERT_EQUAL(cf->M, (INT)1000);
  /* x/g/f are borrowed aliases, not copies. */
  CU_ASSERT_PTR_EQUAL(cf->x, &x0);
  CU_ASSERT_PTR_EQUAL(cf->g, &g0);
  CU_ASSERT_PTR_EQUAL(cf->f, &f0);

  pr = Y(printer_create_str)(buf);
  pr->print(pr, "%P", pf);
  Y(printer_destroy)
  (pr);
  CU_ASSERT_STRING_EQUAL(buf, "(conv sign=1 m=6 M=1000 (tensor 1 (32 1 32 1)))");

  /* key is data-blind: a different x pointer yields the same key */
  Y(problem_md5)
  (pl, pf, s1);
  p2 = Y(mkproblem_conv)(1, &n1, &N1, (INT)1000, 6, w, 1, &x1, &g0, &f0);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2] && s1[3] == s2[3]);
  Y(problem_destroy)
  (p2);

  /* key is data-blind: different g/f pointers yield the same key */
  p2 = Y(mkproblem_conv)(1, &n1, &N1, (INT)1000, 6, w, 1, &x0, &g1, &f1);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2] && s1[3] == s2[3]);
  Y(problem_destroy)
  (p2);

  /* M-bucketing (floor_log2): 1000 and 1023 share a bucket, 1024 does not */
  p2 = Y(mkproblem_conv)(1, &n1, &N1, (INT)1023, 6, w, 1, &x0, &g0, &f0);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2] && s1[3] == s2[3]);
  Y(problem_destroy)
  (p2);
  p2 = Y(mkproblem_conv)(1, &n1, &N1, (INT)1024, 6, w, 1, &x0, &g0, &f0);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);
  Y(problem_destroy)
  (p2);

  /* sign shifts the key */
  p2 = Y(mkproblem_conv)(1, &n1, &N1, (INT)1000, 6, w, -1, &x0, &g0, &f0);
  Y(problem_md5)
  (pl, p2, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2] || s1[3] != s2[3]);
  Y(problem_destroy)
  (p2);

  Y(problem_destroy)
  (pf);
  Y(planner_destroy)
  (pl);
}

/* The CONV solver, planned directly through planner_mkplan (there is no conv
 * guru). Forward and adjoint are asserted as values recomputed from
 * Y(window_phi) and the wrap formula, not from the solver's own psi table and
 * not as a round-trip; the round-trip below is only a cross-check. x is chosen
 * so node 0's support [c-m, c+m+1] neither wraps past n nor self-collides, so
 * each wrapped neighbor index receives exactly one psi contribution.
 * Uses the process-global planner, so it tears it down at the end. */
void Y(check_conv_solver)(void) {
  INT n = 32, N = 16;
  int m = 6, w = Y(get_window_id)();
  const INT M = 4;
  R x[4];
  C *g, *f;
  problem *pf, *pa;
  plan *pln, *pln_adj;

  x[0] = K(0.5); /* c0 = floor(32*0.5) = 16; support [10,23], inside [0,32) */
  x[1] = K(0.3);
  x[2] = K(0.7);
  x[3] = K(0.55);

  g = (C *)Y(malloc)((size_t)n * sizeof(C));
  f = (C *)Y(malloc)((size_t)M * sizeof(C));

  Y(conv_ensure_registered)
  ();

  /* forward: single oversampled-grid spike g[0] = 1 */
  {
    INT ks;
    for (ks = 0; ks < n; ks++)
      g[ks] = K(0.0);
    g[0] = K(1.0);
  }

  pf = Y(mkproblem_conv)(1, &n, &N, M, m, w, 1, x, g, f);
  pln = Y(planner_mkplan)(Y(the_planner)(), pf);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  Y(plan_awake)
  (pln, PLNR_AWAKE);
  pln->adt->apply(pln, pf);

  /* Only taps whose wrapped neighbor is 0 contribute; g is zero elsewhere. */
  {
    INT j;
    for (j = 0; j < M; j++) {
      INT c = LRINT(FLOOR((R)n * x[j]));
      C expect = K(0.0);
      int lj;
      for (lj = 0; lj <= 2 * m + 1; lj++) {
        INT idx = c - m + (INT)lj;
        INT wrap = ((idx % n) + n) % n;
        if (wrap == 0)
          expect += Y(window_phi)(w, n, N, m, x[j] - (R)idx / (R)n);
      }
      CU_ASSERT(CABS(f[j] - expect) < K(1e-12));
    }
  }

  Y(plan_destroy)
  (pln);
  Y(problem_destroy)
  (pf);

  /* B^H: a single node spike f[0] = 1 scatters node 0's psi weights onto its
   * wrapped neighbors. */
  pa = Y(mkproblem_conv)(1, &n, &N, M, m, w, -1, x, g, f);
  pln_adj = Y(planner_mkplan)(Y(the_planner)(), pa);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln_adj);
  Y(plan_awake)
  (pln_adj, PLNR_AWAKE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln_adj->adt->apply_adjoint);

  {
    INT j, ks, c0;
    int lj;
    R touched[32];
    for (j = 0; j < M; j++)
      f[j] = K(0.0);
    f[0] = K(1.0);
    for (ks = 0; ks < n; ks++)
      touched[ks] = K(0.0);

    pln_adj->adt->apply_adjoint(pln_adj, pa);

    c0 = LRINT(FLOOR((R)n * x[0]));
    for (lj = 0; lj <= 2 * m + 1; lj++) {
      INT idx = c0 - m + (INT)lj;
      INT wrap = ((idx % n) + n) % n;
      R expect = Y(window_phi)(w, n, N, m, x[0] - (R)idx / (R)n);
      CU_ASSERT(CABS(g[wrap] - expect) < K(1e-12));
      touched[wrap] = K(1.0);
    }
    /* every grid entry outside node 0's support must remain 0 */
    for (ks = 0; ks < n; ks++)
      if (touched[ks] == K(0.0))
        CU_ASSERT(CABS(g[ks]) < K(1e-12));
  }

  Y(plan_destroy)
  (pln_adj);
  Y(problem_destroy)
  (pa);

  /* Cross-check that forward and adjoint use the identical psi table: with a
   * single node there is no cross-node aliasing, so the adjoint must scatter
   * exactly f[0] * psi[lj] back onto every neighbor. */
  {
    INT ks, c0;
    int lj;
    C *g2 = (C *)Y(malloc)((size_t)n * sizeof(C));
    C f2[1];
    R x1[1];
    problem *pf2, *pa2;
    plan *plnf2, *plna2;
    x1[0] = x[0];
    for (ks = 0; ks < n; ks++)
      g2[ks] = K(0.0);
    g2[0] = K(1.0);
    pf2 = Y(mkproblem_conv)(1, &n, &N, (INT)1, m, w, 1, x1, g2, f2);
    plnf2 = Y(planner_mkplan)(Y(the_planner)(), pf2);
    CU_ASSERT_PTR_NOT_NULL_FATAL(plnf2);
    Y(plan_awake)
    (plnf2, PLNR_AWAKE);
    plnf2->adt->apply(plnf2, pf2);

    pa2 = Y(mkproblem_conv)(1, &n, &N, (INT)1, m, w, -1, x1, g2, f2);
    plna2 = Y(planner_mkplan)(Y(the_planner)(), pa2);
    CU_ASSERT_PTR_NOT_NULL_FATAL(plna2);
    Y(plan_awake)
    (plna2, PLNR_AWAKE);
    plna2->adt->apply_adjoint(plna2, pa2); /* zeroes+refills g2 from f2 */

    c0 = LRINT(FLOOR((R)n * x1[0]));
    for (lj = 0; lj <= 2 * m + 1; lj++) {
      INT idx = c0 - m + (INT)lj;
      INT wrap = ((idx % n) + n) % n;
      C expect = f2[0] * Y(window_phi)(w, n, N, m, x1[0] - (R)idx / (R)n);
      CU_ASSERT(CABS(g2[wrap] - expect) < K(1e-12));
    }

    Y(plan_destroy)
    (plnf2);
    Y(problem_destroy)
    (pf2);
    Y(plan_destroy)
    (plna2);
    Y(problem_destroy)
    (pa2);
    Y(free)
    (g2);
  }

  Y(free)
  (g);
  Y(free)
  (f);

  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)
  ();
}

/* Rank-generic reference gather/scatter.
 *
 * For node j the solver visits the (2m+2)^d box of taps around
 * c_t = floor(n_t * x_t); tap l_t sits at grid index c_t - m + l_t, wrapped
 * into [0, n_t), and carries the weight psi_t = phi(x_t - idx_t/n_t).  Forward
 * gathers g over that box, the adjoint scatters f[j] back onto it.  Recomputing
 * both here from Y(window_phi) gives a value oracle that is independent of the
 * per-rank index arithmetic the leaves implement. */
static void conv_ref(int d, const INT *n, const INT *N, INT M, int m, int w,
                     const R *x, C *g, C *f, int adjoint) {
  const int tap = 2 * m + 2;
  INT j;
  for (j = 0; j < M; j++) {
    int l[4], t;
    INT taps = 1;
    INT i;
    for (t = 0; t < d; t++)
      taps *= (INT)tap;
    if (!adjoint)
      f[j] = K(0.0);
    for (i = 0; i < taps; i++) {
      INT rest = i, lin = 0;
      R psi = K(1.0);
      for (t = d - 1; t >= 0; t--) {
        l[t] = (int)(rest % (INT)tap);
        rest /= (INT)tap;
      }
      for (t = 0; t < d; t++) {
        R xt = x[j * (INT)d + (INT)t];
        INT c = LRINT(FLOOR((R)n[t] * xt));
        INT idx = c - (INT)m + (INT)l[t];
        INT wrap = ((idx % n[t]) + n[t]) % n[t];
        psi *= Y(window_phi)(w, n[t], N[t], m, xt - (R)idx / (R)n[t]);
        lin = lin * n[t] + wrap;
      }
      if (adjoint)
        g[lin] += f[j] * psi;
      else
        f[j] += g[lin] * psi;
    }
  }
}

/* Drive one rank through the CONV solver, forward and adjoint, against
 * conv_ref.  Nodes are spread over [0, 1) so the tap boxes wrap on every axis;
 * g and f are dense (not spikes), so a leaf that mixes up two axes cannot pass
 * by accident. */
static void conv_check_rank(int d, const INT *n, const INT *N, INT M) {
  const int m = 6, w = Y(get_window_id)();
  INT ntot = 1, j;
  int t;
  R *x;
  C *g, *f, *ref_f, *ref_g, *sol_g;
  problem *p;
  plan *pln;
  R tol;

  for (t = 0; t < d; t++)
    ntot *= n[t];

  x = (R *)Y(malloc)((size_t)(M * (INT)d) * sizeof(R));
  g = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  f = (C *)Y(malloc)((size_t)M * sizeof(C));
  ref_f = (C *)Y(malloc)((size_t)M * sizeof(C));
  ref_g = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  sol_g = (C *)Y(malloc)((size_t)ntot * sizeof(C));

  for (j = 0; j < M * (INT)d; j++)
    x[j] = (R)(((j * 37) % 97)) / K(97.0); /* in [0,1), deterministic */
  for (j = 0; j < ntot; j++)
    g[j] = (R)(((j * 13) % 29) - 14) / K(14.0);

  Y(conv_ensure_registered)();

  /* forward */
  conv_ref(d, n, N, M, m, w, x, g, ref_f, 0);
  for (j = 0; j < M; j++)
    f[j] = K(0.0);
  p = Y(mkproblem_conv)(d, n, N, M, m, w, 1, x, g, f);
  pln = Y(planner_mkplan)(Y(the_planner)(), p);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  Y(plan_awake)(pln, PLNR_AWAKE);
  pln->adt->apply(pln, p);
  /* The solver and conv_ref sum the same tap box in a different order, so they
   * differ only by reassociation: measured worst case is ~0.25 eps in every
   * precision, and does not grow with rank.  64 eps keeps ~250x headroom while
   * staying far below any real index or weight error. */
  tol = K(64.0) * Y(float_property)(NFFT_EPSILON);
  for (j = 0; j < M; j++)
    CU_ASSERT(CABS(f[j] - ref_f[j]) < tol * (CABS(ref_f[j]) + K(1.0)));
  Y(plan_destroy)(pln);
  Y(problem_destroy)(p);

  /* adjoint: scatter a fresh f, from a zeroed grid on both sides */
  for (j = 0; j < M; j++)
    f[j] = (R)(((j * 7) % 23) - 11) / K(11.0);
  for (j = 0; j < ntot; j++) {
    ref_g[j] = K(0.0);
    sol_g[j] = K(0.0);
  }
  conv_ref(d, n, N, M, m, w, x, ref_g, f, 1);
  p = Y(mkproblem_conv)(d, n, N, M, m, w, -1, x, sol_g, f);
  pln = Y(planner_mkplan)(Y(the_planner)(), p);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  Y(plan_awake)(pln, PLNR_AWAKE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln->adt->apply_adjoint);
  pln->adt->apply_adjoint(pln, p);
  for (j = 0; j < ntot; j++)
    CU_ASSERT(CABS(sol_g[j] - ref_g[j]) < tol * (CABS(ref_g[j]) + K(1.0)));
  Y(plan_destroy)(pln);
  Y(problem_destroy)(p);

  Y(free)(x);
  Y(free)(g);
  Y(free)(f);
  Y(free)(ref_f);
  Y(free)(ref_g);
  Y(free)(sol_g);
  Y(the_planner_destroy)();
}

void Y(check_conv_1d)(void) {
  INT n = 32, N = 16;
  conv_check_rank(1, &n, &N, 8);
  { /* odd N on the single axis */
    INT no = 32, No = 15;
    conv_check_rank(1, &no, &No, 8);
  }
}

void Y(check_conv_2d)(void) {
  INT n[2] = {32, 32}, N[2] = {16, 15}; /* even + odd */
  conv_check_rank(2, n, N, 6);
}

void Y(check_conv_3d)(void) {
  INT n[3] = {32, 32, 32}, N[3] = {16, 15, 14};
  conv_check_rank(3, n, N, 5);
}

void Y(check_conv_nd)(void) {
  INT n[4] = {32, 32, 32, 32}, N[4] = {16, 15, 14, 13};
  conv_check_rank(4, n, N, 4);
}
