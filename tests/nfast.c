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
#include <stdio.h>
#include <string.h>
#include <strings.h> /* strcasecmp */
#include <CUnit/CUnit.h>
#ifdef _OPENMP
# include <omp.h>
#endif

#include "config.h" /* ABS_SRCDIR */
#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "plan_ng_test.h" /* test-only plan_ng_test_native_parts hook */
#include "nfast.h"
#include "util.h"
#include "alloc_count.h"

/* Native-fast and legacy-fast are the same algorithm, so they agree to
 * round-off: ~1e-12 in double/long-double, but only ~1e-4 in float, where
 * eps ~1e-7 accumulates over the deconv+FFT+conv pipeline. */
#if defined(NFFT_SINGLE)
# define NFAST_LEGACY_REL_BOUND ((R)1e-4)
#else
# define NFAST_LEGACY_REL_BOUND ((R)1e-12)
#endif

/* Select which implemented windows the runtime accuracy loop runs.
 * NFAST_WINDOWS = comma list of kb|kaiserbessel,gaussian,bspline,sinc (case-
 * insensitive); unset or "all" => all four. Unknown token => abort loudly. */
static int win_token(const char *t)
{
  if (!strcasecmp(t, "kb") || !strcasecmp(t, "kaiserbessel"))
    return NFFT_WINDOW_KAISER_BESSEL;
  if (!strcasecmp(t, "gaussian"))
    return NFFT_WINDOW_GAUSSIAN;
  if (!strcasecmp(t, "bspline") || !strcasecmp(t, "b_spline"))
    return NFFT_WINDOW_B_SPLINE;
  if (!strcasecmp(t, "sinc") || !strcasecmp(t, "sinc_power"))
    return NFFT_WINDOW_SINC_POWER;
  return -1;
}
static int windows_from_env(int *out)
{
  const char *env = getenv("NFAST_WINDOWS");
  int n = 0, w;
  if (!env || !*env || !strcasecmp(env, "all")) {
    for (w = NFFT_WINDOW_KAISER_BESSEL; w <= NFFT_WINDOW_SINC_POWER; w++)
      out[n++] = w;
    return n;
  }
  {
    char buf[128], *save = NULL, *tok;
    strncpy(buf, env, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    for (tok = strtok_r(buf, ",", &save); tok;
         tok = strtok_r(NULL, ",", &save)) {
      int id = win_token(tok);
      CU_ASSERT_FATAL(id >= 0);
      CU_ASSERT_FATAL(n < 4);
      out[n++] = id;
    }
  }
  return n;
}

void Y(check_nfast_window_id)(void)
{
  int id = Y(get_window_id)();
  CU_ASSERT(id >= NFFT_WINDOW_KAISER_BESSEL && id <= NFFT_WINDOW_DIRAC_DELTA);
  /* The id must agree with the compile-time window #defines, which both build
   * systems set alike -- unlike WINDOW_NAME's display spelling. */
#if defined(DIRAC_DELTA)
  CU_ASSERT_EQUAL(id, NFFT_WINDOW_DIRAC_DELTA);
#elif defined(GAUSSIAN)
  CU_ASSERT_EQUAL(id, NFFT_WINDOW_GAUSSIAN);
#elif defined(B_SPLINE)
  CU_ASSERT_EQUAL(id, NFFT_WINDOW_B_SPLINE);
#elif defined(SINC_POWER)
  CU_ASSERT_EQUAL(id, NFFT_WINDOW_SINC_POWER);
#else /* Kaiser-Bessel */
  CU_ASSERT_EQUAL(id, NFFT_WINDOW_KAISER_BESSEL);
#endif
}

void Y(check_nfast_deconv_problem)(void)
{
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
  Y(printer_destroy)(pr);
  CU_ASSERT_STRING_EQUAL(buf, "(deconv sign=1 m=6 (tensor 1 (16 1 32 1)) "
                              "variant=0)");

  /* adjoint: tensor_adjoint orientation, distinct key from forward */
  pa = Y(mkproblem_deconv)(1, &N1, 0, &n1, 6, w, -1, &fh, &g);
  CU_ASSERT_EQUAL(Y(problem_deconv_N)(pa, 0), (INT)16); /* direction-aware */
  CU_ASSERT_EQUAL(Y(problem_deconv_n)(pa, 0), (INT)32);

  Y(problem_md5)(pl, pf, s1);
  Y(problem_md5)(pl, pa, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2]
            || s1[3] != s2[3]);

  /* geometry and m both shift the key */
  p2 = Y(mkproblem_deconv)(1, &N1, 0, &n1, 8, w, 1, &fh, &g);
  Y(problem_md5)(pl, p2, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2]
            || s1[3] != s2[3]);
  Y(problem_destroy)(p2);

  /* key is data-blind: a different f_hat/g pointer yields the same key */
  p2 = Y(mkproblem_deconv)(1, &N1, 0, &n1, 6, w, 1, &fh2, &g);
  Y(problem_md5)(pl, p2, s2);
  CU_ASSERT(s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2]
            && s1[3] == s2[3]);
  Y(problem_destroy)(p2);

  /* an even-N type-II axis shifts the key; odd N normalizes type-II -> type-I */
  {
    int v2 = NFFT_NDFT_TYPE_II;
    INT No = 15, no = 32;
    problem *pii = Y(mkproblem_deconv)(1, &N1, &v2, &n1, 6, w, 1, &fh, &g);
    problem *poii = Y(mkproblem_deconv)(1, &No, &v2, &no, 6, w, 1, &fh, &g);
    Y(problem_md5)(pl, pii, s2);
    CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2]
              || s1[3] != s2[3]);
    CU_ASSERT_EQUAL(Y(problem_deconv_variant)(pii, 0), NFFT_NDFT_TYPE_II);
    CU_ASSERT_EQUAL(Y(problem_deconv_variant)(poii, 0), NFFT_NDFT_TYPE_I);
    Y(problem_destroy)(pii);
    Y(problem_destroy)(poii);
  }

  Y(problem_destroy)(pf);
  Y(problem_destroy)(pa);
  Y(planner_destroy)(pl);
}

void Y(check_nfast_conv_problem)(void)
{
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

  pf = Y(mkproblem_conv)(1, &n1, &N1, (INT)1000, 6, w, 1, 0, &x0, &g0, &f0);
  CU_ASSERT_EQUAL(pf->adt->kind, NFFT_PROBLEM_CONV);
  CU_ASSERT_EQUAL(Y(problem_conv_n)(pf, 0), (INT)32);
  CU_ASSERT_EQUAL(Y(problem_conv_N)(pf, 0), (INT)16);
  CU_ASSERT_EQUAL(Y(problem_conv_ntot)(pf), (INT)32);

  /* per-axis N read-back: the scalar problem_conv.N cannot carry d > 1. */
  {
    INT n2[2] = {32, 16}, N2[2] = {16, 8};
    problem *pnd =
         Y(mkproblem_conv)(2, n2, N2, (INT)1000, 6, w, 1, 0, &x0, &g0, &f0);
    CU_ASSERT_EQUAL(Y(problem_conv_N)(pnd, 0), (INT)16);
    CU_ASSERT_EQUAL(Y(problem_conv_N)(pnd, 1), (INT)8);
    Y(problem_destroy)(pnd);
  }

  cf = (problem_conv *)pf;
  CU_ASSERT_EQUAL(cf->M, (INT)1000);
  /* x/g/f are borrowed aliases, not copies. */
  CU_ASSERT_PTR_EQUAL(cf->x, &x0);
  CU_ASSERT_PTR_EQUAL(cf->g, &g0);
  CU_ASSERT_PTR_EQUAL(cf->f, &f0);

  pr = Y(printer_create_str)(buf);
  pr->print(pr, "%P", pf);
  Y(printer_destroy)(pr);
  CU_ASSERT_STRING_EQUAL(buf, "(conv sign=1 m=6 M=1000 (tensor 1 (32 1 32 "
                              "1)))");

  /* key is data-blind: a different x pointer yields the same key */
  Y(problem_md5)(
       pl, pf, s1);
  p2 = Y(mkproblem_conv)(1, &n1, &N1, (INT)1000, 6, w, 1, 0, &x1, &g0, &f0);
  Y(problem_md5)(pl, p2, s2);
  CU_ASSERT(s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2]
            && s1[3] == s2[3]);
  Y(problem_destroy)(p2);

  /* key is data-blind: different g/f pointers yield the same key */
  p2 = Y(mkproblem_conv)(1, &n1, &N1, (INT)1000, 6, w, 1, 0, &x0, &g1, &f1);
  Y(problem_md5)(pl, p2, s2);
  CU_ASSERT(s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2]
            && s1[3] == s2[3]);
  Y(problem_destroy)(p2);

  /* M-bucketing (floor_log2): 1000 and 1023 share a bucket, 1024 does not */
  p2 = Y(mkproblem_conv)(1, &n1, &N1, (INT)1023, 6, w, 1, 0, &x0, &g0, &f0);
  Y(problem_md5)(pl, p2, s2);
  CU_ASSERT(s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2]
            && s1[3] == s2[3]);
  Y(problem_destroy)(p2);
  p2 = Y(mkproblem_conv)(1, &n1, &N1, (INT)1024, 6, w, 1, 0, &x0, &g0, &f0);
  Y(problem_md5)(pl, p2, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2]
            || s1[3] != s2[3]);
  Y(problem_destroy)(p2);

  /* sign shifts the key */
  p2 = Y(mkproblem_conv)(1, &n1, &N1, (INT)1000, 6, w, -1, 0, &x0, &g0, &f0);
  Y(problem_md5)(pl, p2, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2]
            || s1[3] != s2[3]);
  Y(problem_destroy)(p2);

  /* x_sorted shifts the key: it is the one property of x the key records */
  p2 = Y(mkproblem_conv)(1, &n1, &N1, (INT)1000, 6, w, 1, 1, &x0, &g0, &f0);
  CU_ASSERT_EQUAL(((problem_conv *)p2)->x_sorted, 1);
  Y(problem_md5)(pl, p2, s2);
  CU_ASSERT(s1[0] != s2[0] || s1[1] != s2[1] || s1[2] != s2[2]
            || s1[3] != s2[3]);
  Y(problem_destroy)(p2);

  Y(problem_destroy)(pf);
  Y(planner_destroy)(pl);
}

/* The KB vtable is peak-normalized, so these assertions check only positivity
 * and monotone decay, which a uniform positive scale preserves. */
void Y(check_nfast_window_vtable)(void)
{
  INT N = 16, n = 32;
  int m = 6, k;
  /* phi_hut is even and real for |k| <= N/2. */
  for (k = -(int)(N / 2); k < (int)(N / 2); k++) {
    R vh = Y(window_phi_hut)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, (INT)k);
    CU_ASSERT(vh > K(0.0));
  }
  /* phi is peaked at 0 and decays. */
  CU_ASSERT(Y(window_phi)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, K(0.0))
            > Y(window_phi)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, K(0.3)));
  /* Dirac and out-of-range ordinals return 0. */
  CU_ASSERT(Y(window_phi_hut)(NFFT_WINDOW_GAUSSIAN, n, N, m, 0) > K(0.0));
  CU_ASSERT_EQUAL(
       Y(window_phi_hut)(NFFT_WINDOW_DIRAC_DELTA, n, N, m, 0), K(0.0));
  CU_ASSERT_EQUAL(Y(window_phi_hut)(99, n, N, m, 0), K(0.0));
}

/* phi_hut(0) = 1 and the decay is monotone. The scale is uniform, so the
 * ratio between two values equals the ratio of the unscaled I0 values. */
void Y(check_nfast_window_normalized)(void)
{
  const int w = NFFT_WINDOW_KAISER_BESSEL;
  INT n = 32, N = 16;
  int m = 6, k;
  CU_ASSERT(FABS(Y(window_phi_hut)(w, n, N, m, (INT)0) - K(1.0)) <= K(1e-5));
  for (k = -(int)(N / 2); k < (int)(N / 2); k++)
    CU_ASSERT(Y(window_phi_hut)(w, n, N, m, (INT)k) > K(0.0)
              && Y(window_phi_hut)(w, n, N, m, (INT)k) <= K(1.0) + K(1e-5));
  /* uniform-scale invariance: phi_hut(k1)/phi_hut(k2) is scale-free, so it
   * equals I0(a1)/I0(a2) recomputed from the unscaled bessel directly. */
  {
    R b = KPI * (K(2.0) - (R)N / (R)n);
    R t1 = K(2.0) * KPI * K(1.0) / (R)n, t2 = K(2.0) * KPI * K(3.0) / (R)n;
    R a1 = (R)m * SQRT(b * b - t1 * t1), a2 = (R)m * SQRT(b * b - t2 * t2);
    R got = Y(window_phi_hut)(w, n, N, m, (INT)1)
            / Y(window_phi_hut)(w, n, N, m, (INT)3);
    /* I0(a1)/I0(a2) without ever forming I0: exp(a1-a2) times the ratio of
     * the exponentially scaled values. */
    R ref = EXP(a1 - a2) * Y(bessel_i0_exp_scaled)(a1) / Y(bessel_i0_exp_scaled)(a2);
    CU_ASSERT(FABS(got - ref) <= K(1e-4) * (K(1.0) + FABS(ref)));
  }
  CU_ASSERT(
       Y(window_phi)(w, n, N, m, K(0.0)) > Y(window_phi)(w, n, N, m, K(0.3)));
  /* The raw Gaussian also has phi_hut(0) == 1, exactly. */
  CU_ASSERT(FABS(Y(window_phi_hut)(NFFT_WINDOW_GAUSSIAN, n, N, m, 0) - K(1.0))
            <= K(1e-5));
}

/* Each window is checked against an independent recomputation of its infft.h
 * macro math: KB by a scale-free ratio because it is peak-normalized, the
 * other three by relative equality because they are raw. */
void Y(check_nfast_window_all)(void)
{
  const INT n = 32, N = 16;
  const int m = 6;
  const R sigma = (R)n / (R)N;
  const R tol = K(1e-4);
  INT k;

  /* KB: scale-free ratio equals the unscaled I0 ratio. */
  {
    R b = KPI * (K(2.0) - (R)N / (R)n);
    R t1 = K(2.0) * KPI * K(1.0) / (R)n, t2 = K(2.0) * KPI * K(3.0) / (R)n;
    R a1 = (R)m * SQRT(b * b - t1 * t1), a2 = (R)m * SQRT(b * b - t2 * t2);
    R got = Y(window_phi_hut)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, 1)
            / Y(window_phi_hut)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, 3);
    /* I0(a1)/I0(a2) without ever forming I0: exp(a1-a2) times the ratio of
     * the exponentially scaled values. */
    R ref = EXP(a1 - a2) * Y(bessel_i0_exp_scaled)(a1) / Y(bessel_i0_exp_scaled)(a2);
    CU_ASSERT(FABS(got - ref) <= tol * (K(1.0) + FABS(ref)));
  }
  {
    /* WINDOW_STENCIL_REACH is m + 1 for the floor(n x) centring. */
    R b = (K(2.0) * sigma) / (K(2.0) * sigma - K(1.0))
          * (((R)m + K(1.0)) / KPI);
    for (k = 1; k <= 4; k++) {
      R t = KPI * (R)k / (R)n, ref = EXP(-(t * t) * b);
      R got = Y(window_phi_hut)(NFFT_WINDOW_GAUSSIAN, n, N, m, k);
      CU_ASSERT(FABS(got - ref) <= tol * (K(1.0) + FABS(ref)));
    }
    {
      R x = K(0.02), ref = EXP(-POW(x * (R)n, K(2.0)) / b) / SQRT(KPI * b);
      R got = Y(window_phi)(NFFT_WINDOW_GAUSSIAN, n, N, m, x);
      CU_ASSERT(FABS(got - ref) <= tol * (K(1.0) + FABS(ref)));
    }
  }
  {
    R ref0 = K(1.0) / (R)n;
    CU_ASSERT(FABS(Y(window_phi_hut)(NFFT_WINDOW_B_SPLINE, n, N, m, 0) - ref0)
              <= tol * (K(1.0) + FABS(ref0)));
    for (k = 1; k <= 4; k++) {
      R a = (R)k * KPI / (R)n;
      R ref = POW(SIN(a) / a, K(2.0) * (R)m) / (R)n;
      R got = Y(window_phi_hut)(NFFT_WINDOW_B_SPLINE, n, N, m, k);
      CU_ASSERT(FABS(got - ref) <= tol * (K(1.0) + FABS(ref)));
    }
    {
      R x = K(0.02), ref = Y(bsplines)((INT)(2 * m), x * (R)n + (R)m) / (R)n;
      R got = Y(window_phi)(NFFT_WINDOW_B_SPLINE, n, N, m, x);
      CU_ASSERT(FABS(got - ref) <= tol * (K(1.0) + FABS(ref)));
    }
  }
  {
    for (k = 0; k <= 4; k++) {
      R arg =
           (K(2.0) * (R)m * (R)k) / ((K(2.0) * sigma - K(1.0)) * (R)n / sigma)
           + (R)m;
      R ref = Y(bsplines)((INT)(2 * m), arg);
      R got = Y(window_phi_hut)(NFFT_WINDOW_SINC_POWER, n, N, m, k);
      CU_ASSERT(FABS(got - ref) <= tol * (K(1.0) + FABS(ref)));
    }
    {
      R x = K(0.02);
      R ref = ((R)n / sigma) * (K(2.0) * sigma - K(1.0)) / (K(2.0) * (R)m)
              * POW(Y(sinc)(KPI * (R)n / sigma * x * (K(2.0) * sigma - K(1.0))
                         / (K(2.0) * (R)m)),
                    K(2.0) * (R)m)
              / (R)n;
      R got = Y(window_phi)(NFFT_WINDOW_SINC_POWER, n, N, m, x);
      CU_ASSERT(FABS(got - ref) <= tol * (K(1.0) + FABS(ref)));
    }
  }
  /* Dirac and out-of-range ordinals return 0. */
  CU_ASSERT_EQUAL(
       Y(window_phi_hut)(NFFT_WINDOW_DIRAC_DELTA, n, N, m, 0), K(0.0));
  CU_ASSERT_EQUAL(Y(window_phi)(NFFT_WINDOW_DIRAC_DELTA, n, N, m, K(0.0)),
                  K(0.0));
  CU_ASSERT_EQUAL(Y(window_phi_hut)(99, n, N, m, 0), K(0.0));
}

/* window is in the wisdom key: two NFFT problems differing only in window get
 * distinct keys. */
void Y(check_nfast_window_key)(void)
{
  planner *pl = Y(planner_create)();
  INT N = 16, n = 32, M = 50;
  R x = K(0.1);
  C fh = K(1.0), f = K(0.0);
  problem *pkb, *pg;
  md5sig skb, sgb;
  pkb =
       Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u,
                       &x, 0, &fh, &f);
  pg = Y(
       mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_GAUSSIAN, +1, 0u, &x, 0,
                      &fh, &f);
  Y(problem_md5)(pl, pkb, skb);
  Y(problem_md5)(pl, pg, sgb);
  CU_ASSERT(skb[0] != sgb[0] || skb[1] != sgb[1] || skb[2] != sgb[2]
            || skb[3] != sgb[3]);
  Y(problem_destroy)(pkb);
  Y(problem_destroy)(pg);
  Y(planner_destroy)(pl);
}

/* The DECONV solver, planned directly through planner_mkplan (there is no
 * deconv guru). Forward, adjoint and the type-II shift are each asserted as
 * values from a clean input, never as a round-trip: a round-trip also passes
 * for the inverse, which the adjoint of a real diagonal is not.
 * Uses the process-global planner, so it tears it down at the end. */
void Y(check_nfast_deconv_solver)(void)
{
  INT N = 16, n = 32;
  int m = 6, w = Y(get_window_id)();
  INT ks;
  C *f_hat, *g;
  problem *pf, *pa;
  plan *pln, *pln_adj;
  R phi_hut0;

  f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  g = (C *)Y(malloc)((size_t)n * sizeof(C));

  Y(deconv_ensure_registered)();

  /* forward: single spike at ks = N/2 (frequency k = 0) */
  for (ks = 0; ks < N; ks++)
    f_hat[ks] = K(0.0);
  f_hat[N / 2] = K(1.0);
  for (ks = 0; ks < n; ks++)
    g[ks] = K(0.0);

  pf = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, 1, f_hat, g);
  pln = Y(planner_mkplan)(Y(the_planner)(), pf);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  Y(plan_awake)(pln, PLNR_AWAKE);
  pln->adt->apply(pln, pf);

  phi_hut0 = Y(window_phi_hut)(w, n, N, m, (INT)0);
  CU_ASSERT(CABS(g[0] - K(1.0) / phi_hut0)
            <= NFAST_LEGACY_REL_BOUND * (K(1.0) / phi_hut0));
  for (ks = N; ks < n; ks++)
    CU_ASSERT(CABS(g[ks]) < K(1e-12));

  Y(plan_destroy)(pln);
  Y(problem_destroy)(pf);

  /* D is a real diagonal scale-and-pad, so D^H multiplies by the same
   * 1/phi_hut and gathers. */
  pa = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, -1, f_hat, g);
  pln_adj = Y(planner_mkplan)(Y(the_planner)(), pa);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln_adj);
  Y(plan_awake)(pln_adj, PLNR_AWAKE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln_adj->adt->apply_adjoint);

  /* frequency k=0 lives at pos 0: g[0]=1 -> f_hat[N/2] = 1/phi_hut(0). */
  for (ks = 0; ks < N; ks++)
    f_hat[ks] = K(0.0);
  for (ks = 0; ks < n; ks++)
    g[ks] = K(0.0);
  g[0] = K(1.0);
  pln_adj->adt->apply_adjoint(pln_adj, pa);
  CU_ASSERT(CABS(f_hat[N / 2] - K(1.0) / phi_hut0)
            <= NFAST_LEGACY_REL_BOUND * (K(1.0) / phi_hut0));

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
    CU_ASSERT(CABS(f_hat[N / 2 + 1] - K(1.0) / phi_hut1)
              <= NFAST_LEGACY_REL_BOUND * (K(1.0) / phi_hut1));
  }

  Y(plan_destroy)(pln_adj);
  Y(problem_destroy)(pa);

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
    Y(plan_awake)(pln2, PLNR_AWAKE);
    pln2->adt->apply(pln2, p2);
    /* Relative tolerance: 1/phi_hut(0) and 1/phi_hut(1) differ by only ~2.4%,
     * so the bound has to stay well inside that to tell them apart. */
    CU_ASSERT(CABS(g[1] - K(1.0) / phi_hut1)
              <= NFAST_LEGACY_REL_BOUND * (K(1.0) / phi_hut1));
    CU_ASSERT(CABS(g[0]) <= NFAST_LEGACY_REL_BOUND * (K(1.0) / phi_hut1));
    Y(plan_destroy)(pln2);
    Y(problem_destroy)(p2);
  }

  /* AWAKE_ZERO holds placeholder values. SLEEPY -> AWAKE_ZERO -> SLEEPY ->
   * AWAKE must match SLEEPY -> AWAKE. */
  {
    C *gz = (C *)Y(malloc)((size_t)n * sizeof(C));
    problem *pz, *pd;
    plan *plz, *pld;
    for (ks = 0; ks < N; ks++)
      f_hat[ks] = (R)(ks + 1);
    pz = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, 1, f_hat, gz);
    plz = Y(planner_mkplan)(Y(the_planner)(), pz);
    CU_ASSERT_PTR_NOT_NULL_FATAL(plz);
    Y(plan_awake)(plz, PLNR_AWAKE_ZERO);
    plz->adt->apply(plz, pz); /* placeholder tables: runnable, meaningless */
    Y(plan_awake)(plz, PLNR_SLEEPY);
    Y(plan_awake)(plz, PLNR_AWAKE);
    plz->adt->apply(plz, pz);

    pd = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, 1, f_hat, g);
    pld = Y(planner_mkplan)(Y(the_planner)(), pd);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pld);
    Y(plan_awake)(pld, PLNR_AWAKE);
    pld->adt->apply(pld, pd);

    for (ks = 0; ks < n; ks++)
      CU_ASSERT(gz[ks] == g[ks]);

    Y(plan_destroy)(plz);
    Y(problem_destroy)(pz);
    Y(plan_destroy)(pld);
    Y(problem_destroy)(pd);
    Y(free)(gz);
  }

  Y(free)(f_hat);
  Y(free)(g);

  /* Reset the process-global planner so later suites see a fresh generation. */
  Y(the_planner_destroy)();
}

/* Plan a DECONV problem through the process-global planner and wake it. */
static plan *deconv_awake_plan(problem *p)
{
  plan *pln = Y(planner_mkplan)(Y(the_planner)(), p);
  if (pln)
    Y(plan_awake)(pln, PLNR_AWAKE);
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
void Y(check_nfast_deconv_1d_general)(void)
{
  const int m = 6, w = Y(get_window_id)();
  Y(deconv_ensure_registered)();

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
      Y(plan_destroy)(pln);
      Y(problem_destroy)(p);
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
      CU_ASSERT_EQUAL(count_above(f_hat, N, deconv_general_tol() * sc),
                      (INT)1);
      Y(plan_destroy)(pln);
      Y(problem_destroy)(p);
    }
    Y(free)(f_hat);
    Y(free)(g);
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
      Y(plan_destroy)(pln);
      Y(problem_destroy)(p);
    }
    Y(free)(f_hat);
    Y(free)(g);
  }

  /* (c) n < N is declined, not planned: the zero-pad length would wrap. */
  {
    INT N = 32, n = 16;
    C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *g = (C *)Y(malloc)((size_t)n * sizeof(C));
    problem *p = Y(mkproblem_deconv)(1, &N, 0, &n, m, w, 1, f_hat, g);
    CU_ASSERT_PTR_NULL(Y(planner_mkplan)(Y(the_planner)(), p));
    Y(problem_destroy)(p);
    Y(free)(f_hat);
    Y(free)(g);
  }

  Y(the_planner_destroy)();
}

/* Rank-2 DECONV with a type-II axis and an odd axis. Slots are checked as
 * single spikes, so a one-cell index slip shows up as a miss plus a stray. */
void Y(check_nfast_deconv_2d_general)(void)
{
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

  Y(deconv_ensure_registered)();

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
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
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
    CU_ASSERT_EQUAL(count_above(f_hat, Ntot, deconv_general_tol() * sc),
                    (INT)1);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }

  Y(free)(f_hat);
  Y(free)(g);
  Y(the_planner_destroy)();
}

/* Rank-3 DECONV with two type-II axes and one odd axis. */
void Y(check_nfast_deconv_3d_general)(void)
{
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

  Y(deconv_ensure_registered)();

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
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
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
    CU_ASSERT_EQUAL(count_above(f_hat, Ntot, deconv_general_tol() * sc),
                    (INT)1);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }

  Y(free)(f_hat);
  Y(free)(g);
  Y(the_planner_destroy)();
}

/* Rank-4 DECONV with mixed type-II and odd axes, driving the generic carry
 * loop. */
void Y(check_nfast_deconv_nd_general)(void)
{
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

  Y(deconv_ensure_registered)();

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
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
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
    CU_ASSERT_EQUAL(count_above(f_hat, Ntot, deconv_general_tol() * sc),
                    (INT)1);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }

  Y(free)(f_hat);
  Y(free)(g);
  Y(the_planner_destroy)();
}

/* The CONV solver, planned directly through planner_mkplan (there is no conv
 * guru). Forward and adjoint are asserted as values recomputed from
 * Y(window_phi) and the wrap formula, not from the solver's own psi table and
 * not as a round-trip; the round-trip below is only a cross-check. x is chosen
 * so node 0's support [c-m, c+m+1] neither wraps past n nor self-collides, so
 * each wrapped neighbor index receives exactly one psi contribution.
 * Uses the process-global planner, so it tears it down at the end. */
void Y(check_nfast_conv_solver)(void)
{
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

  Y(conv_ensure_registered)();

  /* forward: single oversampled-grid spike g[0] = 1 */
  {
    INT ks;
    for (ks = 0; ks < n; ks++)
      g[ks] = K(0.0);
    g[0] = K(1.0);
  }

  pf = Y(mkproblem_conv)(1, &n, &N, M, m, w, 1, 0, x, g, f);
  pln = Y(planner_mkplan)(Y(the_planner)(), pf);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  Y(plan_awake)(pln, PLNR_AWAKE);
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
      CU_ASSERT(CABS(f[j] - expect)
                <= NFAST_LEGACY_REL_BOUND * FMAX(CABS(f[j]), CABS(expect)));
    }
  }

  Y(plan_destroy)(pln);
  Y(problem_destroy)(pf);

  /* B^H: a single node spike f[0] = 1 scatters node 0's psi weights onto its
   * wrapped neighbors. */
  pa = Y(mkproblem_conv)(1, &n, &N, M, m, w, -1, 0, x, g, f);
  pln_adj = Y(planner_mkplan)(Y(the_planner)(), pa);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln_adj);
  Y(plan_awake)(pln_adj, PLNR_AWAKE);
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
      /* The solver steps the run argument by one grid cell; this recomputes it
       * from x - idx/n, so the two agree to a rounding, not exactly. */
      R expect = Y(window_phi)(w, n, N, m, x[0] - (R)idx / (R)n);
      CU_ASSERT(CABS(g[wrap] - expect)
                <= NFAST_LEGACY_REL_BOUND * FMAX(CABS(g[wrap]), FABS(expect)));
      touched[wrap] = K(1.0);
    }
    /* every grid entry outside node 0's support must remain 0 */
    for (ks = 0; ks < n; ks++)
      if (touched[ks] == K(0.0))
        CU_ASSERT(CABS(g[ks]) < K(1e-12));
  }

  Y(plan_destroy)(pln_adj);
  Y(problem_destroy)(pa);

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
    pf2 = Y(mkproblem_conv)(1, &n, &N, (INT)1, m, w, 1, 0, x1, g2, f2);
    plnf2 = Y(planner_mkplan)(Y(the_planner)(), pf2);
    CU_ASSERT_PTR_NOT_NULL_FATAL(plnf2);
    Y(plan_awake)(plnf2, PLNR_AWAKE);
    plnf2->adt->apply(plnf2, pf2);

    pa2 = Y(mkproblem_conv)(1, &n, &N, (INT)1, m, w, -1, 0, x1, g2, f2);
    plna2 = Y(planner_mkplan)(Y(the_planner)(), pa2);
    CU_ASSERT_PTR_NOT_NULL_FATAL(plna2);
    Y(plan_awake)(plna2, PLNR_AWAKE);
    plna2->adt->apply_adjoint(plna2, pa2); /* zeroes+refills g2 from f2 */

    c0 = LRINT(FLOOR((R)n * x1[0]));
    for (lj = 0; lj <= 2 * m + 1; lj++) {
      INT idx = c0 - m + (INT)lj;
      INT wrap = ((idx % n) + n) % n;
      C expect = f2[0] * Y(window_phi)(w, n, N, m, x1[0] - (R)idx / (R)n);
      CU_ASSERT(CABS(g2[wrap] - expect)
                <= NFAST_LEGACY_REL_BOUND
                        * FMAX(CABS(g2[wrap]), CABS(expect)));
    }

    Y(plan_destroy)(plnf2);
    Y(problem_destroy)(pf2);
    Y(plan_destroy)(plna2);
    Y(problem_destroy)(pa2);
    Y(free)(g2);
  }

  /* AWAKE_ZERO holds placeholder values. SLEEPY -> AWAKE_ZERO -> SLEEPY ->
   * AWAKE must match SLEEPY -> AWAKE. */
  {
    INT ks;
    C *fz = (C *)Y(malloc)((size_t)M * sizeof(C));
    problem *pz, *pd;
    plan *plz, *pld;
    for (ks = 0; ks < n; ks++)
      g[ks] = (R)(ks + 1);
    pz = Y(mkproblem_conv)(1, &n, &N, M, m, w, 1, 0, x, g, fz);
    plz = Y(planner_mkplan)(Y(the_planner)(), pz);
    CU_ASSERT_PTR_NOT_NULL_FATAL(plz);
    Y(plan_awake)(plz, PLNR_AWAKE_ZERO);
    plz->adt->apply(plz, pz); /* placeholder tables: runnable, meaningless */
    Y(plan_awake)(plz, PLNR_SLEEPY);
    Y(plan_awake)(plz, PLNR_AWAKE);
    plz->adt->apply(plz, pz);

    pd = Y(mkproblem_conv)(1, &n, &N, M, m, w, 1, 0, x, g, f);
    pld = Y(planner_mkplan)(Y(the_planner)(), pd);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pld);
    Y(plan_awake)(pld, PLNR_AWAKE);
    pld->adt->apply(pld, pd);

    for (ks = 0; ks < M; ks++)
      CU_ASSERT(fz[ks] == f[ks]);

    Y(plan_destroy)(plz);
    Y(problem_destroy)(pz);
    Y(plan_destroy)(pld);
    Y(problem_destroy)(pd);
    Y(free)(fz);
  }

  Y(free)(g);
  Y(free)(f);

  Y(the_planner_destroy)();
}

/* d == 1 reference case; N is a scalar. */
static int read_1d_case(const char *rel, INT *N, INT *M, R **x, C **f_hat,
                        C **f)
{
  int d;
  INT *Nv, NN;
  if (!Y(test_read_case)(rel, &d, &Nv, &NN, M, x, f_hat, f))
    return 0;
  *N = Nv[0];
  Y(free)(Nv);
  return d == 1;
}

static int read_nd_case(const char *rel, int *d, INT **N, INT *M, R **x,
                        C **f_hat, C **f)
{
  INT NN;
  return Y(test_read_case)(rel, d, N, &NN, M, x, f_hat, f);
}

/* The composed native fast NFFT on the 1D even-N reference case, planned
 * under NFFT_NO_DIRECT so it is the sole survivor. The fast transform is an
 * approximation, so both directions are checked to 1e-5, not the 1e-9 direct
 * tolerance. The adjoint is compared against its own reference file: a
 * round-trip could not distinguish a correct adjoint from one off by a
 * constant that cancels. */
void Y(check_nfast_native_fast_accuracy)(void)
{
  INT N, M, Na, Ma;
  int m = 7;
  R *x, *xa;
  C *f_hat, *f, *f_hat_ref, *fa;
  INT n;
  Y(plan_ng) *p;
  R err;

  CU_ASSERT_TRUE_FATAL(
       read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;

  /* forward */
  {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    p = Y(plan_ng_guru)(1, &N, 0, &n, M, m,
                        Y(get_window_id)(), x, f_hat, got, 0u,
                     NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    Y(precompute)(p);
    Y(execute)(p);
    err = Y(test_rel_max_err)(got, f, M);
    CU_ASSERT(err < (R)1e-5);

    Y(test_assert_plan_names)(p, "nfft_solver_fast_native");
    Y(test_assert_plan_names)(p, "deconv");
    Y(test_assert_plan_names)(p, "conv");

    Y(plan_ng_destroy)(p);
    Y(free)(got);
  }

  /* adjoint: a separately generated reference file, not a round-trip. */
  CU_ASSERT_TRUE_FATAL(read_1d_case("data/nfft_adjoint_1d_20_50.txt", &Na, &Ma,
                                    &xa, &f_hat_ref, &fa));
  CU_ASSERT_EQUAL(Na, N);
  CU_ASSERT_EQUAL(Ma, M);
  {
    C *got_fhat = (C *)Y(malloc)((size_t)Na * sizeof(C));
    Y(plan_ng) *pa;
    pa = Y(plan_ng_guru)(1, &Na, 0, &n, Ma, m,
                         Y(get_window_id)(), xa, got_fhat, fa,
                      0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pa);
    Y(precompute)(pa);
    Y(execute_adjoint)(pa);
    err = Y(test_rel_max_err)(got_fhat, f_hat_ref, Na);
    CU_ASSERT(err < (R)1e-5);
    Y(plan_ng_destroy)(pa);
    Y(free)(got_fhat);
  }

  Y(free)(x);
  Y(free)(f_hat);
  Y(free)(f);
  Y(free)(xa);
  Y(free)(f_hat_ref);
  Y(free)(fa);

  Y(the_planner_destroy)();
}

/* Isolation tests for NFFT_NO_FAST_NATIVE, which gates the composed native
 * fast. Tests that need a specific winner pin it with a disable flag; the
 * rest stay winner-agnostic.
 *
 * With the direct natives excluded, the composed native fast is the sole
 * surviving candidate for a 1D even-N problem. */
void Y(check_nfast_native_tree)(void)
{
  INT N, M;
  int m = 7;
  R *x;
  C *f_hat, *f;
  INT n;
  Y(plan_ng) *p;

  CU_ASSERT_TRUE_FATAL(
       read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;

  {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    p = Y(plan_ng_guru)(1, &N, 0, &n, M, m,
                        Y(get_window_id)(), x, f_hat, got, 0u,
                     NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    Y(test_assert_plan_names)(p, "nfft_solver_fast_native");
    Y(test_assert_plan_names)(p, "deconv");
    Y(test_assert_plan_names)(p, "conv");

    Y(plan_ng_destroy)(p);
    Y(free)(got);
  }

  Y(free)(x);
  Y(free)(f_hat);
  Y(free)(f);
  Y(the_planner_destroy)();
}

/* A window other than the compile-time one makes every NFFT-kind solver
 * decline -- the native fast on its own window gate, the direct natives via
 * NFFT_NO_DIRECT -- so no candidate survives and the guru returns NULL. */
void Y(check_nfast_native_declines_window)(void)
{
  INT N, M;
  int m = 7;
  R *x;
  C *f_hat, *f;
  INT n;
  Y(plan_ng) *p;
  int other;

  CU_ASSERT_TRUE_FATAL(
       read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;

  other = (Y(get_window_id)() == NFFT_WINDOW_KAISER_BESSEL) ?
               NFFT_WINDOW_GAUSSIAN :
               NFFT_WINDOW_KAISER_BESSEL;
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, m, other, x, f_hat, f, 0u,
                   NFFT_ESTIMATE | NFFT_NO_DIRECT | NFFT_NO_FAST_NATIVE);
  CU_ASSERT_PTR_NULL(p);

  Y(free)(x);
  Y(free)(f_hat);
  Y(free)(f);
  Y(the_planner_destroy)();
}

/* NFFT_NO_FAST_NATIVE is selective: with the direct natives still available
 * one of them wins; excluding those too leaves nothing. */
void Y(check_nfast_flag_selective)(void)
{
  INT N, M;
  int m = 7;
  R *x;
  C *f_hat, *f;
  INT n;
  Y(plan_ng) *p;

  CU_ASSERT_TRUE_FATAL(
       read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;

  /* the fast excluded, direct still available */
  {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    p = Y(plan_ng_guru)(1, &N, 0, &n, M, m,
                        Y(get_window_id)(), x, f_hat, got, 0u,
                     NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);
    Y(plan_ng_destroy)(p);
    Y(free)(got);
  }

  /* everything excluded */
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, m, Y(get_window_id)(), x, f_hat, f, 0u,
                   NFFT_ESTIMATE | NFFT_NO_DIRECT | NFFT_NO_FAST_NATIVE);
  CU_ASSERT_PTR_NULL(p);

  Y(free)(x);
  Y(free)(f_hat);
  Y(free)(f);
  Y(the_planner_destroy)();
}

/* The 2D reference cases, the same geometries tests/nfft.c's 2D suite runs,
 * plus the _t210 type-II variants. */
typedef struct {
  const char *file;
  const int *variant; /* NULL = all type-I */
} nfast_2d_case;

static const int nfast_v_ii_i[2] = {NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I};

static const nfast_2d_case files_2d[] = {
     {"data/nfft_2d_10_10_20.txt", 0},
     {"data/nfft_2d_10_10_50.txt", 0},
     {"data/nfft_2d_10_20_20.txt", 0},
     {"data/nfft_2d_10_20_50.txt", 0},
     {"data/nfft_2d_20_10_20.txt", 0},
     {"data/nfft_2d_20_10_50.txt", 0},
     {"data/nfft_2d_20_20_20.txt", 0},
     {"data/nfft_2d_20_20_50.txt", 0},
     {"data/nfft_2d_10_20_50_t210.txt", nfast_v_ii_i},
     {"data/nfft_2d_20_10_50_t210.txt", nfast_v_ii_i},
};
static const nfast_2d_case adjoint_files_2d[] = {
     {"data/nfft_adjoint_2d_10_10_20.txt", 0},
     {"data/nfft_adjoint_2d_10_10_50.txt", 0},
     {"data/nfft_adjoint_2d_10_20_20.txt", 0},
     {"data/nfft_adjoint_2d_10_20_50.txt", 0},
     {"data/nfft_adjoint_2d_20_10_20.txt", 0},
     {"data/nfft_adjoint_2d_20_10_50.txt", 0},
     {"data/nfft_adjoint_2d_20_20_20.txt", 0},
     {"data/nfft_adjoint_2d_20_20_50.txt", 0},
     {"data/nfft_adjoint_2d_10_20_50_t210.txt", nfast_v_ii_i},
     {"data/nfft_adjoint_2d_20_10_50_t210.txt", nfast_v_ii_i},
};
#define NFAST_2D_NFILES (sizeof(files_2d) / sizeof(files_2d[0]))

/* The d=2 forward slice: each reference case planned under NFFT_NO_DIRECT so
 * the composed native fast is the sole survivor, checked against the file to
 * 1e-5 and against an in-test legacy X(trafo_2d) to NFAST_LEGACY_REL_BOUND. */
void Y(check_nfast_native_fast_2d)(void)
{
  size_t fi;
  for (fi = 0; fi < NFAST_2D_NFILES; fi++) {
    int d, m = 7;
    INT *N, M;
    R *x;
    C *f_hat, *f;
    INT n[2];
    Y(plan_ng) *p;
    R err;

    CU_ASSERT_TRUE_FATAL(
         read_nd_case(files_2d[fi].file, &d, &N, &M, &x, &f_hat, &f));
    CU_ASSERT_EQUAL_FATAL(d, 2);
    n[0] = 2 * N[0];
    n[1] = 2 * N[1];

    {
      C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
      p = Y(plan_ng_guru)(2, N, files_2d[fi].variant, n, M, m,
                          Y(get_window_id)(), x,
                       f_hat, got, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL_FATAL(p);

      Y(precompute)(p);
      Y(execute)(p);
      err = Y(test_rel_max_err)(got, f, M);
      CU_ASSERT(err < (R)1e-5);

      Y(test_assert_plan_names)(p, "nfft_solver_fast_native");
      Y(test_assert_plan_names)(p, "deconv");
      Y(test_assert_plan_names)(p, "conv");

      /* in-test legacy reference: the same fast algorithm through the old
       * X(plan) API, which has no type-II concept -- type-I cases only. */
      if (!files_2d[fi].variant) {
        NFFT(plan) lp;
        int Ni[2], ni[2];
        INT j;
        Ni[0] = (int)N[0];
        Ni[1] = (int)N[1];
        ni[0] = (int)n[0];
        ni[1] = (int)n[1];
        NFFT(init_guru)(&lp, 2, Ni, (int)M, ni, m,
                  PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F
                       | FFTW_INIT | FFT_OUT_OF_PLACE,
                  FFTW_ESTIMATE);
        for (j = 0; j < M * 2; j++)
          lp.x[j] = x[j];
        for (j = 0; j < N[0] * N[1]; j++)
          lp.f_hat[j] = f_hat[j];
        NFFT(precompute_one_psi)(&lp);
        NFFT(trafo)(&lp);
        {
          R errl = Y(test_rel_max_err)(got, lp.f, M);
          CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
        }
        NFFT(finalize)(&lp);
      }

      Y(plan_ng_destroy)(p);
      Y(free)(got);
    }

    Y(free)(N);
    Y(free)(x);
    Y(free)(f_hat);
    Y(free)(f);
    Y(the_planner_destroy)();
  }
}

/* The d=2 adjoint slice: same geometries, reading the separately generated
 * nfft_adjoint_2d_*.txt references rather than round-tripping the forward
 * case, and checked against an in-test legacy X(adjoint_2d) as well. */
void Y(check_nfast_native_fast_2d_adjoint)(void)
{
  size_t fi;
  for (fi = 0; fi < NFAST_2D_NFILES; fi++) {
    int d, m = 7;
    INT *N, M;
    R *x;
    C *f_hat_ref, *f;
    INT n[2];
    Y(plan_ng) *p;
    R err;

    CU_ASSERT_TRUE_FATAL(read_nd_case(adjoint_files_2d[fi].file, &d, &N, &M,
                                      &x, &f_hat_ref, &f));
    CU_ASSERT_EQUAL_FATAL(d, 2);
    n[0] = 2 * N[0];
    n[1] = 2 * N[1];

    {
      C *got_fhat = (C *)Y(malloc)((size_t)(N[0] * N[1]) * sizeof(C));
      p = Y(plan_ng_guru)(2, N, adjoint_files_2d[fi].variant, n, M, m,
                          Y(get_window_id)(), x, got_fhat, f, 0u,
                       NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL_FATAL(p);

      Y(precompute)(p);
      Y(execute_adjoint)(p);
      err = Y(test_rel_max_err)(got_fhat, f_hat_ref, N[0] * N[1]);
      CU_ASSERT(err < (R)1e-5);

      if (!adjoint_files_2d[fi].variant) {
        NFFT(plan) lp;
        int Ni[2], ni[2];
        INT j;
        Ni[0] = (int)N[0];
        Ni[1] = (int)N[1];
        ni[0] = (int)n[0];
        ni[1] = (int)n[1];
        NFFT(init_guru)(&lp, 2, Ni, (int)M, ni, m,
                  PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F
                       | FFTW_INIT | FFT_OUT_OF_PLACE,
                  FFTW_ESTIMATE);
        for (j = 0; j < M * 2; j++)
          lp.x[j] = x[j];
        for (j = 0; j < M; j++)
          lp.f[j] = f[j];
        NFFT(precompute_one_psi)(&lp);
        NFFT(adjoint)(&lp);
        {
          R errl = Y(test_rel_max_err)(got_fhat, lp.f_hat, N[0] * N[1]);
          CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
        }
        NFFT(finalize)(&lp);
      }

      Y(plan_ng_destroy)(p);
      Y(free)(got_fhat);
    }

    Y(free)(N);
    Y(free)(x);
    Y(free)(f_hat_ref);
    Y(free)(f);
    Y(the_planner_destroy)();
  }
}

/* The d=3 forward slice, single reference geometry N={10,10,10}, M=10:
 * planned under NFFT_NO_DIRECT so the composed native fast is the sole
 * survivor, checked against the file to 1e-5 and against an in-test legacy
 * X(trafo_3d) to NFAST_LEGACY_REL_BOUND. */
void Y(check_nfast_native_fast_3d)(void)
{
  int d, m = 7;
  INT *N, M;
  R *x;
  C *f_hat, *f;
  INT n[3];
  Y(plan_ng) *p;
  R err;

  CU_ASSERT_TRUE_FATAL(read_nd_case("data/nfft_3d_10_10_10_10.txt", &d, &N, &M,
                                    &x, &f_hat, &f));
  CU_ASSERT_EQUAL_FATAL(d, 3);
  n[0] = 2 * N[0];
  n[1] = 2 * N[1];
  n[2] = 2 * N[2];

  {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    p = Y(plan_ng_guru)(3, N, 0, n, M, m,
                        Y(get_window_id)(), x, f_hat, got, 0u,
                     NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    Y(precompute)(p);
    Y(execute)(p);
    err = Y(test_rel_max_err)(got, f, M);
    /* The double-generated 3D reference is float-ill-conditioned: heavy
     * cancellation in the output, so no fast NFFT (legacy included) reproduces
     * it to 1e-5 in float. Float correctness rests on the native-vs-legacy
     * agreement below and on check_nfast_float_accuracy. */
#ifndef NFFT_SINGLE
    CU_ASSERT(err < (R)1e-5);
#else
    (void)err;
#endif

    Y(test_assert_plan_names)(p, "nfft_solver_fast_native");
    Y(test_assert_plan_names)(p, "deconv");
    Y(test_assert_plan_names)(p, "conv");

    /* in-test legacy reference: the same fast algorithm through the old
     * X(plan) API. */
    {
      NFFT(plan) lp;
      int Ni[3], ni[3];
      INT j;
      Ni[0] = (int)N[0];
      Ni[1] = (int)N[1];
      Ni[2] = (int)N[2];
      ni[0] = (int)n[0];
      ni[1] = (int)n[1];
      ni[2] = (int)n[2];
      NFFT(init_guru)(&lp, 3, Ni, (int)M, ni, m,
                PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F
                     | FFTW_INIT | FFT_OUT_OF_PLACE,
                FFTW_ESTIMATE);
      for (j = 0; j < M * 3; j++)
        lp.x[j] = x[j];
      for (j = 0; j < N[0] * N[1] * N[2]; j++)
        lp.f_hat[j] = f_hat[j];
      NFFT(precompute_one_psi)(&lp);
      NFFT(trafo)(&lp);
      {
        R errl = Y(test_rel_max_err)(got, lp.f, M);
#ifndef NFFT_SINGLE
        CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
#else
        /* legacy collapses in float 3D and diverges from the native */
        (void)errl;
#endif
      }
      NFFT(finalize)(&lp);
    }

    Y(plan_ng_destroy)(p);
    Y(free)(got);
  }

  Y(free)(N);
  Y(free)(x);
  Y(free)(f_hat);
  Y(free)(f);
  Y(the_planner_destroy)();
}

/* The d=3 adjoint slice: same geometry, reading the separately generated
 * nfft_adjoint_3d_10_10_10_10.txt reference rather than round-tripping the
 * forward case, and checked against an in-test legacy X(adjoint_3d) as
 * well. */
void Y(check_nfast_native_fast_3d_adjoint)(void)
{
  int d, m = 7;
  INT *N, M;
  R *x;
  C *f_hat_ref, *f;
  INT n[3];
  Y(plan_ng) *p;
  R err;

  CU_ASSERT_TRUE_FATAL(read_nd_case("data/nfft_adjoint_3d_10_10_10_10.txt", &d,
                                    &N, &M, &x, &f_hat_ref, &f));
  CU_ASSERT_EQUAL_FATAL(d, 3);
  n[0] = 2 * N[0];
  n[1] = 2 * N[1];
  n[2] = 2 * N[2];

  {
    C *got_fhat = (C *)Y(malloc)((size_t)(N[0] * N[1] * N[2]) * sizeof(C));
    p = Y(plan_ng_guru)(3, N, 0, n, M, m,
                        Y(get_window_id)(), x, got_fhat, f, 0u,
                     NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    Y(precompute)(p);
    Y(execute_adjoint)(p);
    err = Y(test_rel_max_err)(got_fhat, f_hat_ref, N[0] * N[1] * N[2]);
    /* Same float-conditioning caveat as the 3D forward case. */
#ifndef NFFT_SINGLE
    CU_ASSERT(err < (R)1e-5);
#else
    (void)err;
#endif

    {
      NFFT(plan) lp;
      int Ni[3], ni[3];
      INT j;
      Ni[0] = (int)N[0];
      Ni[1] = (int)N[1];
      Ni[2] = (int)N[2];
      ni[0] = (int)n[0];
      ni[1] = (int)n[1];
      ni[2] = (int)n[2];
      NFFT(init_guru)(&lp, 3, Ni, (int)M, ni, m,
                PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F
                     | FFTW_INIT | FFT_OUT_OF_PLACE,
                FFTW_ESTIMATE);
      for (j = 0; j < M * 3; j++)
        lp.x[j] = x[j];
      for (j = 0; j < M; j++)
        lp.f[j] = f[j];
      NFFT(precompute_one_psi)(&lp);
      NFFT(adjoint)(&lp);
      {
        R errl = Y(test_rel_max_err)(got_fhat, lp.f_hat, N[0] * N[1] * N[2]);
#ifndef NFFT_SINGLE
        CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
#else
        /* legacy collapses in float 3D and diverges from the native */
        (void)errl;
#endif
      }
      NFFT(finalize)(&lp);
    }

    Y(plan_ng_destroy)(p);
    Y(free)(got_fhat);
  }

  Y(free)(N);
  Y(free)(x);
  Y(free)(f_hat_ref);
  Y(free)(f);
  Y(the_planner_destroy)();
}

/* The generic rnk >= 4 slice, served by the deconv-nd/conv-nd leaves.
 * tests/refgen generates no d >= 4 reference file, so the case is built
 * in-test from a deterministic sequence. Two oracles: an in-test legacy
 * X(plan) to NFAST_LEGACY_REL_BOUND, and the direct NDFT native, forced via
 * NFFT_NO_FAST_NATIVE, to the coarser 1e-5 fast-NFFT bound. */
static R seq_4d(INT k)
{
  /* Deterministic stand-in for a PRNG, so results do not depend on the
   * platform's rand. */
  R v = SIN((R)(k + 1) * K(12.9898)) * K(43758.5453);
  v -= FLOOR(v); /* fractional part, in [0,1) */
  return v;
}

void Y(check_nfast_native_fast_4d)(void)
{
  const int d = 4;
  INT N[4] = {10, 10, 10, 10};
  INT n[4];
  const INT M = 40;
  int m = 7, t;
  INT Ntot, j;
  R *x;
  C *f_hat, *got;
  Y(plan_ng) *p;
  R err;

  for (t = 0; t < d; t++)
    n[t] = 2 * N[t];
  for (t = 0, Ntot = 1; t < d; t++)
    Ntot *= N[t];

  x = (R *)Y(malloc)((size_t)M * (size_t)d * sizeof(R));
  f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  got = (C *)Y(malloc)((size_t)M * sizeof(C));

  for (j = 0; j < M * d; j++)
    x[j] = seq_4d(j) - K(0.5); /* in [-0.5, 0.5) */
  for (j = 0; j < Ntot; j++)
    f_hat[j] = (seq_4d(1000 + 2 * j) - K(0.5))
               + II * (seq_4d(1000 + 2 * j + 1) - K(0.5));

  p = Y(plan_ng_guru)(d, N, 0, n, M, m, Y(get_window_id)(), x, f_hat, got, 0u,
                   NFFT_ESTIMATE | NFFT_NO_DIRECT);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);

  Y(precompute)(p);
  Y(execute)(p);

  Y(test_assert_plan_names)(p, "nfft_solver_fast_native");
  Y(test_assert_plan_names)(p, "deconv");
  Y(test_assert_plan_names)(p, "conv");

  /* in-test legacy reference: the same fast algorithm through the old
   * X(plan) API. */
  {
    NFFT(plan) lp;
    int Ni[4], ni[4];
    INT jj;
    for (t = 0; t < d; t++) {
      Ni[t] = (int)N[t];
      ni[t] = (int)n[t];
    }
    NFFT(init_guru)(&lp, d, Ni, (int)M, ni, m,
              PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F
                   | FFTW_INIT | FFT_OUT_OF_PLACE,
              FFTW_ESTIMATE);
    for (jj = 0; jj < M * d; jj++)
      lp.x[jj] = x[jj];
    for (jj = 0; jj < Ntot; jj++)
      lp.f_hat[jj] = f_hat[jj];
    NFFT(precompute_one_psi)(&lp);
    NFFT(trafo)(&lp);
    {
      R errl = Y(test_rel_max_err)(got, lp.f, M);
      CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
    }
    NFFT(finalize)(&lp);
  }

  /* second oracle: the direct NDFT native, forced via NFFT_NO_FAST_NATIVE. */
  {
    C *got_direct = (C *)Y(malloc)((size_t)M * sizeof(C));
    Y(plan_ng) *pdirect =
         Y(plan_ng_guru)(d, N, 0, n, M, m, NFFT_WINDOW_KAISER_BESSEL, x, f_hat,
                      got_direct, 0u, NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pdirect);
    Y(precompute)(pdirect);
    Y(execute)(pdirect);
    err = Y(test_rel_max_err)(got, got_direct, M);
    CU_ASSERT(err < (R)1e-5);
    Y(plan_ng_destroy)(pdirect);
    Y(free)(got_direct);
  }

  Y(plan_ng_destroy)(p);
  Y(free)(x);
  Y(free)(f_hat);
  Y(free)(got);
  Y(the_planner_destroy)();
}

/* The d=4 adjoint slice: same in-test-generated geometry, oracle is
 * NFFT(adjoint) on a separate legacy plan. f is generated independently, so
 * this is not a round-trip of the forward half. */
void Y(check_nfast_native_fast_4d_adjoint)(void)
{
  const int d = 4;
  INT N[4] = {10, 10, 10, 10};
  INT n[4];
  const INT M = 40;
  int m = 7, t;
  INT Ntot, j;
  R *x;
  C *f, *got_fhat;
  Y(plan_ng) *p;

  for (t = 0; t < d; t++)
    n[t] = 2 * N[t];
  for (t = 0, Ntot = 1; t < d; t++)
    Ntot *= N[t];

  x = (R *)Y(malloc)((size_t)M * (size_t)d * sizeof(R));
  f = (C *)Y(malloc)((size_t)M * sizeof(C));
  got_fhat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));

  for (j = 0; j < M * d; j++)
    x[j] = seq_4d(j + 5000) - K(0.5); /* different offset than forward */
  for (j = 0; j < M; j++)
    f[j] = (seq_4d(9000 + 2 * j) - K(0.5))
           + II * (seq_4d(9000 + 2 * j + 1) - K(0.5));

  p = Y(plan_ng_guru)(d, N, 0, n, M, m, Y(get_window_id)(), x, got_fhat, f, 0u,
                   NFFT_ESTIMATE | NFFT_NO_DIRECT);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);

  Y(precompute)(p);
  Y(execute_adjoint)(p);

  {
    NFFT(plan) lp;
    int Ni[4], ni[4];
    INT jj;
    for (t = 0; t < d; t++) {
      Ni[t] = (int)N[t];
      ni[t] = (int)n[t];
    }
    NFFT(init_guru)(&lp, d, Ni, (int)M, ni, m,
              PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F
                   | FFTW_INIT | FFT_OUT_OF_PLACE,
              FFTW_ESTIMATE);
    for (jj = 0; jj < M * d; jj++)
      lp.x[jj] = x[jj];
    for (jj = 0; jj < M; jj++)
      lp.f[jj] = f[jj];
    NFFT(precompute_one_psi)(&lp);
    NFFT(adjoint)(&lp);
    {
      R errl = Y(test_rel_max_err)(got_fhat, lp.f_hat, Ntot);
      CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
    }
    NFFT(finalize)(&lp);
  }

  Y(plan_ng_destroy)(p);
  Y(free)(x);
  Y(free)(f);
  Y(free)(got_fhat);
  Y(the_planner_destroy)();
}

/* check_nfast_native_fast_3d(_adjoint) gate their native-vs-file check to
 * double. This runs the same comparison on the same reference files in
 * whatever precision is compiled, against Y(test_err_bound). The plans are
 * pinned to KAISER_BESSEL with n = 2*N, so the bound is queried at that
 * window and sigma = 2. */
void Y(check_nfast_float_accuracy)(void)
{
  int m = 7;
  R bound = Y(test_err_bound)(NFFT_WINDOW_KAISER_BESSEL, (R)m, K(2.0));

  /* forward: nfft_3d_10_10_10_10 (true/direct transform, double-generated). */
  {
    int d;
    INT *N, M;
    R *x;
    C *f_hat, *f;
    INT n[3];
    Y(plan_ng) *p;
    R err;

    CU_ASSERT_TRUE_FATAL(read_nd_case("data/nfft_3d_10_10_10_10.txt", &d, &N,
                                      &M, &x, &f_hat, &f));
    CU_ASSERT_EQUAL_FATAL(d, 3);
    n[0] = 2 * N[0];
    n[1] = 2 * N[1];
    n[2] = 2 * N[2];

    {
      C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
      p = Y(
           plan_ng_guru)(3, N, 0, n, M, m, NFFT_WINDOW_KAISER_BESSEL, x, f_hat,
                       got, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL_FATAL(p);

      Y(precompute)(p);
      Y(execute)(p);
      err = Y(test_rel_max_err)(got, f, M);
      CU_ASSERT(err < bound);

      Y(plan_ng_destroy)(p);
      Y(free)(got);
    }

    Y(free)(N);
    Y(free)(x);
    Y(free)(f_hat);
    Y(free)(f);
    Y(the_planner_destroy)();
  }

  /* adjoint: nfft_adjoint_3d_10_10_10_10, a separately generated reference,
   * not a round-trip of the forward case above. */
  {
    int d;
    INT *N, M;
    R *x;
    C *f_hat_ref, *f;
    INT n[3];
    Y(plan_ng) *p;
    R err;

    CU_ASSERT_TRUE_FATAL(read_nd_case("data/nfft_adjoint_3d_10_10_10_10.txt",
                                      &d, &N, &M, &x, &f_hat_ref, &f));
    CU_ASSERT_EQUAL_FATAL(d, 3);
    n[0] = 2 * N[0];
    n[1] = 2 * N[1];
    n[2] = 2 * N[2];

    {
      C *got_fhat = (C *)Y(malloc)((size_t)(N[0] * N[1] * N[2]) * sizeof(C));
      p = Y(plan_ng_guru)(3, N, 0, n, M, m, NFFT_WINDOW_KAISER_BESSEL, x,
                       got_fhat, f, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL_FATAL(p);

      Y(precompute)(p);
      Y(execute_adjoint)(p);
      err = Y(test_rel_max_err)(got_fhat, f_hat_ref, N[0] * N[1] * N[2]);
      CU_ASSERT(err < bound);

      Y(plan_ng_destroy)(p);
      Y(free)(got_fhat);
    }

    Y(free)(N);
    Y(free)(x);
    Y(free)(f_hat_ref);
    Y(free)(f);
    Y(the_planner_destroy)();
  }
}

/* Near the peak, where cancellation bites hardest, window_phi_hut against a
 * stable log-domain reference in the working precision. */
void Y(check_nfast_window_cancellation)(void)
{
  const int w = NFFT_WINDOW_KAISER_BESSEL;
  INT n = 512, N = 256;
  int m = 8, k;
  R b = KPI * (K(2.0) - (R)N / (R)n), xpk = (R)m * b;
  R lt = Y(bessel_i0_logtail)(xpk);
  for (k = 0; k <= 3; k++) /* the near-peak bins where a ~ m*b */
  {
    R t = K(2.0) * KPI * (R)k / (R)n, a = (R)m * SQRT(b * b - t * t);
    R got = Y(window_phi_hut)(w, n, N, m, (INT)k);
    /* stable log-domain reference: log I0(a) - log I0(xpk) is
     * (a - xpk) + logtail(a) - logtail(xpk), with no peak to cancel. */
    R ref = EXP((a - xpk) + Y(bessel_i0_logtail)(a) - lt);
    CU_ASSERT(FABS(got - ref) <= K(20.0) * Y(float_property)(NFFT_EPSILON)
                                      * (K(1.0) + FABS(ref)));
  }
}

void Y(check_nfast_window_apply)(void)
{
  const INT ntab[] = {16, 20, 128};
  const INT Ntab[] = {8, 10, 64};
  const int mtab[] = {4, 7, 8};
  int i;
  for (i = 0; i < 3; i++) {
    INT n = ntab[i], N = Ntab[i], k;
    int m = mtab[i];
    R *buf = (R *)Y(malloc)((size_t)N * sizeof(R));

    /* phi_hut_apply fills the band k in [-N/2, N/2) to match the single-arg
     * form. Under -ffast-math the two paths reassociate differently and differ
     * by ~1 ULP, so the contract is a relative tolerance, not equality. */
    Y(window_phi_hut_apply)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, -N / 2, buf, N);
    CU_ASSERT(FABS(buf[N / 2] - K(1.0)) <= K(1e-5));
    for (k = 0; k < N; k++) {
      R got = buf[k],
        ref = Y(window_phi_hut)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, k - N / 2);
      CU_ASSERT(FABS(got - ref)
                <= NFAST_LEGACY_REL_BOUND * FMAX(FABS(got), FABS(ref)));
      CU_ASSERT(buf[k] > K(0.0) && buf[k] <= K(1.0) + K(1e-5));
    }
    CU_ASSERT(buf[N / 2] > buf[N / 2 + 1]); /* monotone decay */

    /* The per-element fallback serves the raw windows too. */
    {
      R *gbuf = (R *)Y(malloc)((size_t)N * sizeof(R));
      INT kk;
      Y(window_phi_hut_apply)(NFFT_WINDOW_GAUSSIAN, n, N, m, -N / 2, gbuf, N);
      for (kk = 0; kk < N; kk++) {
        R single =
             Y(window_phi_hut)(NFFT_WINDOW_GAUSSIAN, n, N, m, kk - N / 2);
        CU_ASSERT(FABS(gbuf[kk] - single)
                  <= NFAST_LEGACY_REL_BOUND * (K(1.0) + FABS(single)));
      }
      Y(free)(gbuf);
    }

    /* phi_precompute: per-node psi matches single-arg window_phi at every
     * tap, under the same 1-ULP caveat. */
    {
      INT M = 4, j, lj;
      R x[4] = {K(0.1), K(0.37), K(0.5), K(0.83)};
      R *psi = (R *)Y(malloc)((size_t)M * (size_t)(2 * m + 2) * sizeof(R));
      Y(window_phi_precompute)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, x, 1, M, psi,
                            2 * m + 2);
      for (j = 0; j < M; j++) {
        INT c = LRINT(FLOOR((R)n * x[j]));
        for (lj = 0; lj <= 2 * m + 1; lj++) {
          INT idx = c - m + lj;
          R got = psi[j * (2 * m + 2) + lj],
            ref = Y(window_phi)(NFFT_WINDOW_KAISER_BESSEL, n, N, m,
                             x[j] - (R)idx / (R)n);
          CU_ASSERT(FABS(got - ref)
                    <= NFAST_LEGACY_REL_BOUND * FMAX(FABS(got), FABS(ref)));
        }
      }
      Y(free)(psi);
    }
    Y(free)(buf);
  }
}

void Y(check_nfast_native_window_select)(void)
{
  INT N, M, n;
  R *x;
  C *f_hat, *f;
  int w;
  CU_ASSERT_TRUE_FATAL(
       read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;
  for (w = NFFT_WINDOW_KAISER_BESSEL; w <= NFFT_WINDOW_SINC_POWER; w++) {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    Y(plan_ng) *p = Y(plan_ng_guru)(1, &N, 0, &n, M, 7, w, x, f_hat, got, 0u,
                              NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p) {
      Y(plan_ng_destroy)(p);
    }
    Y(free)(got);
  }
  /* Dirac and a garbage ordinal are both declined. */
  {
    int bad[2];
    int i;
    bad[0] = NFFT_WINDOW_DIRAC_DELTA;
    bad[1] = 99;
    for (i = 0; i < 2; i++) {
      Y(plan_ng) *pbad =
           Y(plan_ng_guru)(1, &N, 0, &n, M, 7, bad[i], x, f_hat, f, 0u,
                                   NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NULL(pbad);
    }
  }
  Y(free)(x);
  Y(free)(f_hat);
  Y(free)(f);
  Y(the_planner_destroy)();
}

void Y(check_nfast_window_accuracy)(void)
{
  INT N, M, n;
  int m = 7, wins[4], nw, i;
  R *x, sigma;
  C *f_hat, *f;
  CU_ASSERT_TRUE_FATAL(
       read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;
  sigma = (R)n / (R)N;
  nw = windows_from_env(wins);
  for (i = 0; i < nw; i++) {
    int w = wins[i];
    R bound = Y(test_err_bound)(w, (R)m, sigma);
    /* forward */
    {
      C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
      C *ref = (C *)Y(malloc)((size_t)M * sizeof(C));
      Y(plan_ng) *pn =
           Y(plan_ng_guru)(1, &N, 0, &n, M, m, w, x, f_hat, got, 0u,
                                 NFFT_ESTIMATE | NFFT_NO_DIRECT);
      Y(plan_ng) *pd =
           Y(plan_ng_guru)(1, &N, 0, &n, M, m, w, x, f_hat, ref, 0u,
                                 NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pn);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pd);
      Y(precompute)(pn);
      Y(execute)(pn);
      Y(precompute)(pd);
      Y(execute)(pd);
      CU_ASSERT(Y(test_rel_max_err)(got, ref, M) < bound);
      Y(plan_ng_destroy)(pn);
      Y(plan_ng_destroy)(pd);
      Y(free)(got);
      Y(free)(ref);
    }
    /* adjoint */
    {
      C *got = (C *)Y(malloc)((size_t)N * sizeof(C));
      C *ref = (C *)Y(malloc)((size_t)N * sizeof(C));
      Y(plan_ng) *pn = Y(plan_ng_guru)(1, &N, 0, &n, M, m, w, x, got, f, 0u,
                                 NFFT_ESTIMATE | NFFT_NO_DIRECT);
      Y(plan_ng) *pd = Y(plan_ng_guru)(1, &N, 0, &n, M, m, w, x, ref, f, 0u,
                                 NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pn);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pd);
      Y(precompute)(pn);
      Y(execute_adjoint)(pn);
      Y(precompute)(pd);
      Y(execute_adjoint)(pd);
      CU_ASSERT(Y(test_rel_max_err)(got, ref, N) < bound);
      Y(plan_ng_destroy)(pn);
      Y(plan_ng_destroy)(pd);
      Y(free)(got);
      Y(free)(ref);
    }
  }
  Y(free)(x);
  Y(free)(f_hat);
  Y(free)(f);
  Y(the_planner_destroy)();
}

/* The op counts are what the estimator sees, so they are pinned against the
 * kernel's loops rather than against a formula that could drift with it.
 * mul is an exact per-element count throughout; other is exact for rank 1
 * but for ranks 2, 3 and nd approximates the row-setup index arithmetic at a
 * flat rate per row rather than tracing it. */
void Y(check_nfast_deconv_ops)(void)
{
  planner *pl = Y(planner_create)();
  const int w = Y(get_window_id)();
  C fh = K(1.0), g = K(0.0);

  Y(deconv_solvers_register)(pl);

  { /* rank 1: N=16, n=32 -> mul 32, other (32-16) + 32 = 48 */
    const INT N = 16, n = 32;
    problem *p = Y(mkproblem_deconv)(1, &N, 0, &n, 6, w, 1, &fh, &g);
    plan *pln = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.mul, 32.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.fma, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.other, 48.0, 0.0);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }
  { /* rank 2: N=(4,8) n=(8,16) -> Ntot 32, ntot 128, rows 4
     * mul 96, other 128 + 64 + 2*1*4 = 200 */
    const INT N[2] = {4, 8}, n[2] = {8, 16};
    problem *p = Y(mkproblem_deconv)(2, N, 0, n, 6, w, 1, &fh, &g);
    plan *pln = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.mul, 96.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.fma, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.other, 200.0, 0.0);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }
  { /* rank 3: N=(4,4,8) n=(8,8,16) -> Ntot 128, ntot 1024, rows 16
     * mul 384, other 1024 + 256 + 2*2*16 = 1344 */
    const INT N[3] = {4, 4, 8}, n[3] = {8, 8, 16};
    problem *p = Y(mkproblem_deconv)(3, N, 0, n, 6, w, 1, &fh, &g);
    plan *pln = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.mul, 384.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.fma, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.other, 1344.0, 0.0);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }
  { /* rank 4: N=(2,4,4,8) n=(4,8,8,16) -> Ntot 256, ntot 4096, rows 32
     * mul 768, other 4096 + 512 + 2*3*32 = 4800 */
    const INT N[4] = {2, 4, 4, 8}, n[4] = {4, 8, 8, 16};
    problem *p = Y(mkproblem_deconv)(4, N, 0, n, 6, w, 1, &fh, &g);
    plan *pln = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.mul, 768.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.fma, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.other, 4800.0, 0.0);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }
  Y(planner_destroy)(pl);
}

/* CONV per rank: mul/add trace the per-tap complex-by-real accumulate
 * exactly; other adds a flat per-axis conv_runs term and, for rank >= 2,
 * omits the lower-order axis-combine arithmetic (an estimate, not a full
 * trace -- see the comment at each fill site). M = 10, m = 6 throughout, so
 * len = 14. */
void Y(check_nfast_conv_ops)(void)
{
  planner *pl = Y(planner_create)();
  const int w = Y(get_window_id)();
  R x0 = K(0.0);
  C g0 = K(0.0), f0 = K(0.0);

  Y(conv_solvers_register)(pl);

  /* x_sorted = 1: x is a probe pointer. The sorted solver declines, so the
   * caller-order solver is adopted untimed and nothing reads x. */
  { /* rank 1: taps = 140 -> mul 280, add 280, other 280 + 20 = 300 */
    const INT n = 8, N = 4, M = 10;
    problem *p = Y(mkproblem_conv)(1, &n, &N, M, 6, w, 1, 1, &x0, &g0, &f0);
    plan *pln = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.mul, 280.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 280.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.fma, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.other, 300.0, 0.0);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }
  { /* rank 2: taps = 1960 -> mul 3920, add 3920, other 3920 + 40 = 3960 */
    const INT n[2] = {8, 8}, N[2] = {4, 4}, M = 10;
    problem *p = Y(mkproblem_conv)(2, n, N, M, 6, w, 1, 1, &x0, &g0, &f0);
    plan *pln = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.mul, 3920.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 3920.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.fma, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.other, 3960.0, 0.0);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }
  { /* rank 3: taps = 27440 -> mul 54880, add 54880, other 54880 + 60 = 54940 */
    const INT n[3] = {8, 8, 8}, N[3] = {4, 4, 4}, M = 10;
    problem *p = Y(mkproblem_conv)(3, n, N, M, 6, w, 1, 1, &x0, &g0, &f0);
    plan *pln = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.mul, 54880.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 54880.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.fma, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.other, 54940.0, 0.0);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }
  { /* rank 4 (nd): taps = 384160 -> mul 768320, add 768320,
     * other 768320 + 80 = 768400 */
    const INT n[4] = {8, 8, 8, 8}, N[4] = {4, 4, 4, 4}, M = 10;
    problem *p = Y(mkproblem_conv)(4, n, N, M, 6, w, 1, 1, &x0, &g0, &f0);
    plan *pln = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.mul, 768320.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 768320.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.fma, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.other, 768400.0, 0.0);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }
  Y(planner_destroy)(pl);
}

/* The direct NDFT (rank 1, rank >= 2) and rank-0 solvers. mul/add trace two
 * complex-by-complex multiplies per tap exactly: the phase-weighted
 * accumulate and the w *= dw recurrence step that advances the phase for the
 * next tap. other is a per-tap load plus the block-start trig/reduced_omega
 * setup amortised over NDFT_RECURRENCE_BLOCK -- an estimate, not a full trace
 * (see the comment at the fill site). The rank-0 solver's forward is a pure
 * broadcast: no arithmetic, one store per node. */
void Y(check_nfast_ndft_ops)(void)
{
  planner *pl = Y(planner_create)();
  const int w = Y(get_window_id)();
  R x0 = K(0.0);
  C fh0 = K(0.0), f0 = K(0.0);

  Y(nfft_solver_ndft_1d_register)(pl);
  Y(nfft_solver_ndft_nd_register)(pl);
  Y(nfft_solver_rnk0_register)(pl);

  { /* rank 1: Ntot = 32, M = 10 -> taps = 320
     * mul 2560, add 1920, other 640 + 30 = 670 */
    const INT N = 32, n = 32, M = 10;
    problem *p =
         Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, w, 1, 0u, &x0, 0, &fh0, &f0);
    plan *pln = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.mul, 2560.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 1920.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.fma, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.other, 670.0, 0.0);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }
  { /* rank 2: N = (32, 4) -> Ntot = 128, M = 5 -> taps = 640
     * mul 5120, add 3840, other 1280 + 60 = 1340 */
    const INT N[2] = {32, 4}, n[2] = {32, 4}, M = 5;
    problem *p =
         Y(mkproblem_nfft)(2, N, 0, n, M, 6, w, 1, 0u, &x0, 0, &fh0, &f0);
    plan *pln = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.mul, 5120.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 3840.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.fma, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.other, 1340.0, 0.0);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }
  { /* rank 0: N = (1) elides to sz->rnk == 0; M = 12 -> broadcast only:
     * mul 0, add 0, other 12 */
    const INT N1 = 1, n1 = 1, M = 12;
    problem *p =
         Y(mkproblem_nfft)(1, &N1, 0, &n1, M, 6, w, 1, 0u, &x0, 1, &fh0, &f0);
    plan *pln = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.mul, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.fma, 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(pln->ops.other, 12.0, 0.0);
    Y(plan_destroy)(pln);
    Y(problem_destroy)(p);
  }
  Y(planner_destroy)(pl);
}

/* The composite's estimate is the sum of its three parts and nothing else.
 * Asserting the identity rather than a number keeps the case alive when any
 * child's counts change. */
void Y(check_nfast_composite_ops)(void)
{
  const INT N = 64, n = 128, M = 100;
  R *x = (R *)Y(malloc)((size_t)M * sizeof(R));
  C *f_hat = (C *)Y(malloc)((size_t)N * sizeof(C));
  C *f = (C *)Y(malloc)((size_t)M * sizeof(C));
  Y(plan_ng) *p;
  INT j;

  for (j = 0; j < M; j++)
    x[j] = (R)j / (R)M - K(0.5);

  p = Y(plan_ng_guru)(1, &N, 0, &n, M, 6, Y(get_window_id)(), x, f_hat, f, 0u,
                   NFFT_ESTIMATE);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);
  {
    opcnt top, deconv, conv, fft, sum;
    CU_ASSERT_TRUE_FATAL(
         Y(plan_ng_test_native_parts)(p, &top, &deconv, &conv, &fft));
    Y(ops_add)(&deconv, &conv, &sum);
    Y(ops_add2)(&fft, &sum);
    CU_ASSERT_DOUBLE_EQUAL(top.add, sum.add, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(top.mul, sum.mul, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(top.fma, sum.fma, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(top.other, sum.other, 1e-6);
  }
  Y(plan_ng_destroy)(p);
  Y(free)(f);
  Y(free)(f_hat);
  Y(free)(x);
}

/* 1 when perm is a permutation of 0..M-1 in ascending CONV key, ties in
 * ascending node index. */
static int is_stable_key_order(const R *x, const INT *n, int d, INT M,
                               const INT *perm)
{
  char *seen = (char *)Y(malloc)((size_t)M);
  INT k, prev_key = -1, prev_j = -1;
  int ok = 1;
  memset(seen, 0, (size_t)M);
  for (k = 0; k < M && ok; k++) {
    const INT j = perm[k];
    INT key;
    if (j < 0 || j >= M || seen[j]) {
      ok = 0;
      break;
    }
    seen[j] = 1;
    key = Y(conv_node_key)(x + j * d, n, d);
    if (key < prev_key || (key == prev_key && j < prev_j))
      ok = 0;
    prev_key = key;
    prev_j = j;
  }
  Y(free)(seen);
  return ok;
}

/* The CONV node key is the row-major index of the node's grid cell, each axis
 * shifted by n_t/2; the sort is stable in node index. */
void Y(check_nfast_conv_node_order)(void)
{
  { /* 1D: x = -0.5 is cell 0, x = 0.5 wraps to it */
    const INT n = 16;
    R x = K(-0.5);
    CU_ASSERT_EQUAL(Y(conv_node_key)(&x, &n, 1), (INT)0);
    x = K(0.0);
    CU_ASSERT_EQUAL(Y(conv_node_key)(&x, &n, 1), (INT)8);
    x = K(0.49);
    CU_ASSERT_EQUAL(Y(conv_node_key)(&x, &n, 1), (INT)15);
    x = K(0.5);
    CU_ASSERT_EQUAL(Y(conv_node_key)(&x, &n, 1), (INT)0);
  }
  { /* row-major: axis 0 most significant */
    const INT n2[2] = {8, 4}, n3[3] = {4, 4, 4};
    const R a[2] = {K(0.0), K(-0.25)};
    const R b[3] = {K(-0.5), K(-0.5), K(0.25)};
    const R c[3] = {K(0.25), K(0.0), K(-0.5)};
    CU_ASSERT_EQUAL(Y(conv_node_key)(a, n2, 2), (INT)17);
    CU_ASSERT_EQUAL(Y(conv_node_key)(b, n3, 3), (INT)3);
    CU_ASSERT_EQUAL(Y(conv_node_key)(c, n3, 3), (INT)56);
  }
  { /* sortedness, 1D */
    const INT n = 16;
    const R xs[5] = {K(-0.5), K(-0.2), K(-0.2), K(0.1), K(0.45)};
    const R xu[5] = {K(-0.5), K(0.1), K(-0.2), K(-0.2), K(0.45)};
    const R x1 = K(0.3);
    CU_ASSERT_EQUAL(Y(conv_nodes_sorted)(xs, &n, 1, 5), 1);
    CU_ASSERT_EQUAL(Y(conv_nodes_sorted)(xu, &n, 1, 5), 0);
    CU_ASSERT_EQUAL(Y(conv_nodes_sorted)(&x1, &n, 1, 1), 1);
  }
  { /* sorted by x0 alone is not sorted by cell once two nodes share a row */
    const INT n[2] = {8, 8};
    const R x[4] = {K(0.01), K(0.4), K(0.02), K(-0.4)};
    CU_ASSERT_EQUAL(Y(conv_nodes_sorted)(x, n, 2, 2), 0);
  }
  { /* sort: keys 52, 5, 52, 24, 0, 52 */
    const INT n[2] = {8, 8};
    const R x[12] = {K(0.3),  K(0.1),  K(-0.5), K(0.2),  K(0.3), K(0.1),
                     K(-0.1), K(-0.4), K(-0.5), K(-0.5), K(0.3), K(0.0)};
    const INT expect[6] = {4, 1, 3, 0, 2, 5};
    INT perm[6];
    int k;
    Y(conv_sort_nodes)(x, n, 2, 6, perm);
    for (k = 0; k < 6; k++)
      CU_ASSERT_EQUAL(perm[k], expect[k]);
  }
  { /* ntot = 2^21: three radix passes, and many ties */
    const INT n[3] = {128, 128, 128};
    const INT M = 5000;
    R *x = (R *)Y(malloc)((size_t)(3 * M) * sizeof(R));
    INT *perm = (INT *)Y(malloc)((size_t)M * sizeof(INT));
    INT k;
    for (k = 0; k < 3 * M; k++)
      x[k] = (R)((k * 7919) % 1000) / K(1000.0) - K(0.5);
    Y(conv_sort_nodes)(x, n, 3, M, perm);
    CU_ASSERT_TRUE(is_stable_key_order(x, n, 3, M, perm));
    Y(free)(perm);
    Y(free)(x);
  }
}

/* Calls conv_sort_nodes from thread 0 of an active omp parallel region,
 * nesting disabled and a higher thread count set inside. */
static void sort_in_parallel_call(const R *x, const INT *n, INT M, INT *perm)
{
#ifdef _OPENMP
  const int levels = omp_get_max_active_levels();
  omp_set_max_active_levels(1);
# pragma omp parallel num_threads(2)
  {
    omp_set_num_threads(8);
    if (omp_get_thread_num() == 0)
#endif
      Y(conv_sort_nodes)(x, n, 3, M, perm);
#ifdef _OPENMP
  }
  omp_set_max_active_levels(levels);
#endif
}

/* Same node set as above (three radix passes, many ties), sorted from inside
 * an active OpenMP region. */
void Y(check_nfast_sort_in_parallel)(void)
{
  const INT n[3] = {128, 128, 128};
  const INT M = 5000;
  R *x = (R *)Y(malloc)((size_t)(3 * M) * sizeof(R));
  INT *perm = (INT *)Y(malloc)((size_t)M * sizeof(INT));
  INT k;
  for (k = 0; k < 3 * M; k++)
    x[k] = (R)((k * 7919) % 1000) / K(1000.0) - K(0.5);
  sort_in_parallel_call(x, n, M, perm);
  CU_ASSERT_TRUE(is_stable_key_order(x, n, 3, M, perm));
  Y(free)(perm);
  Y(free)(x);
}

/* The guru's copy of x records whether the live axes are in CONV key order.
 * Paths that never read x record 1. */
void Y(check_nfast_nfft_x_sorted)(void)
{
  const int w = Y(get_window_id)();
  C fh[16], f[4];
  problem *p;

  { /* 1D, copy_x = 1 */
    const INT N = 16, n = 32;
    R xs[4] = {K(-0.5), K(-0.25), K(0.0), K(0.25)};
    R xu[4] = {K(0.25), K(-0.5), K(0.0), K(-0.25)};
    p = Y(mkproblem_nfft)(1, &N, 0, &n, 4, 6, w, 1, 0u, xs, 1, fh, f);
    CU_ASSERT_EQUAL(((problem_nfft *)p)->x_sorted, 1);
    Y(problem_destroy)(p);
    p = Y(mkproblem_nfft)(1, &N, 0, &n, 4, 6, w, 1, 0u, xu, 1, fh, f);
    CU_ASSERT_EQUAL(((problem_nfft *)p)->x_sorted, 0);
    Y(problem_destroy)(p);
    /* borrowed: never read, so 1 */
    p = Y(mkproblem_nfft)(1, &N, 0, &n, 4, 6, w, 1, 0u, xu, 0, fh, f);
    CU_ASSERT_EQUAL(((problem_nfft *)p)->x_sorted, 1);
    Y(problem_destroy)(p);
  }
  { /* a unit axis is elided first; axis 1 is out of order and ignored */
    const INT N[2] = {16, 1}, n[2] = {32, 1};
    R x[8] = {K(-0.5), K(0.4), K(-0.25), K(-0.4),
              K(0.0),  K(0.1), K(0.25),  K(-0.1)};
    p = Y(mkproblem_nfft)(2, N, 0, n, 4, 6, w, 1, 0u, x, 1, fh, f);
    CU_ASSERT_EQUAL(((problem_nfft *)p)->x_sorted, 1);
    Y(problem_destroy)(p);
  }
}

/* One CONV candidate by family: registration id 0 is caller order, 1 is
 * sorted. Each rank's registrar registers both and ranks are disjoint, so the
 * id picks one. NULL when that solver declines. */
static plan *conv_family(planner *pl, const problem *p, int sorted)
{
  plan *cands[8];
  unsigned slv[8];
  plan *got = 0;
  int nc = Y(planner_candidates)(pl, p, cands, slv, 8), i;
  for (i = 0; i < nc; i++) {
    if (got == 0 && pl->slvdescs[slv[i]].reg_id == sorted)
      got = cands[i];
    else
      Y(plan_destroy)(cands[i]);
  }
  return got;
}

static int plan_prints(const plan *pln, const char *needle)
{
  char buf[256];
  printer *pr = Y(printer_create_str)(buf);
  pr->print(pr, "%p", pln);
  Y(printer_destroy)(pr);
  return strstr(buf, needle) != 0;
}

static void conv_sorted_case(int d, const INT *n, const INT *N, INT M)
{
  const int w = Y(get_window_id)(), m = 6;
  planner *pl = Y(planner_create)();
  INT ntot = 1, i;
  int t;
  R tol = K(1.0e-12);
  R *x;
  C *g_in, *g_ref, *g_got, *f_in, *f_ref, *f_got;
  problem *pu, *ps, *pp;
  plan *u, *s, *d0, *d1;

  /* max(1e-12, 1e4*eps): the two orders reach psi through different blocks
   * and accumulate the adjoint in a different node order. */
  if ((R)1.0e4 * EPSILON > tol)
    tol = (R)1.0e4 * EPSILON;
  for (t = 0; t < d; t++)
    ntot *= n[t];
  x = (R *)Y(malloc)((size_t)(M * d) * sizeof(R));
  g_in = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  g_ref = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  g_got = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  f_in = (C *)Y(malloc)((size_t)M * sizeof(C));
  f_ref = (C *)Y(malloc)((size_t)M * sizeof(C));
  f_got = (C *)Y(malloc)((size_t)M * sizeof(C));
  for (i = 0; i < M * d; i++)
    x[i] = (R)((i * 7919) % 1000) / K(1000.0) - K(0.5);
  CU_ASSERT_EQUAL_FATAL(Y(conv_nodes_sorted)(x, n, d, M), 0);
  for (i = 0; i < ntot; i++)
    g_in[i] = (R)(i + 1) / (R)ntot - K(0.25) * II;
  for (i = 0; i < M; i++)
    f_in[i] = (R)(M - i) / (R)M + K(0.75) * II;

  Y(conv_solvers_register)(pl);
  pu = Y(mkproblem_conv)(d, n, N, M, m, w, 1, 0, x, g_in, f_ref);
  ps = Y(mkproblem_conv)(d, n, N, M, m, w, 1, 0, x, g_in, f_got);
  u = conv_family(pl, pu, 0);
  s = conv_family(pl, ps, 1);
  CU_ASSERT_PTR_NOT_NULL_FATAL(u);
  CU_ASSERT_PTR_NOT_NULL_FATAL(s);
  CU_ASSERT_TRUE(plan_prints(s, " sorted"));
  CU_ASSERT_FALSE(plan_prints(u, " sorted"));

  /* forward */
  Y(plan_awake)(u, PLNR_AWAKE);
  Y(plan_awake)(s, PLNR_AWAKE);
  u->adt->apply(u, pu);
  s->adt->apply(s, ps);
  for (i = 0; i < M; i++)
    CU_ASSERT(CABS(f_got[i] - f_ref[i]) <= tol * (CABS(f_ref[i]) + K(1.0)));

  /* adjoint, same plans */
  ((problem_conv *)pu)->f = f_in;
  ((problem_conv *)pu)->g = g_ref;
  ((problem_conv *)ps)->f = f_in;
  ((problem_conv *)ps)->g = g_got;
  u->adt->apply_adjoint(u, pu);
  s->adt->apply_adjoint(s, ps);
  for (i = 0; i < ntot; i++)
    CU_ASSERT(CABS(g_got[i] - g_ref[i]) <= tol * (CABS(g_ref[i]) + K(1.0)));

  /* A full cycle through AWAKE_ZERO rebuilds perm, u and psi: the forward
   * repeats exactly. f_in is scratch from here on. */
  Y(plan_awake)(s, PLNR_SLEEPY);
  Y(plan_awake)(s, PLNR_AWAKE_ZERO);
  Y(plan_awake)(s, PLNR_SLEEPY);
  Y(plan_awake)(s, PLNR_AWAKE);
  ((problem_conv *)ps)->g = g_in;
  ((problem_conv *)ps)->f = f_in;
  s->adt->apply(s, ps);
  for (i = 0; i < M; i++)
    CU_ASSERT(f_in[i] == f_got[i]);

  /* presorted: the sorted solver declines, caller order serves */
  pp = Y(mkproblem_conv)(d, n, N, M, m, w, 1, 1, x, g_in, f_ref);
  d0 = conv_family(pl, pp, 0);
  d1 = conv_family(pl, pp, 1);
  CU_ASSERT_PTR_NOT_NULL(d0);
  CU_ASSERT_PTR_NULL(d1);
  if (d0)
    Y(plan_destroy)(d0);
  if (d1)
    Y(plan_destroy)(d1);

  Y(plan_destroy)(u);
  Y(plan_destroy)(s);
  Y(problem_destroy)(pu);
  Y(problem_destroy)(ps);
  Y(problem_destroy)(pp);
  Y(planner_destroy)(pl);
  Y(free)(f_got);
  Y(free)(f_ref);
  Y(free)(f_in);
  Y(free)(g_got);
  Y(free)(g_ref);
  Y(free)(g_in);
  Y(free)(x);
}

/* The sorted CONV solver per rank: same results as caller order, declines on
 * presorted nodes, survives a full awake cycle. */
void Y(check_nfast_conv_sorted)(void)
{
  {
    const INT n[1] = {128}, N[1] = {64};
    conv_sorted_case(1, n, N, 500);
  }
  {
    const INT n[2] = {48, 64}, N[2] = {16, 24};
    conv_sorted_case(2, n, N, 300);
  }
  {
    const INT n[3] = {32, 32, 32}, N[3] = {8, 10, 12};
    conv_sorted_case(3, n, N, 200);
  }
  {
    const INT n[4] = {24, 24, 24, 24}, N[4] = {6, 6, 8, 8};
    conv_sorted_case(4, n, N, 60);
  }
}

/* Smallest power-of-two axis whose rank-d grid exceeds that rank's sort
 * threshold; 0 when the rank never sorts at ESTIMATE. */
static INT conv_axis_over_threshold(int d)
{
  const double limit = Y(conv_sort_grid_bytes)(d);
  INT a = 2;
  if (limit == HUGE_VAL)
    return 0;
  for (;;) {
    double bytes = (double)sizeof(C);
    int t;
    for (t = 0; t < d; t++)
      bytes *= (double)a;
    if (bytes > limit)
      return a;
    a *= 2;
  }
}

/* Both families' other count on a rank-d grid of axis a, N = 4, M = 10. */
static void conv_order_case(planner *pl, int d, INT a, double u_other,
                            double s_other)
{
  const int w = Y(get_window_id)();
  INT n[4], N[4];
  R x0 = K(0.0);
  C g0 = K(0.0), f0 = K(0.0);
  problem *p;
  plan *u, *s;
  int t;
  for (t = 0; t < d; t++) {
    n[t] = a;
    N[t] = 4;
  }
  p = Y(mkproblem_conv)(d, n, N, 10, 6, w, 1, 0, &x0, &g0, &f0);
  u = conv_family(pl, p, 0);
  s = conv_family(pl, p, 1);
  CU_ASSERT_PTR_NOT_NULL_FATAL(u);
  CU_ASSERT_PTR_NOT_NULL_FATAL(s);
  CU_ASSERT_DOUBLE_EQUAL(u->ops.other, u_other, 0.0);
  CU_ASSERT_DOUBLE_EQUAL(s->ops.other, s_other, 0.0);
  Y(plan_destroy)(u);
  Y(plan_destroy)(s);
  Y(problem_destroy)(p);
}

/* Node order's share of the counts: sorted pays 2M for the permuted f access;
 * caller order pays one cold fetch per tap once the grid exceeds the rank's
 * threshold and the nodes are not presorted. Base counts as in
 * check_nfast_conv_ops: m = 6, len = 14, M = 10. Under is half the smallest
 * axis over the threshold. */
void Y(check_nfast_conv_order_ops)(void)
{
  planner *pl = Y(planner_create)();
  const int w = Y(get_window_id)();
  const INT M = 10;
  const INT a1 = conv_axis_over_threshold(1), a2 = conv_axis_over_threshold(2),
            a3 = conv_axis_over_threshold(3), a4 = conv_axis_over_threshold(4);
  R x0 = K(0.0);
  C g0 = K(0.0), f0 = K(0.0);
  problem *p;
  plan *u, *s;

  Y(conv_solvers_register)(pl);

  { /* rank 1, n = 8: taps = 140, base other 300 */
    const INT n = 8, N = 4;
    p = Y(mkproblem_conv)(1, &n, &N, M, 6, w, 1, 0, &x0, &g0, &f0);
    s = conv_family(pl, p, 1);
    CU_ASSERT_PTR_NOT_NULL_FATAL(s);
    CU_ASSERT_DOUBLE_EQUAL(s->ops.mul, 280.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(s->ops.add, 280.0, 0.0);
    Y(plan_destroy)(s);
    Y(problem_destroy)(p);
    conv_order_case(pl, 1, n, 300.0, 320.0);
  }
  if (a1) {
    const INT N = 4;
    conv_order_case(pl, 1, a1 / 2, 300.0, 320.0);
    conv_order_case(pl, 1, a1, 440.0, 320.0);
    /* over, presorted: no penalty, sorted declines */
    p = Y(mkproblem_conv)(1, &a1, &N, M, 6, w, 1, 1, &x0, &g0, &f0);
    u = conv_family(pl, p, 0);
    CU_ASSERT_PTR_NOT_NULL_FATAL(u);
    CU_ASSERT_DOUBLE_EQUAL(u->ops.other, 300.0, 0.0);
    CU_ASSERT_PTR_NULL(conv_family(pl, p, 1));
    Y(plan_destroy)(u);
    Y(problem_destroy)(p);
  }
  if (a2) { /* taps = 1960, base other 3960 */
    conv_order_case(pl, 2, a2 / 2, 3960.0, 3980.0);
    conv_order_case(pl, 2, a2, 5920.0, 3980.0);
  }
  if (a3) { /* taps = 27440, base other 54940 */
    conv_order_case(pl, 3, a3 / 2, 54940.0, 54960.0);
    conv_order_case(pl, 3, a3, 82380.0, 54960.0);
  }
  /* rank 4 uses rank 3's value */
  CU_ASSERT_DOUBLE_EQUAL(Y(conv_sort_grid_bytes)(4),
                         Y(conv_sort_grid_bytes)(3),
                         0.0);
  if (a4) { /* taps = 384160, base other 768400 */
    conv_order_case(pl, 4, a4 / 2, 768400.0, 768420.0);
    conv_order_case(pl, 4, a4, 1152560.0, 768420.0);
  }
  Y(planner_destroy)(pl);

  /* ESTIMATE picks by these counts; no awake, so the probe x is never read */
  pl = Y(planner_create)();
  pl->flags.u |= PLNR_ESTIMATE;
  Y(conv_solvers_register)(pl);
  {
    const INT n = 8, N = 4;
    p = Y(mkproblem_conv)(1, &n, &N, M, 6, w, 1, 0, &x0, &g0, &f0);
    u = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(u);
    CU_ASSERT_FALSE(plan_prints(u, " sorted"));
    Y(plan_destroy)(u);
    Y(problem_destroy)(p);
  }
  if (a1) {
    const INT N = 4;
    p = Y(mkproblem_conv)(1, &a1, &N, M, 6, w, 1, 0, &x0, &g0, &f0);
    u = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(u);
    CU_ASSERT_TRUE(plan_prints(u, " sorted"));
    Y(plan_destroy)(u);
    Y(problem_destroy)(p);
    p = Y(mkproblem_conv)(1, &a1, &N, M, 6, w, 1, 1, &x0, &g0, &f0);
    u = Y(planner_mkplan)(pl, p);
    CU_ASSERT_PTR_NOT_NULL_FATAL(u);
    CU_ASSERT_FALSE(plan_prints(u, " sorted"));
    Y(plan_destroy)(u);
    Y(problem_destroy)(p);
  }
  Y(planner_destroy)(pl);
}

/* Both CONV families race here, one thread; tables live only while awake, so
 * the race holds one family's tables at a time (ADR-0009). */
void Y(check_nfast_conv_race_memory)(void)
{
  const INT N[2] = {128, 128}, n[2] = {256, 256}, M = 65536;
  const int w = Y(get_window_id)(), d = 2;
  const size_t len = 14;
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
    planner *pl = Y(planner_create)();
    problem *p = Y(mkproblem_conv)(d, n, N, M, 6, w, 1, 0, x, g, f);
    plan *pln;
    Y(conv_solvers_register)(pl);
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
  CU_ASSERT_EQUAL(nmeasured, 2u);
  CU_ASSERT(peak < tables + tables / 2);

  free(f);
  free(g);
  free(x);
}
