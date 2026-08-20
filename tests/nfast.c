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

#include "config.h" /* ABS_SRCDIR */
#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "nfast.h"

/* Native-fast and legacy-fast are the same algorithm, so they agree to
 * round-off: ~1e-12 in double/long-double, but only ~1e-4 in float, where
 * eps ~1e-7 accumulates over the deconv+FFT+conv pipeline. */
#if defined(NFFT_SINGLE)
#define NFAST_LEGACY_REL_BOUND ((R)1e-4)
#else
#define NFAST_LEGACY_REL_BOUND ((R)1e-12)
#endif

/* Runtime per-window NFFT accuracy bound, one case per window; the a/b
 * calibration follows err_trafo in tests/nfft.c. m/sigma are runtime, a/b are
 * precision-fixed. The 56*eps floor covers the fast pipeline's DECONV+CONV
 * round-off accumulation. */
static R err_bound(int window, R m, R s) {
  R eps = Y(float_property)(NFFT_EPSILON), a, b, err;
  switch (window) {
  case NFFT_WINDOW_GAUSSIAN:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    /* Gaussian only: err_trafo's a=0.41 is too tight for the composed
     * pipeline at m=7/sigma=2 (measured ~3.09e-7 against a ~1.76e-7 bound).
     * a=1.0 gives ~1.4x headroom. */
    a = K(1.0);
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
    b = K(3000.0);
#else
    a = K(1.5);
    b = K(50.0);
#endif
    err = KPI * (SQRT(m) + m) * SQRT(SQRT(K(1.0) - K(1.0) / K(2.0))) * EXP(-K2PI * m * SQRT(K(1.0) - K(1.0) / K(2.0)));
    break;
  }
  return FMAX(FMAX(a * err, b * eps), K(56.0) * eps);
}

/* Select which implemented windows the runtime accuracy loop runs.
 * NFAST_WINDOWS = comma list of kb|kaiserbessel,gaussian,bspline,sinc (case-
 * insensitive); unset or "all" => all four. Unknown token => abort loudly. */
static int win_token(const char *t) {
  if (!strcasecmp(t, "kb") || !strcasecmp(t, "kaiserbessel"))
    return NFFT_WINDOW_KAISER_BESSEL;
  if (!strcasecmp(t, "gaussian")) return NFFT_WINDOW_GAUSSIAN;
  if (!strcasecmp(t, "bspline") || !strcasecmp(t, "b_spline"))
    return NFFT_WINDOW_B_SPLINE;
  if (!strcasecmp(t, "sinc") || !strcasecmp(t, "sinc_power"))
    return NFFT_WINDOW_SINC_POWER;
  return -1;
}
static int windows_from_env(int *out) {
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
    for (tok = strtok_r(buf, ",", &save); tok; tok = strtok_r(NULL, ",", &save)) {
      int id = win_token(tok);
      CU_ASSERT_FATAL(id >= 0);
      CU_ASSERT_FATAL(n < 4);
      out[n++] = id;
    }
  }
  return n;
}

void Y(check_nfast_window_id)(void) {
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

void Y(check_nfast_deconv_problem)(void) {
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

void Y(check_nfast_conv_problem)(void) {
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

/* The KB vtable is peak-normalized, so these assertions check only positivity
 * and monotone decay, which a uniform positive scale preserves. Exact
 * normalization is checked in check_nfast_window_normalized. */
void Y(check_nfast_window_vtable)(void) {
  INT N = 16, n = 32;
  int m = 6, k;
  /* phi_hut is even and real for |k| <= N/2. */
  for (k = -(int)(N / 2); k < (int)(N / 2); k++) {
    R vh = Y(window_phi_hut)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, (INT)k);
    CU_ASSERT(vh > K(0.0));
  }
  /* phi is peaked at 0 and decays. */
  CU_ASSERT(Y(window_phi)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, K(0.0)) >
            Y(window_phi)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, K(0.3)));
  /* Dirac and out-of-range ordinals return 0. */
  CU_ASSERT(Y(window_phi_hut)(NFFT_WINDOW_GAUSSIAN, n, N, m, 0) > K(0.0));
  CU_ASSERT_EQUAL(Y(window_phi_hut)(NFFT_WINDOW_DIRAC_DELTA, n, N, m, 0), K(0.0));
  CU_ASSERT_EQUAL(Y(window_phi_hut)(99, n, N, m, 0), K(0.0));
}

/* phi_hut(0) = 1, decay is monotone, and ratios equal the legacy-shape
 * ratios because the scale is uniform. */
void Y(check_nfast_window_normalized)(void) {
  const int w = NFFT_WINDOW_KAISER_BESSEL;
  INT n = 32, N = 16;
  int m = 6, k;
  CU_ASSERT(FABS(Y(window_phi_hut)(w, n, N, m, (INT)0) - K(1.0)) <= K(1e-5));
  for (k = -(int)(N / 2); k < (int)(N / 2); k++)
    CU_ASSERT(Y(window_phi_hut)(w, n, N, m, (INT)k) > K(0.0) && Y(window_phi_hut)(w, n, N, m, (INT)k) <= K(1.0) + K(1e-5));
  /* uniform-scale invariance: phi_hut(k1)/phi_hut(k2) is scale-free, so it
   * equals I0(a1)/I0(a2) recomputed from the unscaled bessel directly. */
  {
    R b = KPI * (K(2.0) - (R)N / (R)n);
    R t1 = K(2.0) * KPI * K(1.0) / (R)n, t2 = K(2.0) * KPI * K(3.0) / (R)n;
    R a1 = (R)m * SQRT(b * b - t1 * t1), a2 = (R)m * SQRT(b * b - t2 * t2);
    R got = Y(window_phi_hut)(w, n, N, m, (INT)1) / Y(window_phi_hut)(w, n, N, m, (INT)3);
    R ref = Y(bessel_i0)(a1) / Y(bessel_i0)(a2);
    CU_ASSERT(FABS(got - ref) <= K(1e-4) * (K(1.0) + FABS(ref)));
  }
  CU_ASSERT(Y(window_phi)(w, n, N, m, K(0.0)) > Y(window_phi)(w, n, N, m, K(0.3)));
  /* Gaussian is raw (not peak-normalized): phi_hut(0) == 1 exactly. */
  CU_ASSERT(FABS(Y(window_phi_hut)(NFFT_WINDOW_GAUSSIAN, n, N, m, 0) - K(1.0)) <= K(1e-5));
}

/* Each window is checked against an independent recomputation of its infft.h
 * macro math: KB by a scale-free ratio because it is peak-normalized, the
 * other three by relative equality because they are raw. */
void Y(check_nfast_window_all)(void) {
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
    R got = Y(window_phi_hut)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, 1) / Y(window_phi_hut)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, 3);
    R ref = Y(bessel_i0)(a1) / Y(bessel_i0)(a2);
    CU_ASSERT(FABS(got - ref) <= tol * (K(1.0) + FABS(ref)));
  }
  {
    R b = (K(2.0) * sigma) / (K(2.0) * sigma - K(1.0)) * ((R)m / KPI);
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
    CU_ASSERT(FABS(Y(window_phi_hut)(NFFT_WINDOW_B_SPLINE, n, N, m, 0) - ref0) <= tol * (K(1.0) + FABS(ref0)));
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
      R arg = (K(2.0) * (R)m * (R)k) / ((K(2.0) * sigma - K(1.0)) * (R)n / sigma) + (R)m;
      R ref = Y(bsplines)((INT)(2 * m), arg);
      R got = Y(window_phi_hut)(NFFT_WINDOW_SINC_POWER, n, N, m, k);
      CU_ASSERT(FABS(got - ref) <= tol * (K(1.0) + FABS(ref)));
    }
    {
      R x = K(0.02);
      R ref = ((R)n / sigma) * (K(2.0) * sigma - K(1.0)) / (K(2.0) * (R)m) * POW(Y(sinc)(KPI * (R)n / sigma * x * (K(2.0) * sigma - K(1.0)) / (K(2.0) * (R)m)), K(2.0) * (R)m) / (R)n;
      R got = Y(window_phi)(NFFT_WINDOW_SINC_POWER, n, N, m, x);
      CU_ASSERT(FABS(got - ref) <= tol * (K(1.0) + FABS(ref)));
    }
  }
  /* Dirac and out-of-range ordinals return 0. */
  CU_ASSERT_EQUAL(Y(window_phi_hut)(NFFT_WINDOW_DIRAC_DELTA, n, N, m, 0), K(0.0));
  CU_ASSERT_EQUAL(Y(window_phi)(NFFT_WINDOW_DIRAC_DELTA, n, N, m, K(0.0)), K(0.0));
  CU_ASSERT_EQUAL(Y(window_phi_hut)(99, n, N, m, 0), K(0.0));
}

/* window is in the wisdom key: two NFFT problems differing only in window get
 * distinct keys. */
void Y(check_nfast_window_key)(void) {
  planner *pl = Y(planner_create)();
  INT N = 16, n = 32, M = 50;
  R x = K(0.1);
  C fh = K(1.0), f = K(0.0);
  problem *pkb, *pg;
  md5sig skb, sgb;
  pkb = Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_KAISER_BESSEL, +1, 0u,
                          &x, 0, &fh, &f);
  pg = Y(mkproblem_nfft)(1, &N, 0, &n, M, 6, NFFT_WINDOW_GAUSSIAN, +1, 0u,
                         &x, 0, &fh, &f);
  Y(problem_md5)
  (pl, pkb, skb);
  Y(problem_md5)
  (pl, pg, sgb);
  CU_ASSERT(skb[0] != sgb[0] || skb[1] != sgb[1] || skb[2] != sgb[2] ||
            skb[3] != sgb[3]);
  Y(problem_destroy)
  (pkb);
  Y(problem_destroy)
  (pg);
  Y(planner_destroy)
  (pl);
}

/* The DECONV solver, planned directly through planner_mkplan (there is no
 * deconv guru). Forward, adjoint and the type-II shift are each asserted as
 * values from a clean input, never as a round-trip: a round-trip also passes
 * for the inverse, which the adjoint of a real diagonal is not.
 * Uses the process-global planner, so it tears it down at the end. */
void Y(check_nfast_deconv_solver)(void) {
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
void Y(check_nfast_deconv_1d_general)(void)
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
void Y(check_nfast_deconv_2d_general)(void) {
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
void Y(check_nfast_deconv_3d_general)(void) {
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
void Y(check_nfast_deconv_nd_general)(void) {
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

/* The CONV solver, planned directly through planner_mkplan (there is no conv
 * guru). Forward and adjoint are asserted as values recomputed from
 * Y(window_phi) and the wrap formula, not from the solver's own psi table and
 * not as a round-trip; the round-trip below is only a cross-check. x is chosen
 * so node 0's support [c-m, c+m+1] neither wraps past n nor self-collides, so
 * each wrapped neighbor index receives exactly one psi contribution.
 * Uses the process-global planner, so it tears it down at the end. */
void Y(check_nfast_conv_solver)(void) {
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

/* Read a d == 1 case file: d, N, M, x[M], f_hat[N] as "re im" pairs, f[M] as
 * "re im" pairs. */
static int read_1d_case(const char *rel, INT *N, INT *M, R **x, C **f_hat,
                        C **f) {
  char path[4096];
  FILE *fp;
  int d;
  long v;
  INT j;
  snprintf(path, sizeof path, "%s/tests/%s", ABS_SRCDIR, rel);
  fp = fopen(path, "r");
  if (!fp) return 0;
  if (fscanf(fp, "%d", &d) != 1 || d != 1) {
    fclose(fp);
    return 0;
  }
  if (fscanf(fp, "%ld", &v) != 1) {
    fclose(fp);
    return 0;
  }
  *N = (INT)v;
  if (fscanf(fp, "%ld", &v) != 1) {
    fclose(fp);
    return 0;
  }
  *M = (INT)v;
  *x = (R *)Y(malloc)((size_t)*M * sizeof(R));
  for (j = 0; j < *M; j++) {
    double dv;
    if (fscanf(fp, "%lf", &dv) != 1) {
      fclose(fp);
      return 0;
    }
    (*x)[j] = (R)dv;
  }
  *f_hat = (C *)Y(malloc)((size_t)*N * sizeof(C));
  for (j = 0; j < *N; j++) {
    double re, im;
    if (fscanf(fp, "%lf %lf", &re, &im) != 2) {
      fclose(fp);
      return 0;
    }
    (*f_hat)[j] = (R)re + II * (R)im;
  }
  *f = (C *)Y(malloc)((size_t)*M * sizeof(C));
  for (j = 0; j < *M; j++) {
    double re, im;
    if (fscanf(fp, "%lf %lf", &re, &im) != 2) {
      fclose(fp);
      return 0;
    }
    (*f)[j] = (R)re + II * (R)im;
  }
  fclose(fp);
  return 1;
}

/* Read a case file of any rank: d, N[0..d-1], M, x[d*M] flat node-major
 * (x[j*d+t]), f_hat[prod N] as "re im" pairs, f[M] as "re im" pairs. */
static int read_nd_case(const char *rel, int *d, INT **N, INT *M, R **x,
                        C **f_hat, C **f) {
  char path[4096];
  FILE *fp;
  int t;
  long v;
  INT j, Ntot;
  snprintf(path, sizeof path, "%s/tests/%s", ABS_SRCDIR, rel);
  fp = fopen(path, "r");
  if (!fp) return 0;
  if (fscanf(fp, "%d", d) != 1 || *d < 1) {
    fclose(fp);
    return 0;
  }
  *N = (INT *)Y(malloc)((size_t)*d * sizeof(INT));
  for (t = 0; t < *d; t++) {
    if (fscanf(fp, "%ld", &v) != 1) {
      fclose(fp);
      return 0;
    }
    (*N)[t] = (INT)v;
  }
  if (fscanf(fp, "%ld", &v) != 1) {
    fclose(fp);
    return 0;
  }
  *M = (INT)v;
  *x = (R *)Y(malloc)((size_t)(*M) * (size_t)(*d) * sizeof(R));
  for (j = 0; j < (*M) * (*d); j++) {
    double dv;
    if (fscanf(fp, "%lf", &dv) != 1) {
      fclose(fp);
      return 0;
    }
    (*x)[j] = (R)dv;
  }
  for (t = 0, Ntot = 1; t < *d; t++)
    Ntot *= (*N)[t];
  *f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  for (j = 0; j < Ntot; j++) {
    double re, im;
    if (fscanf(fp, "%lf %lf", &re, &im) != 2) {
      fclose(fp);
      return 0;
    }
    (*f_hat)[j] = (R)re + II * (R)im;
  }
  *f = (C *)Y(malloc)((size_t)(*M) * sizeof(C));
  for (j = 0; j < *M; j++) {
    double re, im;
    if (fscanf(fp, "%lf %lf", &re, &im) != 2) {
      fclose(fp);
      return 0;
    }
    (*f)[j] = (R)re + II * (R)im;
  }
  fclose(fp);
  return 1;
}

static R rel_max_errC(const C *a, const C *b, INT len) {
  R num = (R)0, den = (R)0;
  INT j;
  for (j = 0; j < len; j++) {
    R e = CABS(a[j] - b[j]);
    if (e > num) num = e;
    if (CABS(b[j]) > den) den = CABS(b[j]);
  }
  return den > (R)0 ? num / den : num;
}

/* The composed native fast NFFT on the 1D even-N reference case, planned
 * under NFFT_NO_DIRECT so it is the sole survivor. The fast transform is an
 * approximation, so both directions are checked to 1e-5, not the 1e-9 direct
 * tolerance. The adjoint is compared against its own reference file: a
 * round-trip could not distinguish a correct adjoint from one off by a
 * constant that cancels. */
void Y(check_nfast_native_fast_accuracy)(void) {
  INT N, M, Na, Ma;
  int m = 7;
  R *x, *xa;
  C *f_hat, *f, *f_hat_ref, *fa;
  INT n;
  Y(plan_ng) * p;
  R err;

  CU_ASSERT_TRUE_FATAL(
      read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;

  /* forward */
  {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    p = Y(plan_ng_guru)(1, &N, 0, &n, M, m, Y(get_window_id)(), x, f_hat,
                        got, 0u,
                        NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    Y(precompute)
    (p);
    Y(execute)
    (p);
    err = rel_max_errC(got, f, M);
    CU_ASSERT(err < (R)1e-5);

    {
      FILE *tf = tmpfile();
      CU_ASSERT_PTR_NOT_NULL_FATAL(tf);
      Y(fprint_plan)
      (p, tf);
      {
        long len;
        char *buf;
        rewind(tf);
        fseek(tf, 0, SEEK_END);
        len = ftell(tf);
        rewind(tf);
        buf = (char *)Y(malloc)((size_t)len + 1);
        if (len > 0)
          CU_ASSERT_EQUAL(fread(buf, 1, (size_t)len, tf), (size_t)len);
        buf[len] = '\0';
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "nfft_solver_fast_native"));
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "deconv"));
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "conv"));
        Y(free)
        (buf);
      }
      fclose(tf);
    }

    Y(plan_ng_destroy)
    (p);
    Y(free)
    (got);
  }

  /* adjoint: a separately generated reference file, not a round-trip. */
  CU_ASSERT_TRUE_FATAL(read_1d_case("data/nfft_adjoint_1d_20_50.txt", &Na,
                                    &Ma, &xa, &f_hat_ref, &fa));
  CU_ASSERT_EQUAL(Na, N);
  CU_ASSERT_EQUAL(Ma, M);
  {
    C *got_fhat = (C *)Y(malloc)((size_t)Na * sizeof(C));
    Y(plan_ng) * pa;
    pa = Y(plan_ng_guru)(1, &Na, 0, &n, Ma, m, Y(get_window_id)(), xa,
                         got_fhat, fa, 0u,
                         NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pa);
    Y(precompute)
    (pa);
    Y(execute_adjoint)
    (pa);
    err = rel_max_errC(got_fhat, f_hat_ref, Na);
    CU_ASSERT(err < (R)1e-5);
    Y(plan_ng_destroy)
    (pa);
    Y(free)
    (got_fhat);
  }

  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(free)
  (xa);
  Y(free)
  (f_hat_ref);
  Y(free)
  (fa);

  Y(the_planner_destroy)
  ();
}

/* Isolation tests for NFFT_NO_FAST_NATIVE, which gates the composed native
 * fast. Tests that need a specific winner pin it with a disable flag; the
 * rest stay winner-agnostic.
 *
 * With the direct natives excluded, the composed native fast is the sole
 * surviving candidate for a 1D even-N problem. */
void Y(check_nfast_native_tree)(void) {
  INT N, M;
  int m = 7;
  R *x;
  C *f_hat, *f;
  INT n;
  Y(plan_ng) * p;

  CU_ASSERT_TRUE_FATAL(
      read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;

  {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    p = Y(plan_ng_guru)(1, &N, 0, &n, M, m, Y(get_window_id)(), x, f_hat,
                        got, 0u,
                        NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    {
      FILE *tf = tmpfile();
      CU_ASSERT_PTR_NOT_NULL_FATAL(tf);
      Y(fprint_plan)
      (p, tf);
      {
        long len;
        char *buf;
        rewind(tf);
        fseek(tf, 0, SEEK_END);
        len = ftell(tf);
        rewind(tf);
        buf = (char *)Y(malloc)((size_t)len + 1);
        if (len > 0)
          CU_ASSERT_EQUAL(fread(buf, 1, (size_t)len, tf), (size_t)len);
        buf[len] = '\0';
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "nfft_solver_fast_native"));
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "deconv"));
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "conv"));
        Y(free)
        (buf);
      }
      fclose(tf);
    }

    Y(plan_ng_destroy)
    (p);
    Y(free)
    (got);
  }

  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* A window other than the compile-time one makes every NFFT-kind solver
 * decline -- the native fast on its own window gate, the direct natives via
 * NFFT_NO_DIRECT -- so no candidate survives and the guru returns NULL. */
void Y(check_nfast_native_declines_window)(void) {
  INT N, M;
  int m = 7;
  R *x;
  C *f_hat, *f;
  INT n;
  Y(plan_ng) * p;
  int other;

  CU_ASSERT_TRUE_FATAL(
      read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;

  other = (Y(get_window_id)() == NFFT_WINDOW_KAISER_BESSEL)
              ? NFFT_WINDOW_GAUSSIAN
              : NFFT_WINDOW_KAISER_BESSEL;
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, m, other, x, f_hat, f, 0u,
                      NFFT_ESTIMATE | NFFT_NO_DIRECT | NFFT_NO_FAST_NATIVE);
  CU_ASSERT_PTR_NULL(p);

  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* NFFT_NO_FAST_NATIVE is selective: with the direct natives still available
 * one of them wins; excluding those too leaves nothing. */
void Y(check_nfast_flag_selective)(void) {
  INT N, M;
  int m = 7;
  R *x;
  C *f_hat, *f;
  INT n;
  Y(plan_ng) * p;

  CU_ASSERT_TRUE_FATAL(
      read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;

  /* the fast excluded, direct still available */
  {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    p = Y(plan_ng_guru)(1, &N, 0, &n, M, m, Y(get_window_id)(), x, f_hat,
                        got, 0u,
                        NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);
    Y(plan_ng_destroy)
    (p);
    Y(free)
    (got);
  }

  /* everything excluded */
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, m, Y(get_window_id)(), x, f_hat, f,
                      0u,
                      NFFT_ESTIMATE | NFFT_NO_DIRECT | NFFT_NO_FAST_NATIVE);
  CU_ASSERT_PTR_NULL(p);

  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* The 2D reference cases, the same geometries tests/nfft.c's 2D suite runs,
 * plus the _t210 type-II variants. */
typedef struct
{
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
void Y(check_nfast_native_fast_2d)(void) {
  size_t fi;
  for (fi = 0; fi < NFAST_2D_NFILES; fi++) {
    int d, m = 7;
    INT *N, M;
    R *x;
    C *f_hat, *f;
    INT n[2];
    Y(plan_ng) * p;
    R err;

    CU_ASSERT_TRUE_FATAL(
        read_nd_case(files_2d[fi].file, &d, &N, &M, &x, &f_hat, &f));
    CU_ASSERT_EQUAL_FATAL(d, 2);
    n[0] = 2 * N[0];
    n[1] = 2 * N[1];

    {
      C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
      p = Y(plan_ng_guru)(2, N, files_2d[fi].variant, n, M, m,
                          Y(get_window_id)(), x, f_hat, got, 0u,
                          NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL_FATAL(p);

      Y(precompute)
      (p);
      Y(execute)
      (p);
      err = rel_max_errC(got, f, M);
      CU_ASSERT(err < (R)1e-5);

      {
        FILE *tf = tmpfile();
        CU_ASSERT_PTR_NOT_NULL_FATAL(tf);
        Y(fprint_plan)
        (p, tf);
        {
          long len;
          char *buf;
          rewind(tf);
          fseek(tf, 0, SEEK_END);
          len = ftell(tf);
          rewind(tf);
          buf = (char *)Y(malloc)((size_t)len + 1);
          if (len > 0)
            CU_ASSERT_EQUAL(fread(buf, 1, (size_t)len, tf), (size_t)len);
          buf[len] = '\0';
          CU_ASSERT_PTR_NOT_NULL(strstr(buf, "nfft_solver_fast_native"));
          CU_ASSERT_PTR_NOT_NULL(strstr(buf, "deconv"));
          CU_ASSERT_PTR_NOT_NULL(strstr(buf, "conv"));
          Y(free)
          (buf);
        }
        fclose(tf);
      }

      /* in-test legacy reference: the same fast algorithm through the old
       * X(plan) API, which has no type-II concept -- type-I cases only. */
      if (!files_2d[fi].variant) {
        NFFT(plan)
        lp;
        int Ni[2], ni[2];
        INT j;
        Ni[0] = (int)N[0];
        Ni[1] = (int)N[1];
        ni[0] = (int)n[0];
        ni[1] = (int)n[1];
        NFFT(init_guru)
        (&lp, 2, Ni, (int)M, ni, m,
         PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
             FFTW_INIT | FFT_OUT_OF_PLACE,
         FFTW_ESTIMATE);
        for (j = 0; j < M * 2; j++)
          lp.x[j] = x[j];
        for (j = 0; j < N[0] * N[1]; j++)
          lp.f_hat[j] = f_hat[j];
        NFFT(precompute_one_psi)
        (&lp);
        NFFT(trafo)
        (&lp);
        {
          R errl = rel_max_errC(got, lp.f, M);
          CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
        }
        NFFT(finalize)
        (&lp);
      }

      Y(plan_ng_destroy)
      (p);
      Y(free)
      (got);
    }

    Y(free)
    (N);
    Y(free)
    (x);
    Y(free)
    (f_hat);
    Y(free)
    (f);
    Y(the_planner_destroy)
    ();
  }
}

/* The d=2 adjoint slice: same geometries, reading the separately generated
 * nfft_adjoint_2d_*.txt references rather than round-tripping the forward
 * case, and checked against an in-test legacy X(adjoint_2d) as well. */
void Y(check_nfast_native_fast_2d_adjoint)(void) {
  size_t fi;
  for (fi = 0; fi < NFAST_2D_NFILES; fi++) {
    int d, m = 7;
    INT *N, M;
    R *x;
    C *f_hat_ref, *f;
    INT n[2];
    Y(plan_ng) * p;
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

      Y(precompute)
      (p);
      Y(execute_adjoint)
      (p);
      err = rel_max_errC(got_fhat, f_hat_ref, N[0] * N[1]);
      CU_ASSERT(err < (R)1e-5);

      if (!adjoint_files_2d[fi].variant) {
        NFFT(plan)
        lp;
        int Ni[2], ni[2];
        INT j;
        Ni[0] = (int)N[0];
        Ni[1] = (int)N[1];
        ni[0] = (int)n[0];
        ni[1] = (int)n[1];
        NFFT(init_guru)
        (&lp, 2, Ni, (int)M, ni, m,
         PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
             FFTW_INIT | FFT_OUT_OF_PLACE,
         FFTW_ESTIMATE);
        for (j = 0; j < M * 2; j++)
          lp.x[j] = x[j];
        for (j = 0; j < M; j++)
          lp.f[j] = f[j];
        NFFT(precompute_one_psi)
        (&lp);
        NFFT(adjoint)
        (&lp);
        {
          R errl = rel_max_errC(got_fhat, lp.f_hat, N[0] * N[1]);
          CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
        }
        NFFT(finalize)
        (&lp);
      }

      Y(plan_ng_destroy)
      (p);
      Y(free)
      (got_fhat);
    }

    Y(free)
    (N);
    Y(free)
    (x);
    Y(free)
    (f_hat_ref);
    Y(free)
    (f);
    Y(the_planner_destroy)
    ();
  }
}

/* The d=3 forward slice, single reference geometry N={10,10,10}, M=10:
 * planned under NFFT_NO_DIRECT so the composed native fast is the sole
 * survivor, checked against the file to 1e-5 and against an in-test legacy
 * X(trafo_3d) to NFAST_LEGACY_REL_BOUND. */
void Y(check_nfast_native_fast_3d)(void) {
  int d, m = 7;
  INT *N, M;
  R *x;
  C *f_hat, *f;
  INT n[3];
  Y(plan_ng) * p;
  R err;

  CU_ASSERT_TRUE_FATAL(
      read_nd_case("data/nfft_3d_10_10_10_10.txt", &d, &N, &M, &x, &f_hat, &f));
  CU_ASSERT_EQUAL_FATAL(d, 3);
  n[0] = 2 * N[0];
  n[1] = 2 * N[1];
  n[2] = 2 * N[2];

  {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    p = Y(plan_ng_guru)(3, N, 0, n, M, m, Y(get_window_id)(), x, f_hat,
                        got, 0u,
                        NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    Y(precompute)
    (p);
    Y(execute)
    (p);
    err = rel_max_errC(got, f, M);
    /* The double-generated 3D reference is float-ill-conditioned: heavy
     * cancellation in the output, so no fast NFFT (legacy included) reproduces
     * it to 1e-5 in float. Float correctness rests on the native-vs-legacy
     * agreement below and on check_nfast_float_accuracy. */
#ifndef NFFT_SINGLE
    CU_ASSERT(err < (R)1e-5);
#else
    (void)err;
#endif

    {
      FILE *tf = tmpfile();
      CU_ASSERT_PTR_NOT_NULL_FATAL(tf);
      Y(fprint_plan)
      (p, tf);
      {
        long len;
        char *buf;
        rewind(tf);
        fseek(tf, 0, SEEK_END);
        len = ftell(tf);
        rewind(tf);
        buf = (char *)Y(malloc)((size_t)len + 1);
        if (len > 0)
          CU_ASSERT_EQUAL(fread(buf, 1, (size_t)len, tf), (size_t)len);
        buf[len] = '\0';
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "nfft_solver_fast_native"));
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "deconv"));
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "conv"));
        Y(free)
        (buf);
      }
      fclose(tf);
    }

    /* in-test legacy reference: the same fast algorithm through the old
     * X(plan) API. */
    {
      NFFT(plan)
      lp;
      int Ni[3], ni[3];
      INT j;
      Ni[0] = (int)N[0];
      Ni[1] = (int)N[1];
      Ni[2] = (int)N[2];
      ni[0] = (int)n[0];
      ni[1] = (int)n[1];
      ni[2] = (int)n[2];
      NFFT(init_guru)
      (&lp, 3, Ni, (int)M, ni, m,
       PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
           FFTW_INIT | FFT_OUT_OF_PLACE,
       FFTW_ESTIMATE);
      for (j = 0; j < M * 3; j++)
        lp.x[j] = x[j];
      for (j = 0; j < N[0] * N[1] * N[2]; j++)
        lp.f_hat[j] = f_hat[j];
      NFFT(precompute_one_psi)
      (&lp);
      NFFT(trafo)
      (&lp);
      {
        R errl = rel_max_errC(got, lp.f, M);
#ifndef NFFT_SINGLE
        CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
#else
        /* legacy collapses in float 3D and diverges from the native */
        (void)errl;
#endif
      }
      NFFT(finalize)
      (&lp);
    }

    Y(plan_ng_destroy)
    (p);
    Y(free)
    (got);
  }

  Y(free)
  (N);
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* The d=3 adjoint slice: same geometry, reading the separately generated
 * nfft_adjoint_3d_10_10_10_10.txt reference rather than round-tripping the
 * forward case, and checked against an in-test legacy X(adjoint_3d) as
 * well. */
void Y(check_nfast_native_fast_3d_adjoint)(void) {
  int d, m = 7;
  INT *N, M;
  R *x;
  C *f_hat_ref, *f;
  INT n[3];
  Y(plan_ng) * p;
  R err;

  CU_ASSERT_TRUE_FATAL(read_nd_case("data/nfft_adjoint_3d_10_10_10_10.txt",
                                    &d, &N, &M, &x, &f_hat_ref, &f));
  CU_ASSERT_EQUAL_FATAL(d, 3);
  n[0] = 2 * N[0];
  n[1] = 2 * N[1];
  n[2] = 2 * N[2];

  {
    C *got_fhat = (C *)Y(malloc)((size_t)(N[0] * N[1] * N[2]) * sizeof(C));
    p = Y(plan_ng_guru)(3, N, 0, n, M, m, Y(get_window_id)(), x,
                        got_fhat, f, 0u,
                        NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    Y(precompute)
    (p);
    Y(execute_adjoint)
    (p);
    err = rel_max_errC(got_fhat, f_hat_ref, N[0] * N[1] * N[2]);
    /* Same float-conditioning caveat as the 3D forward case. */
#ifndef NFFT_SINGLE
    CU_ASSERT(err < (R)1e-5);
#else
    (void)err;
#endif

    {
      NFFT(plan)
      lp;
      int Ni[3], ni[3];
      INT j;
      Ni[0] = (int)N[0];
      Ni[1] = (int)N[1];
      Ni[2] = (int)N[2];
      ni[0] = (int)n[0];
      ni[1] = (int)n[1];
      ni[2] = (int)n[2];
      NFFT(init_guru)
      (&lp, 3, Ni, (int)M, ni, m,
       PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
           FFTW_INIT | FFT_OUT_OF_PLACE,
       FFTW_ESTIMATE);
      for (j = 0; j < M * 3; j++)
        lp.x[j] = x[j];
      for (j = 0; j < M; j++)
        lp.f[j] = f[j];
      NFFT(precompute_one_psi)
      (&lp);
      NFFT(adjoint)
      (&lp);
      {
        R errl = rel_max_errC(got_fhat, lp.f_hat, N[0] * N[1] * N[2]);
#ifndef NFFT_SINGLE
        CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
#else
        /* legacy collapses in float 3D and diverges from the native */
        (void)errl;
#endif
      }
      NFFT(finalize)
      (&lp);
    }

    Y(plan_ng_destroy)
    (p);
    Y(free)
    (got_fhat);
  }

  Y(free)
  (N);
  Y(free)
  (x);
  Y(free)
  (f_hat_ref);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* The generic rnk >= 4 slice, served by the deconv-nd/conv-nd carry-odometer
 * leaves. tests/refgen generates no d >= 4 reference file, so the case is
 * built in-test from a deterministic sequence. Two oracles: an in-test legacy
 * X(plan), whose d==4 hand-unrolled branch is mathematically the generic
 * odometer, to NFAST_LEGACY_REL_BOUND; and the direct NDFT native, forced via
 * NFFT_NO_FAST_NATIVE, to the coarser 1e-5 fast-NFFT bound. */
static R seq_4d(INT k) {
  /* Deterministic stand-in for a PRNG, so results do not depend on the
   * platform's rand. */
  R v = SIN((R)(k + 1) * K(12.9898)) * K(43758.5453);
  v -= FLOOR(v); /* fractional part, in [0,1) */
  return v;
}

void Y(check_nfast_native_fast_4d)(void) {
  const int d = 4;
  INT N[4] = {10, 10, 10, 10};
  INT n[4];
  const INT M = 40;
  int m = 7, t;
  INT Ntot, j;
  R *x;
  C *f_hat, *got;
  Y(plan_ng) * p;
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
    f_hat[j] = (seq_4d(1000 + 2 * j) - K(0.5)) +
               II * (seq_4d(1000 + 2 * j + 1) - K(0.5));

  p = Y(plan_ng_guru)(d, N, 0, n, M, m, Y(get_window_id)(), x, f_hat,
                      got, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);

  Y(precompute)
  (p);
  Y(execute)
  (p);

  {
    FILE *tf = tmpfile();
    CU_ASSERT_PTR_NOT_NULL_FATAL(tf);
    Y(fprint_plan)
    (p, tf);
    {
      long len;
      char *buf;
      rewind(tf);
      fseek(tf, 0, SEEK_END);
      len = ftell(tf);
      rewind(tf);
      buf = (char *)Y(malloc)((size_t)len + 1);
      if (len > 0)
        CU_ASSERT_EQUAL(fread(buf, 1, (size_t)len, tf), (size_t)len);
      buf[len] = '\0';
      CU_ASSERT_PTR_NOT_NULL(strstr(buf, "nfft_solver_fast_native"));
      CU_ASSERT_PTR_NOT_NULL(strstr(buf, "deconv"));
      CU_ASSERT_PTR_NOT_NULL(strstr(buf, "conv"));
      Y(free)
      (buf);
    }
    fclose(tf);
  }

  /* in-test legacy reference: the same fast algorithm through the old
   * X(plan) API. */
  {
    NFFT(plan)
    lp;
    int Ni[4], ni[4];
    INT jj;
    for (t = 0; t < d; t++) {
      Ni[t] = (int)N[t];
      ni[t] = (int)n[t];
    }
    NFFT(init_guru)
    (&lp, d, Ni, (int)M, ni, m,
     PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
         FFTW_INIT | FFT_OUT_OF_PLACE,
     FFTW_ESTIMATE);
    for (jj = 0; jj < M * d; jj++)
      lp.x[jj] = x[jj];
    for (jj = 0; jj < Ntot; jj++)
      lp.f_hat[jj] = f_hat[jj];
    NFFT(precompute_one_psi)
    (&lp);
    NFFT(trafo)
    (&lp);
    {
      R errl = rel_max_errC(got, lp.f, M);
      CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
    }
    NFFT(finalize)
    (&lp);
  }

  /* second oracle: the direct NDFT native, forced via NFFT_NO_FAST_NATIVE. */
  {
    C *got_direct = (C *)Y(malloc)((size_t)M * sizeof(C));
    Y(plan_ng) *pdirect = Y(plan_ng_guru)(
        d, N, 0, n, M, m, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, got_direct, 0u,
        NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pdirect);
    Y(precompute)
    (pdirect);
    Y(execute)
    (pdirect);
    err = rel_max_errC(got, got_direct, M);
    CU_ASSERT(err < (R)1e-5);
    Y(plan_ng_destroy)
    (pdirect);
    Y(free)
    (got_direct);
  }

  Y(plan_ng_destroy)
  (p);
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (got);
  Y(the_planner_destroy)
  ();
}

/* The d=4 adjoint slice: same in-test-generated geometry, oracle is
 * NFFT(adjoint) on a separate legacy plan. f is generated independently, so
 * this is not a round-trip of the forward half. */
void Y(check_nfast_native_fast_4d_adjoint)(void) {
  const int d = 4;
  INT N[4] = {10, 10, 10, 10};
  INT n[4];
  const INT M = 40;
  int m = 7, t;
  INT Ntot, j;
  R *x;
  C *f, *got_fhat;
  Y(plan_ng) * p;

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
    f[j] = (seq_4d(9000 + 2 * j) - K(0.5)) +
           II * (seq_4d(9000 + 2 * j + 1) - K(0.5));

  p = Y(plan_ng_guru)(d, N, 0, n, M, m, Y(get_window_id)(), x, got_fhat,
                      f, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);

  Y(precompute)
  (p);
  Y(execute_adjoint)
  (p);

  {
    NFFT(plan)
    lp;
    int Ni[4], ni[4];
    INT jj;
    for (t = 0; t < d; t++) {
      Ni[t] = (int)N[t];
      ni[t] = (int)n[t];
    }
    NFFT(init_guru)
    (&lp, d, Ni, (int)M, ni, m,
     PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
         FFTW_INIT | FFT_OUT_OF_PLACE,
     FFTW_ESTIMATE);
    for (jj = 0; jj < M * d; jj++)
      lp.x[jj] = x[jj];
    for (jj = 0; jj < M; jj++)
      lp.f[jj] = f[jj];
    NFFT(precompute_one_psi)
    (&lp);
    NFFT(adjoint)
    (&lp);
    {
      R errl = rel_max_errC(got_fhat, lp.f_hat, Ntot);
      CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
    }
    NFFT(finalize)
    (&lp);
  }

  Y(plan_ng_destroy)
  (p);
  Y(free)
  (x);
  Y(free)
  (f);
  Y(free)
  (got_fhat);
  Y(the_planner_destroy)
  ();
}

/* Robust log(I0) and scaled I0 helpers. */
void Y(check_nfast_bessel_log_scaled)(void) {
  int i;
  for (i = 0; i < 20; i++) /* small branch: log I0 == LOG(bessel_i0) */
  {
    R x = K(0.5) * (R)i;
    CU_ASSERT(FABS(Y(bessel_i0_log)(x) - LOG(Y(bessel_i0)(x))) <= K(1e-6) * (K(1.0) + FABS(LOG(Y(bessel_i0)(x)))));
  }
  { /* asymptotic: log I0(x) ~ x - 0.5 log x - 0.5 log(2 pi) */
    R x = K(35.0);
    CU_ASSERT(FABS(Y(bessel_i0_log)(x) - (x - K(0.5) * LOG(x) - K(0.5) * LOG(K(2.0) * KPI))) <= K(0.05));
  }
  { /* scaled(x, log I0(x)) == 1 at the peak; O(<=1) below it; no overflow. */
    R x = K(33.0), lg = Y(bessel_i0_log)(x);
    CU_ASSERT(FABS(Y(bessel_i0_scaled)(x, lg) - K(1.0)) <= K(1e-5));
    {
      R s = Y(bessel_i0_scaled)(K(20.0), lg);
      CU_ASSERT(s > K(0.0) && s <= K(1.0));
    }
  }
  {
    R big = K(200.0), lg = Y(bessel_i0_log)(big), s = Y(bessel_i0_scaled)(big, lg);
    CU_ASSERT(s == s && s > K(0.0) && s < Y(float_property)(NFFT_R_MAX));
  } /* finite */
}

/* The single-precision regression gate for the self-normalizing KB window.
 * check_nfast_native_fast_3d(_adjoint) gate their native-vs-file check to
 * double; this runs the same comparison on the same reference files in
 * whatever precision is compiled, against the precision-aware err_bound
 * instead of a fixed 1e-5. The error was ~1.0 in float before the window
 * self-normalized, ~2-3e-6 after. The plans are pinned to KAISER_BESSEL with
 * n = 2*N, so the bound is queried at that window and sigma = 2. */
void Y(check_nfast_float_accuracy)(void) {
  int m = 7;
  R bound = err_bound(NFFT_WINDOW_KAISER_BESSEL, (R)m, K(2.0));

  /* forward: nfft_3d_10_10_10_10 (true/direct transform, double-generated). */
  {
    int d;
    INT *N, M;
    R *x;
    C *f_hat, *f;
    INT n[3];
    Y(plan_ng) * p;
    R err;

    CU_ASSERT_TRUE_FATAL(read_nd_case("data/nfft_3d_10_10_10_10.txt", &d,
                                      &N, &M, &x, &f_hat, &f));
    CU_ASSERT_EQUAL_FATAL(d, 3);
    n[0] = 2 * N[0];
    n[1] = 2 * N[1];
    n[2] = 2 * N[2];

    {
      C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
      p = Y(plan_ng_guru)(3, N, 0, n, M, m, NFFT_WINDOW_KAISER_BESSEL, x, f_hat,
                          got, 0u,
                          NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL_FATAL(p);

      Y(precompute)
      (p);
      Y(execute)
      (p);
      err = rel_max_errC(got, f, M);
      CU_ASSERT(err < bound);

      Y(plan_ng_destroy)
      (p);
      Y(free)
      (got);
    }

    Y(free)
    (N);
    Y(free)
    (x);
    Y(free)
    (f_hat);
    Y(free)
    (f);
    Y(the_planner_destroy)
    ();
  }

  /* adjoint: nfft_adjoint_3d_10_10_10_10, a separately generated reference,
   * not a round-trip of the forward case above. */
  {
    int d;
    INT *N, M;
    R *x;
    C *f_hat_ref, *f;
    INT n[3];
    Y(plan_ng) * p;
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
                          got_fhat, f, 0u,
                          NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL_FATAL(p);

      Y(precompute)
      (p);
      Y(execute_adjoint)
      (p);
      err = rel_max_errC(got_fhat, f_hat_ref, N[0] * N[1] * N[2]);
      CU_ASSERT(err < bound);

      Y(plan_ng_destroy)
      (p);
      Y(free)
      (got_fhat);
    }

    Y(free)
    (N);
    Y(free)
    (x);
    Y(free)
    (f_hat_ref);
    Y(free)
    (f);
    Y(the_planner_destroy)
    ();
  }
}

/* Near the peak, where cancellation bites hardest, the rationalized KB
 * evaluator beats a naive exp(a - lg_peak) fold. Checked against a
 * log-domain reference in the working precision. */
void Y(check_nfast_window_cancellation)(void) {
  const int w = NFFT_WINDOW_KAISER_BESSEL;
  INT n = 512, N = 256;
  int m = 8, k;
  R b = KPI * (K(2.0) - (R)N / (R)n), xpk = (R)m * b, lg = Y(bessel_i0_log)(xpk);
  for (k = 0; k <= 3; k++) /* the near-peak bins where a ~ m*b */
  {
    R t = K(2.0) * KPI * (R)k / (R)n, a = (R)m * SQRT(b * b - t * t);
    R got = Y(window_phi_hut)(w, n, N, m, (INT)k);
    R ref = EXP(Y(bessel_i0_log)(a) - lg); /* stable log-domain reference */
    CU_ASSERT(FABS(got - ref) <= K(20.0) * Y(float_property)(NFFT_EPSILON) * (K(1.0) + FABS(ref)));
  }
  /* logtail identity: log I0(x) == x + logtail(x). */
  {
    R x = K(40.0);
    CU_ASSERT(FABS(Y(bessel_i0_log)(x) - (x + Y(bessel_i0_logtail)(x))) <= K(1e-5));
  }
}

void Y(check_nfast_window_apply)(void) {
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
    Y(window_phi_hut_apply)
    (NFFT_WINDOW_KAISER_BESSEL, n, N, m, -N / 2, buf, N);
    CU_ASSERT(FABS(buf[N / 2] - K(1.0)) <= K(1e-5));
    for (k = 0; k < N; k++) {
      R got = buf[k], ref =
                          Y(window_phi_hut)(NFFT_WINDOW_KAISER_BESSEL, n, N, m, k - N / 2);
      CU_ASSERT(FABS(got - ref) <= NFAST_LEGACY_REL_BOUND * FMAX(FABS(got), FABS(ref)));
      CU_ASSERT(buf[k] > K(0.0) && buf[k] <= K(1.0) + K(1e-5));
    }
    CU_ASSERT(buf[N / 2] > buf[N / 2 + 1]); /* monotone decay */

    /* The per-element fallback serves the raw windows too. */
    {
      R *gbuf = (R *)Y(malloc)((size_t)N * sizeof(R));
      INT kk;
      Y(window_phi_hut_apply)
      (NFFT_WINDOW_GAUSSIAN, n, N, m, -N / 2, gbuf, N);
      for (kk = 0; kk < N; kk++) {
        R single = Y(window_phi_hut)(NFFT_WINDOW_GAUSSIAN, n, N, m, kk - N / 2);
        CU_ASSERT(FABS(gbuf[kk] - single) <= NFAST_LEGACY_REL_BOUND * (K(1.0) + FABS(single)));
      }
      Y(free)
      (gbuf);
    }

    /* phi_precompute: per-node psi matches single-arg window_phi at every
     * tap, under the same 1-ULP caveat. */
    {
      INT M = 4, j, lj;
      R x[4] = {K(0.1), K(0.37), K(0.5), K(0.83)};
      R *psi = (R *)Y(malloc)((size_t)M * (size_t)(2 * m + 2) * sizeof(R));
      Y(window_phi_precompute)
      (NFFT_WINDOW_KAISER_BESSEL, n, N, m, x, 1, M, psi,
       2 * m + 2);
      for (j = 0; j < M; j++) {
        INT c = LRINT(FLOOR((R)n * x[j]));
        for (lj = 0; lj <= 2 * m + 1; lj++) {
          INT idx = c - m + lj;
          R got = psi[j * (2 * m + 2) + lj], ref =
                                                 Y(window_phi)(NFFT_WINDOW_KAISER_BESSEL, n, N, m,
                                                               x[j] - (R)idx / (R)n);
          CU_ASSERT(FABS(got - ref) <= NFAST_LEGACY_REL_BOUND * FMAX(FABS(got), FABS(ref)));
        }
      }
      Y(free)
      (psi);
    }
    Y(free)
    (buf);
  }
}

/* The composed native fast and its DECONV/CONV children accept KB, Gaussian,
 * B-spline and sinc at runtime; Dirac and out-of-range ordinals are declined
 * with a NULL plan. */
void Y(check_nfast_native_window_select)(void) {
  INT N, M, n;
  R *x;
  C *f_hat, *f;
  int w;
  CU_ASSERT_TRUE_FATAL(
      read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;
  for (w = NFFT_WINDOW_KAISER_BESSEL; w <= NFFT_WINDOW_SINC_POWER; w++) {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    Y(plan_ng) *p = Y(plan_ng_guru)(
        1, &N, 0, &n, M, 7, w, x, f_hat, got, 0u,
        NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p) {
      Y(plan_ng_destroy)
      (p);
    }
    Y(free)
    (got);
  }
  /* Dirac and a garbage ordinal are both declined. */
  {
    int bad[2];
    int i;
    bad[0] = NFFT_WINDOW_DIRAC_DELTA;
    bad[1] = 99;
    for (i = 0; i < 2; i++) {
      Y(plan_ng) *pbad = Y(plan_ng_guru)(
          1, &N, 0, &n, M, 7, bad[i], x, f_hat, f, 0u,
          NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NULL(pbad);
    }
  }
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* Each runtime-selected window's native fast NFFT against the direct NDFT
 * native, which is window-independent, so one build covers all four windows.
 * Selected via NFAST_WINDOWS. */
void Y(check_nfast_window_accuracy)(void) {
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
    R bound = err_bound(w, (R)m, sigma);
    /* forward */
    {
      C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
      C *ref = (C *)Y(malloc)((size_t)M * sizeof(C));
      Y(plan_ng) *pn = Y(plan_ng_guru)(
          1, &N, 0, &n, M, m, w, x, f_hat, got, 0u,
          NFFT_ESTIMATE | NFFT_NO_DIRECT);
      Y(plan_ng) *pd = Y(plan_ng_guru)(
          1, &N, 0, &n, M, m, w, x, f_hat, ref, 0u,
          NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pn);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pd);
      Y(precompute)
      (pn);
      Y(execute)
      (pn);
      Y(precompute)
      (pd);
      Y(execute)
      (pd);
      CU_ASSERT(rel_max_errC(got, ref, M) < bound);
      Y(plan_ng_destroy)
      (pn);
      Y(plan_ng_destroy)
      (pd);
      Y(free)
      (got);
      Y(free)
      (ref);
    }
    /* adjoint */
    {
      C *got = (C *)Y(malloc)((size_t)N * sizeof(C));
      C *ref = (C *)Y(malloc)((size_t)N * sizeof(C));
      Y(plan_ng) *pn = Y(plan_ng_guru)(
          1, &N, 0, &n, M, m, w, x, got, f, 0u,
          NFFT_ESTIMATE | NFFT_NO_DIRECT);
      Y(plan_ng) *pd = Y(plan_ng_guru)(
          1, &N, 0, &n, M, m, w, x, ref, f, 0u,
          NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pn);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pd);
      Y(precompute)
      (pn);
      Y(execute_adjoint)
      (pn);
      Y(precompute)
      (pd);
      Y(execute_adjoint)
      (pd);
      CU_ASSERT(rel_max_errC(got, ref, N) < bound);
      Y(plan_ng_destroy)
      (pn);
      Y(plan_ng_destroy)
      (pd);
      Y(free)
      (got);
      Y(free)
      (ref);
    }
  }
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}
