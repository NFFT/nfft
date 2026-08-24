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

/* Unit tests for the runtime window vtable (kernel/util/window.c): the
 * evaluator must reproduce the compile-time PHI/PHI_HUT macros for every
 * implemented window, stay numerically stable in the log domain, and be
 * discriminated by the wisdom key.  Nothing here runs a transform. */

#include <complex.h> /* before nfft3.h so fftw_complex is C-compatible */
#include <stdio.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "config.h" /* ABS_SRCDIR */
#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "window.h"

/* Two evaluator paths for the same window value reassociate differently under
 * -ffast-math and differ by ~1 ULP, so they are compared relatively. */
#if defined(NFFT_SINGLE)
#define WINDOW_REL_TOL ((R)1e-4)
#else
#define WINDOW_REL_TOL ((R)1e-12)
#endif

void Y(check_window_id)(void) {
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

/* The KB vtable is peak-normalized, so these assertions check only positivity
 * and monotone decay, which a uniform positive scale preserves. Exact
 * normalization is checked in check_window_normalized. */
void Y(check_window_vtable)(void) {
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
void Y(check_window_normalized)(void) {
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
void Y(check_window_all)(void) {
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
void Y(check_window_key)(void) {
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

/* Robust log(I0) and scaled I0 helpers. */
void Y(check_window_bessel_log_scaled)(void) {
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

/* Near the peak, where cancellation bites hardest, the rationalized KB
 * evaluator beats a naive exp(a - lg_peak) fold. Checked against a
 * log-domain reference in the working precision. */
void Y(check_window_cancellation)(void) {
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

void Y(check_window_apply)(void) {
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
      CU_ASSERT(FABS(got - ref) <= WINDOW_REL_TOL * FMAX(FABS(got), FABS(ref)));
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
        CU_ASSERT(FABS(gbuf[kk] - single) <= WINDOW_REL_TOL * (K(1.0) + FABS(single)));
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
          CU_ASSERT(FABS(got - ref) <= WINDOW_REL_TOL * FMAX(FABS(got), FABS(ref)));
        }
      }
      Y(free)
      (psi);
    }
    Y(free)
    (buf);
  }
}
