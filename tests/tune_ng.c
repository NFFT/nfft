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

#include <CUnit/CUnit.h>

#include "infft.h"
#include "nfft3.h"
#include "tune_ng.h"

#include <string.h>

#define TUNE_FWD 0
#define TUNE_ADJ 1

/* Y(tune) is stated for the Kaiser-Bessel window, so the plans built here ask
 * for it by ordinal rather than taking the compile-time window: every window
 * build cell then exercises the same claim. */
#define TUNE_WINDOW NFFT_WINDOW_KAISER_BESSEL

/* Measure the error of one Kaiser-Bessel plan at cut-off m, in the measure
 * Y(tune) predicts: the max node error over the l_1 norm of the input. The
 * oracle is the planner's own direct NDFT on the same nodes and data. */
static R measure(INT N, INT n, INT M, int m, int adjoint)
{
  R *x;
  C *in, *out, *ref;
  Y(plan_ng) * fast;
  Y(plan_ng) * direct;
  const INT in_len = adjoint ? M : N;
  const INT out_len = adjoint ? N : M;
  R num = K(0.0), den = K(0.0);
  INT j;

  x = (R *)Y(malloc)((size_t)M * sizeof(R));
  in = (C *)Y(malloc)((size_t)in_len * sizeof(C));
  out = (C *)Y(malloc)((size_t)out_len * sizeof(C));
  ref = (C *)Y(malloc)((size_t)out_len * sizeof(C));

  Y(srand48)(1723);
  Y(vrand_shifted_unit_double)(x, M);
  Y(vrand_unit_complex)(in, in_len);

  /* Two scratch arrays per plan: the guru binds f_hat and f, and the _on
   * variants then run on the arrays passed here. */
  {
    C *fh = (C *)Y(malloc)((size_t)N * sizeof(C));
    C *ff = (C *)Y(malloc)((size_t)M * sizeof(C));
    memset(fh, 0, (size_t)N * sizeof(C));
    memset(ff, 0, (size_t)M * sizeof(C));

    fast = Y(plan_ng_guru)(1, &N, 0, &n, M, m, TUNE_WINDOW, x, fh, ff, 0u,
                           NFFT_ESTIMATE | NFFT_NO_DIRECT);
    direct = Y(plan_ng_guru)(1, &N, 0, &n, M, m, TUNE_WINDOW, x, fh, ff, 0u,
                             NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
    if (!fast || !direct) {
      CU_FAIL("plan_ng_guru declined a geometry tune reported as usable");
      Y(plan_ng_destroy)(fast);
      Y(plan_ng_destroy)(direct);
      Y(free)(fh);
      Y(free)(ff);
      Y(free)(x);
      Y(free)(in);
      Y(free)(out);
      Y(free)(ref);
      return K(-1.0);
    }
    Y(precompute)(fast);
    Y(precompute)(direct);

    if (adjoint) {
      Y(execute_adjoint_on)(fast, out, in);
      Y(execute_adjoint_on)(direct, ref, in);
    } else {
      Y(execute_on)(fast, in, out);
      Y(execute_on)(direct, in, ref);
    }

    Y(plan_ng_destroy)(fast);
    Y(plan_ng_destroy)(direct);
    Y(free)(fh);
    Y(free)(ff);
  }

  for (j = 0; j < out_len; j++) {
    R e = CABS(out[j] - ref[j]);
    if (e > num)
      num = e;
  }
  for (j = 0; j < in_len; j++)
    den += CABS(in[j]);

  Y(free)(x);
  Y(free)(in);
  Y(free)(out);
  Y(free)(ref);

  return den == K(0.0) ? K(0.0) : num / den;
}

void Y(check_tune_meets_goal)(void)
{
  static const INT Ns[] = {64, 128, 256};
  static const R sigmas[] = {K(1.25), K(1.5), K(2.0), K(3.0)};
  const R eps = Y(float_property)(NFFT_EPSILON);
  /* Stay well above the direct NDFT's own roundoff, so the oracle is clean. */
  const R goal_floor = K(1.0e4) * eps;
  size_t iN, is, ig, dir;
  static const R goals[] = {K(1.0e-2), K(1.0e-4),  K(1.0e-6),
                            K(1.0e-8), K(1.0e-10), K(1.0e-12)};

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (is = 0; is < sizeof(sigmas) / sizeof(sigmas[0]); is++) {
      const INT N = Ns[iN];
      INT n = (INT)(sigmas[is] * (R)N + K(0.5));
      if (n % 2)
        n++;
      if (n <= N)
        n = N + 2;

      for (ig = 0; ig < sizeof(goals) / sizeof(goals[0]); ig++) {
        const R goal = goals[ig];
        if (goal < goal_floor)
          continue;
        for (dir = 0; dir <= 1; dir++) {
          int m = 0;
          R attained = K(0.0);
          const int rc = Y(tune)(N, n, (int)dir, goal, &m, &attained);
          R err;
          CU_ASSERT(rc >= 0);
          if (rc != 1)
            continue; /* goal unattainable here: covered by the other test */
          CU_ASSERT(m >= 1);
          CU_ASSERT(m <= (int)(n / 2 - 1));
          CU_ASSERT(attained <= goal);
          err = measure(N, n, 2 * N, m, (int)dir);
          CU_ASSERT(err >= K(0.0));
          if (err > goal)
            printf("\ntune(N=" __D__ ", n=" __D__
                   ", %s) -> m=%d but err=" __FE__ " > goal=" __FE__ "\n",
                   N, n, dir ? "adjoint" : "forward", m, err, goal);
          CU_ASSERT(err <= goal);
        }
      }
    }
}

void Y(check_tune_unreachable)(void)
{
  static const INT Ns[] = {64, 256};
  static const R sigmas[] = {K(1.25), K(1.5), K(2.0)};
  const R eps = Y(float_property)(NFFT_EPSILON);
  size_t iN, is, dir;

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (is = 0; is < sizeof(sigmas) / sizeof(sigmas[0]); is++) {
      const INT N = Ns[iN];
      INT n = (INT)(sigmas[is] * (R)N + K(0.5));
      if (n % 2)
        n++;
      if (n <= N)
        n = N + 2;

      for (dir = 0; dir <= 1; dir++) {
        int m = -1;
        R attained = K(-1.0);
        /* Below eps no cut-off can help, whatever the oversampling. */
        const int rc = Y(tune)(N, n, (int)dir, eps / K(1.0e3), &m, &attained);
        CU_ASSERT_EQUAL(rc, 0);
        CU_ASSERT(m >= 1);
        CU_ASSERT(m <= (int)(n / 2 - 1));
        CU_ASSERT(attained > eps / K(1.0e3));

        /* Asking for what is on offer needs no more than the argmin cut-off.
         * The goal is nudged up by an ulp or two: the two calls recompute the
         * same expression at different sites, and -ffast-math does not
         * guarantee those agree to the last bit. */
        {
          int m_lo = 0;
          R a_lo = K(0.0);
          const R loosened = attained * (K(1.0) + K(1.0e-6));
          CU_ASSERT_EQUAL(Y(tune)(N, n, (int)dir, loosened, &m_lo, &a_lo), 1);
          CU_ASSERT(m_lo <= m);
          CU_ASSERT(a_lo <= loosened);
        }
      }
    }
}

void Y(check_tune_geometries)(void)
{
  static const INT Ns[] = {32, 64, 128, 256, 512};
  static const R sigmas[] = {K(1.25), K(1.4), K(1.5), K(1.75),
                             K(2.0),  K(2.5), K(3.0), K(4.0)};
  size_t iN, is, ig, dir;
  static const R goals[] = {K(1.0e-2), K(1.0e-4), K(1.0e-6), K(1.0e-8)};

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (dir = 0; dir <= 1; dir++) {
      const INT N = Ns[iN];

      /* Walking the goals from tightest to loosest, m never grows. */
      for (is = 0; is < sizeof(sigmas) / sizeof(sigmas[0]); is++) {
        INT n = (INT)(sigmas[is] * (R)N + K(0.5));
        int prev = -1;
        if (n % 2)
          n++;
        if (n <= N)
          n = N + 2;

        for (ig = sizeof(goals) / sizeof(goals[0]); ig-- > 0;) {
          int m = 0;
          R attained = K(0.0);
          const int rc = Y(tune)(N, n, (int)dir, goals[ig], &m, &attained);
          CU_ASSERT(rc >= 0);
          CU_ASSERT(m >= 1);
          CU_ASSERT(m <= (int)(n / 2 - 1));
          CU_ASSERT(attained > K(0.0));
          if (rc == 1) {
            CU_ASSERT(attained <= goals[ig]);
            if (prev >= 0)
              CU_ASSERT(m <= prev);
            prev = m;
          }
        }
      }

      /* More oversampling never needs a larger m for the same goal. */
      for (ig = 0; ig < sizeof(goals) / sizeof(goals[0]); ig++) {
        int prev = -1;
        for (is = 0; is < sizeof(sigmas) / sizeof(sigmas[0]); is++) {
          INT n = (INT)(sigmas[is] * (R)N + K(0.5));
          int m = 0;
          R attained = K(0.0);
          int rc;
          if (n % 2)
            n++;
          if (n <= N)
            n = N + 2;
          rc = Y(tune)(N, n, (int)dir, goals[ig], &m, &attained);
          if (rc != 1)
            continue;
          if (prev >= 0)
            CU_ASSERT(m <= prev);
          prev = m;
        }
      }
    }
}

void Y(check_tune_bad_args)(void)
{
  int m = 4321;
  R attained = K(-7.0);

  /* n must exceed N. */
  CU_ASSERT_EQUAL(Y(tune)(64, 64, TUNE_FWD, K(1.0e-6), &m, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune)(64, 32, TUNE_FWD, K(1.0e-6), &m, &attained), -1);
  /* sigma below the 5/4 floor is refused, exactly at the boundary. */
  CU_ASSERT_EQUAL(Y(tune)(64, 66, TUNE_FWD, K(1.0e-6), &m, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune)(64, 79, TUNE_FWD, K(1.0e-6), &m, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune)(4, 4, TUNE_FWD, K(1.0e-6), &m, &attained), -1);
  /* Non-positive geometry. */
  CU_ASSERT_EQUAL(Y(tune)(0, 128, TUNE_FWD, K(1.0e-6), &m, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune)(-8, 128, TUNE_FWD, K(1.0e-6), &m, &attained), -1);
  /* Non-positive goal. */
  CU_ASSERT_EQUAL(Y(tune)(64, 128, TUNE_FWD, K(0.0), &m, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune)(64, 128, TUNE_FWD, K(-1.0), &m, &attained), -1);
  /* No output slot. */
  CU_ASSERT_EQUAL(Y(tune)(64, 128, TUNE_FWD, K(1.0e-6), 0, &attained), -1);
  /* Refusals leave the caller's outputs untouched. */
  CU_ASSERT_EQUAL(m, 4321);
  CU_ASSERT(attained == K(-7.0));

  /* sigma exactly 5/4 is accepted. */
  CU_ASSERT(Y(tune)(64, 80, TUNE_FWD, K(1.0e-6), &m, &attained) >= 0);

  /* attained is optional. */
  CU_ASSERT(Y(tune)(64, 128, TUNE_FWD, K(1.0e-6), &m, 0) >= 0);
  CU_ASSERT(m >= 1);
  CU_ASSERT_EQUAL(Y(tune)(64, 128, TUNE_ADJ, K(1.0e-6), &m, 0) >= 0, 1);
}

void Y(check_tune_sigma_agrees)(void)
{
  static const INT Ns[] = {32, 64, 100, 128, 256, 512};
  static const R goals[] = {K(1.0e-2),  K(1.0e-4),  K(1.0e-6), K(1.0e-8),
                            K(1.0e-10), K(1.0e-12), K(1.0e-14)};
  size_t iN, ig, dir;

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (ig = 0; ig < sizeof(goals) / sizeof(goals[0]); ig++)
      for (dir = 0; dir <= 1; dir++) {
        const INT N = Ns[iN];
        const R goal = goals[ig];
        INT n = 0;
        R attained = K(0.0);
        const int rc = Y(tune_sigma)(N, (int)dir, goal, &n, &attained);
        int m = 0;
        R m_att = K(0.0);

        CU_ASSERT(rc >= 0);
        CU_ASSERT(n % 2 == 0);
        CU_ASSERT((INT)4 * n >= (INT)5 * N);
        CU_ASSERT(n <= (INT)4 * N);

        if (rc == 0) {
          /* Out of reach even at the top of the band. */
          CU_ASSERT_EQUAL(n, (INT)4 * N);
          CU_ASSERT(attained > goal);
          continue;
        }

        CU_ASSERT(attained <= goal);
        /* The whole point: tune must agree that this n works. */
        CU_ASSERT_EQUAL(Y(tune)(N, n, (int)dir, goal, &m, &m_att), 1);
        CU_ASSERT(m >= 1);
        CU_ASSERT(m_att <= goal);

        /* And it is the smallest such n on the even grid. */
        if (n - 2 > N && (INT)4 * (n - 2) >= (INT)5 * N)
          CU_ASSERT_EQUAL(Y(tune)(N, n - 2, (int)dir, goal, &m, &m_att), 0);
      }
}

void Y(check_tune_sigma_limits)(void)
{
  const R eps = Y(float_property)(NFFT_EPSILON);
  const INT N = 128;
  INT n = 0;
  R attained = K(0.0);
  int m = 0;
  size_t dir;

  /* Bad arguments. */
  CU_ASSERT_EQUAL(Y(tune_sigma)(0, TUNE_FWD, K(1.0e-6), &n, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune_sigma)(N, TUNE_FWD, K(0.0), &n, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune_sigma)(N, TUNE_FWD, K(1.0e-6), 0, &attained), -1);

  for (dir = 0; dir <= 1; dir++) {
    /* Below eps nothing reaches the goal, at any oversampling. */
    CU_ASSERT_EQUAL(Y(tune_sigma)(N, (int)dir, eps / K(1.0e3), &n, &attained),
                    0);
    CU_ASSERT_EQUAL(n, (INT)4 * N);

    /* A loose goal is met at the bottom of the band. */
    CU_ASSERT_EQUAL(Y(tune_sigma)(N, (int)dir, K(1.0e-2), &n, &attained), 1);
    CU_ASSERT(n <= (INT)2 * N);
  }

  /* A real transform at the recommended geometry meets the goal. */
  for (dir = 0; dir <= 1; dir++) {
    const R goal = K(1.0e-6);
    R err;
    if (goal < K(1.0e4) * eps)
      continue;
    CU_ASSERT_EQUAL(Y(tune_sigma)(N, (int)dir, goal, &n, &attained), 1);
    CU_ASSERT_EQUAL(Y(tune)(N, n, (int)dir, goal, &m, 0), 1);
    err = measure(N, n, 2 * N, m, (int)dir);
    CU_ASSERT(err <= goal);
  }
}

/* n must be even with no prime factor above 5. */
static int is_even_smooth5(INT n)
{
  if (n < 2 || n % 2 != 0)
    return 0;
  while (n % 2 == 0)
    n /= 2;
  while (n % 3 == 0)
    n /= 3;
  while (n % 5 == 0)
    n /= 5;
  return n == 1;
}

void Y(check_tune_plan)(void)
{
  static const INT Ns[] = {32, 64, 100, 128, 251, 256, 500, 512, 517};
  static const R goals[] = {K(1.0e-2), K(1.0e-4),  K(1.0e-6),
                            K(1.0e-8), K(1.0e-10), K(1.0e-12)};
  const R eps = Y(float_property)(NFFT_EPSILON);
  size_t iN, ig, dir;

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (ig = 0; ig < sizeof(goals) / sizeof(goals[0]); ig++)
      for (dir = 0; dir <= 1; dir++) {
        const INT N = Ns[iN];
        const R goal = goals[ig];
        INT n = 0;
        int m = 0;
        R attained = K(0.0);
        const int rc = Y(tune_plan)(N, (int)dir, goal, &n, &m, &attained);
        int m_chk = 0;

        CU_ASSERT(rc >= 0);
        CU_ASSERT(is_even_smooth5(n));
        CU_ASSERT(n > N);
        CU_ASSERT(m >= 1);
        CU_ASSERT(m <= (int)(n / 2 - 1));
        CU_ASSERT(attained > K(0.0));

        if (rc != 1)
          continue; /* below the reachable floor: the capped test covers it */

        CU_ASSERT(attained <= goal);
        /* tune must agree that this cut-off works at this size. */
        CU_ASSERT_EQUAL(Y(tune)(N, n, (int)dir, goal, &m_chk, 0), 1);
        CU_ASSERT(m_chk <= m);

        /* And the pair has to survive a real transform. The oracle is an
         * O(N*M) direct NDFT, so only the small geometries are run: what is
         * under test is the pair, not the bandwidth. */
        if (goal >= K(1.0e4) * eps && N <= (INT)128) {
          const R err = measure(N, n, 2 * N, m, (int)dir);
          if (err > goal)
            printf("\ntune_plan(N=" __D__ ", %s, goal=" __FE__ ") -> n=" __D__
                   " m=%d but err=" __FE__ "\n",
                   N, dir ? "adjoint" : "forward", goal, n, m, err);
          CU_ASSERT(err <= goal);
        }
      }
}

void Y(check_tune_plan_capped)(void)
{
  static const INT Ns[] = {64, 256};
  const R eps = Y(float_property)(NFFT_EPSILON);
  size_t iN, dir;
  int m = 0;
  INT n = 0;
  R attained = K(0.0);

  /* Bad arguments. */
  CU_ASSERT_EQUAL(Y(tune_plan)(0, TUNE_FWD, K(1.0e-6), &n, &m, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune_plan)(64, TUNE_FWD, K(0.0), &n, &m, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune_plan)(64, TUNE_FWD, K(1.0e-6), 0, &m, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune_plan)(64, TUNE_FWD, K(1.0e-6), &n, 0, &attained), -1);

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (dir = 0; dir <= 1; dir++) {
      const INT N = Ns[iN];
      INT n_cap = 0, n_loose = 0;
      int m_cap = 0, m_loose = 0;
      R a_cap = K(0.0), a_loose = K(0.0);

      /* Far below anything reachable: capped, not refused. */
      CU_ASSERT_EQUAL(Y(tune_plan)(N, (int)dir, eps / K(1.0e6), &n_cap, &m_cap,
                                   &a_cap),
                      0);
      CU_ASSERT(is_even_smooth5(n_cap));
      CU_ASSERT(m_cap >= 1);
      CU_ASSERT(a_cap > eps / K(1.0e6));

      /* The cap is honoured: asking for exactly it succeeds. */
      {
        INT n2 = 0;
        int m2 = 0;
        R a2 = K(0.0);
        CU_ASSERT_EQUAL(Y(tune_plan)(N, (int)dir, a_cap, &n2, &m2, &a2), 1);
        CU_ASSERT(a2 <= a_cap);
      }

      /* A loose goal costs no more oversampling than a tight one. */
      CU_ASSERT_EQUAL(Y(tune_plan)(N, (int)dir, K(1.0e-2), &n_loose, &m_loose,
                                   &a_loose),
                      1);
      CU_ASSERT(n_loose <= n_cap);
      CU_ASSERT(m_loose <= m_cap);
    }
}
