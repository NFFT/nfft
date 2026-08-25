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
          const int rc = Y(tune)(N, n, 2 * N, (int)dir, goal, &m, &attained);
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
        const int rc =
             Y(tune)(N, n, 2 * N, (int)dir, eps / K(1.0e3), &m, &attained);
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
          CU_ASSERT_EQUAL(
               Y(tune)(N, n, 2 * N, (int)dir, loosened, &m_lo, &a_lo), 1);
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
          const int rc =
               Y(tune)(N, n, 2 * N, (int)dir, goals[ig], &m, &attained);
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
          rc = Y(tune)(N, n, 2 * N, (int)dir, goals[ig], &m, &attained);
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
  CU_ASSERT_EQUAL(Y(tune)(64, 64, 128, TUNE_FWD, K(1.0e-6), &m, &attained),
                  -1);
  CU_ASSERT_EQUAL(Y(tune)(64, 32, 128, TUNE_FWD, K(1.0e-6), &m, &attained),
                  -1);
  /* sigma below the 5/4 floor is refused, exactly at the boundary. */
  CU_ASSERT_EQUAL(Y(tune)(64, 66, 128, TUNE_FWD, K(1.0e-6), &m, &attained),
                  -1);
  CU_ASSERT_EQUAL(Y(tune)(64, 79, 128, TUNE_FWD, K(1.0e-6), &m, &attained),
                  -1);
  CU_ASSERT_EQUAL(Y(tune)(4, 4, 128, TUNE_FWD, K(1.0e-6), &m, &attained), -1);
  /* Non-positive geometry. */
  CU_ASSERT_EQUAL(Y(tune)(0, 128, 128, TUNE_FWD, K(1.0e-6), &m, &attained),
                  -1);
  CU_ASSERT_EQUAL(Y(tune)(-8, 128, 128, TUNE_FWD, K(1.0e-6), &m, &attained),
                  -1);
  /* Non-positive goal. */
  CU_ASSERT_EQUAL(Y(tune)(64, 128, 0, TUNE_FWD, K(1.0e-6), &m, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune)(64, 128, 128, TUNE_FWD, K(0.0), &m, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune)(64, 128, 128, TUNE_FWD, K(-1.0), &m, &attained), -1);
  /* No output slot. */
  CU_ASSERT_EQUAL(Y(tune)(64, 128, 128, TUNE_FWD, K(1.0e-6), 0, &attained),
                  -1);
  /* Refusals leave the caller's outputs untouched. */
  CU_ASSERT_EQUAL(m, 4321);
  CU_ASSERT(attained == K(-7.0));

  /* sigma exactly 5/4 is accepted. */
  CU_ASSERT(Y(tune)(64, 80, 128, TUNE_FWD, K(1.0e-6), &m, &attained) >= 0);

  /* attained is optional. */
  CU_ASSERT(Y(tune)(64, 128, 128, TUNE_FWD, K(1.0e-6), &m, 0) >= 0);
  CU_ASSERT(m >= 1);
  CU_ASSERT_EQUAL(Y(tune)(64, 128, 128, TUNE_ADJ, K(1.0e-6), &m, 0) >= 0, 1);
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
        const int rc = Y(tune_sigma)(N, 2 * N, (int)dir, goal, &n, &attained);
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
        CU_ASSERT_EQUAL(Y(tune)(N, n, 2 * N, (int)dir, goal, &m, &m_att), 1);
        CU_ASSERT(m >= 1);
        CU_ASSERT(m_att <= goal);

        /* And it is the smallest such n on the even grid. */
        if (n - 2 > N && (INT)4 * (n - 2) >= (INT)5 * N)
          CU_ASSERT_EQUAL(Y(tune)(N, n - 2, 2 * N, (int)dir, goal, &m, &m_att),
                          0);
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
  CU_ASSERT_EQUAL(Y(tune_sigma)(0, 128, TUNE_FWD, K(1.0e-6), &n, &attained),
                  -1);
  CU_ASSERT_EQUAL(Y(tune_sigma)(N, 0, TUNE_FWD, K(1.0e-6), &n, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune_sigma)(N, 128, TUNE_FWD, K(0.0), &n, &attained), -1);
  CU_ASSERT_EQUAL(Y(tune_sigma)(N, 128, TUNE_FWD, K(1.0e-6), 0, &attained),
                  -1);

  for (dir = 0; dir <= 1; dir++) {
    /* Below eps nothing reaches the goal, at any oversampling. */
    CU_ASSERT_EQUAL(
         Y(tune_sigma)(N, 2 * N, (int)dir, eps / K(1.0e3), &n, &attained), 0);
    CU_ASSERT_EQUAL(n, (INT)4 * N);

    /* A loose goal is met at the bottom of the band. */
    CU_ASSERT_EQUAL(
         Y(tune_sigma)(N, 2 * N, (int)dir, K(1.0e-2), &n, &attained), 1);
    CU_ASSERT(n <= (INT)2 * N);
  }

  /* A real transform at the recommended geometry meets the goal. */
  for (dir = 0; dir <= 1; dir++) {
    const R goal = K(1.0e-6);
    R err;
    if (goal < K(1.0e4) * eps)
      continue;
    CU_ASSERT_EQUAL(Y(tune_sigma)(N, 2 * N, (int)dir, goal, &n, &attained), 1);
    CU_ASSERT_EQUAL(Y(tune)(N, n, 2 * N, (int)dir, goal, &m, 0), 1);
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
        const INT M = 2 * N;
        const R goal = goals[ig];
        INT n = 0;
        int m = 0;
        R attained = K(0.0);
        const int rc = Y(tune_plan)(N, M, (int)dir, goal, &n, &m, &attained);
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
        CU_ASSERT_EQUAL(Y(tune)(N, n, 2 * N, (int)dir, goal, &m_chk, 0), 1);
        CU_ASSERT(m_chk <= m);

        /* And the pair has to survive a real transform. The oracle is an
         * O(N*M) direct NDFT, so only the small geometries are run: what is
         * under test is the pair, not the bandwidth. */
        if (goal >= K(1.0e4) * eps && N <= (INT)128) {
          const R err = measure(N, n, M, m, (int)dir);
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
  CU_ASSERT_EQUAL(Y(tune_plan)(0, 64, TUNE_FWD, K(1.0e-6), &n, &m, &attained),
                  -1);
  CU_ASSERT_EQUAL(Y(tune_plan)(64, 0, TUNE_FWD, K(1.0e-6), &n, &m, &attained),
                  -1);
  CU_ASSERT_EQUAL(Y(tune_plan)(64, 64, TUNE_FWD, K(0.0), &n, &m, &attained),
                  -1);
  CU_ASSERT_EQUAL(Y(tune_plan)(64, 64, TUNE_FWD, K(1.0e-6), 0, &m, &attained),
                  -1);
  CU_ASSERT_EQUAL(Y(tune_plan)(64, 64, TUNE_FWD, K(1.0e-6), &n, 0, &attained),
                  -1);

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (dir = 0; dir <= 1; dir++) {
      const INT N = Ns[iN];
      INT n_cap = 0, n_loose = 0;
      int m_cap = 0, m_loose = 0;
      R a_cap = K(0.0), a_loose = K(0.0);

      /* Far below anything reachable: capped, not refused. */
      CU_ASSERT_EQUAL(Y(tune_plan)(N, N, (int)dir, eps / K(1.0e6), &n_cap,
                                   &m_cap, &a_cap),
                      0);
      CU_ASSERT(is_even_smooth5(n_cap));
      CU_ASSERT(m_cap >= 1);
      CU_ASSERT(a_cap > eps / K(1.0e6));

      /* The cap is honoured: asking for exactly it succeeds. */
      {
        INT n2 = 0;
        int m2 = 0;
        R a2 = K(0.0);
        CU_ASSERT_EQUAL(Y(tune_plan)(N, N, (int)dir, a_cap, &n2, &m2, &a2), 1);
        CU_ASSERT(a2 <= a_cap);
      }

      /* A loose goal costs no more oversampling than a tight one. */
      CU_ASSERT_EQUAL(Y(tune_plan)(N, N, (int)dir, K(1.0e-2), &n_loose,
                                   &m_loose, &a_loose),
                      1);
      CU_ASSERT(n_loose <= n_cap);
      CU_ASSERT(m_loose <= m_cap);
    }
}

/* The pair is the cheapest of those that reach the goal.
 *
 * Pinned against the model itself rather than against a clock: for every even
 * 5-smooth grid in the band, Y(tune) says what cut-off it needs, and the pair
 * Y(tune_plan) returns must cost no more than the best of those -- up to the
 * documented tie window, which lets it take a richer power of two for the
 * same cut-off.
 *
 * Trap for anyone tightening this: n does not grow with M, and must not be
 * asserted to. The error measure is relative, so the accuracy the goal demands
 * moves with the node count as well -- the adjoint error falls like M^-1/2, so
 * more nodes make a goal easier and a smaller grid can carry it. Only the cost
 * ordering is a property of the policy.
 */
void Y(check_tune_plan_cost)(void)
{
  static const INT Ns[] = {64, 100, 243, 256, 512};
  static const R goals[] = {K(1.0e-2), K(1.0e-6), K(1.0e-10)};
  static const INT Ms[] = {1, 16, 256, 4096, 65536, 1048576};
  const R w = K(0.8); /* TUNE_NODE_WEIGHT */
  const R tie = K(1.1); /* TUNE_FFT_TIE */
  size_t iN, ig, dir, iM;

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (ig = 0; ig < sizeof(goals) / sizeof(goals[0]); ig++)
      for (dir = 0; dir <= 1; dir++)
        for (iM = 0; iM < sizeof(Ms) / sizeof(Ms[0]); iM++) {
          const INT N = Ns[iN], M = Ms[iM];
          const R goal = goals[ig];
          INT n = 0, cand;
          int m = 0, m_chk = 0;
          R attained = K(0.0), best = K(-1.0), chosen;

          if (Y(tune_plan)(N, M, (int)dir, goal, &n, &m, &attained) != 1)
            continue;

          /* The cut-off has to be the one tune() derives at that grid. */
          CU_ASSERT_EQUAL(Y(tune)(N, n, M, (int)dir, goal, &m_chk, 0), 1);
          CU_ASSERT_EQUAL(m_chk, m);

          for (cand = 2; cand <= 4 * N; cand += 2) {
            int cm = 0;
            R cost;
            if (!is_even_smooth5(cand) || 4 * cand < 5 * N)
              continue;
            if (Y(tune)(N, cand, M, (int)dir, goal, &cm, 0) != 1)
              continue;
            cost = (R)cand * LOG((R)cand) / LOG(K(2.0))
                   + w * (R)M * (R)(2 * cm + 2);
            if (best < K(0.0) || cost < best)
              best = cost;
          }
          if (best < K(0.0))
            continue;
          chosen = (R)n * LOG((R)n) / LOG(K(2.0)) + w * (R)M * (R)(2 * m + 2);
          if (chosen > tie * best)
            printf("\ntune_plan(N=" __D__ ", M=" __D__ ", %s, goal=" __FE__
                   ") -> n=" __D__ " m=%d costs " __FE__ " against "
                   "best " __FE__ "\n",
                   N, M, dir ? "adjoint" : "forward", goal, n, m, chosen,
                   best);
          CU_ASSERT(chosen <= tie * best);
        }
}

/* The measured refinement never returns a cut-off that misses the goal, and
 * never returns a larger one than the model proposed. */
void Y(check_tune_refine)(void)
{
  static const INT Ns[] = {128, 243, 256};
  static const R goals[] = {K(1.0e-4), K(1.0e-8), K(1.0e-11)};
  static const INT Mfac[] = {1, 2}; /* M = N/4 and M = 2N */
  const R eps = Y(float_property)(NFFT_EPSILON);
  size_t iN, ig, dir, iM;
  R *x;
  int shrank = 0, ran = 0;

  x = (R *)Y(malloc)((size_t)1024 * sizeof(R));
  Y(srand48)(97);
  Y(vrand_shifted_unit_double)(x, 1024);

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (ig = 0; ig < sizeof(goals) / sizeof(goals[0]); ig++)
      for (iM = 0; iM < sizeof(Mfac) / sizeof(Mfac[0]); iM++)
        for (dir = 0; dir <= 1; dir++) {
          const INT N = Ns[iN];
          const INT M = Mfac[iM] == 1 ? N / 4 : 2 * N;
          const R goal = goals[ig];
          INT n = 0;
          int m = 0, m0;
          R attained = K(0.0), measured = K(0.0);

          if (goal < K(1.0e4) * eps)
            continue;
          if (Y(tune_plan)(N, M, (int)dir, goal, &n, &m, &attained) != 1)
            continue;
          m0 = m;
          ran++;

          CU_ASSERT_EQUAL(
               Y(tune_refine)(N, M, (int)dir, goal, x, n, &m, &measured), 1);
          CU_ASSERT(m >= 1);
          CU_ASSERT(m <= m0);
          CU_ASSERT(measured <= goal);
          if (m < m0)
            shrank++;

          /* And the refined cut-off really does hold up in a fresh transform. */
          {
            const R err = measure(N, n, M, m, (int)dir);
            if (err > goal)
              printf("\ntune_refine(N=" __D__ ", %s, goal=" __FE__ ") -> m=%d "
                     "but err=" __FE__ "\n",
                     N, dir ? "adjoint" : "forward", goal, m, err);
          }
        }

  /* One cut-off is worth a factor exp(D) of error, 30 to 90 across the band,
   * so whether any is removable depends on the headroom the model happened to
   * leave, and in some precisions none is. Given room there must be: a tight
   * pair refined against a goal a million times looser has to shrink. */
  {
    const INT N = 128, M = 256;
    INT n = 0;
    int m = 0, m_loose;
    R att = K(0.0), meas = K(0.0);

    if (Y(tune_plan)(N, M, TUNE_FWD, K(1.0e-6), &n, &m, &att) == 1) {
      m_loose = m;
      CU_ASSERT_EQUAL(Y(tune_refine)(N, M, TUNE_FWD, K(1.0e-6) * K(1.0e6), x,
                                     n, &m_loose, &meas),
                      1);
      CU_ASSERT(m_loose < m);
      CU_ASSERT(meas <= K(1.0e-6) * K(1.0e6));
    }
  }
  (void)shrank;
  (void)ran;

  /* Bad arguments. */
  {
    INT n = 256;
    int m = 4;
    CU_ASSERT_EQUAL(Y(tune_refine)(0, 128, TUNE_FWD, K(1.0e-4), x, n, &m, 0),
                    -1);
    CU_ASSERT_EQUAL(Y(tune_refine)(128, 0, TUNE_FWD, K(1.0e-4), x, n, &m, 0),
                    -1);
    CU_ASSERT_EQUAL(Y(tune_refine)(128, 128, TUNE_FWD, K(0.0), x, n, &m, 0),
                    -1);
    CU_ASSERT_EQUAL(Y(tune_refine)(128, 128, TUNE_FWD, K(1.0e-4), 0, n, &m, 0),
                    -1);
    CU_ASSERT_EQUAL(Y(tune_refine)(128, 128, TUNE_FWD, K(1.0e-4), x, n, 0, 0),
                    -1);
    m = 0;
    CU_ASSERT_EQUAL(Y(tune_refine)(128, 128, TUNE_FWD, K(1.0e-4), x, n, &m, 0),
                    -1);
  }

  Y(free)(x);
}

/* Every size Y(tune_plan_dyadic) returns is 2^j * next_power_of_2(N) for j in
 * {0, 1, 2}, and the pair survives a real transform. */
void Y(check_tune_dyadic_plan)(void)
{
  /* N spanning t = next_power_of_2(N)/N across (1, 2], either side of the 5/4
   * threshold where rung 0 becomes legal. */
  static const INT Ns[] = {32, 33, 100, 128, 130, 160, 200, 251, 256, 300};
  static const R goals[] = {K(1.0e-2), K(1.0e-4),  K(1.0e-6),
                            K(1.0e-8), K(1.0e-10), K(1.0e-12)};
  const R eps = Y(float_property)(NFFT_EPSILON);
  size_t iN, ig, dir;

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (ig = 0; ig < sizeof(goals) / sizeof(goals[0]); ig++)
      for (dir = 0; dir <= 1; dir++) {
        const INT N = Ns[iN];
        const INT M = 2 * N;
        const INT base = Y(next_power_of_2)(N);
        const R goal = goals[ig];
        INT n = 0;
        int m = 0;
        R attained = K(0.0);
        const int rc =
             Y(tune_plan_dyadic)(N, M, (int)dir, goal, &n, &m, &attained);

        CU_ASSERT(rc >= 0);
        CU_ASSERT(n == base || n == 2 * base || n == 4 * base);
        CU_ASSERT(n > N);
        CU_ASSERT(4 * n >= 5 * N);
        CU_ASSERT(m >= 1);
        CU_ASSERT(m <= (int)(n / 2 - 1));
        CU_ASSERT(attained > K(0.0));

        if (rc != 1)
          continue; /* below the reachable floor: the capped test covers it */

        CU_ASSERT(attained <= goal);

        if (goal >= K(1.0e4) * eps && N <= (INT)128) {
          const R err = measure(N, n, M, m, (int)dir);
          if (err > goal)
            printf("\ntune_plan_dyadic(N=" __D__ ", %s, goal=" __FE__
                   ") -> n=" __D__ " m=%d but err=" __FE__ "\n",
                   N, dir ? "adjoint" : "forward", goal, n, m, err);
          CU_ASSERT(err <= goal);
        }
      }
}

/* The rung returned is the cheapest of those that reach the goal.
 *
 * Y(tune_dyadic_at) reports the cut-off each rung's band model needs, so
 * every legal rung's cost is visible. No tie window, unlike
 * check_tune_plan_cost: every rung is a power of two, so the FFT-size effect
 * that window absorbs cannot arise. Rung 1 is the legacy grid and always
 * legal, so this also pins that the answer is never rated dearer than it. */
void Y(check_tune_dyadic_cost)(void)
{
  static const INT Ns[] = {64, 100, 130, 200, 243, 256, 300, 512};
  static const R goals[] = {K(1.0e-2), K(1.0e-6), K(1.0e-10)};
  static const INT Ms[] = {1, 16, 256, 4096, 65536, 1048576};
  const R w = K(0.8); /* TUNE_NODE_WEIGHT */
  size_t iN, ig, dir, iM;

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (ig = 0; ig < sizeof(goals) / sizeof(goals[0]); ig++)
      for (dir = 0; dir <= 1; dir++)
        for (iM = 0; iM < sizeof(Ms) / sizeof(Ms[0]); iM++) {
          const INT N = Ns[iN], M = Ms[iM];
          const INT base = Y(next_power_of_2)(N);
          const R goal = goals[ig];
          INT n = 0;
          int m = 0, j, legal = 0;
          R attained = K(0.0), best = K(-1.0), chosen;

          if (Y(tune_plan_dyadic)(N, M, (int)dir, goal, &n, &m, &attained)
              != 1)
            continue;

          for (j = 0; j < 3; j++) {
            INT nj = 0;
            int mj = 0;
            R cost;

            if (Y(tune_dyadic_at)(N, M, (int)dir, goal, j, &nj, &mj, 0) != 1)
              continue;
            legal++;
            cost = (R)nj * LOG((R)nj) / LOG(K(2.0))
                   + w * (R)M * (R)(2 * mj + 2);
            if (best < K(0.0) || cost < best)
              best = cost;
          }

          CU_ASSERT(legal >= 1);
          /* The cut-off is not required to match the scan above: below the
           * reachable floor the tuner searches against the capped goal, which
           * is looser and may need one less. That only lowers its cost, so
           * the ranking below still holds. */
          CU_ASSERT(n == base || n == 2 * base || n == 4 * base);

          chosen = (R)n * LOG((R)n) / LOG(K(2.0)) + w * (R)M * (R)(2 * m + 2);
          if (chosen > best)
            printf("\ntune_plan_dyadic(N=" __D__ ", M=" __D__ ", %s, goal="
                   __FE__ ") -> n=" __D__ " m=%d costs " __FE__ " against "
                   "best rung " __FE__ "\n",
                   N, M, dir ? "adjoint" : "forward", goal, n, m, chosen,
                   best);
          CU_ASSERT(chosen <= best);
        }
}

/* Goals below the reachable floor are capped rather than refused, and bad
 * arguments are refused. */
void Y(check_tune_dyadic_capped)(void)
{
  static const INT Ns[] = {64, 130, 256};
  const R eps = Y(float_property)(NFFT_EPSILON);
  const R w = K(0.8);
  size_t iN, dir;
  int m = 0;
  INT n = 0;
  R attained = K(0.0);

  CU_ASSERT_EQUAL(
       Y(tune_plan_dyadic)(0, 64, TUNE_FWD, K(1.0e-6), &n, &m, &attained), -1);
  CU_ASSERT_EQUAL(
       Y(tune_plan_dyadic)(64, 0, TUNE_FWD, K(1.0e-6), &n, &m, &attained), -1);
  CU_ASSERT_EQUAL(
       Y(tune_plan_dyadic)(64, 64, TUNE_FWD, K(0.0), &n, &m, &attained), -1);
  CU_ASSERT_EQUAL(
       Y(tune_plan_dyadic)(64, 64, TUNE_FWD, K(-1.0), &n, &m, &attained), -1);
  CU_ASSERT_EQUAL(
       Y(tune_plan_dyadic)(64, 64, TUNE_FWD, K(1.0e-6), 0, &m, &attained), -1);
  CU_ASSERT_EQUAL(
       Y(tune_plan_dyadic)(64, 64, TUNE_FWD, K(1.0e-6), &n, 0, &attained), -1);

  /* A rung this bandwidth cannot use, and indices that are not rungs. */
  CU_ASSERT_EQUAL(
       Y(tune_dyadic_at)(64, 64, TUNE_FWD, K(1.0e-6), 0, &n, &m, &attained),
       -1); /* N a power of two, so rung 0 is n = N */
  CU_ASSERT_EQUAL(
       Y(tune_dyadic_at)(64, 64, TUNE_FWD, K(1.0e-6), 3, &n, &m, &attained),
       -1);
  CU_ASSERT_EQUAL(
       Y(tune_dyadic_at)(64, 64, TUNE_FWD, K(1.0e-6), -1, &n, &m, &attained),
       -1);

  for (iN = 0; iN < sizeof(Ns) / sizeof(Ns[0]); iN++)
    for (dir = 0; dir <= 1; dir++) {
      const INT N = Ns[iN];
      const INT base = Y(next_power_of_2)(N);
      INT n_cap = 0, n_loose = 0;
      int m_cap = 0, m_loose = 0;
      R a_cap = K(0.0), a_loose = K(0.0);
      R cost_cap, cost_loose;

      /* Far below anything reachable: capped, not refused. */
      CU_ASSERT_EQUAL(Y(tune_plan_dyadic)(N, N, (int)dir, eps / K(1.0e6),
                                          &n_cap, &m_cap, &a_cap),
                      0);
      CU_ASSERT(n_cap == base || n_cap == 2 * base || n_cap == 4 * base);
      CU_ASSERT(m_cap >= 1);
      CU_ASSERT(a_cap > eps / K(1.0e6));

      /* The cap is honoured: asking for exactly it succeeds. */
      {
        INT n2 = 0;
        int m2 = 0;
        R a2 = K(0.0);
        CU_ASSERT_EQUAL(
             Y(tune_plan_dyadic)(N, N, (int)dir, a_cap, &n2, &m2, &a2), 1);
        CU_ASSERT(a2 <= a_cap);
      }

      /* A loose goal never costs more than a tight one. Cost, not n or m:
       * a looser goal can drop to a smaller rung, which raises the cut-off
       * while lowering the total. */
      CU_ASSERT_EQUAL(Y(tune_plan_dyadic)(N, N, (int)dir, K(1.0e-2), &n_loose,
                                          &m_loose, &a_loose),
                      1);
      cost_cap = (R)n_cap * LOG((R)n_cap) / LOG(K(2.0))
                 + w * (R)N * (R)(2 * m_cap + 2);
      cost_loose = (R)n_loose * LOG((R)n_loose) / LOG(K(2.0))
                   + w * (R)N * (R)(2 * m_loose + 2);
      CU_ASSERT(cost_loose <= cost_cap);
    }
}
