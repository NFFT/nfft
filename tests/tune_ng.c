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
