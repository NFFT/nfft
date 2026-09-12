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

#include <complex.h>

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

#include <CUnit/CUnit.h>

#include "ops.h"

/* The helpers are pure arithmetic on four doubles; they exist so that solvers
 * never write the composition by hand. Every one is exercised here because
 * they are too small to debug through a solver. */
void Y(check_ops_helpers)(void)
{
  opcnt a, b, c;

  Y(ops_zero)(&a);
  CU_ASSERT_DOUBLE_EQUAL(a.add, 0.0, 0.0);
  CU_ASSERT_DOUBLE_EQUAL(a.mul, 0.0, 0.0);
  CU_ASSERT_DOUBLE_EQUAL(a.fma, 0.0, 0.0);
  CU_ASSERT_DOUBLE_EQUAL(a.other, 0.0, 0.0);

  Y(ops_other)(7, &a);
  CU_ASSERT_DOUBLE_EQUAL(a.other, 7.0, 0.0);
  CU_ASSERT_DOUBLE_EQUAL(a.add, 0.0, 0.0);

  a.add = 1.0;
  a.mul = 2.0;
  a.fma = 3.0;
  a.other = 4.0;
  Y(ops_cpy)(&a, &b);
  CU_ASSERT_DOUBLE_EQUAL(b.mul, 2.0, 0.0);

  Y(ops_add)(&a, &b, &c);
  CU_ASSERT_DOUBLE_EQUAL(c.add, 2.0, 0.0);
  CU_ASSERT_DOUBLE_EQUAL(c.other, 8.0, 0.0);

  Y(ops_zero)(&c);
  Y(ops_add2)(&a, &c);
  Y(ops_add2)(&a, &c);
  CU_ASSERT_DOUBLE_EQUAL(c.fma, 6.0, 0.0);

  Y(ops_madd)(3, &a, &b, &c); /* 3*a + b */
  CU_ASSERT_DOUBLE_EQUAL(c.add, 4.0, 0.0);
  CU_ASSERT_DOUBLE_EQUAL(c.other, 16.0, 0.0);

  Y(ops_zero)(&c);
  Y(ops_madd2)(5, &a, &c);
  CU_ASSERT_DOUBLE_EQUAL(c.mul, 10.0, 0.0);
}

/* The estimator's weights are FFTW's, unchanged: add + mul + other, with fma
 * counted once where the hardware fuses and twice where it does not. */
void Y(check_ops_estimate)(void)
{
  planner *pl = Y(planner_create)();
  plan *pln = Y(plan_create)(sizeof(plan), 0);
  double got;

  CU_ASSERT_DOUBLE_EQUAL(pln->ops.add, 0.0, 0.0);
  CU_ASSERT_DOUBLE_EQUAL(pln->pcost, 0.0, 0.0);

  pln->ops.add = 10.0;
  pln->ops.mul = 20.0;
  pln->ops.fma = 5.0;
  pln->ops.other = 100.0;
  got = Y(iestimate_cost)(pl, pln, 0);
#if HAVE_FMA
  CU_ASSERT_DOUBLE_EQUAL(got, 135.0, 1e-9);
#else
  CU_ASSERT_DOUBLE_EQUAL(got, 140.0, 1e-9);
#endif

  Y(free)(pln);
  Y(planner_destroy)(pl);
}

static void noop_apply(const plan *ego, const problem *p)
{
  (void)ego;
  (void)p;
}
static void noop_awake(plan *ego, int wakefulness)
{
  (void)ego;
  (void)wakefulness;
}
static void noop_print(const plan *ego, printer *pr)
{
  (void)ego;
  pr->print(pr, "(noop)");
}
static void noop_destroy(plan *ego)
{
  (void)ego;
}
static const plan_adt test_noop_plan_adt = {noop_apply, noop_awake, noop_print,
                                            noop_destroy, noop_apply};

/* FFTW's condition (kernel/planner.c:446) and its per-candidate fallback
 * (462-465): a negative reading means no usable clock, so that one candidate
 * falls back to its estimate and the race carries on. */
void Y(check_ops_evaluate_plan)(void)
{
  planner *pl = Y(planner_create)();
  plan *pln = Y(plan_create)(sizeof(plan), &test_noop_plan_adt);

  pln->ops.mul = 42.0;
  CU_ASSERT_DOUBLE_EQUAL(pln->pcost, 0.0, 0.0);

  /* Estimate mode: pcost is the estimate, whatever the clock does. */
  pl->flags.u = PLNR_ESTIMATE;
  Y(evaluate_plan)(pl, pln, 0);
  CU_ASSERT_DOUBLE_EQUAL(pln->pcost, 42.0, 1e-9);

  /* A plan that already carries a cost is not re-evaluated when
   * BELIEVE_PCOST is set. */
  pl->flags.u = PLNR_BELIEVE_PCOST;
  pln->pcost = 7.0;
  Y(evaluate_plan)(pl, pln, 0);
  CU_ASSERT_DOUBLE_EQUAL(pln->pcost, 7.0, 1e-9);

  Y(free)(pln);
  Y(planner_destroy)(pl);
}

/* No solver may set pcost any more; the contract is ops. A fresh plan from any
 * solver carries a zero pcost and a non-zero estimate. */
void Y(check_ops_solvers_leave_pcost_zero)(void)
{
  planner *pl = Y(planner_create)();
  const INT N = 16, n = 32;
  C fh = K(1.0), g = K(0.0);
  problem *p;
  plan *pln;

  Y(deconv_solvers_register)(pl);
  p = Y(mkproblem_deconv)(1, &N, 0, &n, 6, Y(get_window_id)(), 1, &fh, &g);
  pln = Y(planner_mkplan)(pl, p);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  CU_ASSERT_DOUBLE_EQUAL(pln->pcost, 0.0, 0.0);
  CU_ASSERT_TRUE(Y(iestimate_cost)(pl, pln, p) > 0.0);
  Y(plan_destroy)(pln);
  Y(problem_destroy)(p);
  Y(planner_destroy)(pl);
}
