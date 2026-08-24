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

/* Two test-only NFFT solvers -- a permuting one and a deliberately slow one --
 * registered only through their own register functions, never by the library
 * roster.  See tests/test_solvers.h for what each proves.
 *
 * A coreless plan that borrows problem_nfft.x, the problem's owned copy:
 * SLEEPY -> AWAKE_ZERO reverses the M node-blocks of x in place, and
 * -> SLEEPY reverses them again to restore x exactly. apply/apply_adjoint are
 * no-ops -- this plan exercises the wakefulness and restore machinery, and
 * only its shape matters, not the transform it would compute. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "test_solvers.h"

#undef X
#define X(name) NFFT(name)

/* When set, the SLEEPY restore is skipped and the x-copy is left permuted, so
 * the A()-gated restore guard Y(nfft_x_verify) must fire. */
int Y(nfft_perm_break_restore) = 0;

typedef struct
{
  plan super;
  const problem_nfft *pn;
  INT M, d;
  int permuted;
} perm_plan;

/* Its own inverse: reversing twice restores the original order. */
static void reverse(R *x, INT d, INT M) {
  INT i, t;
  for (i = 0; i < M / 2; i++)
    for (t = 0; t < d; t++) {
      R z = x[i * d + t];
      x[i * d + t] = x[(M - 1 - i) * d + t];
      x[(M - 1 - i) * d + t] = z;
    }
}

static void awake(plan *ego, int w) {
  perm_plan *q = (perm_plan *)ego;
  R *x = (R *)q->pn->x;
  if (w >= PLNR_AWAKE_ZERO && ego->awake_state == PLNR_SLEEPY) {
    reverse(x, q->d, q->M);
    q->permuted = 1;
  } else if (w == PLNR_SLEEPY && ego->awake_state >= PLNR_AWAKE_ZERO) {
    if (q->permuted && !Y(nfft_perm_break_restore))
      reverse(x, q->d, q->M);
    q->permuted = 0;
  }
}

static void apply(const plan *e, const problem *p) {
  (void)e;
  (void)p;
}

static void print(const plan *e, printer *pr) {
  (void)e;
  pr->print(pr, "(nfft_solver_perm_test pcost=0)");
}

static const plan_adt perm_plan_adt = {apply, awake, print, 0, apply};

/* Non-static: direct-drive tests build the plan without a planner. */
plan *Y(nfft_perm_test_mkplan)(const problem *p) {
  const problem_nfft *pn = (const problem_nfft *)p;
  perm_plan *q = (perm_plan *)Y(plan_create)(sizeof(perm_plan), &perm_plan_adt);
  q->super.pcost = 1.0;
  q->pn = pn;
  q->M = pn->M;
  q->d = pn->sz->rnk;
  q->permuted = 0;
  return &q->super;
}

static plan *solver_mkplan(const solver *s, const problem *p, planner *pl) {
  (void)s;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_NFFT)
    return 0;
  return Y(nfft_perm_test_mkplan)(p);
}

static const solver_adt perm_solver_adt = {NFFT_PROBLEM_NFFT, 0, solver_mkplan};

void Y(nfft_solver_perm_test_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &perm_solver_adt));
}

/* A test-only NFFT solver with a huge analytic cost that counts its applies.
 * The measured race must prune it by estimate and never time it. */
INT Y(nfft_slow_test_applies) = 0;

static void slow_apply(const plan *e, const problem *p) {
  (void)e;
  (void)p;
  Y(nfft_slow_test_applies)++;
}

static void slow_print(const plan *e, printer *pr) {
  (void)e;
  pr->print(pr, "(nfft_solver_slow_test)");
}

static const plan_adt slow_plan_adt = {slow_apply, 0, slow_print, 0, slow_apply};

static plan *slow_mkplan(const solver *s, const problem *p, planner *pl) {
  plan *q;
  (void)s;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_NFFT)
    return 0;
  q = Y(plan_create)(sizeof(plan), &slow_plan_adt);
  q->pcost = 1.0e18;
  return q;
}

static const solver_adt slow_solver_adt = {NFFT_PROBLEM_NFFT, 0, slow_mkplan};

void Y(nfft_solver_slow_test_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &slow_solver_adt));
}
