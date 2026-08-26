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

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

/* The rank-0 solver. An all-unit-axis NFFT problem (every N_t == 1) compresses
 * to sz->rnk == 0: one coefficient k=0 whose phase exp(-2pi i * 0 * x) == 1 at
 * every node. Forward broadcasts f_hat[0] to all M outputs, adjoint reduces
 * sum_j f[j] into f_hat[0]. */

typedef struct
{
  plan super;
  const char *reg_nam;
} rnk0_plan;

static void apply(const plan *ego, const problem *p) /* forward broadcast */
{
  const problem_nfft *pn = (const problem_nfft *)p;
  const C v = pn->f_hat[0];
  const INT M = pn->M;
  INT j;
  (void)ego;
  for (j = 0; j < M; j++)
    pn->f[j] = v;
}

static void apply_adjoint(const plan *ego, const problem *p) /* adjoint reduce */
{
  const problem_nfft *pn = (const problem_nfft *)p;
  const INT M = pn->M;
  C acc = K(0.0);
  INT j;
  (void)ego;
  for (j = 0; j < M; j++)
    acc += pn->f[j];
  pn->f_hat[0] = acc;
}

static void print(const plan *ego, printer *pr) {
  const rnk0_plan *pln = (const rnk0_plan *)ego;
  pr->print(pr, "(%s pcost=%D)", pln->reg_nam, (INT)pln->super.pcost);
}

/* Owns only the base allocation: no awake hook, no destroy. */
static const plan_adt rnk0_plan_adt = {apply, 0, print, 0,
                                          apply_adjoint};

static plan *mkplan_rnk0(const solver *ego, const problem *p, planner *pl) {
  const problem_nfft *pn = (const problem_nfft *)p;
  rnk0_plan *pln;
  (void)ego;
  (void)pl; /* the sole rank-0 solver: ungated, never ties */
  if (p->adt->kind != NFFT_PROBLEM_NFFT)
    return 0;
  if (pn->sz->rnk != 0)
    return 0;
  pln = (rnk0_plan *)Y(plan_create)(sizeof(rnk0_plan), &rnk0_plan_adt);
  pln->reg_nam = "nfft_solver_const_0d";
  pln->super.pcost = 2.0 * (double)pn->M; /* one broadcast + one reduce pass */
  return &pln->super;
}

static const solver_adt rnk0_solver_adt = {NFFT_PROBLEM_NFFT, 0,
                                              mkplan_rnk0};

void Y(nfft_solver_rnk0_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &rnk0_solver_adt));
}
