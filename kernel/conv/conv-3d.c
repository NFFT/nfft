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

/* 3D CONV solver: step C of the fast NFFT decomposition, the node convolution
 * (matrix B). Forward sums the oversampled grid g against the window psi at
 * each nonequispaced node x_j; the adjoint scatter-adds f onto g with the same
 * psi weights. psi and the wrapped window start u depend on x/window/n/N/m and
 * are built at awake, so apply does no window evaluation and no FLOOR/LRINT. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "conv.h"

#include <math.h> /* pow */

typedef struct {
  plan super;
  INT n0, n1, n2, N0, N1, N2, M; /* geometry captured at mkplan */
  int m, window, analytic;
  const R *x; /* borrowed alias of the problem's nodes */
  R *psi; /* length M*3*(2m+2): psi[(j*3+t)*(2m+2)+lj]; NULL while SLEEPY */
  INT *u; /* length M*3: wrapped window start, u[j*3+t]; NULL while SLEEPY */
  INT *perm; /* length M, node of row k in key order; NULL while SLEEPY and
              * for caller order */
  int sorted; /* this plan builds perm */
} conv_3d_plan;

static void fill_u(conv_3d_plan *pln)
{
  const INT nn[3] = {pln->n0, pln->n1, pln->n2};
  Y(conv_fill_u)(pln->x, pln->perm, nn, 3, pln->m, pln->M, pln->u);
}

static void fill_psi(conv_3d_plan *pln)
{
  const INT nn[3] = {pln->n0, pln->n1, pln->n2};
  const INT NN[3] = {pln->N0, pln->N1, pln->N2};
  Y(conv_fill_psi)(pln->window, pln->analytic, pln->x, pln->perm, nn, NN, 3,
                pln->m, pln->M, pln->psi);
}

/* The tables live between a wake and the return to SLEEPY (ADR-0009). Node
 * order and window starts are indices: they steer which grid cells apply
 * touches, so they hold true values at AWAKE_ZERO too, and the measured race
 * sees the locality it is choosing between. psi costs a window evaluation per
 * tap and changes no address, so AWAKE_ZERO gets zeros. */
static void awake(plan *ego_, int wakefulness)
{
  conv_3d_plan *pln = (conv_3d_plan *)ego_;
  const size_t M = (size_t)pln->M, len = (size_t)(2 * pln->m + 2);
  if (wakefulness == PLNR_SLEEPY) {
    Y(free)(pln->psi);
    Y(free)(pln->u);
    Y(free)(pln->perm);
    pln->psi = 0;
    pln->u = 0;
    pln->perm = 0;
    return;
  }
  A(pln->psi == 0 && pln->u == 0 && pln->perm == 0);
  if (pln->sorted) {
    const INT nn[3] = {pln->n0, pln->n1, pln->n2};
    pln->perm = (INT *)Y(malloc)(M * sizeof(INT));
    Y(conv_sort_nodes)(pln->x, nn, 3, pln->M, pln->perm);
  }
  pln->u = (INT *)Y(malloc)(M * 3 * sizeof(INT));
  pln->psi = (R *)Y(malloc)(M * 3 * len * sizeof(R));
  fill_u(pln);
  if (wakefulness == PLNR_AWAKE)
    fill_psi(pln);
  else
    memset(pln->psi, 0, M * 3 * len * sizeof(R));
}

/* Forward B (g -> f) / adjoint B^H (f -> g, scatter-add). Each axis splits into
 * at most two contiguous runs, so the tap nest is rectangular. */
static void run(const conv_3d_plan *pln, const problem_conv *pc, int forward)
{
  const INT n0 = pln->n0, n1 = pln->n1, n2 = pln->n2, M = pln->M;
  const INT len = 2 * (INT)pln->m + 2;
  const INT *perm = pln->perm;
  C *g = pc->g;
  C *f = pc->f;
  INT k;
  if (!forward)
    /* The scatter accumulates (+=) into an overlapping, node-dependent set that
     * does not cover the grid, so the whole grid must start zeroed. */
    memset(g, 0, (size_t)(n0 * n1 * n2) * sizeof(C));
  /* Rows stream in order; only the f access goes through perm. */
  for (k = 0; k < M; k++) {
    const INT j = perm ? perm[k] : k;
    C r =
         Y(conv_node_3d)(g, pln->psi + (k * 3 + 0) * len,
                       pln->psi + (k * 3 + 1) * len,
                       pln->psi + (k * 3 + 2) * len, pln->u[k * 3 + 0],
                       pln->u[k * 3 + 1], pln->u[k * 3 + 2], n0, n1, n2, len,
                       0, n0, forward ? K(0.0) : f[j], forward);
    if (forward)
      f[j] = r;
  }
}

static void apply(const plan *ego_, const problem *p)
{
  run((const conv_3d_plan *)ego_, (const problem_conv *)p, 1);
}

static void apply_adjoint(const plan *ego_, const problem *p)
{
  run((const conv_3d_plan *)ego_, (const problem_conv *)p, 0);
}

static void print(const plan *ego_, printer *pr)
{
  const conv_3d_plan *pln = (const conv_3d_plan *)ego_;
  pr->print(pr, "(conv_solver_3d%s pcost=%D)", pln->sorted ? " sorted" : "",
            (INT)pln->super.pcost);
}
static void destroy(plan *ego_)
{
  (void)ego_; /* Y(plan_destroy) sleeps the plan first; x/g/f are borrowed */
}
static const plan_adt conv_3d_plan_adt = {apply, awake, print, destroy,
                                          apply_adjoint};

/* d == 3 only */
static plan *mkplan_conv_3d(const solver *ego_, const problem *p, planner *pl)
{
  const conv_solver *ego = (const conv_solver *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  conv_3d_plan *pln;
  if (p->adt->kind != NFFT_PROBLEM_CONV)
    return 0;
  if (NO_NONTHREADEDP(pl))
    return 0; /* serial: prefer a threaded solver */
  if (pc->sz->rnk != 3)
    return 0;
  if (pc->window < NFFT_WINDOW_KAISER_BESSEL
      || pc->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */
  if (ego->sorted && pc->x_sorted)
    return 0; /* already in key order: nothing to sort */

  pln = (conv_3d_plan *)Y(
       plan_create)(sizeof(conv_3d_plan), &conv_3d_plan_adt);
  pln->n0 = Y(problem_conv_n)(p, 0);
  pln->n1 = Y(problem_conv_n)(p, 1);
  pln->n2 = Y(problem_conv_n)(p, 2);
  pln->N0 = Y(problem_conv_N)(p, 0);
  pln->N1 = Y(problem_conv_N)(p, 1);
  pln->N2 = Y(problem_conv_N)(p, 2);
  pln->M = pc->M;
  pln->m = pc->m;
  pln->window = pc->window;
  pln->analytic = pc->analytic;
  pln->x = pc->x; /* borrowed */
  pln->psi = 0;
  pln->u = 0;
  pln->perm = 0;
  pln->sorted = ego->sorted;

  /* Forward counts: taps = M * len^3, one complex-by-real accumulate per tap
   * (two loads: g, psi); other adds the per-axis conv_runs bookkeeping.
   * Estimate, not a full trace: the per-row axis-combine step (the hoisted
   * inner-axis run reused across the two outer axes) is a lower-order term
   * this does not count. */
  Y(ops_zero)(&pln->super.ops);
  {
    const double taps = (double)pln->M * pow((double)(2 * pln->m + 2), 3.0);
    pln->super.ops.mul = 2.0 * taps;
    pln->super.ops.add = 2.0 * taps;
    pln->super.ops.other = 2.0 * taps + 2.0 * (double)pln->M * 3.0;
  }
  Y(conv_order_ops)(&pln->super.ops, pln->M, 3, pln->m,
                 pln->n0 * pln->n1 * pln->n2, ego->sorted, pc->x_sorted);
  return &pln->super;
}

static const solver_adt conv_3d_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_3d};
void Y(conv_solver_3d_register)(planner *pl)
{
  Y(conv_register_both)(pl, &conv_3d_adt);
}
