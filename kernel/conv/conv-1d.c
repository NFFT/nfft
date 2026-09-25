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

/* 1D CONV solver: step C of the fast NFFT decomposition, the node convolution
 * (matrix B). Forward sums the oversampled grid g against the window psi at
 * each nonequispaced node x_j; the adjoint scatter-adds f onto g with the same
 * psi weights. psi and the wrapped window start u depend on x/window/n/N/m and
 * are built at awake, so apply does no window evaluation and no FLOOR/LRINT. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "conv.h"

typedef struct {
  plan super;
  INT n, N, M; /* geometry captured at mkplan */
  int m, window, analytic;
  const R *x; /* borrowed alias of the problem's nodes */
  R *psi; /* length M*(2m+2): psi[j*(2m+2)+lj]; NULL while SLEEPY */
  INT *u; /* length M: wrapped window start per node; NULL while SLEEPY */
  INT *perm; /* length M, node of row k in key order; NULL while SLEEPY and
              * for caller order */
  int sorted; /* this plan builds perm */
} conv_plan;

static void fill_u(conv_plan *pln)
{
  Y(conv_fill_u)(pln->x, pln->perm, &pln->n, 1, pln->m, pln->M, pln->u);
}

static void fill_psi(conv_plan *pln)
{
  Y(conv_fill_psi)(pln->window, pln->analytic, pln->x, pln->perm, &pln->n,
                &pln->N, 1, pln->m, pln->M, pln->psi);
}

/* The tables live between a wake and the return to SLEEPY (ADR-0009). Node
 * order and window starts are indices: they steer which grid cells apply
 * touches, so they hold true values at AWAKE_ZERO too, and the measured race
 * sees the locality it is choosing between. psi costs a window evaluation per
 * tap and changes no address, so AWAKE_ZERO gets zeros. */
static void awake(plan *ego_, int wakefulness)
{
  conv_plan *pln = (conv_plan *)ego_;
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
    pln->perm = (INT *)Y(malloc)(M * sizeof(INT));
    Y(conv_sort_nodes)(pln->x, &pln->n, 1, pln->M, pln->perm);
  }
  pln->u = (INT *)Y(malloc)(M * sizeof(INT));
  pln->psi = (R *)Y(malloc)(M * len * sizeof(R));
  fill_u(pln);
  if (wakefulness == PLNR_AWAKE)
    fill_psi(pln);
  else
    memset(pln->psi, 0, M * len * sizeof(R));
}

static void apply(const plan *ego_, const problem *p)
{
  const conv_plan *pln = (const conv_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  const INT n = pln->n, M = pln->M, len = 2 * (INT)pln->m + 2;
  const INT *perm = pln->perm;
  const C *g = pc->g;
  C *f = pc->f;
  INT k;
  /* Rows stream in order; only the f write goes through perm. */
  for (k = 0; k < M; k++)
    f[perm ? perm[k] : k] =
         Y(conv_node_1d)(g, pln->psi + k * len, pln->u[k], n,
                                         len, 0, n, K(0.0), 1);
}

static void apply_adjoint(const plan *ego_, const problem *p)
{
  const conv_plan *pln = (const conv_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  const INT n = pln->n, M = pln->M, len = 2 * (INT)pln->m + 2;
  const INT *perm = pln->perm;
  const C *f = pc->f;
  C *g = pc->g;
  INT k;
  /* The scatter accumulates (+=) into an overlapping, node-dependent set that
   * does not cover the grid, so the whole grid must start zeroed. */
  memset(g, 0, (size_t)n * sizeof(C));
  for (k = 0; k < M; k++)
    (void)Y(conv_node_1d)(g, pln->psi + k * len, pln->u[k], n, len, 0, n,
                       f[perm ? perm[k] : k], 0);
}

static void print(const plan *ego_, printer *pr)
{
  const conv_plan *pln = (const conv_plan *)ego_;
  pr->print(pr, "(conv_solver_1d%s pcost=%D)", pln->sorted ? " sorted" : "",
            (INT)pln->super.pcost);
}
static void destroy(plan *ego_)
{
  (void)ego_; /* Y(plan_destroy) sleeps the plan first; x/g/f are borrowed */
}
static const plan_adt conv_plan_adt = {apply, awake, print, destroy,
                                       apply_adjoint};

/* d == 1 only */
static plan *mkplan_conv_1d(const solver *ego_, const problem *p, planner *pl)
{
  const conv_solver *ego = (const conv_solver *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  conv_plan *pln;
  if (p->adt->kind != NFFT_PROBLEM_CONV)
    return 0;
  if (NO_NONTHREADEDP(pl))
    return 0; /* serial: prefer a threaded solver */
  if (pc->sz->rnk != 1)
    return 0;
  if (pc->window < NFFT_WINDOW_KAISER_BESSEL
      || pc->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */
  if (ego->sorted && pc->x_sorted)
    return 0; /* already in key order: nothing to sort */

  pln = (conv_plan *)Y(plan_create)(sizeof(conv_plan), &conv_plan_adt);
  pln->n = Y(problem_conv_n)(p, 0);
  pln->N = Y(problem_conv_N)(p, 0);
  pln->M = pc->M;
  pln->m = pc->m;
  pln->window = pc->window;
  pln->analytic = pc->analytic;
  pln->x = pc->x; /* borrowed */
  pln->psi = 0;
  pln->u = 0;
  pln->perm = 0;
  pln->sorted = ego->sorted;

  /* Forward counts: taps = M * len^1, one complex-by-real accumulate per tap
   * (two loads: g, psi); other adds the per-axis conv_runs bookkeeping. */
  Y(ops_zero)(&pln->super.ops);
  {
    const double taps = (double)pln->M * (double)(2 * pln->m + 2);
    pln->super.ops.mul = 2.0 * taps;
    pln->super.ops.add = 2.0 * taps;
    pln->super.ops.other = 2.0 * taps + 2.0 * (double)pln->M;
  }
  Y(conv_order_ops)(&pln->super.ops, pln->M, 1, pln->m, pln->n, ego->sorted,
                 pc->x_sorted);
  return &pln->super;
}

static const solver_adt conv_1d_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_1d};
void Y(conv_solver_1d_register)(planner *pl)
{
  Y(conv_register_both)(pl, &conv_1d_adt);
}
