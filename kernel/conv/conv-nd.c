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

/* nD CONV solver: step C of the fast NFFT decomposition, the node convolution
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
  int d; /* rank, >= 4 */
  INT *n, *N; /* owned, length d: geometry captured at mkplan */
  INT M; /* node count */
  INT ntot; /* owned product of n[], captured at mkplan */
  int m, window;
  const R *x; /* borrowed alias of the problem's nodes, length d*M */
  R *psi; /* length M*d*(2m+2): psi[(j*d+t)*(2m+2)+lj]; NULL while SLEEPY */
  INT *u; /* length M*d: wrapped window start, u[j*d+t]; NULL while SLEEPY */
  INT *perm; /* length M, node of row k in key order; NULL while SLEEPY and
              * for caller order */
  int sorted; /* this plan builds perm */
} conv_nd_plan;

static void fill_u(conv_nd_plan *pln)
{
  Y(conv_fill_u)(pln->x, pln->perm, pln->n, pln->d, pln->m, pln->M, pln->u);
}

static void fill_psi(conv_nd_plan *pln)
{
  Y(conv_fill_psi)(pln->window, pln->x, pln->perm, pln->n, pln->N, pln->d, pln->m,
                pln->M, pln->psi);
}

/* The tables live between a wake and the return to SLEEPY (ADR-0009). Node
 * order and window starts are indices: they steer which grid cells apply
 * touches, so they hold true values at AWAKE_ZERO too, and the measured race
 * sees the locality it is choosing between. psi costs a window evaluation per
 * tap and changes no address, so AWAKE_ZERO gets zeros. */
static void awake(plan *ego_, int wakefulness)
{
  conv_nd_plan *pln = (conv_nd_plan *)ego_;
  const size_t M = (size_t)pln->M, d = (size_t)pln->d;
  const size_t len = (size_t)(2 * pln->m + 2);
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
    Y(conv_sort_nodes)(pln->x, pln->n, pln->d, pln->M, pln->perm);
  }
  pln->u = (INT *)Y(malloc)(M * d * sizeof(INT));
  pln->psi = (R *)Y(malloc)(M * d * len * sizeof(R));
  fill_u(pln);
  if (wakefulness == PLNR_AWAKE)
    fill_psi(pln);
  else
    memset(pln->psi, 0, M * d * len * sizeof(R));
}

/* Forward B (g -> f) / adjoint B^H (f -> g, scatter-add). */
static void run(const conv_nd_plan *pln, const problem_conv *pc, int forward)
{
  const int d = pln->d;
  const INT *n = pln->n;
  const INT M = pln->M;
  const int m = pln->m;
  const R *psi = pln->psi;
  C *f, *g;
  INT k;

  if (forward) {
    f = pc->f;
    g = (C *)pc->g;
  } else {
    f = (C *)pc->f;
    g = pc->g;
    /* The scatter accumulates (+=) into an overlapping, node-dependent set that
     * does not cover the grid, so the whole grid must start zeroed. */
    memset(g, 0, (size_t)pln->ntot * sizeof(C));
  }

  /* Rows stream in order; only the f access goes through perm. */
  for (k = 0; k < M; k++) {
    const INT j = pln->perm ? pln->perm[k] : k;
    C r = Y(
         conv_node_nd)(g, pln->psi + k * d * (2 * m + 2), pln->u + k * d, n, d,
                       2 * m + 2, 0, n[0], forward ? K(0.0) : f[j], forward);
    if (forward)
      f[j] = r;
  }
}

static void apply(const plan *ego_, const problem *p)
{
  run((const conv_nd_plan *)ego_, (const problem_conv *)p, 1);
}

static void apply_adjoint(const plan *ego_, const problem *p)
{
  run((const conv_nd_plan *)ego_, (const problem_conv *)p, 0);
}

static void print(const plan *ego_, printer *pr)
{
  const conv_nd_plan *pln = (const conv_nd_plan *)ego_;
  pr->print(pr, "(conv_solver_nd%s pcost=%D)", pln->sorted ? " sorted" : "",
            (INT)pln->super.pcost);
}
static void destroy(plan *ego_)
{
  conv_nd_plan *pln = (conv_nd_plan *)ego_;
  Y(free)(pln->N);
  Y(free)(pln->n);
  /* Y(plan_destroy) sleeps the plan first; x/g/f are borrowed */
}
static const plan_adt conv_nd_plan_adt = {apply, awake, print, destroy,
                                          apply_adjoint};

/* d >= 4 only. */
static plan *mkplan_conv_nd(const solver *ego_, const problem *p, planner *pl)
{
  const conv_solver *ego = (const conv_solver *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  conv_nd_plan *pln;
  int t, d;
  if (p->adt->kind != NFFT_PROBLEM_CONV)
    return 0;
  if (NO_NONTHREADEDP(pl))
    return 0; /* serial: prefer a threaded solver */
  if (pc->sz->rnk < 4)
    return 0;
  if (pc->window < NFFT_WINDOW_KAISER_BESSEL
      || pc->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */
  if (ego->sorted && pc->x_sorted)
    return 0; /* already in key order: nothing to sort */

  d = pc->sz->rnk;
  pln = (conv_nd_plan *)Y(
       plan_create)(sizeof(conv_nd_plan), &conv_nd_plan_adt);
  pln->d = d;
  pln->n = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->N = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->ntot = Y(problem_conv_ntot)(p);
  for (t = 0; t < d; t++) {
    pln->n[t] = Y(problem_conv_n)(p, t);
    pln->N[t] = Y(problem_conv_N)(p, t);
  }
  pln->M = pc->M;
  pln->m = pc->m;
  /* run holds stack VLAs sized d*(2m+2); bound the budget so a pathological
   * rank/m cannot overflow the stack. */
  A((size_t)d * (size_t)(2 * pln->m + 2) * sizeof(INT) <= (size_t)(64 * 1024));
  pln->window = pc->window;
  pln->x = pc->x; /* borrowed */
  pln->psi = 0;
  pln->u = 0;
  pln->perm = 0;
  pln->sorted = ego->sorted;

  /* Forward counts: taps = M * len^d, one complex-by-real accumulate per tap
   * (two loads: g, psi); other adds the per-axis wrapped-index setup.
   * Estimate, not a full trace: the odometer's cached prefix-product recompute
   * (run()'s phi_prod/ll_plain cascade) is a lower-order term this does not
   * count. */
  Y(ops_zero)(&pln->super.ops);
  {
    const double taps =
         (double)pln->M * pow((double)(2 * pln->m + 2), (double)d);
    pln->super.ops.mul = 2.0 * taps;
    pln->super.ops.add = 2.0 * taps;
    pln->super.ops.other = 2.0 * taps + 2.0 * (double)pln->M * (double)d;
  }
  Y(conv_order_ops)(&pln->super.ops, pln->M, d, pln->m, pln->ntot, ego->sorted,
                 pc->x_sorted);
  return &pln->super;
}

static const solver_adt conv_nd_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_nd};
void Y(conv_solver_nd_register)(planner *pl)
{
  Y(conv_register_both)(pl, &conv_nd_adt);
}
