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

typedef struct
{
  plan super;
  INT n, N, M; /* geometry captured at mkplan */
  int m, window;
  const R *x; /* borrowed alias of the problem's nodes */
  R *psi;     /* length M*(2m+2): psi[j*(2m+2)+lj] */
  INT *u;     /* length M: wrapped window start per node */
  int level;  /* content of psi/u: SLEEPY (stale), AWAKE_ZERO or AWAKE */
} conv_plan;

static void fill(conv_plan *pln) {
  INT n = pln->n, M = pln->M, j;
  int m = pln->m;
  Y(window_phi_precompute)
  (pln->window, n, pln->N, m, pln->x, 1, M, pln->psi, 2 * m + 2);
  for (j = 0; j < M; j++) {
    INT c = LRINT(FLOOR((R)n * pln->x[j]));
    pln->u[j] = (((c - m) % n) + n) % n;
  }
}

/* AWAKE_ZERO must cost no window evaluation, so the tables get placeholder
 * zeros: u == 0 keeps every apply index in range and psi == 0 keeps every
 * apply flop finite. AWAKE_ZERO reached by downgrade only drops the level, so
 * a later upgrade refills. */
static void awake(plan *ego_, int wakefulness) {
  conv_plan *pln = (conv_plan *)ego_;
  if (wakefulness == PLNR_AWAKE) {
    if (pln->level != PLNR_AWAKE)
      fill(pln);
  } else if (wakefulness == PLNR_AWAKE_ZERO && pln->level == PLNR_SLEEPY) {
    memset(pln->psi, 0, (size_t)pln->M * (size_t)(2 * pln->m + 2) * sizeof(R));
    Y(conv_spread_u)(pln->u, pln->M, 1, &pln->n);
  }
  pln->level = wakefulness;
}

static void apply(const plan *ego_, const problem *p) {
  const conv_plan *pln = (const conv_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  const INT n = pln->n, M = pln->M, len = 2 * (INT)pln->m + 2;
  const C *g = pc->g;
  C *f = pc->f;
  INT j;
  for (j = 0; j < M; j++) {
    const R *psij = pln->psi + j * len;
    INT tof[2], gof[2], rl[2], k;
    C acc = K(0.0);
    int r;
    Y(conv_runs)(pln->u[j], n, len, tof, gof, rl);
    for (r = 0; r < 2; r++) {
      const R *ps = psij + tof[r];
      const C *gr = g + gof[r];
      for (k = 0; k < rl[r]; k++)
        acc += gr[k] * ps[k];
    }
    f[j] = acc;
  }
}

static void apply_adjoint(const plan *ego_, const problem *p) {
  const conv_plan *pln = (const conv_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  const INT n = pln->n, M = pln->M, len = 2 * (INT)pln->m + 2;
  const C *f = pc->f;
  C *g = pc->g;
  INT j;
  /* The scatter accumulates (+=) into an overlapping, node-dependent set that
   * does not cover the grid, so the whole grid must start zeroed. */
  memset(g, 0, (size_t)n * sizeof(C));
  for (j = 0; j < M; j++) {
    const R *psij = pln->psi + j * len;
    C fj = f[j];
    INT tof[2], gof[2], rl[2], k;
    int r;
    Y(conv_runs)(pln->u[j], n, len, tof, gof, rl);
    for (r = 0; r < 2; r++) {
      const R *ps = psij + tof[r];
      C *gr = g + gof[r];
      for (k = 0; k < rl[r]; k++)
        gr[k] += fj * ps[k];
    }
  }
}

static void print(const plan *ego_, printer *pr) {
  const conv_plan *pln = (const conv_plan *)ego_;
  pr->print(pr, "(conv_solver_1d pcost=%D)", (INT)pln->super.pcost);
}
static void destroy(plan *ego_) {
  conv_plan *pln = (conv_plan *)ego_;
  Y(free)
  (pln->psi);
  Y(free)
  (pln->u);
  /* x/g/f are borrowed caller arrays. */
}
static const plan_adt conv_plan_adt = {apply, awake, print, destroy,
                                       apply_adjoint};

/* d == 1 only */
static plan *mkplan_conv_1d(const solver *ego, const problem *p, planner *pl) {
  const problem_conv *pc = (const problem_conv *)p;
  conv_plan *pln;
  (void)ego;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_CONV)
    return 0;
  if (pc->sz->rnk != 1)
    return 0;
  if (pc->window < NFFT_WINDOW_KAISER_BESSEL ||
      pc->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */

  pln = (conv_plan *)Y(plan_create)(sizeof(conv_plan), &conv_plan_adt);
  pln->n = Y(problem_conv_n)(p, 0);
  pln->N = Y(problem_conv_N)(p, 0);
  pln->M = pc->M;
  pln->m = pc->m;
  pln->window = pc->window;
  pln->x = pc->x; /* borrowed */
  pln->psi = (R *)Y(malloc)((size_t)pln->M * (size_t)(2 * pln->m + 2) * sizeof(R));
  pln->u = (INT *)Y(malloc)((size_t)pln->M * sizeof(INT));
  pln->level = PLNR_SLEEPY;
  pln->super.pcost = Y(conv_b_pcost)(p);
  return &pln->super;
}

static const solver_adt conv_1d_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_1d};
void Y(conv_solver_1d_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &conv_1d_adt));
}
