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

/* 2D CONV solver: step C of the fast NFFT decomposition, the node convolution
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
  INT n0, n1, N0, N1, M; /* geometry captured at mkplan */
  int m, window;
  const R *x; /* borrowed alias of the problem's nodes */
  R *psi; /* length M*2*(2m+2): psi[(j*2+t)*(2m+2)+lj] */
  INT *u; /* length M*2: wrapped window start, u[j*2+t] */
  int level; /* content of psi/u: SLEEPY (stale), AWAKE_ZERO or AWAKE */
} conv_2d_plan;

/* Window starts are indices: they steer which grid cells apply touches, so the
 * true values are built for AWAKE_ZERO as well. Only the psi weights, which
 * cost a window evaluation each and change no address, get placeholders. */
static void fill_u(conv_2d_plan *pln)
{
  const INT nn[2] = {pln->n0, pln->n1};
  const INT M = pln->M;
  const int m = pln->m;
  INT j;
  int t;
  for (t = 0; t < 2; t++)
    for (j = 0; j < M; j++) {
      INT c = LRINT(FLOOR(pln->x[j * 2 + t] * (R)nn[t]));
      pln->u[j * 2 + t] = (((c - m) % nn[t]) + nn[t]) % nn[t];
    }
}

static void fill_psi(conv_2d_plan *pln)
{
  const INT nn[2] = {pln->n0, pln->n1};
  const INT NN[2] = {pln->N0, pln->N1};
  const int m = pln->m;
  int t;
  for (t = 0; t < 2; t++)
    Y(window_phi_precompute)(pln->window, nn[t], NN[t], m, pln->x + t, 2, pln->M,
                          pln->psi + t * (2 * m + 2), 2 * (2 * m + 2));
}

/* AWAKE_ZERO must cost no window evaluation, so the tables get placeholder
 * zeros: u == 0 keeps every apply index in range and psi == 0 keeps every
 * apply flop finite. AWAKE_ZERO reached by downgrade only drops the level, so
 * a later upgrade refills. */
static void awake(plan *ego_, int wakefulness)
{
  conv_2d_plan *pln = (conv_2d_plan *)ego_;
  if (wakefulness >= PLNR_AWAKE_ZERO && pln->level == PLNR_SLEEPY)
    fill_u(pln);
  if (wakefulness == PLNR_AWAKE) {
    if (pln->level != PLNR_AWAKE)
      fill_psi(pln);
  } else if (wakefulness == PLNR_AWAKE_ZERO && pln->level == PLNR_SLEEPY) {
    memset(pln->psi, 0,
           (size_t)pln->M * 2 * (size_t)(2 * pln->m + 2) * sizeof(R));
  }
  pln->level = wakefulness;
}

/* Forward B (g -> f) / adjoint B^H (f -> g, scatter-add). Each axis splits into
 * at most two contiguous runs, so the tap nest is rectangular. */
static void run(const conv_2d_plan *pln, const problem_conv *pc, int forward)
{
  const INT n0 = pln->n0, n1 = pln->n1, M = pln->M;
  const INT len = 2 * (INT)pln->m + 2;
  C *g = pc->g;
  C *f = pc->f;
  INT j;
  if (!forward)
    /* The scatter accumulates (+=) into an overlapping, node-dependent set that
     * does not cover the grid, so the whole grid must start zeroed. */
    memset(g, 0, (size_t)(n0 * n1) * sizeof(C));
  for (j = 0; j < M; j++) {
    const R *psi0 = pln->psi + (j * 2 + 0) * len;
    const R *psi1 = pln->psi + (j * 2 + 1) * len;
    INT tof0[2], gof0[2], rl0[2], tof1[2], gof1[2], rl1[2], i, k;
    C acc = K(0.0);
    C fj = forward ? K(0.0) : f[j];
    int a;
    Y(conv_runs)(pln->u[j * 2 + 0], n0, len, tof0, gof0, rl0);
    Y(conv_runs)(pln->u[j * 2 + 1], n1, len, tof1, gof1, rl1);
    { /* the inner axis' runs do not vary over the outer axes; hoist them */
      const INT ka = rl1[0], kb = rl1[1], goa = gof1[0], gob = gof1[1];
      const R *p1a = psi1 + tof1[0], *p1b = psi1 + tof1[1];
      for (a = 0; a < 2; a++)
        for (i = 0; i < rl0[a]; i++) {
          const R p0 = psi0[tof0[a] + i];
          C *grow = g + (gof0[a] + i) * n1;
          C *ga = grow + goa, *gb = grow + gob;
          if (forward) {
            C sub = K(0.0);
            for (k = 0; k < ka; k++)
              sub += ga[k] * p1a[k];
            for (k = 0; k < kb; k++)
              sub += gb[k] * p1b[k];
            acc += sub * p0;
          } else {
            const C fp = fj * p0;
            for (k = 0; k < ka; k++)
              ga[k] += fp * p1a[k];
            for (k = 0; k < kb; k++)
              gb[k] += fp * p1b[k];
          }
        }
    }
    if (forward)
      f[j] = acc;
  }
}

static void apply(const plan *ego_, const problem *p)
{
  run((const conv_2d_plan *)ego_, (const problem_conv *)p, 1);
}

static void apply_adjoint(const plan *ego_, const problem *p)
{
  run((const conv_2d_plan *)ego_, (const problem_conv *)p, 0);
}

static void print(const plan *ego_, printer *pr)
{
  const conv_2d_plan *pln = (const conv_2d_plan *)ego_;
  pr->print(pr, "(conv_solver_2d pcost=%D)", (INT)pln->super.pcost);
}
static void destroy(plan *ego_)
{
  conv_2d_plan *pln = (conv_2d_plan *)ego_;
  Y(free)(pln->psi);
  Y(free)(pln->u);
  /* x/g/f are borrowed caller arrays. */
}
static const plan_adt conv_2d_plan_adt = {apply, awake, print, destroy,
                                          apply_adjoint};

/* d == 2 only */
static plan *mkplan_conv_2d(const solver *ego, const problem *p, planner *pl)
{
  const problem_conv *pc = (const problem_conv *)p;
  conv_2d_plan *pln;
  (void)ego;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_CONV)
    return 0;
  if (pc->sz->rnk != 2)
    return 0;
  if (pc->window < NFFT_WINDOW_KAISER_BESSEL
      || pc->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */

  pln = (conv_2d_plan *)Y(
       plan_create)(sizeof(conv_2d_plan), &conv_2d_plan_adt);
  pln->n0 = Y(problem_conv_n)(p, 0);
  pln->n1 = Y(problem_conv_n)(p, 1);
  pln->N0 = Y(problem_conv_N)(p, 0);
  pln->N1 = Y(problem_conv_N)(p, 1);
  pln->M = pc->M;
  pln->m = pc->m;
  pln->window = pc->window;
  pln->x = pc->x; /* borrowed */
  pln->psi =
       (R *)Y(
            malloc)((size_t)pln->M * 2 * (size_t)(2 * pln->m + 2) * sizeof(R));
  pln->u = (INT *)Y(malloc)((size_t)pln->M * 2 * sizeof(INT));
  pln->level = PLNR_SLEEPY;
  pln->super.pcost = Y(conv_b_pcost)(p);
  return &pln->super;
}

static const solver_adt conv_2d_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_2d};
void Y(conv_solver_2d_register)(planner *pl)
{
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &conv_2d_adt));
}
