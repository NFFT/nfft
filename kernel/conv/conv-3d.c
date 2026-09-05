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

typedef struct {
  plan super;
  INT n0, n1, n2, N0, N1, N2, M; /* geometry captured at mkplan */
  int m, window;
  const R *x; /* borrowed alias of the problem's nodes */
  R *psi; /* length M*3*(2m+2): psi[(j*3+t)*(2m+2)+lj] */
  INT *u; /* length M*3: wrapped window start, u[j*3+t] */
  int level; /* content of psi/u: SLEEPY (stale), AWAKE_ZERO or AWAKE */
} conv_3d_plan;

/* Window starts are indices: they steer which grid cells apply touches, so the
 * true values are built for AWAKE_ZERO as well. Only the psi weights, which
 * cost a window evaluation each and change no address, get placeholders. */
static void fill_u(conv_3d_plan *pln)
{
  const INT nn[3] = {pln->n0, pln->n1, pln->n2};
  const INT M = pln->M;
  const int m = pln->m;
  INT j;
  int t;
  for (t = 0; t < 3; t++)
    for (j = 0; j < M; j++) {
      INT c = LRINT(FLOOR(pln->x[j * 3 + t] * (R)nn[t]));
      pln->u[j * 3 + t] = (((c - m) % nn[t]) + nn[t]) % nn[t];
    }
}

static void fill_psi(conv_3d_plan *pln)
{
  const INT nn[3] = {pln->n0, pln->n1, pln->n2};
  const INT NN[3] = {pln->N0, pln->N1, pln->N2};
  const int m = pln->m;
  int t;
  for (t = 0; t < 3; t++)
    Y(window_phi_precompute)(pln->window, nn[t], NN[t], m, pln->x + t, 3, pln->M,
                          pln->psi + t * (2 * m + 2), 3 * (2 * m + 2));
}

/* AWAKE_ZERO must cost no window evaluation, so the tables get placeholder
 * zeros: u == 0 keeps every apply index in range and psi == 0 keeps every
 * apply flop finite. AWAKE_ZERO reached by downgrade only drops the level, so
 * a later upgrade refills. */
static void awake(plan *ego_, int wakefulness)
{
  conv_3d_plan *pln = (conv_3d_plan *)ego_;
  if (wakefulness >= PLNR_AWAKE_ZERO && pln->level == PLNR_SLEEPY)
    fill_u(pln);
  if (wakefulness == PLNR_AWAKE) {
    if (pln->level != PLNR_AWAKE)
      fill_psi(pln);
  } else if (wakefulness == PLNR_AWAKE_ZERO && pln->level == PLNR_SLEEPY) {
    memset(pln->psi, 0,
           (size_t)pln->M * 3 * (size_t)(2 * pln->m + 2) * sizeof(R));
  }
  pln->level = wakefulness;
}

/* Forward B (g -> f) / adjoint B^H (f -> g, scatter-add). Each axis splits into
 * at most two contiguous runs, so the tap nest is rectangular. */
static void run(const conv_3d_plan *pln, const problem_conv *pc, int forward)
{
  const INT n0 = pln->n0, n1 = pln->n1, n2 = pln->n2, M = pln->M;
  const INT len = 2 * (INT)pln->m + 2;
  C *g = pc->g;
  C *f = pc->f;
  INT j;
  if (!forward)
    /* The scatter accumulates (+=) into an overlapping, node-dependent set that
     * does not cover the grid, so the whole grid must start zeroed. */
    memset(g, 0, (size_t)(n0 * n1 * n2) * sizeof(C));
  for (j = 0; j < M; j++) {
    const R *psi0 = pln->psi + (j * 3 + 0) * len;
    const R *psi1 = pln->psi + (j * 3 + 1) * len;
    const R *psi2 = pln->psi + (j * 3 + 2) * len;
    INT tof0[2], gof0[2], rl0[2], tof1[2], gof1[2], rl1[2];
    INT tof2[2], gof2[2], rl2[2], i, k, l;
    C acc = K(0.0);
    C fj = forward ? K(0.0) : f[j];
    int a, b;
    Y(conv_runs)(pln->u[j * 3 + 0], n0, len, tof0, gof0, rl0);
    Y(conv_runs)(pln->u[j * 3 + 1], n1, len, tof1, gof1, rl1);
    Y(conv_runs)(pln->u[j * 3 + 2], n2, len, tof2, gof2, rl2);
    { /* the inner axis' runs do not vary over the outer axes; hoist them */
      const INT la = rl2[0], lb = rl2[1], goa = gof2[0], gob = gof2[1];
      const R *p2a = psi2 + tof2[0], *p2b = psi2 + tof2[1];
      for (a = 0; a < 2; a++)
        for (i = 0; i < rl0[a]; i++) {
          const R p0 = psi0[tof0[a] + i];
          C *gplane = g + (gof0[a] + i) * n1 * n2;
          for (b = 0; b < 2; b++) {
            const INT nk = rl1[b], to1 = tof1[b], go1 = gof1[b];
            for (k = 0; k < nk; k++) {
              const R p01 = p0 * psi1[to1 + k];
              C *grow = gplane + (go1 + k) * n2;
              C *ga = grow + goa, *gb = grow + gob;
              if (forward) {
                C sub = K(0.0);
                for (l = 0; l < la; l++)
                  sub += ga[l] * p2a[l];
                for (l = 0; l < lb; l++)
                  sub += gb[l] * p2b[l];
                acc += sub * p01;
              } else {
                const C fp = fj * p01;
                for (l = 0; l < la; l++)
                  ga[l] += fp * p2a[l];
                for (l = 0; l < lb; l++)
                  gb[l] += fp * p2b[l];
              }
            }
          }
        }
    }
    if (forward)
      f[j] = acc;
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
  pr->print(pr, "(conv_solver_3d pcost=%D)", (INT)pln->super.pcost);
}
static void destroy(plan *ego_)
{
  conv_3d_plan *pln = (conv_3d_plan *)ego_;
  Y(free)(pln->psi);
  Y(free)(pln->u);
  /* x/g/f are borrowed caller arrays. */
}
static const plan_adt conv_3d_plan_adt = {apply, awake, print, destroy,
                                          apply_adjoint};

/* d == 3 only */
static plan *mkplan_conv_3d(const solver *ego, const problem *p, planner *pl)
{
  const problem_conv *pc = (const problem_conv *)p;
  conv_3d_plan *pln;
  (void)ego;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_CONV)
    return 0;
  if (pc->sz->rnk != 3)
    return 0;
  if (pc->window < NFFT_WINDOW_KAISER_BESSEL
      || pc->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */

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
  pln->x = pc->x; /* borrowed */
  pln->psi =
       (R *)Y(
            malloc)((size_t)pln->M * 3 * (size_t)(2 * pln->m + 2) * sizeof(R));
  pln->u = (INT *)Y(malloc)((size_t)pln->M * 3 * sizeof(INT));
  pln->level = PLNR_SLEEPY;
  pln->super.pcost = Y(conv_b_pcost)(p);
  return &pln->super;
}

static const solver_adt conv_3d_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_3d};
void Y(conv_solver_3d_register)(planner *pl)
{
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &conv_3d_adt));
}
