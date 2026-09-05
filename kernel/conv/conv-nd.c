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

typedef struct {
  plan super;
  int d; /* rank, >= 4 */
  INT *n, *N; /* owned, length d: geometry captured at mkplan */
  INT M; /* node count */
  INT ntot; /* owned product of n[], captured at mkplan */
  int m, window;
  const R *x; /* borrowed alias of the problem's nodes, length d*M */
  R *psi; /* length M*d*(2m+2): psi[(j*d+t)*(2m+2)+lj] */
  INT *u; /* length M*d: wrapped window start, u[j*d+t] */
  int level; /* content of psi/u: SLEEPY (stale), AWAKE_ZERO or AWAKE */
} conv_nd_plan;

static void fill(conv_nd_plan *pln)
{
  const int d = pln->d, m = pln->m;
  const INT M = pln->M;
  INT j;
  int t;
  for (t = 0; t < d; t++) {
    const INT nt = pln->n[t];
    Y(window_phi_precompute)(pln->window, nt, pln->N[t], m, pln->x + t, d, M,
                          pln->psi + t * (2 * m + 2), d * (2 * m + 2));
    for (j = 0; j < M; j++) {
      INT c = LRINT(FLOOR(pln->x[j * d + t] * (R)nt));
      pln->u[j * d + t] = (((c - m) % nt) + nt) % nt;
    }
  }
}

/* AWAKE_ZERO must cost no window evaluation, so the tables get placeholder
 * zeros: u == 0 keeps every apply index in range and psi == 0 keeps every
 * apply flop finite. AWAKE_ZERO reached by downgrade only drops the level, so
 * a later upgrade refills. */
static void awake(plan *ego_, int wakefulness)
{
  conv_nd_plan *pln = (conv_nd_plan *)ego_;
  if (wakefulness == PLNR_AWAKE) {
    if (pln->level != PLNR_AWAKE)
      fill(pln);
  } else if (wakefulness == PLNR_AWAKE_ZERO && pln->level == PLNR_SLEEPY) {
    memset(pln->psi, 0,
           (size_t)pln->M * (size_t)pln->d * (size_t)(2 * pln->m + 2)
                * sizeof(R));
    Y(conv_spread_u)(pln->u, pln->M, pln->d, pln->n);
  }
  pln->level = wakefulness;
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
  INT lprod;
  INT lj[d]; /* multi index over the taps, 0 <= lj <= 2m+1 */
  INT ll_plain[d + 1]; /* postfix plain index in g */
  R phi_prod[d + 1]; /* postfix product of psi */
  INT l_all[d * (2 * m + 2)]; /* wrapped grid indices per axis/tap */
  INT j, t, t2, l_L;

  if (forward) {
    f = pc->f;
    g = (C *)pc->g;
    memset(f, 0, (size_t)M * sizeof(C));
  } else {
    f = (C *)pc->f;
    g = pc->g;
    /* The scatter accumulates (+=) into an overlapping, node-dependent set that
     * does not cover the grid, so the whole grid must start zeroed. */
    memset(g, 0, (size_t)pln->ntot * sizeof(C));
  }

  for (t = 0, lprod = 1; t < d; t++)
    lprod *= (2 * m + 2);

  for (j = 0; j < M; j++) {
    for (t = d - 1; t >= 0; t--) {
      INT lj_t;
      for (lj_t = 0; lj_t < 2 * m + 2; lj_t++)
        l_all[t * (2 * m + 2) + lj_t] = (pln->u[j * d + t] + lj_t) % n[t];
      lj[t] = 0;
    }
    t++;

    phi_prod[0] = K(1.0);
    ll_plain[0] = 0;

    for (l_L = 0; l_L < lprod; l_L++) {
      for (t2 = t; t2 < d; t2++) {
        phi_prod[t2 + 1] =
             phi_prod[t2] * psi[(j * d + t2) * (2 * m + 2) + lj[t2]];
        ll_plain[t2 + 1] =
             ll_plain[t2] * n[t2] + l_all[t2 * (2 * m + 2) + lj[t2]];
      }

      if (forward)
        f[j] += phi_prod[d] * g[ll_plain[d]];
      else
        g[ll_plain[d]] += phi_prod[d] * f[j];

      for (t = d - 1; (t > 0) && (lj[t] == 2 * m + 1); t--)
        lj[t] = 0;
      lj[t]++;
    }
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
  pr->print(pr, "(conv_solver_nd pcost=%D)", (INT)pln->super.pcost);
}
static void destroy(plan *ego_)
{
  conv_nd_plan *pln = (conv_nd_plan *)ego_;
  Y(free)(pln->psi);
  Y(free)(pln->u);
  Y(free)(pln->N);
  Y(free)(pln->n);
  /* x/g/f are borrowed caller arrays. */
}
static const plan_adt conv_nd_plan_adt = {apply, awake, print, destroy,
                                          apply_adjoint};

/* d >= 4 only. */
static plan *mkplan_conv_nd(const solver *ego, const problem *p, planner *pl)
{
  const problem_conv *pc = (const problem_conv *)p;
  conv_nd_plan *pln;
  int t, d;
  (void)ego;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_CONV)
    return 0;
  if (pc->sz->rnk < 4)
    return 0;
  if (pc->window < NFFT_WINDOW_KAISER_BESSEL
      || pc->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */

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
  pln->psi = (R *)Y(
       malloc)((size_t)pln->M * (size_t)d * (size_t)(2 * pln->m + 2)
                         * sizeof(R));
  pln->u = (INT *)Y(malloc)((size_t)pln->M * (size_t)d * sizeof(INT));
  pln->level = PLNR_SLEEPY;
  pln->super.pcost = Y(conv_b_pcost)(p);
  return &pln->super;
}

static const solver_adt conv_nd_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_nd};
void Y(conv_solver_nd_register)(planner *pl)
{
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &conv_nd_adt));
}
