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

/* 1D DECONV solver: step A of the fast NFFT decomposition. Forward divides f_hat
 * by the window's phi_hut factors and zero-pads onto the oversampled grid g;
 * the adjoint gathers g -> f_hat through the same 1/phi_hut. phi_hut is
 * node-independent, so it is built at awake. Slot ks carries frequency
 * k = ks - Nneg, with Nneg = N/2 for type-I and N/2 - 1 for type-II; odd N
 * normalizes to type-I in mkproblem_deconv, so type-II implies even N. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "deconv.h"

typedef struct {
  plan super;
  INT n, N; /* geometry captured at mkplan */
  INT Nneg, Npos; /* slot split: k(ks) = ks - Nneg, Npos = N - Nneg */
  int m, window;
  R *phi_hut_inv; /* length N: 1/phi_hut(ks - Nneg) */
  int level; /* content of phi_hut_inv: SLEEPY (stale), AWAKE_ZERO or AWAKE */
} deconv_plan;

/* AWAKE_ZERO must cost no window evaluation, so the table gets placeholder
 * zeros. AWAKE_ZERO reached by downgrade only drops the level, so a later
 * upgrade refills. */
static void awake(plan *ego_, int wakefulness)
{
  deconv_plan *pln = (deconv_plan *)ego_;
  INT N = pln->N, ks;
  if (wakefulness == PLNR_AWAKE) {
    if (pln->level != PLNR_AWAKE) {
      Y(window_phi_hut_apply)(pln->window, pln->n, N, pln->m, -pln->Nneg,
                           pln->phi_hut_inv, N);
      for (ks = 0; ks < N; ks++)
        pln->phi_hut_inv[ks] = K(1.0) / pln->phi_hut_inv[ks];
    }
  } else if (wakefulness == PLNR_AWAKE_ZERO && pln->level == PLNR_SLEEPY)
    memset(pln->phi_hut_inv, 0, (size_t)N * sizeof(R));
  pln->level = wakefulness;
}

static void apply(const plan *ego_, const problem *p)
{
  const deconv_plan *pln = (const deconv_plan *)ego_;
  const problem_deconv *pd = (const problem_deconv *)p;
  INT N = pln->N, n = pln->n, Nneg = pln->Nneg, Npos = pln->Npos;
  const C *f_hat = pd->f_hat;
  C *g = pd->g;
  INT ks;
  /* The scatter below writes every touched cell exactly once, so they don't need
   * pre-zeroing. The untouched cells are the single contiguous run
   * [Npos .. n - Nneg), of length n - N. */
  memset(g + Npos, 0, (size_t)(n - N) * sizeof(C));
  for (ks = 0; ks < N; ks++) {
    INT pos = (ks < Nneg) ? n - Nneg + ks : ks - Nneg;
    g[pos] = f_hat[ks] * pln->phi_hut_inv[ks];
  }
}

static void apply_adjoint(const plan *ego_, const problem *p)
{
  const deconv_plan *pln = (const deconv_plan *)ego_;
  const problem_deconv *pd = (const problem_deconv *)p;
  INT N = pln->N, n = pln->n, Nneg = pln->Nneg;
  const C *g = pd->g;
  C *f_hat = pd->f_hat;
  INT ks;
  for (ks = 0; ks < N; ks++) {
    INT pos = (ks < Nneg) ? n - Nneg + ks : ks - Nneg;
    f_hat[ks] = g[pos] * pln->phi_hut_inv[ks];
  }
}

static void print(const plan *ego_, printer *pr)
{
  const deconv_plan *pln = (const deconv_plan *)ego_;
  pr->print(pr, "(deconv_solver_1d pcost=%D)", (INT)pln->super.pcost);
}

static void destroy(plan *ego_)
{
  deconv_plan *pln = (deconv_plan *)ego_;
  Y(free)(pln->phi_hut_inv);
}
static const plan_adt deconv_plan_adt = {apply, awake, print, destroy,
                                         apply_adjoint};

/* d == 1 only */
static plan *mkplan_deconv_1d(const solver *ego, const problem *p, planner *pl)
{
  const problem_deconv *pd = (const problem_deconv *)p;
  deconv_plan *pln;
  (void)ego;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_DECONV)
    return 0;
  if (pd->sz->rnk != 1)
    return 0;
  if (pd->window < NFFT_WINDOW_KAISER_BESSEL
      || pd->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */
  if (Y(problem_deconv_n)(p, 0) < Y(problem_deconv_N)(p, 0))
    return 0; /* n < N aliases grid cells */

  pln = (deconv_plan *)Y(plan_create)(sizeof(deconv_plan), &deconv_plan_adt);
  pln->n = Y(problem_deconv_ntot)(p);
  pln->N = Y(problem_deconv_Ntot)(p);
  pln->m = pd->m;
  pln->window = pd->window;
  pln->Nneg =
       pln->N / 2
       - (Y(problem_deconv_variant)(p, 0) == NFFT_NDFT_TYPE_II ? (INT)1 : (INT)0);
  pln->Npos = pln->N - pln->Nneg;
  pln->phi_hut_inv = (R *)Y(malloc)((size_t)pln->N * sizeof(R));
  pln->level = PLNR_SLEEPY;
  pln->super.pcost = Y(deconv_d_pcost)(p);
  return &pln->super;
}

static const solver_adt deconv_1d_adt = {NFFT_PROBLEM_DECONV, 0,
                                         mkplan_deconv_1d};
void Y(deconv_solver_1d_register)(planner *pl)
{
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &deconv_1d_adt));
}
