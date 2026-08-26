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

/* 2D DECONV solver: Step A of the fast NFFT decomposition -- deconvolve f_hat by
 * the window's phi_hut factors and zero-pad onto the oversampled grid g (forward),
 * or the adjoint gather (g -> f_hat, multiplying by the same 1/phi_hut. phi_hut
 * depends only on (n, N, m, window), so it is precomputed once at awake,
 * node-independent. Slot ks carries frequency k = ks - Nneg, with Nneg = N/2
 * for type-I and N/2 - 1 for type-II; odd N normalizes to type-I in
 * mkproblem_deconv, so type-II implies even N. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "deconv.h"

double Y(deconv_d_pcost)(const problem *p) {
  return 2.0 * (double)Y(problem_deconv_Ntot)(p);
}

typedef struct
{
  plan super;
  INT N0, N1, n0, n1;         /* geometry captured at mkplan */
  INT Nneg0, Npos0, Nneg1, Npos1; /* per-axis slot split */
  int m, window;
  R *phi_hut_inv0; /* length N0: 1/phi_hut(n0,N0,m, ks0 - Nneg0), at awake */
  R *phi_hut_inv1; /* length N1: 1/phi_hut(n1,N1,m, ks1 - Nneg1), at awake */
  int precomputed;
} deconv_2d_plan;

/* precompute 1/phi_hut */
static void deconv_2d_awake(plan *ego_, int wakefulness) {
  deconv_2d_plan *pln = (deconv_2d_plan *)ego_;
  if (wakefulness >= PLNR_AWAKE_ZERO) {
    if (!pln->precomputed) {
      INT ks;
      Y(window_phi_hut_apply)
      (pln->window, pln->n0, pln->N0, pln->m, -pln->Nneg0, pln->phi_hut_inv0,
       pln->N0);
      for (ks = 0; ks < pln->N0; ks++)
        pln->phi_hut_inv0[ks] = K(1.0) / pln->phi_hut_inv0[ks];
      Y(window_phi_hut_apply)
      (pln->window, pln->n1, pln->N1, pln->m, -pln->Nneg1, pln->phi_hut_inv1,
       pln->N1);
      for (ks = 0; ks < pln->N1; ks++)
        pln->phi_hut_inv1[ks] = K(1.0) / pln->phi_hut_inv1[ks];
      pln->precomputed = 1;
    }
  } else
    pln->precomputed = 0;
}

/* Each axis splits into two contiguous runs of slots: the negative half maps to
 * the grid tail, the non-negative half to the grid head. Even type-I makes the
 * two runs equal (Nneg == Npos == N/2). */
static void deconv_2d_run(const deconv_2d_plan *pln, const problem_deconv *pd,
                          int forward) {
  const INT N1 = pln->N1, n0 = pln->n0, n1 = pln->n1;
  INT len0[2], sof0[2], gof0[2]; /* run length, f_hat slot offset, grid offset */
  INT len1[2], sof1[2], gof1[2];
  C *f_hat = pd->f_hat;
  C *g_hat = pd->g;
  int a, b;

  len0[0] = pln->Nneg0;
  sof0[0] = 0;
  gof0[0] = n0 - pln->Nneg0;
  len0[1] = pln->Npos0;
  sof0[1] = pln->Nneg0;
  gof0[1] = 0;
  len1[0] = pln->Nneg1;
  sof1[0] = 0;
  gof1[0] = n1 - pln->Nneg1;
  len1[1] = pln->Npos1;
  sof1[1] = pln->Nneg1;
  gof1[1] = 0;

  /* In nD the touched frequencies form a box per grid corner, so the zero-pad is
   * fragmented; clearing the whole grid is typically more efficient. The adjoint
   * writes every f_hat slot, so it needs no clear. */
  if (forward)
    memset(g_hat, 0, (size_t)(n0 * n1) * sizeof(C));

  for (a = 0; a < 2; a++)
    for (b = 0; b < 2; b++) {
      const R *p0 = pln->phi_hut_inv0 + sof0[a];
      const R *p1 = pln->phi_hut_inv1 + sof1[b];
      INT i, j;
      for (i = 0; i < len0[a]; i++) {
        R ck0 = p0[i];
        C *gh = g_hat + (gof0[a] + i) * n1 + gof1[b];
        C *fh = f_hat + (sof0[a] + i) * N1 + sof1[b];
        /* (C * R) * R is four real multiplies; C * (R * R) is three. */
        if (forward)
          for (j = 0; j < len1[b]; j++)
            gh[j] = fh[j] * (ck0 * p1[j]);
        else
          for (j = 0; j < len1[b]; j++)
            fh[j] = gh[j] * (ck0 * p1[j]);
      }
    }
}

/* apply the real diagonal scale-and-pad map (f_hat -> g). */
static void deconv_2d_apply(const plan *ego_, const problem *p) {
  deconv_2d_run((const deconv_2d_plan *)ego_, (const problem_deconv *)p, 1);
}

/* apply the adjoint real diagonal scale-and-pad map (g -> f_hat). The adjoint only swaps scatter->gather. */
static void deconv_2d_apply_adjoint(const plan *ego_, const problem *p) {
  deconv_2d_run((const deconv_2d_plan *)ego_, (const problem_deconv *)p, 0);
}

static void deconv_2d_print(const plan *ego_, printer *pr) {
  const deconv_2d_plan *pln = (const deconv_2d_plan *)ego_;
  pr->print(pr, "(deconv_solver_2d pcost=%D)", (INT)pln->super.pcost);
}

static void deconv_2d_destroy(plan *ego_) {
  deconv_2d_plan *pln = (deconv_2d_plan *)ego_;
  Y(free)
  (pln->phi_hut_inv1);
  Y(free)
  (pln->phi_hut_inv0);
}

static const plan_adt deconv_2d_plan_adt = {deconv_2d_apply, deconv_2d_awake,
                                            deconv_2d_print, deconv_2d_destroy,
                                            deconv_2d_apply_adjoint};

/* d == 2 only */
static plan *mkplan_deconv_2d(const solver *ego, const problem *p, planner *pl) {
  const problem_deconv *pd = (const problem_deconv *)p;
  deconv_2d_plan *pln;
  (void)ego;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_DECONV)
    return 0;
  if (pd->sz->rnk != 2)
    return 0;
  if (pd->window < NFFT_WINDOW_KAISER_BESSEL ||
      pd->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */
  if (Y(problem_deconv_n)(p, 0) < Y(problem_deconv_N)(p, 0) ||
      Y(problem_deconv_n)(p, 1) < Y(problem_deconv_N)(p, 1))
    return 0; /* n < N aliases grid cells */

  pln = (deconv_2d_plan *)Y(plan_create)(sizeof(deconv_2d_plan), &deconv_2d_plan_adt);
  pln->N0 = Y(problem_deconv_N)(p, 0);
  pln->N1 = Y(problem_deconv_N)(p, 1);
  pln->n0 = Y(problem_deconv_n)(p, 0);
  pln->n1 = Y(problem_deconv_n)(p, 1);
  pln->Nneg0 = pln->N0 / 2
               - (Y(problem_deconv_variant)(p, 0) == NFFT_NDFT_TYPE_II ?
                       (INT)1 :
                       (INT)0);
  pln->Npos0 = pln->N0 - pln->Nneg0;
  pln->Nneg1 = pln->N1 / 2
               - (Y(problem_deconv_variant)(p, 1) == NFFT_NDFT_TYPE_II ?
                       (INT)1 :
                       (INT)0);
  pln->Npos1 = pln->N1 - pln->Nneg1;
  pln->m = pd->m;
  pln->window = pd->window;
  pln->phi_hut_inv0 = (R *)Y(malloc)((size_t)pln->N0 * sizeof(R));
  pln->phi_hut_inv1 = (R *)Y(malloc)((size_t)pln->N1 * sizeof(R));
  pln->precomputed = 0;
  pln->super.pcost = Y(deconv_d_pcost)(p);
  return &pln->super;
}

static const solver_adt deconv_2d_adt = {NFFT_PROBLEM_DECONV, 0, mkplan_deconv_2d};
void Y(deconv_solver_2d_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &deconv_2d_adt));
}
