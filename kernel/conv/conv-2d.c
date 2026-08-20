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

/* 2D CONV solver: Step C of the fast NFFT decomposition -- the node convolution 
 * (matrix B). Sums the oversampled grid g against the window psi at each 
 * nonequispaced node x_j (forward), or scatter-adds f onto g with the same psi 
 * weights (adjoint). psi depends on x/window/n/N/m), precomputed once at awake
 * (sparse PRE_PSI strategy in legacy code). */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "conv.h"

double Y(conv_b_pcost)(const problem *p) {
  const problem_conv *pc = (const problem_conv *)p;
  double s = (double)(2 * pc->m + 2);
  return 2.0 * (double)pc->M * s * s;
}

typedef struct
{
  plan super;
  INT n0, n1, N0, N1, M; /* geometry captured at mkplan */
  int m, window;
  const R *x; /* borrowed alias of the problem's nodes */
  R *psi;     /* length M*2*(2m+2): psi[(j*2+t)*(2m+2)+lj], built at awake */
  int precomputed;
} conv_2d_plan;

/* uo2: neighbor window start/end on axis of size n, wrapped mod n. */
static void uo2(INT *u, INT *o, const R x, const INT n, const INT m) {
  INT c = LRINT(FLOOR(x * (R)n));
  *u = (c - m + n) % n;
  *o = (c + 1 + m + n) % n;
}

/* precompute the psi table from ego->x. */
static void conv_2d_awake(plan *ego_, int wakefulness) {
  conv_2d_plan *pln = (conv_2d_plan *)ego_;
  if (wakefulness >= PLNR_AWAKE_ZERO) {
    if (!pln->precomputed) {
      int t;
      INT nn[2];
      nn[0] = pln->n0;
      nn[1] = pln->n1;
      INT NN[2];
      NN[0] = pln->N0;
      NN[1] = pln->N1;
      for (t = 0; t < 2; t++)
        Y(window_phi_precompute)
      (pln->window, nn[t], NN[t], pln->m,
       pln->x + t, 2, pln->M,
       pln->psi + t * (2 * pln->m + 2), 2 * (2 * pln->m + 2));
      pln->precomputed = 1;
    }
  } else
    pln->precomputed = 0; /* -> SLEEPY: psi values now stale */
}

/* Forward B (g -> f[j]) */
static void conv_trafo_2d_compute(C *fj, const C *g, const R *psij_const0,
                                  const R *psij_const1, const R *xj0, const R *xj1, const INT n0,
                                  const INT n1, const int m) {
  INT u0, o0, l0, u1, o1, l1;
  const C *gj;
  const R *psij0, *psij1;

  psij0 = psij_const0;
  psij1 = psij_const1;

  uo2(&u0, &o0, *xj0, n0, m);
  uo2(&u1, &o1, *xj1, n1, m);

  *fj = K(0.0);

  if (u0 < o0)
    if (u1 < o1)
      for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
        psij1 = psij_const1;
        gj = g + (u0 + l0) * n1 + u1;
        for (l1 = 0; l1 <= 2 * m + 1; l1++)
          (*fj) += (*psij0) * (*psij1++) * (*gj++);
      }
    else
      for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
        psij1 = psij_const1;
        gj = g + (u0 + l0) * n1 + u1;
        for (l1 = 0; l1 < 2 * m + 1 - o1; l1++)
          (*fj) += (*psij0) * (*psij1++) * (*gj++);
        gj = g + (u0 + l0) * n1;
        for (l1 = 0; l1 <= o1; l1++)
          (*fj) += (*psij0) * (*psij1++) * (*gj++);
      }
  else if (u1 < o1) {
    for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
      psij1 = psij_const1;
      gj = g + (u0 + l0) * n1 + u1;
      for (l1 = 0; l1 <= 2 * m + 1; l1++)
        (*fj) += (*psij0) * (*psij1++) * (*gj++);
    }
    for (l0 = 0; l0 <= o0; l0++, psij0++) {
      psij1 = psij_const1;
      gj = g + l0 * n1 + u1;
      for (l1 = 0; l1 <= 2 * m + 1; l1++)
        (*fj) += (*psij0) * (*psij1++) * (*gj++);
    }
  } else {
    for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
      psij1 = psij_const1;
      gj = g + (u0 + l0) * n1 + u1;
      for (l1 = 0; l1 < 2 * m + 1 - o1; l1++)
        (*fj) += (*psij0) * (*psij1++) * (*gj++);
      gj = g + (u0 + l0) * n1;
      for (l1 = 0; l1 <= o1; l1++)
        (*fj) += (*psij0) * (*psij1++) * (*gj++);
    }
    for (l0 = 0; l0 <= o0; l0++, psij0++) {
      psij1 = psij_const1;
      gj = g + l0 * n1 + u1;
      for (l1 = 0; l1 < 2 * m + 1 - o1; l1++)
        (*fj) += (*psij0) * (*psij1++) * (*gj++);
      gj = g + l0 * n1;
      for (l1 = 0; l1 <= o1; l1++)
        (*fj) += (*psij0) * (*psij1++) * (*gj++);
    }
  }
}

/* Adjoint B^H (f[j] -> g, scatter-add) */
static void conv_adjoint_2d_compute_serial(const C *fj, C *g,
                                           const R *psij_const0, const R *psij_const1, const R *xj0, const R *xj1,
                                           const INT n0, const INT n1, const int m) {
  INT u0, o0, l0, u1, o1, l1;
  C *gj;
  const R *psij0, *psij1;

  psij0 = psij_const0;
  psij1 = psij_const1;

  uo2(&u0, &o0, *xj0, n0, m);
  uo2(&u1, &o1, *xj1, n1, m);

  if (u0 < o0)
    if (u1 < o1)
      for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
        psij1 = psij_const1;
        gj = g + (u0 + l0) * n1 + u1;
        for (l1 = 0; l1 <= 2 * m + 1; l1++)
          (*gj++) += (*psij0) * (*psij1++) * (*fj);
      }
    else
      for (l0 = 0; l0 <= 2 * m + 1; l0++, psij0++) {
        psij1 = psij_const1;
        gj = g + (u0 + l0) * n1 + u1;
        for (l1 = 0; l1 < 2 * m + 1 - o1; l1++)
          (*gj++) += (*psij0) * (*psij1++) * (*fj);
        gj = g + (u0 + l0) * n1;
        for (l1 = 0; l1 <= o1; l1++)
          (*gj++) += (*psij0) * (*psij1++) * (*fj);
      }
  else if (u1 < o1) {
    for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
      psij1 = psij_const1;
      gj = g + (u0 + l0) * n1 + u1;
      for (l1 = 0; l1 <= 2 * m + 1; l1++)
        (*gj++) += (*psij0) * (*psij1++) * (*fj);
    }
    for (l0 = 0; l0 <= o0; l0++, psij0++) {
      psij1 = psij_const1;
      gj = g + l0 * n1 + u1;
      for (l1 = 0; l1 <= 2 * m + 1; l1++)
        (*gj++) += (*psij0) * (*psij1++) * (*fj);
    }
  } else {
    for (l0 = 0; l0 < 2 * m + 1 - o0; l0++, psij0++) {
      psij1 = psij_const1;
      gj = g + (u0 + l0) * n1 + u1;
      for (l1 = 0; l1 < 2 * m + 1 - o1; l1++)
        (*gj++) += (*psij0) * (*psij1++) * (*fj);
      gj = g + (u0 + l0) * n1;
      for (l1 = 0; l1 <= o1; l1++)
        (*gj++) += (*psij0) * (*psij1++) * (*fj);
    }
    for (l0 = 0; l0 <= o0; l0++, psij0++) {
      psij1 = psij_const1;
      gj = g + l0 * n1 + u1;
      for (l1 = 0; l1 < 2 * m + 1 - o1; l1++)
        (*gj++) += (*psij0) * (*psij1++) * (*fj);
      gj = g + l0 * n1;
      for (l1 = 0; l1 <= o1; l1++)
        (*gj++) += (*psij0) * (*psij1++) * (*fj);
    }
  }
}

static void conv_2d_apply(const plan *ego_, const problem *p) {
  const conv_2d_plan *pln = (const conv_2d_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  INT n0 = pln->n0, n1 = pln->n1, M = pln->M;
  int m = pln->m;
  const C *g = pc->g;
  C *f = pc->f;
  INT j;
  for (j = 0; j < M; j++) {
    const R *psij0 = &pln->psi[(j * 2 + 0) * (2 * m + 2)];
    const R *psij1 = &pln->psi[(j * 2 + 1) * (2 * m + 2)];
    const R *xj0 = &pln->x[j * 2 + 0];
    const R *xj1 = &pln->x[j * 2 + 1];
    conv_trafo_2d_compute(&f[j], g, psij0, psij1, xj0, xj1, n0, n1, m);
  }
}

static void conv_2d_apply_adjoint(const plan *ego_, const problem *p) {
  const conv_2d_plan *pln = (const conv_2d_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  INT n0 = pln->n0, n1 = pln->n1, M = pln->M, ntot = n0 * n1;
  int m = pln->m;
  const C *f = pc->f;
  C *g = pc->g;
  INT j;
  /* The scatter accumulates (+=) into an overlapping, node-dependent set that
   * does not cover the grid, so the whole grid must start zeroed. */
  memset(g, 0, (size_t)ntot * sizeof(C)); /* zero the oversampled grid */
  for (j = 0; j < M; j++) {
    const R *psij0 = &pln->psi[(j * 2 + 0) * (2 * m + 2)];
    const R *psij1 = &pln->psi[(j * 2 + 1) * (2 * m + 2)];
    const R *xj0 = &pln->x[j * 2 + 0];
    const R *xj1 = &pln->x[j * 2 + 1];
    conv_adjoint_2d_compute_serial(&f[j], g, psij0, psij1, xj0, xj1, n0, n1, m);
  }
}

static void conv_2d_print(const plan *ego_, printer *pr) {
  const conv_2d_plan *pln = (const conv_2d_plan *)ego_;
  pr->print(pr, "(conv_solver_2d pcost=%D)", (INT)pln->super.pcost);
}
static void conv_2d_destroy(plan *ego_) {
  conv_2d_plan *pln = (conv_2d_plan *)ego_;
  Y(free)
  (pln->psi);
  /* x/g/f are borrowed caller arrays. */
}
static const plan_adt conv_2d_plan_adt = {conv_2d_apply, conv_2d_awake,
                                          conv_2d_print, conv_2d_destroy,
                                          conv_2d_apply_adjoint};

/* d == 2 only */
static plan *mkplan_conv_2d(const solver *ego, const problem *p, planner *pl) {
  const problem_conv *pc = (const problem_conv *)p;
  conv_2d_plan *pln;
  (void)ego;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_CONV)
    return 0;
  if (pc->sz->rnk != 2)
    return 0;
  if (pc->window < NFFT_WINDOW_KAISER_BESSEL ||
      pc->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */


  pln = (conv_2d_plan *)Y(plan_create)(sizeof(conv_2d_plan), &conv_2d_plan_adt);
  pln->n0 = Y(problem_conv_n)(p, 0);
  pln->n1 = Y(problem_conv_n)(p, 1);
  pln->N0 = Y(problem_conv_N)(p, 0);
  pln->N1 = Y(problem_conv_N)(p, 1);
  pln->M = pc->M;
  pln->m = pc->m;
  pln->window = pc->window;
  pln->x = pc->x; /* borrowed */
  pln->psi = (R *)Y(malloc)((size_t)pln->M * 2 * (size_t)(2 * pln->m + 2) * sizeof(R));
  pln->precomputed = 0;
  pln->super.pcost = Y(conv_b_pcost)(p);
  return &pln->super;
}

static const solver_adt conv_2d_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_2d};
void Y(conv_solver_2d_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &conv_2d_adt));
}
