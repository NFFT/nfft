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

/* nD (n>=4) CONV solver: Step C of the fast NFFT decomposition -- the node convolution 
 * (matrix B). Sums the oversampled grid g against the window psi at each 
 * nonequispaced node x_j (forward), or scatter-adds f onto g with the same psi 
 * weights (adjoint). psi depends on x/window/n/N/m), precomputed once at awake
 * (sparse PRE_PSI strategy in legacy code). */

/* not ported: the d==4 and d==5 hand-unrolled branches of MACRO_B_COMPUTE_ONE_NODE */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "conv.h"

typedef struct
{
  plan super;
  int d;      /* rank, >= 4 */
  INT *n, *N; /* owned, length d: geometry captured at mkplan */
  INT M;      /* node count */
  INT ntot;   /* owned product of n[], captured at mkplan */
  int m, window;
  const R *x; /* borrowed alias of the problem's nodes (length d*M), never freed */
  R *psi;     /* length M*d*(2m+2): psi[(j*d+t)*(2m+2)+lj], built at awake */
  int precomputed;
} conv_nd_plan;

/* uo: unwrapped neighbor window start/end for axis t at node j. */
static void uo(const R *x, int d, INT j, const INT *n, int m, int t, INT *up,
               INT *op) {
  const R xj = x[j * d + t];
  INT c = LRINT(FLOOR(xj * (R)n[t]));
  *up = c - m;
  *op = c + 1 + m;
}

/* precompute the psi table from ego->x. */
static void conv_nd_awake(plan *ego_, int wakefulness) {
  conv_nd_plan *pln = (conv_nd_plan *)ego_;
  if (wakefulness >= PLNR_AWAKE_ZERO) {
    if (!pln->precomputed) {
      int t;
      for (t = 0; t < pln->d; t++)
        Y(window_phi_precompute)
      (pln->window, pln->n[t], pln->N[t], pln->m,
       pln->x + t, pln->d, pln->M,
       pln->psi + t * (2 * pln->m + 2), pln->d * (2 * pln->m + 2));
      pln->precomputed = 1;
    }
  } else
    pln->precomputed = 0; /* -> SLEEPY: psi values now stale */
}

/* Forward B (g -> f[j]) / adjoint B^H (f[j] -> g, scatter-add). */
static void conv_nd_run(const conv_nd_plan *pln, const problem_conv *pc,
                        int forward) {
  const int d = pln->d;
  const INT *n = pln->n;
  const INT M = pln->M;
  const int m = pln->m;
  const R *x = pln->x;
  const R *psi = pln->psi;
  C *f, *g;
  INT lprod;
  INT u[d], o[d];             /* unwrapped multi band w.r.t. x_j */
  INT lj[d];                  /* multi index 0 <= lj <= o-u */
  INT ll_plain[d + 1];        /* postfix plain index in g */
  R phi_prod[d + 1];          /* postfix product of psi */
  INT l_all[d * (2 * m + 2)]; /* wrapped grid indices per axis/tap */
  INT j, t, t2, l_L;

  if (forward) {
    f = pc->f;
    g = (C *)pc->g;
    memset(f, 0, (size_t)M * sizeof(C)); /* MACRO_B_init_result_A */
  } else {
    f = (C *)pc->f;
    g = pc->g;
    memset(g, 0, (size_t)pln->ntot * sizeof(C)); /* MACRO_B_init_result_T */
  }

  for (t = 0, lprod = 1; t < d; t++)
    lprod *= (2 * m + 2);

  for (j = 0; j < M; j++) {
    /* MACRO_init_uo_l_lj_t */
    for (t = d - 1; t >= 0; t--) {
      INT lj_t;
      uo(x, d, j, n, m, t, &u[t], &o[t]);
      for (lj_t = 0; lj_t < 2 * m + 2; lj_t++)
        l_all[t * (2 * m + 2) + lj_t] = (u[t] + lj_t + n[t]) % n[t];
      lj[t] = 0;
    }
    t++;

    phi_prod[0] = K(1.0);
    ll_plain[0] = 0;

    for (l_L = 0; l_L < lprod; l_L++) {
      /* MACRO_update_phi_prod_ll_plain(with_PRE_PSI) */
      for (t2 = t; t2 < d; t2++) {
        phi_prod[t2 + 1] = phi_prod[t2] * psi[(j * d + t2) * (2 * m + 2) + lj[t2]];
        ll_plain[t2 + 1] = ll_plain[t2] * n[t2] + l_all[t2 * (2 * m + 2) + lj[t2]];
      }

      /* MACRO_B_compute_A / MACRO_B_compute_T */
      if (forward)
        f[j] += phi_prod[d] * g[ll_plain[d]];
      else
        g[ll_plain[d]] += phi_prod[d] * f[j];

      /* MACRO_count_uo_l_lj_t */
      for (t = d - 1; (t > 0) && (lj[t] == o[t] - u[t]); t--)
        lj[t] = 0;
      lj[t]++;
    }
  }
}

static void conv_nd_apply(const plan *ego_, const problem *p) {
  conv_nd_run((const conv_nd_plan *)ego_, (const problem_conv *)p, 1);
}

static void conv_nd_apply_adjoint(const plan *ego_, const problem *p) {
  conv_nd_run((const conv_nd_plan *)ego_, (const problem_conv *)p, 0);
}

static void conv_nd_print(const plan *ego_, printer *pr) {
  const conv_nd_plan *pln = (const conv_nd_plan *)ego_;
  pr->print(pr, "(conv_solver_nd pcost=%D)", (INT)pln->super.pcost);
}
static void conv_nd_destroy(plan *ego_) {
  conv_nd_plan *pln = (conv_nd_plan *)ego_;
  Y(free)
  (pln->psi);
  Y(free)
  (pln->N);
  Y(free)
  (pln->n);
  /* x/g/f are borrowed caller arrays. */
}
static const plan_adt conv_nd_plan_adt = {conv_nd_apply, conv_nd_awake,
                                          conv_nd_print, conv_nd_destroy,
                                          conv_nd_apply_adjoint};

/* d >= 4 only. */
static plan *mkplan_conv_nd(const solver *ego, const problem *p, planner *pl) {
  const problem_conv *pc = (const problem_conv *)p;
  conv_nd_plan *pln;
  int t, d;
  (void)ego;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_CONV)
    return 0;
  if (pc->sz->rnk < 4)
    return 0;
  if (pc->window < NFFT_WINDOW_KAISER_BESSEL ||
      pc->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */

  d = pc->sz->rnk;
  pln = (conv_nd_plan *)Y(plan_create)(sizeof(conv_nd_plan), &conv_nd_plan_adt);
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
  /* conv_nd_run holds a per-run stack VLA l_all[d*(2m+2)] plus
   * sibling VLAs.  The composed solver never presents a pathological rank/m,
   * but bound the stack budget so a future direct problem_conv cannot overflow
   * the stack. */
  A((size_t)d * (size_t)(2 * pln->m + 2) * sizeof(INT) <= (size_t)(64 * 1024));
  pln->window = pc->window;
  pln->x = pc->x; /* borrowed alias; never freed here */
  pln->psi = (R *)Y(malloc)(
      (size_t)pln->M * (size_t)d * (size_t)(2 * pln->m + 2) * sizeof(R));
  pln->precomputed = 0;
  pln->super.pcost = Y(conv_b_pcost)(p);
  return &pln->super;
}

static const solver_adt conv_nd_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_nd};
void Y(conv_solver_nd_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &conv_nd_adt));
}
