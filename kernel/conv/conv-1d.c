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

/* 1D CONV solver: Step C of the fast NFFT decomposition -- the node convolution 
 * (matrix B). Sums the oversampled grid g against the window psi at each 
 * nonequispaced node x_j (forward), or scatter-adds f onto g with the same psi 
 * weights (adjoint). psi depends on x/window/n/N/m), precomputed once at awake
 * (sparse PRE_PSI strategy in legacy code). */

/* Performance toggle.  Define CONV_OPTIMIZED for the fast path:
 *   #1  the inner loop walks two contiguous, vectorizable runs instead of the
 *       reference's per-element ((idx % n) + n) % n wrap, and
 *   #2  the wrapped window start u[j] is precomputed once at awake, so apply
 *       has no FLOOR/LRINT or modulo at all.
 * Comment the line out to fall back to the reference implementation; both are
 * exact math equivalents (identical index sequence). */
#define CONV_OPTIMIZED

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
  R *psi;     /* length M*(2m+2): psi[j*(2m+2)+lj], built at awake */
#ifdef CONV_OPTIMIZED
  INT *u;     /* length M: wrapped window start per node, built at awake */
#endif
  int precomputed;
} conv_plan;

/* precompute the psi table */
static void conv_awake(plan *ego_, int wakefulness) {
  conv_plan *pln = (conv_plan *)ego_;
  if (wakefulness >= PLNR_AWAKE_ZERO) {
    if (!pln->precomputed) {
      Y(window_phi_precompute)
      (pln->window, pln->n, pln->N, pln->m,
       pln->x + 0, 1, pln->M,
       pln->psi + 0 * (2 * pln->m + 2), 2 * pln->m + 2);
#ifdef CONV_OPTIMIZED
      { /* precompute the wrapped window start once */
        INT n = pln->n, M = pln->M, j;
        int m = pln->m;
        for (j = 0; j < M; j++) {
          INT c = LRINT(FLOOR((R)n * pln->x[j]));
          pln->u[j] = (((c - m) % n) + n) % n;
        }
      }
#endif
      pln->precomputed = 1;
    }
  } else
    pln->precomputed = 0; /* -> SLEEPY: psi values now stale */
}

#ifdef CONV_OPTIMIZED

static void conv_apply(const plan *ego_, const problem *p) {
  const conv_plan *pln = (const conv_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  INT n = pln->n, M = pln->M;
  int m = pln->m;
  const C *g = pc->g;
  C *f = pc->f;
  const INT len = 2 * (INT)m + 2; /* window support: c-m .. c+m+1 */
  INT j;
  for (j = 0; j < M; j++) {
    INT u = pln->u[j]; /* wrapped start, precomputed at awake */
    const R *psij = pln->psi + j * len;
    C acc = K(0.0);
    INT lj;
    if (u + len <= n) { /* no wrap: single contiguous run */
      for (lj = 0; lj < len; lj++)
        acc += g[u + lj] * psij[lj];
    } else { /* wraps once: tail of g, then head of g */
      INT k = 0;
      for (; u + k < n; k++)
        acc += g[u + k] * psij[k];
      for (; k < len; k++)
        acc += g[k - (n - u)] * psij[k];
    }
    f[j] = acc;
  }
}

static void conv_apply_adjoint(const plan *ego_, const problem *p) {
  const conv_plan *pln = (const conv_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  INT n = pln->n, M = pln->M;
  int m = pln->m;
  const C *f = pc->f;
  C *g = pc->g;
  const INT len = 2 * (INT)m + 2; /* window support: c-m .. c+m+1 */
  INT j;
  /* The scatter accumulates (+=) into an overlapping, node-dependent set that
   * does not cover the grid, so the whole grid must start zeroed. */
  memset(g, 0, (size_t)n * sizeof(C));
  for (j = 0; j < M; j++) {
    INT u = pln->u[j]; /* wrapped start, precomputed at awake */
    const R *psij = pln->psi + j * len;
    C fj = f[j];
    INT lj;
    if (u + len <= n) { /* no wrap: single contiguous run */
      for (lj = 0; lj < len; lj++)
        g[u + lj] += fj * psij[lj];
    } else { /* wraps once: tail of g, then head of g */
      INT k = 0;
      for (; u + k < n; k++)
        g[u + k] += fj * psij[k];
      for (; k < len; k++)
        g[k - (n - u)] += fj * psij[k];
    }
  }
}

#else /* reference implementation: per-element modulo wrap */

static void conv_apply(const plan *ego_, const problem *p) {
  const conv_plan *pln = (const conv_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  INT n = pln->n, M = pln->M;
  int m = pln->m;
  const C *g = pc->g;
  C *f = pc->f;
  INT j;
  for (j = 0; j < M; j++) {
    INT c = LRINT(FLOOR((R)n * pln->x[j]));
    C acc = K(0.0);
    int lj;
    for (lj = 0; lj <= 2 * m + 1; lj++) {
      INT idx = c - m + (INT)lj;
      INT wrap = ((idx % n) + n) % n;
      acc += g[wrap] * pln->psi[j * (2 * m + 2) + lj];
    }
    f[j] = acc;
  }
}

static void conv_apply_adjoint(const plan *ego_, const problem *p) {
  const conv_plan *pln = (const conv_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  INT n = pln->n, M = pln->M;
  int m = pln->m;
  const C *f = pc->f;
  C *g = pc->g;
  INT j;
  memset(g, 0, (size_t)n * sizeof(C)); /* zero the oversampled grid */
  for (j = 0; j < M; j++) {
    INT c = LRINT(FLOOR((R)n * pln->x[j]));
    int lj;
    for (lj = 0; lj <= 2 * m + 1; lj++) {
      INT idx = c - m + (INT)lj;
      INT wrap = ((idx % n) + n) % n;
      g[wrap] += f[j] * pln->psi[j * (2 * m + 2) + lj];
    }
  }
}

#endif /* CONV_OPTIMIZED */

static void conv_print(const plan *ego_, printer *pr) {
  const conv_plan *pln = (const conv_plan *)ego_;
  pr->print(pr, "(conv_solver_1d pcost=%D)", (INT)pln->super.pcost);
}
static void conv_destroy(plan *ego_) {
  conv_plan *pln = (conv_plan *)ego_;
  Y(free)
  (pln->psi);
#ifdef CONV_OPTIMIZED
  Y(free)
  (pln->u);
#endif
  /* x/g/f are borrowed caller arrays. */
}
static const plan_adt conv_plan_adt = {conv_apply, conv_awake, conv_print,
                                       conv_destroy, conv_apply_adjoint};

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
#ifdef CONV_OPTIMIZED
  pln->u = (INT *)Y(malloc)((size_t)pln->M * sizeof(INT));
#endif
  pln->precomputed = 0;
  pln->super.pcost = 2.0 * (double)pln->M * (double)(2 * pln->m + 2);
  return &pln->super;
}

static const solver_adt conv_1d_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_1d};
void Y(conv_solver_1d_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &conv_1d_adt));
}
