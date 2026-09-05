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

/* nD DECONV solver: step A of the fast NFFT decomposition. Forward divides f_hat
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
  int d; /* rank, >= 4 */
  INT *N, *n; /* owned, length d: geometry captured at mkplan */
  INT *Nneg; /* owned, length d: slot split, k(ks) = ks - Nneg[t] */
  INT Ntot, ntot; /* owned products, captured at mkplan */
  int m, window;
  R **phi_hut_inv; /* owned array of d owned tables; phi_hut_inv[t] has
                    * length N[t]: 1/phi_hut(n[t],N[t],m,ks-Nneg[t]) */
  int level; /* content of the tables: SLEEPY (stale), AWAKE_ZERO or AWAKE */
} deconv_nd_plan;

/* AWAKE_ZERO must cost no window evaluation, so the tables get placeholder
 * zeros. AWAKE_ZERO reached by downgrade only drops the level, so a later
 * upgrade refills. */
static void awake(plan *ego_, int wakefulness)
{
  deconv_nd_plan *pln = (deconv_nd_plan *)ego_;
  int t;
  if (wakefulness == PLNR_AWAKE) {
    if (pln->level != PLNR_AWAKE)
      for (t = 0; t < pln->d; t++) {
        INT Nt = pln->N[t], ks;
        Y(window_phi_hut_apply)(pln->window, pln->n[t], Nt, pln->m, -pln->Nneg[t],
                             pln->phi_hut_inv[t], Nt);
        for (ks = 0; ks < Nt; ks++)
          pln->phi_hut_inv[t][ks] = K(1.0) / pln->phi_hut_inv[t][ks];
      }
  } else if (wakefulness == PLNR_AWAKE_ZERO && pln->level == PLNR_SLEEPY)
    for (t = 0; t < pln->d; t++)
      memset(pln->phi_hut_inv[t], 0, (size_t)pln->N[t] * sizeof(R));
  pln->level = wakefulness;
}

static void run(const deconv_nd_plan *pln, const problem_deconv *pd,
                int forward)
{
  const int d = pln->d;
  const INT *N = pln->N;
  const INT *n = pln->n;
  const INT *Nneg = pln->Nneg;
  R *const *phi_hut_inv = (R *const *)pln->phi_hut_inv;
  C *f_hat, *g_hat;
  R c_phi_inv_k[d + 1]; /* postfix product of phi_hut_inv */
  INT t, t2; /* index dimensions */
  INT k_L; /* plain index */
  INT kp[d]; /* multi index (simple) */
  INT k[d]; /* multi index in g_hat */
  INT ks[d]; /* multi index in f_hat, phi_hut_inv */
  INT k_plain[d + 1]; /* postfix plain index */
  INT ks_plain[d + 1]; /* postfix plain index */

  if (forward) {
    f_hat = (C *)pd->f_hat;
    g_hat = pd->g;
    memset(g_hat, 0, (size_t)pln->ntot * sizeof(C));
  } else {
    f_hat = pd->f_hat;
    g_hat = (C *)pd->g;
    memset(f_hat, 0, (size_t)pln->Ntot * sizeof(C));
  }

  c_phi_inv_k[0] = K(1.0);
  k_plain[0] = 0;
  ks_plain[0] = 0;

  for (t = d - 1; 0 <= t; t--) {
    kp[t] = k[t] = 0;
    ks[t] = Nneg[t];
  }
  t++;

  for (k_L = 0; k_L < pln->Ntot; k_L++) {
    for (t2 = t; t2 < d; t2++) {
      c_phi_inv_k[t2 + 1] = c_phi_inv_k[t2] * phi_hut_inv[t2][ks[t2]];
      ks_plain[t2 + 1] = ks_plain[t2] * N[t2] + ks[t2];
      k_plain[t2 + 1] = k_plain[t2] * n[t2] + k[t2];
    }

    if (forward)
      g_hat[k_plain[d]] = f_hat[ks_plain[d]] * c_phi_inv_k[d];
    else
      f_hat[ks_plain[d]] = g_hat[k_plain[d]] * c_phi_inv_k[d];

    for (t = d - 1; (t > 0) && (kp[t] == N[t] - 1); t--) {
      kp[t] = k[t] = 0;
      ks[t] = Nneg[t];
    }
    kp[t]++;
    k[t]++;
    ks[t]++;
    /* The non-negative run has Npos = N - Nneg slots; after it, wrap to the
     * grid tail and restart at f_hat slot 0. */
    if (kp[t] == N[t] - Nneg[t]) {
      k[t] = n[t] - Nneg[t];
      ks[t] = 0;
    }
  }
}

static void apply(const plan *ego_, const problem *p)
{
  run((const deconv_nd_plan *)ego_, (const problem_deconv *)p, 1);
}

static void apply_adjoint(const plan *ego_, const problem *p)
{
  run((const deconv_nd_plan *)ego_, (const problem_deconv *)p, 0);
}

static void print(const plan *ego_, printer *pr)
{
  const deconv_nd_plan *pln = (const deconv_nd_plan *)ego_;
  pr->print(pr, "(deconv_solver_nd pcost=%D)", (INT)pln->super.pcost);
}
static void destroy(plan *ego_)
{
  deconv_nd_plan *pln = (deconv_nd_plan *)ego_;
  int t;
  for (t = 0; t < pln->d; t++)
    Y(free)(pln->phi_hut_inv[t]);
  Y(free)(pln->phi_hut_inv);
  Y(free)(pln->Nneg);
  Y(free)(pln->n);
  Y(free)(pln->N);
}
static const plan_adt deconv_nd_plan_adt = {apply, awake, print, destroy,
                                            apply_adjoint};

/* d >= 4 only */
static plan *mkplan_deconv_nd(const solver *ego, const problem *p, planner *pl)
{
  const problem_deconv *pd = (const problem_deconv *)p;
  deconv_nd_plan *pln;
  int t, d;
  (void)ego;
  (void)pl;
  if (p->adt->kind != NFFT_PROBLEM_DECONV)
    return 0;
  if (pd->sz->rnk < 4)
    return 0;
  if (pd->window < NFFT_WINDOW_KAISER_BESSEL
      || pd->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */

  d = pd->sz->rnk;
  for (t = 0; t < d; t++)
    if (Y(problem_deconv_n)(p, t) < Y(problem_deconv_N)(p, t))
      return 0; /* n < N aliases grid cells */

  pln = (deconv_nd_plan *)Y(plan_create)(sizeof(deconv_nd_plan),
                                      &deconv_nd_plan_adt);
  pln->d = d;
  pln->N = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->n = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->Nneg = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->phi_hut_inv = (R **)Y(malloc)((size_t)d * sizeof(R *));
  pln->Ntot = Y(problem_deconv_Ntot)(p);
  pln->ntot = Y(problem_deconv_ntot)(p);
  for (t = 0; t < d; t++) {
    pln->N[t] = Y(problem_deconv_N)(p, t);
    pln->n[t] = Y(problem_deconv_n)(p, t);
    pln->Nneg[t] = pln->N[t] / 2
                   - (Y(problem_deconv_variant)(p, t) == NFFT_NDFT_TYPE_II ?
                           (INT)1 :
                           (INT)0);
    pln->phi_hut_inv[t] = (R *)Y(malloc)((size_t)pln->N[t] * sizeof(R));
  }
  pln->m = pd->m;
  pln->window = pd->window;
  pln->level = PLNR_SLEEPY;
  pln->super.pcost = Y(deconv_d_pcost)(p);
  return &pln->super;
}

static const solver_adt deconv_nd_adt = {NFFT_PROBLEM_DECONV, 0,
                                         mkplan_deconv_nd};
void Y(deconv_solver_nd_register)(planner *pl)
{
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &deconv_nd_adt));
}
