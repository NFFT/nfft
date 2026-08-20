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

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "ndft.h"

#include <string.h> /* memset — adjoint zeroes f_hat before accumulating */

/* 1D Direct NDFT. 
 *
 * Two changes compared to legacy code:
 *
 * (1) Correctly handle unit axes: the legacy code carries on `== N[t]/2 - 1`. For a
 * unit axis (N[t] == 1) the range [-N/2, N/2-1] is [0, -1], so k[t] starts at
 * 0 and `== -1` never fires. `>=` makes a unit axis carry transparently.
 *
 * (2) Support for odd bandwidth: threshold `(N[t]-1)/2` instead of `N[t]/2 - 1`: 
 * the two are identical for even N (integer division), so the even path is
 * unchanged. For odd N, `(N[t]-1)/2` is the correct upper bound. */
/* Block size for the phase recurrence. */

/* Block size for blocked recurrence variant. */
#define NDFT_RECURRENCE_BLOCK 32

/* Accurate phase for exp(+-i 2pi k x): reduce k*x to ~[-1/2,1/2) with a single-
 * rounding FMA so COS/SIN see a small argument and error does not grow with N. */
static inline R ndft_reduced_omega(const R k, const R x) {
  const R nrnd = RINT(k * x);
  return K2PI * FFMA(k, x, -nrnd);
}

/* NDFT plan: Variant (plain/blocked) is baked into the apply fns at mkplan time. */
typedef struct
{
  plan super;
  void (*apply_native_fwd)(const problem_nfft *pn);
  void (*apply_native_adj)(const problem_nfft *pn);
  const char *reg_nam;
} ndft_plan;

static void np_apply(const plan *ego, const problem *p) {
  const ndft_plan *pln = (const ndft_plan *)ego;
  pln->apply_native_fwd((const problem_nfft *)p);
}
static void np_apply_adjoint(const plan *ego, const problem *p) {
  const ndft_plan *pln = (const ndft_plan *)ego;
  pln->apply_native_adj((const problem_nfft *)p);
}

static void np_print(const plan *ego, printer *pr) {
  const ndft_plan *pln = (const ndft_plan *)ego;
  pr->print(pr, "(%s pcost=%D)", pln->reg_nam, (INT)pln->super.pcost);
}

/* No awake hook (nothing to precompute), no destroy hook. */
static const plan_adt ndft_plan_adt = {np_apply, 0, np_print, 0,
                                       np_apply_adjoint};

static void apply_plain(const problem_nfft *pn) {
  const INT Ntot = Y(problem_nfft_Ntot)((const problem *)pn);
  const INT M = pn->M;
  const R *x = pn->x;
  const C *f_hat = pn->f_hat;
  C *f = pn->f;
  const int type_ii = (pn->variant[0] == NFFT_NDFT_TYPE_II);
  INT j;
  for (j = 0; j < M; j++) {
    C v = K(0.0);
    INT k_L;
    for (k_L = 0; k_L < Ntot; k_L++) {
      /* ascending type-II: uniform +1 shift (odd N normalized to type-I -> 0) */
      R omega = K2PI * ((R)(k_L - Ntot / 2 + type_ii)) * x[j];
      v += f_hat[k_L] * (COS(omega) - II * SIN(omega));
    }
    f[j] = v;
  }
}

static void apply_plain_adjoint(const problem_nfft *pn) {
  const INT Ntot = Y(problem_nfft_Ntot)((const problem *)pn);
  const INT M = pn->M;
  const R *x = pn->x;
  const C *f = pn->f;
  C *f_hat = pn->f_hat;
  const int type_ii = (pn->variant[0] == NFFT_NDFT_TYPE_II);
  INT j;
  memset(f_hat, 0, (size_t)Ntot * sizeof(C));
  for (j = 0; j < M; j++) {
    INT k_L;
    for (k_L = 0; k_L < Ntot; k_L++) {
      /* ascending type-II: uniform +1 shift (odd N normalized to type-I -> 0) */
      R omega = K2PI * ((R)(k_L - Ntot / 2 + type_ii)) * x[j];
      f_hat[k_L] += f[j] * (COS(omega) + II * SIN(omega));
    }
  }
}

static void apply_blocked(const problem_nfft *pn) {
  const INT Ntot = Y(problem_nfft_Ntot)((const problem *)pn);
  const INT M = pn->M;
  const R *x_arr = pn->x;
  const C *f_hat = pn->f_hat;
  C *f = pn->f;
  const int type_ii = (pn->variant[0] == NFFT_NDFT_TYPE_II);
  const INT B = NDFT_RECURRENCE_BLOCK;
  INT j;
  for (j = 0; j < M; j++) {
    C v = K(0.0);
    const R x = x_arr[j];
    const R dphi = K2PI * x;                 /* |dphi| <= pi */
    const C dw = COS(dphi) - II * SIN(dphi); /* per-step phase exp(-i 2pi x) */
    INT k_L = 0;
    while (k_L < Ntot) {
      const R omega = ndft_reduced_omega((R)(k_L - Ntot / 2 + type_ii), x);
      C w = COS(omega) - II * SIN(omega);
      INT kend = k_L + B;
      if (kend > Ntot)
        kend = Ntot;
      for (; k_L < kend; k_L++) {
        v += f_hat[k_L] * w;
        w *= dw;
      }
    }
    f[j] = v;
  }
}

static void apply_blocked_adjoint(const problem_nfft *pn) {
  const INT Ntot = Y(problem_nfft_Ntot)((const problem *)pn);
  const INT M = pn->M;
  const R *x_arr = pn->x;
  const C *f = pn->f;
  C *f_hat = pn->f_hat;
  const int type_ii = (pn->variant[0] == NFFT_NDFT_TYPE_II);
  const INT B = NDFT_RECURRENCE_BLOCK;
  memset(f_hat, 0, (size_t)Ntot * sizeof(C));
  {
    INT j;
    for (j = 0; j < M; j++) {
      const R x = x_arr[j];
      const R dphi = K2PI * x;
      const C dw = COS(dphi) + II * SIN(dphi);
      INT k_L = 0;
      while (k_L < Ntot) {
        const R omega = ndft_reduced_omega((R)(k_L - Ntot / 2 + type_ii), x);
        C w = COS(omega) + II * SIN(omega);
        INT kend = k_L + B;
        if (kend > Ntot)
          kend = Ntot;
        for (; k_L < kend; k_L++) {
          f_hat[k_L] += f[j] * w;
          w *= dw;
        }
      }
    }
  }
}

/* Cost model (1D, Ntot = N_total).  C_TRIG matches the legacy direct estimate
 * scale (two transcendentals per term); blocked replaces the per-term pair with
 * ~2/B transcendentals plus one complex multiply per term, so it is strictly
 * below plain for B = 32 and wins on estimate. */
#define NDFT_COST_TRIG 50
#define NDFT_COST_MUL 6

double Y(nfft_ndft_plain_pcost)(const problem *p) {
  const problem_nfft *ego = (const problem_nfft *)p;
  double Ntot = (double)Y(problem_nfft_Ntot)(p);
  double M = (double)ego->M;
  return NDFT_COST_TRIG * Ntot * M;
}

static double pcost_blocked(const problem *p) {
  const problem_nfft *ego = (const problem_nfft *)p;
  double Ntot = (double)Y(problem_nfft_Ntot)(p);
  double M = (double)ego->M;
  return (NDFT_COST_TRIG / (double)NDFT_RECURRENCE_BLOCK + NDFT_COST_MUL) *
         Ntot * M;
}

plan *Y(nfft_ndft_make_plan)(double pcost,
                             void (*apply_native_fwd)(const problem_nfft *),
                             void (*apply_native_adj)(const problem_nfft *),
                             const char *reg_nam) {
  ndft_plan *pln =
      (ndft_plan *)Y(plan_create)(sizeof(ndft_plan), &ndft_plan_adt);
  pln->super.pcost = pcost;
  pln->apply_native_fwd = apply_native_fwd;
  pln->apply_native_adj = apply_native_adj;
  pln->reg_nam = reg_nam;
  return &pln->super;
}

/* d == 1 only; */
static plan *mkplan_ndft_plain(const solver *ego, const problem *p, planner *pl) {
  const problem_nfft *pn = (const problem_nfft *)p;
  (void)ego;
  if (p->adt->kind != NFFT_PROBLEM_NFFT)
    return 0;
  if (pn->sz->rnk != 1)
    return 0;
  if (Y(problem_nfft_has_unit_axis)(p))
    return 0;
  if (PLNR_L(pl) & PLNR_NO_DIRECT)
    return 0;
  if (PLNR_L(pl) & PLNR_NO_NDFT_PLAIN)
    return 0;
  return Y(nfft_ndft_make_plan)(Y(nfft_ndft_plain_pcost)(p),
                                apply_plain, apply_plain_adjoint,
                                "nfft_solver_ndft_1d");
}

static plan *mkplan_ndft_blocked(const solver *ego, const problem *p, planner *pl) {
  const problem_nfft *pn = (const problem_nfft *)p;
  (void)ego;
  if (p->adt->kind != NFFT_PROBLEM_NFFT)
    return 0;
  if (pn->sz->rnk != 1)
    return 0;
  if (Y(problem_nfft_has_unit_axis)(p))
    return 0;
  if (PLNR_L(pl) & PLNR_NO_DIRECT)
    return 0;
  if (PLNR_L(pl) & PLNR_NO_NDFT_BLOCKED)
    return 0;
  return Y(nfft_ndft_make_plan)(pcost_blocked(p),
                                apply_blocked, apply_blocked_adjoint,
                                "nfft_solver_ndft_1d_blocked");
}

static const solver_adt ndft_plain_adt = {NFFT_PROBLEM_NFFT, 0,
                                          mkplan_ndft_plain};
static const solver_adt ndft_blocked_adt = {NFFT_PROBLEM_NFFT, 0,
                                            mkplan_ndft_blocked};

void Y(nfft_solver_ndft_1d_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &ndft_plain_adt));
}
void Y(nfft_solver_ndft_1d_blocked_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &ndft_blocked_adt));
}
