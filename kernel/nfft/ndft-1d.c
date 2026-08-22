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

#include <string.h> /* memset - adjoint zeroes f_hat before accumulating */

/* Direct NDFT, rank 1. The Ntot frequencies run from k = -Ntot/2 upwards: the
 * symmetric range for odd Ntot, the type-I range for even Ntot. A type-II axis
 * shifts every frequency by +1, which turns that into the ascending type-II
 * range (+Ntot/2 at the last slot).
 *
 * Holds the plan type both direct NDFT solvers share. */

typedef struct
{
  plan super;
  void (*fwd)(const problem_nfft *pn);
  void (*adj)(const problem_nfft *pn);
  const char *reg_nam;
} ndft_plan;

static void apply(const plan *ego, const problem *p) {
  const ndft_plan *pln = (const ndft_plan *)ego;
  pln->fwd((const problem_nfft *)p);
}
static void apply_adjoint(const plan *ego, const problem *p) {
  const ndft_plan *pln = (const ndft_plan *)ego;
  pln->adj((const problem_nfft *)p);
}

static void print(const plan *ego, printer *pr) {
  const ndft_plan *pln = (const ndft_plan *)ego;
  pr->print(pr, "(%s pcost=%D)", pln->reg_nam, (INT)pln->super.pcost);
}

/* Owns only the base allocation: no awake hook, no destroy. */
static const plan_adt ndft_plan_adt = {apply, 0, print, 0, apply_adjoint};

static void sum(const problem_nfft *pn) {
  const INT Ntot = Y(problem_nfft_Ntot)((const problem *)pn);
  const INT M = pn->M;
  const R *x_arr = pn->x;
  const C *f_hat = pn->f_hat;
  C *f = pn->f;
  const R k0 = (R)(-Ntot / 2 + (pn->variant[0] == NFFT_NDFT_TYPE_II));
  const INT B = NDFT_RECURRENCE_BLOCK;
  INT j;
  for (j = 0; j < M; j++) {
    C v = K(0.0);
    const R x = x_arr[j];
    const R dphi = K2PI * x; /* |dphi| <= pi */
    const C dw = COS(dphi) - II * SIN(dphi);
    INT k_L = 0;
    while (k_L < Ntot) {
      const R omega = Y(nfft_ndft_reduced_omega)(k0 + (R)k_L, x, K(0.0));
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

static void sum_adjoint(const problem_nfft *pn) {
  const INT Ntot = Y(problem_nfft_Ntot)((const problem *)pn);
  const INT M = pn->M;
  const R *x_arr = pn->x;
  const C *f = pn->f;
  C *f_hat = pn->f_hat;
  const R k0 = (R)(-Ntot / 2 + (pn->variant[0] == NFFT_NDFT_TYPE_II));
  const INT B = NDFT_RECURRENCE_BLOCK;
  INT j;
  memset(f_hat, 0, (size_t)Ntot * sizeof(C));
  for (j = 0; j < M; j++) {
    const R x = x_arr[j];
    const R dphi = K2PI * x;
    const C dw = COS(dphi) + II * SIN(dphi);
    INT k_L = 0;
    while (k_L < Ntot) {
      const R omega = Y(nfft_ndft_reduced_omega)(k0 + (R)k_L, x, K(0.0));
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

/* NDFT_COST_TRIG is the legacy direct estimate scale, two transcendentals per
 * term; the recurrence spends ~2/B of those plus one complex multiply. */
#define NDFT_COST_TRIG 50
#define NDFT_COST_MUL 6

double Y(nfft_ndft_pcost)(const problem *p) {
  const problem_nfft *ego = (const problem_nfft *)p;
  double Ntot = (double)Y(problem_nfft_Ntot)(p);
  double M = (double)ego->M;
  return (NDFT_COST_TRIG / (double)NDFT_RECURRENCE_BLOCK + NDFT_COST_MUL) *
         Ntot * M;
}

plan *Y(nfft_ndft_make_plan)(double pcost,
                             void (*fwd)(const problem_nfft *),
                             void (*adj)(const problem_nfft *),
                             const char *reg_nam) {
  ndft_plan *pln =
      (ndft_plan *)Y(plan_create)(sizeof(ndft_plan), &ndft_plan_adt);
  pln->super.pcost = pcost;
  pln->fwd = fwd;
  pln->adj = adj;
  pln->reg_nam = reg_nam;
  return &pln->super;
}

/* rnk 1 only */
static plan *mkplan(const solver *ego, const problem *p, planner *pl) {
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
  return Y(nfft_ndft_make_plan)(Y(nfft_ndft_pcost)(p), sum, sum_adjoint,
                                "nfft_solver_ndft_1d");
}

static const solver_adt ndft_1d_adt = {NFFT_PROBLEM_NFFT, 0, mkplan};

void Y(nfft_solver_ndft_1d_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &ndft_1d_adt));
}
