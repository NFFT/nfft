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

/* Direct NDFT, generic rank (d >= 2). An odometer over the leading d-1 axes
 * supplies each row's base phase; the innermost axis varies fastest and is
 * blocked (ndft.h). Axis t runs k[t] in [-N[t]/2, (N[t]-1)/2]; a type-II axis
 * contributes k[t] + 1, the uniform shift that turns the type-I range into the
 * ascending type-II range.
 *
 * The carry differs from the legacy odometer in two places, both needed for
 * the geometries this API admits: `>=` rather than `==`, so a unit axis
 * (range [0, -1]) carries at all, and the upper bound `(N[t]-1)/2` rather than
 * `N[t]/2 - 1`, which is the same for even N and correct for odd N. */

static void sum(const problem_nfft *pn) {
  const problem *p = (const problem *)pn;
  const int d = pn->sz->rnk;
  const INT Ntot = Y(problem_nfft_Ntot)(p);
  const INT M = pn->M;
  const R *xin = pn->x;
  const C *f_hat = pn->f_hat;
  C *f = pn->f;
  const INT B = NDFT_RECURRENCE_BLOCK;
  INT N[d];
  int isII[d];
  int tt;
  for (tt = 0; tt < d; tt++) {
    N[tt] = Y(problem_nfft_N)(p, tt);
    isII[tt] = (pn->variant[tt] == NFFT_NDFT_TYPE_II); /* even only; odd->0 */
  }
  {
    const INT Nin = N[d - 1]; /* innermost axis */
    const R kin0 = (R)(-Nin / 2 + isII[d - 1]); /* its first frequency */
    const INT rows = Ntot / Nin;
    INT j;
    for (j = 0; j < M; j++) {
      C v = K(0.0);
      R x_[d], Omega[d], base;
      INT t, t2, r, k_L, k[d];
      const R x_in = xin[j * d + d - 1];
      const R dphi = K2PI * x_in; /* |dphi| <= pi */
      const C dw = COS(dphi) - II * SIN(dphi);
      Omega[0] = K(0.0);
      for (t = 0; t < d - 1; t++) {
        x_[t] = xin[j * d + t];
        k[t] = -N[t] / 2;
        Omega[t + 1] = ((R)(k[t] + isII[t])) * x_[t] + Omega[t];
      }
      base = Omega[d - 1];
      k_L = 0;
      for (r = 0; r < rows; r++) {
        const R bred = base - RINT(base);
        INT i = 0;
        while (i < Nin) {
          const R omega = Y(nfft_ndft_reduced_omega)(kin0 + (R)i, x_in, bred);
          C w = COS(omega) - II * SIN(omega);
          INT iend = i + B;
          if (iend > Nin)
            iend = Nin;
          for (; i < iend; i++) {
            v += f_hat[k_L++] * w;
            w *= dw;
          }
        }
        for (t = d - 2; (t >= 1) && (k[t] >= (N[t] - 1) / 2); t--)
          k[t] -= N[t] - 1;
        k[t]++;
        for (t2 = t; t2 < d - 1; t2++)
          Omega[t2 + 1] = ((R)(k[t2] + isII[t2])) * x_[t2] + Omega[t2];
        base = Omega[d - 1];
      }
      f[j] = v;
    }
  }
}

static void sum_adjoint(const problem_nfft *pn) {
  const problem *p = (const problem *)pn;
  const int d = pn->sz->rnk;
  const INT Ntot = Y(problem_nfft_Ntot)(p);
  const INT M = pn->M;
  const R *xin = pn->x;
  const C *f = pn->f;
  C *f_hat = pn->f_hat;
  const INT B = NDFT_RECURRENCE_BLOCK;
  INT N[d];
  int isII[d];
  int tt;
  for (tt = 0; tt < d; tt++) {
    N[tt] = Y(problem_nfft_N)(p, tt);
    isII[tt] = (pn->variant[tt] == NFFT_NDFT_TYPE_II);
  }
  memset(f_hat, 0, (size_t)Ntot * sizeof(C));
  {
    const INT Nin = N[d - 1];
    const R kin0 = (R)(-Nin / 2 + isII[d - 1]);
    const INT rows = Ntot / Nin;
    INT j;
    for (j = 0; j < M; j++) {
      R x_[d], Omega[d], base;
      INT t, t2, r, k_L, k[d];
      const C fj = f[j];
      const R x_in = xin[j * d + d - 1];
      const R dphi = K2PI * x_in;
      const C dw = COS(dphi) + II * SIN(dphi);
      Omega[0] = K(0.0);
      for (t = 0; t < d - 1; t++) {
        x_[t] = xin[j * d + t];
        k[t] = -N[t] / 2;
        Omega[t + 1] = ((R)(k[t] + isII[t])) * x_[t] + Omega[t];
      }
      base = Omega[d - 1];
      k_L = 0;
      for (r = 0; r < rows; r++) {
        const R bred = base - RINT(base);
        INT i = 0;
        while (i < Nin) {
          const R omega = Y(nfft_ndft_reduced_omega)(kin0 + (R)i, x_in, bred);
          C w = COS(omega) + II * SIN(omega);
          INT iend = i + B;
          if (iend > Nin)
            iend = Nin;
          for (; i < iend; i++) {
            f_hat[k_L++] += fj * w;
            w *= dw;
          }
        }
        for (t = d - 2; (t >= 1) && (k[t] >= (N[t] - 1) / 2); t--)
          k[t] -= N[t] - 1;
        k[t]++;
        for (t2 = t; t2 < d - 1; t2++)
          Omega[t2 + 1] = ((R)(k[t2] + isII[t2])) * x_[t2] + Omega[t2];
        base = Omega[d - 1];
      }
    }
  }
}

/* rnk >= 2 only */
static plan *mkplan(const solver *ego, const problem *p, planner *pl) {
  const problem_nfft *pn = (const problem_nfft *)p;
  (void)ego;
  if (p->adt->kind != NFFT_PROBLEM_NFFT)
    return 0;
  if (pn->sz->rnk < 2)
    return 0;
  if (Y(problem_nfft_has_unit_axis)(p))
    return 0;
  if (PLNR_L(pl) & PLNR_NO_DIRECT)
    return 0;
  return Y(nfft_ndft_make_plan)(Y(nfft_ndft_pcost)(p), sum, sum_adjoint,
                                "nfft_solver_ndft_nd");
}

static const solver_adt ndft_nd_adt = {NFFT_PROBLEM_NFFT, 0, mkplan};

void Y(nfft_solver_ndft_nd_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &ndft_nd_adt));
}
