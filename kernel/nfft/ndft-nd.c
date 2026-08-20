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

#include <string.h> /* memset -- adjoint zeroes f_hat before accumulating */

/* Direct NDFT, generic rank (d >= 2). 
 *
 * Two changes compared to legacy code:
 *
 * (1) Correctly handle unit axes: the legacy code carries on `== N[t]/2 - 1`. For a
 * unit axis (N[t] == 1) the range [-N/2, N/2-1] is [0, -1], so k[t] starts at
 * 0 and `== -1` never fires. `>=` makes a unit axis carry transparently.
 *
 * (2) Support for odd bandwidth: threshold `(N[t]-1)/2` instead of `N[t]/2 - 1`: 
 * the two are identical for even N (integer division), so the even path is
 * unchanged. For odd N, `(N[t]-1)/2` is the correct upper bound.
 *
 * Per-axis type-II: the odometer still iterates the raw index range
 * k[t] in [-N[t]/2, (N[t]-1)/2] (the carry logic above is untouched); a
 * type-II axis instead contributes the effective frequency `k[t] + isII[t]`
 * to the phase, a uniform +1 shift that turns the type-I range into the
 * ascending type-II range (+N[t]/2 at the last slot).  `isII[t]` is read once
 * per call from `pn->variant[t]`; it is 0 for type-I axes and for odd-N axes. */

static void apply(const problem_nfft *pn) {
  const problem *p = (const problem *)pn;
  const int d = pn->sz->rnk;
  const INT Ntot = Y(problem_nfft_Ntot)(p);
  const INT M = pn->M;
  const R *xin = pn->x;
  const C *f_hat = pn->f_hat;
  C *f = pn->f;
  INT N[d];
  int isII[d];
  int tt;
  for (tt = 0; tt < d; tt++) {
    N[tt] = Y(problem_nfft_N)(p, tt);
    isII[tt] = (pn->variant[tt] == NFFT_NDFT_TYPE_II); /* even only; odd->0 */
  }
  {
    INT j;
    for (j = 0; j < M; j++) {
      C v = K(0.0);
      R x_[d], omega, Omega[d + 1];
      INT t, t2, k_L, k[d];
      Omega[0] = K(0.0);
      for (t = 0; t < d; t++) {
        k[t] = -N[t] / 2;
        x_[t] = K2PI * xin[j * d + t];
        Omega[t + 1] = ((R)(k[t] + isII[t])) * x_[t] + Omega[t];
      }
      omega = Omega[d];
      for (k_L = 0; k_L < Ntot; k_L++) {
        v += f_hat[k_L] * (COS(omega) - II * SIN(omega));
        for (t = d - 1; (t >= 1) && (k[t] >= (N[t] - 1) / 2); t--)
          k[t] -= N[t] - 1;
        k[t]++;
        for (t2 = t; t2 < d; t2++)
          Omega[t2 + 1] = ((R)(k[t2] + isII[t2])) * x_[t2] + Omega[t2];
        omega = Omega[d];
      }
      f[j] = v;
    }
  }
}

static void apply_adjoint(const problem_nfft *pn) {
  const problem *p = (const problem *)pn;
  const int d = pn->sz->rnk;
  const INT Ntot = Y(problem_nfft_Ntot)(p);
  const INT M = pn->M;
  const R *xin = pn->x;
  const C *f = pn->f;
  C *f_hat = pn->f_hat;
  INT N[d];
  int isII[d];
  int tt;
  for (tt = 0; tt < d; tt++) {
    N[tt] = Y(problem_nfft_N)(p, tt);
    isII[tt] = (pn->variant[tt] == NFFT_NDFT_TYPE_II); /* even only; odd->0 */
  }
  memset(f_hat, 0, (size_t)Ntot * sizeof(C));
  {
    INT j;
    for (j = 0; j < M; j++) {
      R x_[d], omega, Omega[d + 1];
      INT t, t2, k_L, k[d];
      const C fj = f[j];
      Omega[0] = K(0.0);
      for (t = 0; t < d; t++) {
        k[t] = -N[t] / 2;
        x_[t] = K2PI * xin[j * d + t];
        Omega[t + 1] = ((R)(k[t] + isII[t])) * x_[t] + Omega[t];
      }
      omega = Omega[d];
      for (k_L = 0; k_L < Ntot; k_L++) {
        f_hat[k_L] += fj * (COS(omega) + II * SIN(omega));
        for (t = d - 1; (t >= 1) && (k[t] >= (N[t] - 1) / 2); t--)
          k[t] -= N[t] - 1;
        k[t]++;
        for (t2 = t; t2 < d; t2++)
          Omega[t2 + 1] = ((R)(k[t2] + isII[t2])) * x_[t2] + Omega[t2];
        omega = Omega[d];
      }
    }
  }
}

/* rnk >= 2 only */
static plan *mkplan_ndft_nd(const solver *ego, const problem *p, planner *pl) {
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
  if (PLNR_L(pl) & PLNR_NO_NDFT_PLAIN)
    return 0;
  return Y(nfft_ndft_make_plan)(Y(nfft_ndft_plain_pcost)(p), apply,
                                apply_adjoint, "nfft_solver_ndft_nd");
}

static const solver_adt ndft_nd_adt = {NFFT_PROBLEM_NFFT, 0, mkplan_ndft_nd};

void Y(nfft_solver_ndft_nd_register)(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &ndft_nd_adt));
}
