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

/* Threaded DECONV solver (Step A). It lives in the add-on library, so a
 * program that does not link it has no threaded solver and cannot raise the
 * planner's thread count.
 *
 * Every rank >= 1: rank 1 distributes its single row's frequencies directly,
 * higher ranks distribute whole rows (one contiguous run along the last
 * axis each). Threads write disjoint grid cells, so there is no reduction
 * and no atomic, and results do not depend on thread count. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "deconv.h"

/* The row loop holds the leading axes' multi index in a fixed-size automatic
 * array, so that each thread's scratch is a few words on its own stack rather
 * than a variable-length array inside a parallel loop. 64 INTs is a few
 * hundred bytes of thread-local stack, cheap at any team size; module-local,
 * not NFFT_-prefixed, since it names no public API. */
#define THREADS_MAX_RANK 64

typedef struct {
  plan super;
  int d; /* rank, >= 1 */
  INT *N, *n, *Nneg; /* owned, length d: geometry and per-axis slot split */
  INT ntot; /* product of n[t]: size of the oversampled grid */
  INT rows; /* product of N[t] over t < d-1; 1 when d == 1 */
  int m, window;
  int nthr; /* team size: the budget, clamped to the row count */
  R **phi_hut_inv; /* owned array of d owned tables, phi_hut_inv[t][0..N[t]) */
  int level; /* content of the tables: SLEEPY, AWAKE_ZERO or AWAKE */
} deconv_omp_plan;

/* AWAKE_ZERO must cost no window evaluation, so the tables get placeholder
 * zeros. AWAKE_ZERO reached by downgrade only drops the level, so a later
 * upgrade refills. */
static void awake(plan *ego_, int wakefulness)
{
  deconv_omp_plan *pln = (deconv_omp_plan *)ego_;
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

static void run(const deconv_omp_plan *pln, const problem_deconv *pd,
                int forward)
{
  const int d = pln->d;
  const INT *N = pln->N, *n = pln->n, *Nneg = pln->Nneg;
  R *const *phi = (R *const *)pln->phi_hut_inv;
  const INT Nl = N[d - 1], nl = n[d - 1], Nnegl = Nneg[d - 1];
  const INT Nposl = Nl - Nnegl;
  const R *pw = phi[d - 1];
  const INT rows = pln->rows, ntot = pln->ntot;
  C *f_hat = pd->f_hat;
  C *g = pd->g;

#pragma omp parallel num_threads(pln->nthr)
  {
    INT i;

    /* omp for's implicit barrier orders this clear before the write below --
     * do not add nowait. */
    if (forward) {
      if (d == 1) {
#pragma omp for schedule(static)
        for (i = Nposl; i < nl - Nnegl; i++)
          g[i] = K(0.0);
      } else {
#pragma omp for schedule(static)
        for (i = 0; i < ntot; i++)
          g[i] = K(0.0);
      }
    }

    if (d == 1) {
      /* A single row: distribute its frequencies instead. */
#pragma omp for schedule(static)
      for (i = 0; i < Nl; i++) {
        const INT pos = (i < Nnegl) ? nl - Nnegl + i : i - Nnegl;
        if (forward)
          g[pos] = f_hat[i] * pw[i];
        else
          f_hat[i] = g[pos] * pw[i];
      }
    } else {
#pragma omp for schedule(static)
      for (i = 0; i < rows; i++) {
        INT ks[THREADS_MAX_RANK];
        R pre = K(1.0);
        INT fh_off = 0, g_off = 0, q = i, k;
        int t;
        /* Recover the row's multi index over axes 0..d-2. d-1 divisions per
         * row, against the legacy D step's d per element. */
        for (t = d - 2; t >= 0; t--) {
          ks[t] = q % N[t];
          q /= N[t];
        }
        /* Fold the leading axes' window factors and offsets once for the whole
         * row. A leading slot below Nneg[t] belongs to the grid's tail. */
        for (t = 0; t < d - 1; t++) {
          const INT kt = ks[t] < Nneg[t] ? n[t] - Nneg[t] + ks[t] :
                                           ks[t] - Nneg[t];
          pre *= phi[t][ks[t]];
          fh_off = fh_off * N[t] + ks[t];
          g_off = g_off * n[t] + kt;
        }
        fh_off *= Nl;
        g_off *= nl;
        {
          /* The last axis splits into two contiguous runs: slots below Nnegl
           * land at the grid's tail, the rest at its head. */
          C *ga = g + g_off + nl - Nnegl, *gb = g + g_off;
          C *fa = f_hat + fh_off, *fb = f_hat + fh_off + Nnegl;
          if (forward) {
            for (k = 0; k < Nnegl; k++)
              ga[k] = fa[k] * (pre * pw[k]);
            for (k = 0; k < Nposl; k++)
              gb[k] = fb[k] * (pre * pw[Nnegl + k]);
          } else {
            for (k = 0; k < Nnegl; k++)
              fa[k] = ga[k] * (pre * pw[k]);
            for (k = 0; k < Nposl; k++)
              fb[k] = gb[k] * (pre * pw[Nnegl + k]);
          }
        }
      }
    }
  }
}

static void apply(const plan *ego_, const problem *p)
{
  run((const deconv_omp_plan *)ego_, (const problem_deconv *)p, 1);
}

static void apply_adjoint(const plan *ego_, const problem *p)
{
  run((const deconv_omp_plan *)ego_, (const problem_deconv *)p, 0);
}

static void print(const plan *ego_, printer *pr)
{
  const deconv_omp_plan *pln = (const deconv_omp_plan *)ego_;
  pr->print(pr, "(deconv_solver_threaded nthr=%d)", pln->nthr);
}

static void destroy(plan *ego_)
{
  deconv_omp_plan *pln = (deconv_omp_plan *)ego_;
  int t;
  for (t = 0; t < pln->d; t++)
    Y(free)(pln->phi_hut_inv[t]);
  Y(free)(pln->phi_hut_inv);
  Y(free)(pln->Nneg);
  Y(free)(pln->n);
  Y(free)(pln->N);
}

static const plan_adt deconv_omp_plan_adt = {apply, awake, print, destroy,
                                             apply_adjoint};

static plan *mkplan_deconv_omp(const solver *ego, const problem *p,
                               planner *pl)
{
  const problem_deconv *pd = (const problem_deconv *)p;
  deconv_omp_plan *pln;
  int t, d;
  INT rows = 1;
  INT Ntot = 1;
  (void)ego;

  if (p->adt->kind != NFFT_PROBLEM_DECONV)
    return 0;
  /* A one-thread budget has nothing to distribute. Declining keeps linking the
   * add-on from changing any plan until the caller explicitly asks for
   * threads. */
  if (pl->nthr < 2)
    return 0;
  d = pd->sz->rnk;
  if (d < 1 || d > THREADS_MAX_RANK)
    return 0;
  if (pd->window < NFFT_WINDOW_KAISER_BESSEL
      || pd->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* reject Dirac or other invalid ordinals */
  for (t = 0; t < d; t++)
    if (Y(problem_deconv_n)(p, t) < Y(problem_deconv_N)(p, t))
      return 0; /* n < N aliases grid cells */

  pln = (deconv_omp_plan *)Y(plan_create)(sizeof(deconv_omp_plan),
                                       &deconv_omp_plan_adt);
  pln->d = d;
  pln->N = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->n = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->Nneg = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->phi_hut_inv = (R **)Y(malloc)((size_t)d * sizeof(R *));
  pln->ntot = Y(problem_deconv_ntot)(p);
  for (t = 0; t < d; t++) {
    pln->N[t] = Y(problem_deconv_N)(p, t);
    pln->n[t] = Y(problem_deconv_n)(p, t);
    pln->Nneg[t] = pln->N[t] / 2
                   - (Y(problem_deconv_variant)(p, t) == NFFT_NDFT_TYPE_II ?
                           (INT)1 :
                           (INT)0);
    pln->phi_hut_inv[t] = (R *)Y(malloc)((size_t)pln->N[t] * sizeof(R));
    Ntot *= pln->N[t];
    if (t < d - 1)
      rows *= pln->N[t];
  }
  pln->rows = rows;
  pln->m = pd->m;
  pln->window = pd->window;
  /* plan_with_nthreads sets a maximum, not a target: a team larger than the
   * work available would only add idle members. */
  pln->nthr = pl->nthr;
  if (d > 1 && (INT)pln->nthr > rows)
    pln->nthr = (int)rows;
  pln->level = PLNR_SLEEPY;
  /* Copy of the serial mkplan's ops formula (kernel/deconv/deconv-1d.c,
   * deconv-nd.c): an op count cannot express a thread count. */
  Y(ops_zero)(&pln->super.ops);
  pln->super.ops.mul = (d == 1) ? 2.0 * (double)Ntot : 3.0 * (double)Ntot;
  pln->super.ops.other =
       (d == 1) ? (double)(pln->n[0] - pln->N[0]) + 2.0 * (double)Ntot :
                  (double)pln->ntot + 2.0 * (double)Ntot
                       + 2.0 * (double)(d - 1) * (double)rows;
  return &pln->super;
}

static const solver_adt deconv_omp_adt = {NFFT_PROBLEM_DECONV, 0,
                                          mkplan_deconv_omp};

void Y(deconv_solver_threaded_register)(planner *pl)
{
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &deconv_omp_adt));
}
