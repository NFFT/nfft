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

/* Threaded CONV solver (Step C), covering every rank.
 *
 * Forward is a gather: nodes are independent, so it is an omp for over j.
 * Adjoint is a scatter, and the tap blocks of different nodes overlap; it is
 * distributed over contiguous bands of axis 0 of the grid instead, so that
 * every grid cell is written by exactly one band and no atomic is needed.
 * Which band a node belongs to comes from a counting sort built beside u in
 * awake. Legacy reaches the same node list with a radix sort of the full
 * linearised window start plus a binary search; the counting sort by band is
 * the part of that the adjoint actually needs.
 *
 * The arithmetic is the serial solvers' own, via the shared kernels in
 * kernel/conv/conv.h, so both directions reproduce the serial result to the
 * last bit at any thread count.
 *
 * The sorted instance builds u and psi in key order and reaches f through
 * perm. Its band lists then hold rows, and the merge still gives each cell
 * ascending rows, which is the serial sorted solver's order. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "conv.h"

#include <math.h> /* pow */

/* No overhead constant; see the DECONV solver's note. The straddle is counted
 * honestly in `other` below. */

typedef struct {
  plan super;
  int d;
  INT *n, *N; /* owned, length d */
  INT M, ntot;
  int m, window, analytic;
  int nthr; /* team size */
  int nb; /* bands along axis 0, 1 <= nb <= nthr */
  INT band; /* band width in grid rows; band >= 2m+2 whenever nb > 1 */
  const R *x; /* borrowed alias of the problem's nodes */
  R *psi; /* M*d*(2m+2); NULL while SLEEPY */
  INT *u; /* M*d; NULL while SLEEPY */
  INT *bstart; /* nb+1 offsets into border; NULL while SLEEPY */
  INT *bcur; /* nb scratch cursors for the counting sort; NULL while SLEEPY */
  INT *border; /* M rows, ascending within each band; NULL while SLEEPY */
  INT *perm; /* M, node of row k in key order; NULL while SLEEPY and for
              * caller order */
  int sorted; /* this plan builds perm */
} conv_omp_plan;

/* ponytail: the rank switch is per node, not hoisted out of the loop. At rank 1
 * that is a few percent; above it, nothing. Hoist it into four loop bodies if a
 * measurement says it matters. */
static inline C node(const conv_omp_plan *pln, C *g, INT j, INT lo, INT hi,
                     C fj, int forward)
{
  const int d = pln->d;
  const INT len = 2 * (INT)pln->m + 2;
  const R *ps = pln->psi + (size_t)j * (size_t)d * (size_t)len;
  const INT *uu = pln->u + (size_t)j * (size_t)d;
  switch (d) {
    case 1:
      return Y(
           conv_node_1d)(g, ps, uu[0], pln->n[0], len, lo, hi, fj, forward);
    case 2:
      return Y(
           conv_node_2d)(g, ps, ps + len, uu[0], uu[1], pln->n[0], pln->n[1],
                          len, lo, hi, fj, forward);
    case 3:
      return Y(
           conv_node_3d)(g, ps, ps + len, ps + 2 * len, uu[0], uu[1], uu[2],
                          pln->n[0], pln->n[1], pln->n[2], len, lo, hi, fj,
                          forward);
    default:
      return Y(conv_node_nd)(g, ps, uu, pln->n, d, len, lo, hi, fj, forward);
  }
}

static void fill_u(conv_omp_plan *pln)
{
  Y(conv_fill_u)(pln->x, pln->perm, pln->n, pln->d, pln->m, pln->M, pln->u);
}

static void fill_psi(conv_omp_plan *pln)
{
  Y(conv_fill_psi)(pln->window, pln->analytic, pln->x, pln->perm, pln->n, pln->N,
                pln->d, pln->m, pln->M, pln->psi);
}

/* Counting sort of the rows by the band their window start falls in. Stable,
 * so each band's list is ascending in row -- which is what makes the
 * adjoint's accumulation order independent of the band count. */
static void fill_bands(conv_omp_plan *pln)
{
  const int d = pln->d, nb = pln->nb;
  const INT M = pln->M, band = pln->band;
  INT j, s;
  int b;

  for (b = 0; b <= nb; b++)
    pln->bstart[b] = 0;
  for (j = 0; j < M; j++) {
    b = (int)(pln->u[j * d] / band);
    if (b >= nb)
      b = nb - 1; /* the last band absorbs the rows the division leaves over */
    pln->bstart[b + 1]++;
  }
  for (b = 0, s = 0; b < nb; b++) {
    s += pln->bstart[b + 1];
    pln->bstart[b + 1] = s;
    pln->bcur[b] = pln->bstart[b];
  }
  for (j = 0; j < M; j++) {
    b = (int)(pln->u[j * d] / band);
    if (b >= nb)
      b = nb - 1;
    pln->border[pln->bcur[b]++] = j;
  }
}

/* The tables live between a wake and the return to SLEEPY (ADR-0009). Window
 * starts and bands are indices: they steer which grid cells apply touches, so
 * they hold true values at AWAKE_ZERO too. psi costs a window evaluation per
 * tap and changes no address, so AWAKE_ZERO gets zeros. */
static void awake(plan *ego_, int wakefulness)
{
  conv_omp_plan *pln = (conv_omp_plan *)ego_;
  const size_t M = (size_t)pln->M, Md = M * (size_t)pln->d;
  const size_t len = (size_t)(2 * pln->m + 2), nb = (size_t)pln->nb;
  if (wakefulness == PLNR_SLEEPY) {
    Y(free)(pln->border);
    Y(free)(pln->bcur);
    Y(free)(pln->bstart);
    Y(free)(pln->u);
    Y(free)(pln->psi);
    Y(free)(pln->perm);
    pln->border = pln->bcur = pln->bstart = pln->u = 0;
    pln->psi = 0;
    pln->perm = 0;
    return;
  }
  A(pln->psi == 0 && pln->u == 0 && pln->border == 0 && pln->perm == 0);
  if (pln->sorted) {
    pln->perm = (INT *)Y(malloc)(M * sizeof(INT));
    Y(conv_sort_nodes)(pln->x, pln->n, pln->d, pln->M, pln->perm);
  }
  pln->psi = (R *)Y(malloc)(Md * len * sizeof(R));
  pln->u = (INT *)Y(malloc)(Md * sizeof(INT));
  pln->bstart = (INT *)Y(malloc)((nb + 1) * sizeof(INT));
  pln->bcur = (INT *)Y(malloc)(nb * sizeof(INT));
  pln->border = (INT *)Y(malloc)(M * sizeof(INT));
  fill_u(pln);
  fill_bands(pln);
  if (wakefulness == PLNR_AWAKE)
    fill_psi(pln);
  else
    memset(pln->psi, 0, Md * len * sizeof(R));
}

static void apply(const plan *ego_, const problem *p)
{
  const conv_omp_plan *pln = (const conv_omp_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  C *g = (C *)pc->g;
  C *f = pc->f;
  const INT M = pln->M, n0 = pln->n[0];
  const INT *perm = pln->perm;
  INT k;
  /* perm is a permutation, so no two iterations write the same f slot. */
#pragma omp parallel for num_threads(pln->nthr) schedule(static)
  for (k = 0; k < M; k++)
    f[perm ? perm[k] : k] = node(pln, g, k, 0, n0, K(0.0), 1);
}

/* Bands, not nodes. Thread ownership follows the omp for, not the thread id,
 * so a team smaller than nb is still correct: every band is visited exactly
 * once and bands write disjoint rows.
 *
 * A node's support is at most len rows wide and every band is at least len
 * rows wide, so a node reaching band b starts in band b or in its cyclic
 * predecessor. Merging those two lists in ascending row gives each cell the
 * same accumulation order the serial solver produces, whatever nb is. */
static void apply_adjoint(const plan *ego_, const problem *p)
{
  const conv_omp_plan *pln = (const conv_omp_plan *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  C *g = pc->g;
  const C *f = pc->f;
  const INT ntot = pln->ntot, n0 = pln->n[0], band = pln->band;
  const int nb = pln->nb;

#pragma omp parallel num_threads(pln->nthr)
  {
    INT k;
    /* The scatter accumulates into an overlapping, node-dependent set that
     * does not cover the grid, so the whole grid must start zeroed. The
     * implicit barrier is what orders the clear before the writes. */
#pragma omp for schedule(static)
    for (k = 0; k < ntot; k++)
      g[k] = K(0.0);

#pragma omp for schedule(dynamic, 1)
    for (k = 0; k < (INT)nb; k++) {
      const int b = (int)k;
      const int prev = (b == 0) ? nb - 1 : b - 1;
      const INT lo = (INT)b * band;
      const INT hi = (b == nb - 1) ? n0 : lo + band;
      INT ia = pln->bstart[prev], ea = pln->bstart[prev + 1];
      INT ib = pln->bstart[b], eb = pln->bstart[b + 1];
      if (nb == 1) /* one band: its own predecessor, so take the list once */
        ia = ea;
      /* Two ascending lists merged into one ascending pass over rows. */
      while (ia < ea || ib < eb) {
        INT j;
        if (ib >= eb || (ia < ea && pln->border[ia] < pln->border[ib]))
          j = pln->border[ia++];
        else
          j = pln->border[ib++];
        (void)node(pln, g, j, lo, hi, f[pln->perm ? pln->perm[j] : j], 0);
      }
    }
  }
}

static void print(const plan *ego_, printer *pr)
{
  const conv_omp_plan *pln = (const conv_omp_plan *)ego_;
  pr->print(pr, "(conv_solver_threaded%s nthr=%d nb=%d)",
            pln->sorted ? " sorted" : "", pln->nthr, pln->nb);
}

static void destroy(plan *ego_)
{
  conv_omp_plan *pln = (conv_omp_plan *)ego_;
  Y(free)(pln->N);
  Y(free)(pln->n);
}

static const plan_adt conv_omp_plan_adt = {apply, awake, print, destroy,
                                           apply_adjoint};

static plan *mkplan_conv_omp(const solver *ego_, const problem *p, planner *pl)
{
  const conv_solver *ego = (const conv_solver *)ego_;
  const problem_conv *pc = (const problem_conv *)p;
  conv_omp_plan *pln;
  const INT len = 2 * (INT)pc->m + 2;
  int t, d;

  if (p->adt->kind != NFFT_PROBLEM_CONV)
    return 0;
  if (pl->nthr < 2)
    return 0;
  d = pc->sz->rnk;
  if (d < 1)
    return 0;
  if (pc->window < NFFT_WINDOW_KAISER_BESSEL
      || pc->window > NFFT_WINDOW_SINC_POWER)
    return 0;
  if (ego->sorted && pc->x_sorted)
    return 0; /* already in key order: nothing to sort */

  pln =
       (conv_omp_plan *)Y(
            plan_create)(sizeof(conv_omp_plan), &conv_omp_plan_adt);
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
  pln->window = pc->window;
  pln->analytic = pc->analytic;
  pln->x = pc->x;
  pln->nthr = pl->nthr;
  /* A node's support is len rows wide, so a band narrower than that could be
   * reached from more than one band away and the adjoint's two-list merge
   * would not be exhaustive. Below two bands the adjoint runs on one thread;
   * the forward still uses the whole team. */
  /* Four bands per thread over a dynamic schedule. One band per thread starves
   * whenever nodes cluster -- polar and Radon geometries concentrate near the
   * origin, MRI trajectories along spokes -- and a band holding most of the
   * nodes serialises the adjoint whatever the team size. Oversubscription does
   * not touch the exactness argument, which never mentions the band count. */
  pln->nb = (int)(pln->n[0] / len);
  if (pln->nb > 4 * pln->nthr)
    pln->nb = 4 * pln->nthr;
  if (pln->nb < 1)
    pln->nb = 1;
  pln->band = pln->n[0] / pln->nb;
  pln->psi = 0;
  pln->u = pln->bstart = pln->bcur = pln->border = 0;
  pln->perm = 0;
  pln->sorted = ego->sorted;
  /* One plan serves both directions and they scale differently -- the forward
   * over the whole team, the adjoint over nb bands -- so the model is the mean
   * of the two. */
  /* The serial solver's counts, plus the straddle: a node crossing a band
   * boundary pays its per-node run setup twice. Its tap work is split, not
   * repeated, so only `other` moves. */
  {
    const double taps = (double)pln->M * pow((double)len, (double)d);
    const double straddle = (double)pln->nb * (double)len / (double)pln->n[0];
    Y(ops_zero)(&pln->super.ops);
    pln->super.ops.mul = 2.0 * taps;
    pln->super.ops.add = 2.0 * taps;
    pln->super.ops.other =
         2.0 * taps + 2.0 * (double)pln->M * (double)d * (1.0 + straddle);
  }
  Y(conv_order_ops)(&pln->super.ops, pln->M, d, pln->m, pln->ntot, ego->sorted,
                 pc->x_sorted);
  return &pln->super;
}

static const solver_adt conv_omp_adt = {NFFT_PROBLEM_CONV, 0, mkplan_conv_omp};

void Y(conv_solver_threaded_register)(planner *pl)
{
  Y(conv_register_both)(pl, &conv_omp_adt);
}
