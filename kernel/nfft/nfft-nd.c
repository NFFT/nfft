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

/* Fast NFFT solver, composes from DECONV, DFT, and CONV child plans.
 *
 * apply_fwd:    DECONV(f_hat -> g1) -> FFTW forward (g1 -> g2) -> CONV(g2 -> f)
 * apply_adjoint: CONV^H(f -> g2) -> FFTW backward (g2 -> g1) -> DECONV^H(g1 -> f_hat)
 *
 * Ownership: the DECONV/CONV child problems borrow the parent's f_hat/f/x. This
 * plan owns and frees g1, g2, the two FFTW plans, and the two child plans/problems. */

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

/* Applicability guard: every axis of sz must satisfy N_t > m, n_t > 2m + 2 and
 * n_t > N_t. The last is oversampling sigma = n/N strictly above 1; at
 * sigma == 1 there is no zero-pad to deconvolve into, and for even N
 * sinc-power's phi_hut is exactly zero at the band edge. Y(plan_ng_guru)
 * applies the same predicate up front unless the caller passed
 * NFFT_NO_FAST_NATIVE. Unit axes are declined at the mkplan gate below. */
int Y(nfft_fast_guards_ok)(const problem *p, int m)
{
  const problem_nfft *ego = (const problem_nfft *)p;
  int t;
  for (t = 0; t < ego->sz->rnk; ++t) {
    INT Nt = Y(problem_nfft_N)(p, t);
    INT nt = Y(problem_nfft_n)(p, t);
    if (!(Nt > (INT)m))
      return 0;
    if (!(nt > (INT)(2 * m + 2)))
      return 0;
    if (!(nt > Nt)) /* sigma = n/N > 1 */
      return 0;
  }
  return 1;
}

typedef struct {
  plan super;
  int d; /* rank (>=1) */
  INT *Nc, *nc; /* per-axis bandwidth / oversampled size (len d) */
  int *narr; /* per-axis nc as int, for FFTW (len d) */
  INT ntot; /* product of nc[t]; size of g1/g2 */
  C *g1, *g2; /* scratch, size ntot each (FFTW out-of-place) */
  FFTW(plan) pfwd, pback; /* child FFTW plans */
  plan *deconv_child, *conv_child; /* child DECONV and CONV plans. */
  problem *deconv_prob, *conv_prob; /* child problems (borrow parent data) */
} native_fast_plan;

static void apply(const plan *ego_, const problem *p)
{
  native_fast_plan *pln = (native_fast_plan *)ego_;
  const problem_nfft *pn = (const problem_nfft *)p;
  ((problem_deconv *)pln->deconv_prob)->f_hat =
       pn->f_hat; /* DECONV reads f_hat */
  ((problem_conv *)pln->conv_prob)->f = pn->f; /* CONV writes f */
  pln->deconv_child->adt->apply(pln->deconv_child,
                                pln->deconv_prob); /* f_hat -> g1 */
  FFTW(execute)(pln->pfwd); /* g1 -> g2 */
  pln->conv_child->adt->apply(pln->conv_child, pln->conv_prob); /* g2 -> f */
}
static void apply_adjoint(const plan *ego_, const problem *p)
{
  native_fast_plan *pln = (native_fast_plan *)ego_;
  const problem_nfft *pn = (const problem_nfft *)p;
  ((problem_conv *)pln->conv_prob)->f = pn->f; /* CONV^H reads f */
  ((problem_deconv *)pln->deconv_prob)->f_hat =
       pn->f_hat; /* DECONV^H writes f_hat */
  pln->conv_child->adt->apply_adjoint(pln->conv_child,
                                      pln->conv_prob); /* f -> g2 */
  FFTW(execute)(pln->pback); /* g2 -> g1 */
  pln->deconv_child->adt->apply_adjoint(pln->deconv_child,
                                        pln->deconv_prob); /* g1 -> f_hat */
}
static void awake(plan *ego_, int wakefulness)
{
  native_fast_plan *pln = (native_fast_plan *)ego_;
  Y(plan_awake)(pln->deconv_child, wakefulness); /* builds phi_hut; idempotent */
  Y(plan_awake)(pln->conv_child, wakefulness); /* builds psi from owned x */
}
static void destroy(plan *ego_)
{
  native_fast_plan *pln = (native_fast_plan *)ego_;
  if (pln->deconv_child)
    Y(plan_destroy)(pln->deconv_child);
  if (pln->conv_child)
    Y(plan_destroy)(pln->conv_child);
  if (pln->deconv_prob)
    Y(problem_destroy)(pln->deconv_prob);
  if (pln->conv_prob)
    Y(problem_destroy)(pln->conv_prob);
  if (pln->pback)
    FFTW(destroy_plan)(pln->pback);
  if (pln->pfwd)
    FFTW(destroy_plan)(pln->pfwd);
  Y(free)(pln->g2);
  Y(free)(pln->g1);
  Y(free)(pln->narr);
  Y(free)(pln->nc);
  Y(free)(pln->Nc);
}
static void print(const plan *ego_, printer *pr)
{
  const native_fast_plan *pln = (const native_fast_plan *)ego_;
  /* %p prints a plan of this planner, so splice in FFTW's own description. */
  char *fftw_desc = pln->pfwd ? FFTW(sprint_plan)(pln->pfwd) : 0;
  pr->print(pr,
            "(nfft_solver_fast_native pcost=%D"
            "%((deconv %p)%)"
            "%((fftw %s)%)"
            "%((conv %p)%))",
            (INT)pln->super.pcost, (void *)pln->deconv_child, fftw_desc,
            (void *)pln->conv_child);
  if (fftw_desc)
    FFTW(free)(fftw_desc);
}
static const plan_adt native_fast_adt = {apply, awake, print, destroy,
                                         apply_adjoint};

/* Test-only (prototype in iplanner.h): native_fast_plan is file-static, so
 * plan.c reaches its children through this accessor rather than the struct. */
int Y(nfft_fast_native_test_parts)(const plan *top, opcnt *deconv, opcnt *conv,
                                opcnt *fft)
{
  const native_fast_plan *pln = (const native_fast_plan *)top;
  double fa, fm, ff;
  if (!top || top->adt != &native_fast_adt)
    return 0;
  Y(ops_cpy)(&pln->deconv_child->ops, deconv);
  Y(ops_cpy)(&pln->conv_child->ops, conv);
  FFTW(flops)(pln->pfwd, &fa, &fm, &ff);
  Y(ops_zero)(fft);
  fft->add = fa;
  fft->mul = fm;
  fft->fma = ff;
  return 1;
}

/* only (rnk >= 1). This solver runs no loop of its own: apply delegates to the
 * DECONV child, the child FFTW plan and the CONV child. A composer is not a
 * serial solver, so it does not test NO_NONTHREADEDP -- the children do, and
 * a child with no threaded solver is what makes the composition decline. */
static plan *mkplan_native_fast(const solver *ego, const problem *p,
                                planner *pl)
{
  const problem_nfft *pn = (const problem_nfft *)p;
  INT M, ntot;
  int d, t, m, window;
  native_fast_plan *pln;
  (void)ego;
  if (p->adt->kind != NFFT_PROBLEM_NFFT)
    return 0;
  if (pn->sz->rnk < 1)
    return 0; /* rnk >= 1 */
  if (Y(problem_nfft_has_unit_axis)(p))
    return 0;
  if (!Y(nfft_fast_guards_ok)(p, pn->m))
    return 0;
  if (pn->window < NFFT_WINDOW_KAISER_BESSEL
      || pn->window > NFFT_WINDOW_SINC_POWER)
    return 0; /* Kaiser-Bessel/Gaussian/B-spline/sinc only; */
  if (PLNR_L(pl) & PLNR_NO_FAST_NATIVE)
    return 0;

  d = pn->sz->rnk;
  M = pn->M;
  m = pn->m;
  window = pn->window;

  pln = (native_fast_plan *)Y(plan_create)(sizeof(native_fast_plan),
                                        &native_fast_adt);
  pln->d = d;
  pln->Nc = 0;
  pln->nc = 0;
  pln->narr = 0;
  pln->deconv_child = 0;
  pln->conv_child = 0;
  pln->deconv_prob = 0;
  pln->conv_prob = 0;
  pln->pfwd = 0;
  pln->pback = 0;
  /* Per-axis bandwidths / oversampled sizes / int-cast n for FFTW; total size
   * of g1/g2 is the product of the per-axis nc[t]. */
  pln->Nc = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->nc = (INT *)Y(malloc)((size_t)d * sizeof(INT));
  pln->narr = (int *)Y(malloc)((size_t)d * sizeof(int));
  ntot = 1;
  for (t = 0; t < d; t++) {
    pln->Nc[t] = Y(problem_nfft_N)(p, t);
    pln->nc[t] = Y(problem_nfft_n)(p, t);
    pln->narr[t] = (int)pln->nc[t];
    ntot *= pln->nc[t];
  }
  pln->ntot = ntot;
  pln->g1 = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  pln->g2 = (C *)Y(malloc)((size_t)ntot * sizeof(C));
  /* Internal FFTW plans. The user's flags follow FFTW's own vocabulary (0 is
   * FFTW_MEASURE); the FFTs touch only the owned g1/g2, so input preservation
   * is stripped and destruction forced. With FFTW_MEASURE (the default) this
   * planning step actually executes transforms over g1/g2 while timing
   * candidates, so they must already be initialised by the time this runs --
   * see the memset below, placed after this block for exactly that reason. */
  {
    unsigned ff = (pn->fftw_flags & ~(unsigned)FFTW_PRESERVE_INPUT)
                  | (unsigned)FFTW_DESTROY_INPUT;
    /* The child FFT gets the same thread budget as the NFFT plan. FFTW's
     * count is planner-global, so it is restored afterwards; pl->nthr is
     * already part of the wisdom key, so the two plans stay distinguishable.
     * Without the add-on threading library both hooks are null and the child
     * FFT is serial, as before. */
    int fftw_nthr_saved = Y(fftw_nthreads_hook) ? Y(fftw_nthreads_hook)() : 1;
    if (Y(fftw_plan_with_nthreads_hook))
      Y(fftw_plan_with_nthreads_hook)(pl->nthr);
    pln->pfwd = FFTW(plan_dft)((int)d, pln->narr, (FC *)pln->g1, (FC *)pln->g2,
                         FFTW_FORWARD, ff);
    pln->pback =
         FFTW(plan_dft)((int)d, pln->narr, (FC *)pln->g2, (FC *)pln->g1,
                          FFTW_BACKWARD, ff);
    if (Y(fftw_plan_with_nthreads_hook))
      Y(fftw_plan_with_nthreads_hook)(fftw_nthr_saved);
  }

  /* Children are measured on their own, so the forward's CONV child would
   * gather from a g2 the DECONV step never wrote. Uninitialised bytes read as
   * doubles are frequently denormal, and a denormal operand can put the
   * multiply on a microcode path an order of magnitude slower -- a plausible
   * wrong number rather than a crash. This must run before the children are
   * planned below, not merely before they are applied: FFTW_MEASURE plans
   * (just created above) already execute real transforms over g1/g2. */
  memset(pln->g1, 0, (size_t)ntot * sizeof(C));
  memset(pln->g2, 0, (size_t)ntot * sizeof(C));

  /* Recurse into the planner for the children. */
  {
    /* DECONV child: f_hat (parent) -> g1; sign +1 (forward orientation). */
    pln->deconv_prob = Y(mkproblem_deconv)(d, pln->Nc, pn->variant, pln->nc, m,
                                        window, +1, pn->f_hat, pln->g1);
    pln->deconv_child = Y(planner_mkplan)(pl, pln->deconv_prob);
    /* CONV child: g2 -> f (parent); x borrowed from parent's owned copy. */
    pln->conv_prob = Y(mkproblem_conv)(d, pln->nc, pln->Nc, M, m, window, +1,
                                    pn->x_sorted, pn->x, pln->g2, pn->f);
    ((problem_conv *)pln->conv_prob)->analytic = pn->analytic;
    pln->conv_child = Y(planner_mkplan)(pl, pln->conv_prob);
  }
  /* FFTW_WISDOM_ONLY can leave an FFTW plan NULL; a child solver can decline. */
  if (!pln->pfwd || !pln->pback || !pln->deconv_child || !pln->conv_child) {
    destroy(&pln->super);
    Y(free)(pln);
    return 0;
  }
  /* The estimate is the sum of the three steps it actually got, which is what
   * the TODO this replaces asked for. FFTW reports its own plan's counts, so
   * the FFT term tracks whatever FFTW chose -- SIMD codelets today, a threaded
   * FFT later -- instead of a fitted n log n. */
  {
    double fa, fm, ff;
    Y(ops_add)(&pln->deconv_child->ops, &pln->conv_child->ops, &pln->super.ops);
    FFTW(flops)(pln->pfwd, &fa, &fm, &ff);
    pln->super.ops.add += fa;
    pln->super.ops.mul += fm;
    pln->super.ops.fma += ff;
  }
  pln->super.uses_x = 1; /* reads x via the CONV child */
  return &pln->super;
}
static const solver_adt native_fast_solver_adt = {NFFT_PROBLEM_NFFT, 0,
                                                  mkplan_native_fast};
void Y(nfft_solver_fast_native_register)(planner *pl)
{
  REGISTER_SOLVER(pl,
                  Y(solver_create)(sizeof(solver), &native_fast_solver_adt));
}
