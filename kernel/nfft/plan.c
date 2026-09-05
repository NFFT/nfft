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
#include "plan_ng_test.h"

/* Direction indices into the bundle's per-direction arrays. */
#define FWD 0 /* sign = +1 */
#define ADJ 1 /* sign = -1 */

struct Y(plan_ng_s) {
  problem *prob[2]; /* [FWD] built and raced; [ADJ] reserved for a separate
                     * adjoint plan, unpopulated (NULL) */
  plan *dir[2];     /* [FWD] the winning plan; [ADJ] NULL */
};

/* Strip the FFTW preservation bits (FFTW_DESTROY_INPUT / FFTW_PRESERVE_INPUT)
 * out of fftw_flags before it reaches the problem's wisdom key: no planner-
 * native candidate mutates its input in place, so the two spellings must not
 * key distinct entries. The remaining bits do affect measured cost and belong
 * in the key. */
static unsigned keyable_fftw_flags(unsigned fftw_flags) {
  return fftw_flags & ~(unsigned)(FFTW_DESTROY_INPUT | FFTW_PRESERVE_INPUT);
}

/* Map the public NFFT_* planning gates to their PLNR_* images. */
static unsigned map_planning_flags(unsigned planning) {
  unsigned F = 0;
  if (planning & NFFT_NO_DIRECT)
    F |= PLNR_NO_DIRECT;
  if (planning & NFFT_NO_FAST_NATIVE)
    F |= PLNR_NO_FAST_NATIVE;
  return F;
}

/* Estimate-mode selection under the planner's current bounds (the caller sets
 * l/u). Fills p->dir[]; returns 1 if at least one direction succeeded. */
static int select_estimate(Y(plan_ng) * p, planner *pl) {
  p->dir[FWD] = p->prob[FWD] ? Y(planner_mkplan)(pl, p->prob[FWD]) : 0;
  p->dir[ADJ] = p->prob[ADJ] ? Y(planner_mkplan)(pl, p->prob[ADJ]) : 0;
  return (p->dir[FWD] != 0 || p->dir[ADJ] != 0);
}

Y(plan_ng) * Y(plan_ng_guru)(int d, const INT *N, const int *variant,
                             const INT *n, INT M, int m, int window, R *x,
                             FC *f_hat, FC *f,
                             unsigned fftw_flags, unsigned planning) {
  Y(plan_ng) * p;
  planner *pl;
  unsigned saved_l, saved_u;
  unsigned F; /* PLNR_* image of the public gate flags */
  int is_estimate;

  Y(nfft_ensure_registered)
  ();
  pl = Y(the_planner)();

  /* Refresh thread count before any keying. */
  pl->nthr = (int)Y(get_num_threads)();

  F = map_planning_flags(planning);
  is_estimate = (planning & NFFT_ESTIMATE) ? 1 : 0;

  /* x/f_hat/f are unconditionally required, in both estimate and measured
   * mode: a plan is always bound to real caller-owned pointers, mirroring
   * FFTW's guru contract. */
  if (d <= 0 || N == 0 || n == 0 || x == 0 || f_hat == 0 || f == 0)
    return 0;
  if (M < (INT)1 || m < 1)
    return 0;
  if (window < NFFT_WINDOW_KAISER_BESSEL || window > NFFT_WINDOW_DIRAC_DELTA)
    return 0;
  {
    int t;
    for (t = 0; t < d; t++)
      if (N[t] <= 0 || n[t] <= 0)
        return 0;
  }

  /* Save the planner's bounds; will be restored on every path. */
  saved_l = pl->flags.l;
  saved_u = pl->flags.u;

  p = (Y(plan_ng) *)Y(malloc)(sizeof(*p));
  p->dir[FWD] = 0;
  p->dir[ADJ] = 0;
  p->prob[FWD] = Y(mkproblem_nfft)(d, N, variant, n, M, m, window, +1,
                                   keyable_fftw_flags(fftw_flags), x,
                                   /*copy_x=*/1, (C *)f_hat, (C *)f);
  p->prob[ADJ] = 0;

  /* With the fast solver in play the geometry and the window must satisfy its
   * applicability predicate. Refused here rather than left to the search,
   * which would serve an O(N*M) direct plan without telling the caller.
   * NFFT_NO_FAST_NATIVE lifts both; unit axes are elided at construction, so
   * the predicate sees only the surviving axes. */
  if (!(planning & NFFT_NO_FAST_NATIVE) &&
      (window > NFFT_WINDOW_SINC_POWER ||
       !Y(nfft_fast_guards_ok)(p->prob[FWD], m))) {
    Y(problem_destroy)
    (p->prob[FWD]);
    Y(free)
    (p);
    return 0;
  }

  if (is_estimate) {
    /* Bounds {l = F, u = PLNR_ESTIMATE | F}; memos stay unblessed. */
    pl->flags.l = F;
    pl->flags.u = PLNR_ESTIMATE | F;

    if (!select_estimate(p, pl)) {
      pl->flags.l = saved_l;
      pl->flags.u = saved_u;
      if (p->prob[FWD])
        Y(problem_destroy)
      (p->prob[FWD]);
      if (p->prob[ADJ])
        Y(problem_destroy)
      (p->prob[ADJ]);
      Y(free)
      (p);
      return 0;
    }

    pl->flags.l = saved_l;
    pl->flags.u = saved_u;

    return p;
  }

  /* Measured mode. Bounds {l = u = F} for queries and bless inserts. */
  pl->flags.l = F;
  pl->flags.u = F;

  /* Candidate arrays (cap 8 per direction). ncands[dirn] counts the plans
   * still owned by cands[dirn]; it drops back to 0 once they are destroyed or
   * a winner has moved into dir[dirn], so no cleanup path can free a plan
   * twice. */
  plan *cands[2][8];
  unsigned cslvndx[2][8];
  int ncands[2];
  int dirn;

  for (dirn = FWD; dirn <= ADJ; dirn++) {
    md5sig sig;
    flags_t q;
    solution *sol;

    ncands[dirn] = 0;

    if (!p->prob[dirn])
      continue;

    Y(problem_md5)
    (pl, p->prob[dirn], sig);
    q.l = F;
    q.u = F;
    q.timelimit_imp = 0;
    q.info = 0;
    q.slvndx = 0;

    sol = Y(planner_hlookup)(pl, sig, &q);
    if (sol != 0 && sol->flags.slvndx == INFEASIBLE_SLVNDX)
      continue; /* infeasible hit: the direction is absent */
    if (sol != 0) {
      /* Feasible hit: copy slvndx before any store mutation, re-run
       * solver's mkplan to adopt without a race. */
      unsigned hit_slvndx = sol->flags.slvndx;
      solver *s;
      plan *pln;
      A(hit_slvndx < pl->nslvdesc);
      s = pl->slvdescs[hit_slvndx].slv;
      pln = s->adt->mkplan(s, p->prob[dirn], pl);
      if (pln != 0) {
        p->dir[dirn] = pln;
        continue;
      }
      /* Stale hit: fall through to the race (bless will replace it). */
    }

    /* Enumerate candidates under current measured bounds. */
    ncands[dirn] = Y(planner_candidates)(pl, p->prob[dirn], cands[dirn],
                                         cslvndx[dirn], 8);
  }

  /* A direction is absent when it has neither a plan nor a candidate. */
  if (p->dir[FWD] == 0 && ncands[FWD] == 0 && p->dir[ADJ] == 0 &&
      ncands[ADJ] == 0) {
    pl->flags.l = saved_l;
    pl->flags.u = saved_u;
    if (p->prob[FWD])
      Y(problem_destroy)
    (p->prob[FWD]);
    if (p->prob[ADJ])
      Y(problem_destroy)
    (p->prob[ADJ]);
    Y(free)
    (p);
    return 0;
  }

  /* Measured race: a candidate is awoken SLEEPY -> AWAKE_ZERO, measured, and
   * put back to SLEEPY. The running best stays AWAKE_ZERO, so the adopted
   * winner's psi is built once and Y(precompute) only upgrades it to AWAKE. */
  {
    double tl = Y(planner_timelimit)(pl);
    double t_start = (tl >= 0.0) ? Y(planner_clock_now)() : 0.0;

    int need_restart = 0;

#ifdef NFFT_DEBUG
    /* x-restore guard: snapshot the top-level problem's owned x copy before
     * any candidate is awoken, re-verified whenever one goes back to SLEEPY. */
    md5sig x_sig_entry;
    /* x is the COMPRESSED copy (rnk*M reals, not d*M). Sizing the guard from
     * the original d would over-read out of bounds in debug builds. */
    INT dM_guard = (INT)((const problem_nfft *)p->prob[FWD])->sz->rnk * M;
    Y(nfft_x_md5)
    (((const problem_nfft *)p->prob[FWD])->x, dM_guard, x_sig_entry);
#define VERIFY_X                                                              \
  A(Y(nfft_x_verify)(((const problem_nfft *)p->prob[FWD])->x, dM_guard,       \
                     x_sig_entry))
#else
#define VERIFY_X ((void)0)
#endif

    for (dirn = FWD; dirn <= ADJ; dirn++) {
      int ci;
      int winner_idx;
      double winner_cost;
      double prune_bound; /* PLNR_PRUNE_RATIO * the cheapest pcost */
      int nc;
      int timed_out = 0;

      if (ncands[dirn] == 0)
        continue; /* hit-adopted or absent: nothing to race */

      nc = ncands[dirn]; /* snapshot: ncands drops to 0 on adoption */

      /* Estimate gate (PLNR_PRUNE_RATIO, iplanner.h). Survivors race on
       * measured cost. */
      {
        double min_pcost = cands[dirn][0]->pcost;
        int nsurvive = 0, lone = 0;
        for (ci = 1; ci < nc; ci++)
          if (cands[dirn][ci]->pcost < min_pcost)
            min_pcost = cands[dirn][ci]->pcost;
        prune_bound = PLNR_PRUNE_RATIO * min_pcost;
        for (ci = 0; ci < nc; ci++)
          if (cands[dirn][ci]->pcost <= prune_bound) {
            nsurvive++;
            lone = ci;
          }
        if (nsurvive == 1) {
          /* Lone survivor: adopt untimed and bless immediately. */
          for (ci = 0; ci < nc; ci++)
            if (ci != lone)
              Y(plan_destroy)(cands[dirn][ci]);
          p->dir[dirn] = cands[dirn][lone];
          ncands[dirn] = 0;
          Y(planner_bless)(pl, p->prob[dirn], cslvndx[dirn][lone]);
          continue;
        }
      }

      /* The candidates run on the caller's arrays, which may hold anything --
       * f_hat is filled after planning -- and NaNs or denormals there would
       * skew every timing. Zero both once, before the first measurement. */
      {
        problem_nfft *pz = (problem_nfft *)p->prob[dirn];
        memset(pz->f_hat, 0,
               (size_t)Y(problem_nfft_Ntot)(p->prob[dirn]) * sizeof(C));
        memset(pz->f, 0, (size_t)pz->M * sizeof(C));
      }

      /* Two or more survivors: measure each. Strict-smallest cost wins; ties
       * keep the earlier-encountered candidate. */
      winner_idx = -1;
      winner_cost = -1.0;

      for (ci = 0; ci < nc; ci++) {
        double cost;
        if (cands[dirn][ci]->pcost > prune_bound)
          continue;
        if (tl >= 0.0 && Y(planner_elapsed_seconds)(t_start) >= tl) {
          timed_out = 1;
          break;
        }

        Y(plan_awake)
        (cands[dirn][ci], PLNR_AWAKE_ZERO);
        cost = Y(plan_measure_cost)(cands[dirn][ci], p->prob[dirn]);

        if (cost < 0.0) {
          need_restart = 1;
          break;
        }

        cands[dirn][ci]->pcost = cost;

        /* Exactly one candidate is left awake: the running best. */
        if (winner_idx < 0 || cost < winner_cost) {
          if (winner_idx >= 0) {
            Y(plan_awake)
            (cands[dirn][winner_idx], PLNR_SLEEPY);
            VERIFY_X;
          }
          winner_idx = ci;
          winner_cost = cost;
        } else {
          Y(plan_awake)
          (cands[dirn][ci], PLNR_SLEEPY);
          VERIFY_X;
        }
      }

      if (timed_out && winner_idx < 0)
        need_restart = 1;

      if (need_restart)
        break;

      /* The winner is never returned to SLEEPY, so verify x here instead. */
      VERIFY_X;

      /* Destroy losers, adopt the winner. */
      for (ci = 0; ci < nc; ci++) {
        if (ci != winner_idx)
          Y(plan_destroy)
        (cands[dirn][ci]);
      }
      p->dir[dirn] = cands[dirn][winner_idx];
      ncands[dirn] = 0;

      if (timed_out) {
        /* Partial race: an untimed survivor may be faster, so this winner is
         * memoised unblessed instead of exported as complete evidence. */
        md5sig sig;
        flags_t fl;
        Y(problem_md5)
        (pl, p->prob[dirn], sig);
        fl.l = PLNR_L(pl);
        fl.u = PLNR_U(pl);
        fl.timelimit_imp = 0;
        fl.info = 0;
        fl.slvndx = cslvndx[dirn][winner_idx];
        Y(planner_hinsert)
        (pl, sig, &fl, cslvndx[dirn][winner_idx]);
      } else
        Y(planner_bless)
        (pl, p->prob[dirn], cslvndx[dirn][winner_idx]);
    }

    /* Shared estimate-grade restart. */
    if (need_restart) {
      int dj, cj;
      for (dj = FWD; dj <= ADJ; dj++) {
        int nc_dj = ncands[dj];
        ncands[dj] = 0;
        for (cj = 0; cj < nc_dj; cj++)
          Y(plan_destroy)
        (cands[dj][cj]);
      }
      if (p->dir[FWD]) {
        Y(plan_destroy)
        (p->dir[FWD]);
        p->dir[FWD] = 0;
      }
      if (p->dir[ADJ]) {
        Y(plan_destroy)
        (p->dir[ADJ]);
        p->dir[ADJ] = 0;
      }

      pl->flags.l = F;
      pl->flags.u = PLNR_ESTIMATE | F;

      if (!select_estimate(p, pl)) {
        pl->flags.l = saved_l;
        pl->flags.u = saved_u;
        if (p->prob[FWD])
          Y(problem_destroy)
        (p->prob[FWD]);
        if (p->prob[ADJ])
          Y(problem_destroy)
        (p->prob[ADJ]);
        Y(free)
        (p);
        return 0;
      }

      pl->flags.l = saved_l;
      pl->flags.u = saved_u;

      return p;
    }
  }

  pl->flags.l = saved_l;
  pl->flags.u = saved_u;

  return p;

#undef VERIFY_X
}

/* Awaken every present direction plan to PLNR_AWAKE, which is what makes the
 * measured winner's placeholder tables hold true values. Mandatory lifecycle
 * step: execute asserts AWAKE. */
void Y(precompute)(Y(plan_ng) * p) {
  int i;
  A(p != 0);

  for (i = 0; i < 2; ++i) {
    if (p->dir[i])
      Y(plan_awake)
    (p->dir[i], PLNR_AWAKE);
  }
}

void Y(execute)(Y(plan_ng) * p) {
  A(p != 0);
  A(p->dir[FWD] != 0);
  A(p->dir[FWD]->awake_state == PLNR_AWAKE);
  p->dir[FWD]->adt->apply(p->dir[FWD], p->prob[FWD]);
}

void Y(execute_adjoint)(Y(plan_ng) * p) {
  A(p != 0);
  A(p->dir[FWD] != 0);
  A(p->dir[FWD]->awake_state == PLNR_AWAKE);
  p->dir[FWD]->adt->apply_adjoint(p->dir[FWD], p->prob[FWD]);
}

/* New-array execute. Plans read their data from the problem, so the problem
 * pointers are swapped and restored around the apply; x and the psi tables are
 * unaffected, they depend on x alone. */
void Y(execute_on)(Y(plan_ng) * p, FC *f_hat, FC *f) {
  problem_nfft *pf;
  C *save_pf_hat, *save_pf;
  A(p != 0);
  A(p->dir[FWD] != 0);
  A(p->dir[FWD]->awake_state == PLNR_AWAKE);
  A(f_hat != 0);
  A(f != 0);

  pf = (problem_nfft *)p->prob[FWD];
  save_pf_hat = pf->f_hat;
  save_pf = pf->f;
  pf->f_hat = (C *)f_hat;
  pf->f = (C *)f;
  p->dir[FWD]->adt->apply(p->dir[FWD], p->prob[FWD]);
  pf->f_hat = save_pf_hat;
  pf->f = save_pf;
}

void Y(execute_adjoint_on)(Y(plan_ng) * p, FC *f_hat,
                           FC *f) {
  problem_nfft *pf;
  C *save_pf_hat, *save_pf;
  A(p != 0);
  A(p->dir[FWD] != 0);
  A(p->dir[FWD]->awake_state == PLNR_AWAKE);
  A(f_hat != 0);
  A(f != 0);

  pf = (problem_nfft *)p->prob[FWD];
  save_pf_hat = pf->f_hat;
  save_pf = pf->f;
  pf->f_hat = (C *)f_hat;
  pf->f = (C *)f;
  p->dir[FWD]->adt->apply_adjoint(p->dir[FWD], p->prob[FWD]);
  pf->f_hat = save_pf_hat;
  pf->f = save_pf;
}

/* Format: (nfft-plan-ng (fwd %p) (adj %p)), each %p the direction plan's own
 * self-description or (null) when the direction is absent. */
void Y(plan_ng_print)(Y(plan_ng) * p, printer *pr) {
  A(p != 0);
  pr->print(pr, "(nfft-plan-ng (fwd %p) (adj %p))", p->dir[FWD], p->dir[ADJ]);
}

/* Test-only: the winning forward plan's wakefulness (prototype in
 * include/plan_ng_test.h). */
int Y(plan_ng_test_awake_state)(const Y(plan_ng) * p) {
  A(p != 0);
  A(p->dir[FWD] != 0);
  return p->dir[FWD]->awake_state;
}

void Y(plan_ng_destroy)(Y(plan_ng) * p) {
  if (!p)
    return;

  /* Y(plan_destroy) asserts a non-NULL plan, so guard each direction. It
   * sleeps the plan first, which is what releases an awake winner's tables. */
  if (p->dir[FWD])
    Y(plan_destroy)
  (p->dir[FWD]);
  if (p->dir[ADJ])
    Y(plan_destroy)
  (p->dir[ADJ]);

  if (p->prob[FWD])
    Y(problem_destroy)
  (p->prob[FWD]);
  if (p->prob[ADJ])
    Y(problem_destroy)
  (p->prob[ADJ]);

  Y(free)
  (p);
}
