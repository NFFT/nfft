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
#include "plan_ng_test.h" /* prototype for the test-only awake_state hook */

/* NFFT plan.
 *
 * One plan ownes two problems/plans ([FWD] sign=+1, [ADJ] sign=-1) of which only [FWD] 
 * is populated and raced today: the adjoint reuses dir[FWD] via apply_adjoint. [ADJ] 
 * is dormant for a potential option with two plans. */
#undef X
#define X(name) NFFT(name)

/* Direction indices into the bundle's per-direction arrays. */
#define FWD 0 /* sign = +1 */
#define ADJ 1 /* sign = -1 */

struct Y(plan_ng_s) {
  problem *prob[2]; /* [FWD] sign=+1 built+raced; [ADJ] dormant (future 2-plan) */
  plan *dir[2];     /* [FWD] the winning plan; [ADJ] dormant (NULL).            */
};

/* Wisdom-key flag normalisation */

/* Strip the FFTW preservation bits (FFTW_DESTROY_INPUT / FFTW_PRESERVE_INPUT)
 * out of fftw_flags before it reaches the problem's wisdom key. Those bits
 * are irrelevant to every planner-native candidate (none of them mutate their
 * input in place), so two problems differing only in a caller's preservation-
 * bit spelling would otherwise get distinct wisdom-key entries for plans that
 * behave identically. Others (FFTW_ESTIMATE/MEASURE/PATIENT/...) flow through 
 * unchanged — those do affect measured cost and belong in the key. */
static unsigned keyable_fftw_flags(unsigned fftw_flags) {
  return fftw_flags & ~(unsigned)(FFTW_DESTROY_INPUT | FFTW_PRESERVE_INPUT);
}

/* Flag mapping: public NFFT_* -> internal PLNR_* (F) */

/* Map the public NFFT_* planning gates to their PLNR_* images. */
static unsigned map_planning_flags(unsigned planning) {
  unsigned F = 0;
  if (planning & NFFT_NO_DIRECT)
    F |= PLNR_NO_DIRECT;
  if (planning & NFFT_NO_NDFT_PLAIN)
    F |= PLNR_NO_NDFT_PLAIN;
  if (planning & NFFT_NO_NDFT_BLOCKED)
    F |= PLNR_NO_NDFT_BLOCKED;
  if (planning & NFFT_NO_FAST_NATIVE)
    F |= PLNR_NO_FAST_NATIVE;
  return F;
}

/* Estimate-mode selection */
/* Run the estimate-mode flow (two planner_mkplan calls) under the planner's
 * current bounds (caller must have set l/u appropriately).  Fills p->dir[].
 * Returns 1 if at least one direction succeeded, 0 if both failed. */
static int select_estimate(Y(plan_ng) * p, planner *pl) {
  p->dir[FWD] = p->prob[FWD] ? Y(planner_mkplan)(pl, p->prob[FWD]) : 0;
  p->dir[ADJ] = p->prob[ADJ] ? Y(planner_mkplan)(pl, p->prob[ADJ]) : 0;
  return (p->dir[FWD] != 0 || p->dir[ADJ] != 0);
}

/* Guru constructor */
Y(plan_ng) * Y(plan_ng_guru)(int d, const INT *N, const int *variant,
                             const INT *n, INT M, int m, int window, R *x,
                             FC *f_hat, FC *f,
                             unsigned fftw_flags, unsigned planning) {
  Y(plan_ng) * p;
  planner *pl;
  unsigned saved_l, saved_u;
  unsigned F; /* PLNR_* image of the public gate flags */
  int is_estimate;

  /* Ensure the roster is registered into the process-global planner. */
  Y(nfft_ensure_registered)
  ();
  pl = Y(the_planner)();

  /* Refresh thread count before any keying. */
  pl->nthr = (int)Y(get_num_threads)();

  /* Map public flags to PLNR_* image F and determine mode. */
  F = map_planning_flags(planning);
  is_estimate = (planning & NFFT_ESTIMATE) ? 1 : 0;

  /* x/f_hat/f are unconditionally required, in both estimate and measured
   * mode: a plan is always bound to real caller-owned pointers, mirroring
   * FFTW's guru contract. */
  if (d <= 0 || N == 0 || n == 0 || x == 0 || f_hat == 0 || f == 0)
    return 0;

  /* Per-axis geometry must be positive; with the fast solver in play, every
   * axis surviving unit-axis elision also needs sigma = n/N > 1. Refused here
   * rather than left to the search, which would serve an O(N*M) direct plan
   * without telling the caller. NFFT_NO_FAST_NATIVE makes sigma <= 1 legal
   * again; unit axes are elided and stay legal with n[t] == 1. */
  {
    int t;
    int wants_fast = !(planning & NFFT_NO_FAST_NATIVE);
    for (t = 0; t < d; t++) {
      if (N[t] <= 0 || n[t] <= 0)
        return 0;
      if (wants_fast && N[t] > (INT)1 && !(n[t] > N[t]))
        return 0;
    }
  }

  /* Save the planner's bounds; will be restored on every path. */
  saved_l = pl->flags.l;
  saved_u = pl->flags.u;

  /* Allocate the bundle. */
  p = (Y(plan_ng) *)Y(malloc)(sizeof(*p));
  p->dir[FWD] = 0;
  p->dir[ADJ] = 0;

  /* Build the forward problem only. */
  {
    unsigned key_fftw_flags = keyable_fftw_flags(fftw_flags);
    p->prob[FWD] = Y(mkproblem_nfft)(d, N, variant, n, M, m, window, +1,
                                     key_fftw_flags, x, /*copy_x=*/1, (C *)f_hat,
                                     (C *)f);
    p->prob[ADJ] = 0;
  }

  /* ------------------------------------------------------------------ */
  if (is_estimate) {
    /* ESTIMATE mode.  Bounds
     *   {l = F, u = PLNR_ESTIMATE | F}
     * Two planner_mkplan calls; memos stay unblessed. */
    pl->flags.l = F;
    pl->flags.u = PLNR_ESTIMATE | F;

    if (!select_estimate(p, pl)) {
      /* Both directions absent (or one excluded, one infeasible). */
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

  /* MEASURED mode. */

  /* Bounds {l = u = F} for measured queries and bless inserts. */
  pl->flags.l = F;
  pl->flags.u = F;

  /* Candidate arrays (cap 8 per direction).
   *
   * ncands[dirn] sentinels:
   *   > 0  : candidates enumerated, awaiting race
   *   -1   : infeasible-absent (infeasible wisdom hit)
   *   -2   : hit-adopted (feasible wisdom hit; plan is in dir[dirn])
   *   -3   : no-candidates-absent (planner_candidates returned 0)
   *   -4   : raced/lone-adopted (winner chosen; plan is in dir[dirn])
   *   -5   : not-requested-absent (under the forward-only race,
   *          prob[ADJ] is never built; a future two-internal-plan option
   *          would populate it and race it here)
   *
   * The -4 sentinel is written immediately after adopting a winner or lone
   * candidate so that the no-clock restart cleanup loop never re-walks that
   * direction's candidate array and double-frees the plan already in dir[]. */
  plan *cands[2][8];
  unsigned cslvndx[2][8];
  int ncands[2];
  int dirn;

  /* First pass: enumerate all candidates for both directions. */
  for (dirn = FWD; dirn <= ADJ; dirn++) {
    /* Wisdom lookup for this direction. */
    md5sig sig;
    flags_t q;
    solution *sol;

    if (!p->prob[dirn]) {
      /* Not requested (forward-only race): prob[ADJ] is never built,
       * so there is nothing to look up, race, or bless. */
      ncands[dirn] = -5; /* sentinel: not-requested-absent */
      continue;
    }

    Y(problem_md5)
    (pl, p->prob[dirn], sig);
    q.l = F;
    q.u = F;
    q.timelimit_imp = 0;
    q.info = 0;
    q.slvndx = 0;

    sol = Y(planner_hlookup)(pl, sig, &q);
    if (sol != 0 && sol->flags.slvndx == INFEASIBLE_SLVNDX) {
      /* Infeasible hit: direction is absent; no candidates. */
      ncands[dirn] = -1; /* sentinel: infeasible-absent */
      continue;
    }
    if (sol != 0) {
      /* Feasible hit: copy slvndx before any store mutation, re-run
       * solver's mkplan to adopt without a race. */
      unsigned hit_slvndx = sol->flags.slvndx; /* copy before mutation */
      solver *s;
      plan *pln;
      A(hit_slvndx < pl->nslvdesc);
      s = pl->slvdescs[hit_slvndx].slv;
      pln = s->adt->mkplan(s, p->prob[dirn], pl);
      if (pln != 0) {
        /* Adopt the hit plan: no race needed, skip to next direction. */
        p->dir[dirn] = pln;
        ncands[dirn] = -2; /* sentinel: hit-adopted */
        continue;
      }
      /* Stale hit: fall through to the race (bless will replace it). */
    }

    /* Enumerate candidates under current measured bounds. */
    ncands[dirn] = Y(planner_candidates)(pl, p->prob[dirn], cands[dirn],
                                         cslvndx[dirn], 8);

    if (ncands[dirn] == 0) {
      /* No candidates: direction absent (user-forced bits). */
      ncands[dirn] = -3; /* sentinel: no-candidates-absent */
      continue;
    }
  }

  /* Absent check: a direction is absent iff dir[dirn]==NULL and ncands is one
   * of the absent sentinels (-1 infeasible, -3 no-candidates).
   * ncands == -2 (hit-adopted) and -4 (raced-adopted) have dir[dirn] set. */
#define DIR_IS_ABSENT(d_) \
  (p->dir[d_] == 0 && (ncands[d_] == -1 || ncands[d_] == -3 || ncands[d_] == -5))

  if (DIR_IS_ABSENT(FWD) && DIR_IS_ABSENT(ADJ)) {
    /* Destroy any live candidate arrays (ncands > 0 = un-raced candidates).
     * All negative sentinels have no live candidate arrays. */
    for (dirn = FWD; dirn <= ADJ; dirn++) {
      int ci;
      if (ncands[dirn] > 0)
        for (ci = 0; ci < ncands[dirn]; ci++)
          Y(plan_destroy)
      (cands[dirn][ci]);
    }
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

  /* Measured race. For each candidate the  */
  /* race awakens it SLEEPY -> AWAKE_ZERO, measures, then AWAKE_ZERO -> SLEEPY */
  {
    /* Capture the global timelimit budget at race entry. */
    double tl = Y(planner_timelimit)(pl);
    double t_start = (tl >= 0.0) ? Y(planner_clock_now)() : 0.0;

    int timed_out = 0;
    int need_restart = 0;

#ifdef NFFT_DEBUG
    /* x-restore byte-identity guard: snapshot the top-level
     * problem's owned x copy once, before any candidate is awoken.  Re-
     * verified after each candidate's AWAKE_ZERO -> SLEEPY restore below.
     * Today's candidates are read-only over x, so this is a tautology --
     * intentional substrate for a future permuting solver. */
    md5sig x_sig_entry;
    /* x is the COMPRESSED copy (rnk*M reals, not d*M). Sizing the guard from
     * the original d would over-read out of bounds in debug builds. */
    INT dM_guard = (INT)((const problem_nfft *)p->prob[FWD])->sz->rnk * M;
    Y(nfft_x_md5)
    (((const problem_nfft *)p->prob[FWD])->x, dM_guard, x_sig_entry);
#endif

    for (dirn = FWD; dirn <= ADJ; dirn++) {
      int ci;
      int winner_idx;
      double winner_cost;
      double prune_bound; /* PLNR_PRUNE_RATIO * the cheapest pcost */
      int nc; /* snapshot of ncands[dirn] at loop entry */

      if (ncands[dirn] <= 0)
        continue; /* hit-adopted or absent — nothing to race */

      nc = ncands[dirn]; /* snapshot: ncands may be overwritten below */

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
          ncands[dirn] = -4; /* sentinel: raced/lone-adopted */
          Y(planner_bless)(pl, p->prob[dirn], cslvndx[dirn][lone]);
          continue;
        }
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

        /* Awaken the candidate to AWAKE_ZERO. */
        Y(plan_awake)
        (cands[dirn][ci], PLNR_AWAKE_ZERO);
        cost = Y(plan_measure_cost)(cands[dirn][ci], p->prob[dirn]);
        Y(plan_awake)
        (cands[dirn][ci], PLNR_SLEEPY);

#ifdef NFFT_DEBUG
        A(Y(nfft_x_verify)(((const problem_nfft *)p->prob[FWD])->x, dM_guard,
                           x_sig_entry)); /* byte-identical */
#endif

        if (cost < 0.0) {
          need_restart = 1;
          break;
        }

        cands[dirn][ci]->pcost = cost;

        if (winner_idx < 0 || cost < winner_cost) {
          winner_idx = ci;
          winner_cost = cost;
        }
      }

      if (timed_out && winner_idx < 0)
        need_restart = 1;

      if (need_restart)
        break;

      /* Winner found: destroy losers, adopt winner (left SLEEPY), bless. */
      for (ci = 0; ci < nc; ci++) {
        if (ci != winner_idx)
          Y(plan_destroy)
        (cands[dirn][ci]);
      }
      p->dir[dirn] = cands[dirn][winner_idx];
      ncands[dirn] = -4; /* sentinel: raced/lone-adopted */
      Y(planner_bless)
      (pl, p->prob[dirn], cslvndx[dirn][winner_idx]);
    }

    /* Shared estimate-grade restart. */
    if (need_restart) {
      int dj, cj;
      for (dj = FWD; dj <= ADJ; dj++) {
        if (ncands[dj] > 0) {
          int nc_dj = ncands[dj];
          ncands[dj] = 0;
          for (cj = 0; cj < nc_dj; cj++)
            Y(plan_destroy)
          (cands[dj][cj]);
        }
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

#undef DIR_IS_ABSENT
}

/* Precompute: awaken every present direction plan (idempotent via state) */

/* The awakening: awakens every present direction plan via
 * Y(plan_awake) to PLNR_AWAKE.  For a fast native winner this triggers its
 * awake hook, which builds the psi tables once per awake period (guarded by
 * the plan's own `precomputed` flag).  A winner left SLEEPY by the race
 * rebuilds psi here (the accepted one-extra build).  For an NDFT native
 * winner plan_awake records AWAKE; its awake hook is NULL so no precompute
 * runs. This is a mandatory lifecycle step, execute asserts AWAKE. */
void Y(precompute)(Y(plan_ng) * p) {
  int i;
  A(p != 0);

  for (i = 0; i < 2; ++i) {
    if (p->dir[i])
      Y(plan_awake)
    (p->dir[i], PLNR_AWAKE);
  }
}

/* Execute: apply the planner-chosen plan for the requested direction     */

void Y(execute)(Y(plan_ng) * p) {
  A(p != 0);
  A(p->dir[FWD] != 0); /* the single forward-winning plan */
  A(p->dir[FWD]->awake_state == PLNR_AWAKE);
  p->dir[FWD]->adt->apply(p->dir[FWD], p->prob[FWD]);
}

void Y(execute_adjoint)(Y(plan_ng) * p) {
  A(p != 0);
  /* Forward-only race: the single plan serves adjoint too. */
  A(p->dir[FWD] != 0);
  A(p->dir[FWD]->awake_state == PLNR_AWAKE);
  p->dir[FWD]->adt->apply_adjoint(p->dir[FWD], p->prob[FWD]);
}

/* New-array execute: run on caller-supplied f_hat/f, x fixed at plan time */
/* Plans read data from the problem; swap the problem pointers so the
 * winner runs on the caller-supplied arrays, then restore (mirrors
 * fftw_execute_dft leaving the plan's arrays intact). x and the psi tables
 * are unaffected -- they depend on x alone. */
void Y(execute_on)(Y(plan_ng) * p, FC *f_hat, FC *f) {
  problem_nfft *pf;
  C *save_pf_hat, *save_pf;
  A(p != 0);
  A(p->dir[FWD] != 0); /* the single plan serves both directions */
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
  A(p->dir[FWD] != 0); /* the single plan serves both directions */
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

/* Format: (nfft-plan-ng (fwd %p) (adj %p))
 * where each %p prints the direction plan's own self-description, or (null)
 * when a direction is absent. */
void Y(plan_ng_print)(Y(plan_ng) * p, printer *pr) {
  A(p != 0);
  /* Forward-only race: one internal plan (dir[FWD]); the adjoint
   * reuses it, so (adj ...) prints (null) today. */
  pr->print(pr, "(nfft-plan-ng (fwd %p) (adj %p))", p->dir[FWD], p->dir[ADJ]);
}

/* Test-only: the winning forward plan's wakefulness (prototype in
 * include/plan_ng_test.h). Lets acceptance tests assert a fresh guru's winner
 * is returned SLEEPY and that precompute awakens it to AWAKE. */
int Y(plan_ng_test_awake_state)(const Y(plan_ng) * p) {
  A(p != 0);
  A(p->dir[FWD] != 0);
  return p->dir[FWD]->awake_state;
}

void Y(plan_ng_destroy)(Y(plan_ng) * p) {
  if (!p)
    return;

  /* Y(plan_destroy) is not documented NULL-safe (it asserts a non-NULL plan),
   * so guard each direction ourselves. Destroying a plan first returns it to
   * SLEEPY  then runs its destroy hook. */
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
