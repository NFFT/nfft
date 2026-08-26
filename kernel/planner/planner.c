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

#include <string.h>

#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"

/* Planner: the solver registry plus the wisdom store. The store is a pair of
 * open-addressed, double-hashed tables (blessed and unblessed) keyed by a
 * 128-bit problem signature, with subsumption-aware lookup and insert. */

/* Must be prime and >= 2: the probe arithmetic divides by size and size - 1,
 * and only a prime size makes the probe sequence cover every slot. */
#define PLNR_INITIAL_TABLE_SIZE 7u

/* The subsumption law. Solution a answers query q, its key already matched,
 * iff
 *   a feasible:   LEQ(a.u, q.u) && LEQ(q.l, a.l)
 *   a infeasible: LEQ(a.l, q.l) && a.timelimit_imp <= q.timelimit_imp */
static int subsumes(const flags_t *a, const flags_t *q) {
  if (a->slvndx != INFEASIBLE_SLVNDX) {
    /* Feasible solutions always have timelimit_imp == 0. */
    A(a->timelimit_imp == 0);
    return LEQ(a->u, q->u) && LEQ(q->l, a->l);
  }
  return LEQ(a->l, q->l) && a->timelimit_imp <= q->timelimit_imp;
}

static int sig_eq(const md5sig a, const md5sig b) {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
}

/* A prime capacity >= 2, every slot invalid. */
static void fresh_table(hashtab *t, unsigned size) {
  A(size >= 2);
  t->size = size;
  t->nelem = 0;
  t->nrehash = 0;
  t->entries = (solution *)Y(malloc)((size_t)size * sizeof(solution));
  /* info == 0 => PLNR_H_VALID clear => slot invalid. */
  memset(t->entries, 0, (size_t)size * sizeof(solution));
}

static void free_table(hashtab *t) {
  Y(free)
  (t->entries);
  t->entries = 0;
  t->size = 0;
  t->nelem = 0;
}

/* Double hashing. A kill marks a slot dead but still valid, so a probe chain
 * ends at the first slot that was never written. */
static unsigned probe_start(const md5sig s, unsigned size) {
  return (unsigned)(s[0] % size);
}

static unsigned probe_step(const md5sig s, unsigned size) {
  return 1u + (unsigned)(s[1] % (size - 1u));
}

/* Place sol at the first non-live slot on its probe chain, no subsumption
 * check. Caller must guarantee room (rehash, post-growth insert). */
static void raw_place(hashtab *t, const solution *sol) {
  unsigned h = probe_start(sol->s, t->size);
  unsigned step = probe_step(sol->s, t->size);
  unsigned k;
  for (k = 0; k < t->size; k++) {
    solution *slot = t->entries + h;
    if (!(slot->flags.info & PLNR_H_LIVE)) {
      *slot = *sol;
      return;
    }
    h += step;
    if (h >= t->size)
      h -= t->size;
  }
  A(0 /* table full: growth invariant violated */);
}

/* Grow to newsize, re-inserting the live entries only. */
static void rehash(hashtab *t, unsigned newsize) {
  solution *old = t->entries;
  unsigned oldsize = t->size;
  unsigned live = t->nelem;
  unsigned i;

  A(Y(is_prime)((INT)newsize) && newsize >= 2);
  t->entries = (solution *)Y(malloc)((size_t)newsize * sizeof(solution));
  memset(t->entries, 0, (size_t)newsize * sizeof(solution));
  t->size = newsize;

  for (i = 0; i < oldsize; i++)
    if (old[i].flags.info & PLNR_H_LIVE)
      raw_place(t, old + i);

  t->nelem = live;
  t->nrehash++;
  Y(free)
  (old);
}

/* Keep size > nelem + nelem/8 (load factor below ~8/9) after the pending
 * insert, so at least one non-live slot remains and probe loops terminate. */
static void maybe_grow(hashtab *t) {
  unsigned need = t->nelem + 1u; /* live count after the pending insert */
  if (t->size > need + need / 8u)
    return;
  {
    /* Target the square of the load bound, (9/8)^2, so the grown table does
     * not re-trip maybe_grow at once. */
    INT target = (INT)((81u * need) / 64u) + 2;
    INT ns = Y(next_prime)(target);
    rehash(t, (unsigned)ns);
  }
}

/* Search one table for a live matching-key entry subsuming q; among those,
 * return one minimal in u under the LEQ order. NULL if none. */
static solution *table_lookup(hashtab *t, const md5sig s, const flags_t *q) {
  unsigned h = probe_start(s, t->size);
  unsigned step = probe_step(s, t->size);
  unsigned k;
  solution *best = 0;
  for (k = 0; k < t->size; k++) {
    solution *slot = t->entries + h;
    if (!(slot->flags.info & PLNR_H_VALID))
      break;
    if ((slot->flags.info & PLNR_H_LIVE) && sig_eq(slot->s, s) && subsumes(&slot->flags, q)) {
      if (best == 0 || LEQ(slot->flags.u, best->flags.u))
        best = slot;
    }
    h += step;
    if (h >= t->size)
      h -= t->size;
  }
  return best;
}

/* The returned pointer aims into a table's entries array and dies at the next
 * insert (which may rehash) or forget. Consume it before mutating the store;
 * a search loop must copy what it needs, not hold the pointer. */
solution *Y(planner_hlookup)(planner *pl, const md5sig s, const flags_t *q) {
  solution *sol = table_lookup(&pl->htab_blessed, s, q);
  if (sol != 0)
    return sol;
  return table_lookup(&pl->htab_unblessed, s, q);
}

void Y(planner_hinsert)(planner *pl, const md5sig s, const flags_t *f,
                        unsigned slvndx) {
  hashtab *t = (f->info & PLNR_BLESSING) ? &pl->htab_blessed
                                         : &pl->htab_unblessed;
  solution newsol;
  unsigned h, step, k;
  int first_killed = -1;

  newsol.s[0] = s[0];
  newsol.s[1] = s[1];
  newsol.s[2] = s[2];
  newsol.s[3] = s[3];
  newsol.flags = *f;
  newsol.flags.slvndx = slvndx;
  newsol.flags.info = (f->info & PLNR_BLESSING) | PLNR_H_VALID | PLNR_H_LIVE;
  A(slvndx == INFEASIBLE_SLVNDX || f->timelimit_imp == 0);

  /* Kill every live matching-key entry the newcomer subsumes. Mutual
   * subsumption (equal bounds) counts, which is what heals stale wisdom. */
  h = probe_start(s, t->size);
  step = probe_step(s, t->size);
  for (k = 0; k < t->size; k++) {
    solution *slot = t->entries + h;
    if (!(slot->flags.info & PLNR_H_VALID))
      break;
    if ((slot->flags.info & PLNR_H_LIVE) && sig_eq(slot->s, s)) {
      if (subsumes(&newsol.flags, &slot->flags)) {
        slot->flags.info &= ~(unsigned)PLNR_H_LIVE; /* dead but valid */
        t->nelem--;
        if (first_killed < 0)
          first_killed = (int)h;
      } else {
        /* The slot already answers every query the newcomer would, so the
         * insert is redundant. Drop it: release builds must keep the
         * one-live-entry-per-key invariant rather than duplicate the key. */
        if (subsumes(&slot->flags, &newsol.flags)) {
          A(0);
          return;
        }
      }
    }
    h += step;
    if (h >= t->size)
      h -= t->size;
  }

  if (first_killed >= 0) {
    t->entries[first_killed] = newsol; /* reuse a just-killed slot */
    t->nelem++;
  } else {
    maybe_grow(t);
    raw_place(t, &newsol);
    t->nelem++;
  }
}

static void reset_table(hashtab *t) {
  free_table(t);
  fresh_table(t, PLNR_INITIAL_TABLE_SIZE);
}

void Y(planner_forget)(planner *pl, amnesia a) {
  reset_table(&pl->htab_unblessed);
  if (a == PLNR_FORGET_ALL)
    reset_table(&pl->htab_blessed);
}

/* Preamble token, e.g. "nfft_wisdom-3.5.4alpha" (double), "nfftf_..." (float).
 * Encodes precision (via Y) and release (PACKAGE_VERSION), so wisdom never
 * crosses either boundary. */
#define WISDOM_PREAMBLE (STRINGIZE(Y(wisdom)) "-" PACKAGE_VERSION)

/* MD5 over sizeof(R), then every registered solver's reg_id and registrar name
 * in registration order. A wisdom file is honoured only by a library whose
 * registry reproduces this signature; a mismatch rejects the import instead of
 * yielding a wrong plan. */
static void config_signature(planner *pl, md5sig out) {
  md5 m;
  unsigned i;

  Y(md5_begin)
  (&m);
  Y(md5_put_unsigned)
  (&m, (unsigned)sizeof(R));
  for (i = 0; i < pl->nslvdesc; i++) {
    slvdesc *d = pl->slvdescs + i;
    Y(md5_put_int)
    (&m, d->reg_id);
    Y(md5_put_str)
    (&m, d->reg_nam ? d->reg_nam : "");
  }
  Y(md5_end)
  (&m);
  out[0] = m.s[0];
  out[1] = m.s[1];
  out[2] = m.s[2];
  out[3] = m.s[3];
}

/* Deterministic string hash pre-filtering registrar-name comparisons on wisdom
 * import. Depends on every character of s; no length cap. Quality is
 * uncritical -- it only avoids a strcmp, never decides correctness. */
static unsigned reg_nam_hash(const char *s) {
  unsigned h = 5381u;
  const unsigned char *p = (const unsigned char *)s;
  while (*p != '\0') {
    h = h * 33u + (unsigned)*p;
    p++;
  }
  return h;
}

/* Resolve an imported (registrar-name, reg-id) pair to a descriptor index.
 * Returns INFEASIBLE_SLVNDX on a miss; a real index is always strictly below
 * it (asserted at registration), so the importer reads that as "unknown
 * solver" unambiguously. */
static unsigned resolve_slvndx(planner *pl, const char *name, int id) {
  unsigned h = reg_nam_hash(name);
  unsigned i;

  for (i = 0; i < pl->nslvdesc; i++) {
    slvdesc *d = pl->slvdescs + i;
    if (d->nam_hash == h && d->reg_nam != 0 && strcmp(d->reg_nam, name) == 0 && d->reg_id == id)
      return i;
  }
  return INFEASIBLE_SLVNDX;
}

void Y(planner_export)(planner *pl, printer *p) {
  hashtab *t = &pl->htab_blessed;
  md5sig cfg;
  unsigned i;

  config_signature(pl, cfg);
  p->print(p, "(%s #x%w #x%w #x%w #x%w", WISDOM_PREAMBLE, cfg[0], cfg[1], cfg[2],
           cfg[3]);

  /* One line per live blessed entry, in table-slot order. */
  for (i = 0; i < t->size; i++) {
    solution *e = t->entries + i;
    const char *name;
    int id;

    if (!(e->flags.info & PLNR_H_LIVE))
      continue;

    if (e->flags.slvndx == INFEASIBLE_SLVNDX) {
      name = "!"; /* reserved: cannot collide with a C-identifier name */
      id = 0;
    } else {
      slvdesc *d = pl->slvdescs + e->flags.slvndx;
      name = d->reg_nam ? d->reg_nam : "";
      id = d->reg_id;
    }

    p->print(p, "\n  (%s %d #x%x #x%x #x%x #x%w #x%w #x%w #x%w)", name, id,
             (unsigned)e->flags.l, (unsigned)e->flags.u,
             (unsigned)e->flags.timelimit_imp, e->s[0], e->s[1], e->s[2],
             e->s[3]);
  }
  p->print(p, "\n)");
}

int Y(planner_import)(planner *pl, scanner *sc) {
  hashtab *t = &pl->htab_blessed;
  hashtab snap;
  char token[64];
  md5uint w0, w1, w2, w3;
  md5sig cfg;

  /* Reject before touching the store: wisdom is a cache, never authoritative,
   * so a mismatch must degrade to a miss. */
  if (!sc->scan(sc, "(%*s #x%w #x%w #x%w #x%w", (int)sizeof(token) - 1, token,
                &w0, &w1, &w2, &w3))
    return 0;
  if (strcmp(token, WISDOM_PREAMBLE) != 0)
    return 0;
  config_signature(pl, cfg);
  if (cfg[0] != w0 || cfg[1] != w1 || cfg[2] != w2 || cfg[3] != w3)
    return 0;

  /* Deep copy so a malformed entry rolls the whole import back: import is
   * atomic. */
  snap = *t;
  snap.entries = (solution *)Y(malloc)((size_t)t->size * sizeof(solution));
  memcpy(snap.entries, t->entries, (size_t)t->size * sizeof(solution));

  for (;;) {
    char name[64];
    int id;
    unsigned l, u, tli;
    md5uint s0, s1, s2, s3;
    md5sig s;
    unsigned slvndx;
    flags_t f, q;

    /* Trying ')' then '(' is safe: the scanner has single-char pushback. */
    if (sc->scan(sc, ")"))
      break;
    if (!sc->scan(sc, "("))
      goto malformed;

    if (!sc->scan(sc, "%*s %d #x%x #x%x #x%x #x%w #x%w #x%w #x%w)",
                  (int)sizeof(name) - 1, name, &id, &l, &u, &tli, &s0, &s1, &s2,
                  &s3))
      goto malformed;

    if (strcmp(name, "!") == 0) {
      slvndx = INFEASIBLE_SLVNDX; /* the only legitimate infeasible name */
    } else {
      slvndx = resolve_slvndx(pl, name, id);
      if (slvndx == INFEASIBLE_SLVNDX)
        goto malformed; /* unknown solver */
      if (tli != 0)
        goto malformed; /* feasible entries always have timelimit_imp == 0 */
    }

    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
    s[3] = s3;
    f.l = (unsigned)(l & 0xFFFFu);
    f.u = (unsigned)(u & 0xFFFFu);
    f.timelimit_imp = (unsigned)(tli & 0xFFFu);
    f.info = PLNR_BLESSING;
    f.slvndx = slvndx;

    /* Every entry must satisfy LEQ(l, u). A crafted file could otherwise
     * plant a live blessed entry that strictly subsumes a later legitimate
     * blessing insert. */
    if (!LEQ(f.l, f.u))
      goto malformed;

    /* Idempotent re-import. Blessed-only lookup, so an unblessed session
     * entry never suppresses a persisted blessed insert. */
    q = f;
    if (table_lookup(&pl->htab_blessed, s, &q) == 0)
      Y(planner_hinsert)
    (pl, s, &f, slvndx);
  }

  Y(free)
  (snap.entries); /* commit */
  return 1;

malformed:
  /* Roll back to the snapshot. */
  Y(free)
  (t->entries);
  *t = snap;
  return 0;
}

planner *Y(planner_create)(void) {
  planner *pl = (planner *)Y(malloc)(sizeof(planner));
  int i;

  pl->slvdescs = 0;
  pl->nslvdesc = 0;
  pl->slvdescsiz = 0;
  pl->cur_reg_nam = 0;
  pl->cur_reg_id = 0;
  for (i = 0; i < NFFT_PROBLEM_LAST; i++)
    pl->kind_head[i] = -1;

  fresh_table(&pl->htab_blessed, PLNR_INITIAL_TABLE_SIZE);
  fresh_table(&pl->htab_unblessed, PLNR_INITIAL_TABLE_SIZE);

  pl->nthr = 1;
  memset(&pl->flags, 0, sizeof(pl->flags));
  pl->timelimit_seconds = -1.0; /* unlimited default */
  return pl;
}

void Y(planner_register_solver)(planner *pl, solver *s) {
  slvdesc *d;
  int kind;

  if (s == 0)
    return;

  Y(solver_use)
  (s);

  if (pl->nslvdesc == pl->slvdescsiz) {
    unsigned newsiz = pl->slvdescsiz ? 2u * pl->slvdescsiz : 8u;
    slvdesc *nd = (slvdesc *)Y(malloc)((size_t)newsiz * sizeof(slvdesc));
    if (pl->slvdescs != 0) {
      memcpy(nd, pl->slvdescs, (size_t)pl->nslvdesc * sizeof(slvdesc));
      Y(free)
      (pl->slvdescs);
    }
    pl->slvdescs = nd;
    pl->slvdescsiz = newsiz;
  }

  d = pl->slvdescs + pl->nslvdesc;
  d->slv = s;
  d->reg_nam = pl->cur_reg_nam;
  d->nam_hash = pl->cur_reg_nam ? reg_nam_hash(pl->cur_reg_nam) : 0u;
  d->reg_id = pl->cur_reg_id++;

  kind = s->adt->problem_kind;
  d->next_same_kind = pl->kind_head[kind];
  pl->kind_head[kind] = (int)pl->nslvdesc;

  pl->nslvdesc++;
  A(pl->nslvdesc < INFEASIBLE_SLVNDX);
}

void Y(planner_destroy)(planner *pl) {
  unsigned i;

  for (i = 0; i < pl->nslvdesc; i++)
    Y(solver_destroy)
  (pl->slvdescs[i].slv);
  Y(free)
  (pl->slvdescs);

  free_table(&pl->htab_blessed);
  free_table(&pl->htab_unblessed);

  Y(free)
  (pl);
}

/* -1.0 is the "unlimited" sentinel. */
double Y(planner_timelimit)(const planner *pl) {
  return pl->timelimit_seconds;
}

void Y(planner_set_timelimit)(planner *pl, double seconds) {
  pl->timelimit_seconds = seconds < 0.0 ? -1.0 : seconds;
}

/* The wisdom key: one md5 context fed sizeof(R) and the thread count, then the
 * problem's own hash. This is the only place sizeof(R) and nthr enter a key,
 * so a concrete problem hash must not repeat them. */
void Y(problem_md5)(planner *pl, const problem *p, md5sig out) {
  md5 m;

  Y(md5_begin)
  (&m);
  Y(md5_put_unsigned)
  (&m, (unsigned)sizeof(R));
  Y(md5_put_int)
  (&m, pl->nthr);
  p->adt->hash(p, &m);
  Y(md5_end)
  (&m);
  out[0] = m.s[0];
  out[1] = m.s[1];
  out[2] = m.s[2];
  out[3] = m.s[3];
}

/* Query flags from the planner's current impatience bounds; timelimit_imp and
 * info stay zero (unblessed session memo). */
static flags_t search_flags(planner *pl) {
  flags_t q;
  q.l = PLNR_L(pl);
  q.u = PLNR_U(pl);
  q.timelimit_imp = 0;
  q.info = 0;
  q.slvndx = 0;
  return q;
}

/* Consult wisdom, else try every registered solver of the problem's kind and
 * keep the cheapest by pcost, memoising the outcome unblessed. Estimate-only:
 * pcost is the solver's analytic number, set at mkplan time. */
plan *Y(planner_mkplan)(planner *pl, const problem *p) {
  md5sig sig;
  flags_t q;
  solution *sol;
  plan *best = 0;
  slvdesc *best_desc = 0;

  Y(problem_md5)
  (pl, p, sig);
  q = search_flags(pl);

  /* A memoised infeasibility short-circuits to NULL; a feasible hit reruns
   * that solver's mkplan. A stale hit (mkplan returns NULL) degrades to a
   * miss and falls through to the full search, never to a failure. */
  sol = Y(planner_hlookup)(pl, sig, &q);
  if (sol != 0) {
    unsigned slvndx = sol->flags.slvndx; /* copy before any store mutation */
    if (slvndx == INFEASIBLE_SLVNDX)
      return 0;
    {
      solver *s = pl->slvdescs[slvndx].slv;
      plan *pln = s->adt->mkplan(s, p, pl);
      if (pln != 0)
        return pln;
    }
  }

  /* The strict compare is the determinism contract: ties keep the
   * earlier-encountered candidate. */
  FORALL_SOLVERS_OF_KIND(p->adt->kind, pl, s, d, {
    plan *pln = s->adt->mkplan(s, p, pl);
    if (pln != 0) {
      if (best == 0 || pln->pcost < best->pcost) {
        if (best != 0)
          Y(plan_destroy)
        (best);
        best = pln;
        best_desc = d;
      } else {
        Y(plan_destroy)
        (pln);
      }
    }
  });

  /* Memoise unblessed. Reached only when the lookup did not answer, so the
   * insert never duplicates a live hit. */
  {
    unsigned slvndx = best_desc != 0 ? (unsigned)(best_desc - pl->slvdescs)
                                     : INFEASIBLE_SLVNDX;
    Y(planner_hinsert)
    (pl, sig, &q, slvndx);
  }

  return best;
}

/* Every applicable solver of p's kind, in the same reverse-registration order
 * as planner_mkplan, storing each non-NULL plan and its descriptor index. No
 * wisdom lookup and no store mutation: the caller owns every returned plan. */
int Y(planner_candidates)(planner *pl, const problem *p, plan **plans,
                          unsigned *slvndx, int cap) {
  int count = 0;

  FORALL_SOLVERS_OF_KIND(p->adt->kind, pl, s, d, {
    plan *pln = s->adt->mkplan(s, p, pl);
    if (pln != 0) {
      if (count >= cap) {
        A(0 /* more candidates than cap: caller must pass a larger array */);
        Y(plan_destroy)
        (pln);
      } else {
        plans[count] = pln;
        slvndx[count] = (unsigned)(d - pl->slvdescs);
        count++;
      }
    }
  });

  return count;
}

/* Bless the winner of a measured race. An equal-bounds entry with a different
 * winner is replaced by the mutual-subsumption kill in planner_hinsert. */
void Y(planner_bless)(planner *pl, const problem *p, unsigned slvndx) {
  md5sig sig;
  flags_t f;
  solution *existing;

  A(slvndx != INFEASIBLE_SLVNDX);
  A(slvndx < pl->nslvdesc); /* wild index must not be persisted */
  A(!(PLNR_U(pl) & PLNR_ESTIMATE)); /* estimate-grade evidence never blessed */

  Y(problem_md5)
  (pl, p, sig);

  f.l = PLNR_L(pl);
  f.u = PLNR_U(pl);
  f.timelimit_imp = 0;
  f.info = PLNR_BLESSING;
  f.slvndx = slvndx;

  existing = table_lookup(&pl->htab_blessed, sig, &f);
  if (existing != 0 && existing->flags.slvndx == slvndx)
    return;

  Y(planner_hinsert)
  (pl, sig, &f, slvndx);
}
