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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "nfft3.h" /* Y(malloc)/Y(free) declarations */
#include "infft.h"
#include "iplanner.h"
#include "planner.h"

/* Convert an RFC-1321 digest (16 bytes, as printed by md5sum) into the four
 * little-endian words of md5.s (word convention). */
static void digest_to_words(const unsigned char d[16], md5uint w[4]) {
  int i;
  for (i = 0; i < 4; i++)
    w[i] = (md5uint)d[4 * i] | ((md5uint)d[4 * i + 1] << 8) | ((md5uint)d[4 * i + 2] << 16) | ((md5uint)d[4 * i + 3] << 24);
}

static int md5_matches(const char *msg, const unsigned char digest[16]) {
  md5 m;
  md5uint w[4];
  Y(md5_begin)
  (&m);
  Y(md5_put_bytes)
  (&m, msg, strlen(msg));
  Y(md5_end)
  (&m);
  digest_to_words(digest, w);
  return m.s[0] == w[0] && m.s[1] == w[1] && m.s[2] == w[2] && m.s[3] == w[3];
}

void Y(check_planner_md5_rfc_vectors)(void) {
  /* RFC 1321, appendix A.5 test suite */
  static const unsigned char d_empty[16] = {0xd4, 0x1d, 0x8c, 0xd9, 0x8f, 0x00,
                                            0xb2, 0x04, 0xe9, 0x80, 0x09, 0x98, 0xec, 0xf8, 0x42, 0x7e};
  static const unsigned char d_a[16] = {0x0c, 0xc1, 0x75, 0xb9, 0xc0, 0xf1,
                                        0xb6, 0xa8, 0x31, 0xc3, 0x99, 0xe2, 0x69, 0x77, 0x26, 0x61};
  static const unsigned char d_abc[16] = {0x90, 0x01, 0x50, 0x98, 0x3c, 0xd2,
                                          0x4f, 0xb0, 0xd6, 0x96, 0x3f, 0x7d, 0x28, 0xe1, 0x7f, 0x72};
  static const unsigned char d_md[16] = {0xf9, 0x6b, 0x69, 0x7d, 0x7c, 0xb7,
                                         0x93, 0x8d, 0x52, 0x5a, 0x2f, 0x31, 0xaa, 0xf1, 0x61, 0xd0};
  static const unsigned char d_az[16] = {0xc3, 0xfc, 0xd3, 0xd7, 0x61, 0x92,
                                         0xe4, 0x00, 0x7d, 0xfb, 0x49, 0x6c, 0xca, 0x67, 0xe1, 0x3b};
  static const unsigned char d_alnum[16] = {0xd1, 0x74, 0xab, 0x98, 0xd2, 0x77,
                                            0xd9, 0xf5, 0xa5, 0x61, 0x1c, 0x2c, 0x9f, 0x41, 0x9d, 0x9f};
  static const unsigned char d_num[16] = {0x57, 0xed, 0xf4, 0xa2, 0x2b, 0xe3,
                                          0xc9, 0x55, 0xac, 0x49, 0xda, 0x2e, 0x21, 0x07, 0xb6, 0x7a};

  CU_ASSERT(md5_matches("", d_empty));
  CU_ASSERT(md5_matches("a", d_a));
  CU_ASSERT(md5_matches("abc", d_abc));
  CU_ASSERT(md5_matches("message digest", d_md));
  CU_ASSERT(md5_matches("abcdefghijklmnopqrstuvwxyz", d_az));
  /* 62 bytes: block boundary crossed only by md5_end padding */
  CU_ASSERT(md5_matches(
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
      d_alnum));
  /* 80 bytes: exercises mid-stream block compression */
  CU_ASSERT(md5_matches("1234567890123456789012345678901234567890"
                        "1234567890123456789012345678901234567890",
                        d_num));
}

void Y(check_planner_md5_feeders)(void) {
  md5 m1, m2;

  /* md5_put_str also hashes the trailing '\0' */
  Y(md5_begin)
  (&m1);
  Y(md5_put_str)
  (&m1, "abc");
  Y(md5_end)
  (&m1);
  Y(md5_begin)
  (&m2);
  Y(md5_put_bytes)
  (&m2, "abc", 4); /* include '\0' */
  Y(md5_end)
  (&m2);
  CU_ASSERT(m1.s[0] == m2.s[0] && m1.s[1] == m2.s[1] && m1.s[2] == m2.s[2] && m1.s[3] == m2.s[3]);

  /* determinism + sensitivity of the typed feeders */
  Y(md5_begin)
  (&m1);
  Y(md5_put_int)
  (&m1, 42);
  Y(md5_put_INT)
  (&m1, (INT)7);
  Y(md5_put_unsigned)
  (&m1, 3u);
  Y(md5_end)
  (&m1);
  Y(md5_begin)
  (&m2);
  Y(md5_put_int)
  (&m2, 42);
  Y(md5_put_INT)
  (&m2, (INT)7);
  Y(md5_put_unsigned)
  (&m2, 3u);
  Y(md5_end)
  (&m2);
  CU_ASSERT(m1.s[0] == m2.s[0] && m1.s[3] == m2.s[3]);

  Y(md5_begin)
  (&m2);
  Y(md5_put_int)
  (&m2, 43); /* single-bit input difference */
  Y(md5_put_INT)
  (&m2, (INT)7);
  Y(md5_put_unsigned)
  (&m2, 3u);
  Y(md5_end)
  (&m2);
  CU_ASSERT(m1.s[0] != m2.s[0] || m1.s[1] != m2.s[1] || m1.s[2] != m2.s[2] || m1.s[3] != m2.s[3]);
}

void Y(check_planner_printer)(void) {
  char buf[256];
  size_t cnt;
  printer *p;

  /* string backend + the directives wisdom export uses */
  p = Y(printer_create_str)(buf);
  p->print(p, "(%s %d #x%x #x%w)", "solver-a", -42, 0xbeefu, (md5uint)0x1234u);
  Y(printer_destroy)
  (p);
  CU_ASSERT_STRING_EQUAL(buf, "(solver-a -42 #xbeef #x00001234)");

  /* %D (INT), %u, %c, literal chars */
  p = Y(printer_create_str)(buf);
  p->print(p, "%D/%u%c", (INT)123456789, 77u, '!');
  Y(printer_destroy)
  (p);
  CU_ASSERT_STRING_EQUAL(buf, "123456789/77!");

  /* NULL string */
  p = Y(printer_create_str)(buf);
  p->print(p, "%s", (char *)0);
  Y(printer_destroy)
  (p);
  CU_ASSERT_STRING_EQUAL(buf, "(null)");

  /* %x of zero prints one digit */
  p = Y(printer_create_str)(buf);
  p->print(p, "%x", 0u);
  Y(printer_destroy)
  (p);
  CU_ASSERT_STRING_EQUAL(buf, "0");

  /* indentation: %( newline+indent, %) dedent */
  p = Y(printer_create_str)(buf);
  p->print(p, "a%(b%)");
  Y(printer_destroy)
  (p);
  CU_ASSERT_STRING_EQUAL(buf, "a\n  b");

  /* counter backend agrees with string backend */
  p = Y(printer_create_cnt)(&cnt);
  p->print(p, "(%s %d #x%x #x%w)", "solver-a", -42, 0xbeefu, (md5uint)0x1234u);
  Y(printer_destroy)
  (p);
  CU_ASSERT_EQUAL(cnt, strlen("(solver-a -42 #xbeef #x00001234)"));

  /* file backend: flushed on destroy, content identical */
  {
    FILE *f = tmpfile();
    long n;
    CU_ASSERT_FATAL(f != NULL);
    p = Y(printer_create_file)(f);
    p->print(p, "hello %d%(world%)", 7);
    Y(printer_destroy)
    (p); /* must flush */
    n = ftell(f);
    rewind(f);
    CU_ASSERT_EQUAL(fread(buf, 1, (size_t)n, f), (size_t)n);
    buf[n] = 0;
    CU_ASSERT_STRING_EQUAL(buf, "hello 7\n  world");
    fclose(f);
  }
}

void Y(check_planner_scanner)(void) {
  scanner *sc;
  char name[64];
  int d;
  unsigned x;
  md5uint m;

  /* parse a wisdom-entry-shaped line */
  sc = Y(scanner_create_str)("(solver-a -42 #xbeef #x00001234)");
  CU_ASSERT(sc->scan(sc, "(%*s %d #x%x #x%w)", (int)sizeof(name) - 1, name,
                     &d, &x, &m));
  CU_ASSERT_STRING_EQUAL(name, "solver-a");
  CU_ASSERT_EQUAL(d, -42);
  CU_ASSERT_EQUAL(x, 0xbeefu);
  CU_ASSERT(m == (md5uint)0x1234u);
  Y(scanner_destroy)
  (sc);

  /* whitespace-tolerant around ( ) and tokens */
  sc = Y(scanner_create_str)("  (\n  hello   7 )");
  CU_ASSERT(sc->scan(sc, "(%*s %d)", (int)sizeof(name) - 1, name, &d));
  CU_ASSERT_STRING_EQUAL(name, "hello");
  CU_ASSERT_EQUAL(d, 7);
  Y(scanner_destroy)
  (sc);

  /* literal mismatch -> 0 */
  sc = Y(scanner_create_str)("(foo)");
  CU_ASSERT_FALSE(sc->scan(sc, "(bar"));
  Y(scanner_destroy)
  (sc);

  /* %d is strict: no digits -> 0, and letters are not digits */
  sc = Y(scanner_create_str)("()");
  CU_ASSERT_FALSE(sc->scan(sc, "(%d", &d));
  Y(scanner_destroy)
  (sc);
  sc = Y(scanner_create_str)("(x)");
  CU_ASSERT_FALSE(sc->scan(sc, "(%d", &d));
  Y(scanner_destroy)
  (sc);

  /* %*s truncates safely */
  sc = Y(scanner_create_str)("abcdefgh");
  CU_ASSERT(sc->scan(sc, "%*s", 4, name));
  CU_ASSERT_STRING_EQUAL(name, "abcd");
  Y(scanner_destroy)
  (sc);

  /* %d rejects out-of-range integers -> mismatch (no int narrowing, no
   * -(INT_MIN) UB) */
  sc = Y(scanner_create_str)("(4294967296)");
  CU_ASSERT_FALSE(sc->scan(sc, "(%d", &d));
  Y(scanner_destroy)
  (sc);
  sc = Y(scanner_create_str)("(-4294967296)");
  CU_ASSERT_FALSE(sc->scan(sc, "(%d", &d));
  Y(scanner_destroy)
  (sc);
  /* INT_MIN is in range and parses (boundary; negate in the long domain) */
  sc = Y(scanner_create_str)("(-2147483648)");
  CU_ASSERT(sc->scan(sc, "(%d)", &d));
  CU_ASSERT_EQUAL(d, (-2147483647 - 1));
  Y(scanner_destroy)
  (sc);

  /* print -> scan round trip through a file */
  {
    FILE *f = tmpfile();
    printer *p;
    CU_ASSERT_FATAL(f != NULL);
    p = Y(printer_create_file)(f);
    p->print(p, "(entry 3 #x%w)\n", (md5uint)0xdeadbeefu);
    Y(printer_destroy)
    (p);
    rewind(f);
    sc = Y(scanner_create_file)(f);
    CU_ASSERT(sc->scan(sc, "(%*s %d #x%w)", (int)sizeof(name) - 1, name, &d, &m));
    CU_ASSERT_STRING_EQUAL(name, "entry");
    CU_ASSERT_EQUAL(d, 3);
    CU_ASSERT(m == (md5uint)0xdeadbeefu);
    Y(scanner_destroy)
    (sc);
    fclose(f);
  }
}

/* Mock solvers */

static const solver_adt mock_adt_nfft = {NFFT_PROBLEM_NFFT, 0};
static const solver_adt mock_adt_deconv = {NFFT_PROBLEM_DECONV, 0};

static void reg_mock_nfft_a(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &mock_adt_nfft));
}

static void reg_mock_nfft_b(planner *pl) {
  /* registers two solvers under one registrar name -> reg_id 0 and 1 */
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &mock_adt_nfft));
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &mock_adt_nfft));
}

static void reg_mock_deconv(planner *pl) {
  REGISTER_SOLVER(pl, Y(solver_create)(sizeof(solver), &mock_adt_deconv));
}

static struct solvtab_s mock_solvtab[] = {SOLVTAB(reg_mock_nfft_a),
                                          SOLVTAB(reg_mock_nfft_b), SOLVTAB(reg_mock_deconv), SOLVTAB_END};

static planner *mk_test_planner(void) {
  planner *pl = Y(planner_create)();
  Y(solvtab_exec)
  (mock_solvtab, pl);
  return pl;
}

/* Helpers */

static void mksig(int i, md5sig s) {
  md5 m;
  Y(md5_begin)
  (&m);
  Y(md5_put_int)
  (&m, i);
  Y(md5_end)
  (&m);
  s[0] = m.s[0];
  s[1] = m.s[1];
  s[2] = m.s[2];
  s[3] = m.s[3];
}

static flags_t mkflags(unsigned l, unsigned u, unsigned info) {
  flags_t f;
  f.l = l;
  f.u = u;
  f.timelimit_imp = 0;
  f.info = info;
  f.slvndx = 0;
  return f;
}

void Y(check_planner_registry)(void) {
  planner *pl = mk_test_planner();
  int n;

  CU_ASSERT_EQUAL(pl->nslvdesc, 4u); /* a=1, b=2, deconv=1 */

  CU_ASSERT_STRING_EQUAL(pl->slvdescs[0].reg_nam, "reg_mock_nfft_a");
  CU_ASSERT_EQUAL(pl->slvdescs[0].reg_id, 0);
  CU_ASSERT_STRING_EQUAL(pl->slvdescs[1].reg_nam, "reg_mock_nfft_b");
  CU_ASSERT_EQUAL(pl->slvdescs[1].reg_id, 0);
  CU_ASSERT_EQUAL(pl->slvdescs[2].reg_id, 1); /* second from same registrar */

  n = 0;
  FORALL_SOLVERS_OF_KIND(NFFT_PROBLEM_NFFT, pl, s, d, {
    UNUSED(s);
    UNUSED(d);
    ++n;
  });
  CU_ASSERT_EQUAL(n, 3);
  n = 0;
  FORALL_SOLVERS_OF_KIND(NFFT_PROBLEM_DECONV, pl, s, d, {
    UNUSED(s);
    UNUSED(d);
    ++n;
  });
  CU_ASSERT_EQUAL(n, 1);
  n = 0;
  FORALL_SOLVERS_OF_KIND(NFFT_PROBLEM_CONV, pl, s, d, {
    UNUSED(s);
    UNUSED(d);
    ++n;
  });
  CU_ASSERT_EQUAL(n, 0);

  Y(planner_destroy)
  (pl); /* exercises refcounted solver teardown */
}

void Y(check_planner_hashtable)(void) {
  planner *pl = Y(planner_create)();
  md5sig sig;
  flags_t fl, q;
  solution *sol;
  int i;

  /* insert + exact lookup */
  mksig(1, sig);
  fl = mkflags(0, 0, PLNR_BLESSING);
  Y(planner_hinsert)
  (pl, sig, &fl, 1u);
  q = mkflags(0, 0, 0);
  sol = Y(planner_hlookup)(pl, sig, &q);
  CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
  CU_ASSERT_EQUAL(sol->flags.slvndx, 1u);

  /* different sig -> miss */
  mksig(2, sig);
  CU_ASSERT_PTR_NULL(Y(planner_hlookup)(pl, sig, &q));

  /* unblessed insert routes to the other table, still found by lookup */
  fl = mkflags(0, 0, 0);
  Y(planner_hinsert)
  (pl, sig, &fl, 2u);
  sol = Y(planner_hlookup)(pl, sig, &q);
  CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
  CU_ASSERT_EQUAL(sol->flags.slvndx, 2u);
  CU_ASSERT_EQUAL(pl->htab_blessed.nelem, 1u);
  CU_ASSERT_EQUAL(pl->htab_unblessed.nelem, 1u);

  /* growth/rehash: 200 blessed entries, all retrievable afterwards */
  for (i = 100; i < 300; i++) {
    mksig(i, sig);
    fl = mkflags(0, 0, PLNR_BLESSING);
    Y(planner_hinsert)
    (pl, sig, &fl, (unsigned)(i % 50));
  }
  CU_ASSERT_EQUAL(pl->htab_blessed.nelem, 201u); /* 1 + 200 */
  CU_ASSERT(pl->htab_blessed.nrehash > 0);
  CU_ASSERT(Y(is_prime)((INT)pl->htab_blessed.size));
  for (i = 100; i < 300; i++) {
    mksig(i, sig);
    sol = Y(planner_hlookup)(pl, sig, &q);
    CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
    CU_ASSERT_EQUAL(sol->flags.slvndx, (unsigned)(i % 50));
  }

  Y(planner_destroy)
  (pl);
}

void Y(check_planner_subsumption)(void) {
  planner *pl = Y(planner_create)();
  md5sig sig;
  flags_t fl, q;
  solution *sol;

  /* Law: feasible a answers q iff LEQ(a.u,q.u) && LEQ(q.l,a.l).
   * All flags here keep the solution invariant LEQ(l, u). */

  /* patient result (u=0) answers an impatient query (u=PLNR_ESTIMATE) */
  mksig(10, sig);
  fl = mkflags(0, 0, PLNR_BLESSING);
  Y(planner_hinsert)
  (pl, sig, &fl, 3u);
  q = mkflags(0, PLNR_ESTIMATE, 0);
  CU_ASSERT_PTR_NOT_NULL(Y(planner_hlookup)(pl, sig, &q));

  /* matching lower bounds: entry avoiding slow solvers answers a query
   * demanding exactly that */
  mksig(15, sig);
  fl = mkflags(PLNR_NO_SLOW, PLNR_NO_SLOW, PLNR_BLESSING);
  Y(planner_hinsert)
  (pl, sig, &fl, 3u);
  q = mkflags(PLNR_NO_SLOW, PLNR_NO_SLOW | PLNR_ESTIMATE, 0);
  CU_ASSERT_PTR_NOT_NULL(Y(planner_hlookup)(pl, sig, &q));

  /* estimate-mode result cannot answer a patient query */
  mksig(11, sig);
  fl = mkflags(PLNR_ESTIMATE, PLNR_ESTIMATE, PLNR_BLESSING);
  Y(planner_hinsert)
  (pl, sig, &fl, 3u);
  q = mkflags(0, 0, 0);
  CU_ASSERT_PTR_NULL(Y(planner_hlookup)(pl, sig, &q));

  /* entry that may contain "slow" components cannot answer a query whose
   * lower bound forbids them */
  mksig(12, sig);
  fl = mkflags(0, 0, PLNR_BLESSING);
  Y(planner_hinsert)
  (pl, sig, &fl, 3u);
  q = mkflags(PLNR_NO_SLOW, PLNR_NO_SLOW, 0);
  CU_ASSERT_PTR_NULL(Y(planner_hlookup)(pl, sig, &q));

  /* inserting a subsuming entry kills the subsumed one (slot reuse) */
  mksig(13, sig);
  fl = mkflags(0, PLNR_ESTIMATE, PLNR_BLESSING);
  Y(planner_hinsert)
  (pl, sig, &fl, 1u);
  {
    unsigned nelem_before = pl->htab_blessed.nelem;
    fl = mkflags(0, 0, PLNR_BLESSING); /* subsumes the previous entry */
    Y(planner_hinsert)
    (pl, sig, &fl, 2u);
    CU_ASSERT_EQUAL(pl->htab_blessed.nelem, nelem_before); /* replaced */
  }
  q = mkflags(0, PLNR_ESTIMATE, 0);
  sol = Y(planner_hlookup)(pl, sig, &q);
  CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
  CU_ASSERT_EQUAL(sol->flags.slvndx, 2u);
  /* Inserting a subsumed entry after a subsuming live one violates the caller
   * contract: --enable-debug aborts via A(0), a release build drops the
   * redundant insert rather than duplicating the key. */
#ifdef NFFT_DEBUG
  CU_PASS("reverse-order subsumed insert aborts via A(0) under --enable-debug");
#else
  mksig(16, sig);
  fl = mkflags(0, 0, PLNR_BLESSING); /* patient/broad entry: u = 0 */
  Y(planner_hinsert)
  (pl, sig, &fl, 5u);
  {
    unsigned before = pl->htab_blessed.nelem;
    /* u = PLNR_ESTIMATE is subsumed by the broad entry above */
    fl = mkflags(0, PLNR_ESTIMATE, PLNR_BLESSING);
    Y(planner_hinsert)
    (pl, sig, &fl, 6u);
    CU_ASSERT_EQUAL(pl->htab_blessed.nelem, before); /* redundant insert dropped */
  }
  q = mkflags(0, 0, 0);
  sol = Y(planner_hlookup)(pl, sig, &q);
  CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
  CU_ASSERT_EQUAL(sol->flags.slvndx, 5u); /* broad entry retained; no duplicate */
#endif

  /* infeasible entries: "failed under fewer restrictions" answers any
   * more-restricted query */
  mksig(14, sig);
  fl = mkflags(0, 0, PLNR_BLESSING);
  Y(planner_hinsert)
  (pl, sig, &fl, INFEASIBLE_SLVNDX);
  q = mkflags(PLNR_NO_SLOW, PLNR_NO_SLOW, 0);
  sol = Y(planner_hlookup)(pl, sig, &q);
  CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
  CU_ASSERT_EQUAL(sol->flags.slvndx, INFEASIBLE_SLVNDX);

  Y(planner_destroy)
  (pl);
}

void Y(check_planner_forget)(void) {
  planner *pl = Y(planner_create)();
  md5sig sig_b, sig_u;
  flags_t fl, q;

  mksig(20, sig_b);
  fl = mkflags(0, 0, PLNR_BLESSING);
  Y(planner_hinsert)
  (pl, sig_b, &fl, 1u);
  mksig(21, sig_u);
  fl = mkflags(0, 0, 0); /* unblessed */
  Y(planner_hinsert)
  (pl, sig_u, &fl, 2u);

  q = mkflags(0, 0, 0);
  Y(planner_forget)
  (pl, PLNR_FORGET_UNBLESSED);
  CU_ASSERT_PTR_NOT_NULL(Y(planner_hlookup)(pl, sig_b, &q)); /* survives */
  CU_ASSERT_PTR_NULL(Y(planner_hlookup)(pl, sig_u, &q));     /* dropped */

  Y(planner_forget)
  (pl, PLNR_FORGET_ALL);
  CU_ASSERT_PTR_NULL(Y(planner_hlookup)(pl, sig_b, &q));

  Y(planner_destroy)
  (pl);
}

/* export current blessed wisdom into a fresh malloc'd string */
static char *export_to_string(planner *pl) {
  size_t cnt;
  char *s;
  printer *p = Y(printer_create_cnt)(&cnt);
  Y(planner_export)
  (pl, p);
  Y(printer_destroy)
  (p);
  s = (char *)Y(malloc)(cnt + 1);
  p = Y(printer_create_str)(s);
  Y(planner_export)
  (pl, p);
  Y(printer_destroy)
  (p);
  return s;
}

static int import_from_string(planner *pl, const char *s) {
  scanner *sc = Y(scanner_create_str)(s);
  int ret = Y(planner_import)(pl, sc);
  Y(scanner_destroy)
  (sc);
  return ret;
}

void Y(check_planner_wisdom_roundtrip)(void) {
  planner *p1 = mk_test_planner();
  planner *p2 = mk_test_planner();
  md5sig s1, s2, s3;
  flags_t fl, q;
  solution *sol;
  char *wis;

  /* empty round trip */
  wis = export_to_string(p1);
  CU_ASSERT(import_from_string(p2, wis));
  Y(free)
  (wis);

  /* three blessed entries: two feasible (slvndx 0 and 2), one infeasible */
  mksig(30, s1);
  fl = mkflags(0, 0, PLNR_BLESSING);
  Y(planner_hinsert)
  (p1, s1, &fl, 0u);
  mksig(31, s2);
  fl = mkflags(PLNR_NO_SLOW, PLNR_NO_SLOW, PLNR_BLESSING);
  Y(planner_hinsert)
  (p1, s2, &fl, 2u);
  mksig(32, s3);
  fl = mkflags(0, 0, PLNR_BLESSING);
  Y(planner_hinsert)
  (p1, s3, &fl, INFEASIBLE_SLVNDX);
  /* plus one unblessed entry, which must not be exported */
  {
    md5sig s4;
    mksig(33, s4);
    fl = mkflags(0, 0, 0);
    Y(planner_hinsert)
    (p1, s4, &fl, 1u);
  }

  wis = export_to_string(p1);
  CU_ASSERT(strstr(wis, "(! 0 ") != NULL); /* infeasible sentinel */

  CU_ASSERT(import_from_string(p2, wis));
  q = mkflags(0, 0, 0);
  sol = Y(planner_hlookup)(p2, s1, &q);
  CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
  CU_ASSERT_EQUAL(sol->flags.slvndx, 0u);
  q = mkflags(PLNR_NO_SLOW, PLNR_NO_SLOW, 0);
  sol = Y(planner_hlookup)(p2, s2, &q);
  CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
  CU_ASSERT_EQUAL(sol->flags.slvndx, 2u);
  q = mkflags(0, 0, 0);
  sol = Y(planner_hlookup)(p2, s3, &q);
  CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
  CU_ASSERT_EQUAL(sol->flags.slvndx, INFEASIBLE_SLVNDX);
  {
    md5sig s4;
    mksig(33, s4);
    CU_ASSERT_PTR_NULL(Y(planner_hlookup)(p2, s4, &q)); /* unblessed absent */
  }

  /* importing the same wisdom twice is idempotent */
  {
    unsigned nelem = p2->htab_blessed.nelem;
    CU_ASSERT(import_from_string(p2, wis));
    CU_ASSERT_EQUAL(p2->htab_blessed.nelem, nelem);
  }

  /* an unblessed session entry that subsumes an incoming line must not
   * suppress the blessed insert: lookup-first is blessed-only */
  {
    planner *p5 = mk_test_planner();
    flags_t uf = mkflags(0, 0, 0); /* unblessed, same key as s1 */
    Y(planner_hinsert)
    (p5, s1, &uf, 3u);
    CU_ASSERT(import_from_string(p5, wis));
    Y(planner_forget)
    (p5, PLNR_FORGET_UNBLESSED);
    q = mkflags(0, 0, 0);
    sol = Y(planner_hlookup)(p5, s1, &q);
    CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
    CU_ASSERT_EQUAL(sol->flags.slvndx, 0u); /* persisted via blessed table */
    Y(planner_destroy)
    (p5);
  }

  /* file round trip */
  {
    FILE *f = tmpfile();
    planner *p3 = mk_test_planner();
    printer *pr;
    scanner *sc;
    CU_ASSERT_FATAL(f != NULL);
    pr = Y(printer_create_file)(f);
    Y(planner_export)
    (p1, pr);
    Y(printer_destroy)
    (pr);
    rewind(f);
    sc = Y(scanner_create_file)(f);
    CU_ASSERT(Y(planner_import)(p3, sc));
    Y(scanner_destroy)
    (sc);
    fclose(f);
    CU_ASSERT_PTR_NOT_NULL(Y(planner_hlookup)(p3, s1, &q));
    Y(planner_destroy)
    (p3);
  }

  Y(free)
  (wis);
  Y(planner_destroy)
  (p1);
  Y(planner_destroy)
  (p2);
}

void Y(check_planner_wisdom_rejects)(void) {
  planner *p1 = mk_test_planner();
  planner *p2 = mk_test_planner();
  md5sig s1, spre;
  flags_t fl, q;
  char *wis;

  /* one blessed entry pointing at slvndx 1 (= reg_mock_nfft_b, id 0) so the
   * exported wisdom names that registrar -- test (b) relies on it */
  mksig(40, s1);
  fl = mkflags(0, 0, PLNR_BLESSING);
  Y(planner_hinsert)
  (p1, s1, &fl, 1u);
  wis = export_to_string(p1);

  /* (a) truncated input: import fails and pre-existing wisdom survives */
  mksig(41, spre);
  fl = mkflags(0, 0, PLNR_BLESSING);
  Y(planner_hinsert)
  (p2, spre, &fl, 2u);
  {
    size_t len = strlen(wis);
    char *bad = (char *)Y(malloc)(len + 1);
    strcpy(bad, wis);
    bad[len - 5] = '\0'; /* cut inside the last entry */
    CU_ASSERT_FALSE(import_from_string(p2, bad));
    Y(free)
    (bad);
  }
  q = mkflags(0, 0, 0);
  CU_ASSERT_PTR_NOT_NULL(Y(planner_hlookup)(p2, spre, &q)); /* restored */

  /* (b) unknown solver name -> reject (same length keeps the shape valid) */
  {
    char *bad = (char *)Y(malloc)(strlen(wis) + 1);
    char *at;
    strcpy(bad, wis);
    at = strstr(bad, "reg_mock_nfft_b");
    CU_ASSERT_PTR_NOT_NULL_FATAL(at);
    at[strlen("reg_mock_nfft_")] = 'z'; /* b -> z */
    CU_ASSERT_FALSE(import_from_string(p2, bad));
    Y(free)
    (bad);
  }

  /* (c) config-signature mismatch: planner with fewer solvers refuses */
  {
    planner *p4 = Y(planner_create)();
    static struct solvtab_s small_tab[] = {SOLVTAB(reg_mock_nfft_a), SOLVTAB_END};
    Y(solvtab_exec)
    (small_tab, p4);
    CU_ASSERT_FALSE(import_from_string(p4, wis));
    Y(planner_destroy)
    (p4);
  }

  /* (d) garbage preamble */
  CU_ASSERT_FALSE(import_from_string(p2, "(bogus-1.0 junk)"));

  /* (g) preamble naming a different-precision store: rejected, and existing
   * wisdom survives. The precision prefix precedes the "_wisdom" token, so
   * flipping its last char breaks the match with this build's token. */
  {
    char *bad = (char *)Y(malloc)(strlen(wis) + 1);
    char *w;
    strcpy(bad, wis);
    w = strstr(bad, "_wisdom");
    CU_ASSERT_PTR_NOT_NULL_FATAL(w);
    w[-1] = (w[-1] == 'f') ? 'l' : 'f';
    CU_ASSERT_FALSE(import_from_string(p2, bad));
    Y(free)
    (bad);
    q = mkflags(0, 0, 0);
    CU_ASSERT_PTR_NOT_NULL(Y(planner_hlookup)(p2, spre, &q)); /* restored */
  }

  /* (e) registrar name longer than 63 bytes: the parse must stay aligned and
   * fail on name resolution. */
  {
    char *bad = (char *)Y(malloc)(strlen(wis) + 128);
    char longname[71];
    const char *nl = strchr(wis, '\n');
    size_t plen;
    int i;
    CU_ASSERT_PTR_NOT_NULL_FATAL(nl);
    plen = (size_t)(nl - wis);
    for (i = 0; i < 70; i++)
      longname[i] = 'a';
    longname[70] = '\0';
    memcpy(bad, wis, plen);
    sprintf(bad + plen, "\n  (%s 0 #x0 #x0 #x0 #x0 #x0 #x0 #x0)\n)", longname);
    CU_ASSERT_FALSE(import_from_string(p2, bad));
    Y(free)
    (bad);
  }

  /* (f) bounds-invariant violation: a feasible entry with !LEQ(l, u) is
   * malformed -- a crafted file must not plant an entry that later strictly
   * subsumes a legitimate blessing */
  {
    char *bad = (char *)Y(malloc)(strlen(wis) + 128);
    const char *nl = strchr(wis, '\n');
    size_t plen;
    CU_ASSERT_PTR_NOT_NULL_FATAL(nl);
    plen = (size_t)(nl - wis);
    memcpy(bad, wis, plen);
    /* l = 0x10, u = 0x0: LEQ(0x10, 0x0) is false */
    sprintf(bad + plen,
            "\n  (reg_mock_nfft_b 0 #x10 #x0 #x0 #x1 #x2 #x3 #x4)\n)");
    CU_ASSERT_FALSE(import_from_string(p2, bad));
    q = mkflags(0, 0, 0);
    CU_ASSERT_PTR_NOT_NULL(Y(planner_hlookup)(p2, spre, &q)); /* restored */
    Y(free)
    (bad);
  }

  Y(free)
  (wis);
  Y(planner_destroy)
  (p1);
  Y(planner_destroy)
  (p2);
}

void Y(check_planner_tensor_basic)(void) {
  tensor *t, *a, *c;

  /* rank 0: the scalar identity */
  t = Y(tensor_create)(0);
  CU_ASSERT_EQUAL(t->rnk, 0);
  CU_ASSERT_EQUAL(Y(tensor_sz_in)(t), (INT)1);
  CU_ASSERT_EQUAL(Y(tensor_sz_out)(t), (INT)1);
  CU_ASSERT(Y(tensor_squarep)(t));
  CU_ASSERT(Y(tensor_kosherp)(t));
  Y(tensor_destroy)
  (t);

  /* a 2-factor rectangular operator: (5x4) kron (7x6) */
  t = Y(tensor_create)(2);
  t->dims[0].n_in = 4;
  t->dims[0].is = 1;
  t->dims[0].n_out = 5;
  t->dims[0].os = 1;
  t->dims[1].n_in = 6;
  t->dims[1].is = 4;
  t->dims[1].n_out = 7;
  t->dims[1].os = 5;
  CU_ASSERT(Y(tensor_kosherp)(t));
  CU_ASSERT_FALSE(Y(tensor_squarep)(t));
  CU_ASSERT_EQUAL(Y(tensor_sz_in)(t), (INT)24);
  CU_ASSERT_EQUAL(Y(tensor_sz_out)(t), (INT)35);

  /* copy is deep and equal; equal is fieldwise */
  c = Y(tensor_copy)(t);
  CU_ASSERT(Y(tensor_equal)(t, c));
  c->dims[1].n_out = 8;
  CU_ASSERT_FALSE(Y(tensor_equal)(t, c));
  Y(tensor_destroy)
  (c);

  /* adjoint swaps sizes and strides per factor; involution restores */
  a = Y(tensor_adjoint)(t);
  CU_ASSERT_EQUAL(a->dims[0].n_in, (INT)5);
  CU_ASSERT_EQUAL(a->dims[0].n_out, (INT)4);
  CU_ASSERT_EQUAL(Y(tensor_sz_in)(a), (INT)35);
  CU_ASSERT_EQUAL(Y(tensor_sz_out)(a), (INT)24);
  c = Y(tensor_adjoint)(a);
  CU_ASSERT(Y(tensor_equal)(c, t));
  Y(tensor_destroy)
  (c);
  Y(tensor_destroy)
  (a);

  /* the iodim case as a state: square constructor */
  {
    mvdim d = Y(mvdim_square)(9, 2, 3);
    CU_ASSERT_EQUAL(d.n_in, (INT)9);
    CU_ASSERT_EQUAL(d.n_out, (INT)9);
    CU_ASSERT_EQUAL(d.is, (INT)2);
    CU_ASSERT_EQUAL(d.os, (INT)3);
  }

  /* kosherp rejects a zero-sized factor */
  t->dims[0].n_in = 0;
  CU_ASSERT_FALSE(Y(tensor_kosherp)(t));
  t->dims[0].n_in = 4;

  Y(tensor_destroy)
  (t);
}

void Y(check_planner_tensor_canonical)(void) {
  tensor *t, *u, *ct, *cu;
  md5 m1, m2;
  char buf[128];
  printer *p;

  /* permuted factors + an interleaved 1x1 unit compress to the same form */
  t = Y(tensor_create)(3);
  t->dims[0].n_in = 1;
  t->dims[0].is = 99;
  t->dims[0].n_out = 1;
  t->dims[0].os = 42; /* unit factor: dropped */
  t->dims[1].n_in = 6;
  t->dims[1].is = 4;
  t->dims[1].n_out = 7;
  t->dims[1].os = 5;
  t->dims[2].n_in = 4;
  t->dims[2].is = 1;
  t->dims[2].n_out = 5;
  t->dims[2].os = 1;
  u = Y(tensor_create)(2);
  u->dims[0].n_in = 4;
  u->dims[0].is = 1;
  u->dims[0].n_out = 5;
  u->dims[0].os = 1;
  u->dims[1].n_in = 6;
  u->dims[1].is = 4;
  u->dims[1].n_out = 7;
  u->dims[1].os = 5;

  ct = Y(tensor_compress)(t);
  cu = Y(tensor_compress)(u);
  CU_ASSERT_EQUAL(ct->rnk, 2);
  CU_ASSERT(Y(tensor_equal)(ct, cu));

  /* a 1xn factor is not unit and survives compression */
  {
    tensor *v = Y(tensor_create)(1);
    tensor *cv;
    v->dims[0].n_in = 1;
    v->dims[0].is = 1;
    v->dims[0].n_out = 8;
    v->dims[0].os = 1;
    cv = Y(tensor_compress)(v);
    CU_ASSERT_EQUAL(cv->rnk, 1);
    CU_ASSERT_EQUAL(cv->dims[0].n_out, (INT)8);
    Y(tensor_destroy)
    (cv);
    Y(tensor_destroy)
    (v);
  }

  /* contiguous merge (FFTW parity, rectangular rule): both sides nest ->
   * sizes multiply per side, the inner factor's strides survive */
  {
    tensor *v = Y(tensor_create)(2);
    tensor *cv;
    v->dims[0].n_in = 6;
    v->dims[0].is = 1;
    v->dims[0].n_out = 4;
    v->dims[0].os = 1;
    v->dims[1].n_in = 2;
    v->dims[1].is = 6; /* = is_b * n_in_b  = 1 * 6 */
    v->dims[1].n_out = 3;
    v->dims[1].os = 4; /* = os_b * n_out_b = 1 * 4 */
    cv = Y(tensor_compress_contiguous)(v);
    CU_ASSERT_EQUAL(cv->rnk, 1);
    CU_ASSERT_EQUAL(cv->dims[0].n_in, (INT)12);
    CU_ASSERT_EQUAL(cv->dims[0].is, (INT)1);
    CU_ASSERT_EQUAL(cv->dims[0].n_out, (INT)12);
    CU_ASSERT_EQUAL(cv->dims[0].os, (INT)1);
    Y(tensor_destroy)
    (cv);

    /* breaking the input-side nesting forbids the merge */
    v->dims[1].is = 8;
    cv = Y(tensor_compress_contiguous)(v);
    CU_ASSERT_EQUAL(cv->rnk, 2);
    Y(tensor_destroy)
    (cv);
    Y(tensor_destroy)
    (v);
  }

  /* total order: sign-differing factors that tie on all absolute keys must
   * still canonicalise deterministically */
  {
    tensor *a = Y(tensor_create)(2);
    tensor *b = Y(tensor_create)(2);
    tensor *ca, *cb;
    md5 ma, mb;
    /* P = {2, +3, 2, +3}, Q = {2, -3, 2, +3}: equal under every absolute
     * key; only the signed-is tiebreaker separates them */
    a->dims[0].n_in = 2;
    a->dims[0].is = 3;
    a->dims[0].n_out = 2;
    a->dims[0].os = 3;
    a->dims[1].n_in = 2;
    a->dims[1].is = -3;
    a->dims[1].n_out = 2;
    a->dims[1].os = 3;
    b->dims[0] = a->dims[1];
    b->dims[1] = a->dims[0];
    ca = Y(tensor_compress)(a);
    cb = Y(tensor_compress)(b);
    CU_ASSERT(Y(tensor_equal)(ca, cb));
    Y(md5_begin)
    (&ma);
    Y(tensor_md5)
    (&ma, ca);
    Y(md5_end)
    (&ma);
    Y(md5_begin)
    (&mb);
    Y(tensor_md5)
    (&mb, cb);
    Y(md5_end)
    (&mb);
    CU_ASSERT(ma.s[0] == mb.s[0] && ma.s[1] == mb.s[1] && ma.s[2] == mb.s[2] && ma.s[3] == mb.s[3]);
    /* the signed tiebreaker is ascending: the negative-is factor first */
    CU_ASSERT_EQUAL(ca->dims[0].is, (INT)-3);
    Y(tensor_destroy)
    (ca);
    Y(tensor_destroy)
    (cb);
    Y(tensor_destroy)
    (a);
    Y(tensor_destroy)
    (b);
  }

  /* canonical forms hash equal */
  Y(md5_begin)
  (&m1);
  Y(tensor_md5)
  (&m1, ct);
  Y(md5_end)
  (&m1);
  Y(md5_begin)
  (&m2);
  Y(tensor_md5)
  (&m2, cu);
  Y(md5_end)
  (&m2);
  CU_ASSERT(m1.s[0] == m2.s[0] && m1.s[1] == m2.s[1] && m1.s[2] == m2.s[2] && m1.s[3] == m2.s[3]);

  /* the hash is direction-sensitive: transposing one factor changes it */
  cu->dims[0].n_in = 5;
  cu->dims[0].n_out = 4; /* transpose in place */
  Y(md5_begin)
  (&m2);
  Y(tensor_md5)
  (&m2, cu);
  Y(md5_end)
  (&m2);
  CU_ASSERT(m1.s[0] != m2.s[0] || m1.s[1] != m2.s[1] || m1.s[2] != m2.s[2] || m1.s[3] != m2.s[3]);

  /* printable S-expression, canonical order (the
   * min(|is|,|os|)=4 factor precedes the min=1 factor) */
  p = Y(printer_create_str)(buf);
  Y(tensor_print)
  (ct, p);
  Y(printer_destroy)
  (p);
  CU_ASSERT_STRING_EQUAL(buf, "(tensor 2 (6 4 7 5) (4 1 5 1))");

  /* rank 0 prints bare */
  {
    tensor *z = Y(tensor_create)(0);
    p = Y(printer_create_str)(buf);
    Y(tensor_print)
    (z, p);
    Y(printer_destroy)
    (p);
    CU_ASSERT_STRING_EQUAL(buf, "(tensor 0)");
    Y(tensor_destroy)
    (z);
  }

  Y(tensor_destroy)
  (ct);
  Y(tensor_destroy)
  (cu);
  Y(tensor_destroy)
  (t);
  Y(tensor_destroy)
  (u);
}

/* Trinity mocks */

typedef struct
{
  problem super;
  int payload;
} mock_problem;

static void mock_prb_hash(const problem *p, md5 *ctx) {
  Y(md5_put_str)
  (ctx, "mock");
  Y(md5_put_int)
  (ctx, ((const mock_problem *)p)->payload);
}

static void mock_prb_print(const problem *p, printer *pr) {
  pr->print(pr, "(mock-problem %d)", ((const mock_problem *)p)->payload);
}

static const problem_adt mock_prb_adt = {NFFT_PROBLEM_NFFT, mock_prb_hash, mock_prb_print, 0};

static problem *mk_mock_problem(int payload) {
  mock_problem *mp = (mock_problem *)Y(problem_create)(sizeof(mock_problem), &mock_prb_adt);
  mp->payload = payload;
  return &mp->super;
}

static int mock_plans_alive = 0;

typedef struct
{
  plan super;
  const char *nam;
  int *sleepy_calls;
} mock_plan;

static void mock_plan_apply(const plan *ego, const problem *p) {
  UNUSED(ego);
  UNUSED(p);
}
static void mock_plan_awake(plan *ego, int w) {
  mock_plan *mpl = (mock_plan *)ego;
  if (w == PLNR_SLEEPY && mpl->sleepy_calls)
    (*mpl->sleepy_calls)++;
}
static void mock_plan_print(const plan *ego, printer *pr) {
  pr->print(pr, "(%s)", ((const mock_plan *)ego)->nam);
}
static void mock_plan_destroy(plan *ego) {
  UNUSED(ego);
  mock_plans_alive--;
}
static const plan_adt mock_plan_adt = {mock_plan_apply, mock_plan_awake, mock_plan_print, mock_plan_destroy};

typedef struct
{
  solver super;
  double pcost;
  unsigned gate; /* decline if PLNR_L(pl) & gate */
  int decline;   /* always decline */
  int *mkplan_calls;
  const char *nam;
} mock_solver;

static const char *mock_force_decline = 0; /* registrar gone stale */

static plan *mock_mkplan(const solver *ego, const problem *p, planner *pl) {
  const mock_solver *ms = (const mock_solver *)ego;
  mock_plan *mpl;
  UNUSED(p);
  if (ms->mkplan_calls)
    (*ms->mkplan_calls)++;
  if (ms->decline || (PLNR_L(pl) & ms->gate) || (mock_force_decline && strcmp(ms->nam, mock_force_decline) == 0))
    return 0;
  mpl = (mock_plan *)Y(plan_create)(sizeof(mock_plan), &mock_plan_adt);
  mpl->super.pcost = ms->pcost;
  mpl->nam = ms->nam;
  mpl->sleepy_calls = 0;
  mock_plans_alive++;
  return &mpl->super;
}
static const solver_adt mock_slv_adt = {NFFT_PROBLEM_NFFT, 0, mock_mkplan};

static void reg_mock(planner *pl, double pcost, unsigned gate, int decline,
                     int *calls, const char *nam) {
  mock_solver *ms = (mock_solver *)Y(solver_create)(sizeof(mock_solver), &mock_slv_adt);
  ms->pcost = pcost;
  ms->gate = gate;
  ms->decline = decline;
  ms->mkplan_calls = calls;
  ms->nam = nam;
  Y(planner_register_solver)
  (pl, &ms->super);
}

static void plan_name_to(plan *pln, char *buf) {
  printer *pr = Y(printer_create_str)(buf);
  pr->print(pr, "%p", pln);
  Y(printer_destroy)
  (pr);
}

void Y(check_planner_trinity_mkplan)(void) {
  planner *pl = Y(planner_create)();
  int c_exp = 0, c_cheap = 0, c_tie = 0, c_decl = 0, c_gated = 0;
  problem *prb;
  plan *pln;
  char buf[64];

  /* registration order; the kind chain iterates in reverse */
  reg_mock(pl, 100.0, 0, 0, &c_exp, "mock-expensive");
  reg_mock(pl, 10.0, 0, 0, &c_cheap, "mock-cheap");
  reg_mock(pl, 10.0, 0, 0, &c_tie, "mock-tie");
  reg_mock(pl, 1.0, 0, 1, &c_decl, "mock-decline");
  reg_mock(pl, 0.5, PLNR_NO_DIRECT, 0, &c_gated, "mock-gated");

  /* ungated search: cheapest applicable (0.5) wins */
  prb = mk_mock_problem(1);
  pln = Y(planner_mkplan)(pl, prb);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  plan_name_to(pln, buf);
  CU_ASSERT_STRING_EQUAL(buf, "(mock-gated)");
  /* every solver of the kind was consulted exactly once */
  CU_ASSERT_EQUAL(c_exp, 1);
  CU_ASSERT_EQUAL(c_cheap, 1);
  CU_ASSERT_EQUAL(c_tie, 1);
  CU_ASSERT_EQUAL(c_decl, 1);
  CU_ASSERT_EQUAL(c_gated, 1);
  /* losers were destroyed: only the winner lives */
  CU_ASSERT_EQUAL(mock_plans_alive, 1);
  Y(plan_destroy)
  (pln);
  CU_ASSERT_EQUAL(mock_plans_alive, 0);
  Y(problem_destroy)
  (prb);

  /* gated search: NO_DIRECT filters the 0.5 solver; the 10.0 tie keeps
   * the earlier-encountered candidate (reverse order: tie before cheap) */
  pl->flags.l = PLNR_NO_DIRECT;
  pl->flags.u = PLNR_NO_DIRECT;
  prb = mk_mock_problem(2);
  pln = Y(planner_mkplan)(pl, prb);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  plan_name_to(pln, buf);
  CU_ASSERT_STRING_EQUAL(buf, "(mock-tie)");
  Y(plan_destroy)
  (pln);
  Y(problem_destroy)
  (prb);

  Y(planner_destroy)
  (pl);
}

void Y(check_planner_trinity_wisdom_memo)(void) {
  planner *pl = Y(planner_create)();
  int c_a = 0, c_b = 0, c_d = 0;
  problem *prb;
  plan *pln;
  char buf[64];
  unsigned nelem;

  reg_mock(pl, 5.0, 0, 0, &c_a, "mock-a");
  reg_mock(pl, 50.0, 0, 0, &c_b, "mock-b");

  prb = mk_mock_problem(7);
  pln = Y(planner_mkplan)(pl, prb);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  CU_ASSERT_EQUAL(c_a, 1);
  CU_ASSERT_EQUAL(c_b, 1);
  nelem = pl->htab_unblessed.nelem;
  CU_ASSERT_EQUAL(nelem, 1u); /* memoised, unblessed */
  CU_ASSERT_EQUAL(pl->htab_blessed.nelem, 0u);
  Y(plan_destroy)
  (pln);
  Y(problem_destroy)
  (prb);

  /* equal problem: wisdom hit -> only the winner's mkplan runs again */
  prb = mk_mock_problem(7);
  pln = Y(planner_mkplan)(pl, prb);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  plan_name_to(pln, buf);
  CU_ASSERT_STRING_EQUAL(buf, "(mock-a)");
  CU_ASSERT_EQUAL(c_a, 2);
  CU_ASSERT_EQUAL(c_b, 1);                          /* loser not re-consulted */
  CU_ASSERT_EQUAL(pl->htab_unblessed.nelem, nelem); /* no growth */
  Y(plan_destroy)
  (pln);
  Y(problem_destroy)
  (prb);
  Y(planner_destroy)
  (pl);

  /* all-decline kind: infeasibility is memoised */
  pl = Y(planner_create)();
  reg_mock(pl, 1.0, 0, 1, &c_d, "mock-decline-only");
  prb = mk_mock_problem(9);
  CU_ASSERT_PTR_NULL(Y(planner_mkplan)(pl, prb));
  CU_ASSERT_EQUAL(c_d, 1);
  CU_ASSERT_EQUAL(pl->htab_unblessed.nelem, 1u);
  Y(problem_destroy)
  (prb);
  prb = mk_mock_problem(9);
  CU_ASSERT_PTR_NULL(Y(planner_mkplan)(pl, prb));
  CU_ASSERT_EQUAL(c_d, 1); /* memoised: not consulted again */
  Y(problem_destroy)
  (prb);
  Y(planner_destroy)
  (pl);

  /* stale hit heals: the memoised winner stops applying -> fallthrough
   * search wins with the runner-up and replaces the stale entry
   * (mutual-subsumption insert) */
  pl = Y(planner_create)();
  c_a = 0;
  c_b = 0;
  reg_mock(pl, 5.0, 0, 0, &c_a, "mock-a");
  reg_mock(pl, 50.0, 0, 0, &c_b, "mock-b");
  prb = mk_mock_problem(11);
  pln = Y(planner_mkplan)(pl, prb);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln); /* c_a=1 c_b=1, winner mock-a */
  Y(plan_destroy)
  (pln);
  Y(problem_destroy)
  (prb);
  mock_force_decline = "mock-a"; /* the winner goes stale */
  prb = mk_mock_problem(11);
  pln = Y(planner_mkplan)(pl, prb); /* hit->rerun (c_a=2, NULL) -> search
                                       (c_a=3, c_b=2) -> mock-b heals */
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  plan_name_to(pln, buf);
  CU_ASSERT_STRING_EQUAL(buf, "(mock-b)");
  CU_ASSERT_EQUAL(pl->htab_unblessed.nelem, 1u); /* replaced, not grown */
  Y(plan_destroy)
  (pln);
  Y(problem_destroy)
  (prb);
  prb = mk_mock_problem(11); /* healed: direct hit on mock-b, no search */
  pln = Y(planner_mkplan)(pl, prb);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  plan_name_to(pln, buf);
  CU_ASSERT_STRING_EQUAL(buf, "(mock-b)");
  CU_ASSERT_EQUAL(c_a, 3);
  CU_ASSERT_EQUAL(c_b, 3);
  Y(plan_destroy)
  (pln);
  Y(problem_destroy)
  (prb);
  mock_force_decline = 0;
  Y(planner_destroy)
  (pl);
}

void Y(check_planner_trinity_print)(void) {
  planner *pl = Y(planner_create)();
  int c = 0, sleepy = 0;
  problem *prb = mk_mock_problem(4);
  plan *pln;
  char buf[64];
  printer *pr;

  /* %P prints through the problem's adt; NULL prints (null) */
  pr = Y(printer_create_str)(buf);
  pr->print(pr, "%P", prb);
  Y(printer_destroy)
  (pr);
  CU_ASSERT_STRING_EQUAL(buf, "(mock-problem 4)");
  pr = Y(printer_create_str)(buf);
  pr->print(pr, "%p/%P", (plan *)0, (problem *)0);
  Y(printer_destroy)
  (pr);
  CU_ASSERT_STRING_EQUAL(buf, "(null)/(null)");

  /* awake lifecycle: destroy of an AWAKE plan sleeps it exactly once */
  reg_mock(pl, 1.0, 0, 0, &c, "mock-w");
  pln = Y(planner_mkplan)(pl, prb);
  CU_ASSERT_PTR_NOT_NULL_FATAL(pln);
  ((mock_plan *)pln)->sleepy_calls = &sleepy;
  Y(plan_awake)
  (pln, PLNR_AWAKE);
  Y(plan_awake)
  (pln, PLNR_AWAKE); /* idempotent */
  CU_ASSERT_EQUAL(sleepy, 0);
  Y(plan_destroy)
  (pln); /* must awake(SLEEPY) first */
  CU_ASSERT_EQUAL(sleepy, 1);

  Y(problem_destroy)
  (prb);
  Y(planner_destroy)
  (pl);

  /* the global planner: lazy create, destroy, recreate */
  CU_ASSERT_PTR_NOT_NULL(Y(the_planner)());
  CU_ASSERT_PTR_EQUAL(Y(the_planner)(), Y(the_planner)());
  Y(the_planner_destroy)
  ();
  CU_ASSERT_PTR_NOT_NULL(Y(the_planner)());
  Y(the_planner_destroy)
  ();
  Y(the_planner_destroy)
  (); /* safe when absent */
}

/* Measurement, candidates, blessing */

static volatile double spin_sink; /* defeats optimisation of the spin loop */

typedef struct
{
  plan super;
  long work; /* spin iterations per apply */
  int *applies;
  int *awakes;
} spin_plan;

static void spin_apply(const plan *ego, const problem *p) {
  const spin_plan *sp = (const spin_plan *)ego;
  double acc = 0.0;
  long i;
  UNUSED(p);
  if (sp->applies)
    (*sp->applies)++;
  for (i = 0; i < sp->work; i++)
    acc += (double)(i & 7);
  spin_sink = acc;
}
static void spin_awake(plan *ego, int w) {
  spin_plan *sp = (spin_plan *)ego;
  UNUSED(w);
  if (sp->awakes)
    (*sp->awakes)++;
}
static void spin_print(const plan *ego, printer *pr) {
  UNUSED(ego);
  pr->print(pr, "(spin)");
}
static const plan_adt spin_adt = {spin_apply, spin_awake, spin_print, 0};

static plan *mk_spin(long work, int *applies, int *awakes) {
  spin_plan *sp = (spin_plan *)Y(plan_create)(sizeof(spin_plan), &spin_adt);
  sp->work = work;
  sp->applies = applies;
  sp->awakes = awakes;
  return &sp->super;
}

void Y(check_planner_measure)(void) {
  int applies = 0, awakes = 0;
  problem *prb = mk_mock_problem(21);
  plan *light = mk_spin(1000L, &applies, &awakes);
  plan *heavy = mk_spin(400000L, 0, 0);
  double tl, th;

  tl = Y(plan_measure_cost)(light, prb);
  th = Y(plan_measure_cost)(heavy, prb);
  if (tl < 0.0 || th < 0.0) {
    /* no usable clock: both must agree; measured mode degrades */
    CU_ASSERT(tl < 0.0 && th < 0.0);
  } else {
#if defined(HAVE_TICK_COUNTER)
    /* Coverage only: the tick path runs and reads strictly positive on
     * counter-equipped CI. The ratio below is the behavioural assertion. */
    CU_ASSERT(tl > 0.0 && th > 0.0);
#endif
    /* 400x the work must measure at least 2x slower. Cost units differ by
     * mode (raw ticks or slow-seconds) but the ratio does not, and
     * min-of-repeats makes it robust to scheduler noise. */
    CU_ASSERT(th > 2.0 * tl);
    /* the doubling loop ran the light plan more than once */
    CU_ASSERT(applies > 1);
  }
  /* measurement never touches wakefulness: the caller owns awake state */
  CU_ASSERT_EQUAL(awakes, 0);

  Y(plan_destroy)
  (light);
  Y(plan_destroy)
  (heavy);
  Y(problem_destroy)
  (prb);
}

void Y(check_planner_candidates)(void) {
  planner *pl = Y(planner_create)();
  problem *prb = mk_mock_problem(23);
  plan *plans[8];
  unsigned ndx[8];
  char buf[64];
  int n, i;

  reg_mock(pl, 100.0, 0, 0, 0, "cand-a");
  reg_mock(pl, 10.0, 0, 1, 0, "cand-decline");
  reg_mock(pl, 1.0, PLNR_NO_DIRECT, 0, 0, "cand-gated");

  /* ungated: decliner drops out; reverse registration order pins
   * cand-gated first, cand-a second; distinct descriptor indices */
  n = Y(planner_candidates)(pl, prb, plans, ndx, 8);
  CU_ASSERT_EQUAL_FATAL(n, 2);
  plan_name_to(plans[0], buf);
  CU_ASSERT_STRING_EQUAL(buf, "(cand-gated)");
  plan_name_to(plans[1], buf);
  CU_ASSERT_STRING_EQUAL(buf, "(cand-a)");
  CU_ASSERT(ndx[0] < pl->nslvdesc);
  CU_ASSERT(ndx[1] < pl->nslvdesc);
  CU_ASSERT(ndx[0] != ndx[1]);
  /* no store mutation: candidates is enumeration, not search */
  CU_ASSERT_EQUAL(pl->htab_unblessed.nelem, 0u);
  CU_ASSERT_EQUAL(pl->htab_blessed.nelem, 0u);
  for (i = 0; i < n; i++)
    Y(plan_destroy)
  (plans[i]);

  /* gated: PLNR_NO_DIRECT filters cand-gated */
  pl->flags.l = PLNR_NO_DIRECT;
  pl->flags.u = PLNR_NO_DIRECT;
  n = Y(planner_candidates)(pl, prb, plans, ndx, 8);
  CU_ASSERT_EQUAL_FATAL(n, 1);
  plan_name_to(plans[0], buf);
  CU_ASSERT_STRING_EQUAL(buf, "(cand-a)");
  Y(plan_destroy)
  (plans[0]);

  Y(problem_destroy)
  (prb);
  Y(planner_destroy)
  (pl);
}

void Y(check_planner_bless)(void) {
  planner *pl = Y(planner_create)();
  problem *prb = mk_mock_problem(29);
  problem *prb2 = mk_mock_problem(31);
  plan *plans[4];
  unsigned ndx[4];
  md5sig sig;
  flags_t q;
  solution *sol;
  int n, i;

  reg_mock(pl, 5.0, 0, 0, 0, "bless-a");
  reg_mock(pl, 50.0, 0, 0, 0, "bless-b");

  /* measured-mode bounds: no PLNR_ESTIMATE anywhere */
  pl->flags.l = 0;
  pl->flags.u = 0;
  n = Y(planner_candidates)(pl, prb, plans, ndx, 4);
  CU_ASSERT_EQUAL_FATAL(n, 2);
  Y(planner_bless)
  (pl, prb, ndx[1]); /* candidate 1 "won the race" */
  CU_ASSERT_EQUAL(pl->htab_blessed.nelem, 1u);
  CU_ASSERT_EQUAL(pl->htab_unblessed.nelem, 0u);

  /* the measured entry answers an ESTIMATE-mode query
   * (PLNR_ESTIMATE lives in u only) */
  Y(problem_md5)
  (pl, prb, sig);
  q.l = 0;
  q.u = PLNR_ESTIMATE;
  q.timelimit_imp = 0;
  q.info = 0;
  q.slvndx = 0;
  sol = Y(planner_hlookup)(pl, sig, &q);
  CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
  CU_ASSERT_EQUAL(sol->flags.slvndx, ndx[1]);
  CU_ASSERT(sol->flags.info & PLNR_BLESSING);

  /* an estimate memo does not answer a measured query (other problem) */
  {
    md5sig s2;
    flags_t f2;
    Y(problem_md5)
    (pl, prb2, s2);
    f2.l = 0;
    f2.u = PLNR_ESTIMATE;
    f2.timelimit_imp = 0;
    f2.info = 0;
    f2.slvndx = 0;
    Y(planner_hinsert)
    (pl, s2, &f2, 0u); /* estimate memo, unblessed */
    q.l = 0;
    q.u = 0; /* measured query */
    CU_ASSERT_PTR_NULL(Y(planner_hlookup)(pl, s2, &q));
  }

  /* re-bless with the same winner: idempotent, no growth */
  Y(planner_bless)
  (pl, prb, ndx[1]);
  CU_ASSERT_EQUAL(pl->htab_blessed.nelem, 1u);

  /* re-bless with a different winner (re-measurement changed the outcome):
   * equal bounds replace via mutual subsumption, no growth */
  Y(planner_bless)
  (pl, prb, ndx[0]);
  CU_ASSERT_EQUAL(pl->htab_blessed.nelem, 1u);
  q.l = 0;
  q.u = 0;
  sol = Y(planner_hlookup)(pl, sig, &q);
  CU_ASSERT_PTR_NOT_NULL_FATAL(sol);
  CU_ASSERT_EQUAL(sol->flags.slvndx, ndx[0]);

  for (i = 0; i < n; i++)
    Y(plan_destroy)
  (plans[i]);
  Y(problem_destroy)
  (prb);
  Y(problem_destroy)
  (prb2);
  Y(planner_destroy)
  (pl);
}

/* Timelimit storage and the coarse clock reader */

void Y(check_planner_timelimit_default_and_set)(void) {
  planner *pl = Y(the_planner)();
  CK(Y(planner_timelimit)(pl) == -1.0);
  Y(planner_set_timelimit)
  (pl, 0.25);
  CK(Y(planner_timelimit)(pl) == 0.25);
  Y(planner_set_timelimit)
  (pl, -7.0);
  CK(Y(planner_timelimit)(pl) == -1.0);
  Y(planner_set_timelimit)
  (pl, -1.0);
}

void Y(check_planner_clock_now_monotonic)(void) {
  double t0 = Y(planner_clock_now)();
  CK(Y(planner_elapsed_seconds)(t0) >= 0.0);
}
