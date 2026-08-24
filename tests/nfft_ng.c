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

/* Planner-API (plan_ng) NFFT accuracy suite, the counterpart of tests/nfft.c:
 * same reference files, same online grid, same metric, same bounds, on the
 * geometry X(init) would have chosen.  The bound tables here are copies;
 * retune them together with tests/nfft.c.
 *
 * The gate flags pin the algorithm so the plan is deterministic:
 * NFFT_NO_FAST_NATIVE leaves only the direct NDFT, NFFT_NO_DIRECT only the
 * composed fast NFFT, neither lets the planner choose.
 *
 * Three legacy behaviours have no counterpart here: N[t] <= m (guards_ok in
 * nfft-nd.c declines, so NFFT_NO_DIRECT finds no plan and the case is skipped,
 * but only after fast_admits() agrees), odd N (legacy X(check) rejects it), and
 * N[t] == 1 (elided; all-unit is rank 0).  docs/agents/test-methodology.md
 * carries the full account.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "accuracy_log.h"
#include "nfft_ng.h"

#include "data/generated/nfft_native_testcases.h"

/* Same seed as tests/nfft.c: the online grids draw the same inputs. */
#define SEED 1234567890L

/* Largest rank the tables below use. */
#define NG_MAX_D 4

/* trafo / adjoint, matching native_testcase_t::kind. */
#define NG_TRAFO 0
#define NG_ADJOINT 1

/* Which algorithm the planner is allowed to pick. */
#define NG_DIRECT 0
#define NG_FAST 1
#define NG_AUTO 2

static const char *mode_name(int mode) {
  return mode == NG_DIRECT ? "plan_ng (no fast)"
                           : (mode == NG_FAST ? "plan_ng (no direct)"
                                              : "plan_ng (auto)");
}

/* The accuracy report's speed axis keys off this name. */
static const char *trafo_name(int mode, int kind) {
  if (mode == NG_DIRECT)
    return kind == NG_TRAFO ? "trafo_direct" : "adjoint_direct";
  return kind == NG_TRAFO ? "trafo" : "adjoint";
}

static unsigned mode_flags(int mode) {
  switch (mode) {
  case NG_DIRECT:
    return NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE;
  case NG_FAST:
    return NFFT_ESTIMATE | NFFT_NO_DIRECT;
  default:
    return NFFT_ESTIMATE;
  }
}

/* err_trafo_direct: the round-off floor. */
static R bound_direct(void) {
  return K(56.0) * Y(float_property)(NFFT_EPSILON);
}

/* err_trafo, for cut-off m and the smallest oversampling factor s. */
static R bound_fast(R m, R s) {
  R eps = Y(float_property)(NFFT_EPSILON);
  R a;
  R b;
  R err;
#if defined(GAUSSIAN)
#if MANT_DIG == 113
  a = K(0.8);
  b = K(50.0);
#elif MANT_DIG == 64
  a = K(0.95);
  b = K(50.0);
#elif MANT_DIG == 53
  a = K(0.41);
  b = K(50.0);
#elif MANT_DIG == 24
  a = K(0.4);
  b = K(2000.0);
#else
  a = K(0.41);
  b = K(50.0);
#endif
  UNUSED(s);
  err = EXP(-m * KPI * (K(1.0) - K(1.0) / (K(2.0) * K(2.0) - K(1.0))));
#elif defined(B_SPLINE)
#if MANT_DIG == 113
  a = K(0.3);
  b = K(50.0);
#elif MANT_DIG == 64
  a = K(0.3);
  b = K(50.0);
#elif MANT_DIG == 53
  a = K(1.0);
  b = K(2000.0);
#elif MANT_DIG == 24
  a = K(0.4);
  b = K(2000.0);
#else
  a = K(1.0);
  b = K(2000.0);
#endif
  err = K(3000.0) * K(4.0) * POW(K(1.0) / (K(2.0) * s - K(1.0)), K(2.0) * m);
#elif defined(SINC_POWER)
#if MANT_DIG == 113
  a = K(0.3);
  b = K(50.0);
#elif MANT_DIG == 64
  a = K(0.3);
  b = K(50.0);
#elif MANT_DIG == 53
  a = K(1.0);
  b = K(2000.0);
#elif MANT_DIG == 24
  a = K(0.4);
  b = K(2000.0);
#else
  a = K(1.0);
  b = K(2000.0);
#endif
  err = (K(1.0) / (m - K(1.0))) *
        ((K(2.0) / (POW(s, K(2.0) * m))) +
         POW(s / (K(2.0) * s - K(1.0)), K(2.0) * m));
#elif defined(KAISER_BESSEL)
#if MANT_DIG == 113
  a = K(1.5);
  b = K(50.0);
#elif MANT_DIG == 64
  a = K(1.5);
  b = K(50.0);
#elif MANT_DIG == 53
  a = K(0.3);
  b = K(2100.0);
#elif MANT_DIG == 24
  a = K(0.4);
  b = K(2000.0);
#else
  a = K(0.3);
  b = K(2100.0);
#endif
  UNUSED(s);
  err = KPI * (SQRT(m) + m) * SQRT(SQRT(K(1.0) - K(1.0) / K(2.0))) *
        EXP(-K2PI * m * SQRT(K(1.0) - K(1.0) / K(2.0)));
#else
#error Unsupported window function.
#endif

  return FMAX(FMAX(a * err, b * eps), bound_direct());
}

/* n[t] = 2 * next_power_of_2(N[t]), what X(init) uses.  Unit axes keep the
 * legacy value too, so sigma -- and with it the bound -- stays legacy-equal. */
static void fill_n(int d, const INT *N, INT *n) {
  int t;
  for (t = 0; t < d; t++)
    n[t] = 2 * Y(next_power_of_2)(N[t]);
}

static R sigma_min(int d, const INT *N, const INT *n) {
  R s = (R)n[0] / (R)N[0];
  int t;
  for (t = 1; t < d; t++)
    s = FMIN(s, (R)n[t] / (R)N[t]);
  return s;
}

static INT prod(int d, const INT *N) {
  INT p = 1;
  int t;
  for (t = 0; t < d; t++)
    p *= N[t];
  return p;
}

/* guards_ok (nfft-nd.c) plus the window restriction: does the composed fast
 * solver accept this geometry?  Unit axes are elided before the guard sees
 * them; all-unit is rank 0, which no gate flag excludes. */
static int fast_admits(int d, const INT *N, const INT *n, int m, int window) {
  int t;
  if (window < NFFT_WINDOW_KAISER_BESSEL || window > NFFT_WINDOW_SINC_POWER) {
    for (t = 0; t < d; t++)
      if (N[t] > (INT)1)
        return 0; /* a live axis needs a fast solver, and none accepts it */
    return 1;     /* rank 0 */
  }
  for (t = 0; t < d; t++) {
    if (N[t] == (INT)1)
      continue;
    if (!(N[t] > (INT)m))
      return 0;
    if (!(n[t] > (INT)(2 * m + 2)))
      return 0;
    if (!(n[t] > N[t]))
      return 0;
  }
  return 1;
}

/* compare_trafo / compare_adjoint: max-norm output deviation over the 1-norm
 * of the input coefficients. */
static R relerr(const C *got, const C *ref, INT len, const C *in,
                   INT in_len) {
  R num = K(0.0), den = K(0.0);
  INT j;
  for (j = 0; j < len; j++) {
    R e = CABS(got[j] - ref[j]);
    if (e > num)
      num = e;
  }
  for (j = 0; j < in_len; j++)
    den += CABS(in[j]);
  return num == K(0.0) ? K(0.0) : num / den;
}

/* Plan, precompute, execute, copy the result to out.
 *   NG_TRAFO:   in is f_hat (NN), out is f (M)
 *   NG_ADJOINT: in is f (M),      out is f_hat (NN)
 * Planning is destructive, so in is written back after Y(precompute).
 * Returns 0 if no candidate survived the gate flags. */
static int run_plan(int d, const INT *N, const int *variant, INT M, INT NN,
                   unsigned planning, const R *x, const C *in, C *out,
                   int kind) {
  INT n[NG_MAX_D];
  Y(plan_ng) * p;
  R *xw;
  C *fh, *ff;
  const size_t in_bytes =
      (size_t)(kind == NG_TRAFO ? NN : M) * sizeof(C);

  /* Fixed-size geometry buffers here and in grade/report_decline. */
  if (d > NG_MAX_D) {
    CU_FAIL("rank exceeds NG_MAX_D");
    return 0;
  }

  fill_n(d, N, n);

  xw = (R *)Y(malloc)((size_t)d * (size_t)M * sizeof(R));
  fh = (C *)Y(malloc)((size_t)NN * sizeof(C));
  ff = (C *)Y(malloc)((size_t)M * sizeof(C));
  memcpy(xw, x, (size_t)d * (size_t)M * sizeof(R));
  memset(fh, 0, (size_t)NN * sizeof(C));
  memset(ff, 0, (size_t)M * sizeof(C));
  memcpy(kind == NG_TRAFO ? fh : ff, in, in_bytes);

  p = Y(plan_ng_guru)(d, N, variant, n, M, WINDOW_HELP_ESTIMATE_m,
                      Y(get_window_id)(), xw, fh, ff, 0u, planning);
  if (!p) {
    Y(free)(xw);
    Y(free)(fh);
    Y(free)(ff);
    return 0;
  }

  Y(precompute)(p);
  memcpy(kind == NG_TRAFO ? fh : ff, in, in_bytes);
  if (kind == NG_TRAFO)
    Y(execute)(p);
  else
    Y(execute_adjoint)(p);

  memcpy(out, kind == NG_TRAFO ? ff : fh,
         (size_t)(kind == NG_TRAFO ? M : NN) * sizeof(C));

  Y(plan_ng_destroy)(p);
  Y(free)(xw);
  Y(free)(fh);
  Y(free)(ff);
  return 1;
}

static void print_geometry(const char *label, int d, const INT *N, INT M) {
  int t;
  printf("%-34s d = %d, N = [", label, d);
  for (t = 0; t < d; t++)
    printf("%s%-5ld", t > 0 ? ", " : "", (long)N[t]);
  for (t = 0; t < NG_MAX_D - d; t++)
    printf("%s%-5s", "  ", "");
  printf("], M = %-5ld", (long)M);
#ifdef _OPENMP
  printf(", nthreads = %d", (int)Y(get_num_threads)());
#endif
}

static void log_case(const char *oracle, int d, const INT *N, INT M, int mode,
                   int kind, R err, R bound, int ok) {
  int Ni[NG_MAX_D];
  int t;
  for (t = 0; t < d; t++)
    Ni[t] = (int)N[t];
  accuracy_log_append("nfft_ng", oracle, d, Ni, (int)M, mode_name(mode),
                      trafo_name(mode, kind), (long double)err,
                      (long double)bound, ok);
}

/* Grade one computed result: 1 pass, 0 fail. */
static int grade(const char *label, const char *oracle, int d, const INT *N,
                    INT M, int mode, int kind, R err) {
  INT n[NG_MAX_D];
  R bound;
  int ok;
  fill_n(d, N, n);
  bound = (mode == NG_DIRECT) ? bound_direct()
                              : bound_fast((R)WINDOW_HELP_ESTIMATE_m,
                                              sigma_min(d, N, n));
  ok = err < bound;
  print_geometry(label, d, N, M);
  printf(", m = %2d, %-18s, %-14s -> %-4s " __FE__ " (" __FE__ ")\n",
         (int)WINDOW_HELP_ESTIMATE_m, mode_name(mode),
         trafo_name(mode, kind), ok ? "OK" : "FAIL", err, bound);
  log_case(oracle, d, N, M, mode, kind, err, bound, ok);
  return ok;
}

/* A geometry the fast solver declines.  Legal only when the guard predicate
 * agrees; otherwise it is a planner regression. */
static int report_decline(const char *label, int d, const INT *N, INT M,
                             int mode) {
  INT n[NG_MAX_D];
  int admits;
  fill_n(d, N, n);
  admits = fast_admits(d, N, n, WINDOW_HELP_ESTIMATE_m, Y(get_window_id)());
  print_geometry(label, d, N, M);
  printf(", m = %2d, %-18s -> %s\n", (int)WINDOW_HELP_ESTIMATE_m,
         mode_name(mode),
         admits ? "FAIL (no plan although the guard admits the geometry)"
                : "SKIP (fast guard declines: N <= m or n <= 2m+2)");
  CU_ASSERT(!admits);
  return admits ? 0 : -1;
}

int Y(nfft_ng_read_case)(const char *rel, int *d, INT **N, INT *NN, INT *M,
                         R **x, C **f_hat, C **f) {
  char path[4096];
  FILE *fp;
  int t;
  INT j, nn = 1;

  *N = NULL;
  *x = NULL;
  *f_hat = NULL;
  *f = NULL;

  snprintf(path, sizeof path, "%s/tests/%s", ABS_SRCDIR, rel);
  fp = fopen(path, "r");
  if (!fp)
    return 0;
  if (fscanf(fp, "%d", d) != 1)
    goto fail;

  *N = (INT *)Y(malloc)((size_t)*d * sizeof(INT));
  for (t = 0; t < *d; t++) {
    long v;
    if (fscanf(fp, "%ld", &v) != 1)
      goto fail;
    (*N)[t] = (INT)v;
    nn *= (INT)v;
  }
  {
    long v;
    if (fscanf(fp, "%ld", &v) != 1)
      goto fail;
    *M = (INT)v;
  }
  *NN = nn;

  *x = (R *)Y(malloc)((size_t)(*d * *M) * sizeof(R));
  for (j = 0; j < *d * *M; j++)
    if (fscanf(fp, __FI__, &((*x)[j])) != 1)
      goto fail;

  *f_hat = (C *)Y(malloc)((size_t)nn * sizeof(C));
  for (j = 0; j < nn; j++) {
    R re, im;
    if (fscanf(fp, __FI__ " " __FI__, &re, &im) != 2)
      goto fail;
    (*f_hat)[j] = re + II * im;
  }

  *f = (C *)Y(malloc)((size_t)*M * sizeof(C));
  for (j = 0; j < *M; j++) {
    R re, im;
    if (fscanf(fp, __FI__ " " __FI__, &re, &im) != 2)
      goto fail;
    (*f)[j] = re + II * im;
  }

  fclose(fp);
  return 1;

fail:
  fclose(fp);
  Y(free)(*N);
  Y(free)(*x);
  Y(free)(*f_hat);
  Y(free)(*f);
  *N = NULL;
  *x = NULL;
  *f_hat = NULL;
  *f = NULL;
  return 0;
}

static const char *case_name(const char *path) {
  const char *slash = strrchr(path, '/');
  return slash ? slash + 1 : path;
}

/* Every reference case of one dimension and kind, through all three modes. */
static void file_dim(int dim, int kind) {
  int i, ok = 1, ran = 0, declined = 0;

  printf("nfft_ng_%s%dd_file:\n", kind == NG_ADJOINT ? "adjoint_" : "", dim);

  for (i = 0; i < native_testcases_count; i++) {
    const native_testcase_t *tc = &native_testcases[i];
    int d, mode;
    INT *N, NN, M;
    R *x;
    C *f_hat, *f;

    if (tc->d != dim || tc->kind != kind)
      continue;
    if (tc->d > NG_MAX_D) {
      CU_FAIL("roster case of rank > NG_MAX_D");
      ok = 0;
      continue;
    }

    if (!Y(nfft_ng_read_case)(tc->filename, &d, &N, &NN, &M, &x, &f_hat, &f)) {
      printf("%-34s -> FAIL (cannot read)\n", case_name(tc->filename));
      CU_FAIL("nfft_ng_read_case failed");
      ok = 0;
      continue;
    }

    for (mode = NG_DIRECT; mode <= NG_AUTO; mode++) {
      C *out = (C *)Y(malloc)((size_t)(kind == NG_TRAFO ? M : NN) * sizeof(C));
      if (!run_plan(d, N, tc->variant, M, NN, mode_flags(mode), x,
                   kind == NG_TRAFO ? f_hat : f, out, kind)) {
        /* Only the fast-only mode may legally find no plan: the direct NDFT
         * and the rank-0 solver together serve every geometry. */
        if (mode == NG_FAST) {
          if (report_decline(case_name(tc->filename), d, N, M, mode) < 0)
            declined++;
          else
            ok = 0;
        } else {
          print_geometry(case_name(tc->filename), d, N, M);
          printf(", %-18s -> FAIL (no plan)\n", mode_name(mode));
          CU_FAIL("plan_ng_guru declined a case this mode must serve");
          ok = 0;
        }
        Y(free)(out);
        continue;
      }
      {
        R err = (kind == NG_TRAFO)
                    ? relerr(out, f, M, f_hat, NN)
                    : relerr(out, f_hat, NN, f, M);
        int r = grade(case_name(tc->filename), "file", d, N, M, mode, kind,
                         err);
        ok = MIN(ok, r);
        ran++;
      }
      Y(free)(out);
    }

    Y(free)(N);
    Y(free)(x);
    Y(free)(f_hat);
    Y(free)(f);
  }

  printf("  (%d graded, %d declined by the fast guard)\n", ran, declined);
  CU_ASSERT(ran > 0);
  CU_ASSERT(ok);
  Y(the_planner_destroy)();
}

void Y(check_nfft_ng_1d_file)(void) { file_dim(1, NG_TRAFO); }
void Y(check_nfft_ng_2d_file)(void) { file_dim(2, NG_TRAFO); }
void Y(check_nfft_ng_3d_file)(void) { file_dim(3, NG_TRAFO); }
void Y(check_nfft_ng_4d_file)(void) { file_dim(4, NG_TRAFO); }
void Y(check_nfft_ng_adjoint_1d_file)(void) { file_dim(1, NG_ADJOINT); }
void Y(check_nfft_ng_adjoint_2d_file)(void) { file_dim(2, NG_ADJOINT); }
void Y(check_nfft_ng_adjoint_3d_file)(void) { file_dim(3, NG_ADJOINT); }
void Y(check_nfft_ng_adjoint_4d_file)(void) { file_dim(4, NG_ADJOINT); }

/* Nodes and input coefficients as setup_online / setup_adjoint_online draw
 * them. */
static void draw(int d, INT M, INT in_len, R *x, C *in) {
  INT j;
  Y(srand48)(SEED);
  for (j = 0; j < M * (INT)d; j++)
    x[j] = Y(drand48)() - K(0.5);
  for (j = 0; j < in_len; j++)
    in[j] = (Y(drand48)() - K(0.5)) + (Y(drand48)() - K(0.5)) * II;
}

/* Run the direct NDFT for the oracle, then the fast NFFT on the same input.
 * Returns 1 pass, 0 fail, -1 skipped. */
static int online_case(const char *label, int d, const INT *N,
                          const int *variant, INT M, int kind) {
  const INT NN = prod(d, N);
  const INT in_len = (kind == NG_TRAFO) ? NN : M;
  const INT out_len = (kind == NG_TRAFO) ? M : NN;
  R *x = (R *)Y(malloc)((size_t)d * (size_t)M * sizeof(R));
  C *in = (C *)Y(malloc)((size_t)in_len * sizeof(C));
  C *ref = (C *)Y(malloc)((size_t)out_len * sizeof(C));
  C *got = (C *)Y(malloc)((size_t)out_len * sizeof(C));
  int r;

  draw(d, M, in_len, x, in);

  if (!run_plan(d, N, variant, M, NN, mode_flags(NG_DIRECT), x, in, ref,
               kind)) {
    print_geometry(label, d, N, M);
    printf(" -> FAIL (no direct NDFT plan)\n");
    CU_FAIL("the direct NDFT must serve every geometry");
    r = 0;
    goto done;
  }

  if (!run_plan(d, N, variant, M, NN, mode_flags(NG_FAST), x, in, got,
               kind)) {
    r = report_decline(label, d, N, M, NG_FAST);
    goto done;
  }

  r = grade(label, "online", d, N, M, NG_FAST, kind,
               relerr(got, ref, out_len, in, in_len));

done:
  Y(free)(x);
  Y(free)(in);
  Y(free)(ref);
  Y(free)(got);
  return r;
}

static void online_uniform(const char *title, int d, const INT *Nvals,
                              int count, INT M, int kind) {
  int i, ok = 1, ran = 0;
  printf("%s:\n", title);
  for (i = 0; i < count; i++) {
    INT N[NG_MAX_D];
    int t, r;
    for (t = 0; t < d; t++)
      N[t] = Nvals[i];
    r = online_case("nfft_ng_online", d, N, 0, M, kind);
    if (r >= 0) {
      ok = MIN(ok, r);
      ran++;
    }
  }
  CU_ASSERT(ran > 0);
  CU_ASSERT(ok);
  Y(the_planner_destroy)();
}

/* The legacy online grids, exhaustive gating included. */
static const INT online_1d[] = {50, 100, 200, 500
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
                                   ,
                                   1000, 2000, 5000, 10000
#endif
};
static const INT online_2d[] = {50, 100, 200
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
                                   ,
                                   500, 1000
#endif
};
static const INT online_3d[] = {50};
static const INT online_4d[] = {28};

void Y(check_nfft_ng_1d_online)(void) {
  online_uniform("nfft_ng_1d_online", 1, online_1d,
                    (int)(SIZE(online_1d)), 50, NG_TRAFO);
}
void Y(check_nfft_ng_adjoint_1d_online)(void) {
  online_uniform("nfft_ng_adjoint_1d_online", 1, online_1d,
                    (int)(SIZE(online_1d)), 50, NG_ADJOINT);
}
void Y(check_nfft_ng_2d_online)(void) {
  online_uniform("nfft_ng_2d_online", 2, online_2d,
                    (int)(SIZE(online_2d)), 50, NG_TRAFO);
}
void Y(check_nfft_ng_adjoint_2d_online)(void) {
  online_uniform("nfft_ng_adjoint_2d_online", 2, online_2d,
                    (int)(SIZE(online_2d)), 50, NG_ADJOINT);
}
void Y(check_nfft_ng_3d_online)(void) {
  online_uniform("nfft_ng_3d_online", 3, online_3d,
                    (int)(SIZE(online_3d)), 50, NG_TRAFO);
}
void Y(check_nfft_ng_adjoint_3d_online)(void) {
  online_uniform("nfft_ng_adjoint_3d_online", 3, online_3d,
                    (int)(SIZE(online_3d)), 50, NG_ADJOINT);
}
void Y(check_nfft_ng_4d_online)(void) {
  online_uniform("nfft_ng_4d_online", 4, online_4d,
                    (int)(SIZE(online_4d)), 50, NG_TRAFO);
}
void Y(check_nfft_ng_adjoint_4d_online)(void) {
  online_uniform("nfft_ng_adjoint_4d_online", 4, online_4d,
                    (int)(SIZE(online_4d)), 50, NG_ADJOINT);
}

typedef struct {
  const char *label;
  int d;
  INT N[NG_MAX_D];
  const int *variant;
  INT M;
} geom;

#define NG_I NFFT_NDFT_TYPE_I
#define NG_II NFFT_NDFT_TYPE_II

static const int v_ii1[1] = {NG_II};
static const int v_ii2[2] = {NG_II, NG_II};
static const int v_ii3[3] = {NG_II, NG_II, NG_II};
static const int v_ii4[4] = {NG_II, NG_II, NG_II, NG_II};
static const int v_ii_i[2] = {NG_II, NG_I};
static const int v_i_i_ii[3] = {NG_I, NG_I, NG_II};
static const int v_i_i_ii_i[4] = {NG_I, NG_I, NG_II, NG_I};
static const int v_unit_ii[4] = {NG_I, NG_II, NG_I, NG_II};

/* Sizes are on the legacy online scale and clear the fast guard for every
 * WINDOW_HELP_ESTIMATE_m in the window matrix (17 is the largest): every live
 * axis has N >= 20, so n = 2*next_power_of_2(N) >= 64 > 2m + 2. */
static const geom variants_1d[] = {
    {"1d even type-I", 1, {50}, 0, 50},
    {"1d even type-II", 1, {50}, v_ii1, 50},
    {"1d odd", 1, {51}, 0, 50},
    {"1d even type-II (100)", 1, {100}, v_ii1, 50},
    {"1d odd (101)", 1, {101}, 0, 50},
};

static const geom variants_2d[] = {
    {"2d even type-I", 2, {50, 50}, 0, 50},
    {"2d even type-II", 2, {50, 50}, v_ii2, 50},
    {"2d odd", 2, {51, 51}, 0, 50},
    {"2d type-II + odd", 2, {50, 51}, v_ii_i, 50},
};

static const geom variants_3d[] = {
    {"3d even type-I", 3, {20, 20, 20}, 0, 50},
    {"3d even type-II", 3, {20, 20, 20}, v_ii3, 50},
    {"3d odd", 3, {21, 21, 21}, 0, 50},
    {"3d type-I + odd + type-II", 3, {20, 21, 22}, v_i_i_ii, 50},
};

static const geom variants_4d[] = {
    {"4d even type-I", 4, {20, 20, 20, 20}, 0, 50},
    {"4d even type-II", 4, {20, 20, 20, 20}, v_ii4, 50},
    {"4d odd", 4, {21, 21, 21, 21}, 0, 50},
    {"4d type-I + odd + type-II", 4, {20, 21, 22, 23}, v_i_i_ii_i, 50},
};

static const geom units_1d[] = {
    {"1d all-unit (rank 0)", 1, {1}, 0, 50},
};

static const geom units_2d[] = {
    {"2d unit axis 0", 2, {1, 50}, 0, 50},
    {"2d unit axis 1", 2, {50, 1}, 0, 50},
    {"2d all-unit (rank 0)", 2, {1, 1}, 0, 50},
};

static const geom units_3d[] = {
    {"3d unit axis 0", 3, {1, 20, 20}, 0, 50},
    {"3d unit axis 1", 3, {20, 1, 20}, 0, 50},
    {"3d unit axis 2", 3, {20, 20, 1}, 0, 50},
    {"3d two unit axes", 3, {1, 1, 20}, 0, 50},
    {"3d all-unit (rank 0)", 3, {1, 1, 1}, 0, 50},
};

static const geom units_4d[] = {
    {"4d unit axis 0", 4, {1, 20, 20, 20}, 0, 50},
    {"4d unit axis 3", 4, {20, 20, 20, 1}, 0, 50},
    {"4d two unit axes", 4, {1, 20, 1, 20}, 0, 50},
    {"4d two unit axes, type-II live", 4, {1, 20, 1, 20}, v_unit_ii, 50},
    {"4d three unit axes", 4, {1, 1, 20, 1}, 0, 50},
    {"4d all-unit (rank 0)", 4, {1, 1, 1, 1}, 0, 50},
};

/* Rank 0 is exact by construction -- forward broadcasts the single
 * coefficient, adjoint sums the nodes -- so assert the values, not a relative
 * error against an oracle running the same solver. */
static int check_rank0(const geom *g, int kind) {
  const INT M = g->M;
  R *x = (R *)Y(malloc)((size_t)g->d * (size_t)M * sizeof(R));
  C *in = (C *)Y(malloc)((size_t)(kind == NG_TRAFO ? 1 : M) * sizeof(C));
  C *out = (C *)Y(malloc)((size_t)(kind == NG_TRAFO ? M : 1) * sizeof(C));
  R tol = K(56.0) * Y(float_property)(NFFT_EPSILON) * (R)M;
  int ok = 1;
  INT j;

  draw(g->d, M, kind == NG_TRAFO ? 1 : M, x, in);

  if (!run_plan(g->d, g->N, g->variant, M, 1, mode_flags(NG_FAST), x, in, out,
               kind)) {
    CU_FAIL("the rank-0 solver must serve an all-unit problem");
    ok = 0;
  } else if (kind == NG_TRAFO) {
    for (j = 0; j < M; j++)
      if (CABS(out[j] - in[0]) > tol)
        ok = 0;
  } else {
    C acc = K(0.0);
    for (j = 0; j < M; j++)
      acc += in[j];
    if (CABS(out[0] - acc) > tol * (CABS(acc) + K(1.0)))
      ok = 0;
  }

  printf("%-34s d = %d, rank 0, %-9s -> %s\n", g->label, g->d,
         kind == NG_TRAFO ? "trafo" : "adjoint", ok ? "OK" : "FAIL");
  CU_ASSERT(ok);

  Y(free)(x);
  Y(free)(in);
  Y(free)(out);
  return ok;
}

static void run_geoms(const char *title, const geom *g, int count) {
  int i, ok = 1, ran = 0;
  printf("%s:\n", title);
  for (i = 0; i < count; i++) {
    int kind;
    for (kind = NG_TRAFO; kind <= NG_ADJOINT; kind++) {
      int r;
      if (prod(g[i].d, g[i].N) == (INT)1) {
        r = check_rank0(&g[i], kind);
      } else {
        r = online_case(g[i].label, g[i].d, g[i].N, g[i].variant, g[i].M,
                           kind);
      }
      if (r >= 0) {
        ok = MIN(ok, r);
        ran++;
      }
    }
  }
  CU_ASSERT(ran > 0);
  CU_ASSERT(ok);
  Y(the_planner_destroy)();
}

void Y(check_nfft_ng_fast_variants_1d)(void) {
  run_geoms("nfft_ng_fast_variants_1d", variants_1d,
               (int)(SIZE(variants_1d)));
}
void Y(check_nfft_ng_fast_variants_2d)(void) {
  run_geoms("nfft_ng_fast_variants_2d", variants_2d,
               (int)(SIZE(variants_2d)));
}
void Y(check_nfft_ng_fast_variants_3d)(void) {
  run_geoms("nfft_ng_fast_variants_3d", variants_3d,
               (int)(SIZE(variants_3d)));
}
void Y(check_nfft_ng_fast_variants_4d)(void) {
  run_geoms("nfft_ng_fast_variants_4d", variants_4d,
               (int)(SIZE(variants_4d)));
}

void Y(check_nfft_ng_fast_unit_axes_1d)(void) {
  run_geoms("nfft_ng_fast_unit_axes_1d", units_1d,
               (int)(SIZE(units_1d)));
}
void Y(check_nfft_ng_fast_unit_axes_2d)(void) {
  run_geoms("nfft_ng_fast_unit_axes_2d", units_2d,
               (int)(SIZE(units_2d)));
}
void Y(check_nfft_ng_fast_unit_axes_3d)(void) {
  run_geoms("nfft_ng_fast_unit_axes_3d", units_3d,
               (int)(SIZE(units_3d)));
}
void Y(check_nfft_ng_fast_unit_axes_4d)(void) {
  run_geoms("nfft_ng_fast_unit_axes_4d", units_4d,
               (int)(SIZE(units_4d)));
}
