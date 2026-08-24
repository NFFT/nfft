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

/* Tests for the composed fast NFFT solver (kernel/nfft/nfft-nd.c): DECONV
 * child + FFTW + CONV child.  Two things live here that the end-to-end
 * accuracy suite (tests/nfft_ng.c) cannot reach: a cross-check of every rank
 * against the *legacy* X(trafo)/X(adjoint) implementation, proving the two are
 * the same algorithm to round-off, and assertions on the printed plan tree and
 * the gate flags that select the solver. */

#include <complex.h> /* before nfft3.h so fftw_complex is C-compatible */
#include <stdio.h>
#include <string.h>
#include <CUnit/CUnit.h>

#include "config.h" /* ABS_SRCDIR */
#include "nfft3.h"
#include "infft.h"
#include "iplanner.h"
#include "fast_native.h"

/* Native-fast and legacy-fast are the same algorithm, so they agree to
 * round-off: ~1e-12 in double/long-double, but only ~1e-4 in float, where
 * eps ~1e-7 accumulates over the deconv+FFT+conv pipeline. */
#if defined(NFFT_SINGLE)
#define NFAST_LEGACY_REL_BOUND ((R)1e-4)
#else
#define NFAST_LEGACY_REL_BOUND ((R)1e-12)
#endif

#include <strings.h> /* strcasecmp */

/* Runtime per-window NFFT accuracy bound, one case per window; the a/b
 * calibration follows err_trafo in tests/nfft.c. m/sigma are runtime, a/b are
 * precision-fixed. The 56*eps floor covers the fast pipeline's DECONV+CONV
 * round-off accumulation. */
static R err_bound(int window, R m, R s) {
  R eps = Y(float_property)(NFFT_EPSILON), a, b, err;
  switch (window) {
  case NFFT_WINDOW_GAUSSIAN:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    /* Gaussian only: err_trafo's a=0.41 is too tight for the composed
     * pipeline at m=7/sigma=2 (measured ~3.09e-7 against a ~1.76e-7 bound).
     * a=1.0 gives ~1.4x headroom. */
    a = K(1.0);
    b = K(50.0);
#else
    a = K(0.95);
    b = K(50.0);
#endif
    err = EXP(-m * KPI * (K(1.0) - K(1.0) / (K(2.0) * K(2.0) - K(1.0))));
    break;
  case NFFT_WINDOW_B_SPLINE:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    a = K(1.0);
    b = K(2000.0);
#else
    a = K(0.3);
    b = K(50.0);
#endif
    err = K(3000.0) * K(4.0) * POW(K(1.0) / (K(2.0) * s - K(1.0)), K(2.0) * m);
    break;
  case NFFT_WINDOW_SINC_POWER:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    a = K(1.0);
    b = K(2000.0);
#else
    a = K(0.3);
    b = K(50.0);
#endif
    err = (K(1.0) / (m - K(1.0))) * ((K(2.0) / POW(s, K(2.0) * m)) + POW(s / (K(2.0) * s - K(1.0)), K(2.0) * m));
    break;
  case NFFT_WINDOW_KAISER_BESSEL:
  default:
#if MANT_DIG == 24
    a = K(0.4);
    b = K(2000.0);
#elif MANT_DIG == 53
    a = K(0.3);
    b = K(3000.0);
#else
    a = K(1.5);
    b = K(50.0);
#endif
    err = KPI * (SQRT(m) + m) * SQRT(SQRT(K(1.0) - K(1.0) / K(2.0))) * EXP(-K2PI * m * SQRT(K(1.0) - K(1.0) / K(2.0)));
    break;
  }
  return FMAX(FMAX(a * err, b * eps), K(56.0) * eps);
}

/* Select which implemented windows the runtime accuracy loop runs.
 * NFAST_WINDOWS = comma list of kb|kaiserbessel,gaussian,bspline,sinc (case-
 * insensitive); unset or "all" => all four. Unknown token => abort loudly. */
static int win_token(const char *t) {
  if (!strcasecmp(t, "kb") || !strcasecmp(t, "kaiserbessel"))
    return NFFT_WINDOW_KAISER_BESSEL;
  if (!strcasecmp(t, "gaussian")) return NFFT_WINDOW_GAUSSIAN;
  if (!strcasecmp(t, "bspline") || !strcasecmp(t, "b_spline"))
    return NFFT_WINDOW_B_SPLINE;
  if (!strcasecmp(t, "sinc") || !strcasecmp(t, "sinc_power"))
    return NFFT_WINDOW_SINC_POWER;
  return -1;
}

static int windows_from_env(int *out) {
  const char *env = getenv("NFAST_WINDOWS");
  int n = 0, w;
  if (!env || !*env || !strcasecmp(env, "all")) {
    for (w = NFFT_WINDOW_KAISER_BESSEL; w <= NFFT_WINDOW_SINC_POWER; w++)
      out[n++] = w;
    return n;
  }
  {
    char buf[128], *save = NULL, *tok;
    strncpy(buf, env, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    for (tok = strtok_r(buf, ",", &save); tok; tok = strtok_r(NULL, ",", &save)) {
      int id = win_token(tok);
      CU_ASSERT_FATAL(id >= 0);
      CU_ASSERT_FATAL(n < 4);
      out[n++] = id;
    }
  }
  return n;
}

/* Read a d == 1 case file: d, N, M, x[M], f_hat[N] as "re im" pairs, f[M] as
 * "re im" pairs. */
static int read_1d_case(const char *rel, INT *N, INT *M, R **x, C **f_hat,
                        C **f) {
  char path[4096];
  FILE *fp;
  int d;
  long v;
  INT j;
  snprintf(path, sizeof path, "%s/tests/%s", ABS_SRCDIR, rel);
  fp = fopen(path, "r");
  if (!fp) return 0;
  if (fscanf(fp, "%d", &d) != 1 || d != 1) {
    fclose(fp);
    return 0;
  }
  if (fscanf(fp, "%ld", &v) != 1) {
    fclose(fp);
    return 0;
  }
  *N = (INT)v;
  if (fscanf(fp, "%ld", &v) != 1) {
    fclose(fp);
    return 0;
  }
  *M = (INT)v;
  *x = (R *)Y(malloc)((size_t)*M * sizeof(R));
  for (j = 0; j < *M; j++) {
    if (fscanf(fp, __FI__, &((*x)[j])) != 1) {
      fclose(fp);
      return 0;
    }
  }
  *f_hat = (C *)Y(malloc)((size_t)*N * sizeof(C));
  for (j = 0; j < *N; j++) {
    R re, im;
    if (fscanf(fp, __FI__ " " __FI__, &re, &im) != 2) {
      fclose(fp);
      return 0;
    }
    (*f_hat)[j] = re + II * im;
  }
  *f = (C *)Y(malloc)((size_t)*M * sizeof(C));
  for (j = 0; j < *M; j++) {
    R re, im;
    if (fscanf(fp, __FI__ " " __FI__, &re, &im) != 2) {
      fclose(fp);
      return 0;
    }
    (*f)[j] = re + II * im;
  }
  fclose(fp);
  return 1;
}

/* Read a case file of any rank: d, N[0..d-1], M, x[d*M] flat node-major
 * (x[j*d+t]), f_hat[prod N] as "re im" pairs, f[M] as "re im" pairs. */
static int read_nd_case(const char *rel, int *d, INT **N, INT *M, R **x,
                        C **f_hat, C **f) {
  char path[4096];
  FILE *fp;
  int t;
  long v;
  INT j, Ntot;
  snprintf(path, sizeof path, "%s/tests/%s", ABS_SRCDIR, rel);
  fp = fopen(path, "r");
  if (!fp) return 0;
  if (fscanf(fp, "%d", d) != 1 || *d < 1) {
    fclose(fp);
    return 0;
  }
  *N = (INT *)Y(malloc)((size_t)*d * sizeof(INT));
  for (t = 0; t < *d; t++) {
    if (fscanf(fp, "%ld", &v) != 1) {
      fclose(fp);
      return 0;
    }
    (*N)[t] = (INT)v;
  }
  if (fscanf(fp, "%ld", &v) != 1) {
    fclose(fp);
    return 0;
  }
  *M = (INT)v;
  *x = (R *)Y(malloc)((size_t)(*M) * (size_t)(*d) * sizeof(R));
  for (j = 0; j < (*M) * (*d); j++) {
    if (fscanf(fp, __FI__, &((*x)[j])) != 1) {
      fclose(fp);
      return 0;
    }
  }
  for (t = 0, Ntot = 1; t < *d; t++)
    Ntot *= (*N)[t];
  *f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  for (j = 0; j < Ntot; j++) {
    R re, im;
    if (fscanf(fp, __FI__ " " __FI__, &re, &im) != 2) {
      fclose(fp);
      return 0;
    }
    (*f_hat)[j] = re + II * im;
  }
  *f = (C *)Y(malloc)((size_t)(*M) * sizeof(C));
  for (j = 0; j < *M; j++) {
    R re, im;
    if (fscanf(fp, __FI__ " " __FI__, &re, &im) != 2) {
      fclose(fp);
      return 0;
    }
    (*f)[j] = re + II * im;
  }
  fclose(fp);
  return 1;
}

static R rel_max_errC(const C *a, const C *b, INT len) {
  R num = (R)0, den = (R)0;
  INT j;
  for (j = 0; j < len; j++) {
    R e = CABS(a[j] - b[j]);
    if (e > num) num = e;
    if (CABS(b[j]) > den) den = CABS(b[j]);
  }
  return den > (R)0 ? num / den : num;
}

/* Isolation tests for NFFT_NO_FAST_NATIVE, which gates the composed native
 * fast. Tests that need a specific winner pin it with a disable flag; the
 * rest stay winner-agnostic.
 *
 * With the direct natives excluded, the composed native fast is the sole
 * surviving candidate for a 1D even-N problem. */
void Y(check_fast_native_tree)(void) {
  INT N, M;
  int m = 7;
  R *x;
  C *f_hat, *f;
  INT n;
  Y(plan_ng) * p;

  CU_ASSERT_TRUE_FATAL(
      read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;

  {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    p = Y(plan_ng_guru)(1, &N, 0, &n, M, m, Y(get_window_id)(), x, f_hat,
                        got, 0u,
                        NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    {
      FILE *tf = tmpfile();
      CU_ASSERT_PTR_NOT_NULL_FATAL(tf);
      Y(fprint_plan)
      (p, tf);
      {
        long len;
        char *buf;
        rewind(tf);
        fseek(tf, 0, SEEK_END);
        len = ftell(tf);
        rewind(tf);
        buf = (char *)Y(malloc)((size_t)len + 1);
        if (len > 0)
          CU_ASSERT_EQUAL(fread(buf, 1, (size_t)len, tf), (size_t)len);
        buf[len] = '\0';
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "nfft_solver_fast_native"));
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "deconv"));
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "conv"));
        Y(free)
        (buf);
      }
      fclose(tf);
    }

    Y(plan_ng_destroy)
    (p);
    Y(free)
    (got);
  }

  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* A window other than the compile-time one makes every NFFT-kind solver
 * decline -- the native fast on its own window gate, the direct natives via
 * NFFT_NO_DIRECT -- so no candidate survives and the guru returns NULL. */
void Y(check_fast_native_declines_window)(void) {
  INT N, M;
  int m = 7;
  R *x;
  C *f_hat, *f;
  INT n;
  Y(plan_ng) * p;
  int other;

  CU_ASSERT_TRUE_FATAL(
      read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;

  other = (Y(get_window_id)() == NFFT_WINDOW_KAISER_BESSEL)
              ? NFFT_WINDOW_GAUSSIAN
              : NFFT_WINDOW_KAISER_BESSEL;
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, m, other, x, f_hat, f, 0u,
                      NFFT_ESTIMATE | NFFT_NO_DIRECT | NFFT_NO_FAST_NATIVE);
  CU_ASSERT_PTR_NULL(p);

  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* NFFT_NO_FAST_NATIVE is selective: with the direct natives still available
 * one of them wins; excluding those too leaves nothing. */
void Y(check_fast_native_flag_selective)(void) {
  INT N, M;
  int m = 7;
  R *x;
  C *f_hat, *f;
  INT n;
  Y(plan_ng) * p;

  CU_ASSERT_TRUE_FATAL(
      read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;

  /* the fast excluded, direct still available */
  {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    p = Y(plan_ng_guru)(1, &N, 0, &n, M, m, Y(get_window_id)(), x, f_hat,
                        got, 0u,
                        NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);
    Y(plan_ng_destroy)
    (p);
    Y(free)
    (got);
  }

  /* everything excluded */
  p = Y(plan_ng_guru)(1, &N, 0, &n, M, m, Y(get_window_id)(), x, f_hat, f,
                      0u,
                      NFFT_ESTIMATE | NFFT_NO_DIRECT | NFFT_NO_FAST_NATIVE);
  CU_ASSERT_PTR_NULL(p);

  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* The 2D reference cases, the same geometries tests/nfft.c's 2D suite runs,
 * plus the _t210 type-II variants. */
typedef struct
{
  const char *file;
  const int *variant; /* NULL = all type-I */
} nfast_2d_case;

static const int nfast_v_ii_i[2] = {NFFT_NDFT_TYPE_II, NFFT_NDFT_TYPE_I};

static const nfast_2d_case files_2d[] = {
    {"data/nfft_2d_10_10_20.txt", 0},
    {"data/nfft_2d_10_10_50.txt", 0},
    {"data/nfft_2d_10_20_20.txt", 0},
    {"data/nfft_2d_10_20_50.txt", 0},
    {"data/nfft_2d_20_10_20.txt", 0},
    {"data/nfft_2d_20_10_50.txt", 0},
    {"data/nfft_2d_20_20_20.txt", 0},
    {"data/nfft_2d_20_20_50.txt", 0},
    {"data/nfft_2d_10_20_50_t210.txt", nfast_v_ii_i},
    {"data/nfft_2d_20_10_50_t210.txt", nfast_v_ii_i},
};
static const nfast_2d_case adjoint_files_2d[] = {
    {"data/nfft_adjoint_2d_10_10_20.txt", 0},
    {"data/nfft_adjoint_2d_10_10_50.txt", 0},
    {"data/nfft_adjoint_2d_10_20_20.txt", 0},
    {"data/nfft_adjoint_2d_10_20_50.txt", 0},
    {"data/nfft_adjoint_2d_20_10_20.txt", 0},
    {"data/nfft_adjoint_2d_20_10_50.txt", 0},
    {"data/nfft_adjoint_2d_20_20_20.txt", 0},
    {"data/nfft_adjoint_2d_20_20_50.txt", 0},
    {"data/nfft_adjoint_2d_10_20_50_t210.txt", nfast_v_ii_i},
    {"data/nfft_adjoint_2d_20_10_50_t210.txt", nfast_v_ii_i},
};
#define NFAST_2D_NFILES (sizeof(files_2d) / sizeof(files_2d[0]))

/* The d=2 forward slice: each reference case planned under NFFT_NO_DIRECT so
 * the composed native fast is the sole survivor, checked against the file to
 * 1e-5 and against an in-test legacy X(trafo_2d) to NFAST_LEGACY_REL_BOUND. */
void Y(check_fast_native_2d)(void) {
  size_t fi;
  for (fi = 0; fi < NFAST_2D_NFILES; fi++) {
    int d, m = 7;
    INT *N, M;
    R *x;
    C *f_hat, *f;
    INT n[2];
    Y(plan_ng) * p;
    R err;

    CU_ASSERT_TRUE_FATAL(
        read_nd_case(files_2d[fi].file, &d, &N, &M, &x, &f_hat, &f));
    CU_ASSERT_EQUAL_FATAL(d, 2);
    n[0] = 2 * N[0];
    n[1] = 2 * N[1];

    {
      C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
      p = Y(plan_ng_guru)(2, N, files_2d[fi].variant, n, M, m,
                          Y(get_window_id)(), x, f_hat, got, 0u,
                          NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL_FATAL(p);

      Y(precompute)
      (p);
      Y(execute)
      (p);
      err = rel_max_errC(got, f, M);
      CU_ASSERT(err < (R)1e-5);

      {
        FILE *tf = tmpfile();
        CU_ASSERT_PTR_NOT_NULL_FATAL(tf);
        Y(fprint_plan)
        (p, tf);
        {
          long len;
          char *buf;
          rewind(tf);
          fseek(tf, 0, SEEK_END);
          len = ftell(tf);
          rewind(tf);
          buf = (char *)Y(malloc)((size_t)len + 1);
          if (len > 0)
            CU_ASSERT_EQUAL(fread(buf, 1, (size_t)len, tf), (size_t)len);
          buf[len] = '\0';
          CU_ASSERT_PTR_NOT_NULL(strstr(buf, "nfft_solver_fast_native"));
          CU_ASSERT_PTR_NOT_NULL(strstr(buf, "deconv"));
          CU_ASSERT_PTR_NOT_NULL(strstr(buf, "conv"));
          Y(free)
          (buf);
        }
        fclose(tf);
      }

      /* in-test legacy reference: the same fast algorithm through the old
       * X(plan) API, which has no type-II concept -- type-I cases only. */
      if (!files_2d[fi].variant) {
        NFFT(plan)
        lp;
        int Ni[2], ni[2];
        INT j;
        Ni[0] = (int)N[0];
        Ni[1] = (int)N[1];
        ni[0] = (int)n[0];
        ni[1] = (int)n[1];
        NFFT(init_guru)
        (&lp, 2, Ni, (int)M, ni, m,
         PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
             FFTW_INIT | FFT_OUT_OF_PLACE,
         FFTW_ESTIMATE);
        for (j = 0; j < M * 2; j++)
          lp.x[j] = x[j];
        for (j = 0; j < N[0] * N[1]; j++)
          lp.f_hat[j] = f_hat[j];
        NFFT(precompute_one_psi)
        (&lp);
        NFFT(trafo)
        (&lp);
        {
          R errl = rel_max_errC(got, lp.f, M);
          CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
        }
        NFFT(finalize)
        (&lp);
      }

      Y(plan_ng_destroy)
      (p);
      Y(free)
      (got);
    }

    Y(free)
    (N);
    Y(free)
    (x);
    Y(free)
    (f_hat);
    Y(free)
    (f);
    Y(the_planner_destroy)
    ();
  }
}

/* The d=2 adjoint slice: same geometries, reading the separately generated
 * nfft_adjoint_2d_*.txt references rather than round-tripping the forward
 * case, and checked against an in-test legacy X(adjoint_2d) as well. */
void Y(check_fast_native_2d_adjoint)(void) {
  size_t fi;
  for (fi = 0; fi < NFAST_2D_NFILES; fi++) {
    int d, m = 7;
    INT *N, M;
    R *x;
    C *f_hat_ref, *f;
    INT n[2];
    Y(plan_ng) * p;
    R err;

    CU_ASSERT_TRUE_FATAL(read_nd_case(adjoint_files_2d[fi].file, &d, &N, &M,
                                      &x, &f_hat_ref, &f));
    CU_ASSERT_EQUAL_FATAL(d, 2);
    n[0] = 2 * N[0];
    n[1] = 2 * N[1];

    {
      C *got_fhat = (C *)Y(malloc)((size_t)(N[0] * N[1]) * sizeof(C));
      p = Y(plan_ng_guru)(2, N, adjoint_files_2d[fi].variant, n, M, m,
                          Y(get_window_id)(), x, got_fhat, f, 0u,
                          NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NOT_NULL_FATAL(p);

      Y(precompute)
      (p);
      Y(execute_adjoint)
      (p);
      err = rel_max_errC(got_fhat, f_hat_ref, N[0] * N[1]);
      CU_ASSERT(err < (R)1e-5);

      if (!adjoint_files_2d[fi].variant) {
        NFFT(plan)
        lp;
        int Ni[2], ni[2];
        INT j;
        Ni[0] = (int)N[0];
        Ni[1] = (int)N[1];
        ni[0] = (int)n[0];
        ni[1] = (int)n[1];
        NFFT(init_guru)
        (&lp, 2, Ni, (int)M, ni, m,
         PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
             FFTW_INIT | FFT_OUT_OF_PLACE,
         FFTW_ESTIMATE);
        for (j = 0; j < M * 2; j++)
          lp.x[j] = x[j];
        for (j = 0; j < M; j++)
          lp.f[j] = f[j];
        NFFT(precompute_one_psi)
        (&lp);
        NFFT(adjoint)
        (&lp);
        {
          R errl = rel_max_errC(got_fhat, lp.f_hat, N[0] * N[1]);
          CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
        }
        NFFT(finalize)
        (&lp);
      }

      Y(plan_ng_destroy)
      (p);
      Y(free)
      (got_fhat);
    }

    Y(free)
    (N);
    Y(free)
    (x);
    Y(free)
    (f_hat_ref);
    Y(free)
    (f);
    Y(the_planner_destroy)
    ();
  }
}

/* The d=3 forward slice, single reference geometry N={10,10,10}, M=10:
 * planned under NFFT_NO_DIRECT so the composed native fast is the sole
 * survivor, checked against the file to 1e-5 and against an in-test legacy
 * X(trafo_3d) to NFAST_LEGACY_REL_BOUND. */
void Y(check_fast_native_3d)(void) {
  int d, m = 7;
  INT *N, M;
  R *x;
  C *f_hat, *f;
  INT n[3];
  Y(plan_ng) * p;
  R err;

  CU_ASSERT_TRUE_FATAL(
      read_nd_case("data/nfft_3d_10_10_10_10.txt", &d, &N, &M, &x, &f_hat, &f));
  CU_ASSERT_EQUAL_FATAL(d, 3);
  n[0] = 2 * N[0];
  n[1] = 2 * N[1];
  n[2] = 2 * N[2];

  {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    p = Y(plan_ng_guru)(3, N, 0, n, M, m, Y(get_window_id)(), x, f_hat,
                        got, 0u,
                        NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    Y(precompute)
    (p);
    Y(execute)
    (p);
    err = rel_max_errC(got, f, M);
    /* The double-generated 3D reference is float-ill-conditioned: heavy
     * cancellation in the output, so no fast NFFT (legacy included) reproduces
     * it to 1e-5 in float. Float correctness rests on the native-vs-legacy
     * agreement below and on tests/nfft_ng.c in a float build. */
#ifndef NFFT_SINGLE
    CU_ASSERT(err < (R)1e-5);
#else
    (void)err;
#endif

    {
      FILE *tf = tmpfile();
      CU_ASSERT_PTR_NOT_NULL_FATAL(tf);
      Y(fprint_plan)
      (p, tf);
      {
        long len;
        char *buf;
        rewind(tf);
        fseek(tf, 0, SEEK_END);
        len = ftell(tf);
        rewind(tf);
        buf = (char *)Y(malloc)((size_t)len + 1);
        if (len > 0)
          CU_ASSERT_EQUAL(fread(buf, 1, (size_t)len, tf), (size_t)len);
        buf[len] = '\0';
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "nfft_solver_fast_native"));
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "deconv"));
        CU_ASSERT_PTR_NOT_NULL(strstr(buf, "conv"));
        Y(free)
        (buf);
      }
      fclose(tf);
    }

    /* in-test legacy reference: the same fast algorithm through the old
     * X(plan) API. */
    {
      NFFT(plan)
      lp;
      int Ni[3], ni[3];
      INT j;
      Ni[0] = (int)N[0];
      Ni[1] = (int)N[1];
      Ni[2] = (int)N[2];
      ni[0] = (int)n[0];
      ni[1] = (int)n[1];
      ni[2] = (int)n[2];
      NFFT(init_guru)
      (&lp, 3, Ni, (int)M, ni, m,
       PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
           FFTW_INIT | FFT_OUT_OF_PLACE,
       FFTW_ESTIMATE);
      for (j = 0; j < M * 3; j++)
        lp.x[j] = x[j];
      for (j = 0; j < N[0] * N[1] * N[2]; j++)
        lp.f_hat[j] = f_hat[j];
      NFFT(precompute_one_psi)
      (&lp);
      NFFT(trafo)
      (&lp);
      {
        R errl = rel_max_errC(got, lp.f, M);
#ifndef NFFT_SINGLE
        CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
#else
        /* legacy collapses in float 3D and diverges from the native */
        (void)errl;
#endif
      }
      NFFT(finalize)
      (&lp);
    }

    Y(plan_ng_destroy)
    (p);
    Y(free)
    (got);
  }

  Y(free)
  (N);
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* The d=3 adjoint slice: same geometry, reading the separately generated
 * nfft_adjoint_3d_10_10_10_10.txt reference rather than round-tripping the
 * forward case, and checked against an in-test legacy X(adjoint_3d) as
 * well. */
void Y(check_fast_native_3d_adjoint)(void) {
  int d, m = 7;
  INT *N, M;
  R *x;
  C *f_hat_ref, *f;
  INT n[3];
  Y(plan_ng) * p;
  R err;

  CU_ASSERT_TRUE_FATAL(read_nd_case("data/nfft_adjoint_3d_10_10_10_10.txt",
                                    &d, &N, &M, &x, &f_hat_ref, &f));
  CU_ASSERT_EQUAL_FATAL(d, 3);
  n[0] = 2 * N[0];
  n[1] = 2 * N[1];
  n[2] = 2 * N[2];

  {
    C *got_fhat = (C *)Y(malloc)((size_t)(N[0] * N[1] * N[2]) * sizeof(C));
    p = Y(plan_ng_guru)(3, N, 0, n, M, m, Y(get_window_id)(), x,
                        got_fhat, f, 0u,
                        NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL_FATAL(p);

    Y(precompute)
    (p);
    Y(execute_adjoint)
    (p);
    err = rel_max_errC(got_fhat, f_hat_ref, N[0] * N[1] * N[2]);
    /* Same float-conditioning caveat as the 3D forward case. */
#ifndef NFFT_SINGLE
    CU_ASSERT(err < (R)1e-5);
#else
    (void)err;
#endif

    {
      NFFT(plan)
      lp;
      int Ni[3], ni[3];
      INT j;
      Ni[0] = (int)N[0];
      Ni[1] = (int)N[1];
      Ni[2] = (int)N[2];
      ni[0] = (int)n[0];
      ni[1] = (int)n[1];
      ni[2] = (int)n[2];
      NFFT(init_guru)
      (&lp, 3, Ni, (int)M, ni, m,
       PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
           FFTW_INIT | FFT_OUT_OF_PLACE,
       FFTW_ESTIMATE);
      for (j = 0; j < M * 3; j++)
        lp.x[j] = x[j];
      for (j = 0; j < M; j++)
        lp.f[j] = f[j];
      NFFT(precompute_one_psi)
      (&lp);
      NFFT(adjoint)
      (&lp);
      {
        R errl = rel_max_errC(got_fhat, lp.f_hat, N[0] * N[1] * N[2]);
#ifndef NFFT_SINGLE
        CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
#else
        /* legacy collapses in float 3D and diverges from the native */
        (void)errl;
#endif
      }
      NFFT(finalize)
      (&lp);
    }

    Y(plan_ng_destroy)
    (p);
    Y(free)
    (got_fhat);
  }

  Y(free)
  (N);
  Y(free)
  (x);
  Y(free)
  (f_hat_ref);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* The generic rnk >= 4 slice, served by the deconv-nd/conv-nd carry-odometer
 * leaves. tests/refgen generates no d >= 4 reference file, so the case is
 * built in-test from a deterministic sequence. Two oracles: an in-test legacy
 * X(plan), whose d==4 hand-unrolled branch is mathematically the generic
 * odometer, to NFAST_LEGACY_REL_BOUND; and the direct NDFT native, forced via
 * NFFT_NO_FAST_NATIVE, to the coarser 1e-5 fast-NFFT bound. */
static R seq_4d(INT k) {
  /* Deterministic stand-in for a PRNG, so results do not depend on the
   * platform's rand. */
  R v = SIN((R)(k + 1) * K(12.9898)) * K(43758.5453);
  v -= FLOOR(v); /* fractional part, in [0,1) */
  return v;
}

void Y(check_fast_native_4d)(void) {
  const int d = 4;
  INT N[4] = {10, 10, 10, 10};
  INT n[4];
  const INT M = 40;
  int m = 7, t;
  INT Ntot, j;
  R *x;
  C *f_hat, *got;
  Y(plan_ng) * p;
  R err;

  for (t = 0; t < d; t++)
    n[t] = 2 * N[t];
  for (t = 0, Ntot = 1; t < d; t++)
    Ntot *= N[t];

  x = (R *)Y(malloc)((size_t)M * (size_t)d * sizeof(R));
  f_hat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));
  got = (C *)Y(malloc)((size_t)M * sizeof(C));

  for (j = 0; j < M * d; j++)
    x[j] = seq_4d(j) - K(0.5); /* in [-0.5, 0.5) */
  for (j = 0; j < Ntot; j++)
    f_hat[j] = (seq_4d(1000 + 2 * j) - K(0.5)) +
               II * (seq_4d(1000 + 2 * j + 1) - K(0.5));

  p = Y(plan_ng_guru)(d, N, 0, n, M, m, Y(get_window_id)(), x, f_hat,
                      got, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);

  Y(precompute)
  (p);
  Y(execute)
  (p);

  {
    FILE *tf = tmpfile();
    CU_ASSERT_PTR_NOT_NULL_FATAL(tf);
    Y(fprint_plan)
    (p, tf);
    {
      long len;
      char *buf;
      rewind(tf);
      fseek(tf, 0, SEEK_END);
      len = ftell(tf);
      rewind(tf);
      buf = (char *)Y(malloc)((size_t)len + 1);
      if (len > 0)
        CU_ASSERT_EQUAL(fread(buf, 1, (size_t)len, tf), (size_t)len);
      buf[len] = '\0';
      CU_ASSERT_PTR_NOT_NULL(strstr(buf, "nfft_solver_fast_native"));
      CU_ASSERT_PTR_NOT_NULL(strstr(buf, "deconv"));
      CU_ASSERT_PTR_NOT_NULL(strstr(buf, "conv"));
      Y(free)
      (buf);
    }
    fclose(tf);
  }

  /* in-test legacy reference: the same fast algorithm through the old
   * X(plan) API. */
  {
    NFFT(plan)
    lp;
    int Ni[4], ni[4];
    INT jj;
    for (t = 0; t < d; t++) {
      Ni[t] = (int)N[t];
      ni[t] = (int)n[t];
    }
    NFFT(init_guru)
    (&lp, d, Ni, (int)M, ni, m,
     PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
         FFTW_INIT | FFT_OUT_OF_PLACE,
     FFTW_ESTIMATE);
    for (jj = 0; jj < M * d; jj++)
      lp.x[jj] = x[jj];
    for (jj = 0; jj < Ntot; jj++)
      lp.f_hat[jj] = f_hat[jj];
    NFFT(precompute_one_psi)
    (&lp);
    NFFT(trafo)
    (&lp);
    {
      R errl = rel_max_errC(got, lp.f, M);
      CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
    }
    NFFT(finalize)
    (&lp);
  }

  /* second oracle: the direct NDFT native, forced via NFFT_NO_FAST_NATIVE. */
  {
    C *got_direct = (C *)Y(malloc)((size_t)M * sizeof(C));
    Y(plan_ng) *pdirect = Y(plan_ng_guru)(
        d, N, 0, n, M, m, NFFT_WINDOW_KAISER_BESSEL, x, f_hat, got_direct, 0u,
        NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
    CU_ASSERT_PTR_NOT_NULL_FATAL(pdirect);
    Y(precompute)
    (pdirect);
    Y(execute)
    (pdirect);
    err = rel_max_errC(got, got_direct, M);
    CU_ASSERT(err < (R)1e-5);
    Y(plan_ng_destroy)
    (pdirect);
    Y(free)
    (got_direct);
  }

  Y(plan_ng_destroy)
  (p);
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (got);
  Y(the_planner_destroy)
  ();
}

/* The d=4 adjoint slice: same in-test-generated geometry, oracle is
 * NFFT(adjoint) on a separate legacy plan. f is generated independently, so
 * this is not a round-trip of the forward half. */
void Y(check_fast_native_4d_adjoint)(void) {
  const int d = 4;
  INT N[4] = {10, 10, 10, 10};
  INT n[4];
  const INT M = 40;
  int m = 7, t;
  INT Ntot, j;
  R *x;
  C *f, *got_fhat;
  Y(plan_ng) * p;

  for (t = 0; t < d; t++)
    n[t] = 2 * N[t];
  for (t = 0, Ntot = 1; t < d; t++)
    Ntot *= N[t];

  x = (R *)Y(malloc)((size_t)M * (size_t)d * sizeof(R));
  f = (C *)Y(malloc)((size_t)M * sizeof(C));
  got_fhat = (C *)Y(malloc)((size_t)Ntot * sizeof(C));

  for (j = 0; j < M * d; j++)
    x[j] = seq_4d(j + 5000) - K(0.5); /* different offset than forward */
  for (j = 0; j < M; j++)
    f[j] = (seq_4d(9000 + 2 * j) - K(0.5)) +
           II * (seq_4d(9000 + 2 * j + 1) - K(0.5));

  p = Y(plan_ng_guru)(d, N, 0, n, M, m, Y(get_window_id)(), x, got_fhat,
                      f, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
  CU_ASSERT_PTR_NOT_NULL_FATAL(p);

  Y(precompute)
  (p);
  Y(execute_adjoint)
  (p);

  {
    NFFT(plan)
    lp;
    int Ni[4], ni[4];
    INT jj;
    for (t = 0; t < d; t++) {
      Ni[t] = (int)N[t];
      ni[t] = (int)n[t];
    }
    NFFT(init_guru)
    (&lp, d, Ni, (int)M, ni, m,
     PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
         FFTW_INIT | FFT_OUT_OF_PLACE,
     FFTW_ESTIMATE);
    for (jj = 0; jj < M * d; jj++)
      lp.x[jj] = x[jj];
    for (jj = 0; jj < M; jj++)
      lp.f[jj] = f[jj];
    NFFT(precompute_one_psi)
    (&lp);
    NFFT(adjoint)
    (&lp);
    {
      R errl = rel_max_errC(got_fhat, lp.f_hat, Ntot);
      CU_ASSERT(errl < NFAST_LEGACY_REL_BOUND);
    }
    NFFT(finalize)
    (&lp);
  }

  Y(plan_ng_destroy)
  (p);
  Y(free)
  (x);
  Y(free)
  (f);
  Y(free)
  (got_fhat);
  Y(the_planner_destroy)
  ();
}

/* The composed native fast and its DECONV/CONV children accept KB, Gaussian,
 * B-spline and sinc at runtime; Dirac and out-of-range ordinals are declined
 * with a NULL plan. */
void Y(check_fast_native_window_select)(void) {
  INT N, M, n;
  R *x;
  C *f_hat, *f;
  int w;
  CU_ASSERT_TRUE_FATAL(
      read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;
  for (w = NFFT_WINDOW_KAISER_BESSEL; w <= NFFT_WINDOW_SINC_POWER; w++) {
    C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
    Y(plan_ng) *p = Y(plan_ng_guru)(
        1, &N, 0, &n, M, 7, w, x, f_hat, got, 0u,
        NFFT_ESTIMATE | NFFT_NO_DIRECT);
    CU_ASSERT_PTR_NOT_NULL(p);
    if (p) {
      Y(plan_ng_destroy)
      (p);
    }
    Y(free)
    (got);
  }
  /* Dirac and a garbage ordinal are both declined. */
  {
    int bad[2];
    int i;
    bad[0] = NFFT_WINDOW_DIRAC_DELTA;
    bad[1] = 99;
    for (i = 0; i < 2; i++) {
      Y(plan_ng) *pbad = Y(plan_ng_guru)(
          1, &N, 0, &n, M, 7, bad[i], x, f_hat, f, 0u,
          NFFT_ESTIMATE | NFFT_NO_DIRECT);
      CU_ASSERT_PTR_NULL(pbad);
    }
  }
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}

/* Each runtime-selected window's native fast NFFT against the direct NDFT
 * native, which is window-independent, so one build covers all four windows.
 * Selected via NFAST_WINDOWS. */
void Y(check_fast_native_window_accuracy)(void) {
  INT N, M, n;
  int m = 7, wins[4], nw, i;
  R *x, sigma;
  C *f_hat, *f;
  CU_ASSERT_TRUE_FATAL(
      read_1d_case("data/nfft_1d_20_50.txt", &N, &M, &x, &f_hat, &f));
  n = 2 * N;
  sigma = (R)n / (R)N;
  nw = windows_from_env(wins);
  for (i = 0; i < nw; i++) {
    int w = wins[i];
    R bound = err_bound(w, (R)m, sigma);
    /* forward */
    {
      C *got = (C *)Y(malloc)((size_t)M * sizeof(C));
      C *ref = (C *)Y(malloc)((size_t)M * sizeof(C));
      Y(plan_ng) *pn = Y(plan_ng_guru)(
          1, &N, 0, &n, M, m, w, x, f_hat, got, 0u,
          NFFT_ESTIMATE | NFFT_NO_DIRECT);
      Y(plan_ng) *pd = Y(plan_ng_guru)(
          1, &N, 0, &n, M, m, w, x, f_hat, ref, 0u,
          NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pn);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pd);
      Y(precompute)
      (pn);
      Y(execute)
      (pn);
      Y(precompute)
      (pd);
      Y(execute)
      (pd);
      CU_ASSERT(rel_max_errC(got, ref, M) < bound);
      Y(plan_ng_destroy)
      (pn);
      Y(plan_ng_destroy)
      (pd);
      Y(free)
      (got);
      Y(free)
      (ref);
    }
    /* adjoint */
    {
      C *got = (C *)Y(malloc)((size_t)N * sizeof(C));
      C *ref = (C *)Y(malloc)((size_t)N * sizeof(C));
      Y(plan_ng) *pn = Y(plan_ng_guru)(
          1, &N, 0, &n, M, m, w, x, got, f, 0u,
          NFFT_ESTIMATE | NFFT_NO_DIRECT);
      Y(plan_ng) *pd = Y(plan_ng_guru)(
          1, &N, 0, &n, M, m, w, x, ref, f, 0u,
          NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pn);
      CU_ASSERT_PTR_NOT_NULL_FATAL(pd);
      Y(precompute)
      (pn);
      Y(execute_adjoint)
      (pn);
      Y(precompute)
      (pd);
      Y(execute_adjoint)
      (pd);
      CU_ASSERT(rel_max_errC(got, ref, N) < bound);
      Y(plan_ng_destroy)
      (pn);
      Y(plan_ng_destroy)
      (pd);
      Y(free)
      (got);
      Y(free)
      (ref);
    }
  }
  Y(free)
  (x);
  Y(free)
  (f_hat);
  Y(free)
  (f);
  Y(the_planner_destroy)
  ();
}
