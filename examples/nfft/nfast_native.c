/*
 * Copyright (c) 2002, 2026 Jens Keiner, Stefan Kunis, Daniel Potts
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

/*! \file nfast_native.c
 *
 * Six-way check of the composed planner-native fast NFFT: legacy direct NDFT,
 * legacy fast NFFT, planner direct NDFT, and planner native fast NFFT, all
 * built over the same problem (data/nfft_1d_8192_128.txt, an offline
 * exact-precision reference) and the same compile-time-selected window. The two
 * fast NFFTs (legacy + native) are each built twice, once with FFTW_ESTIMATE
 * and once with FFTW_MEASURE for the internal DFT plan, so the effect of the
 * FFTW planner flag is visible side by side. Each of the two planner plans
 * prints its plan tree. For every plan, legacy and planner alike,
 * precompute/forward/adjoint are timed separately; forward-transform error is
 * checked against the file reference and PASS/FAIL'd.
 *
 * The two FFTW_MEASURE plans are constructed after every FFTW_ESTIMATE plan:
 * FFTW caches wisdom process-globally, so a measured plan of a given size would
 * otherwise bless a later estimate plan of the same size and spoil the
 * comparison. The NFFT-level planner mode stays NFFT_ESTIMATE throughout, so
 * only the FFTW DFT plan quality changes between the est/meas rows.
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#if defined(NFFT_SINGLE)
#define NFFT_PRECISION_SINGLE
#elif defined(NFFT_LDOUBLE)
#define NFFT_PRECISION_LONG_DOUBLE
#else
#define NFFT_PRECISION_DOUBLE
#endif

#include "nfft3mp.h"
#include "ticks.h"

#if defined(NFFT_SINGLE)
#define NFAST_CABS cabsf
#elif defined(NFFT_LDOUBLE)
#define NFAST_CABS cabsl
#else
#define NFAST_CABS cabs
#endif

#ifndef NFAST_NATIVE_DATA
#define NFAST_NATIVE_DATA "data/nfft_1d_8192_128.txt"
#endif

/* Loose PASS/FAIL gate; the point of the table is the printed magnitudes
 * (the legacy direct worst), not the threshold. Scaled by build precision: at
 * N=8192 the legacy per-term direct accumulates ~N*eps, ~1e-3 in float but
 * ~1e-12/~1e-15 in double/long-double, so one bound cannot fit all three. */
#if defined(NFFT_SINGLE)
#define NFAST_NATIVE_BOUND NFFT_K(1e-2)
#else
#define NFAST_NATIVE_BOUND NFFT_K(1e-5)
#endif


#define WINDOW_M 12

/* Reads d, N[0], M, x[d*M], f_hat[NN], f[M] (same token layout as
 * tests/nplan_data.c:read_case); d is expected to be 1. */
static int read_reference(const char *path, int *d_out, NFFT_INT *N_out,
                          NFFT_INT *M_out, NFFT_R **x_out, NFFT_C **f_hat_out, NFFT_C **f_out) {
  FILE *fp;
  int d, t;
  long v;
  NFFT_INT N, M, j, nn;
  NFFT_R *x;
  NFFT_C *f_hat, *f;

  fp = fopen(path, "r");
  if (!fp) {
    fprintf(stderr, "nfast_native: could not open reference file '%s'\n", path);
    return 0;
  }

  if (fscanf(fp, "%d", &d) != 1) {
    fclose(fp);
    return 0;
  }

  nn = 1;
  N = 0;
  for (t = 0; t < d; t++) {
    if (fscanf(fp, "%ld", &v) != 1) {
      fclose(fp);
      return 0;
    }
    N = (NFFT_INT)v;
    nn *= N;
  }

  if (fscanf(fp, "%ld", &v) != 1) {
    fclose(fp);
    return 0;
  }
  M = (NFFT_INT)v;

  x = (NFFT_R *)malloc((size_t)(d * M) * sizeof(NFFT_R));
  for (j = 0; j < (NFFT_INT)d * M; j++) {
    double xv;
    if (fscanf(fp, "%lf", &xv) != 1) {
      fclose(fp);
      free(x);
      return 0;
    }
    x[j] = (NFFT_R)xv;
  }

  f_hat = (NFFT_C *)malloc((size_t)nn * sizeof(NFFT_C));
  for (j = 0; j < nn; j++) {
    double re, im;
    if (fscanf(fp, "%lf %lf", &re, &im) != 2) {
      fclose(fp);
      free(x);
      free(f_hat);
      return 0;
    }
    f_hat[j] = (NFFT_R)re + (NFFT_R)im * I;
  }

  f = (NFFT_C *)malloc((size_t)M * sizeof(NFFT_C));
  for (j = 0; j < M; j++) {
    double re, im;
    if (fscanf(fp, "%lf %lf", &re, &im) != 2) {
      fclose(fp);
      free(x);
      free(f_hat);
      free(f);
      return 0;
    }
    f[j] = (NFFT_R)re + (NFFT_R)im * I;
  }

  fclose(fp);

  *d_out = d;
  *N_out = N;
  *M_out = M;
  *x_out = x;
  *f_hat_out = f_hat;
  *f_out = f;
  return 1;
}

/* max|a-b| / max|b| over a length-len complex vector pair. */
static NFFT_R rel_max_err(const NFFT_C *a, const NFFT_C *b, NFFT_INT len) {
  NFFT_R num = NFFT_K(0.0), den = NFFT_K(0.0);
  NFFT_INT j;
  for (j = 0; j < len; j++) {
    NFFT_R d = (NFFT_R)NFAST_CABS(a[j] - b[j]);
    NFFT_R m = (NFFT_R)NFAST_CABS(b[j]);
    if (d > num)
      num = d;
    if (m > den)
      den = m;
  }
  if (den == NFFT_K(0.0))
    return num;
  return num / den;
}

/* Cumulative wall-time budget for each measured step: every precompute /
 * forward / adjoint below is re-run until at least this many seconds have
 * elapsed, and the per-run figures are reduced to a mean +/- standard
 * deviation (one run is too noisy). Plan creation is not measured. */
#define NFAST_MEASURE_SECONDS 2.0

/* Aggregated timing over many runs: mean and population standard deviation of
 * the per-run wall-clock seconds (clock_gettime) and CPU ticks (getticks cycle
 * counter), plus the run count that fit in the budget. */
typedef struct
{
  double secs_mean, secs_std;
  double tks_mean, tks_std;
  long runs;
} nfast_timing;

/* Legacy (NFFT(plan)) and planner (NFFT(plan_ng)) transforms/precompute
 * steps behind one void* signature, so time_run() can measure all of them. */
static void run_legacy_precompute(void *ctx) {
  NFFT(precompute_one_psi)
  ((NFFT(plan) *)ctx);
}

static void run_legacy_direct(void *ctx) {
  NFFT(trafo_direct)
  ((NFFT(plan) *)ctx);
}

static void run_legacy_fast(void *ctx) {
  NFFT(trafo)
  ((NFFT(plan) *)ctx);
}

static void run_legacy_adjoint_direct(void *ctx) {
  NFFT(adjoint_direct)
  ((NFFT(plan) *)ctx);
}

static void run_legacy_adjoint_fast(void *ctx) {
  NFFT(adjoint)
  ((NFFT(plan) *)ctx);
}

static void run_precompute(void *ctx) {
  NFFT(precompute)
  ((NFFT(plan_ng) *)ctx);
}

static void run_plan_ng(void *ctx) {
  NFFT(execute)
  ((NFFT(plan_ng) *)ctx);
}

/* execute_adjoint_on takes explicit f_hat/f buffers rather than the plan's
 * own bound arrays, so a context struct carries them through time_run(). */
typedef struct
{
  NFFT(plan_ng) * p;
  NFFT_C *f_hat;
  NFFT_C *f;
} plan_ng_adjoint_ctx;

static void run_plan_ng_adjoint(void *ctx) {
  plan_ng_adjoint_ctx *c = (plan_ng_adjoint_ctx *)ctx;
  NFFT(execute_adjoint_on)
  (c->p, c->f_hat, c->f);
}

/* Wall-clock seconds as a double, independent of the build precision. The
 * public NFFT(clock_gettime_seconds)() returns the precision real R, which in
 * the float build cannot resolve sub-second intervals (the epoch ~1.7e9 has
 * only ~200 s resolution in float), so read CLOCK_MONOTONIC directly. */
static double wall_seconds(void) {
#if defined(HAVE_CLOCK_GETTIME)
  struct timespec ts;
  if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0)
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#endif
  return 0.0;
}

static nfast_timing time_run(void (*fn)(void *), void *ctx) {
  /* Welford online mean/variance for wall-seconds and CPU ticks. */
  double s_mean = 0.0, s_m2 = 0.0, t_mean = 0.0, t_m2 = 0.0;
  long n = 0;
  double probe, start;
  nfast_timing r;

  /* No usable wall clock (wall_seconds() == 0): a budgeted loop would never
   * terminate, so degrade to a single tick-only measurement. */
  probe = wall_seconds();
  if (probe == 0.0) {
    ticks t0, t1;
    t0 = getticks();
    fn(ctx);
    t1 = getticks();
    r.secs_mean = 0.0;
    r.secs_std = 0.0;
    r.tks_mean = elapsed(t1, t0);
    r.tks_std = 0.0;
    r.runs = 1;
    return r;
  }

  /* Re-run until NFAST_MEASURE_SECONDS of wall time accumulate. The loop
   * condition reads the same monotonic clock, so it terminates even when a
   * single run is below the clock resolution. */
  start = probe;
  do {
    ticks t0, t1;
    double s0, s1, ds, dt, d, d2;
    s0 = wall_seconds();
    t0 = getticks();
    fn(ctx);
    t1 = getticks();
    s1 = wall_seconds();
    ds = s1 - s0;
    dt = elapsed(t1, t0);
    n++;
    d = ds - s_mean;
    s_mean += d / (double)n;
    d2 = ds - s_mean;
    s_m2 += d * d2;
    d = dt - t_mean;
    t_mean += d / (double)n;
    d2 = dt - t_mean;
    t_m2 += d * d2;
  } while (wall_seconds() - start < NFAST_MEASURE_SECONDS);

  r.secs_mean = s_mean;
  r.secs_std = n > 1 ? sqrt(s_m2 / (double)n) : 0.0;
  r.tks_mean = t_mean;
  r.tks_std = n > 1 ? sqrt(t_m2 / (double)n) : 0.0;
  r.runs = n;
  return r;
}

/* One timing row. has_err=1 appends the forward accuracy figure on the same
 * line (used only on the "fwd" row, so accuracy sits beside its own timing). */
static void print_timing_row(const char *label, const char *op, nfast_timing s,
                             int has_err, NFFT_R err, int ok) {
  printf("  %-26s %-4s %12.6e +/- %10.4e s   %14.0f +/- %12.4e tk   n=%ld",
         label, op, s.secs_mean, s.secs_std, s.tks_mean, s.tks_std, s.runs);
  if (has_err)
    printf("   err %12.4e bound %g %s", (double)err, (double)NFAST_NATIVE_BOUND,
           ok ? "PASS" : "FAIL");
  printf("\n");
}

static void print_timing(const char *label, nfast_timing pre, nfast_timing fwd,
                         nfast_timing adj, NFFT_R err, int ok) {
  print_timing_row(label, "pre", pre, 0, NFFT_K(0.0), 0);
  print_timing_row("", "fwd", fwd, 1, err, ok);
  print_timing_row("", "adj", adj, 0, NFFT_K(0.0), 0);
}

int main(void) {
  int d;
  NFFT_INT N, M, j;
  NFFT_R *x;
  NFFT_C *f_hat, *f_ref;

  int Narr[1];
  NFFT_INT n[1];
  int m = WINDOW_M;
  int window;

  NFFT(plan)
  lp, lp_m;

  NFFT(plan_ng) * p_dir, *p_native, *p_native_m;
  NFFT_C *f_dir, *f_native, *f_native_m, *f_hat_adj;

  nfast_timing t_pre_legacy, t_pre_none, t_ld, t_lf, t_adj_ld, t_adj_lf;
  nfast_timing t_pre_dir, t_dir, t_adj_dir;
  nfast_timing t_pre_native, t_native, t_adj_native;
  /* FFTW_MEASURE variants of the two fast NFFTs (legacy + native). */
  nfast_timing t_pre_legacy_m, t_lf_m, t_adj_lf_m;
  nfast_timing t_pre_native_m, t_native_m, t_adj_native_m;
  NFFT_R err_legacy_direct_file, err_legacy_file, err_native_file, err_dir_file,
      err_legacy_meas, err_native_meas;
  int ok_legacy_direct_file, ok_legacy_file, ok_native_file, ok_dir_file,
      ok_legacy_meas, ok_native_meas, all_ok;
  plan_ng_adjoint_ctx adj_ctx;

  if (!read_reference(NFAST_NATIVE_DATA, &d, &N, &M, &x, &f_hat, &f_ref)) {
    fprintf(stderr, "nfast_native: failed to read reference data\n");
    return EXIT_FAILURE;
  }

  if (d != 1) {
    fprintf(stderr, "nfast_native: expected d=1, got d=%d\n", d);
    free(x);
    free(f_hat);
    free(f_ref);
    return EXIT_FAILURE;
  }

  Narr[0] = (int)N;
  n[0] = 2 * N;
  window = NFFT(get_window_id)();

  printf("nfast_native: d=%d N=%td M=%td m=%d n=%td window=%s\n", d,
         (ptrdiff_t)N, (ptrdiff_t)M, m, (ptrdiff_t)n[0], NFFT(get_window_name)());

  /* --- legacy NFFT: one plan shared by the direct and fast runs, both
   * directions -------------------------------------------------------- */
  NFFT(init_guru)
  (&lp, 1, Narr, (int)M, (int *)n, m,
   PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
       FFTW_INIT | FFT_OUT_OF_PLACE,
   FFTW_ESTIMATE);
  for (j = 0; j < (NFFT_INT)d * M; j++)
    lp.x[j] = x[j];
  for (j = 0; j < N; j++)
    lp.f_hat[j] = f_hat[j];

  t_pre_legacy = time_run(run_legacy_precompute, &lp);
  t_pre_none.secs_mean = 0.0;
  t_pre_none.secs_std = 0.0;
  t_pre_none.tks_mean = 0.0;
  t_pre_none.tks_std = 0.0;
  t_pre_none.runs = 0;

  /* trafo_direct ignores the PRE_PSI table (exact NDFT); capture before
   * trafo overwrites lp.f. */
  t_ld = time_run(run_legacy_direct, &lp);
  err_legacy_direct_file = rel_max_err(lp.f, f_ref, M);

  t_lf = time_run(run_legacy_fast, &lp);
  err_legacy_file = rel_max_err(lp.f, f_ref, M);

  /* Adjoint direction: feed the file's f[] reference back in as input (the
   * PRE_PSI table above already covers both directions). */
  for (j = 0; j < M; j++)
    lp.f[j] = f_ref[j];
  t_adj_ld = time_run(run_legacy_adjoint_direct, &lp);
  t_adj_lf = time_run(run_legacy_adjoint_fast, &lp);

  /* --- planner-native direct NDFT ------------------------------------- */
  f_dir = (NFFT_C *)malloc((size_t)M * sizeof(NFFT_C));
  p_dir = NFFT(plan_ng_guru)(1, &N, NULL, n, M, m, window, x, f_hat, f_dir, 0u,
                             NFFT_ESTIMATE | NFFT_NO_FAST_NATIVE);
  /* "(adj (null))" is expected here: the adjoint direction reuses the
   * forward winner rather than racing its own plan (Y(plan_ng_print)). */
  printf("\nplan: planner direct\n");
  NFFT(fprint_plan)
  (p_dir, stdout);
  printf("\n");
  t_pre_dir = time_run(run_precompute, p_dir);
  t_dir = time_run(run_plan_ng, p_dir);

  /* --- composed planner-native fast NFFT solver ----------------------- */
  f_native = (NFFT_C *)malloc((size_t)M * sizeof(NFFT_C));
  p_native = NFFT(plan_ng_guru)(1, &N, NULL, n, M, m, window, x, f_hat,
                                f_native, 0u, NFFT_ESTIMATE | NFFT_NO_DIRECT);
  printf("\nplan: planner native fast\n");
  NFFT(fprint_plan)
  (p_native, stdout);
  printf("\n");
  t_pre_native = time_run(run_precompute, p_native);
  t_native = time_run(run_plan_ng, p_native);

  /* --- FFTW_MEASURE variants of the two fast NFFTs -------------------------
   * Same two fast transforms as above, but with FFTW_MEASURE (instead of
   * FFTW_ESTIMATE) for the internal DFT plan.  Built here, AFTER every
   * FFTW_ESTIMATE plan (lp, p_native) is already constructed, so the estimate
   * plans cannot inherit measured FFTW wisdom (FFTW's planner caches wisdom
   * process-globally).  Plan construction -- where FFTW_MEASURE actually spends
   * its time -- is not timed; only precompute/forward/adjoint are. */

  /* legacy fast, FFTW_MEASURE.  FFTW_MEASURE clobbers the internal g1/g2 during
   * planning (allocated inside init_guru), never f_hat/f, so fill f_hat after. */
  NFFT(init_guru)
  (&lp_m, 1, Narr, (int)M, (int *)n, m,
   PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F |
       FFTW_INIT | FFT_OUT_OF_PLACE,
   FFTW_MEASURE);
  for (j = 0; j < (NFFT_INT)d * M; j++)
    lp_m.x[j] = x[j];
  for (j = 0; j < N; j++)
    lp_m.f_hat[j] = f_hat[j];
  t_pre_legacy_m = time_run(run_legacy_precompute, &lp_m);
  t_lf_m = time_run(run_legacy_fast, &lp_m);
  err_legacy_meas = rel_max_err(lp_m.f, f_ref, M);
  for (j = 0; j < M; j++)
    lp_m.f[j] = f_ref[j];
  t_adj_lf_m = time_run(run_legacy_adjoint_fast, &lp_m);

  /* native fast, FFTW_MEASURE (nfft-nd.c turns FFTW_MEASURE into
   * FFTW_MEASURE | FFTW_DESTROY_INPUT, mirroring the estimate path). */
  f_native_m = (NFFT_C *)malloc((size_t)M * sizeof(NFFT_C));
  p_native_m = NFFT(plan_ng_guru)(1, &N, NULL, n, M, m, window, x, f_hat,
                                  f_native_m, (unsigned)FFTW_MEASURE,
                                  NFFT_ESTIMATE | NFFT_NO_DIRECT);
  t_pre_native_m = time_run(run_precompute, p_native_m);
  t_native_m = time_run(run_plan_ng, p_native_m);

  /* Adjoint direction: execute_adjoint_on takes explicit buffers, so the
   * plan's own forward-bound f_hat/f_dir/f_native are left untouched. */
  f_hat_adj = (NFFT_C *)malloc((size_t)N * sizeof(NFFT_C));
  adj_ctx.f = f_ref;
  adj_ctx.f_hat = f_hat_adj;
  adj_ctx.p = p_dir;
  t_adj_dir = time_run(run_plan_ng_adjoint, &adj_ctx);
  adj_ctx.p = p_native;
  t_adj_native = time_run(run_plan_ng_adjoint, &adj_ctx);
  adj_ctx.p = p_native_m;
  t_adj_native_m = time_run(run_plan_ng_adjoint, &adj_ctx);

  /* --- Compare (forward error only, every figure vs the file reference) */
  err_dir_file = rel_max_err(f_dir, f_ref, M);
  err_native_file = rel_max_err(f_native, f_ref, M);
  err_native_meas = rel_max_err(f_native_m, f_ref, M);

  ok_legacy_direct_file = err_legacy_direct_file <= NFAST_NATIVE_BOUND;
  ok_legacy_file = err_legacy_file <= NFAST_NATIVE_BOUND;
  ok_legacy_meas = err_legacy_meas <= NFAST_NATIVE_BOUND;
  ok_dir_file = err_dir_file <= NFAST_NATIVE_BOUND;
  ok_native_file = err_native_file <= NFAST_NATIVE_BOUND;
  ok_native_meas = err_native_meas <= NFAST_NATIVE_BOUND;
  all_ok = ok_legacy_direct_file && ok_legacy_file && ok_legacy_meas &&
           ok_dir_file && ok_native_file && ok_native_meas;

  /* Timing and accuracy on the same rows (fwd row carries the forward error).
   * Ordered so the two directs sit together (legacy, planner) and the two fast
   * NFFTs are adjacent (legacy, planner native), each fast NFFT showing its
   * FFTW_ESTIMATE row immediately above its FFTW_MEASURE row. */
  printf("\ntiming (mean +/- std over ~%g s per step; wall seconds / CPU ticks)"
         " with forward accuracy vs reference (max|got-ref|/max|ref|):\n",
         NFAST_MEASURE_SECONDS);
  print_timing("legacy direct", t_pre_none, t_ld, t_adj_ld,
               err_legacy_direct_file, ok_legacy_direct_file);
  print_timing("planner direct", t_pre_dir, t_dir, t_adj_dir, err_dir_file,
               ok_dir_file);
  print_timing("legacy fast est", t_pre_legacy, t_lf, t_adj_lf, err_legacy_file,
               ok_legacy_file);
  print_timing("legacy fast meas", t_pre_legacy_m, t_lf_m, t_adj_lf_m,
               err_legacy_meas, ok_legacy_meas);
  print_timing("planner native fast est", t_pre_native, t_native, t_adj_native,
               err_native_file, ok_native_file);
  print_timing("planner native fast meas", t_pre_native_m, t_native_m,
               t_adj_native_m, err_native_meas, ok_native_meas);

  printf("%s\n", all_ok ? "PASS" : "FAIL");

  /* lp.f/lp.f_hat are plan-owned (MALLOC_F/MALLOC_F_HAT above), freed by
   * NFFT(finalize). */
  NFFT(finalize)
  (&lp);
  NFFT(finalize)
  (&lp_m);
  NFFT(plan_ng_destroy)
  (p_dir);
  NFFT(plan_ng_destroy)
  (p_native);
  NFFT(plan_ng_destroy)
  (p_native_m);

  free(f_dir);
  free(f_native);
  free(f_native_m);
  free(f_hat_adj);
  free(x);
  free(f_hat);
  free(f_ref);

  return all_ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
