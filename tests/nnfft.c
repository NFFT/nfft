/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "nnfft.h"
#include "accuracy_log.h"
#include "util.h"

#define ABSPATH(x) ABS_SRCDIR "/tests/" x
#define SEED 1234567890L

/* Per-axis truncation term and round-off amplification phi_hat(0) /
 * phi_hat(N/2) on grid n of one NFFT with oversampling sigma. Truncation
 * terms are those of the NFFT tests, with the actual sigma. */
#if defined(KAISER_BESSEL)
static R trunc_w(const int m, const R sigma)
{
  const R z = K2PI * ((R)m + K(1.0)) * SQRT(K(1.0) - K(1.0) / sigma);
  return K(3.5) * EXP(-z) / Y(bessel_i0_exp_scaled)(z);
}

static R amp_w(const int m, const R sigma, const R N, const R n)
{
  const R b = KPI * (K(2.0) - K(1.0) / sigma), r = (R)m + K(1.0);
  return Y(kb_phi_hut)(b, K(1.0), r, n, K(0.0)) / Y(kb_phi_hut)(b, K(1.0), r, n, N / K(2.0));
}
#elif defined(GAUSSIAN)
static R trunc_w(const int m, const R sigma)
{
  return K(0.41) * EXP(-(R)m * KPI * (K(1.0) - K(1.0) / (K(2.0) * sigma - K(1.0))));
}

static R amp_w(const int m, const R sigma, const R N, const R n)
{
  const R b = K(2.0) * sigma / (K(2.0) * sigma - K(1.0)) * ((R)m + K(1.0)) / KPI;
  const R x = KPI * N / (K(2.0) * n);
  return EXP(b * x * x);
}
#elif defined(B_SPLINE)
static R trunc_w(const int m, const R sigma)
{
  return K(12000.0) * POW(K(2.0) * sigma - K(1.0), -K(2.0) * (R)m);
}

static R amp_w(const int m, const R sigma, const R N, const R n)
{
  const R x = KPI * N / (K(2.0) * n);
  UNUSED(sigma);
  return POW(x / SIN(x), K(2.0) * (R)m + K(2.0));
}
#elif defined(SINC_POWER)
static R trunc_w(const int m, const R sigma)
{
  return K(4e-3) / ((R)m - K(1.0)) * (K(2.0) * POW(sigma, -K(2.0) * (R)m)
      + POW(sigma / (K(2.0) * sigma - K(1.0)), K(2.0) * (R)m));
}

/* phi_hat is the cardinal B-spline of order 2m, centred, scaled by w. */
static R amp_w(const int m, const R sigma, const R N, const R n)
{
  const R w = (sigma - SINCPOW_A) / ((R)m * sigma);
  return Y(bsplines)(2 * m, (R)m) / Y(bsplines)(2 * m, (R)m + N / (K(2.0) * w * n));
}
#endif

/* Round-off of the direct transform: it forms the phase v x N in double
 * without reduction, so the error grows with N. */
static R err_direct(const nnfft_plan *p)
{
  R sum_N = K(0.0);
  int t;

  for (t = 0; t < p->d; t++)
    sum_N += (R)p->N[t];

  return (K(56.0) + KPI / K(2.0) * sum_N) * Y(float_property)(NFFT_EPSILON);
}

/* NNFFT = outer window convolution onto the aN1 grid, inner NFFT on that
 * grid, deconvolution at x. Outer truncation plus inner truncation amplified
 * by the outer deconvolution (Kircheis, Potts, Tasche 2023, Thm 3.5, with
 * the sharper factor of Elbel, Steidl 1998), plus the round-off of both
 * deconvolutions and of the direct transform. */
static R err_trafo(const nnfft_plan *p)
{
  const nfft_plan *in = p->direct_plan;
  R t1 = K(0.0), t2 = K(0.0), a1 = K(1.0), a2 = K(1.0);
  int t;

  for (t = 0; t < p->d; t++)
  {
    const R s1 = (R)p->N1[t] / (R)p->N[t], s2 = (R)in->n[t] / (R)p->aN1[t];
    t1 += trunc_w(p->m, s1);
    t2 += trunc_w(p->m, s2);
    a1 *= amp_w(p->m, s1, (R)p->N[t], (R)p->N1[t]);
    a2 *= amp_w(p->m, s2, (R)p->aN1[t], (R)in->n[t]);
  }

  return t1 + a1 * t2 + K(3.5) * (a1 + a2) * Y(float_property)(NFFT_EPSILON)
      + err_direct(p);
}

/* The Gaussian edge error stays below the bound from m = 9 on, with or
 * without the full run on the grid. B-spline and sinc-power windows have
 * compact support and show no edge error at any m. */
#if defined(GAUSSIAN)
#define EDGE_M 4
#else
#define EDGE_M WINDOW_HELP_ESTIMATE_m
#endif

/* Test case delegate. setup allocates N, x, v, f_hat, f, destroy frees
 * them; fill writes input and reference once the plan exists. */
typedef struct testcase_delegate_s testcase_delegate_t;
typedef void (*setup_t)(const testcase_delegate_t *ego_, int *d, int **N, int *M, int *NV, R **x, R **v, C **f_hat, C **f);
typedef void (*destroy_t)(int *N, R *x, R *v, C *f_hat, C *f);
typedef void (*fill_t)(const testcase_delegate_t *ego_, nnfft_plan *p, const int adjoint, C *f_hat, C *f);

struct testcase_delegate_s
{
  setup_t setup;
  destroy_t destroy;
  fill_t fill;
};

typedef struct testcase_delegate_file_s
{
  setup_t setup;
  destroy_t destroy;
  fill_t fill;
  const char *filename;
} testcase_delegate_file_t;

typedef struct testcase_delegate_online_s
{
  setup_t setup;
  destroy_t destroy;
  fill_t fill;
  const int d;
  const int N;
  const int M; /* also the number of frequency nodes NV */
  /* Worst-case input: the first frequency node at a band-edge corner, v0 on
   * every axis, with a unit coefficient there, or a unit sample at the first
   * node for the adjoint. With l1 norm 1 the error is the sup norm the bound
   * describes. 0: random input. */
  const R v0;
} testcase_delegate_online_t;

typedef struct init_delegate_s init_delegate_t;
typedef void (*init_t)(const init_delegate_t *ego, nnfft_plan *p, const int d, int *N, const int M, const int NV);

struct init_delegate_s
{
  const char *name;
  init_t init;
  const int m;
  const unsigned flags;
};

typedef struct trafo_delegate_s
{
  const char *name;
  void (*trafo)(nnfft_plan *p);
  R (*acc)(const nnfft_plan *p);
  const int adjoint;
} trafo_delegate_t;

static void print_size(const char *name, const int d, const int *N, const int M)
{
  int j;

  printf("%-31s d = %-1d, N = [", name, d);
  for (j = 0; j < d; j++)
    printf("%s%-5d", IF(j > 0, ", ", ""), N[j]);
  for (j = 0; j < (4 - d); j++)
    printf("%s%-5s", "  ", "");
  printf("], M = %-5d", M);
#ifdef _OPENMP
  printf(" nthreads = " __D__, NFFT(get_num_threads)());
#endif
}

/* File layout: tests/refgen io_format.py, NNFFT variant. */
static void setup_file(const testcase_delegate_t *ego_, int *d, int **N, int *M, int *NV, R **x, R **v, C **f_hat, C **f)
{
  const testcase_delegate_file_t *ego = (const testcase_delegate_file_t*)ego_;
  FILE *file = fopen(ego->filename, "r");
  int j;

  fscanf(file, "%d", d);
  *N = Y(malloc)((size_t)(*d) * sizeof(int));
  for (j = 0; j < *d; j++)
    fscanf(file, "%d", &((*N)[j]));
  fscanf(file, "%d", M);
  fscanf(file, "%d", NV);

  print_size(strrchr(ego->filename, '/') + 1, *d, *N, *M);

  *x = Y(malloc)((size_t)(*M * *d) * sizeof(R));
  for (j = 0; j < *M * *d; j++)
    fscanf(file, __FI__, &((*x)[j]));
  *v = Y(malloc)((size_t)(*NV * *d) * sizeof(R));
  for (j = 0; j < *NV * *d; j++)
    fscanf(file, __FI__, &((*v)[j]));

  *f_hat = Y(malloc)((size_t)(*NV) * sizeof(C));
  for (j = 0; j < *NV; j++)
  {
    R re, im;
    fscanf(file, __FI__ " " __FI__, &re, &im);
    (*f_hat)[j] = re + im * I;
  }
  *f = Y(malloc)((size_t)(*M) * sizeof(C));
  for (j = 0; j < *M; j++)
  {
    R re, im;
    fscanf(file, __FI__ " " __FI__, &re, &im);
    (*f)[j] = re + im * I;
  }

  fclose(file);
}

/* Frees what setup_file or setup_online allocated. */
static void destroy_file(int *N, R *x, R *v, C *f_hat, C *f)
{
  Y(free)(N);
  Y(free)(x);
  Y(free)(v);
  Y(free)(f_hat);
  Y(free)(f);
}

/* The file holds input and reference, read in setup_file. */
static void fill_file(const testcase_delegate_t *ego_, nnfft_plan *p, const int adjoint, C *f_hat, C *f)
{
  UNUSED(ego_); UNUSED(p); UNUSED(adjoint); UNUSED(f_hat); UNUSED(f);
}

static void setup_online(const testcase_delegate_t *ego_, int *d, int **N, int *M, int *NV, R **x, R **v, C **f_hat, C **f)
{
  const testcase_delegate_online_t *ego = (const testcase_delegate_online_t*)ego_;
  int j;

  Y(srand48)(SEED);

  *d = ego->d;
  *N = Y(malloc)((size_t)(*d) * sizeof(int));
  for (j = 0; j < *d; j++)
    (*N)[j] = ego->N;
  *M = ego->M;
  *NV = ego->M;

  print_size(IF(ego->v0 != K(0.0), "nnfft_online_edge", "nnfft_online"), *d, *N, *M);

  *x = Y(malloc)((size_t)(*M * *d) * sizeof(R));
  for (j = 0; j < *M * *d; j++)
    (*x)[j] = Y(drand48)() - K(0.5);
  *v = Y(malloc)((size_t)(*NV * *d) * sizeof(R));
  for (j = 0; j < *NV * *d; j++)
    (*v)[j] = Y(drand48)() - K(0.5);
  if (ego->v0 != K(0.0))
    for (j = 0; j < *d; j++)
      (*v)[j] = ego->v0;

  *f_hat = Y(malloc)((size_t)(*NV) * sizeof(C));
  *f = Y(malloc)((size_t)(*M) * sizeof(C));
}

/* Random or unit input, reference from the direct transform. */
static void fill_online(const testcase_delegate_t *ego_, nnfft_plan *p, const int adjoint, C *f_hat, C *f)
{
  const testcase_delegate_online_t *ego = (const testcase_delegate_online_t*)ego_;
  C *in = adjoint ? f : f_hat, *out = adjoint ? f_hat : f;
  C *p_in = adjoint ? p->f : p->f_hat, *p_out = adjoint ? p->f_hat : p->f;
  const int n_in = adjoint ? p->M_total : p->N_total, n_out = adjoint ? p->N_total : p->M_total;
  int j;

  for (j = 0; j < n_in; j++)
    in[j] = IF(ego->v0 != K(0.0), IF(j == 0, K(1.0), K(0.0)),
        (Y(drand48)() - K(0.5)) + (Y(drand48)() - K(0.5)) * I);
  for (j = 0; j < n_in; j++)
    p_in[j] = in[j];

  if (adjoint)
    nnfft_adjoint_direct(p);
  else
    nnfft_trafo_direct(p);

  for (j = 0; j < n_out; j++)
    out[j] = p_out[j];
}

static void init_(const init_delegate_t *ego, nnfft_plan *p, const int d, int *N, const int M, const int NV)
{
  UNUSED(ego);
  nnfft_init(p, d, NV, M, N);
}

static void init_guru_(const init_delegate_t *ego, nnfft_plan *p, const int d, int *N, const int M, const int NV)
{
  int t, N1[d];

  for (t = 0; t < d; t++)
    N1[t] = 2 * N[t];
  nnfft_init_guru(p, d, NV, M, N, N1, ego->m, ego->flags | MALLOC_X | MALLOC_V | MALLOC_F | MALLOC_F_HAT);
}

static const init_delegate_t init = {"init", init_, 0, 0};
static const init_delegate_t init_guru_pre_psi = {"init_guru (PRE PSI)", init_guru_, WINDOW_HELP_ESTIMATE_m, PRE_PSI | PRE_PHI_HUT};
static const init_delegate_t init_guru_pre_psi_m_half = {"init_guru (PRE PSI, m/2)", init_guru_, WINDOW_HELP_ESTIMATE_m / 2, PRE_PSI | PRE_PHI_HUT};
static const init_delegate_t init_guru_no_pre = {"init_guru (no PRE_PHI_HUT, no PRE_PSI)", init_guru_, WINDOW_HELP_ESTIMATE_m, 0};
static const init_delegate_t init_guru_edge = {"init_guru (PRE PSI)", init_guru_, EDGE_M, PRE_PSI | PRE_PHI_HUT};

static const init_delegate_t *initializers_direct[] = {&init_guru_pre_psi};
static const init_delegate_t *initializers[] = {&init, &init_guru_pre_psi, &init_guru_pre_psi_m_half, &init_guru_no_pre};
static const init_delegate_t *initializers_edge[] = {&init_guru_edge};

static void trafo_direct_(nnfft_plan *p) { nnfft_trafo_direct(p); }
static void adjoint_direct_(nnfft_plan *p) { nnfft_adjoint_direct(p); }

static const trafo_delegate_t trafo_direct = {"trafo_direct", trafo_direct_, err_direct, 0};
static const trafo_delegate_t trafo = {"trafo", nnfft_trafo, err_trafo, 0};
static const trafo_delegate_t adjoint_direct = {"adjoint_direct", adjoint_direct_, err_direct, 1};
static const trafo_delegate_t adjoint = {"adjoint", nnfft_adjoint, err_trafo, 1};

static int check_single(const testcase_delegate_t *testcase, const init_delegate_t *init_delegate,
    const trafo_delegate_t *trafo_delegate)
{
  nnfft_plan p;
  int ok, d, j, *N, M, NV;
  R *x, *v, err = K(0.0), norm = K(0.0), bound;
  C *f_hat, *f, *in, *out, *p_out;
  int n_in, n_out;
  const int adj = trafo_delegate->adjoint;

  testcase->setup(testcase, &d, &N, &M, &NV, &x, &v, &f_hat, &f);
  printf(", %-28s", init_delegate->name);
  init_delegate->init(init_delegate, &p, d, N, M, NV);
  printf(", m = %2d, %-14s", p.m, trafo_delegate->name);

  for (j = 0; j < M * d; j++)
    p.x[j] = x[j];
  for (j = 0; j < NV * d; j++)
    p.v[j] = v[j];
  nnfft_precompute_one_psi(&p);

  testcase->fill(testcase, &p, adj, f_hat, f);

  in = adj ? f : f_hat;
  out = adj ? f_hat : f;
  n_in = adj ? M : NV;
  n_out = adj ? NV : M;
  for (j = 0; j < n_in; j++)
  {
    (adj ? p.f : p.f_hat)[j] = in[j];
    norm += CABS(in[j]);
  }

  trafo_delegate->trafo(&p);

  p_out = adj ? p.f_hat : p.f;
  for (j = 0; j < n_out; j++)
    err = FMAX(err, CABS(out[j] - p_out[j]));
  err = IF(err == K(0.0), K(0.0), err / norm);

  bound = trafo_delegate->acc(&p);
  ok = IF(err < bound, 1, 0);
  printf(" -> %-4s " __FE__ " (" __FE__ ")\n", IF(ok == 0, "FAIL", "OK"), err, bound);
  accuracy_log_append("nnfft",
    testcase->setup == setup_online
      ? IF(((const testcase_delegate_online_t*)testcase)->v0 != K(0.0), "online_edge", "online") : "file",
    d, N, M, init_delegate->name, trafo_delegate->name, (long double)err, (long double)bound, ok);

  testcase->destroy(N, x, v, f_hat, f);
  nnfft_finalize(&p);
  return ok;
}

static void check_many(const size_t nf, const size_t ni, const testcase_delegate_t **testcases,
    const init_delegate_t **initializers_, const trafo_delegate_t *trafo_delegate)
{
  size_t i, j;
  int ok = 1, r;

  for (i = 0; i < nf; i++)
  {
    for (j = 0; j < ni; j++)
    {
      r = check_single(testcases[i], initializers_[j], trafo_delegate);
      ok = MIN(ok, r);
    }
  }
  CU_ASSERT(ok);
}

#include "data/generated/nnfft_testcases.h"

#define CHECK_FILE(D) \
void nnfft_check_##D##d_direct_file(void) \
{ \
  printf("check_" #D "d_direct_file:\n"); \
  check_many(SIZE(testcases_##D##d_file), SIZE(initializers_direct), testcases_##D##d_file_, initializers_direct, &trafo_direct); \
} \
void nnfft_check_##D##d_fast_file(void) \
{ \
  printf("check_" #D "d_fast_file:\n"); \
  check_many(SIZE(testcases_##D##d_file), SIZE(initializers), testcases_##D##d_file_, initializers, &trafo); \
} \
void nnfft_check_adjoint_##D##d_direct_file(void) \
{ \
  printf("check_adjoint_" #D "d_direct_file:\n"); \
  check_many(SIZE(testcases_adjoint_##D##d_file), SIZE(initializers_direct), testcases_adjoint_##D##d_file_, initializers_direct, &adjoint_direct); \
} \
void nnfft_check_adjoint_##D##d_fast_file(void) \
{ \
  printf("check_adjoint_" #D "d_fast_file:\n"); \
  check_many(SIZE(testcases_adjoint_##D##d_file), SIZE(initializers), testcases_adjoint_##D##d_file_, initializers, &adjoint); \
}

CHECK_FILE(1)
CHECK_FILE(2)
CHECK_FILE(3)

#define ONLINE(name, d, N, M, v0) \
static const testcase_delegate_online_t name = {setup_online, destroy_file, fill_online, d, N, M, v0};

#define V0_UPPER (K(0.5) - K(1.0) / K(1048576.0))
#define V0_LOWER K(-0.5)

ONLINE(online_1d_4, 1, 4, 200, K(0.0))
ONLINE(online_1d_32, 1, 32, 200, K(0.0))
ONLINE(online_1d_200, 1, 200, 200, K(0.0))
ONLINE(online_2d_8, 2, 8, 200, K(0.0))
ONLINE(online_2d_32, 2, 32, 200, K(0.0))
ONLINE(online_3d_8, 3, 8, 200, K(0.0))
ONLINE(online_4d_4, 4, 4, 200, K(0.0))

ONLINE(edge_1d_upper, 1, 32, 50, V0_UPPER)
ONLINE(edge_1d_lower, 1, 32, 50, V0_LOWER)
ONLINE(edge_2d_upper, 2, 16, 50, V0_UPPER)
ONLINE(edge_2d_lower, 2, 16, 50, V0_LOWER)
ONLINE(edge_3d_upper, 3, 8, 50, V0_UPPER)
ONLINE(edge_3d_lower, 3, 8, 50, V0_LOWER)
ONLINE(edge_4d_upper, 4, 4, 50, V0_UPPER)
ONLINE(edge_4d_lower, 4, 4, 50, V0_LOWER)

static const testcase_delegate_online_t *testcases_1d_online[] = {&online_1d_4, &online_1d_32, &online_1d_200};
static const testcase_delegate_online_t *testcases_2d_online[] = {&online_2d_8, &online_2d_32};
static const testcase_delegate_online_t *testcases_3d_online[] = {&online_3d_8};
static const testcase_delegate_online_t *testcases_4d_online[] = {&online_4d_4};
static const testcase_delegate_online_t *testcases_1d_online_edge[] = {&edge_1d_upper, &edge_1d_lower};
static const testcase_delegate_online_t *testcases_2d_online_edge[] = {&edge_2d_upper, &edge_2d_lower};
static const testcase_delegate_online_t *testcases_3d_online_edge[] = {&edge_3d_upper, &edge_3d_lower};
static const testcase_delegate_online_t *testcases_4d_online_edge[] = {&edge_4d_upper, &edge_4d_lower};

#define CHECK_ONLINE(D, kind, inits) \
void nnfft_check_##D##d_online##kind(void) \
{ \
  printf("check_" #D "d_online" #kind ":\n"); \
  check_many(SIZE(testcases_##D##d_online##kind), SIZE(inits), \
    (const testcase_delegate_t**)testcases_##D##d_online##kind, inits, &trafo); \
} \
void nnfft_check_adjoint_##D##d_online##kind(void) \
{ \
  printf("check_adjoint_" #D "d_online" #kind ":\n"); \
  check_many(SIZE(testcases_##D##d_online##kind), SIZE(inits), \
    (const testcase_delegate_t**)testcases_##D##d_online##kind, inits, &adjoint); \
}

CHECK_ONLINE(1, , initializers)
CHECK_ONLINE(2, , initializers)
CHECK_ONLINE(1, _edge, initializers_edge)
CHECK_ONLINE(2, _edge, initializers_edge)
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
CHECK_ONLINE(3, , initializers)
CHECK_ONLINE(4, , initializers)
CHECK_ONLINE(3, _edge, initializers_edge)
CHECK_ONLINE(4, _edge, initializers_edge)
#endif
