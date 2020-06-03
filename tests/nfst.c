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
#include <math.h>
#include <complex.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3.h"
#include "infft.h"
#include "cycle.h"
#include "nfst.h"

#define ABSPATH(x) ABS_SRCDIR "/tests/" x

/* Testcase delegate. */
typedef struct testcase_delegate_s testcase_delegate_t;

typedef void (*setup_t)(const testcase_delegate_t *ego_, int *d, int **N, int *NN, int *M, R **x, R **f_hat, R **f);
typedef void (*destroy_t)(const testcase_delegate_t *ego_, int *N, R *x, R *f_hat, R *f);

struct testcase_delegate_s
{
  setup_t setup;
  destroy_t destroy;
};

typedef struct testcase_delegate_file_s
{
  setup_t setup;
  destroy_t destroy;
  const char *filename;
} testcase_delegate_file_t;

static void setup_file(const testcase_delegate_t *ego_, int *d, int **N, int *NN, int *M, R **x, R **f_hat, R **f);
static void destroy_file(const testcase_delegate_t *ego_, int *N, R *x, R *f_hat, R *f);

typedef struct testcase_delegate_online_s
{
  setup_t setup;
  destroy_t destroy;
  const int d;
  const int N;
  const int M;
} testcase_delegate_online_t;

static void setup_online(const testcase_delegate_t *ego_, int *d, int **N, int *NN, int *M, R **x, R **f_hat, R **f);
static void destroy_online(const testcase_delegate_t *ego_, int *N, R *x, R *f_hat, R *f);

/* Initialization delegate. */
typedef struct init_delegate_s init_delegate_t;
typedef void (*init_t)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);

struct init_delegate_s
{
  const char *name;
  init_t init;
  const int m;
  const unsigned flags;
  const unsigned fftw_flags;
};

/* Prepare delegate. */
typedef struct check_delegate_s check_delegate_t;
typedef void (*prepare_t)(check_delegate_t *ego, X(plan) *p, const int NN, const int M, const R *f, const R *f_hat);
typedef R (*compare_t)(check_delegate_t *ego, X(plan) *p, const int NN, const int M, const R *f, const R *f_hat);

struct check_delegate_s
{
  prepare_t prepare;
  compare_t compare;
};

/* Trafo delegate. */
typedef void (*trafo_t)(X(plan) *p);
typedef R (*cost_t)(X(plan) *p);
typedef const char* (*check_t)(X(plan) *p);
typedef R (*acc_t)(X(plan) *p);

typedef struct trafo_delegate_s
{
  const char *name;
  trafo_t trafo;
  check_t check;
  cost_t cost;
  acc_t acc;

} trafo_delegate_t;

static R trafo_direct_cost(X(plan) *p);

static R err_trafo(X(plan) *p);
static R err_trafo_direct(X(plan) *p);

/* Check single test case.*/
static int check_single(const testcase_delegate_t *testcase,
    init_delegate_t *init_delegate, check_delegate_t *check_delegate,
    trafo_delegate_t *trafo_delegate);

/* Check multiple test cases.*/
static void check_many(const size_t nf, const size_t ni, const size_t nt,
  const testcase_delegate_t **testcases, init_delegate_t **initializers,
  check_delegate_t *check_delegate, trafo_delegate_t **trafos);

/* Initializers. */
static void init_1d_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);
static void init_2d_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);
static void init_3d_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);
static void init_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);
static void init_advanced_pre_psi_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M);

#define DEFAULT_NFFT_FLAGS MALLOC_X | MALLOC_F | MALLOC_F_HAT | FFTW_INIT | FFT_OUT_OF_PLACE
#define DEFAULT_FFTW_FLAGS FFTW_ESTIMATE | FFTW_DESTROY_INPUT

static init_delegate_t init_1d;
static init_delegate_t init_2d;
static init_delegate_t init_3d;
static init_delegate_t init;
static init_delegate_t init_advanced_pre_psi;
static init_delegate_t init_advanced_pre_full_psi;
static init_delegate_t init_advanced_pre_lin_psi;
#if defined(GAUSSIAN)
static init_delegate_t init_advanced_pre_fg_psi;
#endif

static check_delegate_t check_trafo;
static check_delegate_t check_adjoint;

static trafo_delegate_t trafo_direct;
static trafo_delegate_t trafo;

static R trafo_direct_cost_factor = K(1.0E-6);

static R trafo_direct_cost(X(plan) *p)
{
  if (trafo_direct_cost_factor == 0.0)
  {
    int M, d, Nd, x = 0;
    for (d = 1; d <= 4; d++)
    {
      for (Nd = 4; Nd < 128; Nd *= 2)
      {
        for (M = 4; M <= 128; M *= 2)
        {
          X(plan) p2;
          int *N = Y(malloc)((size_t)(d) * sizeof(int)), i;
          for (i = 0; i < d; i++)
          {
            N[i] = Nd;
          }
          X(init)(&p2, d, N, M);
          for (i = 0; i < M; i++)
            p2.x[i] = K(0.0);
          if(p2.flags & PRE_ONE_PSI)
            X(precompute_one_psi)(&p2);
          for (i = 0; i < d * (Nd - 1); i++)
          {
            p2.f_hat[i] = K(0.0);
          }
          {
            R r;
            ticks t0, t1;
            t0 = getticks();
            X(trafo_direct)(&p2);
            t1 = getticks();
            r = Y(elapsed_seconds)(t1, t0) / (R)(M);
            for (i = 0; i < d; i++)
              r = r / Nd;
            trafo_direct_cost_factor += r;
            printf(__FE__ "\n", r);
            x += 1;
          }
          X(finalize)(&p2);
          Y(free)(N);
        }
      }
    }
    trafo_direct_cost_factor = trafo_direct_cost_factor/((double)x);
    printf("--> " __FE__ "\n", trafo_direct_cost_factor);
  }

  {
    INT c = p->M_total, i;

    for (i = 0; i < p->d; i++)
      c *= p->N[i];

    return trafo_direct_cost_factor * (R)(c);
  }
}

static R err_trafo_direct(X(plan) *p)
{
  UNUSED(p);
  return K(130.0) * Y(float_property)(NFFT_EPSILON);
}

static R err_trafo(X(plan) *p)
{
  const R m = ((R)p->m);
  R s; /* oversampling factor */
  R a;
  R b;
  R eps = Y(float_property)(NFFT_EPSILON);
  R err;
  int i;
  for (i = 0, s = ((R)p->sigma[0]); i < p->d; i++)
    s = FMIN(s, ((R)p->sigma[i]));
#if defined(GAUSSIAN)
#if defined(NFFT_LDOUBLE)
    a = K(3.4);
    b = K(50.0);
#elif defined(NFFT_SINGLE)
    a = K(0.4);
    b = K(11700.0);
#else
    a = K(1.95);
    b = K(50.0);
#endif
    err = EXP(-m*KPI*(K(1.0)-K(1.0)/(K(2.0)*K(2.0) - K(1.0))));
#elif defined(B_SPLINE)
    //printf("m = %E, s = %E, a1 = %E, a2 = %E, z = %E\n", m, s, K(1.0)/(K(2.0)*s-K(1.0)), K(2.0)*m, K(4.0) * POW(K(1.0)/(K(2.0)*s-K(1.0)),K(2.0)*m));
    //printf("\n<s = %E>\n", s);
    //fflush(stdout);
#if defined(NFFT_LDOUBLE)
    a = K(0.3);
    b = K(50.0);
#elif defined(NFFT_SINGLE)
    a = K(0.4);
    b = K(4800.0);
#else
    a = K(1.0);
    b = K(4100.0);
#endif
    err = K(3000.0) * K(4.0) * POW(K(1.0)/(K(2.0)*s-K(1.0)),K(2.0)*m);
  #elif defined(SINC_POWER)
#if defined(NFFT_LDOUBLE)
    a = K(0.3);
    b = K(50.0);
#elif defined(NFFT_SINGLE)
    a = K(0.4);
    b = K(4800.0);
#else
    a = K(1.0);
    b = K(4100.0);
#endif
    err = (K(1.0)/(m-K(1.0))) * ((K(2.0)/(POW(s,K(2.0)*m))) + POW(s/(K(2.0)*s-K(1.0)),K(2.0)*m));
  #elif defined(KAISER_BESSEL)
#if defined(NFFT_LDOUBLE)
    a = K(2.9);
    b = K(50.0);
#elif defined(NFFT_SINGLE)
    a = K(0.95);
    b = K(4800.0);
#else
    a = K(0.7);
    b = K(5000.0);
#endif
    err = KPI * (SQRT(m) + m) * SQRT(SQRT(K(1.0) - K(1.0)/K(2.0))) * EXP(-K2PI * m * SQRT(K(1.0) - K(1.0) / K(2.0)));
  #else
    #error Unsupported window function.
  #endif

  return FMAX(FMAX(a * err, b * eps), err_trafo_direct(p));
}

#define MAX_SECONDS 0.1

static int check_single(const testcase_delegate_t *testcase,
  init_delegate_t *init_delegate, check_delegate_t *check_delegate,
  trafo_delegate_t *trafo_delegate)
{
  int ok = 0;
  X(plan) p;
  int d, j, *N, NN, M;
  R *x;
  R *f_hat, *f;

  testcase->setup(testcase, &d, &N, &NN, &M, &x, &f_hat, &f);

  /* Init plan. */
  printf(", %-28s", init_delegate->name);
  init_delegate->init(init_delegate, &p, d, N, M);

  printf(", m = %2d", (int)p.m);
  printf(", %-14s", trafo_delegate->name);

  /* Nodes. */
  for (j = 0; j < M*d; j++)
  {
    p.x[j] = x[j];
  }

  if (trafo_delegate->check)
  {
    const char* check = trafo_delegate->check(&p);
    if (check != 0)
    {
      printf(" -> %-4s (","OK");
      printf("%s", check);
      printf(")\n");
      ok = 1;
      goto cleanup;
    }
  }
  else if (trafo_delegate->cost)
  {
    const R cost = trafo_delegate->cost(&p);
    if (cost > MAX_SECONDS)
    {
      printf(" -> %-4s (cost too high)\n","OK");
      ok = 1;
      goto cleanup;
    }
  }

  /* Pre-compute Psi, maybe. */
  if(p.flags & PRE_ONE_PSI)
    X(precompute_one_psi)(&p);

  check_delegate->prepare(check_delegate, &p, NN, M, f, f_hat);

  trafo_delegate->trafo(&p);

  /* debug */
  /*fprintf(stderr, "\n");
  for (j = 0; j < M; j++)
    fprintf(stderr, "f[%2d] = " __FE__ ", f[%2d] = " __FE__ ", err = " __FE__ "\n", j,
      f[j], j, p.f[j], ABS(f[j] - p.f[j]) / ABS(f[j]));*/

  /* Standard NFFT error measure. */
  {
    R err = check_delegate->compare(check_delegate, &p, NN, M, f, f_hat);
    R bound = trafo_delegate->acc(&p);
    ok = IF(err < bound, 1, 0);
    printf(" -> %-4s " __FE__ " (" __FE__ ")\n", IF(ok == 0, "FAIL", "OK"), err, bound);
  }

cleanup:
  testcase->destroy(testcase, N, x, f_hat, f);
  X(finalize)(&p);

  return ok;
}

static void check_many(const size_t nf, const size_t ni, const size_t nt,
  const testcase_delegate_t **testcases, init_delegate_t **initializers,
  check_delegate_t *check_delegate, trafo_delegate_t **trafos)
{
  size_t i, j, k;
  int ok = 1, r;
  for (k = 0; k < nt; k++)
  {
    for (i = 0; i < nf; i++)
    {
      for (j = 0; j < ni; j++)
      {
         r = check_single(testcases[i], initializers[j], check_delegate, trafos[k]);
         ok = MIN(ok, r);
      }
    }
  }
  CU_ASSERT(ok);
}

static void setup_file(const testcase_delegate_t *ego_, int *d, int **N, int *NN, int *M, R **x, R **f_hat, R **f)
{
  const testcase_delegate_file_t *ego = (const testcase_delegate_file_t*)ego_;
  int j;
  char filename[200];
  char* c = strrchr(ego->filename, '/');
  FILE *file = fopen(ego->filename, "r");

  filename[0] = (char) 0;
  strncpy(filename, &c[1], 200);
  filename[199] = (char) 0;
  printf("%-31s", filename);

  /* Dimensions. */
  fscanf(file, "%d", d);
  /* Bandwidths. */
  *N = Y(malloc)((size_t)(*d) * sizeof(int));
  for (j = 0; j < *d; j++)
    fscanf(file, "%d", &((*N)[j]));
  /* Number of nodes. */
  fscanf(file, "%d", M);

  printf(" d = %-1d, N = [", *d);
  {
    for (j = 0; j < *d; j++)
    {
      printf("%s%-5d", IF(j > 0,", ", ""), (*N)[j]);
    }
    for (j = 0; j < (4-*d); j++)
    {
      printf("%s%-5s", "  ", "");
    }
  }
  printf("],");
  printf(" M = %-5d", *M);

  for (j = 0, *NN = 1; j < *d; j++)
    *NN *= ((*N)[j] - 1);

  /* Nodes. */
  *x = Y(malloc)((size_t)(M[0]*d[0])*sizeof(R));
  for (j = 0; j < M[0]*d[0]; j++)
  {
    fscanf(file, __FI__, &((*x)[j]));
  }

  /* Fourier coefficients. */
  *f_hat = Y(malloc)((size_t)(NN[0])*sizeof(R));
  for (j = 0; j < NN[0]; j++)
  {
    R re;
    fscanf(file, __FI__, &re);
    (*f_hat)[j] = re;
  }

  /* Reference function values. */
  *f = Y(malloc)((size_t)(M[0]) * sizeof(R));
  for (j = 0; j < M[0]; j++)
  {
    R re;
    fscanf(file, __FI__, &re);
    (*f)[j] = re;
  }

  fclose(file);
}

static void destroy_file(const testcase_delegate_t *ego_, int *N, R *x, R *f_hat, R *f)
{
  UNUSED(ego_);
  Y(free)(N);
  Y(free)(x);
  Y(free)(f_hat);
  Y(free)(f);
}

static void setup_online(const testcase_delegate_t *ego_, int *d, int **N, int *NN, int *M, R **x, R **f_hat, R **f)
{
  const testcase_delegate_online_t *ego = (const testcase_delegate_online_t*)ego_;
  int j;

  /* Dimensions. */
  *d = ego->d;
  /* Bandwidths. */
  *N = Y(malloc)((size_t)(*d) * sizeof(int));
  for (j = 0; j < *d; j++)
    (*N)[j] = ego->N;
  /* Number of nodes. */
  *M = ego->M;

  printf("%-31s", "nfst_online");

  printf(" d = %-1d, N = [", *d);
  {
    for (j = 0; j < *d; j++)
    {
      printf("%s%-5d", IF(j > 0,", ", ""), (*N)[j]);
    }
    for (j = 0; j < (4-*d); j++)
    {
      printf("%s%-5s", "  ", "");
    }
  }
  printf("],");
  printf(" M = %-5d", *M);

  for (j = 0, *NN = 1; j < *d; j++)
    *NN *= ((*N)[j] - 1);

  /* Nodes. */
  *x = Y(malloc)((size_t)(M[0]*d[0])*sizeof(R));
  for (j = 0; j < M[0]*d[0]; j++)
  {
    (*x)[j] = K(0.5) * Y(drand48)();
  }

  /* Fourier coefficients. */
  *f_hat = Y(malloc)((size_t)(NN[0])*sizeof(R));
  for (j = 0; j < NN[0]; j++)
  {
    (*f_hat)[j] = Y(drand48)() - K(0.5);
  }

  {
    X(plan) p;

    X(init)(&p, *d, *N, *M);

    /* Nodes. */
    for (j = 0; j < M[0]*d[0]; j++)
    {
      p.x[j] = (*x)[j];
    }

    /* Pre-compute Psi, maybe. */
    if(p.flags & PRE_ONE_PSI)
      X(precompute_one_psi)(&p);

    /* Fourier coefficients. */
    for (j = 0; j < *NN; j++)
    {
      p.f_hat[j] = (*f_hat)[j];
    }

    X(trafo_direct)(&p);

    /* Reference function values. */
    *f = Y(malloc)((size_t)(M[0]) * sizeof(R));
    for (j = 0; j < M[0]; j++)
    {
      (*f)[j] = p.f[j];
    }

    X(finalize)(&p);
  }
}

static void setup_adjoint_online(const testcase_delegate_t *ego_, int *d, int **N, int *NN, int *M, R **x, R **f_hat, R **f)
{
  const testcase_delegate_online_t *ego = (const testcase_delegate_online_t*)ego_;
  int j;

  /* Dimensions. */
  *d = ego->d;

  /* Bandwidths. */
  *N = Y(malloc)((size_t)(*d) * sizeof(int));

  for (j = 0; j < *d; j++)
    (*N)[j] = ego->N;

  /* Number of nodes. */
  *M = ego->M;

  printf("%-31s", "nfst_online");

  printf(" d = %-1d, N = [", *d);
  {
    for (j = 0; j < *d; j++)
    {
      printf("%s%-5d", IF(j > 0,", ", ""), (*N)[j]);
    }
    for (j = 0; j < (4-*d); j++)
    {
      printf("%s%-5s", "  ", "");
    }
  }
  printf("],");
  printf(" M = %-5d", *M);

  for (j = 0, *NN = 1; j < *d; j++)
    *NN *= ((*N)[j] - 1);

  /* Nodes. */
  *x = Y(malloc)((size_t)(M[0]*d[0])*sizeof(R));

  for (j = 0; j < M[0]*d[0]; j++)
  {
    (*x)[j] = K(0.5) * Y(drand48)();
  }

  /* Function values. */
  *f = Y(malloc)((size_t)(M[0]) * sizeof(C));
  for (j = 0; j < M[0]; j++)
  {
    (*f)[j] = Y(drand48)() - K(0.5);
  }

  {
    X(plan) p;

    X(init)(&p, *d, *N, *M);

    /* Nodes. */
    for (j = 0; j < M[0]*d[0]; j++)
    {
      p.x[j] = (*x)[j];
    }

    /* Pre-compute Psi, maybe. */
    if(p.flags & PRE_ONE_PSI)
      X(precompute_one_psi)(&p);

    /* Function values. */
    for (j = 0; j < M[0]; j++)
    {
      p.f[j] = (*f)[j];
    }

    X(adjoint_direct)(&p);

    /* Reference pseudo Fourier coefficients. */
    *f_hat = Y(malloc)((size_t)(NN[0])*sizeof(R));

    for (j = 0; j < NN[0]; j++)
    {
      (*f_hat)[j] = p.f_hat[j];
    }

    X(finalize)(&p);
  }
}

static void destroy_online(const testcase_delegate_t *ego_, int *N, R *x, R *f_hat, R *f)
{
  UNUSED(ego_);
  Y(free)(N);
  Y(free)(x);
  Y(free)(f_hat);
  Y(free)(f);
}

/* Initializers. */
static void init_1d_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  UNUSED(ego);
  UNUSED(d);
  X(init_1d)(p, N[0], M);
}

static void init_2d_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  UNUSED(ego);
  UNUSED(d);
  X(init_2d)(p, N[0], N[1], M);
}

static void init_3d_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  UNUSED(ego);
  UNUSED(d);
  X(init_3d)(p, N[0], N[1], N[2], M);
}

static void init_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  UNUSED(ego);
  X(init)(p, d, N, M);
}

static void init_advanced_pre_psi_(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  int *n = Y(malloc)((size_t)(d)*sizeof(int));
  int i;
  for (i = 0; i < d; i++)
    n[i] = 2 * (int)(Y(next_power_of_2)(N[i]));
  X(init_guru)(p, d, N, M, n, ego->m, ego->flags, ego->fftw_flags);
  Y(free)(n);
}

static init_delegate_t init_direct = {"init_guru ()", init_advanced_pre_psi_, WINDOW_HELP_ESTIMATE_m, (DEFAULT_NFFT_FLAGS ^ PRE_PSI), DEFAULT_FFTW_FLAGS};
static init_delegate_t init_1d = {"init_1d", init_1d_, 0, 0, 0};
static init_delegate_t init_2d = {"init_2d", init_2d_, 0, 0, 0};
static init_delegate_t init_3d = {"init_3d", init_3d_, 0, 0, 0};
static init_delegate_t init = {"init", init_, 0, 0, 0};
static init_delegate_t init_advanced_pre_psi = {"init_guru (PRE PSI)", init_advanced_pre_psi_, WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | PRE_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
static init_delegate_t init_advanced_pre_full_psi = {"init_guru (PRE FULL PSI)", init_advanced_pre_psi_, WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | PRE_FULL_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
static init_delegate_t init_advanced_pre_lin_psi = {"init_guru (PRE LIN PSI)", init_advanced_pre_psi_, WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | PRE_LIN_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
#if defined(GAUSSIAN)
static init_delegate_t init_advanced_pre_fg_psi = {"init_guru (PRE FG PSI)", init_advanced_pre_psi_, WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | FG_PSI | PRE_FG_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
#endif

/* Check routines. */
static void prepare_trafo(check_delegate_t *ego, X(plan) *p, const int NN, const int M, const R *f, const R *f_hat)
{
  UNUSED(ego);
  UNUSED(M);
  UNUSED(f);
  int j;

  /* Fourier coefficients. */
  for (j = 0; j < NN; j++)
  {
    p->f_hat[j] = f_hat[j];
  }
}

static void prepare_adjoint(check_delegate_t *ego, X(plan) *p, const int NN, const int M, const R *f, const R *f_hat)
{
  UNUSED(ego);
  UNUSED(NN);
  UNUSED(f_hat);
  int j;

  /* Fourier coefficients. */
  for (j = 0; j < M; j++)
  {
    p->f[j] = f[j];
  }
}

static R compare_trafo(check_delegate_t *ego, X(plan) *p, const int NN, const int M, const R *f, const R *f_hat)
{
  UNUSED(ego);
  UNUSED(f_hat);
  int j;
  R numerator = K(0.0), denominator = K(0.0);

  /* debug */
//  fprintf(stderr, "\n");
//  for (j = 0; j < M; j++)
//    fprintf(stderr, "f[%2d] = " __FE__ ", f[%2d] = " __FE__ ", err = " __FE__ "\n", j,
//      f[j], j, p->f[j], ABS(f[j] - p->f[j]) / ABS(f[j]));

  for (j = 0; j < M; j++)
    numerator = MAX(numerator, ABS(f[j] - p->f[j]));

  for (j = 0; j < NN; j++)
    denominator += ABS(p->f_hat[j]);

  return numerator == K(0.0) ? K(0.0) : numerator/denominator;
}

static R compare_adjoint(check_delegate_t *ego, X(plan) *p, const int NN, const int M, const R *f, const R *f_hat)
{
  UNUSED(ego);
  UNUSED(f);
  int j;
  R numerator = K(0.0), denominator = K(0.0);

  /* debug */
//  fprintf(stderr, "\n");
//  for (j = 0; j < NN; j++)
//    fprintf(stderr, "f_hat[%2d] = " __FE__ ", f_hat[%2d] = " __FE__ ", err = " __FE__ "\n", j,
//      f_hat[j], j, p->f_hat[j], ABS(f_hat[j] - p->f_hat[j]) / ABS(f_hat[j]));

  for (j = 0; j < NN; j++)
    numerator = MAX(numerator, ABS(f_hat[j] - p->f_hat[j]));

  for (j = 0; j < M; j++)
    denominator += ABS(p->f[j]);

  return numerator == K(0.0) ? K(0.0) : numerator/denominator;
}

static check_delegate_t check_trafo = {prepare_trafo, compare_trafo};
static check_delegate_t check_adjoint = {prepare_adjoint, compare_adjoint};

static trafo_delegate_t trafo_direct = {"trafo_direct", X(trafo_direct), 0, trafo_direct_cost, err_trafo_direct};
static trafo_delegate_t trafo = {"trafo", X(trafo), X(check), 0, err_trafo};
//static trafo_delegate_t trafo_1d = {"trafo_1d", X(trafo_1d), X(check), 0, err_trafo};
//static trafo_delegate_t trafo_2d = {"trafo_2d", X(trafo_2d), X(check), 0, err_trafo};
//static trafo_delegate_t trafo_3d = {"trafo_3d", X(trafo_3d), X(check), 0, err_trafo};

static trafo_delegate_t adjoint_direct = {"adjoint_direct", X(adjoint_direct), 0, trafo_direct_cost, err_trafo_direct};
static trafo_delegate_t adjoint = {"adjoint", X(adjoint), X(check), 0, err_trafo};
//static trafo_delegate_t adjoint_1d = {"adjoint_1d", adjoint_1d, X(check), 0, err_trafo};
//static trafo_delegate_t adjoint_2d = {"adjoint_2d", adjoint_2d, X(check), 0, err_trafo};
//static trafo_delegate_t adjoint_3d = {"adjoint_3d", adjoint_3d, X(check), 0, err_trafo};

/* 1D */

/* Initializers. */
static const init_delegate_t* initializers_direct[] =
{
  &init_direct,
};

static const init_delegate_t* initializers_1d[] =
{
  &init_1d,
  &init,
  &init_advanced_pre_psi,
  &init_advanced_pre_full_psi,
//  &init_advanced_pre_lin_psi,
#if defined(GAUSSIAN)
  &init_advanced_pre_fg_psi,
#endif
};

static const testcase_delegate_file_t nfst_1d_2_1 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_2_1.txt")};
static const testcase_delegate_file_t nfst_1d_2_10 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_2_10.txt")};
static const testcase_delegate_file_t nfst_1d_2_25 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_2_25.txt")};
static const testcase_delegate_file_t nfst_1d_2_50 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_2_50.txt")};
static const testcase_delegate_file_t nfst_1d_4_1 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_4_1.txt")};
static const testcase_delegate_file_t nfst_1d_4_10 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_4_10.txt")};
static const testcase_delegate_file_t nfst_1d_4_25 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_4_25.txt")};
static const testcase_delegate_file_t nfst_1d_4_50 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_4_50.txt")};
static const testcase_delegate_file_t nfst_1d_10_1 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_10_1.txt")};
static const testcase_delegate_file_t nfst_1d_10_10 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_10_10.txt")};
static const testcase_delegate_file_t nfst_1d_10_25 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_10_25.txt")};
static const testcase_delegate_file_t nfst_1d_10_50 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_10_50.txt")};
static const testcase_delegate_file_t nfst_1d_25_1 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_25_1.txt")};
static const testcase_delegate_file_t nfst_1d_25_10 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_25_10.txt")};
static const testcase_delegate_file_t nfst_1d_25_25 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_25_25.txt")};
static const testcase_delegate_file_t nfst_1d_25_50 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_25_50.txt")};
static const testcase_delegate_file_t nfst_1d_50_1 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_50_1.txt")};
static const testcase_delegate_file_t nfst_1d_50_10 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_50_10.txt")};
static const testcase_delegate_file_t nfst_1d_50_25 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_50_25.txt")};
static const testcase_delegate_file_t nfst_1d_50_50 = {setup_file,destroy_file,ABSPATH("data/nfst_1d_50_50.txt")};

static const testcase_delegate_file_t *testcases_1d_file[] =
{
    &nfst_1d_2_1,
    &nfst_1d_2_10,
    &nfst_1d_2_25,
    &nfst_1d_2_50,
    &nfst_1d_4_1,
    &nfst_1d_4_10,
    &nfst_1d_4_25,
    &nfst_1d_4_50,
    &nfst_1d_10_1,
    &nfst_1d_10_10,
    &nfst_1d_10_25,
    &nfst_1d_10_50,
    &nfst_1d_25_1,
    &nfst_1d_25_10,
    &nfst_1d_25_25,
    &nfst_1d_25_50,
    &nfst_1d_50_1,
    &nfst_1d_50_10,
    &nfst_1d_50_25,
    &nfst_1d_50_50,
};

static const trafo_delegate_t* trafos_1d_direct_file[] = {&trafo_direct};

void X(check_1d_direct_file)(void)
{
  printf("check_1d_direct_file:\n");
  check_many(SIZE(testcases_1d_file), SIZE(initializers_direct), SIZE(trafos_1d_direct_file),
    testcases_1d_file, initializers_direct, &check_trafo, trafos_1d_direct_file);
}

static const trafo_delegate_t* trafos_1d_fast_file[] = {&trafo/*, &trafo_1d*/};

void X(check_1d_fast_file)(void)
{
  printf("check_1d_fast_file:\n");
  check_many(SIZE(testcases_1d_file), SIZE(initializers_1d), SIZE(trafos_1d_fast_file),
    testcases_1d_file, initializers_1d, &check_trafo, trafos_1d_fast_file);
}

static const testcase_delegate_file_t nfst_adjoint_1d_2_1 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_2_1.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_2_10 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_2_10.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_2_25 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_2_25.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_2_50 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_2_50.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_4_1 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_4_1.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_4_10 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_4_10.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_4_25 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_4_25.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_4_50 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_4_50.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_10_1 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_10_1.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_10_10 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_10_10.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_10_25 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_10_25.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_10_50 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_10_50.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_25_1 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_25_1.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_25_10 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_25_10.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_25_25 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_25_25.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_25_50 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_25_50.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_50_1 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_50_1.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_50_10 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_50_10.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_50_25 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_50_25.txt")};
static const testcase_delegate_file_t nfst_adjoint_1d_50_50 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_1d_50_50.txt")};

static const testcase_delegate_file_t *testcases_adjoint_1d_file[] =
{
  &nfst_adjoint_1d_2_1,
  &nfst_adjoint_1d_2_10,
  &nfst_adjoint_1d_2_25,
  &nfst_adjoint_1d_2_50,
  &nfst_adjoint_1d_4_1,
  &nfst_adjoint_1d_4_10,
  &nfst_adjoint_1d_4_25,
  &nfst_adjoint_1d_4_50,
  &nfst_adjoint_1d_10_1,
  &nfst_adjoint_1d_10_10,
  &nfst_adjoint_1d_10_25,
  &nfst_adjoint_1d_10_50,
  &nfst_adjoint_1d_25_1,
  &nfst_adjoint_1d_25_10,
  &nfst_adjoint_1d_25_25,
  &nfst_adjoint_1d_25_50,
  &nfst_adjoint_1d_50_1,
  &nfst_adjoint_1d_50_10,
  &nfst_adjoint_1d_50_25,
  &nfst_adjoint_1d_50_50,
};

static const trafo_delegate_t* trafos_adjoint_direct_1d_file[] = {&adjoint_direct};

void X(check_adjoint_1d_direct_file)(void)
{
  printf("check_adjoint_1d_direct_file:\n");
  check_many(SIZE(testcases_adjoint_1d_file), SIZE(initializers_direct), SIZE(trafos_adjoint_direct_1d_file),
    testcases_adjoint_1d_file, initializers_direct, &check_adjoint, trafos_adjoint_direct_1d_file);
}

static const trafo_delegate_t* trafos_adjoint_fast_1d_file[] = {&adjoint/*, &adjoint_1d*/};

void X(check_adjoint_1d_fast_file)(void)
{
  printf("check_adjoint_1d_fast_file:\n");
  check_many(SIZE(testcases_adjoint_1d_file), SIZE(initializers_1d), SIZE(trafos_adjoint_fast_1d_file),
    testcases_adjoint_1d_file, initializers_1d, &check_adjoint, trafos_adjoint_fast_1d_file);
}

static const testcase_delegate_online_t nfst_online_1d_50_50 = {setup_online, destroy_online, 1, 50 ,50};
static const testcase_delegate_online_t nfst_online_1d_100_50 = {setup_online, destroy_online, 1, 100 ,50};
static const testcase_delegate_online_t nfst_online_1d_200_50 = {setup_online, destroy_online, 1, 200 ,50};
static const testcase_delegate_online_t nfst_online_1d_500_50 = {setup_online, destroy_online, 1, 500 ,50};
static const testcase_delegate_online_t nfst_online_1d_1000_50 = {setup_online, destroy_online, 1, 1000 ,50};
static const testcase_delegate_online_t nfst_online_1d_2000_50 = {setup_online, destroy_online, 1, 2000 ,50};
static const testcase_delegate_online_t nfst_online_1d_5000_50 = {setup_online, destroy_online, 1, 5000 ,50};
static const testcase_delegate_online_t nfst_online_1d_10000_50 = {setup_online, destroy_online, 1, 10000 ,50};

static const testcase_delegate_online_t *testcases_1d_online[] =
{
  &nfst_online_1d_50_50,
  &nfst_online_1d_100_50,
  &nfst_online_1d_200_50,
  &nfst_online_1d_500_50,
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
  &nfst_online_1d_1000_50,
  &nfst_online_1d_2000_50,
  &nfst_online_1d_5000_50,
  &nfst_online_1d_10000_50,
#endif
};

static const trafo_delegate_t* trafos_1d_online[] = {&trafo/*, &trafo_1d*/};

void X(check_1d_online)(void)
{
  check_many(SIZE(testcases_1d_online), SIZE(initializers_1d), SIZE(trafos_1d_online),
    testcases_1d_online, initializers_1d, &check_trafo, trafos_1d_online);
}

static const testcase_delegate_online_t nfst_adjoint_online_1d_50_50 = {setup_adjoint_online, destroy_online, 1, 50 ,50};
static const testcase_delegate_online_t nfst_adjoint_online_1d_100_50 = {setup_adjoint_online, destroy_online, 1, 100 ,50};
static const testcase_delegate_online_t nfst_adjoint_online_1d_200_50 = {setup_adjoint_online, destroy_online, 1, 200 ,50};
static const testcase_delegate_online_t nfst_adjoint_online_1d_500_50 = {setup_adjoint_online, destroy_online, 1, 500 ,50};
static const testcase_delegate_online_t nfst_adjoint_online_1d_1000_50 = {setup_adjoint_online, destroy_online, 1, 1000 ,50};
static const testcase_delegate_online_t nfst_adjoint_online_1d_2000_50 = {setup_adjoint_online, destroy_online, 1, 2000 ,50};
static const testcase_delegate_online_t nfst_adjoint_online_1d_5000_50 = {setup_adjoint_online, destroy_online, 1, 5000 ,50};
static const testcase_delegate_online_t nfst_adjoint_online_1d_10000_50 = {setup_adjoint_online, destroy_online, 1, 10000 ,50};

static const testcase_delegate_online_t *testcases_adjoint_1d_online[] =
{
  &nfst_adjoint_online_1d_50_50,
  &nfst_adjoint_online_1d_100_50,
  &nfst_adjoint_online_1d_200_50,
  &nfst_adjoint_online_1d_500_50,
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
  &nfst_adjoint_online_1d_1000_50,
  &nfst_adjoint_online_1d_2000_50,
  &nfst_adjoint_online_1d_5000_50,
  &nfst_adjoint_online_1d_10000_50,
#endif
};

static const trafo_delegate_t* trafos_adjoint_1d_online[] = {&adjoint/*, &adjoint_1d*/};

void X(check_adjoint_1d_online)(void)
{
  check_many(SIZE(testcases_adjoint_1d_online), SIZE(initializers_1d), SIZE(trafos_adjoint_1d_online),
    testcases_adjoint_1d_online, initializers_1d, &check_adjoint, trafos_adjoint_1d_online);
}

/* 2D */

/* Initializers. */
static const init_delegate_t* initializers_2d[] =
{
  &init_2d,
  &init,
  &init_advanced_pre_psi,
  &init_advanced_pre_full_psi,
//  &init_advanced_pre_lin_psi,
#if defined(GAUSSIAN)
  &init_advanced_pre_fg_psi,
#endif
};

static const testcase_delegate_file_t nfst_2d_10_10_25 = {setup_file,destroy_file,ABSPATH("data/nfst_2d_10_10_25.txt")};
static const testcase_delegate_file_t nfst_2d_10_10_50 = {setup_file,destroy_file,ABSPATH("data/nfst_2d_10_10_50.txt")};
static const testcase_delegate_file_t nfst_2d_10_25_25 = {setup_file,destroy_file,ABSPATH("data/nfst_2d_10_25_25.txt")};
static const testcase_delegate_file_t nfst_2d_10_25_50 = {setup_file,destroy_file,ABSPATH("data/nfst_2d_10_25_50.txt")};
static const testcase_delegate_file_t nfst_2d_25_10_25 = {setup_file,destroy_file,ABSPATH("data/nfst_2d_25_10_25.txt")};
static const testcase_delegate_file_t nfst_2d_25_10_50 = {setup_file,destroy_file,ABSPATH("data/nfst_2d_25_10_50.txt")};
static const testcase_delegate_file_t nfst_2d_25_25_25 = {setup_file,destroy_file,ABSPATH("data/nfst_2d_25_25_25.txt")};
static const testcase_delegate_file_t nfst_2d_25_25_50 = {setup_file,destroy_file,ABSPATH("data/nfst_2d_25_25_50.txt")};

static const testcase_delegate_file_t *testcases_2d_file[] =
{
  &nfst_2d_10_10_25,
  &nfst_2d_10_10_50,
  &nfst_2d_10_25_25,
  &nfst_2d_10_25_50,
  &nfst_2d_25_10_25,
  &nfst_2d_25_10_50,
  &nfst_2d_25_25_25,
  &nfst_2d_25_25_50,
};

static const trafo_delegate_t* trafos_2d_direct_file[] = {&trafo_direct};

void X(check_2d_direct_file)(void)
{
  printf("check_2d_direct_file:\n");
  check_many(SIZE(testcases_2d_file), SIZE(initializers_direct), SIZE(trafos_2d_direct_file),
    testcases_2d_file, initializers_direct, &check_trafo, trafos_2d_direct_file);
}

static const trafo_delegate_t* trafos_2d_fast_file[] = {&trafo/*, &trafo_2d*/};

void X(check_2d_fast_file)(void)
{
  printf("check_2d_fast_file:\n");
  check_many(SIZE(testcases_2d_file), SIZE(initializers_2d), SIZE(trafos_2d_fast_file),
    testcases_2d_file, initializers_2d, &check_trafo, trafos_2d_fast_file);
}

static const testcase_delegate_file_t nfst_adjoint_2d_10_10_25 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_2d_10_10_25.txt")};
static const testcase_delegate_file_t nfst_adjoint_2d_10_10_50 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_2d_10_10_50.txt")};
static const testcase_delegate_file_t nfst_adjoint_2d_10_25_25 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_2d_10_25_25.txt")};
static const testcase_delegate_file_t nfst_adjoint_2d_10_25_50 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_2d_10_25_50.txt")};
static const testcase_delegate_file_t nfst_adjoint_2d_25_10_25 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_2d_25_10_25.txt")};
static const testcase_delegate_file_t nfst_adjoint_2d_25_10_50 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_2d_25_10_50.txt")};
static const testcase_delegate_file_t nfst_adjoint_2d_25_25_25 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_2d_25_25_25.txt")};
static const testcase_delegate_file_t nfst_adjoint_2d_25_25_50 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_2d_25_25_50.txt")};

static const testcase_delegate_file_t *testcases_adjoint_2d_file[] =
{
  &nfst_adjoint_2d_10_10_25,
  &nfst_adjoint_2d_10_10_50,
  &nfst_adjoint_2d_10_25_25,
  &nfst_adjoint_2d_10_25_50,
  &nfst_adjoint_2d_25_10_25,
  &nfst_adjoint_2d_25_10_50,
  &nfst_adjoint_2d_25_25_25,
  &nfst_adjoint_2d_25_25_50,
};

static const trafo_delegate_t* trafos_adjoint_2d_direct_file[] = {&adjoint_direct};

void X(check_adjoint_2d_direct_file)(void)
{
  printf("check_adjoint_2d_direct_file:\n");
  check_many(SIZE(testcases_adjoint_2d_file), SIZE(initializers_direct), SIZE(trafos_adjoint_2d_direct_file),
    testcases_adjoint_2d_file, initializers_direct, &check_adjoint, trafos_adjoint_2d_direct_file);
}

static const trafo_delegate_t* trafos_adjoint_2d_fast_file[] = {&adjoint/*, &adjoint_2d*/};

void X(check_adjoint_2d_fast_file)(void)
{
  printf("check_adjoint_2d_fast_file:\n");
  check_many(SIZE(testcases_adjoint_2d_file), SIZE(initializers_2d), SIZE(trafos_adjoint_2d_fast_file),
    testcases_adjoint_2d_file, initializers_2d, &check_adjoint, trafos_adjoint_2d_fast_file);
}

static const testcase_delegate_online_t nfst_online_2d_50_50 = {setup_online, destroy_online, 2, 50 ,50};
static const testcase_delegate_online_t nfst_online_2d_100_50 = {setup_online, destroy_online, 2, 100 ,50};
static const testcase_delegate_online_t nfst_online_2d_200_50 = {setup_online, destroy_online, 2, 200 ,50};
static const testcase_delegate_online_t nfst_online_2d_500_50 = {setup_online, destroy_online, 2, 500 ,50};
static const testcase_delegate_online_t nfst_online_2d_1000_50 = {setup_online, destroy_online, 2, 1000 ,50};

static const testcase_delegate_online_t *testcases_2d_online[] =
{
  &nfst_online_2d_50_50,
  &nfst_online_2d_100_50,
  &nfst_online_2d_200_50,
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
  &nfst_online_2d_500_50,
  &nfst_online_2d_1000_50,
#endif
};

static const trafo_delegate_t* trafos_2d_online[] = {&trafo/*, &trafo_2d*/};

void X(check_2d_online)(void)
{
  check_many(SIZE(testcases_2d_online), SIZE(initializers_2d), SIZE(trafos_2d_online),
    testcases_2d_online, initializers_2d, &check_trafo, trafos_2d_online);
}

static const testcase_delegate_online_t nfst_adjoint_online_2d_50_50 = {setup_adjoint_online, destroy_online, 2, 50 ,50};
static const testcase_delegate_online_t nfst_adjoint_online_2d_100_50 = {setup_adjoint_online, destroy_online, 2, 100 ,50};
static const testcase_delegate_online_t nfst_adjoint_online_2d_200_50 = {setup_adjoint_online, destroy_online, 2, 200 ,50};
static const testcase_delegate_online_t nfst_adjoint_online_2d_500_50 = {setup_adjoint_online, destroy_online, 2, 500 ,50};
static const testcase_delegate_online_t nfst_adjoint_online_2d_1000_50 = {setup_adjoint_online, destroy_online, 2, 1000 ,50};

static const testcase_delegate_online_t *testcases_adjoint_2d_online[] =
{
  &nfst_adjoint_online_2d_50_50,
  &nfst_adjoint_online_2d_100_50,
  &nfst_adjoint_online_2d_200_50,
#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
  &nfst_adjoint_online_2d_500_50,
  &nfst_adjoint_online_2d_1000_50,
#endif
};

static const trafo_delegate_t* trafos_adjoint_2d_online[] = {&adjoint/*, &adjoint_2d*/};

void X(check_adjoint_2d_online)(void)
{
  check_many(SIZE(testcases_adjoint_2d_online), SIZE(initializers_2d), SIZE(trafos_adjoint_2d_online),
    testcases_adjoint_2d_online, initializers_2d, &check_adjoint, trafos_adjoint_2d_online);
}

/* 3D */

/* Initializers. */
static const init_delegate_t* initializers_3d[] =
{
  &init_3d,
  &init,
  &init_advanced_pre_psi,
  &init_advanced_pre_full_psi,
//  &init_advanced_pre_lin_psi,
#if defined(GAUSSIAN)
  &init_advanced_pre_fg_psi,
#endif
};

static const testcase_delegate_file_t nfst_3d_10_10_10_10 = {setup_file,destroy_file,ABSPATH("data/nfst_3d_10_10_10_10.txt")};

static const testcase_delegate_file_t *testcases_3d_file[] =
{
  &nfst_3d_10_10_10_10,
};

static const trafo_delegate_t* trafos_3d_direct_file[] = {&trafo_direct};

void X(check_3d_direct_file)(void)
{
  printf("check_3d_direct_file:\n");
  check_many(SIZE(testcases_3d_file), SIZE(initializers_direct), SIZE(trafos_3d_direct_file),
    testcases_3d_file, initializers_direct, &check_trafo, trafos_3d_direct_file);
}

static const trafo_delegate_t* trafos_3d_fast_file[] = {&trafo/*, &trafo_3d*/};

void X(check_3d_fast_file)(void)
{
  printf("check_3d_fast_file:\n");
  check_many(SIZE(testcases_3d_file), SIZE(initializers_3d), SIZE(trafos_3d_fast_file),
    testcases_3d_file, initializers_3d, &check_trafo, trafos_3d_fast_file);
}

static const testcase_delegate_file_t nfst_adjoint_3d_10_10_10_10 = {setup_file,destroy_file,ABSPATH("data/nfst_adjoint_3d_10_10_10_10.txt")};

static const testcase_delegate_file_t *testcases_adjoint_3d_file[] =
{
  &nfst_adjoint_3d_10_10_10_10,
};

static const trafo_delegate_t* trafos_adjoint_3d_direct_file[] = {&adjoint_direct};

void X(check_adjoint_3d_direct_file)(void)
{
  printf("check_adjoint_3d_direct_file:\n");
  check_many(SIZE(testcases_adjoint_3d_file), SIZE(initializers_direct), SIZE(trafos_adjoint_3d_direct_file),
    testcases_adjoint_3d_file, initializers_direct, &check_adjoint, trafos_adjoint_3d_direct_file);
}

static const trafo_delegate_t* trafos_adjoint_3d_fast_file[] = {&adjoint/*, &adjoint_3d*/};

void X(check_adjoint_3d_fast_file)(void)
{
  printf("check_adjoint_3d_fast_file:\n");
  check_many(SIZE(testcases_adjoint_3d_file), SIZE(initializers_3d), SIZE(trafos_adjoint_3d_fast_file),
    testcases_adjoint_3d_file, initializers_3d, &check_adjoint, trafos_adjoint_3d_fast_file);
}

#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
static const testcase_delegate_online_t nfst_online_3d_50_50 = {setup_online, destroy_online, 3, 50 ,50};

static const testcase_delegate_online_t *testcases_3d_online[] =
{
  &nfst_online_3d_50_50,
};

static const trafo_delegate_t* trafos_3d_online[] = {&trafo/*, &trafo_3d*/};

void X(check_3d_online)(void)
{
  check_many(SIZE(testcases_3d_online), SIZE(initializers_3d), SIZE(trafos_3d_online),
    testcases_3d_online, initializers_3d, &check_trafo, trafos_3d_online);
}

static const testcase_delegate_online_t nfst_adjoint_online_3d_50_50 = {setup_adjoint_online, destroy_online, 3, 50 ,50};

static const testcase_delegate_online_t *testcases_adjoint_3d_online[] =
{
  &nfst_adjoint_online_3d_50_50,
};

static const trafo_delegate_t* trafos_adjoint_3d_online[] = {&adjoint/*, &adjoint_3d*/};

void X(check_adjoint_3d_online)(void)
{
  check_many(SIZE(testcases_adjoint_3d_online), SIZE(initializers_3d), SIZE(trafos_adjoint_3d_online),
    testcases_adjoint_3d_online, initializers_3d, &check_adjoint, trafos_adjoint_3d_online);
}
#endif

/* 4D. */

/* Initializers. */
static const init_delegate_t* initializers_4d[] =
{
  &init,
  &init_advanced_pre_psi,
  &init_advanced_pre_full_psi,
//  &init_advanced_pre_lin_psi,
#if defined(GAUSSIAN)
  &init_advanced_pre_fg_psi,
#endif
};

#ifdef NFFT_EXHAUSTIVE_UNIT_TESTS
static const testcase_delegate_online_t nfst_online_4d_28_50 = {setup_online, destroy_online, 4, 28 ,50};

static const testcase_delegate_online_t *testcases_4d_online[] =
{
  &nfst_online_4d_28_50,
};

static const trafo_delegate_t* trafos_4d_online[] = {&trafo};

void X(check_4d_online)(void)
{
  check_many(SIZE(testcases_4d_online), SIZE(initializers_4d), SIZE(trafos_4d_online),
    testcases_4d_online, initializers_4d, &check_trafo, trafos_4d_online);
}

static const testcase_delegate_online_t nfst_adjoint_online_4d_28_50 = {setup_adjoint_online, destroy_online, 4, 28 ,50};

static const testcase_delegate_online_t *testcases_adjoint_4d_online[] =
{
  &nfst_adjoint_online_4d_28_50,
};

static const trafo_delegate_t* trafos_adjoint_4d_online[] = {&adjoint};

void X(check_adjoint_4d_online)(void)
{
  check_many(SIZE(testcases_adjoint_4d_online), SIZE(initializers_4d), SIZE(trafos_adjoint_4d_online),
    testcases_adjoint_4d_online, initializers_4d, &check_adjoint, trafos_adjoint_4d_online);
}
#endif
