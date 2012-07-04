/*
 * Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* $Id: simple_test.c 3509 2010-05-25 19:00:59Z keiner $ */

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <CUnit/CUnit.h>

#include "config.h"
#include "nfft3util.h"
#include "nfft3.h"
#include "infft.h"
#include "cycle.h"
#include "check_nfft.h"

static double trafo_direct_cost_factor = 1.0E-6;

double X(trafo_direct_cost)(X(plan) *p)
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
          int *N = malloc(d*sizeof(int)), i;
          for (i = 0; i < d; i++)
          {
            N[i] = Nd;
          }
          X(init)(&p2, d, N, M);
          for (i = 0; i < M; i++)
            p2.x[i] = K(0.0);
          if(p2.nfft_flags & PRE_ONE_PSI)
            X(precompute_one_psi)(&p2);
          for (i = 0; i < d*Nd; i++)
          {
            p2.f_hat[i] = K(0.0) + K(0.0) * I;
          }
          {
            double r;
            ticks t0, t1;
            t0 = getticks();
            X(trafo_direct)(&p2);
            t1 = getticks();
            r = X(elapsed_seconds)(t1, t0)/M;
            for (i = 0; i < d; i++)
              r = r / Nd;
            trafo_direct_cost_factor += r;
            printf("%E\n", r);
            x += 1;
          }
          X(finalize)(&p2);
          free(N);
        }
      }
    }
    trafo_direct_cost_factor = trafo_direct_cost_factor/((double)x);
    printf("--> %E\n", trafo_direct_cost_factor);
  }

  {
    int c = p->M_total, i;

    for (i = 0; i < p->d; i++)
      c *= p->N[i];

    return trafo_direct_cost_factor * c;
  }
}

R X(err_trafo_direct)(X(plan) *p)
{
  UNUSED(p);
  return K(30.0) * EPSILON;
}

R X(err_trafo)(X(plan) *p)
{
  if (p->nfft_flags & PRE_LIN_PSI)
    return K(3.5)*K(10E-09);
  {
    const R m = ((R)p->m);
    R s, err;
    int i;
    for (i = 0, s = ((R)p->sigma[0]); i < p->d; i++)
      s = FMIN(s, ((R)p->sigma[i]));
#if defined(GAUSSIAN)
    err = K(4.0) * EXP(-m*KPI*(K(1.0)-K(1.0)/(K(2.0)*s-K(1.0))));
#elif defined(B_SPLINE)
    //printf("m = %E, s = %E, a1 = %E, a2 = %E, z = %E\n", m, s, K(1.0)/(K(2.0)*s-K(1.0)), K(2.0)*m, K(4.0) * POW(K(1.0)/(K(2.0)*s-K(1.0)),K(2.0)*m));
    //printf("<s = %E>", s);
    //fflush(stdout);
    err = K(4.0) * POW(K(1.0)/(K(2.0)*s-K(1.0)),K(2.0)*m);
  #elif defined(SINC_POWER)
    err = (K(1.0)/(m-K(1.0))) * ((K(2.0)/(POW(s,K(2.0)*m))) + POW(s/(K(2.0)*s-K(1.0)),K(2.0)*m));
  #elif defined(KAISER_BESSEL)
    if (p->nfft_flags | PRE_LIN_PSI)
    {
      R K = ((R)p->K);
      err = EXP(K2PI * m)/(K(8.0) * K * K);
    }
    else
      err = K(4.0) * KPI * (SQRT(m) + m) * SQRT(SQRT(K(1.0) - K(1.0)/s)) * EXP(-K2PI*m*SQRT(K(1.0)-K(1.0)/s));
  #else
    #error Unsupported window function.
  #endif

    return FMAX(K(70.0) * EPSILON, err);
  }
}

#define MAX_SECONDS 0.1

int X(check_single)(const testcase_delegate_t *testcase,
    init_delegate_t *init_delegate, trafo_delegate_t *trafo_delegate)
{
  int result = EXIT_FAILURE;
  X(plan) p;
  int d, j, *N, NN, M;
  R *x;
  C *f_hat, *f;

  testcase->setup(testcase, &d, &N, &NN, &M, &x, &f_hat, &f);

  /* Init plan. */
  printf(", %-24s", init_delegate->name);
  init_delegate->init(init_delegate, &p, d, N, M);

  /* Nodes. */
  for (j = 0; j < M*d; j++)
  {
    p.x[j] = x[j];
  }

  /* Pre-compute Psi, maybe. */
  if(p.nfft_flags & PRE_ONE_PSI)
    X(precompute_one_psi)(&p);

  /* Fourier coefficients. */
  for (j = 0; j < NN; j++)
  {
    p.f_hat[j] = f_hat[j];
  }

  printf(", %-12s", trafo_delegate->name);

  if (trafo_delegate->check)
  {
    const char* check = trafo_delegate->check(&p);
    if (check != 0)
    {
      printf(" -> %-4s (","OK");
      printf("%s", check);
      printf(")\n");
      result = EXIT_SUCCESS;
      goto cleanup;
    }
  }
  else if (trafo_delegate->cost)
  {
    const double cost = trafo_delegate->cost(&p);
    if (cost > MAX_SECONDS)
    {
      printf(" -> %-4s (cost too high)\n","OK");
      result = EXIT_SUCCESS;
      goto cleanup;
    }
  }

  trafo_delegate->trafo(&p);

  /* debug */
  /*for (j = 0; j < M; j++)
    fprintf(stderr, "f[%2d] = " FE_ " + " FE_ "I, f[%2d] = " FE_ " + " FE_ "I, err = " FE_ "\n", j,
      CREAL(f[j]), CIMAG(f[j]), j, CREAL(p.f[j]), CIMAG(p.f[j]), CABS(f[j] - p.f[j]) / CABS(f[j]));*/

  /* Standard NFFT error measure. */
  {
    R numerator = K(0.0), denominator = K(0.0);
    for (j = 0; j < M; j++)
      numerator = MAX(numerator, CABS(f[j] - p.f[j]));
    for (j = 0; j < NN; j++)
      denominator += CABS(p.f_hat[j]);
    {
      R err = numerator/denominator;
      R bound = trafo_delegate->acc(&p);
      result = IF(err < bound, EXIT_SUCCESS, EXIT_FAILURE);
      printf(" -> %-4s " FE_ " (" FE_ ")\n", IF(result == EXIT_FAILURE, "FAIL", "OK"), err, bound);
    }
  }

cleanup:
  testcase->destroy(testcase, x, f_hat, f);
  X(finalize)(&p);

  CU_ASSERT(result == EXIT_SUCCESS);
  return result;
}

void X(check_many)(const size_t nf, const size_t ni, const size_t nt,
  const testcase_delegate_t **testcases, init_delegate_t **initializers,
  trafo_delegate_t **trafos)
{
  size_t i, j, k;
  int result = EXIT_SUCCESS, r;
  for (k = 0; k < nt; k++)
  {
    for (i = 0; i < nf; i++)
    {
      for (j = 0; j < ni; j++)
      {
         r = X(check_single)(testcases[i], initializers[j], trafos[k]);
         result = IF(r == EXIT_FAILURE, EXIT_FAILURE, result);
      }
    }
  }
}

void X(setup_file)(const testcase_delegate_t *ego_, int *d, int **N, int *NN, int *M, R **x, C **f_hat, C **f)
{
  const testcase_delegate_file_t *ego = (const testcase_delegate_file_t*)ego_;
  int j;
  FILE *file = fopen(ego->filename, "r");

  printf("%-28s", ego->filename);

  /* Dimensions. */
  fscanf(file, "%d", d);
  /* Bandwidths. */
  *N = malloc(*d * sizeof(int));
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
    *NN *= (*N)[j];

  /* Nodes. */
  *x = malloc(M[0]*d[0]*sizeof(R));
  for (j = 0; j < M[0]*d[0]; j++)
  {
    fscanf(file, FFI, &((*x)[j]));
  }

  /* Fourier coefficients. */
  *f_hat = malloc(NN[0]*sizeof(C));
  for (j = 0; j < NN[0]; j++)
  {
    R re, im;
    fscanf(file, FFI " " FFI, &re, &im);
    (*f_hat)[j] = re + im * I;
  }

  /* Reference function values. */
  *f = malloc(M[0] * sizeof(C));
  for (j = 0; j < M[0]; j++)
  {
    R re, im;
    fscanf(file, FFI " " FFI, &re, &im);
    (*f)[j] = re + im * I;
  }

  fclose(file);
}

void X(destroy_file)(const testcase_delegate_t *ego_, R *x, C *f_hat, C *f)
{
  UNUSED(ego_);
  free(x);
  free(f_hat);
  free(f);
}

void X(setup_online)(const testcase_delegate_t *ego_, int *d, int **N, int *NN, int *M, R **x, C **f_hat, C **f)
{
  const testcase_delegate_online_t *ego = (const testcase_delegate_online_t*)ego_;
  int j;

  /* Dimensions. */
  *d = ego->d;
  /* Bandwidths. */
  *N = malloc(*d * sizeof(int));
  for (j = 0; j < *d; j++)
    (*N)[j] = ego->N;
  /* Number of nodes. */
  *M = ego->M;

  printf("%-28s", "online");

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
    *NN *= (*N)[j];

  /* Nodes. */
  *x = malloc(M[0]*d[0]*sizeof(R));
  for (j = 0; j < M[0]*d[0]; j++)
  {
    (*x)[j] = nfft_drand48() - K(0.5);
  }

  /* Fourier coefficients. */
  *f_hat = malloc(NN[0]*sizeof(C));
  for (j = 0; j < NN[0]; j++)
  {
    (*f_hat)[j] = (nfft_drand48() - K(0.5)) + (nfft_drand48() - K(0.5)) * I;
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
    if(p.nfft_flags & PRE_ONE_PSI)
      X(precompute_one_psi)(&p);

    /* Fourier coefficients. */
    for (j = 0; j < *NN; j++)
    {
      p.f_hat[j] = (*f_hat)[j];
    }

    X(trafo_direct)(&p);

    /* Reference function values. */
    *f = malloc(M[0] * sizeof(C));
    for (j = 0; j < M[0]; j++)
    {
      (*f)[j] = p.f[j];
    }

    X(finalize)(&p);
  }
}

void X(destroy_online)(const testcase_delegate_t *ego_, R *x, C *f_hat, C *f)
{
  UNUSED(ego_);
  free(x);
  free(f_hat);
  free(f);
}

/* Initializers. */
void X(init_1d_)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  UNUSED(ego);
  UNUSED(d);
  X(init_1d)(p, N[0], M);
}

void X(init_2d_)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  UNUSED(ego);
  UNUSED(d);
  X(init_2d)(p, N[0], N[1], M);
}

void X(init_3d_)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  UNUSED(ego);
  UNUSED(d);
  X(init_3d)(p, N[0], N[1], N[2], M);
}

void X(init_)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  UNUSED(ego);
  X(init)(p, d, N, M);
}

void X(init_advanced_pre_psi_)(init_delegate_t *ego, X(plan) *p, const int d, const int *N, const int M)
{
  int *n = malloc(d*sizeof(int));
  int i;
  for (i = 0; i < d; i++)
    n[i] = 2*X(next_power_of_2)(N[i]);
  X(init_guru)(p, d, N, M, n, ego->m, ego->nfft_flags, ego->fftw_flags);
  free(n);
}

init_delegate_t init_1d = {"init_1d", X(init_1d_), 0, 0, 0};
init_delegate_t init_2d = {"init_2d", X(init_2d_), 0, 0, 0};
init_delegate_t init_3d = {"init_3d", X(init_3d_), 0, 0, 0};
init_delegate_t init = {"init", X(init_), 0, 0, 0};
init_delegate_t init_advanced_pre_psi = {"init_guru (PRE PSI)", X(init_advanced_pre_psi_), WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | PRE_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
init_delegate_t init_advanced_pre_full_psi = {"init_guru (PRE FULL PSI)", X(init_advanced_pre_psi_), WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | PRE_FULL_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
init_delegate_t init_advanced_pre_lin_psi = {"init_guru (PRE LIN PSI)", X(init_advanced_pre_psi_), WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | PRE_LIN_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
#if defined(GAUSSIAN)
init_delegate_t init_advanced_pre_fg_psi = {"init_guru (PRE FG PSI)", X(init_advanced_pre_psi_), WINDOW_HELP_ESTIMATE_m, PRE_PHI_HUT | FG_PSI | PRE_FG_PSI | DEFAULT_NFFT_FLAGS, DEFAULT_FFTW_FLAGS};
#endif

trafo_delegate_t trafo_direct = {"trafo_direct", X(trafo_direct), 0, X(trafo_direct_cost), X(err_trafo_direct)};
trafo_delegate_t trafo = {"trafo", X(trafo), X(check), 0, X(err_trafo)};
trafo_delegate_t trafo_1d = {"trafo_1d", X(trafo_1d), X(check), 0, X(err_trafo)};
trafo_delegate_t trafo_2d = {"trafo_2d", X(trafo_2d), X(check), 0, X(err_trafo)};
trafo_delegate_t trafo_3d = {"trafo_3d", X(trafo_3d), X(check), 0, X(err_trafo)};

/* 1D */

/* Initializers. */
static const init_delegate_t* initializers_1d[] =
{
  &init_1d,
  &init,
  &init_advanced_pre_psi,
  &init_advanced_pre_full_psi,
  &init_advanced_pre_lin_psi,
#if defined(GAUSSIAN)
  &init_advanced_pre_fg_psi,
#endif
};

static const testcase_delegate_file_t nfft_1d_1_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_1.txt"};
static const testcase_delegate_file_t nfft_1d_1_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_10.txt"};
static const testcase_delegate_file_t nfft_1d_1_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_20.txt"};
static const testcase_delegate_file_t nfft_1d_1_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_1_50.txt"};
static const testcase_delegate_file_t nfft_1d_2_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_1.txt"};
static const testcase_delegate_file_t nfft_1d_2_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_10.txt"};
static const testcase_delegate_file_t nfft_1d_2_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_20.txt"};
static const testcase_delegate_file_t nfft_1d_2_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_2_50.txt"};
static const testcase_delegate_file_t nfft_1d_4_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_1.txt"};
static const testcase_delegate_file_t nfft_1d_4_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_10.txt"};
static const testcase_delegate_file_t nfft_1d_4_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_20.txt"};
static const testcase_delegate_file_t nfft_1d_4_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_4_50.txt"};
static const testcase_delegate_file_t nfft_1d_10_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_1.txt"};
static const testcase_delegate_file_t nfft_1d_10_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_10.txt"};
static const testcase_delegate_file_t nfft_1d_10_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_20.txt"};
static const testcase_delegate_file_t nfft_1d_10_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_10_50.txt"};
static const testcase_delegate_file_t nfft_1d_20_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_1.txt"};
static const testcase_delegate_file_t nfft_1d_20_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_10.txt"};
static const testcase_delegate_file_t nfft_1d_20_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_20.txt"};
static const testcase_delegate_file_t nfft_1d_20_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_20_50.txt"};
static const testcase_delegate_file_t nfft_1d_50_1 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_1.txt"};
static const testcase_delegate_file_t nfft_1d_50_10 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_10.txt"};
static const testcase_delegate_file_t nfft_1d_50_20 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_20.txt"};
static const testcase_delegate_file_t nfft_1d_50_50 = {X(setup_file), X(destroy_file), "data/nfft_1d_50_50.txt"};

static const testcase_delegate_file_t *testcases_1d_file[] =
{
  &nfft_1d_1_1, &nfft_1d_1_10, &nfft_1d_1_20, &nfft_1d_1_50,
  &nfft_1d_2_1, &nfft_1d_2_10, &nfft_1d_2_20, &nfft_1d_2_50,
  &nfft_1d_4_1, &nfft_1d_4_10, &nfft_1d_4_20, &nfft_1d_4_50,
  &nfft_1d_10_1, &nfft_1d_10_10, &nfft_1d_10_20, &nfft_1d_10_50,
  &nfft_1d_20_1, &nfft_1d_20_10, &nfft_1d_20_20, &nfft_1d_20_50,
  &nfft_1d_50_1, &nfft_1d_50_10, &nfft_1d_50_20, &nfft_1d_50_50,
};

static const trafo_delegate_t* trafos_1d_file[] = {&trafo_direct, &trafo, &trafo_1d};

static void check_nfft_1d_file(void)
{
  X(check_many)(SIZE(testcases_1d_file), SIZE(initializers_1d), SIZE(trafos_1d_file),
    testcases_1d_file, initializers_1d, trafos_1d_file);
}

static const testcase_delegate_online_t nfft_online_1d_50_50 = {X(setup_online), X(destroy_online), 1, 50 ,50};
static const testcase_delegate_online_t nfft_online_1d_100_50 = {X(setup_online), X(destroy_online), 1, 100 ,50};
static const testcase_delegate_online_t nfft_online_1d_200_50 = {X(setup_online), X(destroy_online), 1, 200 ,50};
static const testcase_delegate_online_t nfft_online_1d_500_50 = {X(setup_online), X(destroy_online), 1, 500 ,50};
static const testcase_delegate_online_t nfft_online_1d_1000_50 = {X(setup_online), X(destroy_online), 1, 1000 ,50};
static const testcase_delegate_online_t nfft_online_1d_2000_50 = {X(setup_online), X(destroy_online), 1, 2000 ,50};
static const testcase_delegate_online_t nfft_online_1d_5000_50 = {X(setup_online), X(destroy_online), 1, 5000 ,50};
static const testcase_delegate_online_t nfft_online_1d_10000_50 = {X(setup_online), X(destroy_online), 1, 10000 ,50};
/*static const testcase_delegate_online_t nfft_online_1d_100000_50 = {X(setup_online), X(destroy_online), 1, 100000 ,50};
static const testcase_delegate_online_t nfft_online_1d_1000000_50 = {X(setup_online), X(destroy_online), 1, 1000000 ,50};*/

static const testcase_delegate_online_t *testcases_1d_online[] =
{
  &nfft_online_1d_50_50,
  &nfft_online_1d_100_50,
  &nfft_online_1d_200_50,
  &nfft_online_1d_500_50,
  &nfft_online_1d_1000_50,
  &nfft_online_1d_2000_50,
  &nfft_online_1d_5000_50,
  &nfft_online_1d_10000_50,
/*  &nfft_online_1d_100000_50,
  &nfft_online_1d_1000000_50,*/
};

static const trafo_delegate_t* trafos_1d_online[] = {&trafo, &trafo_1d};

static void check_nfft_1d_online(void)
{
  X(check_many)(SIZE(testcases_1d_online), SIZE(initializers_1d), SIZE(trafos_1d_online),
    testcases_1d_online, initializers_1d, trafos_1d_online);
}

/* 2D */

/* Initializers. */
static const init_delegate_t* initializers_2d[] =
{
  &init_2d,
  &init,
  &init_advanced_pre_psi,
  &init_advanced_pre_full_psi,
  &init_advanced_pre_lin_psi,
#if defined(GAUSSIAN)
  &init_advanced_pre_fg_psi,
#endif
};

static const testcase_delegate_file_t nfft_2d_10_10_20 = {X(setup_file),X(destroy_file),"data/nfft_2d_10_10_20.txt"};
static const testcase_delegate_file_t nfft_2d_10_10_50 = {X(setup_file),X(destroy_file),"data/nfft_2d_10_10_50.txt"};
static const testcase_delegate_file_t nfft_2d_10_20_20 = {X(setup_file),X(destroy_file),"data/nfft_2d_10_20_20.txt"};
static const testcase_delegate_file_t nfft_2d_10_20_50 = {X(setup_file),X(destroy_file),"data/nfft_2d_10_20_50.txt"};
static const testcase_delegate_file_t nfft_2d_20_10_20 = {X(setup_file),X(destroy_file),"data/nfft_2d_20_10_20.txt"};
static const testcase_delegate_file_t nfft_2d_20_10_50 = {X(setup_file),X(destroy_file),"data/nfft_2d_20_10_50.txt"};
static const testcase_delegate_file_t nfft_2d_20_20_20 = {X(setup_file),X(destroy_file),"data/nfft_2d_20_20_20.txt"};
static const testcase_delegate_file_t nfft_2d_20_20_50 = {X(setup_file),X(destroy_file),"data/nfft_2d_20_20_50.txt"};

static const testcase_delegate_file_t *testcases_2d_file[] =
{
  &nfft_2d_10_10_20,
  &nfft_2d_10_10_50,
  &nfft_2d_10_20_20,
  &nfft_2d_10_20_50,
  &nfft_2d_20_10_20,
  &nfft_2d_20_10_50,
  &nfft_2d_20_20_20,
  &nfft_2d_20_20_50,
};

static const trafo_delegate_t* trafos_2d_file[] = {&trafo_direct, &trafo, &trafo_2d};

static void check_nfft_2d_file(void)
{
  X(check_many)(SIZE(testcases_2d_file), SIZE(initializers_2d), SIZE(trafos_2d_file),
    testcases_2d_file, initializers_2d, trafos_2d_file);
}

static const testcase_delegate_online_t nfft_online_2d_50_50 = {X(setup_online), X(destroy_online), 2, 50 ,50};
static const testcase_delegate_online_t nfft_online_2d_100_50 = {X(setup_online), X(destroy_online), 2, 100 ,50};
static const testcase_delegate_online_t nfft_online_2d_200_50 = {X(setup_online), X(destroy_online), 2, 200 ,50};
static const testcase_delegate_online_t nfft_online_2d_500_50 = {X(setup_online), X(destroy_online), 2, 500 ,50};
static const testcase_delegate_online_t nfft_online_2d_1000_50 = {X(setup_online), X(destroy_online), 2, 1000 ,50};
static const testcase_delegate_online_t nfft_online_2d_2000_50 = {X(setup_online), X(destroy_online), 2, 2000 ,50};

static const testcase_delegate_online_t *testcases_2d_online[] =
{
  &nfft_online_2d_50_50,
  &nfft_online_2d_100_50,
  &nfft_online_2d_200_50,
  &nfft_online_2d_500_50,
  &nfft_online_2d_1000_50,
  &nfft_online_2d_2000_50,
};

static const trafo_delegate_t* trafos_2d_online[] = {&trafo, &trafo_2d};

static void check_nfft_2d_online(void)
{
  X(check_many)(SIZE(testcases_2d_online), SIZE(initializers_2d), SIZE(trafos_2d_online),
    testcases_2d_online, initializers_2d, trafos_2d_online);
}

/* 3D */

/* Initializers. */
static const init_delegate_t* initializers_3d[] =
{
  &init_3d,
  &init,
  &init_advanced_pre_psi,
  &init_advanced_pre_full_psi,
  &init_advanced_pre_lin_psi,
#if defined(GAUSSIAN)
  &init_advanced_pre_fg_psi,
#endif
};

static const testcase_delegate_file_t nfft_3d_10_10_10_10 = {X(setup_file),X(destroy_file),"data/nfft_3d_10_10_10_10.txt"};

static const testcase_delegate_file_t *testcases_3d_file[] =
{
  &nfft_3d_10_10_10_10,
};

static const trafo_delegate_t* trafos_3d_file[] = {&trafo_direct, &trafo, &trafo_3d};

static void check_nfft_3d_file(void)
{
  X(check_many)(SIZE(testcases_3d_file), SIZE(initializers_3d), SIZE(trafos_3d_file),
    testcases_3d_file, initializers_3d, trafos_3d_file);
}

static const testcase_delegate_online_t nfft_online_3d_50_50 = {X(setup_online), X(destroy_online), 3, 50 ,50};
static const testcase_delegate_online_t nfft_online_3d_100_50 = {X(setup_online), X(destroy_online), 3, 100 ,50};

static const testcase_delegate_online_t *testcases_3d_online[] =
{
  &nfft_online_3d_50_50,
  &nfft_online_3d_100_50,
};

static const trafo_delegate_t* trafos_3d_online[] = {&trafo, &trafo_3d};

static void check_nfft_3d_online(void)
{
  X(check_many)(SIZE(testcases_3d_online), SIZE(initializers_3d), SIZE(trafos_3d_online),
    testcases_3d_online, initializers_3d, trafos_3d_online);
}

/* 4D. */

/* Initializers. */
static const init_delegate_t* initializers_4d[] =
{
  &init,
  &init_advanced_pre_psi,
  &init_advanced_pre_full_psi,
  &init_advanced_pre_lin_psi,
#if defined(GAUSSIAN)
  &init_advanced_pre_fg_psi,
#endif
};

static const testcase_delegate_online_t nfft_online_4d_35_50 = {X(setup_online), X(destroy_online), 4, 35 ,50};

static const testcase_delegate_online_t *testcases_4d_online[] =
{
  &nfft_online_4d_35_50,
};

static const trafo_delegate_t* trafos_4d_online[] = {&trafo};

static void check_nfft_4d_online(void)
{
  X(check_many)(SIZE(testcases_4d_online), SIZE(initializers_4d), SIZE(trafos_4d_online),
    testcases_4d_online, initializers_4d, trafos_4d_online);
}

int main(void)
{
  CU_pSuite s;
  CU_initialize_registry();
  CU_set_output_filename("nfft");
  s = CU_add_suite("nfft", 0, 0);
  /*CU_add_test(s, "nfft_1d_file", check_nfft_1d_file);*/
  CU_add_test(s, "nfft_1d_online", check_nfft_1d_online);
  /*CU_add_test(s, "nfft_2d_file", check_nfft_2d_file);*/
  /*CU_add_test(s, "nfft_2d_online", check_nfft_2d_online);*/
  /*CU_add_test(s, "nfft_3d_file", check_nfft_3d_file);*/
  /*CU_add_test(s, "nfft_3d_online", check_nfft_3d_online);*/
  /*CU_add_test(s, "nfft_4d_online", check_nfft_4d_online);*/
  CU_automated_run_tests();
  /*CU_basic_run_tests();*/
  CU_cleanup_registry();
  return 0;
}
