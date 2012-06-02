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

#include "nfft3util.h"
#include "nfft3.h"
#include "infft.h"
#include "cycle.h"

typedef void (*init_t)(X(plan) *p, const int d, const int *N, const int M);

typedef struct init_delegate_s
{
  const char *name;
  init_t init;
} init_delegate_t;

static void init_1d_function(X(plan) *p, const int d, const int *N, const int M)
{X(init_1d)(p, N[0], M);}

static void init_2d_function(X(plan) *p, const int d, const int *N, const int M)
{X(init_2d)(p, N[0], N[1], M);}

static init_delegate_t init_1d = {"init_1d", init_1d_function};
static init_delegate_t init_2d = {"init_2d", init_2d_function};

typedef void (*trafo_t)(X(plan) *p);
typedef double (*cost_t)(X(plan) *p);
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

static double trafo_direct_cost_factor = 1.0E-6;

static double trafo_direct_cost(X(plan) *p)
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
          X(plan) p;
          int *N = malloc(d*sizeof(int)), i;
          for (i = 0; i < d; i++)
          {
            N[i] = Nd;
          }
          X(init)(&p, d, N, M);
          for (i = 0; i < M; i++)
            p.x[i] = K(0.0);
          if(p.nfft_flags & PRE_ONE_PSI)
            X(precompute_one_psi)(&p);
          for (i = 0; i < d*Nd; i++)
          {
            p.f_hat[i] = K(0.0) + K(0.0) * I;
          }
          {
            double r;
            ticks t0, t1;
            t0 = getticks();
            X(trafo_direct)(&p);
            t1 = getticks();
            r = X(elapsed_seconds)(t1, t0)/M;
            for (i = 0; i < d; i++)
              r = r / Nd;
            trafo_direct_cost_factor += r;
            printf("%E\n", r);
            x += 1;
          }
          X(finalize)(&p);
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

static R err_trafo_direct(X(plan) *p)
{
  return K(30.0) * EPSILON;
}

static R err_trafo(X(plan) *p)
{
  R m = ((R)p->m), s = K(0.0);
  int i;
  for (i = 0; i < p->d; i++)
    s = FMAX(s, ((R)p->sigma[i]));
#if defined(GAUSSIAN)
  R err = K(4.0) * EXP(-m*KPI*(K(1.0)-K(1.0)/(K(2.0)*s-K(1.0))));
#elif defined(B_SPLINE)
  R err = K(4.0) * POW(K(1.0)/(K(2.0)*s-K(1.0)),K(2.0)*m);
#elif defined(SINC_POWER)
  R err = (K(1.0)/(m-K(1.0))) * ((K(2.0)/(POW(s,K(2.0)*m))) + POW(s/(K(2.0)*s-K(1.0)),K(2.0)*m));
#elif defined(KAISER_BESSEL)
  R err = K(4.0) * KPI * (SQRT(m) + m) * SQRT(SQRT(K(1.0) - K(1.0)/s)) * EXP(-K2PI*m*SQRT(K(1.0)-K(1.0)/s));
#else
  #error Unsupported window function.
#endif

  return FMAX(K(30.0) * EPSILON, err);
}

static trafo_delegate_t trafo_direct = {"trafo_direct", X(trafo_direct), 0, trafo_direct_cost, err_trafo_direct};
static trafo_delegate_t trafo = {"trafo       ", X(trafo), X(check), 0, err_trafo};

#define MAX_SECONDS 0.1

static void pad(const char *s, const int n)
{
  const int m = n - strlen(s);
  int i;
  for (i = 0; i < m; i++)
    printf(" ");
}

static int check_nfft_xd(const char *filename, init_delegate_t *init_delegate, trafo_delegate_t *trafo_delegate)
{
  int result = EXIT_FAILURE;
  FILE *file = fopen(filename, "r");
  C *f = 0;
  X(plan) p;

  /* Check if input file could be opened. */
  if (!file)
    goto cleanup;

  printf("%-35s:", filename);

  {
    int d, j, *N, M;

    /* Dimensions. */
    fscanf(file, "%d", &d);
    /* Bandwidths. */
    N = malloc(d * sizeof(int));
    for (j = 0; j < d; j++)
      fscanf(file, "%d", &N[j]);
    /* Number of nodes. */
    fscanf(file, "%d", &M);

    printf(" d = %-1d, N = [", d);
    {
      for (j = 0; j < d; j++)
      {
        printf("%s%-4d", IF(j > 0,", ", ""), N[j]);
      }
    }
    printf(" M = %-4d,", M);
    printf(" %-10s,", init_delegate->name);
    printf(" %-10s", trafo_delegate->name);

    {
      int NN;

      /* Total number of Fourier coefficients. */
      for (j = 0, NN = 1; j < d; j++)
        NN *= N[j];

      /* Init plan. */
      //X(init_2d)(&p, N[0], N[1], M);
      init_delegate->init(&p, d, N, M);

      /* Nodes. */
      for (j = 0; j < M*d; j++)
      {
        fscanf(file, FFI, &p.x[j]);
      }

      /* Pre-compute Psi, maybe. */
      if(p.nfft_flags & PRE_ONE_PSI)
        X(precompute_one_psi)(&p);

      /* Fourier coefficients. */
      for (j = 0; j < NN; j++)
      {
        R re, im;
        fscanf(file, FFI " " FFI, &re, &im);
        p.f_hat[j] = re + im * I;
      }

      /* Reference function values. */
      {
        f = malloc(M * sizeof(C));

        for (j = 0; j < M; j++)
        {
          R re, im;
          fscanf(file, FFI " " FFI, &re, &im);
          f[j] = re + im * I;
        }

        if (trafo_delegate->check)
        {
          const char* check = trafo_delegate->check(&p);
          if (check != 0)
          {
            printf(" -> OK (");
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
            printf(" -> OK (cost too high)\n");
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
            result = IF(err < trafo_delegate->acc(&p), EXIT_SUCCESS, EXIT_FAILURE);
            printf(" -> %s " FE_ " (" FE_ ")\n", IF(result == EXIT_FAILURE, "FAIL", "OK"), err, bound);
          }
        }
      }
    }
  }

cleanup:
  if (file)
    fclose(file);
  if (f)
    free(f);
  X(finalize)(&p);
  return result;
}

static const char* filenames[] =
{
  "data/nfft_1d_1_1.txt",
  /*"data/nfft_1d_1_10.txt",
  "data/nfft_1d_1_20.txt",
  "data/nfft_1d_1_50.txt",
  "data/nfft_1d_2_1.txt",
  "data/nfft_1d_2_10.txt",
  "data/nfft_1d_2_20.txt",
  "data/nfft_1d_2_50.txt",
  "data/nfft_1d_4_1.txt",
  "data/nfft_1d_4_10.txt",
  "data/nfft_1d_4_20.txt",
  "data/nfft_1d_4_50.txt",
  "data/nfft_1d_10_1.txt",
  "data/nfft_1d_10_10.txt",
  "data/nfft_1d_10_20.txt",
  "data/nfft_1d_10_50.txt",
  "data/nfft_1d_20_1.txt",
  "data/nfft_1d_20_10.txt",
  "data/nfft_1d_20_20.txt",
  "data/nfft_1d_20_50.txt",
  "data/nfft_1d_50_1.txt",
  "data/nfft_1d_50_10.txt",
  "data/nfft_1d_50_20.txt",
  "data/nfft_1d_50_50.txt",*/
};

static const init_delegate_t* initializers[] = {&init_1d};
static const trafo_delegate_t* trafos[] = {&trafo_direct, &trafo};

int test(const int nf, const int ni, const int nt, const char **files,
  init_delegate_t **initializers, trafo_delegate_t **trafos)
{
  int i, j, k, result = EXIT_SUCCESS, r;
  for (i = 0; i < nf; i++)
    for (j = 0; j < ni; j++)
      for (k = 0; k < nt; k++)
      {
         r = check_nfft_xd(filenames[i], initializers[j], trafos[k]);
         result = IF(r == EXIT_FAILURE, EXIT_FAILURE, result);
      }
  return result;
}

#define SIZE(x) sizeof(x)/sizeof(x[0])

int main(void)
{
  return EXIT_SUCCESS;
  //test(SIZE(filenames), SIZE(initializers), SIZE(trafos), filenames,
  //  initializers, trafos);
}
