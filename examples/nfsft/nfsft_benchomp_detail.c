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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

#include "nfft3.h"
#include "infft.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void bench_openmp_readfile(FILE *infile, int *trafo_adjoint, int *N, int *M, double **x, C **f_hat, C **f)
{
  double re,im;
  int k, n, j, t;
  nfsft_plan plan;

  fscanf(infile, "%d %d %d", trafo_adjoint, N, M);
  *x = (double *)nfft_malloc(2*(*M)*sizeof(double));
  *f_hat = (C*)nfft_malloc((2*(*N)+2) * (2*(*N)+2) * sizeof(C));
  *f = (C*)nfft_malloc((*M)*sizeof(C));

  memset(*f_hat,0U,(2*(*N)+2) * (2*(*N)+2) * sizeof(C));
  memset(*f,0U,(*M)*sizeof(C));

#ifdef _OPENMP
  fftw_import_wisdom_from_filename("nfsft_benchomp_detail_threads.plan");
#else
  fftw_import_wisdom_from_filename("nfsft_benchomp_detail_single.plan");
#endif

  nfsft_init_guru(&plan, *N, *M, NFSFT_MALLOC_X | NFSFT_MALLOC_F |
    NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED | NFSFT_PRESERVE_F_HAT,
    PRE_PHI_HUT | FFTW_INIT | FFT_OUT_OF_PLACE, 6);

#ifdef _OPENMP
  fftw_export_wisdom_to_filename("nfsft_benchomp_detail_threads.plan");
#else
  fftw_export_wisdom_to_filename("nfsft_benchomp_detail_single.plan");
#endif

  for (j=0; j < *M; j++)
    for (t=0; t < 2; t++)
      fscanf(infile, "%lg", (*x)+2*j+t);

  if (trafo_adjoint==0)
  {
    for (k = 0; k <= *N; k++)
      for (n = -k; n <= k; n++)
      {
        fscanf(infile, "%lg %lg", &re, &im);
        (*f_hat)[NFSFT_INDEX(k,n,&plan)] = re + _Complex_I * im;
      }
  }
  else
  {
    for (j=0; j < *M; j++)
    {
      fscanf(infile, "%lg %lg", &re, &im);
      (*f)[j] = re + _Complex_I * im;
    }
  }

  nfsft_finalize(&plan);
}

void bench_openmp(int trafo_adjoint, int N, int M, double *x, C *f_hat, C *f, int m, int nfsft_flags, int psi_flags)
{
  nfsft_plan plan;
  int k, n;
//  int N, M, trafo_adjoint;
  int t, j;
  ticks t0, t1;
  double tt_total, tt_pre;

//  fscanf(infile, "%d %d %d", &trafo_adjoint, &N, &M);

/*#ifdef _OPENMP
  fftw_import_wisdom_from_filename("nfsft_benchomp_detail_threads.plan");
#else
  fftw_import_wisdom_from_filename("nfsft_benchomp_detail_single.plan");
#endif*/

  /* precomputation (for fast polynomial transform) */
//  nfsft_precompute(N,1000.0,0U,0U);

  /* Initialize transform plan using the guru interface. All input and output
   * arrays are allocated by nfsft_init_guru(). Computations are performed with
   * respect to L^2-normalized spherical harmonics Y_k^n. The array of spherical
   * Fourier coefficients is preserved during transformations. The NFFT uses a
   * cut-off parameter m = 6. See the NFFT 3 manual for details.
   */
  nfsft_init_guru(&plan, N, M, nfsft_flags | NFSFT_MALLOC_X | NFSFT_MALLOC_F |
    NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED | NFSFT_PRESERVE_F_HAT,
    PRE_PHI_HUT | psi_flags | FFTW_INIT | FFT_OUT_OF_PLACE, m);

/*#ifdef _OPENMP
  fftw_export_wisdom_to_filename("nfsft_benchomp_detail_threads.plan");
#else
  fftw_export_wisdom_to_filename("nfsft_benchomp_detail_single.plan");
#endif*/

  for (j=0; j < plan.M_total; j++)
  {
    for (t=0; t < 2; t++)
 //     fscanf(infile, "%lg", plan.x+2*j+t);
      plan.x[2*j+t] = x[2*j+t];
  }

  if (trafo_adjoint==0)
  {
    memset(plan.f_hat,0U,plan.N_total*sizeof(double _Complex));
    for (k = 0; k <= plan.N; k++)
      for (n = -k; n <= k; n++)
      {
//        fscanf(infile, "%lg %lg", &re, &im);
//        plan.f_hat[NFSFT_INDEX(k,n,&plan)] = re + _Complex_I * im;
        plan.f_hat[NFSFT_INDEX(k,n,&plan)] = f_hat[NFSFT_INDEX(k,n,&plan)];
      }
  }
  else
  {
    for (j=0; j < plan.M_total; j++)
    {
//      fscanf(infile, "%lg %lg", &re, &im);
//      plan.f[j] = re + _Complex_I * im;
      plan.f[j] = f[j];
    }

    memset(plan.f_hat,0U,plan.N_total*sizeof(double _Complex));
  }

  t0 = getticks();
  /* precomputation (for NFFT, node-dependent) */
  nfsft_precompute_x(&plan);
  t1 = getticks();
  tt_pre = nfft_elapsed_seconds(t1,t0);

  if (trafo_adjoint==0)
    nfsft_trafo(&plan);
  else
    nfsft_adjoint(&plan);
  t1 = getticks();
  tt_total = nfft_elapsed_seconds(t1,t0);

#ifndef MEASURE_TIME
  plan.MEASURE_TIME_t[0] = 0.0;
  plan.MEASURE_TIME_t[2] = 0.0;
#endif

#ifndef MEASURE_TIME_FFTW
  plan.MEASURE_TIME_t[1] = 0.0;
#endif

  printf("%.6e %.6e %6e %.6e %.6e %.6e\n", tt_pre, plan.MEASURE_TIME_t[0], plan.MEASURE_TIME_t[1], plan.MEASURE_TIME_t[2], tt_total-tt_pre-plan.MEASURE_TIME_t[0]-plan.MEASURE_TIME_t[1]-plan.MEASURE_TIME_t[2], tt_total);

  /** finalise the one dimensional plan */
  nfsft_finalize(&plan);
}

int main(int argc, char **argv)
{
  int m, nfsft_flags, psi_flags;
  int nrepeat;
  int trafo_adjoint, N, M, r;
  double *x;
  C *f_hat, *f;
#ifdef _OPENMP
  int nthreads;

  if (argc != 6)
    return 1;

  nthreads = atoi(argv[5]);
  fftw_init_threads();
  omp_set_num_threads(nthreads);
#else
  if (argc != 5)
    return 1;
#endif

  m = atoi(argv[1]);
  nfsft_flags = atoi(argv[2]);
  psi_flags = atoi(argv[3]);
  nrepeat = atoi(argv[4]);

  bench_openmp_readfile(stdin, &trafo_adjoint, &N, &M, &x, &f_hat, &f);

  /* precomputation (for fast polynomial transform) */
  nfsft_precompute(N,1000.0,0U,0U);

  for (r = 0; r < nrepeat; r++)
    bench_openmp(trafo_adjoint, N, M, x, f_hat, f, m, nfsft_flags, psi_flags);

  return 0;
}
