/*
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* $Id: simple_test.c 3372 2009-10-21 06:04:05Z skunis $ */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

#include "nfft3util.h"
#include "nfft3.h"
#include "infft.h"
#include "fastsum.h"
#include "kernels.h"
#ifdef _OPENMP
#include <omp.h>
#endif

int bench_openmp(FILE *infile, int n, int m, int p, double _Complex (*kernel)(double , int , const double *), double c, double eps_I, double eps_B)
{
  fastsum_plan my_fastsum_plan;
  int d, L, M;
  int t, j;
  double re,im;
  double r_max = 0.25 - my_fastsum_plan.eps_B/2.0;
  ticks t0, t1;
  double tt_total;

  fscanf(infile, "%d %d %d", &d, &L, &M);

#ifdef _OPENMP
  fftw_import_wisdom_from_filename("fastsum_benchomp_detail_threads.plan");
#else
  fftw_import_wisdom_from_filename("fastsum_benchomp_detail_single.plan");
#endif

  fastsum_init_guru(&my_fastsum_plan, d, L, M, kernel, &c, 0, n, m, p, eps_I, eps_B);

#ifdef _OPENMP
  fftw_export_wisdom_to_filename("fastsum_benchomp_detail_threads.plan");
#else
  fftw_export_wisdom_to_filename("fastsum_benchomp_detail_single.plan");
#endif

  for (j=0; j < L; j++)
  {
    for (t=0; t < d; t++)
    {
      double v;
      fscanf(infile, "%lg", &v);
      my_fastsum_plan.x[d*j+t] = v * r_max;
    }
  }

  for (j=0; j < L; j++)
  {
    fscanf(infile, "%lg %lg", &re, &im);
    my_fastsum_plan.alpha[j] = re + _Complex_I * im;
  }

  for (j=0; j < M; j++)
  {
    for (t=0; t < d; t++)
    {
      double v;
      fscanf(infile, "%lg", &v);
      my_fastsum_plan.y[d*j+t] = v * r_max;
    }
  }

  /** precomputation */
  t0 = getticks();
  fastsum_precompute(&my_fastsum_plan);

  /** fast computation */
  fastsum_trafo(&my_fastsum_plan);
  t1 = getticks();
  tt_total = nfft_elapsed_seconds(t1,t0);

#ifndef MEASURE_TIME
  my_fastsum_plan.MEASURE_TIME_t[0] = 0.0;
  my_fastsum_plan.MEASURE_TIME_t[1] = 0.0;
  my_fastsum_plan.MEASURE_TIME_t[2] = 0.0;
  my_fastsum_plan.MEASURE_TIME_t[3] = 0.0;
  my_fastsum_plan.MEASURE_TIME_t[4] = 0.0;
  my_fastsum_plan.MEASURE_TIME_t[5] = 0.0;
  my_fastsum_plan.MEASURE_TIME_t[6] = 0.0;
  my_fastsum_plan.MEASURE_TIME_t[7] = 0.0;
  my_fastsum_plan.mv1.MEASURE_TIME_t[0] = 0.0;
  my_fastsum_plan.mv1.MEASURE_TIME_t[2] = 0.0;
  my_fastsum_plan.mv2.MEASURE_TIME_t[0] = 0.0;
  my_fastsum_plan.mv2.MEASURE_TIME_t[2] = 0.0;
#endif
#ifndef MEASURE_TIME_FFTW
  my_fastsum_plan.mv1.MEASURE_TIME_t[1] = 0.0;
  my_fastsum_plan.mv2.MEASURE_TIME_t[1] = 0.0;
#endif

  printf("%.6e %.6e %.6e %6e %.6e %.6e %.6e %.6e %.6e %6e %.6e %.6e %6e %.6e %.6e %6e\n", my_fastsum_plan.MEASURE_TIME_t[0], my_fastsum_plan.MEASURE_TIME_t[1], my_fastsum_plan.MEASURE_TIME_t[2], my_fastsum_plan.MEASURE_TIME_t[3], my_fastsum_plan.MEASURE_TIME_t[4], my_fastsum_plan.MEASURE_TIME_t[5], my_fastsum_plan.MEASURE_TIME_t[6], my_fastsum_plan.MEASURE_TIME_t[7], tt_total-my_fastsum_plan.MEASURE_TIME_t[0]-my_fastsum_plan.MEASURE_TIME_t[1]-my_fastsum_plan.MEASURE_TIME_t[2]-my_fastsum_plan.MEASURE_TIME_t[3]-my_fastsum_plan.MEASURE_TIME_t[4]-my_fastsum_plan.MEASURE_TIME_t[5]-my_fastsum_plan.MEASURE_TIME_t[6]-my_fastsum_plan.MEASURE_TIME_t[7], tt_total, my_fastsum_plan.mv1.MEASURE_TIME_t[0], my_fastsum_plan.mv1.MEASURE_TIME_t[1], my_fastsum_plan.mv1.MEASURE_TIME_t[2], my_fastsum_plan.mv2.MEASURE_TIME_t[0], my_fastsum_plan.mv2.MEASURE_TIME_t[1], my_fastsum_plan.mv2.MEASURE_TIME_t[2]);

  fastsum_finalize(&my_fastsum_plan);

  return 0;
}

int main(int argc, char **argv)
{
  int n;                                             /**< expansion degree        */
  int m;                                             /**< cut-off parameter       */
  int p;                                             /**< degree of smoothness    */
  char *s;                                           /**< name of kernel          */
  double _Complex (*kernel)(double , int , const double *);  /**< kernel function         */
  double c;                                          /**< parameter for kernel    */
  double eps_I;                                      /**< inner boundary          */
  double eps_B;                                      /**< outer boundary          */


#ifdef _OPENMP
  int nthreads;

  if (argc != 9)
    return 1;

  nthreads = atoi(argv[8]);
  fftw_init_threads();
  omp_set_num_threads(nthreads);
#else
  if (argc != 8)
    return 1;
#endif

  n=atoi(argv[1]);
  m=atoi(argv[2]);
  p=atoi(argv[3]);
  s=argv[4];
  c=atof(argv[5]);
  eps_I=atof(argv[6]);
  eps_B=atof(argv[7]);
  if (strcmp(s,"gaussian")==0)
    kernel = gaussian;
  else if (strcmp(s,"multiquadric")==0)
    kernel = multiquadric;
  else if (strcmp(s,"inverse_multiquadric")==0)
    kernel = inverse_multiquadric;
  else if (strcmp(s,"logarithm")==0)
    kernel = logarithm;
  else if (strcmp(s,"thinplate_spline")==0)
    kernel = thinplate_spline;
  else if (strcmp(s,"one_over_square")==0)
    kernel = one_over_square;
  else if (strcmp(s,"one_over_modulus")==0)
    kernel = one_over_modulus;
  else if (strcmp(s,"one_over_x")==0)
    kernel = one_over_x;
  else if (strcmp(s,"inverse_multiquadric3")==0)
    kernel = inverse_multiquadric3;
  else if (strcmp(s,"sinc_kernel")==0)
    kernel = sinc_kernel;
  else if (strcmp(s,"cosc")==0)
    kernel = cosc;
  else if (strcmp(s,"cot")==0)
    kernel = kcot;
  else
  {
    s="multiquadric";
    kernel = multiquadric;
  }

  bench_openmp(stdin, n, m, p, kernel, c, eps_I, eps_B);

  return 0;
}

