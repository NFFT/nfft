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

void bench_openmp(FILE *infile, int m, int psi_flag)
{
  NFFT(plan) p;
  int *N;
  int *n;
  int M, d, trafo_adjoint;
  int t, j;
  double re,im;
  ticks t0, t1;
  double tt_total, tt_preonepsi;

  fscanf(infile, "%d %d", &d, &trafo_adjoint);

  N = malloc(d*sizeof(int));
  n = malloc(d*sizeof(int));

  for (t=0; t<d; t++)
    fscanf(infile, "%d", N+t);

  for (t=0; t<d; t++)
    fscanf(infile, "%d", n+t);

  fscanf(infile, "%d", &M);

#ifdef _OPENMP
  FFTW(import_wisdom_from_filename)("nfft_benchomp_detail_threads.plan");
#else
  FFTW(import_wisdom_from_filename)("nfft_benchomp_detail_single.plan");
#endif

  /** init an d-dimensional plan */
  NFFT(init_guru)(&p, d, N, M, n, m,
                   PRE_PHI_HUT| psi_flag | MALLOC_X | MALLOC_F_HAT| MALLOC_F| FFTW_INIT | FFT_OUT_OF_PLACE,
                   FFTW_MEASURE| FFTW_DESTROY_INPUT);

#ifdef _OPENMP
  FFTW(export_wisdom_to_filename)("nfft_benchomp_detail_threads.plan");
#else
  FFTW(export_wisdom_to_filename)("nfft_benchomp_detail_single.plan");
#endif

  for (j=0; j < p.M_total; j++)
  {
    for (t=0; t < p.d; t++)
      fscanf(infile, "%lg", p.x+p.d*j+t);
  }

  if (trafo_adjoint==0)
  {
    for (j=0; j < p.N_total; j++)
    {
      fscanf(infile, "%lg %lg", &re, &im);
      p.f_hat[j] = re + _Complex_I * im;
    }
  }
  else
  {
    for (j=0; j < p.M_total; j++)
    {
      fscanf(infile, "%lg %lg", &re, &im);
      p.f[j] = re + _Complex_I * im;
    }
  }

  t0 = getticks();
  /** precompute psi, the entries of the matrix B */
  if(p.flags & PRE_ONE_PSI)
      NFFT(precompute_one_psi)(&p);
  t1 = getticks();
  tt_preonepsi = NFFT(elapsed_seconds)(t1,t0);

  if (trafo_adjoint==0)
    NFFT(trafo)(&p);
  else
    NFFT(adjoint)(&p);
  t1 = getticks();
  tt_total = NFFT(elapsed_seconds)(t1,t0);

#ifndef MEASURE_TIME
  p.MEASURE_TIME_t[0] = 0.0;
  p.MEASURE_TIME_t[2] = 0.0;
#endif

#ifndef MEASURE_TIME_FFTW
  p.MEASURE_TIME_t[1] = 0.0;
#endif

  printf("%.6e %.6e %6e %.6e %.6e %.6e\n", tt_preonepsi, p.MEASURE_TIME_t[0], p.MEASURE_TIME_t[1], p.MEASURE_TIME_t[2], tt_total-tt_preonepsi-p.MEASURE_TIME_t[0]-p.MEASURE_TIME_t[1]-p.MEASURE_TIME_t[2], tt_total);
//  printf("%.6e\n", tt);

  free(N);
  free(n);

  /** finalise the one dimensional plan */
  NFFT(finalize)(&p);
}

int main(int argc, char **argv)
{
  int m, psi_flag;
#ifdef _OPENMP
  int nthreads;

  if (argc != 4)
    return 1;

  nthreads = atoi(argv[3]);
  FFTW(init_threads)();
  omp_set_num_threads(nthreads);
#else
  if (argc != 3)
    return 1;
#endif

  m = atoi(argv[1]);
  psi_flag = atoi(argv[2]);

  bench_openmp(stdin, m, psi_flag);

  return 0;
}
