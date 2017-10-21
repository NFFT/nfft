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
#include "config.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3.h"
#include "infft.h"

int global_n;
int global_d;

static int comp1(const void *x, const void *y)
{
  return ((*(const R*) x) < (*(const R*) y) ? -1 : 1);
}

static int comp2(const void *x, const void *y)
{
  int nx0, nx1, ny0, ny1;
  nx0 = global_n * (int)LRINT(*((const R*) x + 0));
  nx1 = global_n * (int)LRINT(*((const R*) x + 1));
  ny0 = global_n * (int)LRINT(*((const R*) y + 0));
  ny1 = global_n * (int)LRINT(*((const R*) y + 1));

  if (nx0 < ny0)
    return -1;
  else if (nx0 == ny0)
    if (nx1 < ny1)
      return -1;
    else
      return 1;
  else
    return 1;
}

static int comp3(const void *x, const void *y)
{
  int nx0, nx1, nx2, ny0, ny1, ny2;
  nx0 = global_n * (int)LRINT(*((const R*) x + 0));
  nx1 = global_n * (int)LRINT(*((const R*) x + 1));
  nx2 = global_n * (int)LRINT(*((const R*) x + 2));
  ny0 = global_n * (int)LRINT(*((const R*) y + 0));
  ny1 = global_n * (int)LRINT(*((const R*) y + 1));
  ny2 = global_n * (int)LRINT(*((const R*) y + 2));

  if (nx0 < ny0)
    return -1;
  else if (nx0 == ny0)
    if (nx1 < ny1)
      return -1;
    else if (nx1 == ny1)
      if (nx2 < ny2)
        return -1;
      else
        return 1;
    else
      return 1;
  else
    return 1;
}

static void measure_time_nfft(int d, int N, unsigned test_ndft)
{
  int r, M, NN[d], nn[d];
  R t, t_fft, t_ndft, t_nfft;
  ticks t0, t1;

  NFFT(plan) p;
  FFTW(plan) p_fft;

  printf("\\verb+%ld+&\t", LRINT(LOG((R)(N)) / LOG((R)(2)) * (R)(d) + K(0.5)));

  for (r = 0, M = 1; r < d; r++)
  {
    M = N * M;
    NN[r] = N;
    nn[r] = 2 * N;
  }

  NFFT(init_guru)(&p, d, NN, M, nn, 2,
  PRE_PHI_HUT |
  PRE_FULL_PSI | MALLOC_F_HAT | MALLOC_X | MALLOC_F |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  p_fft = FFTW(plan_dft)(d, NN, p.f_hat, p.f, FFTW_FORWARD, FFTW_MEASURE);

  /** init pseudo random nodes */
  NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);

  global_n = nn[0];
  global_d = d;

  switch (d)
  {
    case 1:
    {
      qsort(p.x, (size_t)(p.M_total), (size_t)(d) * sizeof(R), comp1);
      break;
    }
    case 2:
    {
      qsort(p.x, (size_t)(p.M_total), (size_t)(d) * sizeof(R), comp2);
      break;
    }
    case 3:
    {
      qsort(p.x, (size_t)(p.M_total), (size_t)(d) * sizeof(R), comp3);
      break;
    }
  }

  NFFT(precompute_one_psi)(&p);

  /* init pseudo random Fourier coefficients */
  NFFT(vrand_unit_complex)(p.f_hat, p.N_total);

  /** FFT */
  t_fft = K(0.0);
  r = 0;

  while (t_fft < K(1.0))
  {
    r++;
    t0 = getticks();
    FFTW(execute)(p_fft);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_fft += t;
  }
  t_fft /= (R)(r);

  //  printf("\\verb+%.1" __FES__ "+ & \\verb+(%.1f)+ &\t",t_fft,t_fft/(p.N_total*(log(N)/log(2)*d))*auxC);
  printf("\\verb+%.1" __FES__ "+ &\t", t_fft);

  /** NDFT */
  if (test_ndft)
  {
    t_ndft = K(0.0);
    r = 0;
    while (t_ndft < K(1.0))
    {
      r++;
      t0 = getticks();
      NFFT(trafo_direct)(&p);
      t1 = getticks();
      t = NFFT(elapsed_seconds)(t1, t0);
      t_ndft += t;
    }
    t_ndft /= (R)(r);
    //printf("\\verb+%.1" __FES__ "+ & \\verb+(%d)+&\t",t_ndft,(int)round(t_ndft/(p.N_total*p.N_total)*auxC));
    printf("\\verb+%.1" __FES__ "+ &\t", t_ndft);
  }
  else
    //    printf("\\verb+*+\t&\t&\t");
    printf("\\verb+*+\t&\t");

  /** NFFT */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(1.0))
  {
    r++;
    t0 = getticks();
    switch (d)
    {
      case 1:
      {
        NFFT(trafo_1d)(&p);
        break;
      }
      case 2:
      {
        NFFT(trafo_2d)(&p);
        break;
      }
      case 3:
      {
        NFFT(trafo_3d)(&p);
        break;
      }
      default:
        NFFT(trafo)(&p);
    }
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);

  //  printf("\\verb+%.1" __FES__ "+ & \\verb+(%d)+ & \\verb+(%.1" __FES__ ")+\\\\\n",t_nfft,(int)round(t_nfft/(p.N_total*(log(N)/log(2)*d))*auxC),t_nfft/t_fft);
  printf("\\verb+%.1" __FES__ "+ & \\verb+(%3.1" __FIS__ ")+\\\\\n", t_nfft, t_nfft / t_fft);

  FFTW(destroy_plan)(p_fft);
  NFFT(finalize)(&p);
}

static void measure_time_nfft_XXX2(int d, int N, unsigned test_ndft)
{
  int r, M, NN[d], nn[d];
  R t, t_fft, t_ndft, t_nfft;
  ticks t0, t1;

  NFFT(plan) p;
  FFTW(plan) p_fft;

  printf("%ld\t", LRINT(LOG((R)(N)) / LOG((R)(2)) * (R)(d) + K(0.5)));
  fflush(stdout);

  for (r = 0, M = 1; r < d; r++)
  {
    M = N * M;
    NN[r] = N;
    nn[r] = 2 * N;
  }

  NFFT(init_guru)(&p, d, NN, M, nn, 2,
  PRE_PHI_HUT |
  PRE_PSI |
  MALLOC_F_HAT | MALLOC_X | MALLOC_F |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  p_fft = FFTW(plan_dft)(d, NN, p.f_hat, p.f, FFTW_FORWARD, FFTW_MEASURE);

  C *swapndft = (C*) NFFT(malloc)((size_t)(p.M_total) * sizeof(C));

  /** init pseudo random nodes */
  NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);

  qsort(p.x, (size_t)(p.M_total), (size_t)(d) * sizeof(R), comp1);
  //nfft_vpr_double(p.x,p.M_total,"nodes x");

  NFFT(precompute_one_psi)(&p);

  /** init pseudo random Fourier coefficients */
  NFFT(vrand_unit_complex)(p.f_hat, p.N_total);

  /** FFT */
  t_fft = K(0.0);
  r = 0;
  while (t_fft < K(0.1))
  {
    r++;
    t0 = getticks();
    FFTW(execute)(p_fft);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_fft += t;
  }
  t_fft /= (R)(r);

  printf("%.1" __FES__ "\t", t_fft);

  /** NDFT */
  if (test_ndft)
  {
    CSWAP(p.f, swapndft);
    t_ndft = K(0.0);
    r = 0;
    while (t_ndft < K(0.1))
    {
      r++;
      t0 = getticks();
      NFFT(trafo_direct)(&p);
      t1 = getticks();
      t = NFFT(elapsed_seconds)(t1, t0);
      t_ndft += t;
    }
    t_ndft /= (R)(r);
    printf("%.1" __FES__ "\t", t_ndft);
    CSWAP(p.f, swapndft);
  }
  else
    printf("\t");

  /** NFFT */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(trafo)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f, p.M_total));

  /** NFFT_1d */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(trafo_1d)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f, p.M_total));

  printf("\n");

  NFFT(free)(swapndft);
  FFTW(destroy_plan)(p_fft);
  NFFT(finalize)(&p);
}

static void measure_time_nfft_XXX3(int d, int N, unsigned test_ndft)
{
  int r, M, NN[d], nn[d];
  R t, t_fft, t_ndft, t_nfft;
  ticks t0, t1;

  NFFT(plan) p;
  FFTW(plan) p_fft;

  printf("%ld\t", LRINT(LOG((R)(N)) / LOG((R)(2)) * (R)(d) + K(0.5)));
  fflush(stdout);

  for (r = 0, M = 1; r < d; r++)
  {
    M = N * M;
    NN[r] = N;
    nn[r] = 2 * N;
  }

  NFFT(init_guru)(&p, d, NN, M, nn, 2,
  PRE_PHI_HUT |
  PRE_PSI |
  MALLOC_F_HAT | MALLOC_X | MALLOC_F |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  p_fft = FFTW(plan_dft)(d, NN, p.f, p.f_hat, FFTW_BACKWARD, FFTW_MEASURE);

  C *swapndft = (C*) NFFT(malloc)((size_t)(p.N_total) * sizeof(C));

  /** init pseudo random nodes */
  NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);

  qsort(p.x, (size_t)(p.M_total), (size_t)(d) * sizeof(R), comp1);
  //nfft_vpr_double(p.x,p.M_total,"nodes x");

  NFFT(precompute_one_psi)(&p);

  /** init pseudo random samples */
  NFFT(vrand_unit_complex)(p.f, p.N_total);

  /** FFT */
  t_fft = K(0.0);
  r = 0;
  while (t_fft < K(0.1))
  {
    r++;
    t0 = getticks();
    FFTW(execute)(p_fft);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_fft += t;
  }
  t_fft /= (R)(r);

  printf("%.1" __FES__ "\t", t_fft);

  /** NDFT */
  if (test_ndft)
  {
    CSWAP(p.f_hat, swapndft);
    t_ndft = K(0.0);
    r = 0;
    while (t_ndft < K(0.1))
    {
      r++;
      t0 = getticks();
      NFFT(adjoint_direct)(&p);
      t1 = getticks();
      t = NFFT(elapsed_seconds)(t1, t0);
      t_ndft += t;
    }
    t_ndft /= (R)(r);
    printf("%.1" __FES__ "\t", t_ndft);
    CSWAP(p.f_hat, swapndft);
  }
  else
    printf("\t");

  /** NFFT */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(adjoint)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f_hat, p.N_total));

  /** NFFT_1d */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(adjoint_1d)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f_hat, p.N_total));

  printf("\n");

  NFFT(free)(swapndft);
  FFTW(destroy_plan)(p_fft);
  NFFT(finalize)(&p);
}

static void measure_time_nfft_XXX4(int d, int N, unsigned test_ndft)
{
  int r, M, NN[d], nn[d];
  R t, t_fft, t_ndft, t_nfft;
  ticks t0, t1;

  NFFT(plan) p;
  FFTW(plan) p_fft;

  printf("%ld\t", LRINT(LOG((R)(N)) / LOG((R)(2)) * (R)(d) + K(0.5)));
  fflush(stdout);

  for (r = 0, M = 1; r < d; r++)
  {
    M = N * M;
    NN[r] = N;
    nn[r] = 2 * N;
  }

  NFFT(init_guru)(&p, d, NN, M, nn, 4,
  PRE_PHI_HUT |
  PRE_PSI |
  MALLOC_F_HAT | MALLOC_X | MALLOC_F |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  p_fft = FFTW(plan_dft)(d, NN, p.f_hat, p.f, FFTW_FORWARD, FFTW_MEASURE);

  C *swapndft = (C*) NFFT(malloc)((size_t)(p.M_total) * sizeof(C));

  /** init pseudo random nodes */
  NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);

  //for(j=0;j<2*M;j+=2)
  //   p.x[j]=0.5*p.x[j]+0.25*(p.x[j]>=0?1:-1);

  //qsort(p.x,p.M_total,d*sizeof(double),comp1);
  //nfft_vpr_double(p.x,p.M_total,"nodes x");

  NFFT(precompute_one_psi)(&p);

  /** init pseudo random Fourier coefficients */
  NFFT(vrand_unit_complex)(p.f_hat, p.N_total);

  /** FFT */
  t_fft = K(0.0);
  r = 0;
  while (t_fft < K(0.1))
  {
    r++;
    t0 = getticks();
    FFTW(execute)(p_fft);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_fft += t;
  }
  t_fft /= (R)(r);

  printf("%.1" __FES__ "\t", t_fft);

  /** init pseudo random Fourier coefficients */
  NFFT(vrand_unit_complex)(p.f_hat, p.N_total);

  /** NDFT */
  if (test_ndft)
  {
    CSWAP(p.f, swapndft);
    t_ndft = K(0.0);
    r = 0;
    while (t_ndft < K(0.1))
    {
      r++;
      t0 = getticks();
      NFFT(trafo_direct)(&p);
      t1 = getticks();
      t = NFFT(elapsed_seconds)(t1, t0);
      t_ndft += t;
    }
    t_ndft /= (R)(r);
    printf("%.1" __FES__ "\t", t_ndft);

    //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

    CSWAP(p.f, swapndft);
  }
  else
    printf("\t");

  /** NFFT */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(trafo)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f, p.M_total));

  //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

  /** NFFT_2d */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(trafo_2d)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f, p.M_total));

  //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

  printf("\n");

  NFFT(free)(swapndft);
  FFTW(destroy_plan)(p_fft);
  NFFT(finalize)(&p);
}

static void measure_time_nfft_XXX5(int d, int N, unsigned test_ndft)
{
  int r, M, NN[d], nn[d];
  R t, t_fft, t_ndft, t_nfft;
  ticks t0, t1;

  NFFT(plan) p;
  FFTW(plan) p_fft;

  printf("%ld\t", LRINT(LOG((R)(N)) / LOG((R)(2)) * (R)(d) + K(0.5)));
  fflush(stdout);

  for (r = 0, M = 1; r < d; r++)
  {
    M = N * M;
    NN[r] = N;
    nn[r] = 2 * N;
  }

  NFFT(init_guru)(&p, d, NN, M, nn, 4,
  PRE_PHI_HUT |
  PRE_PSI |
  MALLOC_F_HAT | MALLOC_X | MALLOC_F |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  p_fft = FFTW(plan_dft)(d, NN, p.f, p.f_hat, FFTW_FORWARD, FFTW_MEASURE);

  C *swapndft = (C*) NFFT(malloc)((size_t)(p.N_total) * sizeof(C));

  /** init pseudo random nodes */
  NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);

  //sort_nodes(p.x,p.d,p.M_total,

  NFFT(precompute_one_psi)(&p);

  /** init pseudo random samples */
  NFFT(vrand_unit_complex)(p.f, p.M_total);

  /** FFT */
  t_fft = K(0.0);
  r = 0;
  while (t_fft < K(0.1))
  {
    r++;
    t0 = getticks();
    FFTW(execute)(p_fft);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_fft += t;
  }
  t_fft /= (R)(r);

  printf("%.1" __FES__ "\t", t_fft);

  /** init pseudo random samples */
  NFFT(vrand_unit_complex)(p.f, p.M_total);

  /** NDFT */
  if (test_ndft)
  {
    CSWAP(p.f_hat, swapndft);
    t_ndft = K(0.0);
    r = 0;
    while (t_ndft < K(0.1))
    {
      r++;
      t0 = getticks();
      NFFT(adjoint_direct)(&p);
      t1 = getticks();
      t = NFFT(elapsed_seconds)(t1, t0);
      t_ndft += t;
    }
    t_ndft /= (R)(r);
    printf("%.1" __FES__ "\t", t_ndft);

    //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

    CSWAP(p.f_hat, swapndft);
  }
  else
    printf("\t");

  /** NFFT */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(adjoint)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f_hat, p.N_total));

  //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

  /** NFFT_2d */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(adjoint_2d)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f_hat, p.N_total));

  //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

  printf("\n");

  NFFT(free)(swapndft);
  FFTW(destroy_plan)(p_fft);
  NFFT(finalize)(&p);
}

static void measure_time_nfft_XXX6(int d, int N, unsigned test_ndft)
{
  int r, M, NN[d], nn[d];
  R t, t_fft, t_ndft, t_nfft;
  ticks t0, t1;

  NFFT(plan) p;
  FFTW(plan) p_fft;

  printf("%ld\t", LRINT(LOG((R)(N)) / LOG((R)(2)) * (R)(d) + K(0.5)));
  fflush(stdout);

  for (r = 0, M = 1; r < d; r++)
  {
    M = N * M;
    NN[r] = N;
    nn[r] = 2 * N;
  }

  NFFT(init_guru)(&p, d, NN, M, nn, 2,
  PRE_PHI_HUT |
  FG_PSI |
  MALLOC_F_HAT | MALLOC_X | MALLOC_F |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  p_fft = FFTW(plan_dft)(d, NN, p.f_hat, p.f, FFTW_FORWARD, FFTW_MEASURE);

  C *swapndft = (C*) NFFT(malloc)((size_t)(p.M_total) * sizeof(C));

  /** init pseudo random nodes */
  NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);

  //qsort(p.x,p.M_total,d*sizeof(double),comp1);
  //nfft_vpr_double(p.x,p.M_total,"nodes x");

  NFFT(precompute_one_psi)(&p);

  /** init pseudo random Fourier coefficients */
  NFFT(vrand_unit_complex)(p.f_hat, p.N_total);

  /** FFT */
  t_fft = K(0.0);
  r = 0;
  while (t_fft < K(0.1))
  {
    r++;
    t0 = getticks();
    FFTW(execute)(p_fft);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_fft += t;
  }
  t_fft /= (R)(r);

  printf("%.1" __FES__ "\t", t_fft);

  /** init pseudo random Fourier coefficients */
  NFFT(vrand_unit_complex)(p.f_hat, p.N_total);

  /** NDFT */
  if (test_ndft)
  {
    CSWAP(p.f, swapndft);
    t_ndft = K(0.0);
    r = 0;
    while (t_ndft < K(0.1))
    {
      r++;
      t0 = getticks();
      NFFT(trafo_direct)(&p);
      t1 = getticks();
      t = NFFT(elapsed_seconds)(t1, t0);
      t_ndft += t;
    }
    t_ndft /= (R)(r);
    printf("%.1" __FES__ "\t", t_ndft);

    //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

    CSWAP(p.f, swapndft);
  }
  else
    printf("\t");

  /** NFFT */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(trafo)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f, p.M_total));

  //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

  /** NFFT_3d */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(trafo_3d)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f, p.M_total));

  //printf("f=%e+i%e\t",creal(p.f[0]),cimag(p.f[0]));

  printf("\n");

  NFFT(free)(swapndft);
  FFTW(destroy_plan)(p_fft);
  NFFT(finalize)(&p);
}

static void measure_time_nfft_XXX7(int d, int N, unsigned test_ndft)
{
  int r, M, NN[d], nn[d];
  R t, t_fft, t_ndft, t_nfft;
  ticks t0, t1;

  NFFT(plan) p;
  FFTW(plan) p_fft;

  printf("%ld\t", LRINT(LOG((R)(N)) / LOG((R)(2)) * (R)(d) + K(0.5)));
  fflush(stdout);

  for (r = 0, M = 1; r < d; r++)
  {
    M = N * M;
    NN[r] = N;
    nn[r] = 2 * N;
  }

  NFFT(init_guru)(&p, d, NN, M, nn, 2,
  PRE_PHI_HUT |
  FG_PSI |
  MALLOC_F_HAT | MALLOC_X | MALLOC_F |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  p_fft = FFTW(plan_dft)(d, NN, p.f, p.f_hat, FFTW_FORWARD, FFTW_MEASURE);

  C *swapndft = (C*) NFFT(malloc)((size_t)(p.N_total) * sizeof(C));

  /** init pseudo random nodes */
  NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);

  //sort_nodes(p.x,p.d,p.M_total,

  NFFT(precompute_one_psi)(&p);

  /** init pseudo random samples */
  NFFT(vrand_unit_complex)(p.f, p.M_total);

  /** FFT */
  t_fft = K(0.0);
  r = 0;
  while (t_fft < K(0.1))
  {
    r++;
    t0 = getticks();
    FFTW(execute)(p_fft);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_fft += t;
  }
  t_fft /= (R)(r);

  printf("%.1" __FES__ "\t", t_fft);

  /** init pseudo random samples */
  NFFT(vrand_unit_complex)(p.f, p.M_total);

  /** NDFT */
  if (test_ndft)
  {
    CSWAP(p.f_hat, swapndft);
    t_ndft = K(0.0);
    r = 0;
    while (t_ndft < K(0.1))
    {
      r++;
      t0 = getticks();
      NFFT(adjoint_direct)(&p);
      t1 = getticks();
      t = NFFT(elapsed_seconds)(t1, t0);
      t_ndft += t;
    }
    t_ndft /= (R)(r);
    printf("%.1" __FES__ "\t", t_ndft);

    //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

    CSWAP(p.f_hat, swapndft);
  }
  else
    printf("\t");

  /** NFFT */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(adjoint)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f_hat, p.N_total));

  //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

  /** NFFT_3d */
  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(adjoint_3d)(&p);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);
  printf("%.1" __FES__ "\t", t_nfft);
  if (test_ndft)
    printf("(%.1" __FES__ ")\t", NFFT(error_l_2_complex)(swapndft, p.f_hat, p.N_total));

  //printf("\nf_hat=%e+i%e\t",creal(p.f_hat[0]),cimag(p.f_hat[0]));

  printf("\n");

  NFFT(free)(swapndft);
  FFTW(destroy_plan)(p_fft);
  NFFT(finalize)(&p);
}

//static int main(void)
//{
//  int l, d, logIN;
//
//  for (l = 3; l <= 6; l++)
//  {
//    d = 3;
//    logIN = d * l;
//    int N = (int)(1U << (logIN / d));
//    measure_time_nfft_XXX6(d, N, logIN <= 15 ? 1 : 0);
//    measure_time_nfft_XXX7(d, N, logIN <= 15 ? 1 : 0);
//  }
//
//  for (l = 7; l <= 12; l++)
//  {
//    d = 2;
//    logIN = d * l;
//    int N = (int)(1U << (logIN / d));
//    measure_time_nfft_XXX4(d, N, logIN <= 15 ? 1 : 0);
//    measure_time_nfft_XXX5(d, N, logIN <= 15 ? 1 : 0);
//  }
//
//  for (l = 3; l <= 12; l++)
//  {
//    logIN = l;
//    int N = (int)(1U << (logIN));
//    measure_time_nfft_XXX2(1, N, logIN <= 15 ? 1 : 0);
//    measure_time_nfft_XXX3(1, N, logIN <= 15 ? 1 : 0);
//  }
//
//  return EXIT_SUCCESS;
//}

int main(void)
{
  int l, d, logIN;

  UNUSED(measure_time_nfft_XXX2);
  UNUSED(measure_time_nfft_XXX3);
  UNUSED(measure_time_nfft_XXX4);
  UNUSED(measure_time_nfft_XXX5);
  UNUSED(measure_time_nfft_XXX6);
  UNUSED(measure_time_nfft_XXX7);

  printf("\\hline $l_N$ & FFT & NDFT & NFFT & NFFT/FFT\\\\\n");
  printf("\\hline \\hline \\multicolumn{5}{|c|}{$d=1$} \\\\ \\hline\n");
  for (l = 3; l <= 22; l++)
  {
    d = 1;
    logIN = l;
    int N = (int)(1U << (logIN / d));
    measure_time_nfft(d, N, logIN <= 15 ? 1 : 0);

    fflush(stdout);
  }

  printf("\\hline $l_N$ & FFT & NDFT & NFFT & NFFT/FFT\\\\\n");
  printf("\\hline \\hline \\multicolumn{5}{|c|}{$d=2$} \\\\ \\hline\n");
  for (l = 3; l <= 11; l++)
  {
    d = 2;
    logIN = d * l;
    int N = (int)(1U << (logIN / d));
    measure_time_nfft(d, N, logIN <= 15 ? 1 : 0);

    fflush(stdout);
  }

  printf("\\hline \\hline \\multicolumn{5}{|c|}{$d=3$} \\\\ \\hline\n");
  for (l = 3; l <= 7; l++)
  {
    d = 3;
    logIN = d * l;
    int N = (int)(1U << (logIN / d));
    measure_time_nfft(d, N, logIN <= 15 ? 1 : 0);

    fflush(stdout);
  }

  return 1;
}
