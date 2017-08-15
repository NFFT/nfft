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

/*! \file ndft_fast.c
 *
 * \brief Testing ndft, Horner-like ndft, and fully precomputed ndft.
 *
 * \author Stefan Kunis
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

static void ndft_horner_trafo(NFFT(plan) *ths)
{
  INT j, k;
  C *f_hat_k, *f_j;
  C exp_omega_0;

  for (j = 0, f_j = ths->f; j < ths->M_total; j++, f_j++)
    (*f_j) = K(0.0);

  for (j = 0, f_j = ths->f; j < ths->M_total; j++, f_j++)
  {
    exp_omega_0 = CEXP(+K2PI * II * ths->x[j]);
    for (k = 0, f_hat_k = ths->f_hat; k < ths->N[0]; k++, f_hat_k++)
    {
      (*f_j) += (*f_hat_k);
      (*f_j) *= exp_omega_0;
    }
    (*f_j) *= CEXP(-KPI * II * (R)(ths->N[0]) * ths->x[j]);
  }
} /* ndft_horner_trafo */

static void ndft_pre_full_trafo(NFFT(plan) *ths, C *A)
{
  INT j, k;
  C *f_hat_k, *f_j;
  C *A_local;

  for (j = 0, f_j = ths->f; j < ths->M_total; j++, f_j++)
    (*f_j) = K(0.0);

  for (j = 0, f_j = ths->f, A_local = A; j < ths->M_total; j++, f_j++)
    for (k = 0, f_hat_k = ths->f_hat; k < ths->N[0]; k++, f_hat_k++, A_local++)
      (*f_j) += (*f_hat_k) * (*A_local);
} /* ndft_pre_full_trafo */

static void ndft_pre_full_init(NFFT(plan) *ths, C *A)
{
  INT j, k;
  C *f_hat_k, *f_j, *A_local;

  for (j = 0, f_j = ths->f, A_local = A; j < ths->M_total; j++, f_j++)
    for (k = 0, f_hat_k = ths->f_hat; k < ths->N[0]; k++, f_hat_k++, A_local++)
      (*A_local) = CEXP(
          -K2PI * II * ((R) (k) - (R) (ths->N[0]) / K(2.0)) * ths->x[j]);

} /* ndft_pre_full_init */

static void ndft_time(int N, int M, unsigned test_ndft, unsigned test_pre_full)
{
  int r;
  R t, t_ndft, t_horner, t_pre_full, t_nfft;
  C *A = NULL;
  ticks t0, t1;

  NFFT(plan) np;

  printf("%d\t%d\t", N, M);

  NFFT(init_1d)(&np, N, M);

  /* Initialize pseudo random nodes. */
  NFFT(vrand_shifted_unit_double)(np.x, np.M_total);

  if (test_pre_full)
  {
    A = (C*) NFFT(malloc)((size_t)(N * M) * sizeof(R));
    ndft_pre_full_init(&np, A);
  }

  /* Initialize pseudo random Fourier coefficients. */
  NFFT(vrand_unit_complex)(np.f_hat, np.N_total);

  /* NDFT */
  if (test_ndft)
  {
    t_ndft = K(0.0);
    r = 0;
    while (t_ndft < K(0.1))
    {
      r++;
      t0 = getticks();
      NFFT(trafo_direct)(&np);
      t1 = getticks();
      t = NFFT(elapsed_seconds)(t1, t0);
      t_ndft += t;
    }
    t_ndft /= (R) (r);

    printf("%.2" __FES__ "\t", t_ndft);
  }
  else
    printf("N/A\t\t");

  /* Horner NDFT */
  t_horner = K(0.0);
  r = 0;
  while (t_horner < K(0.1))
  {
    r++;
    t0 = getticks();
    ndft_horner_trafo(&np);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_horner += t;
  }
  t_horner /= (R)(r);

  printf("%.2" __FES__ "\t", t_horner);

  /* Fully precomputed NDFT */
  if (test_pre_full)
  {
    t_pre_full = K(0.0);
    r = 0;
    while (t_pre_full < K(0.1))
    {
      r++;
      t0 = getticks();
      ndft_pre_full_trafo(&np, A);
      t1 = getticks();
      t = NFFT(elapsed_seconds)(t1, t0);
      t_pre_full += t;
    }
    t_pre_full /= (R)(r);

    printf("%.2" __FES__ "\t", t_pre_full);
  }
  else
    printf("N/A\t\t");

  t_nfft = K(0.0);
  r = 0;
  while (t_nfft < K(0.1))
  {
    r++;
    t0 = getticks();
    NFFT(trafo)(&np);
    t1 = getticks();
    t = NFFT(elapsed_seconds)(t1, t0);
    t_nfft += t;
  }
  t_nfft /= (R)(r);

  printf("%.2" __FES__ "\n", t_nfft);

  fflush(stdout);

  if (test_pre_full)
    NFFT(free)(A);

  NFFT(finalize)(&np);
}

int main(int argc, char **argv)
{
  int l, trial;

  if (argc < 4)
  {
    fprintf(stderr, "ndft_fast type first last trials\n");
    return EXIT_FAILURE;
  }
  else
  {
    int arg2 = (atoi(argv[2]));
    int arg3 = (atoi(argv[3]));
    int arg4 = (atoi(argv[4]));
    fprintf(stderr, "Testing ndft, Horner-like ndft, fully precomputed ndft.\n");
    fprintf(stderr, "Columns: N, M, t_ndft, t_horner, t_pre_full, t_nfft\n\n");

    /* time vs. N=M */
    if (atoi(argv[1]) == 0)
    {
      for (l = arg2; l <= arg3; l++)
      {
        for (trial = 0; trial < arg4; trial++)
        {
          int N = (int)(1U << l);
          int M = (int)(1U << l);
          ndft_time(N, M, l <= 15 ? 1 : 0, l <= 13 ? 1 : 0);
        }
      }
    }
  }

  return EXIT_SUCCESS;
}
