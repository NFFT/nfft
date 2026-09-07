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

/*! \file flags.c
 *
 * \brief Testing the nfft.
 *
 * \author Stefan Kunis
 *
 * References: Time and Memory Requirements of the Nonequispaced FFT
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
#include "nfft3mp.h"
#include "nfft3util.h"

#ifdef GAUSSIAN
  unsigned test_fg=1;
#else
  unsigned test_fg=0;
#endif

#ifdef MEASURE_TIME_FFTW
  unsigned test_fftw=1;
#else
  unsigned test_fftw=0;
#endif

#ifdef MEASURE_TIME
  unsigned test=1;
#else
  unsigned test=0;
#endif

static void flags_cp(NFFT(plan) *dst, NFFT(plan) *src)
{
  dst->x = src->x;
  dst->f_hat = src->f_hat;
  dst->f = src->f;
  dst->g1 = src->g1;
  dst->g2 = src->g2;
  dst->my_fftw_plan1 = src->my_fftw_plan1;
  dst->my_fftw_plan2 = src->my_fftw_plan2;
}

static void time_accuracy(int d, int N, int M, int n, int m, unsigned test_ndft,
    unsigned test_pre_full_psi)
{
  int r, NN[d], nn[d];
  NFFT_R t_ndft, t, e;
  NFFT_C *swapndft = NULL;
  double t0, t1;

  NFFT(plan) p;
  NFFT(plan) p_pre_phi_hut;
  NFFT(plan) p_fg_psi;
  NFFT(plan) p_pre_lin_psi;
  NFFT(plan) p_pre_fg_psi;
  NFFT(plan) p_pre_psi;
  NFFT(plan) p_pre_full_psi;

  printf("%d\t%d\t", d, N);

  for (r = 0; r < d; r++)
  {
    NN[r] = N;
    nn[r] = n;
  }

  /* output vector ndft */
  if (test_ndft)
    swapndft = (NFFT_C*) NFFT(malloc)((size_t)(M) * sizeof(NFFT_C));

  NFFT(init_guru)(&p, d, NN, M, nn, m,
  MALLOC_X | MALLOC_F_HAT | MALLOC_F |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  /** init pseudo random nodes */
  NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);

  NFFT(init_guru)(&p_pre_phi_hut, d, NN, M, nn, m, PRE_PHI_HUT, 0);
  flags_cp(&p_pre_phi_hut, &p);
  NFFT(precompute_one_psi)(&p_pre_phi_hut);

  if (test_fg)
  {
    NFFT(init_guru)(&p_fg_psi, d, NN, M, nn, m, FG_PSI, 0);
    flags_cp(&p_fg_psi, &p);
    NFFT(precompute_one_psi)(&p_fg_psi);
  }

  NFFT(init_guru)(&p_pre_lin_psi, d, NN, M, nn, m, PRE_LIN_PSI, 0);
  flags_cp(&p_pre_lin_psi, &p);
  NFFT(precompute_one_psi)(&p_pre_lin_psi);

  if (test_fg)
  {
    NFFT(init_guru)(&p_pre_fg_psi, d, NN, M, nn, m, PRE_FG_PSI, 0);
    flags_cp(&p_pre_fg_psi, &p);
    NFFT(precompute_one_psi)(&p_pre_fg_psi);
  }

  NFFT(init_guru)(&p_pre_psi, d, NN, M, nn, m, PRE_PSI, 0);
  flags_cp(&p_pre_psi, &p);
  NFFT(precompute_one_psi)(&p_pre_psi);

  if (test_pre_full_psi)
  {
    NFFT(init_guru)(&p_pre_full_psi, d, NN, M, nn, m, PRE_FULL_PSI, 0);
    flags_cp(&p_pre_full_psi, &p);
    NFFT(precompute_one_psi)(&p_pre_full_psi);
  }

  /* init pseudo random Fourier coefficients */
  NFFT(vrand_unit_complex)(p.f_hat, p.N_total);

  /* NDFT */
  if (test_ndft)
  {
    NFFT_CSWAP(p.f, swapndft);

    t_ndft = NFFT_K(0.0);
    r = 0;
    while (t_ndft < NFFT_K(0.01))
    {
      r++;
      t0 = NFFT(clock_gettime_seconds)();
      NFFT(trafo_direct)(&p);
      t1 = NFFT(clock_gettime_seconds)();
      t = t1 - t0;
      t_ndft += t;
    }
    t_ndft /= (NFFT_R)(r);

    NFFT_CSWAP(p.f, swapndft);
  }
  else
    t_ndft = NFFT_MKNAN("");

  /* NFFTs */
  NFFT(trafo)(&p);
  NFFT(trafo)(&p_pre_phi_hut);
  if (test_fg)
    NFFT(trafo)(&p_fg_psi);
  else
    p_fg_psi.MEASURE_TIME_t[2] = NFFT_MKNAN("");
  NFFT(trafo)(&p_pre_lin_psi);
  if (test_fg)
    NFFT(trafo)(&p_pre_fg_psi);
  else
    p_pre_fg_psi.MEASURE_TIME_t[2] = NFFT_MKNAN("");
  NFFT(trafo)(&p_pre_psi);
  if (test_pre_full_psi)
    NFFT(trafo)(&p_pre_full_psi);
  else
    p_pre_full_psi.MEASURE_TIME_t[2] = NFFT_MKNAN("");

  if (test_fftw == 0)
    p.MEASURE_TIME_t[1] = NFFT_MKNAN("");

  if (test_ndft)
    e = NFFT(error_l_2_complex)(swapndft, p.f, p.M_total);
  else
    e = NFFT_MKNAN("");

  printf(
      "%.2" NFFT__FES__ "\t%d\t%.2" NFFT__FES__ "\t%.2" NFFT__FES__ "\t%.2" NFFT__FES__ "\t%.2" NFFT__FES__ "\t%.2" NFFT__FES__ "\t%.2" NFFT__FES__ "\t%.2" NFFT__FES__ "\t%.2" NFFT__FES__ "\t%.2" NFFT__FES__ "\t%.2" NFFT__FES__ "\n",
      t_ndft, m, e, p.MEASURE_TIME_t[0], p_pre_phi_hut.MEASURE_TIME_t[0],
      p.MEASURE_TIME_t[1], p.MEASURE_TIME_t[2], p_fg_psi.MEASURE_TIME_t[2],
      p_pre_lin_psi.MEASURE_TIME_t[2], p_pre_fg_psi.MEASURE_TIME_t[2],
      p_pre_psi.MEASURE_TIME_t[2], p_pre_full_psi.MEASURE_TIME_t[2]);

  fflush(stdout);

  /** finalise */
  if (test_pre_full_psi)
    NFFT(finalize)(&p_pre_full_psi);
  NFFT(finalize)(&p_pre_psi);
  if (test_fg)
    NFFT(finalize)(&p_pre_fg_psi);
  NFFT(finalize)(&p_pre_lin_psi);
  if (test_fg)
    NFFT(finalize)(&p_fg_psi);
  NFFT(finalize)(&p_pre_phi_hut);
  NFFT(finalize)(&p);

  if (test_ndft)
    NFFT(free)(swapndft);
}

static void accuracy_pre_lin_psi(int d, int N, int M, int n, int m, int K)
{
  int r, NN[d], nn[d];
  NFFT_R e;
  NFFT_C *swapndft;

  NFFT(plan) p;

  for (r = 0; r < d; r++)
  {
    NN[r] = N;
    nn[r] = n;
  }

  /* output vector ndft */
  swapndft = (NFFT_C*) NFFT(malloc)((size_t)(M) * sizeof(NFFT_C));

  NFFT(init_guru)(&p, d, NN, M, nn, m,
  MALLOC_X | MALLOC_F_HAT | MALLOC_F |
  PRE_PHI_HUT | PRE_LIN_PSI |
  FFTW_INIT | FFT_OUT_OF_PLACE,
  FFTW_MEASURE | FFTW_DESTROY_INPUT);

  /** realloc psi */
  NFFT(free)(p.psi);
  p.K = K;
  p.psi = (NFFT_R*) NFFT(malloc)((size_t)((p.K + 1) * p.d) * sizeof(NFFT_R));

  /** precomputation can be done before the nodes are initialised */
  NFFT(precompute_one_psi)(&p);

  /** init pseudo random nodes */
  NFFT(vrand_shifted_unit_double)(p.x, p.d * p.M_total);

  /** init pseudo random Fourier coefficients */
  NFFT(vrand_unit_complex)(p.f_hat, p.N_total);

  /** compute exact result */
  NFFT_CSWAP(p.f, swapndft);
  NFFT(trafo_direct)(&p);
  NFFT_CSWAP(p.f, swapndft);

  /** NFFT */
  NFFT(trafo)(&p);
  e = NFFT(error_l_2_complex)(swapndft, p.f, p.M_total);

  //  printf("%d\t%d\t%d\t%d\t%.2e\n",d,N,m,K,e);
  printf("$%.1" NFFT__FES__ "$&\t", e);

  fflush(stdout);

  /** finalise */
  NFFT(finalize)(&p);
  NFFT(free)(swapndft);
}

int main(int argc, char **argv)
{
  int l, trial;

  if (argc <= 2)
  {
    fprintf(stderr, "flags type first last trials d m\n");
    return EXIT_FAILURE;
  }

  if ((test == 0) && (atoi(argv[1]) < 2))
  {
    fprintf(stderr, "configure with --enable-measure-time\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "Testing different precomputation schemes for the nfft.\n");
  fprintf(stderr, "Columns: d, N=M, t_ndft, e_nfft, t_D, t_pre_phi_hut, ");
  fprintf(stderr, "t_fftw, t_B, t_fg_psi, t_pre_lin_psi, t_pre_fg_psi, ");
  fprintf(stderr, "t_pre_psi, t_pre_full_psi\n\n");

  int arg2 = atoi(argv[2]);
  int arg3 = atoi(argv[3]);
  int arg4 = atoi(argv[4]);

  /* time vs. N=M */
  if (atoi(argv[1]) == 0)
  {
    int d = atoi(argv[5]);
    int m = atoi(argv[6]);

    for (l = arg2; l <= arg3; l++)
    {
      int N = (int)(1U << l);
      int M = (int)(1U << (d * l));
      for (trial = 0; trial < arg4; trial++)
      {
        time_accuracy(d, N, M, 2 * N, m, 0, 0);
      }
    }
  }
  else if (atoi(argv[1]) == 1) /* accuracy vs. time */
  {
    int d = atoi(argv[5]);
    int N = atoi(argv[6]);
    int m;

    for (m = arg2; m <= arg3; m++)
    {
      for (trial = 0; trial < arg4; trial++)
      {
        time_accuracy(d, N, (int)(NFFT_LRINT(NFFT_POW((NFFT_R)(N), (NFFT_R)(d)))), 2 * N, m, 1, 1);
      }
    }
  }
  else if (atoi(argv[1]) == 2) /* accuracy vs. K for linear interpolation, assumes (m+1)|K */
  {
    int d = atoi(argv[5]);
    int N = atoi(argv[6]);
    int m = atoi(argv[7]);

    printf("$\\log_2(K/(m+1))$&\t");

    for (l = arg2; l < arg3; l++)
      printf("$%d$&\t", l);

    printf("$%d$\\\\\n", arg3);

    printf("$\\tilde E_2$&\t");
    for (l = arg2; l <= arg3; l++)
    {
      int x = (m + 1) * (int)(1U << l);
      accuracy_pre_lin_psi(d, N, (int)(NFFT_LRINT(NFFT_POW((NFFT_R)(N), (NFFT_R)(d)))), 2 * N, m, x);
    }

    printf("\n");
  }


  return EXIT_SUCCESS;
}
