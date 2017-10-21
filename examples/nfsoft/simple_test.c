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

/* Include standard C headers. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

/* Include NFFT3 library header. */
#include "nfft3.h"

static void simple_test_nfsoft(int bw, int M)
{
  nfsoft_plan plan_nfsoft; /**< Plan for the NFSOFT   */
  nfsoft_plan plan_ndsoft; /**< Plan for the NDSOFT   */

  double t0, t1;
  int j; /** just an index*/
  int k, m; /** the two parameters controlling the accuracy of the NFSOFT*/
  double d1, d2, d3; /** indeces for initializing the Euler angles*/
  double time, error; /**...self-explainatory*/
  unsigned int flags = NFSOFT_MALLOC_X | NFSOFT_MALLOC_F | NFSOFT_MALLOC_F_HAT; /**flags for memory allocation \see nfft3.h*/

  /**set the accuracy controlling parameters*/
  k = 1000; /**  k resembles the FPT kappa */
  m = 5; /**  m the NFFT cut-off parameter*/

  /** Init two transform plans using the guru interface.
   * All arrays for input and
   * output variables are allocated by nfsoft_init_guru(). The
   * array of SO(3) Fourier coefficients is preserved during
   * transformations.
   */

  nfsoft_init_guru(&plan_ndsoft, bw, M, flags | NFSOFT_USE_NDFT
      | NFSOFT_USE_DPT, PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT
      | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE, m, k);

  nfsoft_init_guru(&plan_nfsoft, bw, M, flags, PRE_PHI_HUT | PRE_PSI | MALLOC_X
      | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE, m, k);

  /** Init random nodes (for both plans, the same). */
  for (j = 0; j < plan_nfsoft.M_total; j++)
  {
    d1 = ((double) rand()) / RAND_MAX - 0.5;
    d2 = 0.5 * ((double) rand()) / RAND_MAX;
    d3 = ((double) rand()) / RAND_MAX - 0.5;

    plan_nfsoft.x[3* j ] = d1; /**alpha*/
    plan_nfsoft.x[3* j + 1] = d2; /**beta*/
    plan_nfsoft.x[3* j + 2] = d3; /**gamma*/

    plan_ndsoft.x[3* j ] = d1; /**alpha*/
    plan_ndsoft.x[3* j + 1] = d2; /**beta*/
    plan_ndsoft.x[3* j + 2] = d3; /**gamma*/
  }

  /** init random Fourier coefficients (again the same for both plans) and display them*/
  for (j = 0; j < (bw + 1) * (4* (bw +1)*(bw+1)-1)/3;j++)
  {
    d1=((double)rand())/RAND_MAX - 0.5;
    d2=((double)rand())/RAND_MAX - 0.5;
    plan_nfsoft.f_hat[j]=d1 + I*d2;
    plan_ndsoft.f_hat[j]=d1 + I*d2;
  }

  if ((bw+1)*(4*(bw+1)*(bw+1)-1)/3<=20)
  nfft_vpr_complex(plan_nfsoft.f_hat,(bw+1)*(4*(bw+1)*(bw+1)-1)/3,"randomly generated SO(3) Fourier coefficients");
	else
  nfft_vpr_complex(plan_ndsoft.f_hat,20,"1st-20th randomly generated SO(3) Fourier coefficient");

  printf("\n---------------------------------------------\n");

  /**Do precomputation for all transforms*/
  nfsoft_precompute(&plan_nfsoft);
  nfsoft_precompute(&plan_ndsoft);


  /** Compute NFSOFT and display the time needed. */
  t0 = nfft_clock_gettime_seconds();
  nfsoft_trafo(&plan_nfsoft);
  t1 = nfft_clock_gettime_seconds();
  time = t1-t0;
  if (plan_nfsoft.M_total<=20)
    nfft_vpr_complex(plan_nfsoft.f,plan_nfsoft.M_total,"NFSOFT, function samples");
  else
    nfft_vpr_complex(plan_nfsoft.f,20,"NFSOFT, 1st-20th function sample");
  printf(" computed in %11le seconds\n",time);

  /** Compute NDSOFT and display the time needed. */
  t0 = nfft_clock_gettime_seconds();
  nfsoft_trafo(&plan_ndsoft);
  t1 = nfft_clock_gettime_seconds();
  time = t1-t0;
  if (plan_ndsoft.M_total<=20)
    nfft_vpr_complex(plan_ndsoft.f,plan_ndsoft.M_total,"NDSOFT, function samples");
  else
    nfft_vpr_complex(plan_ndsoft.f,20,"NDSOFT, 1st-20th function sample");
  printf(" computed in %11le seconds\n",time);

  /**compute the error between the NFSOFT and NDSOFT and display it*/
  error= nfft_error_l_infty_complex(plan_ndsoft.f,plan_nfsoft.f, plan_nfsoft.M_total);
  printf("\n The NFSOFT of bandwidth=%d for %d rotations has infty-error %11le \n",bw, M,error);

  printf("\n---------------------------------------------\n");

  plan_nfsoft.f[0]=1.0;
  plan_ndsoft.f[0]=1.0;
  nfft_vpr_complex(plan_ndsoft.f,plan_ndsoft.M_total, "function samples to start adjoint trafo");

  /** Compute the adjoint NFSOFT and display the time needed.*/
  t0 = nfft_clock_gettime_seconds();
  nfsoft_adjoint(&plan_nfsoft);
  t1 = nfft_clock_gettime_seconds();
  time = t1-t0;
  if ((bw+1)*(4*(bw+1)*(bw+1)-1)/3<=20)
     nfft_vpr_complex(plan_nfsoft.f_hat,(bw+1)*(4*(bw+1)*(bw+1)-1)/3,"SO(3) Fourier coefficients");
  else
     nfft_vpr_complex(plan_nfsoft.f_hat,20,"adjoint NFSOFT, 1st-20th Fourier coefficient");
  printf(" computed in %11le seconds\n",time);

  /** Compute adjoint NDSOFT and display the time needed.*/
  t0 = nfft_clock_gettime_seconds();
  nfsoft_adjoint(&plan_ndsoft);
  t1 = nfft_clock_gettime_seconds();
  time = t1-t0;
  if ((bw+1)*(4*(bw+1)*(bw+1)-1)/3<=20)
	nfft_vpr_complex(plan_ndsoft.f_hat,(bw+1)*(4*(bw+1)*(bw+1)-1)/3,"SO(3) Fourier coefficients");
  else
    nfft_vpr_complex(plan_ndsoft.f_hat,20,"adjoint NDSOFT, 1st-20th Fourier coefficient");
  printf(" computed in %11le seconds\n",time);


  /**compute the error between the adjoint NFSOFT and NDSOFT and display it*/
  error=nfft_error_l_infty_complex(plan_ndsoft.f_hat,plan_nfsoft.f_hat, (bw+1)*(4*(bw+1)*(bw+1)-1)/3);
  printf("\n The adjoint NFSOFT of bandwidth=%d for %d rotations has infty-error %11le \n",bw, M,error);

  printf("\n---------------------------------------------\n");

  /**destroy the plans*/
  nfsoft_finalize(&plan_ndsoft);
  nfsoft_finalize(&plan_nfsoft);
}

  /**
   * The main program.
   *
   * computes an NDSOFT and its adjoint as well as an NFSOFT and its adjoint
   *
   * \f[
   *  f(g_q)=\sum^{N}_{l=0}\sum_{m,n=-l}^l \hat f^{mn}_l D_{mn}^l(\alpha_q,\beta_q,\gamma_q)
   * \f]

   * at the desired bandwidth N for all M random nodes (g_q) with q=0,...,M-1
   *
   * \arg N the bandwidth
   * \arg M the number of nodes
   *
   */

int main(int argc, char **argv)
{
  int N; /**< The bandwidth N  */
  int M; /**< The number of nodes M */

  if (argc < 2)
  {
    printf(
        "This test programm computes the NFSOFT with maximum polynomial degree N at M input rotations\n");
    printf("Usage: simple_test N M \n");
    printf("e.g.: simple_test 8 64\n");
    exit(0);
  }

  /**Read in bandwidth and number of nodes*/
  N = atoi(argv[1]);
  M = atoi(argv[2]);

  printf(
      "computing an NDSOFT, an NFSOFT, an adjoint NDSOFT, and an adjoint NFSOFT\n\n");

  simple_test_nfsoft(N, M);



  /* Exit the program. */
  return EXIT_SUCCESS;

}
