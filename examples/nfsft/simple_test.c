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

/* standard headers */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
/* It is important to include complex.h before nfft3.h. */
#include <complex.h>

#include "nfft3.h" /* NFFT3 header */

#define __FES__ "E"
#define K(x) ((double) x)

static void simple_test_nfsft(void)
{
  const int N = 4; /* bandwidth/maximum degree */
  const int M = 8; /* number of nodes */
  nfsft_plan plan; /* transform plan */
  int j, k, n; /* loop variables */

  /* precomputation (for fast polynomial transform) */
  nfsft_precompute(N,1000.0,0U,0U);

  /* Initialize transform plan using the guru interface. All input and output
   * arrays are allocated by nfsft_init_guru(). Computations are performed with
   * respect to L^2-normalized spherical harmonics Y_k^n. The array of spherical
   * Fourier coefficients is preserved during transformations. The NFFT uses a
   * cut-off parameter m = 6. See the NFFT 3 manual for details.
   */
  nfsft_init_guru(&plan, N, M, NFSFT_MALLOC_X | NFSFT_MALLOC_F |
    NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED | NFSFT_PRESERVE_F_HAT,
    PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFT_OUT_OF_PLACE, 6);

  /* pseudo-random nodes */
  for (j = 0; j < plan.M_total; j++)
  {
    plan.x[2*j]= nfft_drand48() - K(0.5);
    plan.x[2*j+1]= K(0.5) * nfft_drand48();
  }

  /* precomputation (for NFFT, node-dependent) */
  nfsft_precompute_x(&plan);

  /* pseudo-random Fourier coefficients */
  for (k = 0; k <= plan.N; k++)
    for (n = -k; n <= k; n++)
      plan.f_hat[NFSFT_INDEX(k,n,&plan)] =
          nfft_drand48() - K(0.5) + _Complex_I*(nfft_drand48() - K(0.5));

  /* Direct transformation, display result. */
  nfsft_trafo_direct(&plan);
  printf("Vector f (NDSFT):\n");
  for (j = 0; j < plan.M_total; j++)
    printf("f[%+2d] = %+5.3" __FES__ " %+5.3" __FES__ "*I\n",j,
      creal(plan.f[j]), cimag(plan.f[j]));

  printf("\n");

  /* Fast approximate transformation, display result. */
  nfsft_trafo(&plan);
  printf("Vector f (NFSFT):\n");
  for (j = 0; j < plan.M_total; j++)
    printf("f[%+2d] = %+5.3" __FES__ " %+5.3" __FES__ "*I\n",j,
      creal(plan.f[j]), cimag(plan.f[j]));

  printf("\n");

  /* Direct adjoint transformation, display result. */
  nfsft_adjoint_direct(&plan);
  printf("Vector f_hat (NDSFT):\n");
  for (k = 0; k <= plan.N; k++)
    for (n = -k; n <= k; n++)
      fprintf(stdout,"f_hat[%+2d,%+2d] = %+5.3" __FES__ " %+5.3" __FES__ "*I\n",k,n,
        creal(plan.f_hat[NFSFT_INDEX(k,n,&plan)]),
        cimag(plan.f_hat[NFSFT_INDEX(k,n,&plan)]));

  printf("\n");

  /* Fast approximate adjoint transformation, display result. */
  nfsft_adjoint(&plan);
  printf("Vector f_hat (NFSFT):\n");
  for (k = 0; k <= plan.N; k++)
  {
    for (n = -k; n <= k; n++)
    {
      fprintf(stdout,"f_hat[%+2d,%+2d] = %+5.3" __FES__ " %+5.3" __FES__ "*I\n",k,n,
        creal(plan.f_hat[NFSFT_INDEX(k,n,&plan)]),
        cimag(plan.f_hat[NFSFT_INDEX(k,n,&plan)]));
    }
  }

  /* Finalize the plan. */
  nfsft_finalize(&plan);

  /* Destroy data precomputed for fast polynomial transform. */
	nfsft_forget();
}

int main(void)
{
  printf("Computing an NDSFT, an NFSFT, an adjoint NDSFT, and an adjoint NFSFT"
    "...\n\n");
  simple_test_nfsft();
  return EXIT_SUCCESS;
}
