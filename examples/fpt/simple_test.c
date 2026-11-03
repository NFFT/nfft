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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* It is important to include complex.h before nfft3.h and fftw3.h. */
#include <complex.h>

#include <fftw3.h>

/* NFFT3 header */
#include "nfft3.h"
/* NFFT_PRECISION_DOUBLE, NFFT_PRECISION_SINGLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

int main(void)
{
  /* This example shows the use of the fast polynomial transform to evaluate a
   * finite expansion in Legendre polynomials,
   *
   *   f(x) = a_0 P_0(x) + a_1 P_1(x) + ... + a_N P_N(x)                     (1)
   *
   * at the Chebyshev nodes x_j = cos(j*pi/N), j=0,1,...,N. */
  const int N = 8;

  /* An fpt_set is a data structure that contains precomputed data for a number
   * of different polynomial transforms. Here, we need only one transform. the
   * second parameter (t) is the exponent of the maximum transform size desired
   * (2^t), i.e., t = 3 means that N in (1) can be at most N = 8. */
  FPT(set) set = FPT(init)(1,lrint(ceil(log2((double)N))),0U);

  /* Three-term recurrence coefficients for Legendre polynomials */
  NFFT_R *alpha = malloc((N+2)*sizeof(NFFT_R)),
    *beta = malloc((N+2)*sizeof(NFFT_R)),
    *gamma = malloc((N+2)*sizeof(NFFT_R));

  /* alpha[0] and beta[0] are not referenced. */
  alpha[0] = beta[0] = NFFT_K(0.0);
  /* gamma[0] contains the value of P_0(x) (which is a constant). */
  gamma[0] = NFFT_K(1.0);

  /* Actual three-term recurrence coefficients for Legendre polynomials */
  {
    int k;
    for (k = 0; k <= N; k++)
    {
      alpha[k+1] = ((NFFT_R)(2*k+1))/((NFFT_R)(k+1));
      beta[k+1] = NFFT_K(0.0);
      gamma[k+1] = -((NFFT_R)(k))/((NFFT_R)(k+1));
    }
  }

  printf(
    "Computing a fast polynomial transform (FPT) and a fast discrete cosine \n"
    "transform (DCT) to evaluate\n\n"
    "  f_j = a_0 P_0(x_j) + a_1 P_1(x_j) + ... + a_N P_N(x_j), j=0,1,...,N,\n\n"
    "with N=%d, x_j = cos(j*pi/N), j=0,1,...N, the Chebyshev nodes, a_k,\n"
    "k=0,1,...,N, random Fourier coefficients in [-1,1]x[-1,1]*I, and P_k,\n"
    "k=0,1,...,N, the Legendre polynomials.",N
  );

  /* Random seed, makes things reproducible. */
  NFFT(srand48)(314);

  /* The function fpt_repcompute actually does the precomputation for a single
   * transform. It needs arrays alpha, beta, and gamma, containing the three-
   * term recurrence coefficients, here of the Legendre polynomials. The format
   * is explained above. The sixth parameter (k_start) is where the index in the
   * linear combination (1) starts, here k_start=0. The seventh parameter
   * (kappa) is the threshold which has an influence on the accuracy of the fast
   * polynomial transform. Usually, kappa = 1000 is a good choice. */
  FPT(precompute)(set,0,alpha,beta,gamma,0,NFFT_K(1000.0));


  {
    /* Arrays for Fourier coefficients and function values. */
    NFFT_C *a = malloc((N+1)*sizeof(NFFT_C));
    NFFT_C *b = malloc((N+1)*sizeof(NFFT_C));
    NFFT_R *f = malloc((N+1)*sizeof(NFFT_C));

    /* Plan for discrete cosine transform */
    const int NP1 = N + 1;
    FFTW(r2r_kind) kind = FFTW_REDFT00;
    FFTW(plan) p = FFTW(plan_many_r2r)(1, &NP1, 1, (NFFT_R*)b, NULL, 2, 1,
      (NFFT_R*)f, NULL, 1, 1, &kind, 0U);

    /* random Fourier coefficients */
    {
      int k;
      printf("\n2) Random Fourier coefficients a_k, k=0,1,...,N:\n");
      for (k = 0; k <= N; k++)
      {
        a[k] = NFFT_K(2.0)*NFFT(drand48)() - NFFT_K(1.0); /* for debugging: use k+1 */
        printf("   a_%-2d = %+5.3" NFFT__FES__ "\n",k,NFFT_CREAL(a[k]));
      }
    }

    /* fast polynomial transform */
    FPT(trafo)(set,0,a,b,N,0U);

    /* Renormalize coefficients b_j, j=1,2,...,N-1 owing to how FFTW defines a
     * DCT-I; see
     * http://www.fftw.org/fftw3_doc/1d-Real_002deven-DFTs-_0028DCTs_0029.html
     * for details */
    {
      int j;
      for (j = 1; j < N; j++)
        b[j] *= NFFT_K(0.5);
    }

    /* discrete cosine transform */
    FFTW(execute)(p);

    {
      int j;
      printf("\n3) Function values f_j, j=1,1,...,M:\n");
      for (j = 0; j <= N; j++)
        printf("   f_%-2d = %+5.3" NFFT__FES__ "\n",j,f[j]);
    }

    /* cleanup */
    free(a);
    free(b);
    free(f);

    /* cleanup */
    FFTW(destroy_plan)(p);
  }

  /* cleanup */
  FPT(finalize)(set);
  free(alpha);
  free(beta);
  free(gamma);

  return EXIT_SUCCESS;
}
