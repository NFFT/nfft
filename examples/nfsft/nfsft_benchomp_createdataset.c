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

#include "config.h"

#include "nfft3.h"
#include "infft.h"

void nfsft_benchomp_createdataset(unsigned int trafo_adjoint, int N, int M)
{
  int t, j, k, n;
  R *x;
  C *f, *f_hat;
  int N_total = (2*N+2) * (2*N+2);
  nfsft_plan ptemp;

  nfsft_init_guru(&ptemp, N, M, NFSFT_MALLOC_X | NFSFT_MALLOC_F |
    NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED | NFSFT_PRESERVE_F_HAT,
    PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFT_OUT_OF_PLACE, 6);

  x = (R*) nfft_malloc(2*M*sizeof(R));
  f = (C*) nfft_malloc(M*sizeof(C));
  f_hat = (C*) nfft_malloc(N_total*sizeof(C));

  /* init pseudo-random nodes */
  for (j = 0; j < M; j++)
  {
    x[2*j]= X(drand48)() - K(0.5);
    x[2*j+1]= K(0.5) * X(drand48)();
  }
 
  if (trafo_adjoint==0)
  {
    for (k = 0; k <= N; k++)
      for (n = -k; n <= k; n++)
        nfft_vrand_unit_complex(f_hat+NFSFT_INDEX(k,n,&ptemp),1);
  }
  else
  {
    nfft_vrand_unit_complex(f,M);
  }

  printf("%d %d %d\n", trafo_adjoint, N, M);

  for (j=0; j < M; j++)
  {
    for (t=0; t < 2; t++)
      printf("%.16e ", x[2*j+t]);
    printf("\n");
  }

  if (trafo_adjoint==0)
  {
    for (k = 0; k <= N; k++)
      for (n = -k; n <= k; n++)
        printf("%.16e %.16e\n", creal(f_hat[NFSFT_INDEX(k,n,&ptemp)]), cimag(f_hat[NFSFT_INDEX(k,n,&ptemp)]));
  }
  else
  {
    for (j=0; j < M; j++)
      printf("%.16e %.16e\n", creal(f[j]), cimag(f[j]));
  }

  nfft_free(x);
  nfft_free(f);
  nfft_free(f_hat);
}

int main(int argc, char **argv)
{
  int trafo_adjoint;
  int N;
  int M;

  if (argc < 4) {
    fprintf(stderr, "usage: tr_adj N M\n");
    return -1;
  }

  trafo_adjoint = atoi(argv[1]);
  if (trafo_adjoint < 0 && trafo_adjoint > 1)
    trafo_adjoint = 1;

  N = atoi(argv[2]);
  M = atoi(argv[3]);
  fprintf(stderr, "tr_adj=%d, N=%d, M=%d\n", trafo_adjoint, N, M);

  nfsft_benchomp_createdataset(trafo_adjoint, N, M);

  return 0;
}
