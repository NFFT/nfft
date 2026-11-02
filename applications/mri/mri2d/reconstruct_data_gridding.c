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
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "nfft3.h"
/* NFFT_PRECISION_SINGLE, NFFT_PRECISION_DOUBLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

/**
 * \defgroup applications_mri2d_construct_data_gridding construct_data_gridding
 * \ingroup applications_mri2d
 * \{
 */

/**
 * reconstruct makes a 2d-adjoint-nfft
 */
static void reconstruct(char* filename, int N, int M, int weight)
{
  int j;                   /* some variables  */
  NFFT_R weights;          /* store one weight temporary */
  NFFT_R real,imag;        /* to read the real and imag part of a complex number */
  NFFT(plan) my_plan;       /* plan for the two dimensional nfft  */
  FILE* fin;               /* input file  */
  FILE* fweight;           /* input file for the weights */
  FILE *fout_real;         /* output file  */
  FILE *fout_imag;         /* output file  */
  int my_N[2],my_n[2];
  int flags = PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
                      MALLOC_F| FFTW_INIT| FFTW_MEASURE;

  /* initialise nfft */
  my_N[0]=N; my_n[0]=NFFT_CEIL(N*NFFT_K(1.2));
  my_N[1]=N; my_n[1]=NFFT_CEIL(N*NFFT_K(1.2));
  NFFT(init_guru)(&my_plan, 2, my_N, M, my_n, 6,flags,
                      FFTW_MEASURE| FFTW_DESTROY_INPUT);

  fin=fopen(filename,"r");

  fweight=fopen("weights.dat","r");
  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fweight,NFFT__FR__ " ",&weights);
    fscanf(fin,NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__,&my_plan.x[2*j+0],&my_plan.x[2*j+1],&real,&imag);
    my_plan.f[j] = real + _Complex_I*imag;
    if (weight)
      my_plan.f[j] = my_plan.f[j] * weights;
  }
  fclose(fweight);

  /* precompute psi */
  if(my_plan.flags & PRE_PSI)
    NFFT(precompute_psi)(&my_plan);

  /* precompute full psi */
  if(my_plan.flags & PRE_FULL_PSI)
    NFFT(precompute_full_psi)(&my_plan);


  /* compute the adjoint nfft */
  NFFT(adjoint)(&my_plan);

  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");

  for (j=0;j<N*N;j++) {
    fprintf(fout_real,NFFT__FR__ " ",NFFT_CREAL(my_plan.f_hat[j]));
    fprintf(fout_imag,NFFT__FR__ " ",NFFT_CIMAG(my_plan.f_hat[j]));
  }

  fclose(fin);
  fclose(fout_real);
  fclose(fout_imag);

  NFFT(finalize)(&my_plan);
}


int main(int argc, char **argv)
{
  if (argc <= 5) {
    printf("usage: ./reconstruct_data_gridding FILENAME N M ITER WEIGHTS\n");
    return 1;
  }
  reconstruct(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[5]));

  return 0;
}
/* \} */
