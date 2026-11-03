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
/* NFFT_PRECISION_DOUBLE, NFFT_PRECISION_SINGLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

/**
 * \defgroup applications_mri3d_reconstruct_data_gridding reconstruct_data_gridding
 * \ingroup applications_mri3d
 * \{
 */

/**
 * reconstruct makes an 2d-adjoint-nfft for every slice
 */
static void reconstruct(char* filename,int N,int M,int Z, int weight ,NFFT_C *mem)
{
  int j,k,z;               /* some variables  */
  NFFT_R weights;          /* store one weight temporary */
  NFFT_R tmp;              /* tmp to read the obsolent z from the input file */
  NFFT_R real,imag;        /* to read the real and imag part of a complex number */
  NFFT(plan) my_plan;       /* plan for the two dimensional nfft  */
  int my_N[2],my_n[2];     /* to init the nfft */
  FILE* fin;               /* input file  */
  FILE* fweight;           /* input file for the weights */

  /* initialise my_plan */
  my_N[0]=N; my_n[0]=NFFT_CEIL(N*NFFT_K(1.2));
  my_N[1]=N; my_n[1]=NFFT_CEIL(N*NFFT_K(1.2));
  NFFT(init_guru)(&my_plan, 2, my_N, M/Z, my_n, 6, PRE_PHI_HUT| PRE_PSI|
                        MALLOC_X| MALLOC_F_HAT| MALLOC_F|
                        FFTW_INIT,
                        FFTW_MEASURE);

  /* precompute lin psi if set */
  if(my_plan.flags & PRE_LIN_PSI)
    NFFT(precompute_lin_psi)(&my_plan);

  fin=fopen(filename,"r");

  for(z=0;z<Z;z++) {
    fweight=fopen("weights.dat","r");
    for(j=0;j<my_plan.M_total;j++)
    {
      fscanf(fweight,NFFT__FR__ " ",&weights);
      fscanf(fin,NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__,
             &my_plan.x[2*j+0],&my_plan.x[2*j+1],&tmp,&real,&imag);
      my_plan.f[j] = real + _Complex_I*imag;
      if(weight)
        my_plan.f[j] = my_plan.f[j] * weights;
    }
    fclose(fweight);

    /* precompute psi if set just one time because the knots equal each slice */
    if(z==0 && my_plan.flags & PRE_PSI)
      NFFT(precompute_psi)(&my_plan);

    /* precompute full psi if set just one time because the knots equal each slice */
    if(z==0 && my_plan.flags & PRE_FULL_PSI)
      NFFT(precompute_full_psi)(&my_plan);

    /* compute the adjoint nfft */
    NFFT(adjoint)(&my_plan);

    for(k=0;k<my_plan.N_total;k++) {
      /* write every slice in the memory.
      here we make an fftshift direct */
      mem[(Z*N*N/2+z*N*N+ k)%(Z*N*N)] = my_plan.f_hat[k];
    }
  }
  fclose(fin);

  NFFT(finalize)(&my_plan);
}

 /**
 * print writes the memory back in a file
 * output_real.dat for the real part and output_imag.dat for the imaginary part
 */
static void print(int N,int M,int Z, NFFT_C *mem)
{
  int i,j;
  FILE* fout_real;
  FILE* fout_imag;
  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");

  for(i=0;i<Z;i++) {
    for (j=0;j<N*N;j++) {
      fprintf(fout_real,NFFT__FR__ " ",NFFT_CREAL(mem[(Z*N*N/2+i*N*N+ j)%(Z*N*N)]) /(NFFT_R)Z);
      fprintf(fout_imag,NFFT__FR__ " ",NFFT_CIMAG(mem[(Z*N*N/2+i*N*N+ j)%(Z*N*N)]) /(NFFT_R)Z);
    }
    fprintf(fout_real,"\n");
    fprintf(fout_imag,"\n");
  }

  fclose(fout_real);
  fclose(fout_imag);
}


int main(int argc, char **argv)
{
  NFFT_C *mem;
  FFTW(plan) plan;
  int N,M,Z;

  if (argc <= 6) {
    printf("usage: ./reconstruct_data_gridding FILENAME N M Z ITER WEIGHTS\n");
    return 1;
  }

  N=atoi(argv[2]);
  M=atoi(argv[3]);
  Z=atoi(argv[4]);

  /* Allocate memory to hold every slice in memory after the
  2D-infft */
  mem = (NFFT_C*) NFFT(malloc)(sizeof(NFFT_C) * atoi(argv[2]) * atoi(argv[2]) * atoi(argv[4]));

  /* Create plan for the 1d-ifft */
  plan = FFTW(plan_many_dft)(1, &Z, N*N,
                                  mem, NULL,
                                  N*N, 1,
                                  mem, NULL,
                                  N*N,1 ,
                                  FFTW_BACKWARD, FFTW_MEASURE);

  /* execute the 2d-nfft's */
  reconstruct(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[6]),mem);

  /* execute the 1d-fft's */
  FFTW(execute)(plan);

  /* write the memory back in files */
  print(N,M,Z, mem);

  /* free memory */
  NFFT(free)(mem);

  return 0;
}
/* \} */
