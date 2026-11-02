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
 * \defgroup applications_mri3d_construct_data_1d2d construct_data_1d2d
 * \ingroup applications_mri3d
 * \{
 */

/**
 * construct makes an 2d-nfft for every slice
 */
static void construct(char * file, int N, int M, int Z, NFFT_C *mem)
{
  int j,z;                /* some variables */
  NFFT_R tmp;             /* a placeholder */
  NFFT(plan) my_plan;      /* plan for the two dimensional nfft  */
  FILE* fp;

  /* initialise my_plan */
  NFFT(init_2d)(&my_plan,N,N,M/Z);

  fp=fopen("knots.dat","r");

  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__,&my_plan.x[2*j+0],&my_plan.x[2*j+1],&tmp);
  }
  fclose(fp);

  fp=fopen(file,"w");

  for(z=0;z<Z;z++) {
    tmp = (NFFT_R) z;

    for(j=0;j<N*N;j++)
      my_plan.f_hat[j] = mem[(z*N*N+N*N*Z/2+j)%(N*N*Z)];

    if(my_plan.flags & PRE_PSI)
      NFFT(precompute_psi)(&my_plan);

    NFFT(trafo)(&my_plan);

    for(j=0;j<my_plan.M_total;j++)
    {
      fprintf(fp,NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ "\n",my_plan.x[2*j+0],my_plan.x[2*j+1],tmp/Z-NFFT_K(0.5),
              NFFT_CREAL(my_plan.f[j]),NFFT_CIMAG(my_plan.f[j]));
    }
  }
  fclose(fp);

  NFFT(finalize)(&my_plan);
}

/**
 * fft makes an 1D-ftt for every knot through
 * all layers
 */
static void fft(int N,int M,int Z, NFFT_C *mem)
{
  FFTW(plan) plan;
  plan = FFTW(plan_many_dft)(1, &Z, N*N,
                                  mem, NULL,
                                  N*N, 1,
                                  mem, NULL,
                                  N*N,1 ,
                                  FFTW_FORWARD, FFTW_ESTIMATE);

  FFTW(execute)(plan); /* execute the fft */
  FFTW(destroy_plan)(plan);
}

/**
 * read fills the memory with the file input_data_f.dat as
 * the real part of f and with zeros for the imag part of f
 */
static void read_data(int N,int M,int Z, NFFT_C *mem)
{
  int i,z;
  NFFT_R real;
  FILE* fin;
  fin=fopen("input_f.dat","r");

  for(z=0;z<Z;z++)
  {
    for(i=0;i<N*N;i++)
    {
      fscanf(fin,NFFT__FR__ " ",&real );
      mem[(z*N*N+N*N*Z/2+i)%(N*N*Z)]=real;
    }
  }
  fclose(fin);
}

int main(int argc, char **argv)
{
  NFFT_C *mem;

  if (argc <= 4) {
    printf("usage: ./construct_data FILENAME N M Z\n");
    return 1;
  }

  mem = (NFFT_C*) NFFT(malloc)(sizeof(NFFT_C) * atoi(argv[2]) * atoi(argv[2]) * atoi(argv[4]));

  read_data(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]), mem);

  fft(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]), mem);

  construct(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]), mem);

  NFFT(free)(mem);

  return 0;
}
/* \} */
