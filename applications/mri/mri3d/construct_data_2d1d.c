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

#include <stdlib.h>
#include <math.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3.h"

/**
 * \defgroup applications_mri3d_construct_data_1d2d construct_data_1d2d
 * \ingroup applications_mri3d
 * \{
 */

/**
 * construct makes an 2d-nfft for every slice
 */
static void construct(char * file, int N, int M, int Z, fftw_complex *mem)
{
  int j,z;                /* some variables */
  double tmp;             /* a placeholder */
  nfft_plan my_plan;      /* plan for the two dimensional nfft  */
  FILE* fp;

  /* initialise my_plan */
  nfft_init_2d(&my_plan,N,N,M/Z);

  fp=fopen("knots.dat","r");

  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,"%le %le %le",&my_plan.x[2*j+0],&my_plan.x[2*j+1],&tmp);
  }
  fclose(fp);

  fp=fopen(file,"w");

  for(z=0;z<Z;z++) {
    tmp = (double) z;

    for(j=0;j<N*N;j++)
      my_plan.f_hat[j] = mem[(z*N*N+N*N*Z/2+j)%(N*N*Z)];

    if(my_plan.flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);

    nfft_trafo(&my_plan);

    for(j=0;j<my_plan.M_total;j++)
    {
      fprintf(fp,"%le %le %le %le %le\n",my_plan.x[2*j+0],my_plan.x[2*j+1],tmp/Z-0.5,
              creal(my_plan.f[j]),cimag(my_plan.f[j]));
    }
  }
  fclose(fp);

  nfft_finalize(&my_plan);
}

/**
 * fft makes an 1D-ftt for every knot through
 * all layers
 */
static void fft(int N,int M,int Z, fftw_complex *mem)
{
  fftw_plan plan;
  plan = fftw_plan_many_dft(1, &Z, N*N,
                                  mem, NULL,
                                  N*N, 1,
                                  mem, NULL,
                                  N*N,1 ,
                                  FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(plan); /* execute the fft */
  fftw_destroy_plan(plan);
}

/**
 * read fills the memory with the file input_data_f.dat as
 * the real part of f and with zeros for the imag part of f
 */
static void read_data(int N,int M,int Z, fftw_complex *mem)
{
  int i,z;
  double real;
  FILE* fin;
  fin=fopen("input_f.dat","r");

  for(z=0;z<Z;z++)
  {
    for(i=0;i<N*N;i++)
    {
      fscanf(fin,"%le ",&real );
      mem[(z*N*N+N*N*Z/2+i)%(N*N*Z)]=real;
    }
  }
  fclose(fin);
}

int main(int argc, char **argv)
{
  fftw_complex *mem;

  if (argc <= 4) {
    printf("usage: ./construct_data FILENAME N M Z\n");
    return 1;
  }

  mem = (fftw_complex*) nfft_malloc(sizeof(fftw_complex) * atoi(argv[2]) * atoi(argv[2]) * atoi(argv[4]));

  read_data(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]), mem);

  fft(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]), mem);

  construct(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]), mem);

  nfft_free(mem);

  return 1;
}
/* \} */
