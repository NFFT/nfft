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
 * \defgroup applications_mri3d_reconstruct_data_1d2d reconstruct_data_1d2d
 * \ingroup applications_mri3d
 * \{
 */

/**
 * reconstruct makes an inverse 2d-nfft for every slice
 */
static void reconstruct(char* filename,int N,int M,int Z,int iteration, int weight, fftw_complex *mem)
{
  int j,k,l,z;                  /* some variables  */
  double real,imag;             /* to read the real and imag part of a complex number */
  nfft_plan my_plan;            /* plan for the two dimensional nfft  */
  solver_plan_complex my_iplan; /* plan for the two dimensional infft */
  FILE* fin;                    /* input file */
  int my_N[2],my_n[2];          /* to init the nfft */
  double tmp, epsilon=0.0000003;/* tmp to read the obsolent z from the input file
                                   epsilon is the break criterium for
                                   the iteration */
  unsigned infft_flags = CGNR | PRECOMPUTE_DAMP; /* flags for the infft */

  /* initialise my_plan */
  my_N[0]=N;my_n[0]=ceil(N*1.2);
  my_N[1]=N; my_n[1]=ceil(N*1.2);
  nfft_init_guru(&my_plan, 2, my_N, M/Z, my_n, 6, PRE_PHI_HUT| PRE_PSI|
                         MALLOC_X| MALLOC_F_HAT| MALLOC_F|
                        FFTW_INIT,
                        FFTW_MEASURE);

  /* precompute lin psi if set */
  if(my_plan.flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&my_plan);

  /* set the flags for the infft*/
  if (weight)
    infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

  /* initialise my_iplan, advanced */
  solver_init_advanced_complex(&my_iplan,(nfft_mv_plan_complex*)(&my_plan), infft_flags );

  /* get the weights */
  if(my_iplan.flags & PRECOMPUTE_WEIGHT)
  {
    fin=fopen("weights.dat","r");
    for(j=0;j<my_plan.M_total;j++)
    {
        fscanf(fin,"%le ",&my_iplan.w[j]);
    }
    fclose(fin);
  }

  /* get the damping factors */
  if(my_iplan.flags & PRECOMPUTE_DAMP)
  {
    for(j=0;j<N;j++){
      for(k=0;k<N;k++) {
        int j2= j-N/2;
        int k2= k-N/2;
        double r=sqrt(j2*j2+k2*k2);
        if(r>(double) N/2)
          my_iplan.w_hat[j*N+k]=0.0;
        else
          my_iplan.w_hat[j*N+k]=1.0;
      }
    }
  }

  /* open the input file */
  fin=fopen(filename,"r");

  /* For every Layer*/
  for(z=0;z<Z;z++) {

    /* read x,y,freal and fimag from the knots */
    for(j=0;j<my_plan.M_total;j++)
    {
      fscanf(fin,"%le %le %le %le %le ",&my_plan.x[2*j+0],&my_plan.x[2*j+1], &tmp,
      &real,&imag);
      my_iplan.y[j] = real + _Complex_I*imag;
    }

    /* precompute psi if set just one time because the knots equal each plane */
    if(z==0 && my_plan.flags & PRE_PSI)
      nfft_precompute_psi(&my_plan);

    /* precompute full psi if set just one time because the knots equal each plane */
    if(z==0 && my_plan.flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&my_plan);

    /* init some guess */
    for(k=0;k<my_plan.N_total;k++)
      my_iplan.f_hat_iter[k]=0.0;

    /* inverse trafo */
    solver_before_loop_complex(&my_iplan);
    for(l=0;l<iteration;l++)
    {
      /* break if dot_r_iter is smaller than epsilon*/
      if(my_iplan.dot_r_iter<epsilon)
      break;
      fprintf(stderr,"%e,  %i of %i\n",sqrt(my_iplan.dot_r_iter),
      iteration*z+l+1,iteration*Z);
      solver_loop_one_step_complex(&my_iplan);
    }
    for(k=0;k<my_plan.N_total;k++) {
      /* write every slice in the memory.
      here we make an fftshift direct */
      mem[(Z*N*N/2+z*N*N+ k)%(Z*N*N)] = my_iplan.f_hat_iter[k];
    }
  }

  fclose(fin);

  /* finalize the infft */
  solver_finalize_complex(&my_iplan);

  /* finalize the nfft */
  nfft_finalize(&my_plan);
}

/**
 * print writes the memory back in a file
 * output_real.dat for the real part and output_imag.dat for the imaginary part
 */
static void print(int N,int M,int Z, fftw_complex *mem)
{
  int i,j;
  FILE* fout_real;
  FILE* fout_imag;
  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");

  for(i=0;i<Z;i++) {
    for (j=0;j<N*N;j++) {
      fprintf(fout_real,"%le ",creal(mem[(Z*N*N/2+i*N*N+ j)%(Z*N*N)]) /Z);
      fprintf(fout_imag,"%le ",cimag(mem[(Z*N*N/2+i*N*N+ j)%(Z*N*N)]) /Z);
    }
    fprintf(fout_real,"\n");
    fprintf(fout_imag,"\n");
  }

  fclose(fout_real);
  fclose(fout_imag);
}

int main(int argc, char **argv)
{
  fftw_complex *mem;
  fftw_plan plan;
  int N,M,Z;

  if (argc <= 6) {
    printf("usage: ./reconstruct FILENAME N M Z ITER WEIGHTS\n");
    return 1;
  }

  N=atoi(argv[2]);
  M=atoi(argv[3]);
  Z=atoi(argv[4]);

  /* Allocate memory to hold every layer in memory after the
  2D-infft */
  mem = (fftw_complex*) nfft_malloc(sizeof(fftw_complex) * atoi(argv[2]) * atoi(argv[2]) * atoi(argv[4]));

  /* Create plan for the 1d-ifft */
  plan = fftw_plan_many_dft(1, &Z, N*N,
                                  mem, NULL,
                                  N*N, 1,
                                  mem, NULL,
                                  N*N,1 ,
                                  FFTW_BACKWARD, FFTW_MEASURE);

  /* execute the 2d-infft's */
  reconstruct(argv[1],N,M,Z,atoi(argv[5]),atoi(argv[6]),mem);

  /* execute the 1d-fft's */
  fftw_execute(plan);

  /* write the memory back in files */
  print(N,M,Z, mem);

  /* free memory */
  nfft_free(mem);
  fftw_destroy_plan(plan);
  return 1;
}
/* \} */
