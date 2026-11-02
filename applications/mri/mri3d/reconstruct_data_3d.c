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
#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include "nfft3.h"
/* NFFT_PRECISION_SINGLE, NFFT_PRECISION_DOUBLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

/**
 * \defgroup applications_mri3d_reconstruct_data_3d reconstruct_data_3d
 * \ingroup applications_mri3d
 * \{
 */

/**
 * reconstruct makes an inverse 3d-nfft
 */
static void reconstruct(char* filename,int N,int M,int Z,int iteration, int weight)
{
  int j,k,z,l;                  /* some variables  */
  NFFT_R real,imag;             /* to read the real and imag part of a complex number */
  NFFT(plan) my_plan;            /* plan for the two dimensional nfft  */
  SOLVER(plan_complex) my_iplan;          /* plan for the two dimensional infft */
  FILE* fin;                    /* input file                         */
  FILE* fout_real;              /* output file (real part) */
  FILE* fout_imag;              /* output file (imag part) */
  int my_N[3],my_n[3];          /* to init the nfft */
  NFFT_R epsilon=NFFT_K(0.0000003);     /* tmp to read the obsolent z from 700.acs
                                   epsilon is a the break criterion for
                                   the iteration */
  unsigned infft_flags = CGNR | PRECOMPUTE_DAMP;  /* flags for the infft */

  /* initialise my_plan, specific.
     we don't precompute psi */
  my_N[0]=Z; my_n[0]=NFFT_CEIL(Z*NFFT_K(1.2));
  my_N[1]=N; my_n[1]=NFFT_CEIL(N*NFFT_K(1.2));
  my_N[2]=N; my_n[2]=NFFT_CEIL(N*NFFT_K(1.2));
  NFFT(init_guru)(&my_plan, 3, my_N, M, my_n, 6,
                      PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
                      MALLOC_F| FFTW_INIT,
                      FFTW_MEASURE);

  /* precompute lin psi */
  if(my_plan.flags & PRE_LIN_PSI)
    NFFT(precompute_lin_psi)(&my_plan);

  if (weight)
    infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

  /* initialise my_iplan, advanced */
  SOLVER(init_advanced_complex)(&my_iplan,(NFFT(mv_plan_complex)*)(&my_plan), infft_flags );

  /* get the weights */
  if(my_iplan.flags & PRECOMPUTE_WEIGHT)
  {
    fin=fopen("weights.dat","r");
    for(j=0;j<M;j++)
    {
      fscanf(fin,NFFT__FR__ " ",&my_iplan.w[j]);
    }
    fclose(fin);
  }

  /* get the damping factors */
  if(my_iplan.flags & PRECOMPUTE_DAMP)
  {
    for(j=0;j<N;j++){
      for(k=0;k<N;k++) {
        for(z=0;z<Z;z++) {
        int j2= j-N/2;
        int k2= k-N/2;
        int z2= z-Z/2;
        NFFT_R r=NFFT_SQRT(j2*j2+k2*k2+z2*z2);
        if(r>(NFFT_R) N/2)
          my_iplan.w_hat[z*N*N+j*N+k]=NFFT_K(0.0);
        else
          my_iplan.w_hat[z*N*N+j*N+k]=NFFT_K(1.0);
        }
      }
    }
  }

  /* open the input file */
  fin=fopen(filename,"r");

  /* open the output files */
  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");

  /* read x,y,freal and fimag from the knots */
  for(j=0;j<M;j++)
  {
    fscanf(fin,NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " ",&my_plan.x[3*j+1],&my_plan.x[3*j+2], &my_plan.x[3*j+0],
    &real,&imag);
    my_iplan.y[j] = real + _Complex_I*imag;
  }

  /* precompute psi */
  if(my_plan.flags & PRE_PSI)
    NFFT(precompute_psi)(&my_plan);

  /* precompute full psi */
  if(my_plan.flags & PRE_FULL_PSI)
    NFFT(precompute_full_psi)(&my_plan);

  /* init some guess */
  for(k=0;k<my_plan.N_total;k++)
    my_iplan.f_hat_iter[k]=NFFT_K(0.0);

  /* inverse trafo */
  SOLVER(before_loop_complex)(&my_iplan);
  for(l=0;l<iteration;l++)
  {
    /* break if dot_r_iter is smaller than epsilon*/
    if(my_iplan.dot_r_iter<epsilon)
      break;
    fprintf(stderr,NFFT__FR__ ",  %i of %i\n",NFFT_SQRT(my_iplan.dot_r_iter),
    l+1,iteration);
    SOLVER(loop_one_step_complex)(&my_iplan);
  }

  for(l=0;l<Z;l++)
  {
    for(k=0;k<N*N;k++)
    {
      /* write every Layer in the files */
      fprintf(fout_real,NFFT__FR__ " ",NFFT_CREAL(my_iplan.f_hat_iter[ k+N*N*l ]));
      fprintf(fout_imag,NFFT__FR__ " ",NFFT_CIMAG(my_iplan.f_hat_iter[ k+N*N*l ]));
    }
    fprintf(fout_real,"\n");
    fprintf(fout_imag,"\n");
  }

  fclose(fout_real);
  fclose(fout_imag);

  SOLVER(finalize_complex)(&my_iplan);
  NFFT(finalize)(&my_plan);
}

int main(int argc, char **argv)
{
  if (argc <= 6) {
    printf("usage: ./reconstruct3D FILENAME N M Z ITER WEIGHTS\n");
    return 1;
  }

  reconstruct(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
  return 0;
}
/* \} */
