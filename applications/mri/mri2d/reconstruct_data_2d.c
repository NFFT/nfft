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
/* NFFT_PRECISION_DOUBLE, NFFT_PRECISION_SINGLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

/**
 * \defgroup applications_mri2d_reconstruct_data_2d reconstruct_data_2d
 * \ingroup applications_mri2d
 * \{
 */

/**
 * reconstruct makes an inverse 2d nfft
 */
static void reconstruct(char* filename,int N,int M,int iteration, int weight)
{
  int j,k,l;                    /* some variables  */
  double t0, t1, t;
  NFFT_R real,imag;             /* to read the real and imag part of a complex number */
  NFFT(plan) my_plan;            /* plan for the two dimensional nfft  */
  SOLVER(plan_complex) my_iplan; /* plan for the two dimensional infft */
  FILE* fin;                    /* input file                         */
  FILE* fout_real;              /* output file                        */
  FILE* fout_imag;              /* output file                        */
  int my_N[2],my_n[2];          /* to init the nfft */
  NFFT_R epsilon=NFFT_K(0.0000003);     /* epsilon is a the break criterium for
                                   the iteration */
  unsigned infft_flags = CGNR | PRECOMPUTE_DAMP;  /* flags for the infft*/
  int m = 6;
  NFFT_R alpha = NFFT_K(2.0);
  /* initialise my_plan */
  my_N[0]=N; my_n[0]=NFFT_CEIL(N*alpha);
  my_N[1]=N; my_n[1]=NFFT_CEIL(N*alpha);
  NFFT(init_guru)(&my_plan, 2, my_N, M, my_n, m, PRE_PHI_HUT| PRE_PSI|
                         MALLOC_X| MALLOC_F_HAT| MALLOC_F|
                         FFTW_INIT,
                         FFTW_MEASURE);

  /* set the flags for the infft*/
  if (weight)
    infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

  /* initialise my_iplan, advanced */
  SOLVER(init_advanced_complex)(&my_iplan,(NFFT(mv_plan_complex)*)&my_plan, infft_flags );

  /* get the weights */
  if(my_iplan.flags & PRECOMPUTE_WEIGHT)
  {
    fin=fopen("weights.dat","r");
    for(j=0;j<my_plan.M_total;j++)
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
        int j2= j-N/2;
        int k2= k-N/2;
        NFFT_R r=NFFT_SQRT(j2*j2+k2*k2);
        if(r>(NFFT_R) N/2)
          my_iplan.w_hat[j*N+k]=NFFT_K(0.0);
        else
          my_iplan.w_hat[j*N+k]=NFFT_K(1.0);
      }
    }
  }

  /* open the input file */
  fin=fopen(filename,"r");

  /* read x,y,freal and fimag from the knots */
  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fin,NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " ",&my_plan.x[2*j+0],&my_plan.x[2*j+1],
    &real,&imag);
    my_iplan.y[j] = real + _Complex_I*imag;
  }

	fclose(fin);

  NFFT(precompute_one_psi)(&my_plan);

  /* init some guess */
  for(k=0;k<my_plan.N_total;k++)
    my_iplan.f_hat_iter[k]=NFFT_K(0.0);

  t0 = NFFT(clock_gettime_seconds)();

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


  t1 = NFFT(clock_gettime_seconds)();
  t=t1-t0;

  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");

  for(k=0;k<my_plan.N_total;k++) {
    fprintf(fout_real,NFFT__FR__ " ", NFFT_CREAL(my_iplan.f_hat_iter[k]));
    fprintf(fout_imag,NFFT__FR__ " ", NFFT_CIMAG(my_iplan.f_hat_iter[k]));
  }

  fclose(fout_real);
  fclose(fout_imag);

  /* finalize the infft */
  SOLVER(finalize_complex)(&my_iplan);

  /* finalize the nfft */
  NFFT(finalize)(&my_plan);
}

int main(int argc, char **argv)
{
  if (argc <= 5) {
    printf("usage: ./reconstruct_data_2d FILENAME N M ITER WEIGHTS\n");
    return 1;
  }

  reconstruct(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));

  return 0;
}
/* \} */
