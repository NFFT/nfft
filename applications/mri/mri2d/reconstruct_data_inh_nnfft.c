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
#include <limits.h>
#include <complex.h>

#include "nfft3.h"
/* NFFT_PRECISION_SINGLE, NFFT_PRECISION_DOUBLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

/**
 * \defgroup applications_mri2d_construct_data_inh_nnfft construct_data_inh_nnfft
 * \ingroup applications_mri2d
 * \{
 */

/**
 * reconstruct
 */
static void reconstruct(char* filename,int N,int M,int iteration, int weight)
{
  int j,k,l;                    /* some variables  */
  NNFFT(plan) my_plan;            /* plan for the two dimensional nfft  */
  SOLVER(plan_complex) my_iplan;          /* plan for the two dimensional infft */
  FILE* fin;                    /* input file                         */
  FILE* finh;
  FILE* ftime;
  FILE* fout_real;              /* output file                        */
  FILE* fout_imag;              /* output file                        */
  int my_N[3],my_n[3];          /* to init the nfft */
  NFFT_R t0, t1;
  NFFT_R t,epsilon=NFFT_K(0.0000003);     /* epsilon is a the break criterium for
                                   the iteration */
  unsigned infft_flags = CGNR | PRECOMPUTE_DAMP; /* flags for the infft*/
  NFFT_R time,min_time,max_time,min_inh,max_inh;
  NFFT_R real,imag;
  NFFT_R *w;

  NFFT_R Ts;
  NFFT_R W;
  int N3;
  int m=2;
  NFFT_R sigma = NFFT_K(1.25);

  w = (NFFT_R*)NFFT(malloc)(sizeof(NFFT_R)*N*N);

  ftime=fopen("readout_time.dat","r");
  finh=fopen("inh.dat","r");

  min_time=INT_MAX; max_time=INT_MIN;
  for(j=0;j<M;j++)
  {
    fscanf(ftime,NFFT__FR__ " ",&time);
    if(time<min_time)
      min_time = time;
    if(time>max_time)
      max_time = time;
  }

  fclose(ftime);

  Ts=(min_time+max_time)/NFFT_K(2.0);

  min_inh=INT_MAX; max_inh=INT_MIN;
  for(j=0;j<N*N;j++)
  {
    fscanf(finh,NFFT__FR__ " ",&w[j]);
    if(w[j]<min_inh)
      min_inh = w[j];
    if(w[j]>max_inh)
      max_inh = w[j];
  }
  fclose(finh);

  N3=NFFT_CEIL((MAX(NFFT_FABS(min_inh),NFFT_FABS(max_inh))*(max_time-min_time)/NFFT_K(2.0))*4);


  W=MAX(NFFT_FABS(min_inh),NFFT_FABS(max_inh))*NFFT_K(2.0);

  fprintf(stderr,"3:  %i " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ "\n",N3,W,min_inh,max_inh,min_time,max_time,Ts);

  /* initialise my_plan */
  my_N[0]=N;my_n[0]=NFFT_CEIL(N*sigma);
  my_N[1]=N; my_n[1]=NFFT_CEIL(N*sigma);
  my_N[2]=N3; my_n[2]=NFFT_CEIL(N3*sigma);
  NNFFT(init_guru)(&my_plan, 3, N*N, M, my_N,my_n,m,
        PRE_PSI| PRE_PHI_HUT| MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F );

  /* precompute lin psi if set */
  if(my_plan.nnfft_flags & PRE_LIN_PSI)
    NNFFT(precompute_lin_psi)(&my_plan);

  /* set the flags for the infft*/
  if (weight)
    infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

  /* initialise my_iplan, advanced */
  SOLVER(init_advanced_complex)(&my_iplan,(NFFT(mv_plan_complex)*)(&my_plan), infft_flags );

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
  ftime=fopen("readout_time.dat","r");

  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fin,NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " ",&my_plan.x[3*j+0],&my_plan.x[3*j+1],&real,&imag);
    my_iplan.y[j]=real+ _Complex_I*imag;
    fscanf(ftime,NFFT__FR__ " ",&my_plan.x[3*j+2]);

    my_plan.x[3*j+2] = (my_plan.x[3*j+2]-Ts)*W/N3;
  }

  for(j=0;j<N;j++)
    {
    for(l=0;l<N;l++)
      {
        my_plan.v[3*(N*j+l)+0]=(((NFFT_R) j) -(((NFFT_R) N)/NFFT_K(2.0)))/((NFFT_R) N);
        my_plan.v[3*(N*j+l)+1]=(((NFFT_R) l) -(((NFFT_R) N)/NFFT_K(2.0)))/((NFFT_R) N);
        my_plan.v[3*(N*j+l)+2] = w[N*j+l]/W ;
      }
    }

  /* precompute psi */
  if(my_plan.nnfft_flags & PRE_PSI) {
    NNFFT(precompute_psi)(&my_plan);
    if(my_plan.nnfft_flags & PRE_FULL_PSI)
      NNFFT(precompute_full_psi)(&my_plan);
  }

  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    NNFFT(precompute_phi_hut)(&my_plan);

  /* init some guess */
  for(k=0;k<my_plan.N_total;k++)
  {
    my_iplan.f_hat_iter[k]=NFFT_K(0.0);
  }

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
  t = t1-t0;

  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");

  for(k=0;k<my_plan.N_total;k++) {

    my_iplan.f_hat_iter[k]*=NFFT_CEXP(-NFFT_K(2.0)*_Complex_I*NFFT_KPI*Ts*w[k]);

    fprintf(fout_real,NFFT__FR__ " ", NFFT_CREAL(my_iplan.f_hat_iter[k]));
    fprintf(fout_imag,NFFT__FR__ " ", NFFT_CIMAG(my_iplan.f_hat_iter[k]));
  }


  fclose(fout_real);
  fclose(fout_imag);


  /* finalize the infft */
  SOLVER(finalize_complex)(&my_iplan);

  /* finalize the nfft */
  NNFFT(finalize)(&my_plan);

  NFFT(free)(w);
}

int main(int argc, char **argv)
{
  if (argc <= 5) {
    printf("usage: ./reconstruct_data_inh_nnfft FILENAME N M ITER WEIGHTS\n");
    return 1;
  }

  reconstruct(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));

  return 0;
}
/* \} */
