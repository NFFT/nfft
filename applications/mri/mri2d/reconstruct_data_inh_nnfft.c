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
  nnfft_plan my_plan;            /* plan for the two dimensional nfft  */
  solver_plan_complex my_iplan;          /* plan for the two dimensional infft */
  FILE* fin;                    /* input file                         */
  FILE* finh;
  FILE* ftime;
  FILE* fout_real;              /* output file                        */
  FILE* fout_imag;              /* output file                        */
  int my_N[3],my_n[3];          /* to init the nfft */
  double t0, t1;
  double t,epsilon=0.0000003;     /* epsilon is a the break criterium for
                                   the iteration */
  unsigned infft_flags = CGNR | PRECOMPUTE_DAMP; /* flags for the infft*/
  double time,min_time,max_time,min_inh,max_inh;
  double real,imag;
  double *w;

  double Ts;
  double W;
  int N3;
  int m=2;
  double sigma = 1.25;

  w = (double*)nfft_malloc(N*N*sizeof(double));

  ftime=fopen("readout_time.dat","r");
  finh=fopen("inh.dat","r");

  min_time=INT_MAX; max_time=INT_MIN;
  for(j=0;j<M;j++)
  {
    fscanf(ftime,"%le ",&time);
    if(time<min_time)
      min_time = time;
    if(time>max_time)
      max_time = time;
  }

  fclose(ftime);

  Ts=(min_time+max_time)/2.0;

  min_inh=INT_MAX; max_inh=INT_MIN;
  for(j=0;j<N*N;j++)
  {
    fscanf(finh,"%le ",&w[j]);
    if(w[j]<min_inh)
      min_inh = w[j];
    if(w[j]>max_inh)
      max_inh = w[j];
  }
  fclose(finh);

  N3=ceil((MAX(fabs(min_inh),fabs(max_inh))*(max_time-min_time)/2.0)*4);


  W=MAX(fabs(min_inh),fabs(max_inh))*2.0;

  fprintf(stderr,"3:  %i %e %e %e %e %e %e\n",N3,W,min_inh,max_inh,min_time,max_time,Ts);

  /* initialise my_plan */
  my_N[0]=N;my_n[0]=ceil(N*sigma);
  my_N[1]=N; my_n[1]=ceil(N*sigma);
  my_N[2]=N3; my_n[2]=ceil(N3*sigma);
  nnfft_init_guru(&my_plan, 3, N*N, M, my_N,my_n,m,
        PRE_PSI| PRE_PHI_HUT| MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F );

  /* precompute lin psi if set */
  if(my_plan.nnfft_flags & PRE_LIN_PSI)
    nnfft_precompute_lin_psi(&my_plan);

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
  ftime=fopen("readout_time.dat","r");

  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fin,"%le %le %le %le ",&my_plan.x[3*j+0],&my_plan.x[3*j+1],&real,&imag);
    my_iplan.y[j]=real+ _Complex_I*imag;
    fscanf(ftime,"%le ",&my_plan.x[3*j+2]);

    my_plan.x[3*j+2] = (my_plan.x[3*j+2]-Ts)*W/N3;
  }

  for(j=0;j<N;j++)
    {
    for(l=0;l<N;l++)
      {
        my_plan.v[3*(N*j+l)+0]=(((double) j) -(((double) N)/2.0))/((double) N);
        my_plan.v[3*(N*j+l)+1]=(((double) l) -(((double) N)/2.0))/((double) N);
        my_plan.v[3*(N*j+l)+2] = w[N*j+l]/W ;
      }
    }

  /* precompute psi */
  if(my_plan.nnfft_flags & PRE_PSI) {
    nnfft_precompute_psi(&my_plan);
    if(my_plan.nnfft_flags & PRE_FULL_PSI)
      nnfft_precompute_full_psi(&my_plan);
  }

  if(my_plan.nnfft_flags & PRE_PHI_HUT)
    nnfft_precompute_phi_hut(&my_plan);

  /* init some guess */
  for(k=0;k<my_plan.N_total;k++)
  {
    my_iplan.f_hat_iter[k]=0.0;
  }

  t0 = nfft_clock_gettime_seconds();

  /* inverse trafo */
  solver_before_loop_complex(&my_iplan);
  for(l=0;l<iteration;l++)
  {
    /* break if dot_r_iter is smaller than epsilon*/
    if(my_iplan.dot_r_iter<epsilon)
    break;
    fprintf(stderr,"%e,  %i of %i\n",sqrt(my_iplan.dot_r_iter),
    l+1,iteration);
    solver_loop_one_step_complex(&my_iplan);
  }

  t1 = nfft_clock_gettime_seconds();
  t = t1-t0;

  fout_real=fopen("output_real.dat","w");
  fout_imag=fopen("output_imag.dat","w");

  for(k=0;k<my_plan.N_total;k++) {

    my_iplan.f_hat_iter[k]*=cexp(2.0*_Complex_I*M_PI*Ts*w[k]);

    fprintf(fout_real,"%le ", creal(my_iplan.f_hat_iter[k]));
    fprintf(fout_imag,"%le ", cimag(my_iplan.f_hat_iter[k]));
  }


  fclose(fout_real);
  fclose(fout_imag);


  /* finalize the infft */
  solver_finalize_complex(&my_iplan);

  /* finalize the nfft */
  nnfft_finalize(&my_plan);

  nfft_free(w);
}

int main(int argc, char **argv)
{
  if (argc <= 5) {
    printf("usage: ./reconstruct_data_inh_nnfft FILENAME N M ITER WEIGHTS\n");
    return 1;
  }

  reconstruct(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));

  return 1;
}
/* \} */
