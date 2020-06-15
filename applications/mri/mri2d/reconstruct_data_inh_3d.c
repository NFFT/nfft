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
 * \defgroup applications_mri2d_reconstruct_data_inh_3d reconstruct_data_inh_3d
 * \ingroup applications_mri2d
 * \{
 */

static void reconstruct(char* filename,int N,int M,int iteration , int weight)
{
  int j,k,l;
  double t0, t1;
  double time,min_time,max_time,min_inh,max_inh;
  double t,real,imag;
  double w,epsilon=0.0000003;     /* epsilon is a the break criterium for
                                   the iteration */;
  mri_inh_3d_plan my_plan;
  solver_plan_complex my_iplan;
  FILE* fp,*fw,*fout_real,*fout_imag,*finh,*ftime;
  int my_N[3],my_n[3];
  int flags = PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
                      MALLOC_F| FFTW_INIT;
  unsigned infft_flags = CGNR | PRECOMPUTE_DAMP;

  double Ts;
  double W;
  int N3;
  int m=2;
  double sigma = 1.25;

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
    fscanf(finh,"%le ",&w);
    if(w<min_inh)
      min_inh = w;
    if(w>max_inh)
      max_inh = w;
  }
  fclose(finh);

  N3=ceil((MAX(fabs(min_inh),fabs(max_inh))*(max_time-min_time)/2.0+m/(2*sigma))*4*sigma);
	/* N3 has to be even */
	if(N3%2!=0)
	  N3++;

  W= MAX(fabs(min_inh),fabs(max_inh))/(0.5-((double) m)/N3);

  my_N[0]=N;my_n[0]=ceil(N*sigma);
  my_N[1]=N; my_n[1]=ceil(N*sigma);
  my_N[2]=N3; my_n[2]=ceil(N3*sigma);

  /* initialise nfft */
  mri_inh_3d_init_guru(&my_plan, my_N, M, my_n, m, sigma, flags,
                      FFTW_MEASURE| FFTW_DESTROY_INPUT);

  if (weight)
    infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

  /* initialise my_iplan, advanced */
  solver_init_advanced_complex(&my_iplan,(nfft_mv_plan_complex*)(&my_plan), infft_flags );

  /* get the weights */
  if(my_iplan.flags & PRECOMPUTE_WEIGHT)
  {
    fw=fopen("weights.dat","r");
    for(j=0;j<my_plan.M_total;j++)
    {
        fscanf(fw,"%le ",&my_iplan.w[j]);
    }
    fclose(fw);
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

  fp=fopen(filename,"r");
  ftime=fopen("readout_time.dat","r");

  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,"%le %le %le %le",&my_plan.plan.x[3*j+0],&my_plan.plan.x[3*j+1],&real,&imag);
    my_iplan.y[j]=real+ _Complex_I*imag;
    fscanf(ftime,"%le ",&my_plan.plan.x[3*j+2]);

    my_plan.plan.x[3*j+2] = (my_plan.plan.x[3*j+2]-Ts)*W/N3;
  }
  fclose(fp);
  fclose(ftime);


  finh=fopen("inh.dat","r");
  for(j=0;j<N*N;j++)
  {
    fscanf(finh,"%le ",&my_plan.w[j]);
    my_plan.w[j]/=W;
  }
  fclose(finh);


  if(my_plan.plan.flags & PRE_PSI) {
    nfft_precompute_psi(&my_plan.plan);
  }
  if(my_plan.plan.flags & PRE_FULL_PSI) {
      nfft_precompute_full_psi(&my_plan.plan);
  }

  /* init some guess */
  for(j=0;j<my_plan.N_total;j++)
  {
    my_iplan.f_hat_iter[j]=0.0;
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

  for (j=0;j<N*N;j++) {
    /* Verschiebung wieder herausrechnen */
    my_iplan.f_hat_iter[j]*=cexp(-2.0*_Complex_I*M_PI*Ts*my_plan.w[j]*W);

    fprintf(fout_real,"%le ",creal(my_iplan.f_hat_iter[j]));
    fprintf(fout_imag,"%le ",cimag(my_iplan.f_hat_iter[j]));
  }

  fclose(fout_real);
  fclose(fout_imag);
  solver_finalize_complex(&my_iplan);
  mri_inh_3d_finalize(&my_plan);
}


int main(int argc, char **argv)
{
  if (argc <= 5) {

    printf("usage: ./reconstruct_data_inh_3d FILENAME N M ITER WEIGHTS\n");
    return 1;
  }

  reconstruct(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));

  return 1;
}
/* \} */
