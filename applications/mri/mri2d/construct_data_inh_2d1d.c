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
#include <limits.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3.h"

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

/**
 * \defgroup applications_mri2d_construct_data_inh_2d1d construct_data__inh_2d1d
 * \ingroup applications_mri2d
 * \{
 */

/**
 * construct
 */
static void construct(char * file, int N, int M)
{
  int j;                  /* some variables */
  double real;
  double w;
  double time,min_time,max_time,min_inh,max_inh;
  mri_inh_2d1d_plan my_plan;
  FILE *fp,*fout,*fi,*finh,*ftime;
  int my_N[3],my_n[3];
  int flags = PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
                      MALLOC_F| FFTW_INIT| FFTW_MEASURE;

  double Ts;
  double W,T;
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
  T=((max_time-min_time)/2.0)/(0.5-((double) m)/N3);
  W=N3/T;

  my_N[0]=N; my_n[0]=ceil(N*sigma);
  my_N[1]=N; my_n[1]=ceil(N*sigma);
  my_N[2]=N3; my_n[2]=N3;

  /* initialise nfft */
  mri_inh_2d1d_init_guru(&my_plan, my_N, M, my_n, m, sigma, flags,
                      FFTW_MEASURE| FFTW_DESTROY_INPUT);

  ftime=fopen("readout_time.dat","r");
  fp=fopen("knots.dat","r");

  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,"%le %le ",&my_plan.plan.x[2*j+0],&my_plan.plan.x[2*j+1]);
    fscanf(ftime,"%le ",&my_plan.t[j]);
    my_plan.t[j] = (my_plan.t[j]-Ts)/T;
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


  fi=fopen("input_f.dat","r");
  for(j=0;j<N*N;j++)
  {
    fscanf(fi,"%le ",&real);
    my_plan.f_hat[j] = real*cexp(2.0*_Complex_I*M_PI*Ts*my_plan.w[j]*W);
  }

  if(my_plan.plan.flags & PRE_PSI)
    nfft_precompute_psi(&my_plan.plan);

  mri_inh_2d1d_trafo(&my_plan);

  fout=fopen(file,"w");

  for(j=0;j<my_plan.M_total;j++)
  {
    fprintf(fout,"%le %le %le %le\n",my_plan.plan.x[2*j+0],my_plan.plan.x[2*j+1],creal(my_plan.f[j]),cimag(my_plan.f[j]));
  }

  fclose(fout);

  mri_inh_2d1d_finalize(&my_plan);
}

int main(int argc, char **argv)
{
  if (argc <= 3) {
    printf("usage: ./construct_data_inh_2d1d FILENAME N M\n");
    return 1;
  }

  construct(argv[1],atoi(argv[2]),atoi(argv[3]));

  return 1;
}
/* \} */
