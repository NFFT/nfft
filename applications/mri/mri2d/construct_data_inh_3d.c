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
/* NFFT_PRECISION_DOUBLE, NFFT_PRECISION_SINGLE, or NFFT_PRECISION_LONG_DOUBLE must be defined. */
#include "nfft3mp.h"

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

/**
 * \defgroup applications_mri2d_construct_data_inh_3d construct_data_inh_3d
 * \ingroup applications_mri2d
 * \{
 */

/**
 * construct
 */
static void construct(char * file, int N, int M)
{
  int j;                  /* some variables */
  NFFT_R real;
  NFFT_R w;
  NFFT_R time,min_time,max_time,min_inh,max_inh;
  MRI(inh_3d_plan) my_plan;
  FILE *fp,*fout,*fi,*finh,*ftime;
  int my_N[3],my_n[3];
  int flags = PRE_PHI_HUT| PRE_PSI |MALLOC_X| MALLOC_F_HAT|
                      MALLOC_F| FFTW_INIT| FFTW_MEASURE;

  NFFT_R Ts;
  NFFT_R W;
  int N3;
  int m=2;
  NFFT_R sigma = NFFT_K(1.25);

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
    fscanf(finh,NFFT__FR__ " ",&w);
    if(w<min_inh)
      min_inh = w;
    if(w>max_inh)
      max_inh = w;
  }
  fclose(finh);

  N3=NFFT_CEIL((MAX(NFFT_FABS(min_inh),NFFT_FABS(max_inh))*(max_time-min_time)/NFFT_K(2.0)+m/(2*sigma))*4*sigma);

  W= MAX(NFFT_FABS(min_inh),NFFT_FABS(max_inh))/(NFFT_K(0.5)-((NFFT_R)m)/N3);

  my_N[0]=N; my_n[0]=NFFT_CEIL(N*sigma);
  my_N[1]=N; my_n[1]=NFFT_CEIL(N*sigma);
  my_N[2]=N3; my_n[2]=NFFT_CEIL(N3*sigma);

  /* initialise nfft */
  MRI(inh_3d_init_guru)(&my_plan, my_N, M, my_n, m, sigma, flags,
                      FFTW_MEASURE| FFTW_DESTROY_INPUT);

  ftime=fopen("readout_time.dat","r");
  fp=fopen("knots.dat","r");

  for(j=0;j<my_plan.M_total;j++)
  {
    fscanf(fp,NFFT__FR__ " " NFFT__FR__,&my_plan.plan.x[3*j+0],&my_plan.plan.x[3*j+1]);
    fscanf(ftime,NFFT__FR__ " ",&my_plan.plan.x[3*j+2]);
    my_plan.plan.x[3*j+2] = (my_plan.plan.x[3*j+2]-Ts)*W/N3;
  }
  fclose(fp);
  fclose(ftime);

  finh=fopen("inh.dat","r");
  for(j=0;j<N*N;j++)
  {
    fscanf(finh,NFFT__FR__ " ",&my_plan.w[j]);
    my_plan.w[j]/=W;
  }
  fclose(finh);


  fi=fopen("input_f.dat","r");
  for(j=0;j<N*N;j++)
  {
    fscanf(fi,NFFT__FR__ " ",&real);
    my_plan.f_hat[j] = real*NFFT_CEXP(NFFT_K(2.0)*_Complex_I*NFFT_KPI*Ts*my_plan.w[j]*W);
  }

  NFFT(precompute_one_psi)(&my_plan.plan);

  MRI(inh_3d_trafo)(&my_plan);

  fout=fopen(file,"w");

  for(j=0;j<my_plan.M_total;j++)
  {
    fprintf(fout,NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ " " NFFT__FR__ "\n",my_plan.plan.x[3*j+0],my_plan.plan.x[3*j+1],NFFT_CREAL(my_plan.f[j]),NFFT_CIMAG(my_plan.f[j]));
  }

  fclose(fout);

  MRI(inh_3d_finalize)(&my_plan);
}

int main(int argc, char **argv)
{
  if (argc <= 3) {
    printf("usage: ./construct_data_inh_3d FILENAME N M\n");
    return 1;
  }

  construct(argv[1],atoi(argv[2]),atoi(argv[3]));

  return 0;
}
/* \} */
