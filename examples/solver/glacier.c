/*
 * Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* $Id$ */
#include "config.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3.h"
#include "infft.h"

/**
 * \defgroup examples_solver_glacier Reconstruction of a glacier from \
 scattered data
 * \ingroup examples_solver
 * \{
 */

/** Generalised Sobolev weight */
static double my_weight(double z,double a,double b,double c)
{
    return pow(0.25-z*z,b)/(c+pow(fabs(z),2*a));
}

/** Reconstruction routine */
static void glacier(int N,int M)
{
  int j,k,k0,k1,l,my_N[2],my_n[2];
  double tmp_y;
  nfft_plan p;
  solver_plan_complex ip;
  FILE* fp;

  /* initialise p */
  my_N[0]=N; my_n[0]=X(next_power_of_2)(N);
  my_N[1]=N; my_n[1]=X(next_power_of_2)(N);
  nfft_init_guru(&p, 2, my_N, M, my_n, 6,
		 PRE_PHI_HUT| PRE_FULL_PSI|
		 MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  /* initialise ip, specific */
  solver_init_advanced_complex(&ip,(nfft_mv_plan_complex*)(&p), CGNE| PRECOMPUTE_DAMP);
  fprintf(stderr,"Using the generic solver!");

  /* init nodes */
  fp=fopen("input_data.dat","r");
  for(j=0;j<p.M_total;j++)
  {
      fscanf(fp,"%le %le %le",&p.x[2*j+0],&p.x[2*j+1],&tmp_y);
      ip.y[j]=tmp_y;
  }
  fclose(fp);

  /* precompute psi */
  if(p.flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p);

  /* initialise damping factors */
  if(ip.flags & PRECOMPUTE_DAMP)
    for(k0=0;k0<p.N[0];k0++)
      for(k1=0;k1<p.N[1];k1++)
        ip.w_hat[k0*p.N[1]+k1]=
	    my_weight(((double)(k0-p.N[0]/2))/p.N[0],0.5,3,0.001)*
	    my_weight(((double)(k1-p.N[1]/2))/p.N[1],0.5,3,0.001);

  /* init some guess */
  for(k=0;k<p.N_total;k++)
      ip.f_hat_iter[k]=0;

  /* inverse trafo */
  solver_before_loop_complex(&ip);
  for(l=0;l<40;l++)
    {
      fprintf(stderr,"Residual ||r||=%e,\n",sqrt(ip.dot_r_iter));
      solver_loop_one_step_complex(&ip);
    }

  for(k=0;k<p.N_total;k++)
    printf("%le %le\n",creal(ip.f_hat_iter[k]),cimag(ip.f_hat_iter[k]));

  solver_finalize_complex(&ip);
  nfft_finalize(&p);
}

/** Reconstruction routine with cross validation */
static void glacier_cv(int N,int M,int M_cv,unsigned solver_flags)
{
  int j,k,k0,k1,l,my_N[2],my_n[2];
  double tmp_y,r;
  nfft_plan p,cp;
  solver_plan_complex ip;
  double _Complex* cp_y;
  FILE* fp;
  int M_re=M-M_cv;

  /* initialise p for reconstruction */
  my_N[0]=N; my_n[0]=X(next_power_of_2)(N);
  my_N[1]=N; my_n[1]=X(next_power_of_2)(N);
  nfft_init_guru(&p, 2, my_N, M_re, my_n, 6,
		 PRE_PHI_HUT| PRE_FULL_PSI|
		 MALLOC_X| MALLOC_F_HAT| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);


  /* initialise ip, specific */
  solver_init_advanced_complex(&ip,(nfft_mv_plan_complex*)(&p), solver_flags);

  /* initialise cp for validation */
  cp_y = (double _Complex*) nfft_malloc(M*sizeof(double _Complex));
  nfft_init_guru(&cp, 2, my_N, M, my_n, 6,
		 PRE_PHI_HUT| PRE_FULL_PSI|
		 MALLOC_X| MALLOC_F|
		 FFTW_INIT| FFT_OUT_OF_PLACE,
		 FFTW_MEASURE| FFTW_DESTROY_INPUT);

  cp.f_hat=ip.f_hat_iter;

  /* set up data in cp and cp_y */
  fp=fopen("input_data.dat","r");
  for(j=0;j<cp.M_total;j++)
    {
      fscanf(fp,"%le %le %le",&cp.x[2*j+0],&cp.x[2*j+1],&tmp_y);
      cp_y[j]=tmp_y;
    }
  fclose(fp);

  /* copy part of the data to p and ip */
  for(j=0;j<p.M_total;j++)
  {
      p.x[2*j+0]=cp.x[2*j+0];
      p.x[2*j+1]=cp.x[2*j+1];
      ip.y[j]=tmp_y;
  }

  /* precompute psi */
  if(p.flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&p);

  /* precompute psi */
  if(cp.flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&cp);

  /* initialise damping factors */
  if(ip.flags & PRECOMPUTE_DAMP)
    for(k0=0;k0<p.N[0];k0++)
      for(k1=0;k1<p.N[1];k1++)
        ip.w_hat[k0*p.N[1]+k1]=
	    my_weight(((double)(k0-p.N[0]/2))/p.N[0],0.5,3,0.001)*
	    my_weight(((double)(k1-p.N[1]/2))/p.N[1],0.5,3,0.001);

  /* init some guess */
  for(k=0;k<p.N_total;k++)
      ip.f_hat_iter[k]=0;

  /* inverse trafo */
  solver_before_loop_complex(&ip);
  //  fprintf(stderr,"iteration starts,\t");
  for(l=0;l<40;l++)
    solver_loop_one_step_complex(&ip);

  //fprintf(stderr,"r=%1.2e, ",sqrt(ip.dot_r_iter)/M_re);

  CSWAP(p.f_hat,ip.f_hat_iter);
  nfft_trafo(&p);
  CSWAP(p.f_hat,ip.f_hat_iter);
  nfft_upd_axpy_complex(p.f,-1,ip.y,M_re);
  r=sqrt(X(dot_complex)(p.f,M_re)/X(dot_complex)(cp_y,M));
  fprintf(stderr,"r=%1.2e, ",r);
  printf("$%1.1e$ & ",r);

  nfft_trafo(&cp);
  nfft_upd_axpy_complex(&cp.f[M_re],-1,&cp_y[M_re],M_cv);
  r=sqrt(X(dot_complex)(&cp.f[M_re],M_cv)/X(dot_complex)(cp_y,M));
  fprintf(stderr,"r_1=%1.2e\t",r);
  printf("$%1.1e$ & ",r);

  nfft_finalize(&cp);
  solver_finalize_complex(&ip);
  nfft_finalize(&p);
}


/** Main routine */
int main(int argc, char **argv)
{
  int M_cv;

  if(argc<3)
    {
      fprintf(stderr,"Call this program from the Matlab script glacier.m!");
      exit(-1);
    }

  if(argc==3)
    glacier(atoi(argv[1]),atoi(argv[2]));
  else
    for(M_cv=atoi(argv[3]);M_cv<=atoi(argv[5]);M_cv+=atoi(argv[4]))
      {
	fprintf(stderr,"\nM_cv=%d,\t",M_cv);
	printf("$%d$ & ",M_cv);
	fprintf(stderr,"cgne+damp: ");
	glacier_cv(atoi(argv[1]),atoi(argv[2]),M_cv,CGNE| PRECOMPUTE_DAMP);
	//fprintf(stderr,"cgne: ");
	//glacier_cv(atoi(argv[1]),atoi(argv[2]),M_cv,CGNE);
	fprintf(stderr,"cgnr: ");
	glacier_cv(atoi(argv[1]),atoi(argv[2]),M_cv,CGNR);
	fprintf(stderr,"cgnr: ");
	glacier_cv(atoi(argv[1])/4,atoi(argv[2]),M_cv,CGNR);
	printf("XXX \\\\\n");
      }

  fprintf(stderr,"\n");

  return 1;
}
/* \} */
