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

#include <string.h>
#include <math.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include "nfft3.h"
#include "infft.h"

/**
 * window_funct_plan is a plan to use the window functions
 * independent of the nfft
 */

typedef struct window_funct_plan_ {
  int d;
	int m;
	int n[1];
	double sigma[1];
	double *b;
} window_funct_plan;

/**
 * init the window_funct_plan
 */
static void window_funct_init(window_funct_plan* ths, int m, int n, double sigma) {
	ths->d=1;
	ths->m=m;
	ths->n[0]=n;
	ths->sigma[0]=sigma;
  WINDOW_HELP_INIT
}

/*
 * mri_inh_2d1d
 */

void mri_inh_2d1d_trafo(mri_inh_2d1d_plan *that) {
  int l,j;
  double _Complex *f = (double _Complex*) nfft_malloc(that->M_total*sizeof(double _Complex));
  double _Complex *f_hat = (double _Complex*) nfft_malloc(that->N_total*sizeof(double _Complex));

  window_funct_plan *ths = (window_funct_plan*) nfft_malloc(sizeof(window_funct_plan));
	window_funct_init(ths,that->plan.m,that->N3,that->sigma3);

	/* the pointers that->f and that->f_hat have been modified by the solver */
	that->plan.f = that->f;
  that->plan.f_hat = that->f_hat;


	memset(f,0,that->M_total*sizeof(double _Complex));
  for(j=0;j<that->N_total;j++)
  {
    f_hat[j]=that->f_hat[j];
  }

  for(l=-ths->n[0]/2;l<=ths->n[0]/2;l++) {
    for(j=0;j<that->N_total;j++)
      that->f_hat[j]*=cexp(-2*KPI*_Complex_I*that->w[j]*((double)l))/PHI_HUT(ths->n[0], ths->n[0]*that->w[j],0);
    nfft_trafo(&that->plan);
    for(j=0;j<that->M_total;j++){
      /* PHI has compact support */
			if(fabs(that->t[j]-((double)l)/((double)ths->n[0]))<that->plan.m/((double)ths->n[0]))
      {
        double phi_val = PHI(ths->n[0],that->t[j]-((double)l)/((double)ths->n[0]),0);
        f[j]+=that->f[j]*phi_val;
// the line below causes internal compiler error for gcc 4.7.1
//        f[j]+=that->f[j]*PHI(ths->n[0],that->t[j]-((double)l)/((double)ths->n[0]),0);
      }
    }
    for(j=0;j<that->N_total;j++)
      that->f_hat[j]=f_hat[j];
  }

  nfft_free(that->plan.f);
  that->f=f;
  that->plan.f = that->f;

  nfft_free(f_hat);

  WINDOW_HELP_FINALIZE
  nfft_free(ths);
}

void mri_inh_2d1d_adjoint(mri_inh_2d1d_plan *that) {
  int l,j;
  double _Complex *f = (double _Complex*) nfft_malloc(that->M_total*sizeof(double _Complex));
  double _Complex *f_hat = (double _Complex*) nfft_malloc(that->N_total*sizeof(double _Complex));

  window_funct_plan *ths = (window_funct_plan*) nfft_malloc(sizeof(window_funct_plan));
	window_funct_init(ths,that->plan.m,that->N3,that->sigma3);

	memset(f_hat,0,that->N_total*sizeof(double _Complex));

	/* the pointers that->f and that->f_hat have been modified by the solver */
	that->plan.f = that->f;
  that->plan.f_hat = that->f_hat;

	for(j=0;j<that->M_total;j++)
  {
    f[j]=that->f[j];
  }



  for(l=-ths->n[0]/2;l<=ths->n[0]/2;l++) {

    for(j=0;j<that->M_total;j++) {
      /* PHI has compact support */
      if(fabs(that->t[j]-((double)l)/((double)ths->n[0]))<that->plan.m/((double)ths->n[0]))
        that->f[j]*=PHI(ths->n[0],that->t[j]-((double)l)/((double)ths->n[0]),0);
      else
      	that->f[j]=0.0;
    }
    nfft_adjoint(&that->plan);
    for(j=0;j<that->N_total;j++)
      f_hat[j]+=that->f_hat[j]*cexp(2*KPI*_Complex_I*that->w[j]*((double)l));
    for(j=0;j<that->M_total;j++)
      that->f[j]=f[j];
  }

  for(j=0;j<that->N_total;j++)
  {
    f_hat[j] /= PHI_HUT(ths->n[0],ths->n[0]*that->w[j],0);
  }

  nfft_free(that->plan.f_hat);
  that->f_hat=f_hat;
  that->plan.f_hat = that->f_hat;

	nfft_free(f);

	WINDOW_HELP_FINALIZE
  nfft_free(ths);
}

void mri_inh_2d1d_init_guru(mri_inh_2d1d_plan *ths, int *N, int M, int *n,
                    int m, double sigma, unsigned nfft_flags, unsigned fftw_flags) {

  nfft_init_guru(&ths->plan,2,N,M,n,m,nfft_flags,fftw_flags);
  ths->N3=N[2];
	ths->sigma3=sigma;
  ths->N_total = ths->plan.N_total;
  ths->M_total = ths->plan.M_total;
  ths->f = ths->plan.f;
  ths->f_hat = ths->plan.f_hat;

  ths->t = (double*) nfft_malloc(ths->M_total*sizeof(double));
  ths->w = (double*) nfft_malloc(ths->N_total*sizeof(double));

  ths->mv_trafo = (void (*) (void* ))mri_inh_2d1d_trafo;
  ths->mv_adjoint = (void (*) (void* ))mri_inh_2d1d_adjoint;
}

void mri_inh_2d1d_finalize(mri_inh_2d1d_plan *ths) {
  nfft_free(ths->t);
  nfft_free(ths->w);

	/* the pointers ths->f and ths->f_hat have been modified by the solver */
	ths->plan.f = ths->f;
  ths->plan.f_hat = ths->f_hat;

	nfft_finalize(&ths->plan);
}

/*
 * mri_inh_3d
 */

void mri_inh_3d_trafo(mri_inh_3d_plan *that) {
  int l,j;
  window_funct_plan *ths = (window_funct_plan*) nfft_malloc(sizeof(window_funct_plan));
	window_funct_init(ths,that->plan.m,that->N3,that->sigma3);

	/* the pointers that->f has been modified by the solver */
  that->plan.f =that->f ;



  for(j=0;j<that->N_total;j++) {
    for(l=-ths->n[0]/2;l<ths->n[0]/2;l++)
    {
      /* PHI has compact support */
      if(fabs(that->w[j]-((double)l)/((double)ths->n[0]))<ths->m/((double)ths->n[0]))
        that->plan.f_hat[j*ths->n[0]+(l+ths->n[0]/2)]= that->f_hat[j]*PHI(ths->n[0],that->w[j]-((double)l)/((double)ths->n[0]),0);
      else
	      that->plan.f_hat[j*ths->n[0]+(l+ths->n[0]/2)]=0.0;
    }
  }

  nfft_trafo(&that->plan);

  for(j=0;j<that->M_total;j++)
  {
    that->f[j] /= PHI_HUT(ths->n[0],ths->n[0]*that->plan.x[3*j+2],0);
  }

	WINDOW_HELP_FINALIZE
  nfft_free(ths);
}

void mri_inh_3d_adjoint(mri_inh_3d_plan *that) {
  int l,j;
  window_funct_plan *ths = (window_funct_plan*) nfft_malloc(sizeof(window_funct_plan));
	window_funct_init(ths,that->plan.m,that->N3,that->sigma3);

	/* the pointers that->f has been modified by the solver */
  that->plan.f =that->f ;

  for(j=0;j<that->M_total;j++)
  {
    that->f[j] /= PHI_HUT(ths->n[0],ths->n[0]*that->plan.x[3*j+2],0);
  }

  nfft_adjoint(&that->plan);

  for(j=0;j<that->N_total;j++) {
    that->f_hat[j]=0.0;
    for(l=-ths->n[0]/2;l<ths->n[0]/2;l++)
    {
      /* PHI has compact support */
      if(fabs(that->w[j]-((double)l)/((double)ths->n[0]))<ths->m/((double)ths->n[0]))
        that->f_hat[j]+= that->plan.f_hat[j*ths->n[0]+(l+ths->n[0]/2)]*PHI(ths->n[0],that->w[j]-((double)l)/((double)ths->n[0]),0);
    }
  }


	WINDOW_HELP_FINALIZE
  nfft_free(ths);
}

void mri_inh_3d_init_guru(mri_inh_3d_plan *ths, int *N, int M, int *n,
                    int m, double sigma, unsigned nfft_flags, unsigned fftw_flags) {
  ths->N3=N[2];
	ths->sigma3=sigma;
  nfft_init_guru(&ths->plan,3,N,M,n,m,nfft_flags,fftw_flags);
  ths->N_total = N[0]*N[1];
  ths->M_total = ths->plan.M_total;
  ths->f = ths->plan.f;
  ths->f_hat = (double _Complex*) nfft_malloc(ths->N_total*sizeof(double _Complex));
  ths->w = (double*) nfft_malloc(ths->N_total*sizeof(double));

  ths->mv_trafo = (void (*) (void* ))mri_inh_3d_trafo;
  ths->mv_adjoint = (void (*) (void* ))mri_inh_3d_adjoint;
}

void mri_inh_3d_finalize(mri_inh_3d_plan *ths) {
  nfft_free(ths->w);
  nfft_free(ths->f_hat);
  nfft_finalize(&ths->plan);
}
