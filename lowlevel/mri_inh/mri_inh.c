#include "string.h"

#include "nfft3.h"
#include "options.h"
#include "util.h"
#include "window_defines.h"
#include "math.h"

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
  double *spline_coeffs;                /**< input for de Boor algorithm, if   
                                             B_SPLINE or SINC_2m is defined   */ 
} window_funct_plan;


/**
 * init the window_funct_plan
 */
void window_funct_init(window_funct_plan* ths, int m, int n, double sigma) {
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
  complex *f = (complex*) fftw_malloc(that->M_total*sizeof(complex));
  complex *f_hat = (complex*) fftw_malloc(that->N_total*sizeof(complex));

  window_funct_plan *ths = (window_funct_plan*) fftw_malloc(sizeof(window_funct_plan));
	window_funct_init(ths,that->plan.m,that->N3,that->sigma3);
  
	/* the pointers that->f and that->f_hat have been modified by the solver */
	that->plan.f = that->f;
  that->plan.f_hat = that->f_hat;

	
	memset(f,0,that->M_total*sizeof(complex));
  for(j=0;j<that->N_total;j++)
  {
    f_hat[j]=that->f_hat[j];
  }

  for(l=-ths->n[0]/2;l<=ths->n[0]/2;l++) {
    for(j=0;j<that->N_total;j++)
      that->f_hat[j]*=cexp(-2*PI*I*that->w[j]*((double)l))/PHI_HUT(ths->n[0]*that->w[j],0);
    nfft_trafo(&that->plan);
    for(j=0;j<that->M_total;j++){
      /* PHI has compact support */
			if(fabs(that->t[j]-((double)l)/((double)ths->n[0]))<that->plan.m/((double)ths->n[0])) 
        f[j]+=that->f[j]*PHI(that->t[j]-((double)l)/((double)ths->n[0]),0);
    }
    for(j=0;j<that->N_total;j++)
      that->f_hat[j]=f_hat[j];
  }

  fftw_free(that->plan.f);
  that->f=f;
	that->plan.f = that->f;
	
	fftw_free(f_hat);
  
	WINDOW_HELP_FINALIZE
  free(ths);
}

void mri_inh_2d1d_adjoint(mri_inh_2d1d_plan *that) {
  int l,j;
  complex *f = (complex*) fftw_malloc(that->M_total*sizeof(complex));
  complex *f_hat = (complex*) fftw_malloc(that->N_total*sizeof(complex));

  window_funct_plan *ths = (window_funct_plan*) fftw_malloc(sizeof(window_funct_plan));
	window_funct_init(ths,that->plan.m,that->N3,that->sigma3);
  
	memset(f_hat,0,that->N_total*sizeof(complex));
	
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
        that->f[j]*=PHI(that->t[j]-((double)l)/((double)ths->n[0]),0);
      else
      	that->f[j]=0.0;
    }
    nfft_adjoint(&that->plan);
    for(j=0;j<that->N_total;j++)
      f_hat[j]+=that->f_hat[j]*cexp(2*PI*I*that->w[j]*((double)l));
    for(j=0;j<that->M_total;j++)
      that->f[j]=f[j];
  }

  for(j=0;j<that->N_total;j++)
  {
    f_hat[j] /= PHI_HUT(ths->n[0]*that->w[j],0);
  }

  fftw_free(that->plan.f_hat);
  that->f_hat=f_hat;
  that->plan.f_hat = that->f_hat;
  
	fftw_free(f);
  
	WINDOW_HELP_FINALIZE
  free(ths);
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
  
  ths->t = (double*) fftw_malloc(ths->M_total*sizeof(double));
  ths->w = (double*) fftw_malloc(ths->N_total*sizeof(double));
}

void mri_inh_2d1d_finalize(mri_inh_2d1d_plan *ths) {
  fftw_free(ths->t);
  fftw_free(ths->w);
	
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
  window_funct_plan *ths = (window_funct_plan*) fftw_malloc(sizeof(window_funct_plan));
	window_funct_init(ths,that->plan.m,that->N3,that->sigma3);

	/* the pointers that->f has been modified by the solver */
  that->plan.f =that->f ;
  


  for(j=0;j<that->N_total;j++) {
    for(l=-ths->n[0]/2;l<ths->n[0]/2;l++)
    {
      /* PHI has compact support */
      if(fabs(that->w[j]-((double)l)/((double)ths->n[0]))<ths->m/((double)ths->n[0])) 
        that->plan.f_hat[j*ths->n[0]+(l+ths->n[0]/2)]= that->f_hat[j]*PHI(that->w[j]-((double)l)/((double)ths->n[0]),0);
      else
	      that->plan.f_hat[j*ths->n[0]+(l+ths->n[0]/2)]=0.0;
    }
  }
  
  nfft_trafo(&that->plan);
  
  for(j=0;j<that->M_total;j++)
  {
    that->f[j] /= PHI_HUT(ths->n[0]*that->plan.x[3*j+2],0);
  }

	WINDOW_HELP_FINALIZE
  free(ths);
}

void mri_inh_3d_adjoint(mri_inh_3d_plan *that) {
  int l,j;
  window_funct_plan *ths = (window_funct_plan*) fftw_malloc(sizeof(window_funct_plan));
	window_funct_init(ths,that->plan.m,that->N3,that->sigma3);

	/* the pointers that->f has been modified by the solver */
  that->plan.f =that->f ;
  
  for(j=0;j<that->M_total;j++)
  {
    that->f[j] /= PHI_HUT(ths->n[0]*that->plan.x[3*j+2],0);
  }
  
  nfft_adjoint(&that->plan);

  for(j=0;j<that->N_total;j++) {
    that->f_hat[j]=0.0;
    for(l=-ths->n[0]/2;l<ths->n[0]/2;l++)
    {
      /* PHI has compact support */
      if(fabs(that->w[j]-((double)l)/((double)ths->n[0]))<ths->m/((double)ths->n[0])) 
        that->f_hat[j]+= that->plan.f_hat[j*ths->n[0]+(l+ths->n[0]/2)]*PHI(that->w[j]-((double)l)/((double)ths->n[0]),0);
    }
  }

  
	WINDOW_HELP_FINALIZE
  free(ths);
}

void mri_inh_3d_init_guru(mri_inh_3d_plan *ths, int *N, int M, int *n,
                    int m, double sigma, unsigned nfft_flags, unsigned fftw_flags) {
  ths->N3=N[2];
	ths->sigma3=sigma;
  nfft_init_guru(&ths->plan,3,N,M,n,m,nfft_flags,fftw_flags);
  ths->N_total = N[0]*N[1];
  ths->M_total = ths->plan.M_total;
  ths->f = ths->plan.f;
  ths->f_hat = (complex*) fftw_malloc(ths->N_total*sizeof(complex));
  ths->w = (double*) fftw_malloc(ths->N_total*sizeof(double));
}

void mri_inh_3d_finalize(mri_inh_3d_plan *ths) {
  fftw_free(ths->w);
  fftw_free(ths->f_hat);
  nfft_finalize(&ths->plan);
}
