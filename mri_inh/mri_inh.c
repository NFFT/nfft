#include "nfft3.h"
#include "options.h"
#include "util.h"
#include "window_defines.h"
#include "math.h"

#define PHI_periodic(x) ((x>0.5)?(PHI(x-1.0,0)):((x<-0.5)?PHI(x+1.0,0):PHI(x,0)))

void mri_inh_2d1d_trafo(mri_inh_2d1d_plan *that) {
  int l,j;
  complex *f = (complex*) fftw_malloc(that->M_total*sizeof(complex));
  complex *f_hat = (complex*) fftw_malloc(that->N_total*sizeof(complex));

  nfft_plan *ths = (nfft_plan*) malloc(sizeof(nfft_plan));
  nfft_init_guru(ths, 1, &that->N3, 1,&that->N3, 6,0,0);
  memset(f,0,that->M_total*sizeof(complex));
  for(j=0;j<that->N_total;j++)
  {
    that->f_hat[j]/=PHI_HUT(ths->n[0]*that->w[j],0);
    f_hat[j]=that->f_hat[j];
  }

  for(l=-ths->n[0]/2;l<ths->n[0]/2;l++) {
    for(j=0;j<that->N_total;j++)
      that->f_hat[j]*=cexp(-2*PI*I*that->w[j]*l);
    nfft_trafo((nfft_plan*)that);
    for(j=0;j<that->M_total;j++)
      f[j]+=that->f[j]*PHI_periodic(that->t[j]-((double)l)/((double)ths->n[0]));
    for(j=0;j<that->N_total;j++)
      that->f_hat[j]=f_hat[j];
  }

  fftw_free(that->f);
  that->f=f;
  fftw_free(f_hat);
  
  nfft_finalize(ths);
  free(ths);
}

void mri_inh_2d1d_adjoint(mri_inh_2d1d_plan *that) {
  int l,j;
  complex *f = (complex*) fftw_malloc(that->M_total*sizeof(complex));
  complex *f_hat = (complex*) fftw_malloc(that->N_total*sizeof(complex));

  nfft_plan *ths = (nfft_plan*) malloc(sizeof(nfft_plan));
  nfft_init_guru(ths, 1, &that->N3, 1,&that->N3, 6,0,0);
  memset(f_hat,0,that->N_total*sizeof(complex));
  for(j=0;j<that->M_total;j++)
  {
    f[j]=that->f[j];
  }


  
  for(l=-ths->n[0]/2;l<ths->n[0]/2;l++) {
    
    for(j=0;j<that->M_total;j++)
      that->f[j]*=PHI_periodic(that->t[j]-((double)l)/((double)ths->n[0]));
    nfft_adjoint((nfft_plan*)that);
    for(j=0;j<that->N_total;j++)
      f_hat[j]+=that->f_hat[j]*cexp(2*PI*I*that->w[j]*l);
    for(j=0;j<that->M_total;j++)
      that->f[j]=f[j];
  }

  for(j=0;j<that->N_total;j++)
  {
    f_hat[j] /= PHI_HUT(ths->n[0]*that->w[j],0);
  }



  fftw_free(that->f_hat);
  that->f_hat=f_hat;
  fftw_free(f);
  
  nfft_finalize(ths);
  free(ths);
}

void mri_inh_2d1d_init_guru(mri_inh_2d1d_plan *ths, int *N, int M, int *n,
                    int m, unsigned nfft_flags, unsigned fftw_flags) {

  nfft_init_guru((nfft_plan*)ths,2,N,M,n,m,nfft_flags,fftw_flags);
  ths->N3=N[2];
  
  ths->t = (double*) fftw_malloc(ths->M_total*sizeof(double));
  ths->w = (double*) fftw_malloc(ths->N_total*sizeof(double));
}

void mri_inh_2d1d_finalize(mri_inh_2d1d_plan *ths) {
  fftw_free(ths->t);
  fftw_free(ths->w);
  nfft_finalize((nfft_plan*)ths);
}























void mri_inh_3d_trafo(mri_inh_3d_plan *that) {
  int l,j;
  nfft_plan *ths = (nfft_plan*) malloc(sizeof(nfft_plan));
  nfft_init_guru(ths, 1, &that->N3, 1,&that->N3, 6,0,0);

  that->plan.f =that->f ;
  


  for(j=0;j<that->N_total;j++) {
    for(l=-ths->n[0]/2;l<ths->n[0]/2;l++)
    {
      that->plan.f_hat[j*ths->n[0]+(l+ths->n[0]/2)]= that->f_hat[j]*PHI_periodic(that->w[j]-((double)l)/((double)ths->n[0]));
    }
  }
  
  nfft_trafo(&that->plan);
  
  for(j=0;j<that->M_total;j++)
  {
    that->f[j] /= PHI_HUT(ths->n[0]*that->plan.x[3*j+2],0);
  }

  nfft_finalize(ths);
  free(ths);
}

void mri_inh_3d_adjoint(mri_inh_3d_plan *that) {
  int l,j;
  nfft_plan *ths = (nfft_plan*) malloc(sizeof(nfft_plan));
  nfft_init_guru(ths, 1, &that->N3, 1,&that->N3, 6,0,0);

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
      that->f_hat[j]+= that->plan.f_hat[j*ths->n[0]+(l+ths->n[0]/2)]*PHI_periodic(that->w[j]-((double)l)/((double)ths->n[0]));
    }
  }

  
  nfft_finalize(ths);
  free(ths);
}

void mri_inh_3d_init_guru(mri_inh_3d_plan *ths, int *N, int M, int *n,
                    int m, unsigned nfft_flags, unsigned fftw_flags) {
  ths->N3=N[2];
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
