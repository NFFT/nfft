#include "nfft3.h"

typedef struct mri_inh_plan_
{
  /** api */
  MACRO_MV_PLAN(complex);          

  int d;                                /**< dimension, rank                 */
  int *N;                               /**< cut-off-frequencies             */
  double *sigma;                        /**< oversampling-factor             */
  int *n;                               /**< fftw-length = sigma*N           */
  int n_total;                          /**< total size of fftw              */
  int m;                                /**< cut-off parameter in time-domain*/
  double *b;                            /**< shape parameters                */
  int K;                                /**< number of precomp. uniform psi  */

  unsigned nfft_flags;                  /**< flags for precomputation, malloc*/
  unsigned fftw_flags;                  /**< flags for the fftw              */

  double *x;                            /**< nodes (in time/spatial domain)  */

  /** internal*/
  fftw_plan  my_fftw_plan1;             /**< fftw_plan forward               */
  fftw_plan  my_fftw_plan2;             /**< fftw_plan backward              */

  double **c_phi_inv;                   /**< precomputed data, matrix D      */
  double *psi;                          /**< precomputed data, matrix B      */
  int *psi_index_g;                     /**< only for PRE_FULL_PSI           */
  int *psi_index_f;                     /**< only for PRE_FULL_PSI           */

  complex *g;
  complex *g_hat;
  complex *g1;                          /**< input of fftw                   */
  complex *g2;                          /**< output of fftw                  */

  double *spline_coeffs;                /**< input for de Boor algorithm, if
                                             B_SPLINE or SINC_2m is defined  */
  int N3;
  double *t;
  double *w;
} mri_inh_plan;


void mri_inh_trafo(mri_inh_plan *ths);

void mri_inh_adjoint(mri_inh_plan *ths);

void mri_inh_init_guru(mri_inh_plan *ths, int d, int *N, int M, int *n,
                    int m, unsigned nfft_flags, unsigned fftw_flags);

void mri_inh_finalize(mri_inh_plan *ths);

void mri_inh_trafo2(mri_inh_plan *ths);

void mri_inh_adjoint2(mri_inh_plan *ths);

void mri_inh_init_guru2(mri_inh_plan *ths, int d, int *N, int M, int *n,
                    int m, unsigned nfft_flags, unsigned fftw_flags);

void mri_inh_finalize2(mri_inh_plan *ths);

MACRO_SOLVER_PLAN(mri_inh, complex)




void mri_inh_trafo2(mri_inh_plan *that) {

}

void mri_inh_adjoint2(mri_inh_plan *that) {
  int l,j;
  complex *f_hat = (complex*) fftw_malloc(that->N_total*sizeof(complex));

  nfft_plan *ths = (nfft_plan*) malloc(sizeof(nfft_plan));
  nfft_init_1d(ths,that->N3,1);
  memset(f_hat,0,that->N_total*sizeof(complex));

  for(j=0;j<that->M_total;j++)
  {
    that->f[j] /= PHI_HUT(ths->n[0]*that->x[3*j+2],0);
  }
  
  nfft_adjoint((nfft_plan*)that);
  

  for(j=0;j<that->N[0]*that->N[1];j++) {
    for(l=-ths->n[0]/2;l<ths->n[0]/2;l++)
    {
      f_hat[j]+= that->f_hat[j*ths->n[0]+(l+ths->n[0]/2)]*PHI_periodic(that->w[j]-((double)l)/((double)ths->n[0]));
    }
  }





  fftw_free(that->f_hat);
  that->f_hat=f_hat;
  
  nfft_finalize(ths);
  free(ths);
}

void mri_inh_init_guru2(mri_inh_plan *ths, int d, int *N, int M, int *n,
                    int m, unsigned nfft_flags, unsigned fftw_flags) {
  ths->N3=N[2];
  N[2]=n[2];
  nfft_init_guru((nfft_plan*)ths,d,N,M,n,m,nfft_flags,fftw_flags);
  ths->w = (double*) fftw_malloc(ths->N_total*sizeof(double));
}

void mri_inh_finalize2(mri_inh_plan *ths) {
  fftw_free(ths->w);
  nfft_finalize((nfft_plan*)ths);
}
