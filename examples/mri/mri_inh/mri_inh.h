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
