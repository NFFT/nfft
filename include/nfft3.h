/*! \file nfft3.h
 *  \brief Header file for the nfft3 library.
 */
#ifndef NFFT3_H
#define NFFT3_H

#include "fftw3.h"

/** 
 * Constant symbols for precomputation and memory usage
 */
#define PRE_PHI_HUT      (1U<< 0)
#define PRE_PSI          (1U<< 1)
#define PRE_FULL_PSI     (1U<< 2)
#define MALLOC_X         (1U<< 3)
#define MALLOC_F_HAT     (1U<< 4)
#define MALLOC_F         (1U<< 5)
#define FFT_OUT_OF_PLACE (1U<< 6)

#define MACRO_MV_PLAN(float_type)			                      \
  int N_total;                          /**< total number of Fourier coeffs.*/\
  int M_total;                          /**< total number of samples        */\
  float_type *f_hat;                    /**< Fourier coefficients           */\
  float_type *f;                        /**< samples                        */\
  
/** @defgroup nfft Group
 * Direct and fast computation of the
 * discrete Fourier transform at nonequispaced knots
 * @{ 
 */
typedef struct nfft_plan_ 
{
  /** api */
  MACRO_MV_PLAN(complex);          

  int d;                                /**< dimension, rank                 */
  int *N;                               /**< cut-off-frequencies             */
  double *sigma;	                /**< oversampling-factor             */
  int *n;                               /**< fftw-length = sigma*N           */
  int n_total;                          /**< total size of fftw              */	
  int m;                                /**< cut-off parameter in time-domain*/
  double *b;                            /**< shape parameters                */

  unsigned nfft_flags;                  /**< flags for precomputation, malloc*/
  unsigned fftw_flags;                  /**< flags for the fftw              */

  double *x;                            /**< nodes (in time/spatial domain)  */

  /** internal*/
  fftw_plan  my_fftw_plan1;             /**< fftw_plan forward               */
  fftw_plan  my_fftw_plan2;		/**< fftw_plan backward              */

  double **c_phi_inv;                   /**< precomputed data, matrix D      */
  double *psi;				/**< precomputed data, matrix B      */
  int *psi_index_g;                     /**< only for PRE_FULL_PSI           */
  int *psi_index_f;                     /**< only for PRE_FULL_PSI           */

  complex *g;
  complex *g_hat;
  complex *g1;                          /**< input of fftw                   */
  complex *g2;                          /**< output of fftw                  */

  double *spline_coeffs;          	/**< input for de Boor algorithm, if  
					     B_SPLINE or SINC_2m is defined  */
} nfft_plan;


/** see formula (1.1), computes
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 */
void ndft_trafo(nfft_plan *ths);

/** see formula (1.1), computes
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 */
void ndft_conjugated(nfft_plan *ths);

/** see formula (1.2), computes
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 */
void ndft_adjoint(nfft_plan *ths);

/** see formula (1.2), computes
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(-2(pi) k x[j])
 */
void ndft_transposed(nfft_plan *ths);

/** see formula (1.1), computes in a fast and approximative way
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 */
void nfft_trafo(nfft_plan *ths);

/** see formula (1.1), computes in a fast and approximative way
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 */
void nfft_conjugated(nfft_plan *ths);

/** see formula (1.2), computes in a fast and approximative way
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 */
void nfft_adjoint(nfft_plan *ths);

/** see formula (1.2), computes in a fast and approximative way
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 */
void nfft_transposed(nfft_plan *ths);

/** wrapper for nfft_init and d=1
 */
void nfft_init_1d(nfft_plan *ths, int N1, int M);

/** wrapper for nfft_init and d=2
 */
void nfft_init_2d(nfft_plan *ths, int N1, int N2, int M);

/** wrapper for nfft_init and d=3
 */
void nfft_init_3d(nfft_plan *ths, int N1, int N2, int N3, int M);

/** initialisation of transform, simple interface
 */
void nfft_init(nfft_plan *ths, int d, int *N, int M);

/** initialisation of transform, advanced interface
 */
void nfft_init_advanced(nfft_plan *ths, int d, int *N, int M,
			unsigned nfft_flags_on, unsigned nfft_flags_off);

/** initialisation of transform, guru interface
 */
void nfft_init_guru(nfft_plan *ths, int d, int *N, int M, int *n,
	            int m, unsigned nfft_flags, unsigned fftw_flags);

/** precomputes equally spaced values of psi
 *  if PRE_LIN_PSI is set the application program has to call this routine
 */
void nfft_precompute_lin_psi(nfft_plan *ths);

/** precomputes the values psi in a tensor product form
 *  if PRE_PSI is set the application program has to call this routine
 *  after setting the nodes x
 */
void nfft_precompute_psi(nfft_plan *ths);

/** precomputes the values psi and their indices in non tensor product form
 *  if PRE_FULL_PSI is set the application program has to call this routine
 */
void nfft_precompute_full_psi(nfft_plan *ths);

/** finalisation of transform
 */
void nfft_finalize(nfft_plan *ths);

/* @} 
 */


#endif
/* nfft.h */
