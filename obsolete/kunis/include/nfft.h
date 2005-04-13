/*! \file nfft.h
 *  \brief Header file for the library.
 *
 * Includes simple and fast computation of the NDFT (direct problem) as well
 * as (iterative) solution to the inverse problem.
 * authors: D. Potts, S. Kunis (c) 2002,2003
 */

#ifndef nfft_h_inc
#define nfft_h_inc

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "utils.h"

/** 
 * Constant symbols for precomputation and memory usage (direct problem)
 */
#define PRE_PHI_HUT      (1U<< 0)
#define PRE_PSI          (1U<< 1)
#define PRE_FULL_PSI     (1U<< 2)
#define MALLOC_X         (1U<< 3)
#define MALLOC_F_HAT     (1U<< 4)
#define MALLOC_F         (1U<< 5)
#define FFT_OUT_OF_PLACE (1U<< 6)

/** 
 * Constant symbols for precomputation and memory usage (inverse problem)
 */
#define LANDWEBER           (1U<< 0)
#define STEEPEST_DESCENT    (1U<< 1)
#define CGNR_E   	    (1U<< 2)
#define CGNE_R 		    (1U<< 3)
#define ITERATE_2nd         (1U<< 4)
#define NORMS_FOR_LANDWEBER (1U<< 5)
#define PRECOMPUTE_WEIGHT   (1U<< 6)
#define PRECOMPUTE_DAMP     (1U<< 7)

/**
 * Structs for transformations 
 */
#ifndef nfft_plan_h
#define nfft_plan_h
typedef struct nfft_plan_ 
{
  /** api */
  int d;                                /**< dimension, rank                  */
  int *N;                               /**< cut-off-frequencies              */
  int N_L;                              /**< total number of frequencies      */
  int M;                                /**< number of nodes                  */
  double *sigma;	                /**< oversampling-factor              */
  int *n;                               /**< fftw-length = sigma*N            */
  int n_L;                              /**< total size of fftw               */	
  int m;                                /**< cut-off parameter in time-domain */
  double *b;                            /**< shape parameters                 */

  unsigned nfft_flags;                  /**< flags for precomputation, malloc */
  unsigned fftw_flags;                  /**< flags for the fftw               */

  double *x;                            /**< nodes (in time/spatial domain)   */

  fftw_complex *f_hat;                  /**< fourier coefficients, equispaced */
  fftw_complex *f;                      /**< samples at nodes x               */

  /** internal */
  fftw_plan  my_fftw_plan1;             /**< fftw_plan forward                */
  fftw_plan  my_fftw_plan2;		/**< fftw_plan backward               */

  double **c_phi_inv;                   /**< precomputed data, matrix D       */
  double *psi;				/**< precomputed data, matrix B       */
  int size_psi;                         /**< only for thin B                  */
  int *psi_index_g;                     /**< only for thin B                  */
  int *psi_index_f;                     /**< only for thin B                  */

  fftw_complex *g;
  fftw_complex *g_hat;
  fftw_complex *g1;                     /**< input of fftw                    */
  fftw_complex *g2;                     /**< output of fftw                   */

  double *spline_coeffs;          	/**< input for de Boor algorithm, if  
					     B_SPLINE or SINC_2m is defined   */
} nfft_plan;

typedef struct infft_plan_ 
{
  nfft_plan *direct_plan;		/**< matrix vector multiplication     */
  unsigned infft_flags;			/**< iteration type, damping, weights */
  
  double *w;                         	/**< weighting factors                */
  double *w_hat;                     	/**< damping factors                  */

  /* right hand side */
  fftw_complex *given_f;                /**< right hand side, sampled values  */

  /* solutions */ 
  fftw_complex *f_hat_iter;             /**< iterative solution               */
  fftw_complex *f_hat_iter_2nd;         /**< iterative solution, 2nd algo     */

  /* vectors */
  fftw_complex *r_iter;                 /**< iterated original residual vector*/ 
  fftw_complex *z_hat_iter;             /**< residual vector of normal eq.    */
  fftw_complex *p_hat_iter;             /**< (non damped) search direction    */
  fftw_complex *v_iter;                 /**< residual vector update in CGNR   */

  /* factors */
  double alpha_iter;                 /**< step size for search direction   */
  double alpha_iter_2nd;             /**< step size for search direction   */
  double beta_iter;                  /**< step size for search correction  */
  double gamma_iter;                 /**< factor in CGNE_R                 */
  double gamma_iter_old;             /**< factor in CGNE_R                 */

  /* dot products */
  double dot_r_iter;                 /**< dotproductc{_w}(r_iter)          */
  double dot_r_iter_old;             /**< old dotproductc{_w}(r_iter)      */
  double dot_z_hat_iter;             /**< dotproductc{_w}(z_hat_iter)      */
  double dot_z_hat_iter_old;         /**< old dotproductc{_w}(z_hat_iter)  */
  double dot_p_hat_iter;             /**< dotproductc{_w}(p_hat_iter)      */
  double dot_v_iter;                 /**< dotproductc{_w}(v_iter)          */  
} infft_plan;
#endif



/** @defgroup ndft Group for direct ndft
 * This group contains routines for 
 * direct discrete Fourier transform at nonequispaced knots
 * @{ 
 */

/** see formula (1.1), computes
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 */
void ndft_trafo(nfft_plan *this_plan);

/** see formula (1.1), computes
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 */
void ndft_conjugated(nfft_plan *this_plan);

/** see formula (1.2), computes
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 */
void ndft_adjoint(nfft_plan *this_plan);

/** see formula (1.2), computes
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(-2(pi) k x[j])
 */
void ndft_transposed(nfft_plan *this_plan);
/** @} 
 */ 

/** @defgroup nfft Group for direct nfft
 * This group contains routines for 
 * direct fast Fourier transform at nonequispaced knots
 * @{ 
 */

/** wrapper for nfft_init and d=1
 */
void nfft_init_1d(nfft_plan *this_plan, int N1, int M);

/** wrapper for nfft_init and d=2
 */
void nfft_init_2d(nfft_plan *this_plan, int N1, int N2, int M);

/** wrapper for nfft_init and d=3
 */
void nfft_init_3d(nfft_plan *this_plan, int N1, int N2, int N3, int M);

/** initialisation for direct transform, simple interface
 */
void nfft_init(nfft_plan *this_plan, int d, int *N, int M);

/** initialisation for direct transform, specific interface
 */
void nfft_init_specific(nfft_plan *this_plan, int d, int *N, int M, int *n,
	                int m, unsigned nfft_flags, unsigned fftw_flags);

/** precomputes the values psi
 *  if the PRE_PSI is set the application program has to call this routine
 *  after setting the nodes this_plan->x
 */
void nfft_precompute_psi(nfft_plan *this_plan);

/** precomputes the values psi and their indices in non tensor form
 *  if the PRE_PSI and PRE_FULL_PSI is set the application program has to call
 *  this routine after calling nfft_precompute_psi
 */
void nfft_full_psi(nfft_plan *this_plan, double eps);

/** see formula (1.1), computes in a fast and approximative way
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(-2 (pi) k x[j])
 */
void nfft_trafo(nfft_plan *this_plan);

/** see formula (1.1), computes in a fast and approximative way
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{k in I_N^d} f_hat[k] * exp(+2 (pi) k x[j])
 */
void nfft_conjugated(nfft_plan *this_plan);

/** see formula (1.2), computes in a fast and approximative way
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 */
void nfft_adjoint(nfft_plan *this_plan);

/** see formula (1.2), computes in a fast and approximative way
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * exp(+2(pi) k x[j])
 */
void nfft_transposed(nfft_plan *this_plan);

/** finalisation for direct transform
 */
void nfft_finalize(nfft_plan *this_plan);			
/** @} 
 */ 

/** @defgroup infft Group for inverse nfft
 * This group contains routines for 
 * inverse fast Fourier transform at nonequispaced knots
 * @{ 
 */

/** initialisation for inverse transform, simple interface
 */
void infft_init(infft_plan *this_iplan, nfft_plan *direct_plan);

/** initialisation for inverse transform, specific interface
 */
void infft_init_specific(infft_plan *this_iplan, nfft_plan *direct_plan,
			 int infft_flags);

/** computes start residuals
 */
void infft_before_loop(infft_plan *this_iplan);

/** iterates one step
 */
void infft_loop_one_step(infft_plan *this_iplan);

/** finalisation for inverse transform
 */
void infft_finalize(infft_plan *this_iplan);
/** @} 
 */ 

/** debug functions */
void plan_info(nfft_plan *this_plan);

#endif
/* nfft.h */
