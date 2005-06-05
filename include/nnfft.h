/*! THIS HEADER WILL BE MERGED WITH NFFT3.H
 */

#ifndef nnfft_h_inc
#define nnfft_h_inc

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "nfft3.h"

/** 
 * Constant symbols for precomputation and memory usage (direct problem)
 */

#define FFTW_INIT        (1U<< 23)
#define MALLOC_V         (1U<< 24)


#define MACRO_MV_PLAN(float_type)			                      \
  int N_total;                          /**< total number of Fourier coeffs.*/\
  int M_total;                          /**< total number of samples        */\
  float_type *f_hat;                    /**< Fourier coefficients           */\
  float_type *f;                        /**< samples                        */\
  
/**
 * Structs for transformations 
 */
#ifndef nnfft_plan_h
#define nnfft_plan_h

typedef struct nnfft_plan_ 
{
  /** api */
  MACRO_MV_PLAN(complex);

  int d;                                /**< dimension, rank                  */
  double *sigma;                        /**< oversampling-factor 1            */
  double *a;                            /**< 1 + 2*m1/N1                      */ 
  int *N;                               /**< cut-off-frequencies              */
  int *N1;
  int *aN1;
  int m;                                /**< cut-off parameter in time-domain */
  double *b;                            /**< shape parameters                 */
  int K;                                /**< K+1 values of phi                */

  int aN1_total;
  
  nfft_plan *direct_plan;               /**< plan for the nfft                */
  unsigned nnfft_flags;                 /**< flags for precomputation, malloc */

  int *n;                               /**<  = a*N                           */
  
  double *x;                            /**< nodes (in time/spatial domain)   */
  
  double *v;                            /**< nodes (in fourier domain)   */

  double *c_phi_inv;                    /**< precomputed data, matrix D       */
  double *psi;                          /**< precomputed data, matrix B       */
  int size_psi;                         /**< only for thin B                  */
  int *psi_index_g;                     /**< only for thin B                  */
  int *psi_index_f;                     /**< only for thin B                  */

  complex *F;

  double *spline_coeffs;                /**< input for de Boor algorithm, if  
                                             B_SPLINE or SINC_2m is defined   */
} nnfft_plan;

/** @defgroup nndft Group for direct nndft
 * ths group contains routines for 
 * direct discrete Fourier transform at nonequispaced knots in time and fouier domain
 * @{ 
 */

/** see formula (???), computes
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{l = 0}^{N-1} f_hat[k_l] * exp(-2 (pi) k_l x[j])
 */
void nndft_trafo(nnfft_plan *ths_plan);

/** see formula (???), computes
 * for j=0,...,M-1                                                             
 *  f[j] = sum_{l = 0}^{N-1} f_hat[k_l] * exp(+2 (pi) k_l x[j])
 */
void nndft_conjugated(nnfft_plan *ths_plan);

void nndft_adjoint(nnfft_plan *ths_plan);

void nndft_transposed(nnfft_plan *ths_plan);

/** @defgroup nfft Group for direct nfft
 * ths group contains routines for 
 * direct fast Fourier transform at nonequispaced knots
 * @{ 
 */
 
 /** initialisation for direct transform, simple interface
 */
void nnfft_init(nnfft_plan *ths_plan, int d, int N_total, int M_total, int *N);

void nnfft_init_specific(nnfft_plan *ths_plan, int d, int N_total, int M_total, int *N, int *N1,
                        int m, unsigned nnfft_flags);
                        
/** finalisation for direct transform
 */
void nnfft_finalize(nnfft_plan *ths_plan);

void nnfft_trafo(nnfft_plan *ths_plan);

void nnfft_conjugated(nnfft_plan *ths_plan);

void nnfft_adjoint(nnfft_plan *ths_plan);

void nnfft_transposed(nnfft_plan *ths_plan);

void nnfft_precompute_lin_psi(nnfft_plan *ths_plan, int K);

void nnfft_precompute_psi(nnfft_plan *ths_plan);

void nnfft_full_psi(nnfft_plan *ths_plan, double eps);

void nnfft_precompute_phi_hut(nnfft_plan *ths_plan);
#endif

#endif
/* nnfft.h */
