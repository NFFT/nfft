/*! \file nnfft.h
 *  \brief Header file for the library.
 *
 * Includes simple and fast computation of the NNDFT (direct problem) as well
 * as (iterative) solution to the inverse problem.
 * authors: D. Potts, S. Kunis (c) 2002,2003
 */

#ifndef nnfft_h_inc
#define nnfft_h_inc

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "nfft.h"

/** 
 * Constant symbols for precomputation and memory usage (direct problem)
 */
#define PRE_PHI_HUT      (1U<< 0)
#define PRE_LIN_PSI      (1U<< 2)
#define PRE_PSI          (1U<< 3)
#define PRE_FULL_PSI     (1U<< 4)
#define MALLOC_X         (1U<< 5)
#define MALLOC_F_HAT     (1U<< 6)
#define MALLOC_F         (1U<< 7)
#define FFTW_INIT        (1U<< 8)
#define FFT_OUT_OF_PLACE (1U<< 9)
#define MALLOC_V         (1U<< 10)

/**
 * Structs for transformations 
 */
#ifndef nnfft_plan_h
#define nnfft_plan_h

typedef struct nnfft_plan_ 
{
  /** api */
  int M_total;                          /**< number of frequencies            */
  int N_total;                          /**< number of nodes                  */
  int d;                                /**< dimension, rank                  */
  double *sigma;                        /**< oversampling-factor 1            */
  double *a;                                /**< 1 + 2*m1/N1                      */ 
  int *N;                               /**< fftw-length = sigma1*N           */ 
  int *N1;
  int *aN1;
  int m;                                /**< cut-off parameter in time-domain */
  double *b;                            /**< shape parameters                 */
  int K;                                /**< K+1 values of phi                */

  int aN1_L;
  
  nfft_plan *direct_plan;                /**< plan for the nfft                */
  unsigned nnfft_flags;                 /**< flags for precomputation, malloc */

  int *n;                               /**< fftw-length = a*N            */
  
  double *x;                            /**< nodes (in time/spatial domain)   */
  
  double *v;                            /**< nodes (in fourier domain)   */

  fftw_complex *f_hat;                  /**< fourier coefficients, equispaced */
  fftw_complex *f;                      /**< samples at nodes x               */

  double *c_phi_inv;                    /**< precomputed data, matrix D       */
  double *psi;                          /**< precomputed data, matrix B       */
  int size_psi;                         /**< only for thin B                  */
  int *psi_index_g;                     /**< only for thin B                  */
  int *psi_index_f;                     /**< only for thin B                  */

  fftw_complex *F;

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
/* nnfft.h */
