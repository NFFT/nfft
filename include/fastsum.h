/*! \file fastsum.h
 *  \brief Header file for the fast NFFT-based summation algorithm.
 *
 *  reference: M. Fenn, G. Steidl,
 *    Fast NFFT based summation of radial functions.
 *    Sampl. Theory Signal Image Process., 3, 1-28, 2004.
 *
 *  \author Markus Fenn
 *  \date 2003-2006
 */

#ifndef fastsum_h_inc
#define fastsum_h_inc

#include <complex.h>
#include "util.h"
#include "nfft3.h"

/**
 * Constant symbols
 */

/**********************************************************************
 * plan for fast summation algorithm
 **********************************************************************/
typedef struct fastsum_plan_
{
  /** api */

  int d;                                /**< number of dimensions            */

  int N_total;                          /**< number of source knots          */
  int M_total;                          /**< number of target knots          */

  complex *alpha;                       /**< source coefficients             */
  complex *f;                           /**< target evaluations              */

  double *x;                            /**< source knots in [-1/4,1/4]      */
  double *y;                            /**< target knots in [-1/4,1/4]      */

  complex (*kernel)(double , int , const double *);  /**< kernel function    */
  double *kernel_param;                 /**< parameters for kernel function  */

  unsigned flags;                       /**< flags precomp. and approx.type  */

  /** internal */

  int Ad;                               /**< number of spline knots for nearfield computation of regularized kernel*/
  double *Add;                          /**< spline values */

  /** DS_PRE - direct summation */
  complex *pre_K;                       /**< precomputed K(x_j-y_l)          */

  /** FS__ - fast summation */
  int n;                                /**< expansion degree                */
  fftw_complex *b;                      /**< expansion coefficients          */

  int p;                                /**< degree of smoothness of regularization */
  double eps_I;                         /**< inner boundary */
  double eps_B;                         /**< outer boundary */

  nfft_plan mv1;                        /**< source nfft plan                */
  nfft_plan mv2;                        /**< target nfft plan                */

  /* near field??? */

  /* things for computing *b - are they used only once?? */
  fftw_plan fft_plan;
} fastsum_plan;

/** initialize fast summation plan */
void fastsum_init_guru(fastsum_plan *ths, int d, int N_total, int M_total, complex (*kernel)(), double *param, unsigned flags, int nn, int p);

/** finalize plan */
void fastsum_finalize(fastsum_plan *ths);

/** direct summation */
void fastsum_exact(fastsum_plan *ths);

/** sort source nodes, precompute Fourier coefficients, etc. */
void fastsum_precompute(fastsum_plan *ths);

/** fast NFFT-based summation algorithm */
void fastsum_trafo(fastsum_plan *ths);

#endif
/* fastsum.h */
