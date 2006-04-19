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

/** Include header for C99 complex datatype. */
#include <complex.h>
/** Include header for utils from NFFT3 library. */
#include "util.h"
/** Include header for NFFT3 library. */
#include "nfft3.h"

/**
 * Constant symbols
 */
#define EXACT_NEARFIELD  (1U<< 0)

/** plan for fast summation algorithm */
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

  unsigned flags;                       /**< flags precomp. and approx.type  */ /* not used so far */

  /** internal */

  int Ad;                               /**< number of spline knots for nearfield computation of regularized kernel*/
  double *Add;                          /**< spline values */

  /** DS_PRE - direct summation */
  complex *pre_K;                       /**< precomputed K(x_j-y_l)          */

  /** FS__ - fast summation */
  int n;                                /**< expansion degree                */
  fftw_complex *b;                      /**< expansion coefficients          */

  int p;                                /**< degree of smoothness of regularization */
  double eps_I;                         /**< inner boundary                  */  /* fixed to p/n so far  */
  double eps_B;                         /**< outer boundary                  */  /* fixed to 1/16 so far */

  nfft_plan mv1;                        /**< source nfft plan                */
  nfft_plan mv2;                        /**< target nfft plan                */

  /* near field??? */

  /* things for computing *b - are they used only once?? */
  fftw_plan fft_plan;
} fastsum_plan;

/** initialize fast summation plan */
void fastsum_init_guru(fastsum_plan *ths, int d, int N_total, int M_total, complex (*kernel)(), double *param, unsigned flags, int nn, int m, int p);

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
