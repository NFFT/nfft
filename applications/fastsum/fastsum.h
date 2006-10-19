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

/** 
 * \defgroup applications_fastsum Fast summation
 * \ingroup applications
 * \{
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

  double *x;                            /**< source knots in d-ball with radius 1/4-eps_b/2 */
  double *y;                            /**< target knots in d-ball with radius 1/4-eps_b/2 */

  complex (*kernel)(double , int , const double *);  /**< kernel function    */
  double *kernel_param;                 /**< parameters for kernel function  */

  unsigned flags;                       /**< flags precomp. and approx.type  */

  /** internal */

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

  /** near field */
  int Ad;                               /**< number of spline knots for nearfield computation of regularized kernel */
  double *Add;                          /**< spline values */

  /* things for computing *b - are they used only once?? */
  fftw_plan fft_plan;
} fastsum_plan;

/** initialize fast summation plan
 *
 * \param ths The pointer to a fastsum plan.
 * \param d The dimension of the problem.
 * \param N_total The number of source knots x.
 * \param M_total The number of target knots y.
 * \param kernel The kernel function.
 * \param param The parameters for the kernel function.
 * \param flags Fastsum flags.
 * \param nn The expansion degree.
 * \param m The cut-off parameter for the NFFT.
 * \param p The degree of smoothness.
 * \param eps_I The inner boundary.
 * \param eps_B the outer boundary.
 *
 */
void fastsum_init_guru(fastsum_plan *ths, int d, int N_total, int M_total, complex (*kernel)(), double *param, unsigned flags, int nn, int m, int p, double eps_I, double eps_B);

/** finalize plan
 *
 * \param ths The pointer to a fastsum plan.
 */
void fastsum_finalize(fastsum_plan *ths);

/** direct summation
 *
 * \param ths The pointer to a fastsum plan.
 */
void fastsum_exact(fastsum_plan *ths);

/** sort source nodes, precompute Fourier coefficients, etc.
 *
 * \param ths The pointer to a fastsum plan.
 */
void fastsum_precompute(fastsum_plan *ths);

/** fast NFFT-based summation algorithm
 *
 * \param ths The pointer to a fastsum plan.
 */
void fastsum_trafo(fastsum_plan *ths);
/* \} */

#endif
/* fastsum.h */
