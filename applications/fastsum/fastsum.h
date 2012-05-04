/*
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* $Id$ */

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

#include "config.h"

/** Include header for C99 complex datatype. */
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
/** Include header for utils from NFFT3 library. */
#include "nfft3util.h"
/** Include header for NFFT3 library. */
#include "nfft3.h"

#if !(defined(NF_LIN) || defined(NF_QUADR) || defined(NF_KUB))
  #define NF_KUB
#endif

#if !(defined(NF_ST) || defined(NF_BO))
  #define NF_ST
#endif

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef double _Complex (*kernel)(double , int , const double *);

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

  double _Complex *alpha;                       /**< source coefficients             */
  double _Complex *f;                           /**< target evaluations              */

  double *x;                            /**< source knots in d-ball with radius 1/4-eps_b/2 */
  double *y;                            /**< target knots in d-ball with radius 1/4-eps_b/2 */

  kernel k;  /**< kernel function    */
  double *kernel_param;                 /**< parameters for kernel function  */

  unsigned flags;                       /**< flags precomp. and approx.type  */

  /** internal */

  /** DS_PRE - direct summation */
  double _Complex *pre_K;                       /**< precomputed K(x_j-y_l)          */

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
  double _Complex *Add;                 /**< spline values */

  /* things for computing *b - are they used only once?? */
  fftw_plan fft_plan;

#ifdef NF_BO
  int box_count;
  int box_count_per_dim;
  int *box_offset;
  double *box_x;
  double _Complex *box_alpha;
#endif

  double MEASURE_TIME_t[8]; /**< Measured time for each step if MEASURE_TIME is set */
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
void fastsum_init_guru(fastsum_plan *ths, int d, int N_total, int M_total, kernel k, double *param, unsigned flags, int nn, int m, int p, double eps_I, double eps_B);

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

double _Complex regkern(kernel k, double xx, int p, const double *param, double a, double b);

/** cubic spline interpolation in near field with even kernels */
double _Complex kubintkern(const double x, const double _Complex *Add,
  const int Ad, const double a);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
/* fastsum.h */
