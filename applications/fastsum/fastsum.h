/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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
 * Direct and fast summation (convolution)
 * 
 * Computes the sums
 * \f[
 *   f(y_j) = \sum_{k=1}^N \alpha_k K(x_k-y_j),\quad   j=1\dots M.
 * \f]
 */

#ifndef fastsum_h_inc
#define fastsum_h_inc

#include "config.h"

/** Include header for C99 complex datatype. */
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
/** Include header for utils from NFFT3 library. */
/** Include header for NFFT3 library. */
#include "nfft3.h"
#include "infft.h"

#undef X
#define X(name) NFFT(name)

#if !(defined(NF_LIN) || defined(NF_QUADR) || defined(NF_KUB))
  #define NF_KUB
#endif

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef C (*kernel)(R , int , const R *);

/**
 * Constant symbols
 */
#define EXACT_NEARFIELD  (1U<< 0)

#define NEARFIELD_BOXES (1U<< 1)

/** If this flag is set, and eps_I > 0.0 and NEARFIELD_BOXES is not set,
 * then the vector permutation_x_alpha is stored. */
#define STORE_PERMUTATION_X_ALPHA (1U<< 2)

/** plan for fast summation algorithm */
typedef struct fastsum_plan_
{
  /** api */

  int d;                                /**< number of dimensions            */

  int N_total;                          /**< number of source knots          */
  int M_total;                          /**< number of target knots          */

  C *alpha;                       /**< source coefficients             */
  C *f;                           /**< target evaluations              */

  R *x;                            /**< source knots in d-ball with radius 1/4-eps_b/2 */
  R *y;                            /**< target knots in d-ball with radius 1/4-eps_b/2 */

  kernel k;  /**< kernel function    */
  R *kernel_param;                 /**< parameters for kernel function  */

  unsigned flags;                       /**< flags precomp. and approx.type  */

  /** internal */

  /** DS_PRE - direct summation */
  C *pre_K;                       /**< precomputed K(x_j-y_l)          */

  /** FS__ - fast summation */
  int n;                                /**< expansion degree                */
  C *b;                      /**< expansion coefficients          */
  C *f_hat;  /**< Fourier coefficients of nfft plans */

  int p;                                /**< degree of smoothness of regularization */
  R eps_I;                         /**< inner boundary                  */  /* fixed to p/n so far  */
  R eps_B;                         /**< outer boundary                  */  /* fixed to 1/16 so far */

  X(plan) mv1;                        /**< source nfft plan                */
  X(plan) mv2;                        /**< target nfft plan                */

  /** near field */
  int Ad;                               /**< number of spline knots for nearfield computation of regularized kernel */
  C *Add;                 /**< spline values */

  /* things for computing *b - are they used only once?? */
  FFTW(plan) fft_plan;

  int box_count;
  int box_count_per_dim;
  int *box_offset;
  R *box_x;
  C *box_alpha;
  
  int *permutation_x_alpha;    /**< permutation vector of source nodes if STORE_PERMUTATION_X_ALPHA is set */

  R MEASURE_TIME_t[8]; /**< Measured time for each step if MEASURE_TIME is set */

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
void fastsum_init_guru(fastsum_plan *ths, int d, int N_total, int M_total, kernel k, R *param, unsigned flags, int nn, int m, int p, R eps_I, R eps_B);

/** initialize node independent part of fast summation plan
 *
 * \param ths The pointer to a fastsum plan.
 * \param d The dimension of the problem.
 * \param kernel The kernel function.
 * \param param The parameters for the kernel function.
 * \param flags Fastsum flags.
 * \param nn The expansion degree.
 * \param p The degree of smoothness.
 * \param eps_I The inner boundary.
 * \param eps_B the outer boundary.
 *
 */
void fastsum_init_guru_kernel(fastsum_plan *ths, int d, kernel k, R *param,
    unsigned flags, int nn, int p, R eps_I, R eps_B);

/** initialize source nodes dependent part of fast summation plan
 *
 * \param ths The pointer to a fastsum plan.
 * \param N_total The number of source knots x.
 * \param nn_oversampled The oversampled expansion degree for nfft.
 * \param m The cut-off parameter for the NFFT.
 *
 */
void fastsum_init_guru_source_nodes(fastsum_plan *ths, int N_total, int nn_oversampled, int m);

/** initialize target nodes dependent part of fast summation plan
 *
 * \param ths The pointer to a fastsum plan.
 * \param M_total The number of target knots y.
 * \param nn_oversampled The oversampled expansion degree for nfft.
 * \param m The cut-off parameter for the NFFT.
 *
 */
void fastsum_init_guru_target_nodes(fastsum_plan *ths, int M_total, int nn_oversampled, int m);

/** finalize plan
 *
 * \param ths The pointer to a fastsum plan.
 */
void fastsum_finalize(fastsum_plan *ths);

/** finalize source nodes dependent part of plan
 *
 * \param ths The pointer to a fastsum plan.
 */
void fastsum_finalize_source_nodes(fastsum_plan *ths);

/** finalize target nodes dependent part of plan
 *
 * \param ths The pointer to a fastsum plan.
 */
void fastsum_finalize_target_nodes(fastsum_plan *ths);

/** finalize node independent part of plan
 *
 * \param ths The pointer to a fastsum plan.
 */
void fastsum_finalize_kernel(fastsum_plan *ths);

/** direct summation
 *
 * \param ths The pointer to a fastsum plan.
 */
void fastsum_exact(fastsum_plan *ths);

/** sort source nodes, precompute nfft source plan.
 *
 * \param ths The pointer to a fastsum plan.
 */
void fastsum_precompute_source_nodes(fastsum_plan *ths);

/** precompute nfft target plan.
 *
 * \param ths The pointer to a fastsum plan.
 */
void fastsum_precompute_target_nodes(fastsum_plan *ths);

/** sort source nodes, precompute nfft plans etc.
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

C regkern(kernel k, R xx, int p, const R *param, R a, R b);

/** cubic spline interpolation in near field with even kernels */
C kubintkern(const R x, const C *Add,
  const int Ad, const R a);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
/* fastsum.h */
