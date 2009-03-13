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

/*! \file nfft3.h
 *  \brief Header file for the nfft3 library.
 */
#ifndef __NFFT3_H__
#define __NFFT3_H__

/** Include header for FFTW3 library for its complex type. */
#include <fftw3.h>

/* config header */
#include "nfftconf.h"

/* Malloc and free functions */
extern void *nfft_malloc(size_t n);
extern void nfft_free(void *p);
extern void nfft_die(const char *s);

/* Malloc and free hooks */
typedef void *(*nfft_malloc_type_function) (size_t n);
typedef void  (*nfft_free_type_function) (void *p);
typedef void  (*nfft_die_type_function) (const char *errString);
extern nfft_malloc_type_function nfft_malloc_hook;
extern nfft_free_type_function nfft_free_hook;
extern nfft_die_type_function nfft_die_hook;

/** Macros for public members inherited by all plan structures. */
#define MACRO_MV_PLAN(float_type)                                             \
  int N_total;                          /**< Total number of Fourier          \
					     coefficients                   */\
  int M_total;                          /**< Total number of samples        */\
                                                                              \
  float_type *f_hat;                    /**< Vector of Fourier coefficients,  \
                                             size is N_total float_types    */\
  float_type *f;                        /**< Vector of samples,               \
				             size is M_total float types    */\
  void (*mv_trafo)(void*);              /**< Pointer to the own transform   */\
  void (*mv_adjoint)(void*);            /**< Pointer to the own adjoint     */\

typedef struct
{
  MACRO_MV_PLAN(fftw_complex)
} mv_plan_complex;

typedef struct
{
  MACRO_MV_PLAN(double)
} mv_plan_double;

/** Macros for window functions. */
#if defined(DIRAC_DELTA)
  #define PHI_HUT(k,d) 1.0
  #define PHI(x,d) (fabs((x))<10e-8)? 1.0 : 0.0
  #define WINDOW_HELP_INIT(d)
  #define WINDOW_HELP_FINALIZE
  #define WINDOW_HELP_ESTIMATE_m {ths->m = 0;}
#elif defined(GAUSSIAN)
  #define PHI_HUT(k,d) ((double)exp(-(pow(PI*(k)/ths->n[d],2.0)*ths->b[d])))
  #define PHI(x,d) ((double)exp(-pow((x)*ths->n[d],2.0)/ ths->b[d])/sqrt(PI*ths->b[d]))
  #define WINDOW_HELP_INIT \
    {                                                                          \
      int WINDOW_idx;                                                          \
      ths->b = (double*) nfft_malloc(ths->d*sizeof(double));                   \
      for(WINDOW_idx=0; WINDOW_idx<ths->d; WINDOW_idx++)                       \
      ths->b[WINDOW_idx]=((double)2*ths->sigma[WINDOW_idx])/                   \
        (2*ths->sigma[WINDOW_idx]-1)*(((double)ths->m) / PI);                  \
      }
  #define WINDOW_HELP_FINALIZE {nfft_free(ths->b);}
  #define WINDOW_HELP_ESTIMATE_m {ths->m =12;}
#elif defined(B_SPLINE)
  #define PHI_HUT(k,d) ((double)(((k)==0)? 1.0/ths->n[(d)] :                   \
    pow(sin((k)*PI/ths->n[(d)])/((k)*PI/ths->n[(d)]),2*ths->m)/ths->n[(d)]))
  #define PHI(x,d) (nfft_bspline(2*ths->m,((x)*ths->n[(d)])+                   \
    (double)ths->m,ths->spline_coeffs)/ths->n[(d)])
  #define WINDOW_HELP_INIT \
    {                                                                          \
      ths->spline_coeffs= (double*)nfft_malloc(2*ths->m*sizeof(double));       \
    }
  #define WINDOW_HELP_FINALIZE {nfft_free(ths->spline_coeffs);}
  #define WINDOW_HELP_ESTIMATE_m {ths->m =11;}
#elif defined(SINC_POWER)
  #define PHI_HUT(k,d) (nfft_bspline(2*ths->m,((double)2*ths->m*(k))/          \
    ((2*ths->sigma[(d)]-1)*ths->n[(d)]/ths->sigma[(d)])+ (double)ths->m,       \
    ths->spline_coeffs))
  #define PHI(x,d) ((double)(ths->n[(d)]/ths->sigma[(d)]*(2*ths->sigma[(d)]-1)/\
    (2*ths->m)*pow(nfft_sinc(PI*ths->n[(d)]/ths->sigma[(d)]*(x)*               \
    (2*ths->sigma[(d)]-1)/(2*ths->m)),2*ths->m)/ths->n[(d)]))
  #define WINDOW_HELP_INIT \
    {                                                                          \
      ths->spline_coeffs= (double*)nfft_malloc(2*ths->m*sizeof(double));       \
    }
  #define WINDOW_HELP_FINALIZE {nfft_free(ths->spline_coeffs);}
  #define WINDOW_HELP_ESTIMATE_m {ths->m = 9;}
#else /* Kaiser-Bessel is the default. */
  #define PHI_HUT(k,d) ((double)nfft_i0( ths->m*sqrt(\
    pow((double)(ths->b[d]),2.0) - pow(2.0*PI*(k)/ths->n[d],2.0))))
  #define PHI(x,d) ((double)((pow((double)(ths->m),2.0)\
    -pow((x)*ths->n[d],2.0))>0)? \
    sinh(ths->b[d]*sqrt(pow((double)(ths->m),2.0)-                             \
    pow((x)*ths->n[d],2.0)))/(PI*sqrt(pow((double)(ths->m),2.0)-               \
    pow((x)*ths->n[d],2.0))): (((pow((double)(ths->m),2.0)-                    \
    pow((x)*ths->n[d],2.0))<0)? sin(ths->b[d]*                                 \
    sqrt(pow(ths->n[d]*(x),2.0)-pow((double)(ths->m),2.0)))/                   \
    (PI*sqrt(pow(ths->n[d]*(x),2.0)-pow((double)(ths->m),2.0))):1.0))
  #define WINDOW_HELP_INIT \
    {                                                                          \
      int WINDOW_idx;                                                          \
      ths->b = (double*) nfft_malloc(ths->d*sizeof(double));                   \
      for(WINDOW_idx=0; WINDOW_idx<ths->d; WINDOW_idx++)                       \
      ths->b[WINDOW_idx] = ((double)PI*(2.0-1.0/ths->sigma[WINDOW_idx]));      \
  }
  #define WINDOW_HELP_FINALIZE {nfft_free(ths->b);}
  #define WINDOW_HELP_ESTIMATE_m {ths->m = 6;}
#endif

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/**
 * @defgroup nfft NFFT - Nonequispaced fast Fourier transform
 * Direct and fast computation of the nonequispaced discrete Fourier transform.
 * @{
 *
 * This module implements the nonequispaced fast Fourier transforms.
 * In the following, we abbreviate the term "nonequispaced fast Fourier
 * transform" by NFFT.
 *
 * We introduce our notation and nomenclature for discrete Fourier transforms.
 * Let the torus
 * \f[
 *   \mathbb{T}^d
 *    := \left\{ \mathbf{x}=\left(x_t\right)_{t=0,\hdots,d-1}\in\mathbb{R}^{d}:
 *    \; - \frac{1}{2} \le x_t < \frac{1}{2},\; t=0,\hdots,d-1 \right\}
 * \f]
 * of dimension \f$d\f$ be given.
 * It will serve as domain from which the nonequispaced nodes \f$\mathbf{x}\f$
 * are taken.
 * The sampling set is given by \f${\cal X}:=\{\mathbf{x}_j \in {\mathbb T}^d:
 * \,j=0,\hdots,M-1\}\f$.
 * Possible frequencies \f$\mathbf{k}\f$ are collected in the multi index set
 * \f[
 *   I_{\mathbf{N}} := \left\{ \mathbf{k}=\left(k_t\right)_{t=0,\hdots,d-1}\in
 *   \mathbb{Z}^d: - \frac{N_t}{2} \le k_t < \frac{N_t}{2} ,\;t=0,\hdots,d-1
 * \right\}.
 * \f]
 *
 * Our concern is the computation of the
 * \e nonequispaced discrete Fourier transform \e (NDFT) \anchor ndft_formula
 * \f[
 * f_j = \sum_{\mathbf{k}\in I_{\mathbf{N}}}
 * \hat{f}_{\mathbf{k}} {\rm e}^{-2\pi{\rm i} \mathbf{k}\mathbf{x}_j}, \qquad
 * j=0,\hdots,M-1.
 * \f]
 * The corresponding adjoint NDFT is the computation of
 * \f[
 *   \hat f_{\mathbf{k}}=\sum_{j=0}^{M-1} f_j {\rm e}^{+2\pi{\rm i}
 *    \mathbf{k}\mathbf{x}_j}, \qquad \mathbf{k}\in I_{\mathbf{N}}.
 * \f]
 * Direct implementations are given by \ref ndft_trafo and \ref ndft_adjoint
 * taking \f${\cal O}(|I_{\mathbf{N}}|M)\f$ floating point operations.
 * Approximative realisations take only
 * \f${\cal O}(|I_{\mathbf{N}}|\log|I_{\mathbf{N}}|+M)\f$ floating point operations.
 * These are provided by \ref nfft_trafo and \ref nfft_adjoint, respectively.
 */

/* Planner flags, i.e. constant symbols for precomputation and memory usage */

/**
 * If this flag is set, the deconvolution step (the multiplication with the
 * diagonal matrix \f$\mathbf{D}\f$) uses precomputed values of the Fourier
 * transformed window function.
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \author Stefan Kunis
 */
#define PRE_PHI_HUT      (1U<< 0)

/**
 * If this flag is set, the convolution step (the multiplication with the
 * sparse matrix \f$\mathbf{B}\f$) uses particular properties of the Gaussian
 * window function to trade multiplications for direct calls to exponential
 * function.
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \author Stefan Kunis
 */
#define FG_PSI           (1U<< 1)

/**
 * If this flag is set, the convolution step (the multiplication with the
 * sparse matrix \f$\mathbf{B}\f$) uses linear interpolation from a lookup
 * table of equispaced samples of the window function instead of exact values
 * of the window function.
 * At the moment a table of size \f$(K+1)d\f$ is used, where
 * \f$K=2^{10}(m+1)\f$.
 * An estimate for the size of the lookup table with respect to the target
 * accuracy should be implemented.
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \author Stefan Kunis
 */
#define PRE_LIN_PSI      (1U<< 2)

/**
 * If this flag is set, the convolution step (the multiplication with the
 * sparse matrix \f$\mathbf{B}\f$) uses particular properties of the Gaussian
 * window function to trade multiplications for direct calls to exponential
 * function (the remaining \f$2dM\f$ direct calls are precomputed).
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \author Stefan Kunis
 */
#define PRE_FG_PSI       (1U<< 3)

/**
 * If this flag is set, the convolution step (the multiplication with the
 * sparse matrix \f$\mathbf{B}\f$) uses \f$(2m+2)dM\f$ precomputed values of
 * the window function.
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \author Stefan Kunis
 */
#define PRE_PSI          (1U<< 4)

/**
 * If this flag is set, the convolution step (the multiplication with the
 * sparse matrix \f$\mathbf{B}\f$) uses \f$(2m+2)^dM\f$ precomputed values of
 * the window function, in addition indices of source and target vectors are
 * stored.
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \author Stefan Kunis
 */
#define PRE_FULL_PSI     (1U<< 5)

/**
 * If this flag is set, (de)allocation of the node vector is done.
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \see nfft_finalize
 * \author Stefan Kunis
 */
#define MALLOC_X         (1U<< 6)

/**
 * If this flag is set, (de)allocation of the vector of Fourier coefficients is
 * done.
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \see nfft_finalize
 * \author Stefan Kunis
 */
#define MALLOC_F_HAT     (1U<< 7)

/**
 * If this flag is set, (de)allocation of the vector of samples is done.
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \see nfft_finalize
 * \author Stefan Kunis
 */
#define MALLOC_F         (1U<< 8)

/**
 * If this flag is set, FFTW uses disjoint input/output vectors.
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \see nfft_finalize
 * \author Stefan Kunis
 */
#define FFT_OUT_OF_PLACE (1U<< 9)

/**
 * If this flag is set, fftw_init/fftw_finalize is called.
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \see nfft_finalize
 * \author Stefan Kunis
 */
#define FFTW_INIT        (1U<< 10)

/**
 * Summarises if precomputation is used within the convolution step (the
 * multiplication with the sparse matrix \f$\mathbf{B}\f$).
 * If testing against this flag is positive, \ref nfft_precompute_one_psi has
 * to be called.
 *
 * \see nfft_init
 * \see nfft_init_advanced
 * \see nfft_init_guru
 * \see nfft_precompute_one_psi
 * \see nfft_finalize
 * \author Stefan Kunis
 */
#define PRE_ONE_PSI (PRE_LIN_PSI| PRE_FG_PSI| PRE_PSI| PRE_FULL_PSI)

/** Structure for a NFFT plan */
typedef struct
{
  /* api */
  MACRO_MV_PLAN(fftw_complex)

  int d;                                /**< Dimension, rank                 */
  int *N;                               /**< Multi bandwidth                 */
  double *sigma;                        /**< Oversampling-factor             */
  int *n;                               /**< FFTW length, equal to sigma*N,
					     default is the power of 2 such
					     that \f$2\le\sigma<4\f$         */
  int n_total;                          /**< Total size of FFTW              */
  int m;                                /**< Cut-off parameter of the window
                                             function, default value is
					     6 (KAISER_BESSEL),
					     9 (SINC_POWER),
					     11 (B_SPLINE),
					     12 (GAUSSIAN)                   */
  double *b;                            /**< Shape parameter of the window
                                             function                        */
  int K;                                /**< Number of equispaced samples of
					     the window function for \ref
					     PRE_LIN_PSI                     */

  unsigned nfft_flags;                  /**< Flags for precomputation,
                                             (de)allocation, and FFTW usage,
					     default setting is
					     PRE_PHI_HUT| PRE_PSI|
					     MALLOC_X| MALLOC_F_HAT| MALLOC_F|
					     FFTW_INIT| FFT_OUT_OF_PLACE     */

  unsigned fftw_flags;                  /**< Flags for the FFTW, default is
					     FFTW_ESTIMATE| FFTW_DESTROY_INPUT
					                                     */

  double *x;                            /**< Nodes in time/spatial domain,
					     size is \f$dM\f$ doubles        */

  double MEASURE_TIME_t[3];             /**< Measured time for each step if
					     MEASURE_TIME is set             */

  /* internal*/
  fftw_plan  my_fftw_plan1;             /**< Forward FFTW plan               */
  fftw_plan  my_fftw_plan2;             /**< Backward FFTW plan              */

  double **c_phi_inv;                   /**< Precomputed data for the
                                             diagonal matrix \f$D\f$, size is
					     \f$N_0+\hdots+N_{d-1}\f$ doubles*/
  double *psi;                          /**< Precomputed data for the sparse
                                             matrix \f$B\f$, size depends on
					     precomputation scheme           */
  int *psi_index_g;                     /**< Indices in source/target vector
					     for \ref PRE_FULL_PSI           */
  int *psi_index_f;                     /**< Indices in source/target vector
					     for \ref PRE_FULL_PSI           */

  fftw_complex *g;                    /**< Oversampled vector of samples,
					     size is \ref n_total double
					     complex                         */
  fftw_complex *g_hat;                /**< Zero-padded vector of Fourier
                                             coefficients, size is
					     \ref n_total fftw_complex     */
  fftw_complex *g1;                   /**< Input of fftw                   */
  fftw_complex *g2;                   /**< Output of fftw                  */

  double *spline_coeffs;                /**< Input for de Boor algorithm if
                                             B_SPLINE or SINC_POWER is
                                             defined                         */
} nfft_plan;


/**
 * Computes a NDFT, see the \ref ndft_formula "definition".
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void ndft_trafo(nfft_plan *ths);

/**
 * Computes an adjoint NDFT, see the \ref ndftH_formula "definition".
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void ndft_adjoint(nfft_plan *ths);

/**
 * Computes a NFFT, see the \ref ndft_formula "definition".
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_trafo(nfft_plan *ths);
void nfft_trafo_1d(nfft_plan *ths);
void nfft_trafo_2d(nfft_plan *ths);
void nfft_trafo_3d(nfft_plan *ths);

/**
 * Computes an adjoint NFFT, see the \ref ndftH_formula "definition".
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_adjoint(nfft_plan *ths);
void nfft_adjoint_1d(nfft_plan *ths);
void nfft_adjoint_2d(nfft_plan *ths);
void nfft_adjoint_3d(nfft_plan *ths);

/**
 * Initialisation of a transform plan, wrapper d=1.
 *
 * \arg ths The pointer to a nfft plan
 * \arg N1 bandwidth
 * \arg M The number of nodes
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init_1d(nfft_plan *ths, int N1, int M);

/**
 * Initialisation of a transform plan, wrapper d=2.
 *
 * \arg ths The pointer to a nfft plan
 * \arg N1 bandwidth
 * \arg N2 bandwidth
 * \arg M The number of nodes
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init_2d(nfft_plan *ths, int N1, int N2, int M);

/**
 * Initialisation of a transform plan, wrapper d=3.
 *
 * \arg ths The pointer to a nfft plan
 * \arg N1 bandwidth
 * \arg N2 bandwidth
 * \arg N3 bandwidth
 * \arg M The number of nodes
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init_3d(nfft_plan *ths, int N1, int N2, int N3, int M);

/**
 * Initialisation of a transform plan, simple.
 *
 * \arg ths The pointer to a nfft plan
 * \arg d The dimension
 * \arg N The multi bandwidth
 * \arg M The number of nodes
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init(nfft_plan *ths, int d, int *N, int M);

/**
 * Initialisation of a transform plan, advanced.
 * NOT YET IMPLEMENTED!!
 *
 * \arg ths The pointer to a nfft plan
 * \arg d The dimension
 * \arg N The multi bandwidth
 * \arg M The number of nodes
 * \arg nfft_flags_on NFFT flags to switch on
 * \arg nfft_flags_off NFFT flags to switch off
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init_advanced(nfft_plan *ths, int d, int *N, int M,
                        unsigned nfft_flags_on, unsigned nfft_flags_off);

/**
 * Initialisation of a transform plan, guru.
 *
 * \arg ths The pointer to a nfft plan
 * \arg d The dimension
 * \arg N The multi bandwidth
 * \arg M The number of nodes
 * \arg n The oversampled multi bandwidth
 * \arg m The spatial cut-off
 * \arg nfft_flags NFFT flags to use
 * \arg fftw_flags_off FFTW flags to use
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_init_guru(nfft_plan *ths, int d, int *N, int M, int *n,
                    int m, unsigned nfft_flags, unsigned fftw_flags);


/**
 * Precomputation for a transform plan.
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis
 *
 * wrapper for precompute*_psi
 *
 * if PRE_*_PSI is set the application program has to call this routine
 * (after) setting the nodes x
 */
void nfft_precompute_one_psi(nfft_plan *ths);

/**
 * Superceded by nfft_precompute_one_psi.
 * \author Stefan Kunis
 */
void nfft_precompute_full_psi(nfft_plan *ths);

/**
 * Superceded by nfft_precompute_one_psi.
 * \author Stefan Kunis
 */
void nfft_precompute_psi(nfft_plan *ths);

/**
 * Superceded by nfft_precompute_one_psi.
 * \author Stefan Kunis
 */
void nfft_precompute_lin_psi(nfft_plan *ths);

/**
 * Checks a transform plan for frequently used bad parameter.
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_check(nfft_plan *ths);

/**
 * Destroys a transform plan.
 *
 * \arg ths The pointer to a nfft plan
 *
 * \author Stefan Kunis, Daniel Potts
 */
void nfft_finalize(nfft_plan *ths);

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup nfsct NFCT/NFST - Nonequispaced fast (co)sine transform
 * Direct and fast computation of the discrete nonequispaced (co)sine
 * transform.
 * @{
 */

#ifdef HAVE_NFCT

/** Structure for a transform plan */
typedef struct
{
  /* api */
  MACRO_MV_PLAN(double)

  int d;                                /**< dimension, rank                  */
  int *N;                               /**< cut-off-frequencies (kernel)     */
  int *n;                               /**< length of dct-i                  */
  double *sigma;                        /**< oversampling-factor              */
  int m;                                /**< cut-off parameter in time-domain */

  double nfct_full_psi_eps;
  double *b;                            /**< shape parameters                 */

  unsigned nfct_flags;                  /**< flags for precomputation, malloc */
  unsigned fftw_flags;                  /**< flags for the fftw               */

  double *x;                            /**< nodes (in time/spatial domain)   */

  double MEASURE_TIME_t[3];             /**< measured time for each step      */

  /** internal */
  fftw_plan  my_fftw_r2r_plan;          /**< fftw_plan                        */
  fftw_r2r_kind *r2r_kind;              /**< r2r transform type (dct-i)       */

  double **c_phi_inv;                   /**< precomputed data, matrix D       */
  double *psi;                          /**< precomputed data, matrix B       */
  int size_psi;                         /**< only for thin B                  */
  int *psi_index_g;                     /**< only for thin B                  */
  int *psi_index_f;                     /**< only for thin B                  */

  double *g;
  double *g_hat;
  double *g1;                           /**< input of fftw                    */
  double *g2;                           /**< output of fftw                   */

  double *spline_coeffs;                /**< input for de Boor algorithm, if
                                             B_SPLINE or SINC_2m is defined   */
} nfct_plan;


/**
 * Creates a 1-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth \f$N\f$
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfct_init_1d( nfct_plan *ths_plan, int N0, int M_total);

/**
 * Creates a 3-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth of dimension 1
 * \arg N1 The bandwidth of dimension 2
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfct_init_2d( nfct_plan *ths_plan, int N0, int N1, int M_total);

/**
 * Creates a 3-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth of dimension 1
 * \arg N1 The bandwidth of dimension 2
 * \arg N2 The bandwidth of dimension 3
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfct_init_3d( nfct_plan *ths_plan, int N0, int N1, int N2, int M_total);

/**
 * Creates a d-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d the dimension
 * \arg N The bandwidths
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfct_init( nfct_plan *ths_plan, int d, int *N, int M_total);

/**
 * Creates a d-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d the dimension
 * \arg N The bandwidths
 * \arg M_total The number of nodes \f$x\f$
 * \arg n The oversampled bandwidths
 * \arg m The cut-off parameter
 * \arg nfct_flags The flags known to nfct
 * \arg fftw_flags The flags known to fftw
 *
 * \author Steffen Klatt
 */
void nfct_init_guru( nfct_plan *ths_plan, int d, int *N, int M_total, int *n,
                         int m, unsigned nfct_flags, unsigned fftw_flags);

/**
 * precomputes the values psi
 * if the PRE_PSI is set the application program has to call this routine
 * after setting the nodes this_plan->x
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfct_precompute_psi( nfct_plan *ths_plan);


/**
 * executes a NFCT (approximate,fast), computes for \f$j=0,...,M\_total-1\f$
 * \f$f_j^C(x_j) = \sum_{k \in I_0^{N,d}} \hat{f}_k^C * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfct_trafo( nfct_plan *ths_plan);

/**
 * executes a NDCT (exact,slow), computes for \f$j=0,...,M\_total-1\f$
 * \f$f_j^C(x_j) = \sum_{k \in I_0^{N,d}} \hat{f}_k^C * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void ndct_trafo( nfct_plan *ths_plan);

/**
 * executes a transposed NFCT (approximate,fast), computes for \f$k \in I_0^{N,d}\f$
 * \f$h^C(k) = \sum_{j \in I_0^{(M\_total,1)}} f_j^C * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfct_adjoint( nfct_plan *ths_plan);

/**
 * executes a direct transposed NDCT (exact,slow), computes for \f$k \in I_0^{N,d}\f$
 * \f$h^C(k) = \sum_{j \in I_0^{(M\_total,1)}} f_j^C * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void ndct_adjoint( nfct_plan *ths_plan);

/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfct_finalize( nfct_plan *ths_plan);

/**
 * do some adjustments (N,n) then compute PHI_HUT
 *
 * \arg ths_plan the plan for the transform
 * \arg k        index of c_phi
 * \arg d        dimension
 *
 * \author Steffen Klatt
 */
double nfct_phi_hut( nfct_plan *ths_plan, int k, int d);

/**
 * do some adjustments (N,n) then compute PHI
 *
 * \arg ths_plan the plan for the transform
 * \arg x        node \f$x\f$
 * \arg d        dimension
 *
 * \author Steffen Klatt
 */
double nfct_phi ( nfct_plan *ths_plan, double x, int d);

/**
 * returns 2(n-1),  fftw related issue
 *
 * \arg n       i.e. length of dct-1
 *
 * \author Steffen Klatt
 */
int nfct_fftw_2N( int n);

/**
 * returns 0.5n+1,  fftw related issue
 *
 * \arg n       i.e. length of dct-1
 *
 * \author Steffen Klatt
 */
int nfct_fftw_2N_rev(int n);

#endif

/*###########################################################################*/

#ifdef HAVE_NFST

/** Structure for a transform plan */
typedef struct
{
  /* api */
  MACRO_MV_PLAN(double)

  int d;                                /**< dimension, rank                  */
  int *N;                               /**< bandwidth                        */
  int *n;                               /**< length of dst-1                  */
  double *sigma;                        /**< oversampling-factor              */
  int m;                                /**< cut-off parameter in time-domain */

  double nfst_full_psi_eps;
  double *b;                            /**< shape parameters                 */

  unsigned nfst_flags;                  /**< flags for precomputation, malloc */
  unsigned fftw_flags;                  /**< flags for the fftw               */

  double *x;                            /**< nodes (in time/spatial domain)   */

  double MEASURE_TIME_t[3];             /**< measured time for each step     */

  /** internal */
  fftw_plan  my_fftw_r2r_plan;         /**< fftw_plan forward                */
  fftw_r2r_kind *r2r_kind;              /**< r2r transform type (dct-i)       */

  double **c_phi_inv;                   /**< precomputed data, matrix D       */
  double *psi;                          /**< precomputed data, matrix B       */
  int size_psi;                         /**< only for thin B                  */
  int *psi_index_g;                     /**< only for thin B                  */
  int *psi_index_f;                     /**< only for thin B                  */


  double *g;
  double *g_hat;
  double *g1;                           /**< input of fftw                    */
  double *g2;                           /**< output of fftw                   */

  double *spline_coeffs;                /**< input for de Boor algorithm, if
                                             B_SPLINE or SINC_2m is defined   */
} nfst_plan;


/**
 * Creates a 1-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth \f$N\f$
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfst_init_1d( nfst_plan *ths_plan, int N0, int M_total);

/**
 * Creates a 3-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth of dimension 1
 * \arg N1 The bandwidth of dimension 2
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfst_init_2d( nfst_plan *ths_plan, int N0, int N1, int M_total);

/**
 * Creates a 3-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N0 The bandwidth of dimension 1
 * \arg N1 The bandwidth of dimension 2
 * \arg N2 The bandwidth of dimension 3
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfst_init_3d( nfst_plan *ths_plan, int N0, int N1, int N2, int M_total);

/**
 * Creates a d-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d the dimension
 * \arg N The bandwidths
 * \arg M_total The number of nodes \f$x\f$
 *
 * \author Steffen Klatt
 */
void nfst_init( nfst_plan *ths_plan, int d, int *N, int M_total);

/**
 * Creates a d-dimensional transform plan with pcific m.
 * (just for convenience)
 *
 * \arg ths_plan The plan for the transform
 * \arg d the dimension
 * \arg N The bandwidths
 * \arg M_total The number of nodes \f$x\f$
 * \arg m cut-off parameter
 *
 * \author Steffen Klatt
 */
void nfst_init_m( nfst_plan *ths_plan, int d, int *N, int M_total, int m);

/**
 * Creates a d-dimensional transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d the dimension
 * \arg N The bandwidths
 * \arg M_total The number of nodes \f$x\f$
 * \arg n The oversampled bandwidths
 * \arg m The cut-off parameter
 * \arg nfst_flags The flags known to nfst
 * \arg fftw_flags The flags known to fftw
 *
 * \author Steffen Klatt
 */
void nfst_init_guru( nfst_plan *ths_plan, int d, int *N, int M_total, int *n,
                     int m, unsigned nfst_flags, unsigned fftw_flags);

/**
 * precomputes the values psi
 * if the PRE_PSI is set the application program has to call this routine
 * after setting the nodes this_plan->x
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfst_precompute_psi( nfst_plan *ths_plan);

/**
 * executes a NFST (approximate,fast), computes for \f$j=0,...,M\_total-1\f$
 * \f$f_j^S(x_j) = \sum_{k \in I_1^{N,d}} \hat{f}_k^S * sin(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfst_trafo( nfst_plan *ths_plan);

/**
 * executes a NDST (exact,slow), computes for \f$j=0,...,M\_total-1\f$
 * \f$f_j^S(x_j) = \sum_{k \in I_1^{N,d}} \hat{f}_k^S * sin(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void ndst_trafo( nfst_plan *ths_plan);



/**
 * executes a transposed NFST (approximate,fast), computes for \f$k \in I_1^{N,d}\f$
 * \f$h^S(k) = \sum_{j \in I_0^{M\_total,1}} f_j^S * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfst_adjoint( nfst_plan *ths_plan);

/**
 * executes a direct transposed NDST (exact,slow), computes for \f$k \in I_1^{N,d}\f$
 * \f$h^S(k) = \sum_{j \in I_0^{M\_total,1}} f_j^S * cos(2 \pi k x_j)\f$
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void ndst_adjoint( nfst_plan *ths_plan);

/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan for the transform
 *
 * \author Steffen Klatt
 */
void nfst_finalize( nfst_plan *ths_plan);

/**
 *
 *
 * \arg ths_plan The plan for the transform
 *
 */
void nfst_full_psi( nfst_plan *ths_plan, double eps);

/**
 * do some adjustments (N,n) then compute PHI_HUT
 *
 * \arg ths_plan the plan for the transform
 * \arg k        index of c_phi
 * \arg d        dimension
 *
 * \author Steffen Klatt
 */
double nfst_phi_hut( nfst_plan *ths_plan, int k, int d);

/**
 * do some adjustments (N,n) then compute PHI
 *
 * \arg ths_plan the plan for the transform
 * \arg x        node \f$x\f$
 * \arg d        dimension
 *
 * \author Steffen Klatt
 */
double nfst_phi ( nfst_plan *ths_plan, double x, int d);

/**
 * returns 2(n+1),  fftw related issue
 *
 * \arg n       i.e. length of dst-1
 *
 * \author Steffen Klatt
 */
int nfst_fftw_2N( int n);

/**
 * returns 0.5n-1,  fftw related issue
 *
 * \arg n       i.e. length of dct-1
 *
 * \author Steffen Klatt
 */
int nfst_fftw_2N_rev( int n);

#endif

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup nnfft NNFFT - Nonequispaced in time and frequency FFT
 * Direct and fast computation of the discrete nonequispaced in time and
 * frequency Fourier transform.
 * @{
 */

#ifdef HAVE_NNFFT

/**
 * If this flag is set, (de)allocation of the frequency node vector is done.
 *
 * \see nnfft_init
 * \see nnfft_init_guru
 * \see nnfft_finalize
 * \author Tobias Knopp
 */
#define MALLOC_V         (1U<< 11)

/** Structure for a transform plan */
typedef struct
{
  /* api */
  MACRO_MV_PLAN(fftw_complex)

  int d;                                /**< dimension, rank                 */
  double *sigma;                        /**< oversampling-factor             */
  double *a;                            /**< 1 + 2*m/N1                      */
  int *N;                               /**< cut-off-frequencies             */
  int *N1;                              /**< sigma*N                         */
  int *aN1;                             /**< sigma*a*N                       */
  int m;                                /**< cut-off parameter in time-domain*/
  double *b;                            /**< shape parameters                */
  int K;                                /**< number of precomp. uniform psi  */

  int aN1_total;                        /**< aN1_total=aN1[0]* ... *aN1[d-1] */

  nfft_plan *direct_plan;               /**< plan for the nfft               */
  unsigned nnfft_flags;                 /**< flags for precomputation, malloc*/
  int *n;                               /**< n=N1, for the window function   */

  double *x;                            /**< nodes (in time/spatial domain)  */
  double *v;                            /**< nodes (in fourier domain)       */

  double *c_phi_inv;                    /**< precomputed data, matrix D      */
  double *psi;                          /**< precomputed data, matrix B      */
  int size_psi;                         /**< only for thin B                 */
  int *psi_index_g;                     /**< only for thin B                 */
  int *psi_index_f;                     /**< only for thin B                 */
  fftw_complex *F;

  double *spline_coeffs;                /**< input for de Boor algorithm, if
                                             B_SPLINE or SINC_2m is defined  */
} nnfft_plan;


/**
 * Creates a transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d The dimension
 * \arg N_total The number of nodes \f$v\f$
 * \arg M_total The number of nodes \f$x\f$
 * \arg N The bandwidth \f$N\f$
 *
 * \author Tobias Knopp
 */
void nnfft_init(nnfft_plan *ths_plan, int d, int N_total, int M_total, int *N);

/**
 * Creates a transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg d The dimension
 * \arg N_total The number of nodes \f$v\f$
 * \arg M_total The number of nodes \f$x\f$
 * \arg N The bandwidth \f$N\f$
 * \arg N1 The oversampled bandwidth \f$N\f$
 * \arg m The cut-off parameter
 * \arg nnfft_flags The flags
 *
 * \author Tobias Knopp
 */
void nnfft_init_guru(nnfft_plan *ths_plan, int d, int N_total, int M_total,
                     int *N, int *N1, int m, unsigned nnfft_flags);

/**
 * Executes a direct NNDFT, i.e. computes for \f$j=0,...,M_{total}-1\f$
 * \f[
 *   f(x_j) = \sum_{k = 0}^{N_{total}-1} \hat{f}(v_k) {\rm e}^{-2 \pi
 *            \mbox{\rm\scriptsize i} v_k x_j \odot N}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void nndft_trafo(nnfft_plan *ths_plan);

/**
 * Executes a direct adjoint NNDFT, i.e. computes for \f$k=0,...,N_{total}-1\f$
 * \f[
 *   \hat{f}(v_k) = \sum_{j = 0}^{M_{total}-1} f(x_j) {\rm e}^{2 \pi
 *                  \mbox{\rm\scriptsize i} v_k x_j \odot N}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void nndft_adjoint(nnfft_plan *ths_plan);

/**
 * Executes a NNFFT, i.e. computes for \f$j=0,...,M_{total}-1\f$
 * \f[
 *   f(x_j) = \sum_{k = 0}^{N_{total}-1} \hat{f}(v_k) {\rm e}^{-2 \pi
 *            \mbox{\rm\scriptsize i} v_k x_j \odot N}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void nnfft_trafo(nnfft_plan *ths_plan);

/**
 * Executes a adjoint NNFFT, i.e. computes for \f$k=0,...,N_{total}-1\f$
 * \f[
 *   \hat{f}(v_k) = \sum_{j = 0}^{M_{tota}l-1} f(x_j) {\rm e}^{2 \pi
 *                  \mbox{\rm\scriptsize i} v_k x_j \odot N}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void nnfft_adjoint(nnfft_plan *ths_plan);

/**
 * Precomputation for a transform plan.
 *
 * \arg ths_plan The pointer to a nfft plan
 *
 * \author Tobias Knopp
 *
 * precomputes equally spaced values of the window function psi
 *
 * if PRE_LIN_PSI is set the application program has to call this routine
 */
void nnfft_precompute_lin_psi(nnfft_plan *ths_plan);

/**
 * Precomputation for a transform plan.
 *
 * \arg ths_plan The pointer to a nfft plan
 *
 * \author Tobias Knopp
 *
 * precomputes the values of the window function psi in a tensor product form
 *
 * if PRE_PSI is set the application program has to call this routine after
 * setting the nodes x
 */
void nnfft_precompute_psi(nnfft_plan *ths_plan);

/**
 * Precomputation for a transform plan.
 *
 * \arg ths_plan The pointer to a nfft plan
 *
 * \author Tobias Knopp
 *
 * precomputes the values of the window function psi and their indices in
 * non tensor product form
 *
 * if PRE_FULL_PSI is set the application program has to call this routine
 * after setting the nodes x
 */
void nnfft_precompute_full_psi(nnfft_plan *ths_plan);

/**
 * Precomputation for a transform plan.
 *
 * \arg ths_plan The pointer to a nfft plan
 *
 * \author Tobias Knopp
 *
 * precomputes the values of the fourier transformed window function, i.e. phi_hut
 *
 * if PRE_PHI_HUT is set the application program has to call this routine
 * after setting the nodes v
 */
void nnfft_precompute_phi_hut(nnfft_plan *ths_plan);

/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void nnfft_finalize(nnfft_plan *ths_plan);

#endif

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup nsfft NSFFT - Nonequispaced sparse FFT
 * Direct and fast computation of the nonequispaced FFT on the hyperbolic
 * cross.
 * @{
 */

#ifdef HAVE_NSFFT

/**
 * If this flag is set, the member \ref index_sparse_to_full is (de)allocated
 * and initialised for the use in the routine \ref nsdft_trafo and
 * \ref nsdft_adjoint.
 *
 * \see nsfft_init
 * \author Stefan Kunis
 */
#define NSDFT            (1U<< 12)

/** Structure for a NFFT plan */
typedef struct
{
  MACRO_MV_PLAN(fftw_complex)

  int d;                                /**< dimension, rank; d=2,3          */
  int J;                                /**< problem size, i.e.,
                                             d=2: N_total=(J+4) 2^(J+1)
                                             d=3: N_total=2^J 6(2^((J+1)/2+1)
                                                          -1)+2^(3(J/2+1))   */
  int sigma;                            /**< oversampling-factor             */

  unsigned flags;                       /**< flags for precomputation, malloc*/

  int *index_sparse_to_full;            /**< index conversation,
                                             overflow for d=3, J=9!          */

  int r_act_nfft_plan;                  /**< index of current nfft block     */
  nfft_plan *act_nfft_plan;             /**< current nfft block              */
  nfft_plan *center_nfft_plan;          /**< central nfft block              */

  fftw_plan* set_fftw_plan1;            /**< fftw plan for the nfft blocks   */
  fftw_plan* set_fftw_plan2;            /**< fftw plan for the nfft blocks   */

  nfft_plan *set_nfft_plan_1d;          /**< nfft plans for short nffts      */
  nfft_plan *set_nfft_plan_2d;          /**< nfft plans for short nffts      */

  double *x_transposed;                 /**< coordinate exchanged nodes, d=2 */
  double *x_102,*x_201,*x_120,*x_021;   /**< coordinate exchanged nodes, d=3 */

} nsfft_plan;

/**
 * Executes an NSDFT, computes for \f$j=0,\hdots,M-1\f$:
 * \f[
 *   f_j = \sum_{k\in H_N^d}\hat f_k {\rm e}^{-2\pi{\rm\scriptsize i}k x_j}
 * \f]
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsdft_trafo(nsfft_plan *ths);

/**
 * Executes an adjoint NSFFT, computes for \f$k\in H_N^d\f$:
 * \f[
 *   \hat f_k = \sum_{j=0,\hdots,M-1} f_j {\rm e}^{+2\pi{\rm\scriptsize i}k x_j}
 * \f]
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Stefan Kunis
 */
void nsdft_adjoint(nsfft_plan *ths);

/**
 * Executes an NSFFT, computes \b fast and \b approximate for
 * \f$j=0,\hdots,M-1\f$:
 * \f[
 *   f_j = \sum_{k\in H_N^d}\hat f_k {\rm e}^{-2\pi{\rm\scriptsize i}k x_j}
 * \f]
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsfft_trafo(nsfft_plan *ths);

/**
 * Executes an adjoint NSFFT, computes \b fast and \b approximate for
 * \f$k\in H_N^d\f$:
 * \f[
 *   \hat f_k = \sum_{j=0,\hdots,M-1} f_j {\rm e}^{+2\pi{\rm\scriptsize i}k x_j}
 * \f]
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Stefan Kunis
 */
void nsfft_adjoint(nsfft_plan *ths);

/**
 * Copy coefficients from nsfft plan to a nfft plan.
 *
 * \arg ths Pointers to a nsfft plan and to a nfft plan
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsfft_cp(nsfft_plan *ths, nfft_plan *ths_nfft);

/**
 * Initialisation of pseudo random nodes and coefficients.
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsfft_init_random_nodes_coeffs(nsfft_plan *ths);

/**
 * Initialisation of a transform plan.
 *
 * \arg ths The pointer to a nsfft plan
 * \arg d The dimension
 * \arg J The problem size
 * \arg M The number of nodes
 * \arg m nfft cut-off parameter
 * \arg flags
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsfft_init(nsfft_plan *ths, int d, int J, int M, int m, unsigned flags);

/**
 * Destroys a transform plan.
 *
 * \arg ths The pointer to a nsfft plan
 *
 * \author Markus Fenn, Stefan Kunis
 */
void nsfft_finalize(nsfft_plan *ths);

#endif

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup mri MRI - Transforms in magnetic resonance imaging
 * @{
 */

#ifdef HAVE_MRI

/**
 * The structure for the transform plan.
 */
typedef struct
{
  /* api */
  MACRO_MV_PLAN(fftw_complex)

  nfft_plan plan;

  int N3;
  double sigma3;
  double *t;
  double *w;
} mri_inh_2d1d_plan;

/**
 * The structure for the transform plan.
 */
typedef struct
{
  /* api */
  MACRO_MV_PLAN(fftw_complex)

  nfft_plan plan;

  int N3;
  double sigma3;
  double *t;
  double *w;
} mri_inh_3d_plan;


/**
 * Executes a mri transformation considering the field inhomogeneity with the 2d1d method,
 * i.e. computes for \f$j=0,...,M_{total}-1\f$
 * \f[
 *   f(x_j) = \sum_{k \in I_N^2} \hat{f}(k) {\rm e}^{\mbox{\rm\scriptsize i} t_j \omega(k)}
 *                      {\rm e}^{-2 \pi \mbox{\rm\scriptsize i} k x_j}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_2d1d_trafo(mri_inh_2d1d_plan *ths);

/**
 * Executes an adjoint mri transformation considering the field inhomogeneity with the 2d1d method,
 * i.e. computes for \f$k \in I_N^2\f$
 * \f[
 *   \hat{f}(k) = \sum_{j=0}^{M_{total}-1} f(x_j) {\rm e}^{\mbox{\rm\scriptsize i} t_j \omega(k)}
 *                      {\rm e}^{-2 \pi \mbox{\rm\scriptsize i} k x_j}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_2d1d_adjoint(mri_inh_2d1d_plan *ths);

/**
 * Creates a transform plan.
 *
 * \arg ths_plan The plan for the transform
 * \arg N The bandwidth \f$N\f$
 * \arg M_total The number of nodes \f$x\f$
 * \arg n The oversampled bandwidth \f$N\f$
 * \arg m The cut-off parameter
 * \arg sigma The oversampling factor
 * \arg nnfft_flags The flags
 *
 * \author Tobias Knopp
 */
void mri_inh_2d1d_init_guru(mri_inh_2d1d_plan *ths, int *N, int M, int *n,
                    int m, double sigma, unsigned nfft_flags, unsigned fftw_flags);

/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_2d1d_finalize(mri_inh_2d1d_plan *ths);

/**
 * Executes a mri transformation considering the field inhomogeneity with the 3d method,
 * i.e. computes for \f$j=0,...,M_{total}-1\f$
 * \f[
 *   f(x_j) = \sum_{k \in I_N^2} \hat{f}(k) {\rm e}^{\mbox{\rm\scriptsize i} t_j \omega(k)}
 *                      {\rm e}^{-2 \pi \mbox{\rm\scriptsize i} k x_j}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_3d_trafo(mri_inh_3d_plan *ths);

/**
 * Executes an adjoint mri transformation considering the field inhomogeneity with the 3d method,
 * i.e. computes for \f$k \in I_N^2\f$
 * \f[
 *   \hat{f}(k) = \sum_{j=0}^{M_{total}-1} f(x_j) {\rm e}^{\mbox{\rm\scriptsize i} t_j \omega(k)}
 *                      {\rm e}^{-2 \pi \mbox{\rm\scriptsize i} k x_j}
 * \f]
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_3d_adjoint(mri_inh_3d_plan *ths);

void mri_inh_3d_init_guru(mri_inh_3d_plan *ths, int *N, int M, int *n,
                    int m, double sigma, unsigned nfft_flags, unsigned fftw_flags);

/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */
void mri_inh_3d_finalize(mri_inh_3d_plan *ths);

#endif

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/**
 * @defgroup nfsft NFSFT - Nonequispaced fast spherical Fourier transform
 * @{
 *
 * This module implements nonuniform fast spherical Fourier transforms. In the
 * following, we abbreviate the term "nonuniform fast spherical Fourier
 * transform" by NFSFT.
 *
 * \section Preliminaries
 * This section summarises basic definitions and properties related to spherical
 * Fourier transforms.
 *
 * \subsection sc Spherical Coordinates
 * Every point in \f$\mathbb{R}^3\f$ can be described in \e spherical \e
 * coordinates by a vector \f$(r,\vartheta,\varphi)^{\mathrm{T}}\f$ with the
 * radius \f$r \in \mathbb{R}^{+}\f$ and two angles \f$\vartheta \in [0,\pi]\f$,
 * \f$\varphi \in [-\pi,\pi)\f$.
 * We denote by \f$\mathbb{S}^2\f$ the two-dimensional unit sphere embedded
 * into \f$\mathbb{R}^3\f$, i.e.
 * \f[
 *   \mathbb{S}^2 := \left\{\mathbf{x} \in \mathbb{R}^{3}:\;
 *   \|\mathbf{x}\|_2=1\right\}
 * \f]
 * and identify a point from \f$\mathbb{S}^2\f$ with the corresponding vector
 * \f$(\vartheta,\varphi)^{\mathrm{T}}\f$. The
 * spherical coordinate system is illustrated in the following figure:
 * \image html sphere.png ""
 * \image latex sphere.pdf "" width=0.45\textwidth
 * For consistency with the other modules and the conventions used there, we
 * also use \e swapped \e scaled \e spherical \e coordinates \f$x_1 :=
 * \frac{\varphi}{2\pi}\f$, \f$x_2 := \frac{\vartheta}{2\pi}\f$ and identify a
 * point from \f$\mathbb{S}^2\f$ with the vector
 * \f$\mathbf{x} := \left(x_1,x_2\right) \in
 *  [-\frac{1}{2}, \frac{1}{2}) \times [0,\frac{1}{2}]\f$.
 *
 * \subsection lp Legendre Polynomials
 * The \e Legendre \e polynomials \f$P_k : [-1,1]
 * \rightarrow \mathbb{R}$, $k \in \mathbb{N}_{0}\f$ as \e classical \e
 * orthogonal \e polynomials are given by their corresponding \e Rodrigues \e
 * formula
 * \f[
 *   P_k(t) := \frac{1}{2^k k!} \frac{\text{d}^k}{\text{d} t^k}
 *   \left(t^2-1\right)^k.
 * \f]
 * The corresponding three-term recurrence relation is
 * \f[
 *   (k+1)P_{k+1}(t) = (2k+1) x P_{k}(t) - k P_{k-1}(t) \quad (k \in
 *   \mathbb{N}_0).
 * \f]
 * With
 * \f[
 *   \left< f,g \right>_{\text{L}^2\left([-1,1]\right)} :=
 *   \int_{-1}^{1} f(t) g(t) \text{d} t
 * \f]
 * being the usual \f$\text{L}^2\left([-1,1]\right)\f$ inner product,
 * the Legendre polynomials obey the orthogonality condition
 * \f[
 *   \left< P_k,P_l \right>_{\text{L}^2\left([-1,1]\right)} = \frac{2}{2k+1}
 *   \delta_{k,l}.
 * \f]
 *
 * \remark The normalisation constant \f$ c_k := \sqrt{\frac{2k+1}{2}}\f$
 * renders the scaled Legendre polynomials \f$c_k P_k\f$ orthonormal with
 * respect to the induced \f$\text{L}^2\left([-1,1]\right)\f$ norm
 * \f[
 *   \|f\|_{\text{L}^2\left([-1,1]\right)} :=
 *   \left(<f,f>_{\text{L}^2\left([-1,1]\right)}\right)^{1/2} =
 *   \left(\int_{-1}^{1} |f(t)|^2 \; \text{d} t\right)^{1/2}.
 * \f]
 *
 * \subsection alf Associated Legendre Functions
 * The \a associated \a Legendre \a functions \f$P_k^n : [-1,1] \rightarrow
 * \mathbb{R} \f$, \f$n \in \mathbb{N}_0\f$, \f$k \ge n\f$ are defined by
 * \f[
 *   P_k^n(t) := \left(\frac{(k-n)!}{(k+n)!}\right)^{1/2}
 *   \left(1-t^2\right)^{n/2} \frac{\text{d}^n}{\text{d} t^n} P_k(t).
 * \f]
 * For \f$n = 0\f$, they coincide with the Legendre polynomials, i.e.
 * \f$P_k^0 = P_k\f$.
 * The associated Legendre functions obey the three-term recurrence relation
 * \f[
 *   P_{k+1}^n(t) = v_{k}^n t P_k^n(t) + w_{k}^n P_{k-1}^n(t) \quad (k \ge n),
 * \f]
 * with \f$P_{n-1}^n(t) := 0\f$, \f$P_{n}^n(t) := \frac{\sqrt{(2n)!}}{2^n n!}
 * \left(1-t^2\right)^{n/2}\f$, and
 * \f[
 *   v_{k}^n := \frac{2k+1}{((k-n+1)(k+n+1))^{1/2}}\; ,\qquad
 *   w_{k}^n := - \frac{((k-n)(k+n))^{1/2}}{((k-n+1)(k+n+1))^{1/2}}.
 * \f]
 * For fixed \f$n\f$, the set \f$\left\{P_k^n:\: k
 * \ge n\right\}\f$ forms a complete set of orthogonal functions in
 * \f$\text{L}^2\left([-1,1]\right)\f$
 * with
 * \f[
 *   \left< P_k^n,P_l^n \right>_{\text{L}^2\left([-1,1]\right)} = \frac{2}{2k+1}
 *   \delta_{k,l} \quad (0 \le n \le k,l).
 * \f]
 *
 * \remark The normalisation constant \f$ c_k = \sqrt{\frac{2k+1}{2}}\f$
 * renders the scaled associated Legendre functions \f$c_k P_k^n\f$ orthonormal
 * with respect to the induced \f$\text{L}^2\left([-1,1]\right)\f$ norm
 * \f[
 *   \|f\|_{\text{L}^2\left([-1,1]\right)} :=
 *   \left(<f,f>_{\text{L}^2\left([-1,1]\right)}\right)^{1/2} =
 *   \left(\int_{-1}^{1} |f(t)|^2 \; \text{d} t\right)^{1/2}.
 * \f]
 *
 * \subsection sh Spherical Harmonics
 * The standard orthogonal basis of spherical harmonics for \f$\text{L}^2
 * \left(\mathbb{S}^2\right)\f$ with yet unnormalised basis functions
 * \f$\tilde{Y}_k^n : \mathbb{S}^2 \rightarrow \mathbb{C}\f$ is given by
 * \f[
 *   \tilde{Y}_k^n(\vartheta,\varphi) := P_k^{|n|}(\cos\vartheta)
 *   \mathrm{e}^{\mathrm{i} n \varphi}
 * \f]
 * with the usual \f$\text{L}^2\left(\mathbb{S}^2\right)\f$ inner product
 * \f[
 *   \left< f,g \right>_{\mathrm{L}^2\left(\mathbb{S}^2\right)} :=
 *   \int_{\mathbb{S}^2} f(\vartheta,\varphi) \overline{g(\vartheta,\varphi)}
 *   \: \mathrm{d} \mathbf{\xi} := \int_{-\pi}^{\pi} \int_{0}^{\pi}
 *   f(\vartheta,\varphi) \overline{g(\vartheta,\varphi)} \sin \vartheta
 *   \; \mathrm{d} \vartheta \; \mathrm{d} \varphi.
 * \f]
 * The normalisation constant \f$c_k^n := \sqrt{\frac{2k+1}{4\pi}}\f$ renders
 * the  scaled basis functions
 * \f[
 *   Y_k^n(\vartheta,\varphi) := c_k^n P_k^{|n|}(\cos\vartheta)
 *   \mathrm{e}^{\mathrm{i} n \varphi}
 * \f]
 * orthonormal with respect to the induced \f$\text{L}^2\left(\mathbb{S}^2
 * \right)\f$ norm
 * \f[
 *   \|f\|_{\text{L}^2\left(\mathbb{S}^2\right)} =
 *   \left(<f,f>_{\text{L}^2\left(\mathbb{S}^2\right)}\right)^{1/2} =
 *   \left(\int_{-\pi}^{\pi} \int_{0}^{\pi} |f(\vartheta,\varphi)|^2 \sin
 *   \vartheta \; \mathrm{d} \vartheta \; \mathrm{d} \varphi\right)^{1/2}.
 * \f]
 * A function \f$f \in \mathrm{L}^2\left(\mathbb{S}^2\right)\f$ has the
 * orthogonal expansion
 * \f[
 *   f = \sum_{k=0}^{\infty} \sum_{n=-k}^{k} \hat{f}(k,n) Y_k^n,
 * \f]
 * where the coefficients \f$\hat{f}(k,n) := \left< f, Y_k^{n}
 * \right>_{\mathrm{L}^2\left(\mathbb{S}^2\right)}\f$ are the \e spherical
 * \e Fourier \e coefficients and the equivalence is understood in the
 * \f$\mathrm{L}^2\f$-sense.
 *
 *
 * \section nfsfts Nonuniform Fast Spherical Fourier Transforms
 *
 * This section describes the input and output relation of the spherical
 * Fourier transform algorithms and the layout of the corresponding plan
 * structure.
 *
 * \subsection ndsft Nonuniform Discrete Spherical Fourier Transform
 * The \e nonuniform \e discrete \e spherical \e Fourier \e transform (\e NDSFT)
 * is defined as follows:
 * \f[
 *     \begin{array}{rcl}
 *       \text{\textbf{Input}} & : & \text{coefficients }
 *         \hat{f}(k,n) \in \mathbb{C} \text{ for } k=0,\ldots,N,\;n=-k,
 *         \ldots,k,\; N \in \mathbb{N}_0,\\[1ex]
 *                             &   & \text{arbitrary nodes } \mathbf{x}(m) \in
 *         [-\frac{1}{2},\frac{1}{2}] \times [0,\frac{1}{2}]
 *         \text{ for } m=0,\ldots,M-1, M \in \mathbb{N}. \\[1ex]
 *       \text{\textbf{Task}}  & : & \text{evaluate } f(m) := f\left(
 *       \mathbf{x}(m)\right) = \sum_{k=0}^N \sum_{n=-k}^k \hat{f}_k^n
 *         Y_k^n\left(\mathbf{x}(m)\right) \text{ for } m=0,\ldots,M-1.
 *         \\[1ex]
 *       \text{\textbf{Output}} & : & \text{coefficients } f(m) \in
 *         \mathbb{C} \text{ for } m=0,\ldots,M-1.\\
 *     \end{array}
 * \f]
 *
 * \subsection andsft Adjoint Nonuniform Discrete Spherical Fourier Transform
 * The \e adjoint \e nonuniform \e discrete \e spherical \e Fourier \e transform
 * (\e adjoint \e NDSFT)
 * is defined as follows:
 * \f[
 *     \begin{array}{rcl}
 *       \text{\textbf{Input}} & : & \text{coefficients } f(m) \in
 *         \mathbb{C} \text{ for } m=0,\ldots,M-1, M \in \mathbb{N},\\
 *                             &   & \text{arbitrary nodes } \mathbf{x}(m) \in
 *         [-\frac{1}{2},\frac{1}{2}] \times [0,\frac{1}{2}] \text{ for }
 *         m=0,\ldots,M-1, N \in \mathbb{N}_0.\\[1ex]
 *       \text{\textbf{Task}}  & : & \text{evaluate } \hat{f}(k,n)
 *         := \sum_{m=0}^{M-1} f(m) \overline{Y_k^n\left(\mathbf{x}(m)\right)}cd Do
 *         \text{ for } k=0,\ldots,N,\;n=-k,\ldots,k.\\[1ex]
 *       \text{\textbf{Output}} & : & \text{coefficients }
 *         \hat{f}(k,n) \in \mathbb{C} \text{ for }
 *         k=0,\ldots,N,\;n=-k,\ldots,k.\\[1ex]
 *     \end{array}
 * \f]
 *
 * \subsection dl Data Layout
 * This section describes the public  layout of the \ref nfsft_plan structure
 * which
 * contains all data for the computation of the aforementioned spherical Fourier
 * transforms. The structure contains private (no read or write allowed), public
 * read-only (only
 * read access permitted), and public read-write (read and write access allowed)
 * members. In the following, we indicate read and write access by \c read and
 * \c write. The public members are structured as follows:
 * \li \c N_total (\c read)
 *        The total number of components in \c f_hat. If the bandwidth is
 *        \f$N \in \mathbb{N}_0\f$, the total number of components in \c f_hat
 *        is \c N_total \f$ = (2N+2)^2\f$.
 * \li \c M_total (\c read)
 *        the total number of samples \f$M\f$
 * \li \c f_hat (\c read-write)
 *        The flattened array of spherical Fourier coefficents. The array
 *        has length \f$(2N+2)^2\f$ such that valid indices \f$i \in
 *        \mathbb{N}_0\f$ for array access \c f_hat \c[ \f$i\f$ \c] are
 *        \f$i=0,1,\ldots,(2N+2)^2-1\f$.
 *        However, only read and write access to indices corresponding to
 *        spherical Fourier coefficients \f$\hat{f}(k,n)\f$ is defined. The index
 *        \f$i\f$ corresponding to the spherical Fourier coefficient
 *        \f$\hat{f}(k,n)\f$ with \f$0 \le k \le M\f$, \f$-k \le n \le k\f$ is
 *        \f$i = (N+2)(N-n+1)+N+k+1\f$. For convenience, the helper macro
 *        \ref NFSFT_INDEX(k,n) provides the necessary index calculations such
 *        that
 *        one can write \c f_hat[ \c NFSFT_INDEX(\f$k,n\f$\c)] \c =
 *        \c ... to access
 *        the component corresponding to \f$\hat{f}(k,n)\f$.
 *        The data layout is due to implementation details.
 * \li \c f (\c read-write)
 *        the array of coefficients \f$f(m)\f$ for \f$m=0,\ldots,M-1\f$ such
 *        that \c f[\f$m\f$\c] = \f$f(m)\f$
 * \li \c N (\c read)
 *        the bandwidth \f$N \in \mathbb{N}_0\f$
 * \li \c x
 *        the array of nodes \f$\mathbf{x}(m) \in
 *        [-\frac{1}{2},\frac{1}{2}] \times [0,\frac{1}{2}]\f$ for \f$m = 0,
 *        \ldots,M-1\f$ such that \c f[\f$2m\f$\c] = \f$x_1\f$ and
 *        \c f[\f$2m+1\f$\c] = \f$x_2\f$
 *
 * \subsection gtn Good to know...
 * When using the routines of this module you should bear in mind the following:
 * \li The bandwidth \f$N_{\text{max}}\f$ up to which precomputation is
 *   performed is always chosen as the next power of two with respect to the
 *   specified maximum bandwidth.
 * \li By default, the NDSFT transforms (see \ref ndsft_trafo, \ref nfsft_trafo)
 *   are allowed to destroy the input \c f_hat while the input \c x is
 *   preserved. On the contrary, the adjoint NDSFT transforms
 *   (see \ref ndsft_adjoint, \ref nfsft_adjoint) do not destroy the input
 *   \c f and \c x by default. The desired behaviour can be assured by using the
 *   \ref NFSFT_PRESERVE_F_HAT, \ref NFSFT_PRESERVE_X, \ref NFSFT_PRESERVE_F and
 *   \ref NFSFT_DESTROY_F_HAT, \ref NFSFT_DESTROY_X, \ref NFSFT_DESTROY_F
 *   flags.
 */

#ifdef HAVE_NFSFT

/* Planner flags */

/**
 * By default, all computations are performed with respect to the
 * unnormalized basis functions
 * \f[
 *   \tilde{Y}_k^n(\vartheta,\varphi) = P_k^{|n|}(\cos\vartheta)
 *   \mathrm{e}^{\mathrm{i} n \varphi}.
 * \f]
 * If this flag is set, all computations are carried out using the \f$L_2\f$-
 * normalized basis functions
 * \f[
 *   Y_k^n(\vartheta,\varphi) = \sqrt{\frac{2k+1}{4\pi}} P_k^{|n|}(\cos\vartheta)
 *   \mathrm{e}^{\mathrm{i} n \varphi}.
 * \f]
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_NORMALIZED    (1U << 0)

/**
 * If this flag is set, the fast NFSFT algorithms (see \ref nfsft_trafo,
 * \ref nfsft_adjoint) will use internally the exact but usually slower direct
 * NDFT algorithm in favor of fast but approximative NFFT algorithm.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_USE_NDFT      (1U << 1)

/**
 * If this flag is set, the fast NFSFT algorithms (see \ref nfsft_trafo,
 * \ref nfsft_adjoint) will use internally the usually slower direct
 * DPT algorithm in favor of the fast FPT algorithm.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 * \warning This feature is not implemented yet!
 */
#define NFSFT_USE_DPT       (1U << 2)

/**
 * If this flag is set, the init methods (see \ref nfsft_init , \ref
 * nfsft_init_advanced , and \ref nfsft_init_guru) will allocate memory and the
 * method \ref nfsft_finalize will free the array \c x for you. Otherwise,
 * you have to assure by yourself that \c x points to an array of
 * proper size before excuting a transform and you are responsible for freeing
 * the corresponding memory before program termination.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_MALLOC_X      (1U << 3)

/**
 * If this flag is set, the init methods (see \ref nfsft_init , \ref
 * nfsft_init_advanced , and \ref nfsft_init_guru) will allocate memory and the
 * method \ref nfsft_finalize will free the array \c f_hat for you. Otherwise,
 * you have to assure by yourself that \c f_hat points to an array of
 * proper size before excuting a transform and you are responsible for freeing
 * the corresponding memory before program termination.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_MALLOC_F_HAT  (1U << 5)

/**
 * If this flag is set, the init methods (see \ref nfsft_init , \ref
 * nfsft_init_advanced , and \ref nfsft_init_guru) will allocate memory and the
 * method \ref nfsft_finalize will free the array \c f for you. Otherwise,
 * you have to assure by yourself that \c f points to an array of
 * proper size before excuting a transform and you are responsible for freeing
 * the corresponding memory before program termination.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_MALLOC_F      (1U << 6)

/**
 * If this flag is set, it is guaranteed that during an execution of
 * \ref ndsft_trafo or \ref nfsft_trafo the content of \c f_hat remains
 * unchanged.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_PRESERVE_F_HAT (1U << 7)

/**
 * If this flag is set, it is guaranteed that during an execution of
 * \ref ndsft_trafo, \ref nfsft_trafo or \ref ndsft_adjoint, \ref nfsft_adjoint
 * the content of \c x remains
 * unchanged.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_PRESERVE_X     (1U << 8)

/**
 * If this flag is set, it is guaranteed that during an execution of
 * \ref ndsft_adjoint or \ref nfsft_adjoint the content of \c f remains
 * unchanged.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_PRESERVE_F     (1U << 9)

/**
 * If this flag is set, it is explicitely allowed that during an execution of
 * \ref ndsft_trafo or \ref nfsft_trafo the content of \c f_hat may be changed.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_DESTROY_F_HAT    (1U << 10)

/**
 * If this flag is set, it is explicitely allowed that during an execution of
 * \ref ndsft_trafo, \ref nfsft_trafo or \ref ndsft_adjoint, \ref nfsft_adjoint
 * the content of \c x may be changed.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_DESTROY_X      (1U << 11)

/**
 * If this flag is set, it is explicitely allowed that during an execution of
 * \ref ndsft_adjoint or \ref nfsft_adjoint the content of \c f may be changed.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_DESTROY_F      (1U << 12)

/* Precomputation flags */

/**
 * If this flag is set, the transforms \ref ndsft_trafo and \ref ndsft_adjoint
 * do not work. Setting this flag saves some memory for precomputed data.
 *
 * \see nfsft_precompute
 * \see ndsft_trafo
 * \see ndsft_adjoint
 * \author Jens Keiner
 */
#define NFSFT_NO_DIRECT_ALGORITHM    (1U << 13)

/**
 * If this flag is set, the transforms \ref nfsft_trafo and \ref nfsft_adjoint
 * do not work. Setting this flag saves memory for precomputed data.
 *
 * \see nfsft_precompute
 * \see nfsft_trafo
 * \see nfsft_adjoint
 * \author Jens Keiner
 */
#define NFSFT_NO_FAST_ALGORITHM      (1U << 14)

/**
 * If this flag is set, the transforms \ref nfsft_adjoint and
 * \ref ndsft_adjoint set all unused entries in \c f_hat not corresponding to
 * spherical Fourier coefficients to zero.
 *
 * \author Jens Keiner
 */
#define NFSFT_ZERO_F_HAT             (1U << 16)

/* */

/**
 * This helper macro expands to the index \f$i\f$
 * corresponding to the spherical Fourier coefficient
 * \f$f_hat(k,n)\f$ for \f$0 \le k \le N\f$, \f$-k \le n \le k\f$ with
 * \f[
 *   (N+2)(N-n+1)+N+k+1
 * \f]
 */
#define NFSFT_INDEX(k,n,plan)        ((2*(plan)->N+2)*((plan)->N-n+1)+(plan)->N+k+1)

/**
 * This helper macro expands to the logical size of a spherical Fourier coefficients
 * array for a bandwidth N.
 */
#define NFSFT_F_HAT_SIZE(N)          ((2*N+2)*(2*N+2))

/** Structure for a NFSFT transform plan */
typedef struct
{
  /** Inherited public members */
  MACRO_MV_PLAN(fftw_complex)

  /* Public members */
  int N;                              /**< the bandwidth \f$N\f$              */
  double *x;                          /**< the nodes \f$\mathbf{x}(m) =       *
                                           \left(x_1,x_2\right) \in           *
                                           [-\frac{1}{2},\frac{1}{2}) \times  *
                                           [0,\frac{1}{2}]\f$ for             *
                                           \f$m=0,\ldots,M-1\f$,\f$M \in      *
                                           \mathbb{N},\f$                     */

  /* Private members */
  /*int NPT;*/                        /**< the next greater power of two with *
                                           respect to \f$N\f$                 */
  int t;                              /**< the logarithm of NPT with           *
                                           respect to the basis 2             */
  unsigned int flags;                 /**< the planner flags                  */
  nfft_plan plan_nfft;                /**< the internal NFFT plan             */
  fftw_complex *f_hat_intern;              /**< Internally used pointer to         *
                                           spherical Fourier coefficients     */
} nfsft_plan;

/**
 * Creates a transform plan.
 *
 * \arg plan a pointer to a \ref nfsft_plan structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 *
 * \author Jens Keiner
 */
void nfsft_init(nfsft_plan *plan, int N, int M);

/**
 * Creates a transform plan.
 *
 * \arg plan a pointer to a \verbatim nfsft_plan \endverbatim structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 * \arg nfsft_flags the NFSFT flags
 *
 * \author Jens Keiner
 */
void nfsft_init_advanced(nfsft_plan* plan, int N, int M, unsigned int
                         nfsft_flags);

/**
 * Creates a transform plan.
 *
 * \arg plan a pointer to a \verbatim nfsft_plan \endverbatim structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 * \arg nfsft_flags the NFSFT flags
 * \arg nfft_cutoff the NFFT cutoff parameter
 *
 * \author Jens Keiner
 */
void nfsft_init_guru(nfsft_plan *plan, int N, int M, unsigned int nfsft_flags,
    unsigned int nfft_flags, int nfft_cutoff);

/**
 * Performes precomputation up to the next power of two with respect to a given
 * bandwidth \f$N \in \mathbb{N}_2\f$. The threshold parameter \f$\kappa \in
 * \mathbb{R}^{+}\f$ determines the number of stabilization steps computed in
 * the discrete polynomial transform and thereby its accuracy.
 *
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg threshold the threshold \f$\kappa \in \mathbb{R}^{+}\f$
 * \arg nfsft_precomputation_flags the NFSFT precomputation flags
 * \arg fpt_precomputation_flags the FPT precomputation flags
 *
 * \author Jens Keiner
 */
void nfsft_precompute(int N, double kappa, unsigned int nfsft_flags,
  unsigned int fpt_flags);

/**
 * Forgets all precomputed data.
 *
 * \author Jens Keiner
 */
void nfsft_forget(void);

/**
 * Executes a direct NDSFT, i.e. computes for \f$m = 0,\ldots,M-1\f$
 * \f[
 *   f(m) = \sum_{k=0}^N \sum_{n=-k}^k \hat{f}(k,n) Y_k^n\left(2\pi x_1(m),
 *   2\pi x_2(m)\right).
 * \f]
 *
 * \arg plan the plan
 *
 * \author Jens Keiner
 */
void ndsft_trafo(nfsft_plan* plan);

/**
 * Executes a direct adjoint NDSFT, i.e. computes for \f$k=0,\ldots,N;
 * n=-k,\ldots,k\f$
 * \f[
 *   \hat{f}(k,n) = \sum_{m = 0}^{M-1} f(m) Y_k^n\left(2\pi x_1(m),
 *   2\pi x_2(m)\right).
 * \f]
 *
 * \arg plan the plan
 *
 * \author Jens Keiner
 */
void ndsft_adjoint(nfsft_plan* plan);

/**
 * Executes a NFSFT, i.e. computes for \f$m = 0,\ldots,M-1\f$
 * \f[
 *   f(m) = \sum_{k=0}^N \sum_{n=-k}^k \hat{f}(k,n) Y_k^n\left(2\pi x_1(m),
 *   2\pi x_2(m)\right).
 * \f]
 *
 * \arg plan the plan
 *
 * \author Jens Keiner
 */
void nfsft_trafo(nfsft_plan* plan);

/**
 * Executes an adjoint NFSFT, i.e. computes for \f$k=0,\ldots,N;
 * n=-k,\ldots,k\f$
 * \f[
 *   \hat{f}(k,n) = \sum_{m = 0}^{M-1} f(m) Y_k^n\left(2\pi x_1(m),
 *   2\pi x_2(m)\right).
 * \f]
 *
 * \arg plan the plan
 *
 * \author Jens Keiner
 */
void nfsft_adjoint(nfsft_plan* plan);

/**
 * Destroys a plan.
 *
 * \arg plan the plan to be destroyed
 *
 * \author Jens Keiner
 */
void nfsft_finalize(nfsft_plan* plan);

void nfsft_precompute_x(nfsft_plan *plan);

#endif

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/**
 * @defgroup fpt FPT - Fast polynomial transform
 * @{
 *
 * This module implements fast polynomial transforms. In the following, we
 * abbreviate the term "fast polynomial transforms" by FPT.
 */

#ifdef HAVE_FPT

/* Flags for fpt_init() */
#define FPT_NO_FAST_ALGORITHM (1U << 2) /**< If set, TODO complete comment.   */
#define FPT_NO_DIRECT_ALGORITHM (1U << 3) /**< If set, TODO complete comment.   */
#define FPT_NO_STABILIZATION  (1U << 0) /**< If set, no stabilization will be
                                             used.                            */

#define FPT_PERSISTENT_DATA   (1U << 4) /**< If set, TODO complete comment.   */

/* Flags for fpt_trafo(), dpt_transposed(), fpt_trafo(), fpt_transposed() */
#define FPT_FUNCTION_VALUES   (1U << 5) /**< If set, the output are function
                                             values at Chebyshev nodes rather
                                             than Chebyshev coefficients.     */
#define FPT_AL_SYMMETRY       (1U << 6) /**< TODO Don't use this flag!        */

/* Data structures */
typedef struct fpt_set_s_ *fpt_set;     /**< A set of precomputed data for a
                                             set of DPT transforms of equal
                                             maximum length.                  */

/**
 * Initializes a set of precomputed data for DPT transforms of equal length.
 *
 * \arg M The maximum DPT transform index \f$M \in \mathbb{N}_0\f$. The
 *        individual transforms are addressed by and index number \f$m \in
 *        \mathbb{N}_0\f$ with range \f$m = 0,\ldots,M\f$. The total number
 *        of transforms is therefore \f$M+1\f$.
 * \arg t The exponent \f$t \in \mathbb{N}, t \ge 2\f$ of the transform length
 *        \f$N = 2^t \in \mathbb{N}, N \ge 4\f$
 * \arg flags A bitwise combination of the flags FPT_NO_STABILIZATION,
 *
 * \author Jens Keiner
 */
fpt_set fpt_init(const int M, const int t, const unsigned int flags);

/**
 * Computes the data required for a single DPT transform.
 *
 * \arg set The set of DPT transform data where the computed data will be stored.
 * \arg m The transform index \f$m \in \mathbb{N}, 0 \le m \le M\f$.
 * \arg alpha The three-term recurrence coefficients \f$\alpha_k \in
 *      \mathbb{R}\f$ for \f$k=0,\ldots,N\f$ such that \verbatim alpha[k]
 *      \endverbatim \f$=\alpha_k\f$.
 * \arg beta The three-term recurrence coefficients \f$\beta_k \in \mathbb{R}\f$
 *            for \f$k=0,\ldots,N\f$ such that \verbatim beta[k] \endverbatim
 *            \f$=\beta_k\f$.
 * \arg gamma The three-term recurrence coefficients \f$\gamma_k \in
 *            \mathbb{R}\f$ for \f$k=0,\ldots,N\f$ such that \verbatim gamma[k]
 *            \endverbatim \f$=\gamma_k\f$.
 * \arg k_start The index \f$k_{\text{start}} \in \mathbb{N}_0,
 *              0 \le k_{\text{start}} \le N\f$
 * \arg threshold The treshold \f$\kappa \in \mathbb{R}, \kappa > 0\f$.
 *
 * \author Jens Keiner
 */
void fpt_precompute(fpt_set set, const int m, double *alpha, double *beta,
  double *gam, int k_start, const double threshold);

/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void dpt_trafo(fpt_set set, const int m, const fftw_complex *x, fftw_complex *y,
  const int k_end, const unsigned int flags);

/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void fpt_trafo(fpt_set set, const int m, const fftw_complex *x, fftw_complex *y,
  const int k_end, const unsigned int flags);

/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void dpt_transposed(fpt_set set, const int m, fftw_complex *x,
  fftw_complex *y, const int k_end, const unsigned int flags);

/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void fpt_transposed(fpt_set set, const int m, fftw_complex *x,
  fftw_complex *y, const int k_end, const unsigned int flags);

void fpt_finalize(fpt_set set);

#endif

/** @}
 */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/**
 * @defgroup nfsoft NFSOFT - Nonequispaced fast SO(3) Fourier transform
 * @{
 *
 * This module implements nonuniform fast SO(3) Fourier transforms. In the
 * following, we abbreviate the term "nonuniform fast SO(3) Fourier
 * transform" by NFSOFT.
 *
 */

#ifdef HAVE_NFSOFT

/* Planner flags */
/**
 * By default, all computations are performed with respect to the
 * unnormalized basis functions
 * \f[
 *   D_{mn}^l(\alpha,\beta,\gamma) = d^{mn}_{l}(\cos\beta)
 *   \mathrm{e}^{-\mathrm{i} m \alpha}\mathrm{e}^{-\mathrm{i} n \gamma}.
 * \f]
 * If this flag is set, all computations are carried out using the \f$L_2\f$-
 * normalized basis functions
 * \f[
 *  \tilde D_{mn}^l(\alpha,\beta,\gamma) = \sqrt{\frac{2l+1}{8\pi^2}}d^{mn}_{l}(\cos\beta)
 *   \mathrm{e}^{-\mathrm{i} m \alpha}\mathrm{e}^{-\mathrm{i} n \gamma}
 * \f]
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_NORMALIZED    (1U << 0)

/**
 * If this flag is set, the fast NFSOFT algorithms (see \ref nfsoft_trafo,
 * \ref nfsoft_adjoint) will use internally the exact but usually slower direct
 * NDFT algorithm in favor of fast but approximative NFFT algorithm.
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_USE_NDFT      (1U << 1)

/**
 * If this flag is set, the fast NFSOFT algorithms (see \ref nfsoft_trafo,
 * \ref nfsoft_adjoint) will use internally the usually slower direct
 * DPT algorithm in favor of the fast FPT algorithm.
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_USE_DPT       (1U << 2)

/**
 * If this flag is set, the init methods (see \ref nfsoft_init ,
 * \ref nfsoft_init_advanced , and \ref nfsoft_init_guru) will allocate memory and the
 * method \ref nfsoft_finalize will free the array \c x for you. Otherwise,
 * you have to assure by yourself that \c x points to an array of
 * proper size before excuting a transform and you are responsible for freeing
 * the corresponding memory before program termination.
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_MALLOC_X      (1U << 3)

/**
 * If this flag is set, the Wigner-D functions will be normed
 * such that they satisfy the representation property of
 * the spherical harmonics as defined in the NFFT software package
 *
 * \author Antje Vollrath
 */
#define NFSOFT_REPRESENT      (1U << 4)


/**
 * If this flag is set, the init methods (see \ref nfsoft_init ,
 * \ref nfsoft_init_advanced , and \ref nfsoft_init_guru) will allocate memory and the
 * method \ref nfsoft_finalize will free the array \c f_hat for you. Otherwise,
 * you have to assure by yourself that \c f_hat points to an array of
 * proper size before excuting a transform and you are responsible for freeing
 * the corresponding memory before program termination.
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_MALLOC_F_HAT  (1U << 5)

/**
 * If this flag is set, the init methods (see \ref nfsoft_init ,
 * \ref nfsoft_init_advanced , and \ref nfsoft_init_guru) will allocate memory and the
 * method \ref nfsoft_finalize will free the array \c f for you. Otherwise,
 * you have to assure by yourself that \c f points to an array of
 * proper size before excuting a transform and you are responsible for freeing
 * the corresponding memory before program termination.
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_MALLOC_F      (1U << 6)

/**
 * If this flag is set, it is guaranteed that during an execution of
 * \ref nfsoft_trafo the content of \c f_hat remains
 * unchanged.
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_PRESERVE_F_HAT (1U << 7)

/**
 * If this flag is set, it is guaranteed that during an execution of
 * \ref nfsoft_trafo or \ref nfsoft_adjoint
 * the content of \c x remains
 * unchanged.
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_PRESERVE_X     (1U << 8)

/**
 * If this flag is set, it is guaranteed that during an execution of
 * \ref ndsoft_adjoint or \ref nfsoft_adjoint the content of \c f remains
 * unchanged.
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_PRESERVE_F     (1U << 9)

/**
 * If this flag is set, it is explicitely allowed that during an execution of
 * \ref nfsoft_trafo the content of \c f_hat may be changed.
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_DESTROY_F_HAT    (1U << 10)

/**
 * If this flag is set, it is explicitely allowed that during an execution of
 * \ref nfsoft_trafo or \ref nfsoft_adjoint
 * the content of \c x may be changed.
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_DESTROY_X      (1U << 11)

/**
 * If this flag is set, it is explicitely allowed that during an execution of
 * \ref ndsoft_adjoint or \ref nfsoft_adjoint the content of \c f may be changed.
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_DESTROY_F      (1U << 12)

/**
 * If this flag is set, the fast NFSOFT algorithms (see \ref nfsoft_trafo,
 * \ref nfsoft_adjoint) will use internally the FPT algorithm without the
 * stabilization scheme and thus making bigger errors for higher
 * bandwidth but becoming significantly faster
 *
 * \author Antje Vollrath
 */
#define NFSOFT_NO_STABILIZATION      (1U << 13)

/**
 * If this flag is set, the fast NFSOFT algorithms (see \ref nfsoft_trafo,
 * \ref nfsoft_adjoint) will decide whether to use the DPT or
 * FPT algorithm depending on which is faster for the chosen orders.
 *
 * not yet included in the checked-in version
 *
 * \author Antje Vollrath
 */
#define NFSOFT_CHOOSE_DPT            (1U << 14)

/**
 * If this flag is set, the fast NFSOFT algorithms (see \ref nfsoft_trafo,
 * \ref nfsoft_adjoint) becomes a SOFT, i.e., we use equispaced nodes.
 * The FFTW will be used instead of the NFFT.-->not included yet
 *
 * \see nfsoft_init
 * \see nfsoft_init_advanced
 * \see nfsoft_init_guru
 * \author Antje Vollrath
 */
#define NFSOFT_SOFT                  (1U << 15)


/**
 * If this flag is set, the transform \ref nfsoft_adjoint
 * sets all unused entries in \c f_hat not corresponding to
 * SO(3) Fourier coefficients to zero.
 *
 * \author Antje Vollrath
 */
#define NFSOFT_ZERO_F_HAT             (1U << 16)


/* Helper macros*/
/**
 * These macro expands to the index \f$i\f$
 * corresponding to the SO(3) Fourier coefficient
 * \f$f_hat^{mn}_l\f$ for \f$l=0,...,B\f$, \f$m,n =-l,...,l\f$ with
 */
#define NFSOFT_INDEX(m,n,l,B)        (((l)+((B)+1))+(2*(B)+2)*(((n)+((B)+1))+(2*(B)+2)*((m)+((B)+1))))
#define NFSOFT_INDEX_TWO(m,n,l,B) ((B+1)*(B+1)+(B+1)*(B+1)*(m+B)-((m-1)*m*(2*m-1)+(B+1)*(B+2)*(2*B+3))/6)+(posN(n,m,B))+(l-MAX(ABS(m),ABS(n)))
int posN(int n,int m, int B);

/**
 * This macro expands to the logical size of a SO(3) Fourier coefficients
 * array for a bandwidth B.
 */
#define NFSOFT_F_HAT_SIZE(B)          (((B)+1)*(4*((B)+1)*((B)+1)-1)/3)

/** Structure for a NFSOFT transform plan */
typedef struct nfsoft_plan_
{
  /** Inherited public members */
  MACRO_MV_PLAN(fftw_complex)

  double *x;                           /**< the  input nodes                    */
  /**some auxillary memory*/
  fftw_complex *wig_coeffs;		       /**< contains a set of SO(3) Fourier coefficients*
                                    for fixed orders m and n*/
  fftw_complex *cheby;		       /**< contains a set of Chebychev coefficients for*
                                    fixed orders m and n*/
  fftw_complex *aux;			       /**< used when converting Chebychev to Fourier*
                                    coeffcients*/

  /** Private members */
  int t;                               /**< the logaritm of NPT with          *
                                          respect to the basis 2              */
  unsigned int flags;                  /**< the planner flags                 */
  nfft_plan nfft_plan;                /**< the internal NFFT plan             */
  fftw_plan fftw_plan;                /**< the optional internal FFTW plan    */

  int fpt_kappa;       /**a parameter controlling the accuracy of the FPT*/

} nfsoft_plan;


/**
 * Creates a transform plan.
 *
 * \arg plan a pointer to a \ref nfsoft_plan structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 *
 * \author Antje Vollrath
 */

/** Functions for NFSOFT plans*/
/**
 * Computes a single transposed FPT transform.
 *
 * \arg coeffs the Chebychev coefficientss that should be transformed
 * \arg l the polynomial degree
 * \arg k the first order
 * \arg m the second order
 * \arg kappa the parameter controlling the accuracy
 * \arg nfsoft_flags
 */

fpt_set SO3_fpt_init(int l, unsigned int flags,int kappa);
void SO3_fpt(fftw_complex *coeffs, fpt_set set, int l, int k, int m, unsigned int nfsoft_flags);
void SO3_fpt_transposed(fftw_complex *coeffs,fpt_set set,int l, int k, int m,unsigned int nfsoft_flags);


/**
 * Creates a NFSOFT transform plan.
 *
 * \arg plan a pointer to a \ref nfsoft_plan structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 *
 * \author Antje Vollrath
 */
void nfsoft_init(nfsoft_plan *plan, int N, int M);
/**
 * Creates a NFSOFT transform plan.
 *
 * \arg plan a pointer to a \ref nfsoft_plan structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 * \arg nfsoft_flags the NFSOFT flags
 *
 * \author Antje Vollrath
 */
void nfsoft_init_advanced(nfsoft_plan *plan, int N, int M,unsigned int nfsoft_flags);
/**
 * Creates a  NFSOFT transform plan.
 *
 * \arg plan a pointer to a \ref nfsoft_plan structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 * \arg nfsoft_flags the NFSFT flags
 * \arg nfft_flags the NFFT flags
 * \arg fpt_kappa a parameter contolling the accuracy of the FPT
 * \arg nfft_cutoff the NFFT cutoff parameter
 *
 * \author Antje Vollrath
 */
void nfsoft_init_guru(nfsoft_plan *plan, int N, int M,unsigned int nfsoft_flags,unsigned int nfft_flags,int nfft_cutoff,int fpt_kappa);

/**
 * Executes a NFSOFT, i.e. computes for \f$m = 0,\ldots,M-1\f$
 * \f[
 *   f(g_m) = \sum_{l=0}^B \sum_{m=-l}^l \sum_{n=-l}^l \hat{f}^{mn}_l
 *            D_l^{mn}\left( \alpha_m,\beta_m,\gamma_m\right).
 * \f]
 *
 * \arg plan_nfsoft the plan
 *
 * \author Antje Vollrath
 */
void nfsoft_trafo(nfsoft_plan *plan_nfsoft);
/**
 * Executes an adjoint NFSOFT, i.e. computes for \f$l=0,\ldots,B;
 * m,n=-l,\ldots,l\f$
 * \f[
 *   \hat{f}^{mn}_l = \sum_{m = 0}^{M-1} f(g_m)
 *                    D_l^{mn}\left( \alpha_m,\beta_m,\gamma_m\right)
 * \f]
 *
 * \arg plan_nfsoft the plan
 *
 * \author Antje Vollrath
 */
void nfsoft_adjoint(nfsoft_plan *plan_nfsoft);
/**
 * Destroys a plan.
 *
 * \arg plan the plan to be destroyed
 *
 * \author Antje Vollrath
 */
void nfsoft_finalize(nfsoft_plan *plan);


#endif

/** @}
 */


/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/** @defgroup solver Solver - Inverse transforms
 * @{
 */

/* Planner flags, i.e. constant symbols for methods */

/**
 * If this flag is set, the Landweber (Richardson) iteration is used to compute
 * an inverse transform.
 *
 * \author Stefan Kunis
 */
#define LANDWEBER             (1U<< 0)

/**
 * If this flag is set, the method of steepest descent (gradient) is used to
 * compute an inverse transform.
 *
 * \author Stefan Kunis
 */
#define STEEPEST_DESCENT      (1U<< 1)

/**
 * If this flag is set, the conjugate gradient method for the normal equation
 * of first kind is used to compute an inverse transform.
 * Each iterate minimises the residual in the current Krylov subspace.
 *
 * \author Stefan Kunis
 */
#define CGNR                  (1U<< 2)

/**
 * If this flag is set, the conjugate gradient method for the normal equation
 * of second kind is used to compute an inverse transform.
 * Each iterate minimises the error in the current Krylov subspace.
 *
 * \author Stefan Kunis
 */
#define CGNE                  (1U<< 3)

/**
 * If this flag is set, the Landweber iteration updates the member
 * dot_r_iter.
 *
 * \author Stefan Kunis
 */
#define NORMS_FOR_LANDWEBER   (1U<< 4)

/**
 * If this flag is set, the samples are weighted, eg to cope with varying
 * sampling density.
 *
 * \author Stefan Kunis
 */
#define PRECOMPUTE_WEIGHT     (1U<< 5)

/**
 * If this flag is set, the Fourier coefficients are damped, eg to favour
 * fast decaying coefficients.
 *
 * \author Stefan Kunis
 */
#define PRECOMPUTE_DAMP       (1U<< 6)


typedef struct
{
  mv_plan_complex *mv;                  /**< matrix vector multiplication   */
  unsigned flags;                       /**< iteration type                 */

  double *w;                            /**< weighting factors              */
  double *w_hat;                        /**< damping factors                */

  fftw_complex *y;                      /**< right hand side, samples       */

  fftw_complex *f_hat_iter;             /**< iterative solution             */

  fftw_complex *r_iter;                 /**< iterated residual vector       */
  fftw_complex *z_hat_iter;             /**< residual of normal equation of
					     first kind                     */
  fftw_complex *p_hat_iter;             /**< search direction               */
  fftw_complex *v_iter;                 /**< residual vector update         */

  double alpha_iter;                    /**< step size for search direction */
  double beta_iter;                     /**< step size for search correction*/

  double dot_r_iter;                    /**< weighted dotproduct of r_iter  */
  double dot_r_iter_old;                /**< previous dot_r_iter            */
  double dot_z_hat_iter;                /**< weighted dotproduct of
					     z_hat_iter                     */
  double dot_z_hat_iter_old;            /**< previous dot_z_hat_iter        */
  double dot_p_hat_iter;                /**< weighted dotproduct of
					     p_hat_iter                     */
  double dot_v_iter;                    /**< weighted dotproduct of v_iter  */
} solver_plan_complex;

void solver_init_advanced_complex(solver_plan_complex* ths, mv_plan_complex *mv, unsigned flags);
void solver_init_complex(solver_plan_complex* ths, mv_plan_complex *mv);
void solver_before_loop_complex(solver_plan_complex* ths);
void solver_loop_one_step_complex(solver_plan_complex *ths);
void solver_finalize_complex(solver_plan_complex *ths);

typedef struct
{
  mv_plan_double *mv;                   /**< matrix vector multiplication   */
  unsigned flags;                       /**< iteration type                 */

  double *w;                            /**< weighting factors              */
  double *w_hat;                        /**< damping factors                */

  double *y;                            /**< right hand side, samples       */

  double *f_hat_iter;                   /**< iterative solution             */

  double *r_iter;                       /**< iterated residual vector       */
  double *z_hat_iter;                   /**< residual of normal equation of
					     first kind                     */
  double *p_hat_iter;                   /**< search direction               */
  double *v_iter;                       /**< residual vector update         */

  double alpha_iter;                    /**< step size for search direction */
  double beta_iter;                     /**< step size for search correction*/

  double dot_r_iter;                    /**< weighted dotproduct of r_iter  */
  double dot_r_iter_old;                /**< previous dot_r_iter            */
  double dot_z_hat_iter;                /**< weighted dotproduct of
					     z_hat_iter                     */
  double dot_z_hat_iter_old;            /**< previous dot_z_hat_iter        */
  double dot_p_hat_iter;                /**< weighted dotproduct of
					     p_hat_iter                     */
  double dot_v_iter;                    /**< weighted dotproduct of v_iter  */
} solver_plan_double;

void solver_init_advanced_double(solver_plan_double* ths, mv_plan_double *mv, unsigned flags);
void solver_init_double(solver_plan_double* ths, mv_plan_double *mv);
void solver_before_loop_double(solver_plan_double* ths);
void solver_loop_one_step_double(solver_plan_double *ths);
void solver_finalize_double(solver_plan_double *ths);

/** @}
 */

#endif
/* nfft3.h */
