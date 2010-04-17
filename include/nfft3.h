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

#ifndef __NFFT3_H__
#define __NFFT3_H__

/* module configuration */
#include "nfft3conf.h"

/* fftw_complex */
#include <fftw3.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#define NFFT_CONCAT(prefix, name) prefix ## name

#define NFFT_MANGLE_DOUBLE(name) NFFT_CONCAT(nfft_, name)
#define NFFT_MANGLE_FLOAT(name) NFFT_CONCAT(nfftf_, name)
#define NFFT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nfftl_, name)

/* IMPORTANT: for Windows compilers, you should add a line
 *   #define FFTW_DLL
 * here and in kernel/infft.h if you are compiling/using NFFT as a DLL, in order
 * to do the proper importing/exporting, or alternatively compile with
 * -DNFFT_DLL or the equivalent command-line flag. This is not necessary under
 * MinGW/Cygwin, where libtool does the imports/exports automatically. */
#if defined(NFFT_DLL) && (defined(_WIN32) || defined(__WIN32__))
  /* annoying Windows syntax for shared-library declarations */
#  if defined(COMPILING_NFFT) /* defined in api.h when compiling NFFT */
#    define NFFT_EXTERN extern __declspec(dllexport)
#  else /* user is calling NFFT; import symbol */
#    define NFFT_EXTERN extern __declspec(dllimport)
#  endif
#else
#  define NFFT_EXTERN extern
#endif

/* Malloc and free functions */
NFFT_EXTERN void *nfft_malloc(size_t n);
NFFT_EXTERN void nfft_free(void *p);
NFFT_EXTERN void nfft_die(const char *s);

/* Malloc and free hooks */
typedef void *(*nfft_malloc_type_function) (size_t n);
typedef void  (*nfft_free_type_function) (void *p);
typedef void  (*nfft_die_type_function) (const char *errString);
NFFT_EXTERN nfft_malloc_type_function nfft_malloc_hook;
NFFT_EXTERN nfft_free_type_function nfft_free_hook;
NFFT_EXTERN nfft_die_type_function nfft_die_hook;

/** Macros for public members inherited by all plan structures. */
#define MACRO_MV_PLAN(RC) \
  int N_total; /**< Total number of Fourier coefficients */\
  int M_total; /**< Total number of samples */\
\
  RC *f_hat; /**< Vector of Fourier coefficients, size is N_total * sizeof(RC) */\
  RC *f; /**< Vector of samples, size is M_total * sizeof(RC) */\
  void (*mv_trafo)(void*); /**< Pointer to the own transform */\
  void (*mv_adjoint)(void*); /**< Pointer to the own adjoint */

/* huge second-order macro that defines prototypes for all API functions. We
 * expand this macro for each supported precision.
 *   X: name-mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define NFFT_DEFINE_API(X,Y,R,C) \
\
typedef struct \
{ \
  MACRO_MV_PLAN(C) \
} X(mv_plan_complex); \
\
typedef struct \
{ \
  MACRO_MV_PLAN(R) \
} X(mv_plan_double); \
\
/** Structure for a NFFT plan */\
typedef struct\
{\
  /* api */\
  MACRO_MV_PLAN(C)\
\
  int d; /**< dimension aka rank */\
  int *N; /**< multi-bandwidth */\
  R *sigma; /**< oversampling-factor */\
  int *n; /**< FFTW length, equal to sigma*N, default is the power of 2 such
               that \f$2\le\sigma<4\f$ */\
  int n_total; /**< Total size of FFTW */\
  int m; /**< Cut-off parameter of the window function, default value is
               6 (KAISER_BESSEL),
               9 (SINC_POWER),
               11 (B_SPLINE),
               12 (GAUSSIAN) */\
  R *b; /**< Shape parameter of the window function */\
  int K; /**< Number of equispaced samples of the window function for \ref
              PRE_LIN_PSI */\
\
  unsigned nfft_flags; /**< Flags for precomputation, (de)allocation, and FFTW
                            usage, default setting is
                            PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT |
                            MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE */\
\
  unsigned fftw_flags; /**< Flags for the FFTW, default is
                            FFTW_ESTIMATE | FFTW_DESTROY_INPUT */\
\
  R *x; /**< Nodes in time/spatial domain, size is \f$dM\f$ doubles */\
\
  double MEASURE_TIME_t[3]; /**< Measured time for each step if MEASURE_TIME is
                                 set (always a double!)*/\
\
  /* internal*/\
  Y(plan)  my_fftw_plan1; /**< Forward FFTW plan */\
  Y(plan)  my_fftw_plan2; /**< Backward FFTW plan */\
\
  R **c_phi_inv; /**< Precomputed data for the diagonal matrix \f$D\f$,
                      size is \f$N_0+\hdots+N_{d-1}\f$ doubles*/\
  R *psi; /**< Precomputed data for the sparse matrix \f$B\f$, size depends
                    on precomputation scheme */\
  int *psi_index_g; /**< Indices in source/target vector for \ref PRE_FULL_PSI */\
  int *psi_index_f; /**< Indices in source/target vector for \ref PRE_FULL_PSI */\
\
  C *g; /**< Oversampled vector of samples, size is \ref n_total
                        double complex */\
  C *g_hat; /**< Zero-padded vector of Fourier coefficients, size is
                            \ref n_total fftw_complex */\
  C *g1; /**< Input of fftw */\
  C *g2; /**< Output of fftw */\
\
  R *spline_coeffs; /**< Input for de Boor algorithm if B_SPLINE or
                              SINC_POWER is defined */\
} X(plan); \
\
NFFT_EXTERN void X(direct_trafo)(X(plan) *ths);\
NFFT_EXTERN void X(direct_adjoint)(X(plan) *ths);\
NFFT_EXTERN void X(trafo)(X(plan) *ths);\
NFFT_EXTERN void X(trafo_1d)(X(plan) *ths);\
NFFT_EXTERN void X(trafo_2d)(X(plan) *ths);\
NFFT_EXTERN void X(trafo_3d)(X(plan) *ths);\
NFFT_EXTERN void X(adjoint)(X(plan) *ths);\
NFFT_EXTERN void X(adjoint_1d)(X(plan) *ths);\
NFFT_EXTERN void X(adjoint_2d)(X(plan) *ths);\
NFFT_EXTERN void X(adjoint_3d)(X(plan) *ths);\
NFFT_EXTERN void X(init_1d)(X(plan) *ths, int N1, int M);\
NFFT_EXTERN void X(init_2d)(X(plan) *ths, int N1, int N2, int M);\
NFFT_EXTERN void X(init_3d)(X(plan) *ths, int N1, int N2, int N3, int M);\
NFFT_EXTERN void X(init)(X(plan) *ths, int d, int *N, int M);\
NFFT_EXTERN void X(init_advanced)(X(plan) *ths, int d, int *N, int M, \
  unsigned nfft_flags_on, unsigned nfft_flags_off);\
NFFT_EXTERN void X(init_guru)(X(plan) *ths, int d, int *N, int M, int *n, \
  int m, unsigned nfft_flags, unsigned fftw_flags);\
NFFT_EXTERN void X(precompute_one_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_full_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_lin_psi)(X(plan) *ths);\
NFFT_EXTERN void X(check)(X(plan) *ths);\
NFFT_EXTERN void X(finalize)(X(plan) *ths);

NFFT_DEFINE_API(NFFT_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,float,fftwf_complex)
NFFT_DEFINE_API(NFFT_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,double,fftw_complex)
NFFT_DEFINE_API(NFFT_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

/* init flags */

#define PRE_PHI_HUT      (1U<< 0)
#define FG_PSI           (1U<< 1)
#define PRE_LIN_PSI      (1U<< 2)
#define PRE_FG_PSI       (1U<< 3)
#define PRE_PSI          (1U<< 4)
#define PRE_FULL_PSI     (1U<< 5)
#define MALLOC_X         (1U<< 6)
#define MALLOC_F_HAT     (1U<< 7)
#define MALLOC_F         (1U<< 8)
#define FFT_OUT_OF_PLACE (1U<< 9)
#define FFTW_INIT        (1U<< 10)
#define PRE_ONE_PSI (PRE_LIN_PSI| PRE_FG_PSI| PRE_PSI| PRE_FULL_PSI)

/*###########################################################################*/

#define NFCT_MANGLE_DOUBLE(name) NFFT_CONCAT(nfct_, name)
#define NFCT_MANGLE_FLOAT(name) NFFT_CONCAT(nfctf_, name)
#define NFCT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nfctl_, name)

#define NFCT_DEFINE_API(X,Y,R,C) \
/** Structure for a transform plan */\
typedef struct\
{\
  /* api */\
  MACRO_MV_PLAN(R)\
\
  int d;                                /**< dimension, rank                  */\
  int *N;                               /**< cut-off-frequencies (kernel)     */\
  int *n;                               /**< length of dct-i                  */\
  R *sigma;                             /**< oversampling-factor              */\
  int m;                                /**< cut-off parameter in time-domain */\
\
  R nfct_full_psi_eps;\
  R *b;                            /**< shape parameters                 */\
\
  unsigned nfct_flags;                  /**< flags for precomputation, malloc */\
  unsigned fftw_flags;                  /**< flags for the fftw               */\
\
  R *x;                            /**< nodes (in time/spatial domain)   */\
\
  double MEASURE_TIME_t[3];             /**< measured time for each step      */\
\
  /** internal */\
  Y(plan)  my_fftw_r2r_plan;          /**< fftw_plan                        */\
  Y(r2r_kind) *r2r_kind;              /**< r2r transform type (dct-i)       */\
\
  R **c_phi_inv;                   /**< precomputed data, matrix D       */\
  R *psi;                          /**< precomputed data, matrix B       */\
  int size_psi;                         /**< only for thin B                  */\
  int *psi_index_g;                     /**< only for thin B                  */\
  int *psi_index_f;                     /**< only for thin B                  */\
\
  R *g;\
  R *g_hat;\
  R *g1;                           /**< input of fftw                    */\
  R *g2;                           /**< output of fftw                   */\
\
  R *spline_coeffs;                /**< input for de Boor algorithm, if\
                                             B_SPLINE or SINC_2m is defined   */\
} X(plan);\
\
NFFT_EXTERN void X(init_1d)(X(plan) *ths_plan, int N0, int M_total); \
NFFT_EXTERN void X(init_2d)(X(plan) *ths_plan, int N0, int N1, int M_total); \
NFFT_EXTERN void X(init_3d)(X(plan) *ths_plan, int N0, int N1, int N2, int M_total); \
NFFT_EXTERN void X(init)(X(plan) *ths_plan, int d, int *N, int M_total); \
NFFT_EXTERN void X(init_guru)(X(plan) *ths_plan, int d, int *N, int M_total, int *n, \
  int m, unsigned nfct_flags, unsigned fftw_flags); \
NFFT_EXTERN void X(precompute_psi)(X(plan) *ths_plan); \
NFFT_EXTERN void X(trafo)(X(plan) *ths_plan); \
NFFT_EXTERN void X(direct_trafo)(X(plan) *ths_plan); \
NFFT_EXTERN void X(adjoint)(X(plan) *ths_plan); \
NFFT_EXTERN void X(direct_adjoint)(X(plan) *ths_plan); \
NFFT_EXTERN void X(finalize)(X(plan) *ths_plan); \
NFFT_EXTERN R X(phi_hut)(X(plan) *ths_plan, int k, int d); \
NFFT_EXTERN R X(phi)(X(plan) *ths_plan, R x, int d); \
NFFT_EXTERN int X(fftw_2N)(int n); \
NFFT_EXTERN int X(fftw_2N_rev)(int n);

#if defined(HAVE_NFCT)
NFCT_DEFINE_API(NFCT_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,float,fftwf_complex)
NFCT_DEFINE_API(NFCT_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,double,fftw_complex)
NFCT_DEFINE_API(NFCT_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,long double,fftwl_complex)
#endif

/*###########################################################################*/

#define NFST_MANGLE_DOUBLE(name) NFFT_CONCAT(nfst_, name)
#define NFST_MANGLE_FLOAT(name) NFFT_CONCAT(nfstf_, name)
#define NFST_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nfstl_, name)

#define NFST_DEFINE_API(X,Y,R,C) \
/** Structure for a transform plan */\
typedef struct\
{\
  /* api */\
  MACRO_MV_PLAN(R)\
\
  int d;                                /**< dimension, rank                  */\
  int *N;                               /**< bandwidth                        */\
  int *n;                               /**< length of dst-1                  */\
  R *sigma;                        /**< oversampling-factor              */\
  int m;                                /**< cut-off parameter in time-domain */\
\
  R nfst_full_psi_eps;\
  R *b;                            /**< shape parameters                 */\
\
  unsigned nfst_flags;                  /**< flags for precomputation, malloc */\
  unsigned fftw_flags;                  /**< flags for the fftw               */\
\
  R *x;                            /**< nodes (in time/spatial domain)   */\
\
  double MEASURE_TIME_t[3];             /**< measured time for each step     */\
\
  /** internal */\
  Y(plan)  my_fftw_r2r_plan;         /**< fftw_plan forward                */\
  Y(r2r_kind) *r2r_kind;              /**< r2r transform type (dct-i)       */\
\
  R **c_phi_inv;                   /**< precomputed data, matrix D       */\
  R *psi;                          /**< precomputed data, matrix B       */\
  int size_psi;                         /**< only for thin B                  */\
  int *psi_index_g;                     /**< only for thin B                  */\
  int *psi_index_f;                     /**< only for thin B                  */\
\
  R *g;\
  R *g_hat;\
  R *g1;                           /**< input of fftw                    */\
  R *g2;                           /**< output of fftw                   */\
\
  R *spline_coeffs;                /**< input for de Boor algorithm, if\
                                             B_SPLINE or SINC_2m is defined   */\
} X(plan);\
\
NFFT_EXTERN void X(init_1d)(X(plan) *ths_plan, int N0, int M_total); \
NFFT_EXTERN void X(init_2d)(X(plan) *ths_plan, int N0, int N1, int M_total); \
NFFT_EXTERN void X(init_3d)(X(plan) *ths_plan, int N0, int N1, int N2, int M_total); \
NFFT_EXTERN void X(init)(X(plan) *ths_plan, int d, int *N, int M_total); \
NFFT_EXTERN void X(init_m)(X(plan) *ths_plan, int d, int *N, int M_total, int m);\
NFFT_EXTERN void X(init_guru)(X(plan) *ths_plan, int d, int *N, int M_total, int *n, \
  int m, unsigned nfst_flags, unsigned fftw_flags); \
NFFT_EXTERN void X(precompute_psi)(X(plan) *ths_plan); \
NFFT_EXTERN void X(trafo)(X(plan) *ths_plan); \
NFFT_EXTERN void X(direct_trafo)(X(plan) *ths_plan); \
NFFT_EXTERN void X(adjoint)(X(plan) *ths_plan); \
NFFT_EXTERN void X(direct_adjoint)(X(plan) *ths_plan); \
NFFT_EXTERN void X(finalize)(X(plan) *ths_plan); \
NFFT_EXTERN void X(full_psi)(X(plan) *ths_plan, R eps); \
NFFT_EXTERN R X(phi_hut)(X(plan) *ths_plan, int k, int d); \
NFFT_EXTERN R X(phi)(X(plan) *ths_plan, R x, int d); \
NFFT_EXTERN int X(fftw_2N)(int n); \
NFFT_EXTERN int X(fftw_2N_rev)(int n);

#ifdef HAVE_NFST
NFST_DEFINE_API(NFST_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,float,fftwf_complex)
NFST_DEFINE_API(NFST_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,double,fftw_complex)
NFST_DEFINE_API(NFST_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,long double,fftwl_complex)
#endif

/** @}
 */

/*###########################################################################*/

#define NNFFT_MANGLE_DOUBLE(name) NFFT_CONCAT(nnfft_, name)
#define NNFFT_MANGLE_FLOAT(name) NFFT_CONCAT(nnfftf_, name)
#define NNFFT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nnfftl_, name)

#define NNFFT_DEFINE_API(X,Y,Z,R,C) \
/** Structure for a transform plan */\
typedef struct\
{\
  /* api */\
  MACRO_MV_PLAN(C)\
\
  int d;                                /**< dimension, rank                 */\
  R *sigma;                        /**< oversampling-factor             */\
  R *a;                            /**< 1 + 2*m/N1                      */\
  int *N;                               /**< cut-off-frequencies             */\
  int *N1;                              /**< sigma*N                         */\
  int *aN1;                             /**< sigma*a*N                       */\
  int m;                                /**< cut-off parameter in time-domain*/\
  R *b;                            /**< shape parameters                */\
  int K;                                /**< number of precomp. uniform psi  */\
\
  int aN1_total;                        /**< aN1_total=aN1[0]* ... *aN1[d-1] */\
\
  Z(plan) *direct_plan;               /**< plan for the nfft               */\
  unsigned nnfft_flags;                 /**< flags for precomputation, malloc*/\
  int *n;                               /**< n=N1, for the window function   */\
\
  R *x;                            /**< nodes (in time/spatial domain)  */\
  R *v;                            /**< nodes (in fourier domain)       */\
\
  R *c_phi_inv;                    /**< precomputed data, matrix D      */\
  R *psi;                          /**< precomputed data, matrix B      */\
  int size_psi;                         /**< only for thin B                 */\
  int *psi_index_g;                     /**< only for thin B                 */\
  int *psi_index_f;                     /**< only for thin B                 */\
  C *F;\
\
  R *spline_coeffs;                /**< input for de Boor algorithm, if\
                                             B_SPLINE or SINC_2m is defined  */\
} X(plan);\
\
NFFT_EXTERN void X(init)(X(plan) *ths_plan, int d, int N_total, int M_total, int *N); \
NFFT_EXTERN void X(init_guru)(X(plan) *ths_plan, int d, int N_total, int M_total, \
  int *N, int *N1, int m, unsigned nnfft_flags); \
NFFT_EXTERN void X(direct_trafo)(X(plan) *ths_plan); \
NFFT_EXTERN void X(direct_adjoint)(X(plan) *ths_plan); \
NFFT_EXTERN void X(trafo)(X(plan) *ths_plan); \
NFFT_EXTERN void X(adjoint)(X(plan) *ths_plan); \
NFFT_EXTERN void X(precompute_lin_psi)(X(plan) *ths_plan); \
NFFT_EXTERN void X(precompute_psi)(X(plan) *ths_plan); \
NFFT_EXTERN void X(precompute_full_psi)(X(plan) *ths_plan); \
NFFT_EXTERN void X(precompute_phi_hut)(X(plan) *ths_plan); \
NFFT_EXTERN void X(finalize)(X(plan) *ths_plan);

#ifdef HAVE_NNFFT
NNFFT_DEFINE_API(NNFFT_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,float,fftwf_complex)
NNFFT_DEFINE_API(NNFFT_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,double,fftw_complex)
NNFFT_DEFINE_API(NNFFT_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)
#endif

#define MALLOC_V         (1U<< 11)

/** @}
 */

/*###########################################################################*/

#define NSFFT_MANGLE_DOUBLE(name) NFFT_CONCAT(nsfft_, name)
#define NSFFT_MANGLE_FLOAT(name) NFFT_CONCAT(nsfftf_, name)
#define NSFFT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nsfftl_, name)

#define NSFFT_DEFINE_API(X,Y,Z,R,C) \
/** Structure for a NFFT plan */\
typedef struct\
{\
  MACRO_MV_PLAN(C)\
\
  int d;                                /**< dimension, rank; d=2,3          */\
  int J;                                /**< problem size, i.e.,\
                                             d=2: N_total=(J+4) 2^(J+1)\
                                             d=3: N_total=2^J 6(2^((J+1)/2+1)\
                                                          -1)+2^(3(J/2+1))   */\
  int sigma;                            /**< oversampling-factor             */\
\
  unsigned flags;                       /**< flags for precomputation, malloc*/\
\
  int *index_sparse_to_full;            /**< index conversation,\
                                             overflow for d=3, J=9!          */\
\
  int r_act_nfft_plan;                  /**< index of current nfft block     */\
  Z(plan) *act_nfft_plan;             /**< current nfft block              */\
  Z(plan) *center_nfft_plan;          /**< central nfft block              */\
\
  Y(plan) *set_fftw_plan1;            /**< fftw plan for the nfft blocks   */\
  Y(plan) *set_fftw_plan2;            /**< fftw plan for the nfft blocks   */\
\
  Z(plan) *set_nfft_plan_1d;          /**< nfft plans for short nffts      */\
  Z(plan) *set_nfft_plan_2d;          /**< nfft plans for short nffts      */\
\
  R *x_transposed;                 /**< coordinate exchanged nodes, d=2 */\
  R *x_102,*x_201,*x_120,*x_021;   /**< coordinate exchanged nodes, d=3 */\
\
} X(plan);\
\
NFFT_EXTERN void X(direct_trafo)(X(plan) *ths); \
NFFT_EXTERN void X(direct_adjoint)(X(plan) *ths); \
NFFT_EXTERN void X(trafo)(X(plan) *ths); \
NFFT_EXTERN void X(adjoint)(X(plan) *ths); \
NFFT_EXTERN void X(cp)(X(plan) *ths, Z(plan) *ths_nfft); \
NFFT_EXTERN void X(init_random_nodes_coeffs)(X(plan) *ths); \
NFFT_EXTERN void X(init)(X(plan) *ths, int d, int J, int M, int m, unsigned flags); \
NFFT_EXTERN void X(finalize)(X(plan) *ths);

#ifdef HAVE_NSFFT
NSFFT_DEFINE_API(NSFFT_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,float,fftwf_complex)
NSFFT_DEFINE_API(NSFFT_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,double,fftw_complex)
NSFFT_DEFINE_API(NSFFT_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)
#endif

#define NSDFT            (1U<< 12)

/** @}
 */

/*###########################################################################*/

#define MRI_MANGLE_DOUBLE(name) NFFT_CONCAT(mri_, name)
#define MRI_MANGLE_FLOAT(name) NFFT_CONCAT(mrif_, name)
#define MRI_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(mril_, name)

#define MRI_DEFINE_API(X,Y,Z,R,C) \
typedef struct\
{\
  /* api */\
  MACRO_MV_PLAN(C)\
\
  Z(plan) plan;\
\
  int N3;\
  R sigma3;\
  R *t;\
  R *w;\
} X(inh_2d1d_plan);\
\
typedef struct\
{\
  /* api */\
  MACRO_MV_PLAN(C)\
\
  Z(plan) plan;\
\
  int N3;\
  R sigma3;\
  R *t;\
  R *w;\
} X(inh_3d_plan);\
\
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
 */ \
void X(inh_2d1d_trafo)(X(inh_2d1d_plan) *ths); \
\
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
 */ \
void X(inh_2d1d_adjoint)(X(inh_2d1d_plan) *ths); \
\
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
 */ \
void X(inh_2d1d_init_guru)(X(inh_2d1d_plan) *ths, int *N, int M, int *n, \
                    int m, R sigma, unsigned nfft_flags, unsigned fftw_flags); \
\
/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */ \
void X(inh_2d1d_finalize)(X(inh_2d1d_plan) *ths); \
\
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
 */ \
void X(inh_3d_trafo)(X(inh_3d_plan) *ths); \
\
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
 */ \
void X(inh_3d_adjoint)(X(inh_3d_plan) *ths); \
\
void X(inh_3d_init_guru)(X(inh_3d_plan) *ths, int *N, int M, int *n, \
                    int m, R sigma, unsigned nfft_flags, unsigned fftw_flags); \
\
/**
 * Destroys a plan.
 *
 * \arg ths_plan The plan
 *
 * \author Tobias Knopp
 */ \
void X(inh_3d_finalize)(X(inh_3d_plan) *ths);

#ifdef HAVE_MRI
MRI_DEFINE_API(MRI_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,float,fftwf_complex)
MRI_DEFINE_API(MRI_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,double,fftw_complex)
MRI_DEFINE_API(MRI_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)
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
 * \li By default, the NDSFT transforms (see \ref nfsft_direct_trafo, \ref nfsft_trafo)
 *   are allowed to destroy the input \c f_hat while the input \c x is
 *   preserved. On the contrary, the adjoint NDSFT transforms
 *   (see \ref nfsft_direct_adjoint, \ref nfsft_adjoint) do not destroy the input
 *   \c f and \c x by default. The desired behaviour can be assured by using the
 *   \ref NFSFT_PRESERVE_F_HAT, \ref NFSFT_PRESERVE_X, \ref NFSFT_PRESERVE_F and
 *   \ref NFSFT_DESTROY_F_HAT, \ref NFSFT_DESTROY_X, \ref NFSFT_DESTROY_F
 *   flags.
 */

#define NFSFT_MANGLE_DOUBLE(name) NFFT_CONCAT(nfsft_, name)
#define NFSFT_MANGLE_FLOAT(name) NFFT_CONCAT(nfsftf_, name)
#define NFSFT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nfsftl_, name)

#define NFSFT_DEFINE_API(X,Y,Z,R,C) \
/** Structure for a NFSFT transform plan */ \
typedef struct\
{\
  /** Inherited public members */\
  MACRO_MV_PLAN(C)\
\
  /* Public members */\
  int N;                              /**< the bandwidth \f$N\f$              */\
  R *x;                          /**< the nodes \f$\mathbf{x}(m) =       *\
                                           \left(x_1,x_2\right) \in           *\
                                           [-\frac{1}{2},\frac{1}{2}) \times  *\
                                           [0,\frac{1}{2}]\f$ for             *\
                                           \f$m=0,\ldots,M-1\f$,\f$M \in      *\
                                           \mathbb{N},\f$                     */\
\
  /* Private members */\
  /*int NPT;*/                        /**< the next greater power of two with *\
                                           respect to \f$N\f$                 */\
  int t;                              /**< the logarithm of NPT with           *\
                                           respect to the basis 2             */\
  unsigned int flags;                 /**< the planner flags                  */\
  Z(plan) plan_nfft;                /**< the internal NFFT plan             */\
  C *f_hat_intern;              /**< Internally used pointer to         *\
                                           spherical Fourier coefficients     */\
} X(plan);\
\
/**
 * Creates a transform plan.
 *
 * \arg plan a pointer to a \ref nfsft_plan structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 *
 * \author Jens Keiner
 */ \
NFFT_EXTERN void X(init)(X(plan) *plan, int N, int M); \
\
/**
 * Creates a transform plan.
 *
 * \arg plan a pointer to a \verbatim nfsft_plan \endverbatim structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 * \arg nfsft_flags the NFSFT flags
 *
 * \author Jens Keiner
 */ \
NFFT_EXTERN void X(init_advanced)(X(plan)* plan, int N, int M, unsigned int \
                         nfsft_flags); \
\
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
 */ \
NFFT_EXTERN void X(init_guru)(X(plan) *plan, int N, int M, unsigned int nfsft_flags, \
    unsigned int nfft_flags, int nfft_cutoff); \
\
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
 */ \
NFFT_EXTERN void X(precompute)(int N, R kappa, unsigned int nfsft_flags, \
  unsigned int fpt_flags); \
\
/**
 * Forgets all precomputed data.
 *
 * \author Jens Keiner
 */ \
NFFT_EXTERN void X(forget)(void); \
\
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
 */ \
NFFT_EXTERN void X(direct_trafo)(X(plan)* plan); \
\
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
 */ \
NFFT_EXTERN void X(direct_adjoint)(X(plan)* plan); \
\
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
 */ \
NFFT_EXTERN void X(trafo)(X(plan)* plan); \
\
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
 */ \
NFFT_EXTERN void X(adjoint)(X(plan)* plan); \
\
/**
 * Destroys a plan.
 *
 * \arg plan the plan to be destroyed
 *
 * \author Jens Keiner
 */ \
NFFT_EXTERN void X(finalize)(X(plan) *plan); \
\
NFFT_EXTERN void X(precompute_x)(X(plan) *plan);

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
 * \ref nfsft_direct_trafo or \ref nfsft_trafo the content of \c f_hat remains
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
 * \ref nfsft_direct_trafo, \ref nfsft_trafo or \ref nfsft_direct_adjoint, \ref nfsft_adjoint
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
 * \ref nfsft_direct_adjoint or \ref nfsft_adjoint the content of \c f remains
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
 * \ref nfsft_direct_trafo or \ref nfsft_trafo the content of \c f_hat may be changed.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_DESTROY_F_HAT    (1U << 10)

/**
 * If this flag is set, it is explicitely allowed that during an execution of
 * \ref nfsft_direct_trafo, \ref nfsft_trafo or \ref nfsft_direct_adjoint, \ref nfsft_adjoint
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
 * \ref nfsft_direct_adjoint or \ref nfsft_adjoint the content of \c f may be changed.
 *
 * \see nfsft_init
 * \see nfsft_init_advanced
 * \see nfsft_init_guru
 * \author Jens Keiner
 */
#define NFSFT_DESTROY_F      (1U << 12)

/* Precomputation flags */

/**
 * If this flag is set, the transforms \ref nfsft_direct_trafo and \ref nfsft_direct_adjoint
 * do not work. Setting this flag saves some memory for precomputed data.
 *
 * \see nfsft_precompute
 * \see nfsft_direct_trafo
 * \see nfsft_direct_adjoint
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
 * \ref nfsft_direct_adjoint set all unused entries in \c f_hat not corresponding to
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

NFSFT_DEFINE_API(NFSFT_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,float,fftwf_complex)
NFSFT_DEFINE_API(NFSFT_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,double,fftw_complex)
NFSFT_DEFINE_API(NFSFT_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)
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

#define FPT_MANGLE_DOUBLE(name) NFFT_CONCAT(fpt_, name)
#define FPT_MANGLE_FLOAT(name) NFFT_CONCAT(fptf_, name)
#define FPT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(fptl_, name)

#define FPT_DEFINE_API(X,Y,R,C) \
/* Data structures */\
typedef struct X(set_s_) *X(set);     /**< A set of precomputed data for a\
                                             set of DPT transforms of equal\
                                             maximum length.                  */\
\
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
 */ \
NFFT_EXTERN X(set) X(init)(const int M, const int t, const unsigned int flags); \
\
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
 */ \
NFFT_EXTERN void X(precompute)(X(set) set, const int m, R *alpha, R *beta, \
  R *gam, int k_start, const R threshold); \
\
/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */ \
NFFT_EXTERN void X(direct_trafo)(X(set) set, const int m, const C *x, C *y, \
  const int k_end, const unsigned int flags); \
\
/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */ \
NFFT_EXTERN void X(trafo)(X(set) set, const int m, const C *x, C *y, \
  const int k_end, const unsigned int flags); \
\
/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */ \
NFFT_EXTERN void X(direct_transposed)(X(set) set, const int m, C *x, \
  C *y, const int k_end, const unsigned int flags); \
\
/**
 * Computes a single DPT transform.
 *
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */ \
NFFT_EXTERN void X(transposed)(X(set) set, const int m, C *x, \
  C *y, const int k_end, const unsigned int flags); \
\
NFFT_EXTERN void X(finalize)(X(set) set);

#ifdef HAVE_FPT

/* Flags for fpt_init() */
#define FPT_NO_FAST_ALGORITHM (1U << 2) /**< If set, TODO complete comment.   */
#define FPT_NO_DIRECT_ALGORITHM (1U << 3) /**< If set, TODO complete comment.   */
#define FPT_NO_STABILIZATION  (1U << 0) /**< If set, no stabilization will be
                                             used.                            */

#define FPT_PERSISTENT_DATA   (1U << 4) /**< If set, TODO complete comment.   */

/* Flags for fpt_direct_trafo(), fpt_direct_transposed(), fpt_trafo(), fpt_transposed() */
#define FPT_FUNCTION_VALUES   (1U << 5) /**< If set, the output are function
                                             values at Chebyshev nodes rather
                                             than Chebyshev coefficients.     */
#define FPT_AL_SYMMETRY       (1U << 6) /**< TODO Don't use this flag!        */

FPT_DEFINE_API(FPT_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,float,fftwf_complex)
FPT_DEFINE_API(FPT_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,double,fftw_complex)
FPT_DEFINE_API(FPT_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,long double,fftwl_complex)
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

#define NFSOFT_MANGLE_DOUBLE(name) NFFT_CONCAT(nfsoft_, name)
#define NFSOFT_MANGLE_FLOAT(name) NFFT_CONCAT(nfsoftf_, name)
#define NFSOFT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nfsoftl_, name)

#define NFSOFT_DEFINE_API(X,Y,Z,R,C) \
typedef struct X(plan_)\
{\
  /** Inherited public members */\
  MACRO_MV_PLAN(C)\
\
  R *x;                           /**< the  input nodes                    */\
  /**some auxillary memory*/\
  C *wig_coeffs;          /**< contains a set of SO(3) Fourier coefficients*\
                                    for fixed orders m and n*/\
  C *cheby;           /**< contains a set of Chebychev coefficients for*\
                                    fixed orders m and n*/\
  C *aux;             /**< used when converting Chebychev to Fourier*\
                                    coeffcients*/\
\
  /** Private members */\
  int t;                               /**< the logaritm of NPT with          *\
                                          respect to the basis 2              */\
  unsigned int flags;                  /**< the planner flags                 */\
  Y(plan) p_nfft;                /**< the internal NFFT plan             */\
  Z(set) internal_fpt_set;                    /**< the internal FPT plan */\
\
  int fpt_kappa;       /**a parameter controlling the accuracy of the FPT*/\
\
} X(plan);\
\
/** Functions for NFSOFT plans*/\
\
/**
 * Does all node-dependent and node-independent precomputations needed for the NFSOFT.
 *
 * \arg plan a pointer to a \ref nfsoft_plan structure
 */\
NFFT_EXTERN void X(precompute)(X(plan) *plan); \
\
/**
 * Computes the FPT transform.
 *
 * \arg coeffs the Chebychev coefficients that should be transformed
 * \arg set the FPT-set containing precomputed data
 * \arg l the polynomial degree
 * \arg k the first order
 * \arg m the second order
 * \arg nfsoft_flags
 */ \
NFFT_EXTERN Z(set) X(SO3_single_fpt_init)(int l, int k, int m, unsigned int flags, int kappa); \
NFFT_EXTERN void X(SO3_fpt)(C *coeffs, Z(set) set, int l, int k, int m, unsigned int nfsoft_flags); \
NFFT_EXTERN void X(SO3_fpt_transposed)(C *coeffs, Z(set) set,int l, int k, int m,unsigned int nfsoft_flags); \
\
/**
 * Creates a NFSOFT transform plan.
 *
 * \arg plan a pointer to a \ref nfsoft_plan structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 *
 * \author Antje Vollrath
 */ \
NFFT_EXTERN void X(init)(X(plan) *plan, int N, int M); \
/**
 * Creates a NFSOFT transform plan.
 *
 * \arg plan a pointer to a \ref nfsoft_plan structure
 * \arg N the bandwidth \f$N \in \mathbb{N}_0\f$
 * \arg M the number of nodes \f$M \in \mathbb{N}\f$
 * \arg nfsoft_flags the NFSOFT flags
 *
 * \author Antje Vollrath
 */ \
NFFT_EXTERN void X(init_advanced)(X(plan) *plan, int N, int M,unsigned int nfsoft_flags); \
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
 */ \
NFFT_EXTERN void X(init_guru)(X(plan) *plan, int N, int M,unsigned int nfsoft_flags,unsigned int nfft_flags,int nfft_cutoff,int fpt_kappa); \
\
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
 */ \
NFFT_EXTERN void X(trafo)(X(plan) *plan_nfsoft); \
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
 */ \
NFFT_EXTERN void X(adjoint)(X(plan) *plan_nfsoft); \
/**
 * Destroys a plan.
 *
 * \arg plan the plan to be destroyed
 *
 * \author Antje Vollrath
 */ \
NFFT_EXTERN void X(finalize)(X(plan) *plan); \
\
NFFT_EXTERN int X(posN)(int n,int m, int B);

#ifdef HAVE_NFSOFT

NFSOFT_DEFINE_API(NFSOFT_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,FPT_MANGLE_FLOAT,float,fftwf_complex)
NFSOFT_DEFINE_API(NFSOFT_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,FPT_MANGLE_DOUBLE,double,fftw_complex)
NFSOFT_DEFINE_API(NFSOFT_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,FPT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

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
 * the spherical harmonics as defined in the NFFT software package, i.e.
 * for every rotation matrix \f$A$\f with Euler angles \f$\alpha, \beta, \gamma$\f
 * and every unit vector \f$x$\f the Wigner-D functions will be normed such that
 *
 * \f[
 *  \sum_{m=-l}^l D_{mn}^l(\alpha,\beta,\gamma) Y_m^l(x) = Y_n^l(A^{-1} x)
 * \f]
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
#define NFSOFT_INDEX(m,n,l,B) (((l)+((B)+1))+(2*(B)+2)*(((n)+((B)+1))+(2*(B)+2)*((m)+((B)+1))))
#define NFSOFT_INDEX_TWO(m,n,l,B) ((B+1)*(B+1)+(B+1)*(B+1)*(m+B)-((m-1)*m*(2*m-1)+(B+1)*(B+2)*(2*B+3))/6)+(posN(n,m,B))+(l-MAX(ABS(m),ABS(n)))

/**
 * This macro expands to the logical size of a SO(3) Fourier coefficients
 * array for a bandwidth B.
 */
#define NFSOFT_F_HAT_SIZE(B)          (((B)+1)*(4*((B)+1)*((B)+1)-1)/3)

/** Structure for a NFSOFT transform plan */

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

#define SOLVER_MANGLE_DOUBLE(name) NFFT_CONCAT(solver_, name)
#define SOLVER_MANGLE_FLOAT(name) NFFT_CONCAT(solverf_, name)
#define SOLVER_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(solverl_, name)

#define SOLVER_DEFINE_API(X,Y,R,C)\
typedef struct\
{\
  Y(mv_plan_complex) *mv;                  /**< matrix vector multiplication   */\
  unsigned flags;                       /**< iteration type                 */\
\
  R *w;                            /**< weighting factors              */\
  R *w_hat;                        /**< damping factors                */\
\
  C *y;                      /**< right hand side, samples       */\
\
  C *f_hat_iter;             /**< iterative solution             */\
\
  C *r_iter;                 /**< iterated residual vector       */\
  C *z_hat_iter;             /**< residual of normal equation of\
               first kind                     */\
  C *p_hat_iter;             /**< search direction               */\
  C *v_iter;                 /**< residual vector update         */\
\
  R alpha_iter;                    /**< step size for search direction */\
  R beta_iter;                     /**< step size for search correction*/\
\
  R dot_r_iter;                    /**< weighted dotproduct of r_iter  */\
  R dot_r_iter_old;                /**< previous dot_r_iter            */\
  R dot_z_hat_iter;                /**< weighted dotproduct of\
               z_hat_iter                     */\
  R dot_z_hat_iter_old;            /**< previous dot_z_hat_iter        */\
  R dot_p_hat_iter;                /**< weighted dotproduct of\
               p_hat_iter                     */\
  R dot_v_iter;                    /**< weighted dotproduct of v_iter  */\
} X(plan_complex);\
\
NFFT_EXTERN void X(init_advanced_complex)(X(plan_complex)* ths, Y(mv_plan_complex) *mv, unsigned flags);\
NFFT_EXTERN void X(init_complex)(X(plan_complex)* ths, Y(mv_plan_complex) *mv);\
NFFT_EXTERN void X(before_loop_complex)(X(plan_complex)* ths);\
NFFT_EXTERN void X(loop_one_step_complex)(X(plan_complex) *ths);\
NFFT_EXTERN void X(finalize_complex)(X(plan_complex) *ths);\
\
typedef struct\
{\
  Y(mv_plan_double) *mv;                   /**< matrix vector multiplication   */\
  unsigned flags;                       /**< iteration type                 */\
\
  R *w;                            /**< weighting factors              */\
  R *w_hat;                        /**< damping factors                */\
\
  R *y;                            /**< right hand side, samples       */\
\
  R *f_hat_iter;                   /**< iterative solution             */\
\
  R *r_iter;                       /**< iterated residual vector       */\
  R *z_hat_iter;                   /**< residual of normal equation of\
               first kind                     */\
  R *p_hat_iter;                   /**< search direction               */\
  R *v_iter;                       /**< residual vector update         */\
\
  R alpha_iter;                    /**< step size for search direction */\
  R beta_iter;                     /**< step size for search correction*/\
\
  R dot_r_iter;                    /**< weighted dotproduct of r_iter  */\
  R dot_r_iter_old;                /**< previous dot_r_iter            */\
  R dot_z_hat_iter;                /**< weighted dotproduct of\
               z_hat_iter                     */\
  R dot_z_hat_iter_old;            /**< previous dot_z_hat_iter        */\
  R dot_p_hat_iter;                /**< weighted dotproduct of\
               p_hat_iter                     */\
  R dot_v_iter;                    /**< weighted dotproduct of v_iter  */\
} X(plan_double);\
\
NFFT_EXTERN void X(init_advanced_double)(X(plan_double)* ths, Y(mv_plan_double) *mv, unsigned flags);\
NFFT_EXTERN void X(solver_init_double)(X(plan_double)* ths, Y(mv_plan_double) *mv);\
NFFT_EXTERN void X(solver_before_loop_double)(X(plan_double)* ths);\
NFFT_EXTERN void X(solver_loop_one_step_double)(X(plan_double) *ths);\
NFFT_EXTERN void X(solver_finalize_double)(X(plan_double) *ths);

SOLVER_DEFINE_API(SOLVER_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,float,fftwf_complex)
SOLVER_DEFINE_API(SOLVER_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,double,fftw_complex)
SOLVER_DEFINE_API(SOLVER_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

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

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

/** @}
 */

#endif
/* nfft3.h */
