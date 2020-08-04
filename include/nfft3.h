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

#ifndef __NFFT3_H__
#define __NFFT3_H__

/* fftw_complex */
#include <fftw3.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#define NFFT_CONCAT(prefix, name) prefix ## name

/* IMPORTANT: for Windows compilers, you should add a line
 *   #define NFFT_DLL
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

/* Integral type large enough to contain a stride (what ``int'' should have been
 * in the first place) */
typedef ptrdiff_t NFFT_INT;

/* Members inherited by all plans. */
#define MACRO_MV_PLAN(RC) \
  NFFT_INT N_total; /**< Total number of Fourier coefficients. */\
  NFFT_INT M_total; /**< Total number of samples. */\
  RC *f_hat; /**< Fourier coefficients. */\
  RC *f; /**< Samples. */\
  void (*mv_trafo)(void*); /**< Transform. */\
  void (*mv_adjoint)(void*); /**< Adjoint transform. */

/* nfft */

/* Name mangling macros. */
#define NFFT_MANGLE_DOUBLE(name) NFFT_CONCAT(nfft_, name)
#define NFFT_MANGLE_FLOAT(name) NFFT_CONCAT(nfftf_, name)
#define NFFT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nfftl_, name)

#define NFFT_DEFINE_MALLOC_API(X) \
/* our own memory allocation and exit functions */ \
NFFT_EXTERN void *X(malloc)(size_t n); \
NFFT_EXTERN void X(free)(void *p); \
NFFT_EXTERN void X(die)(const char *s); \
\
/* You can replace the hooks with your own functions, if necessary. We */ \
/* need this for the Matlab interface. */ \
typedef void *(*X(malloc_type_function)) (size_t n); \
typedef void  (*X(free_type_function)) (void *p); \
typedef void  (*X(die_type_function)) (const char *errString); \
NFFT_EXTERN X(malloc_type_function) X(malloc_hook); \
NFFT_EXTERN X(free_type_function) X(free_hook); \
NFFT_EXTERN X(die_type_function) X(die_hook);

/* Nfft module API. */
NFFT_DEFINE_MALLOC_API(NFFT_MANGLE_FLOAT)
NFFT_DEFINE_MALLOC_API(NFFT_MANGLE_DOUBLE)
NFFT_DEFINE_MALLOC_API(NFFT_MANGLE_LONG_DOUBLE)

/* Macro to define prototypes for all NFFT API functions.
 * We expand this macro for each supported precision.
 *   X: NFFT name-mangling macro
 *   Y: FFTW name-mangling macro
 *   R: Real float type
 *   C: Complex float type
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
/** data structure for an NFFT (nonequispaced fast Fourier transform) plan with R precision */ \
typedef struct\
{\
  MACRO_MV_PLAN(C)\
\
  NFFT_INT d; /**< Dimension (rank). */\
  NFFT_INT *N; /**< Multi-bandwidth. */\
  R *sigma; /**< Oversampling factor. */\
  NFFT_INT *n; /**< Length of FFTW transforms. This is equal to sigma*N. The default
               is to use a power of two that satifies  \f$2\le\sigma<4\f$. */\
  NFFT_INT n_total; /**< Combined total length of FFTW transforms. */\
  NFFT_INT m; /**< Cut-off parameter for window function. Default values for the
              different window functions are
              -  6 (KAISER_BESSEL),
              -  9 (SINC_POWER),
              - 11 (B_SPLINE),
              - 12 (GAUSSIAN) */\
  R *b; /**< Shape parameter for window function */\
  NFFT_INT K; /**< Number of equispaced samples of window function. Used for flag
             PRE_LIN_PSI. */\
\
  unsigned flags; /**< Flags for precomputation, (de)allocation, and FFTW
                       usage, default setting is PRE_PHI_HUT | PRE_PSI
                       | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT
                       | FFT_OUT_OF_PLACE */\
\
  unsigned fftw_flags; /**< Flags for the FFTW, default is
                            FFTW_ESTIMATE | FFTW_DESTROY_INPUT */\
\
  R *x; /**< Nodes in time/spatial domain, size is \f$dM\f$ R ## s */\
\
  R MEASURE_TIME_t[3]; /**< Measured time for each step if MEASURE_TIME is
    set */\
\
  /* internal use only */\
  Y(plan) my_fftw_plan1; /**< Forward FFTW plan */\
  Y(plan) my_fftw_plan2; /**< Backward FFTW plan */\
\
  R **c_phi_inv; /**< Precomputed data for the diagonal matrix \f$D\f$, size \
    is \f$N_0+\dots+N_{d-1}\f$ doubles*/\
  R *psi; /**< Precomputed data for the sparse matrix \f$B\f$, size depends
                    on precomputation scheme */\
  NFFT_INT *psi_index_g; /**< Indices in source/target vector for \ref PRE_FULL_PSI */\
  NFFT_INT *psi_index_f; /**< Indices in source/target vector for \ref PRE_FULL_PSI */\
\
  C *g; /**< Oversampled vector of samples, size is \ref n_total double complex */\
  C *g_hat; /**< Zero-padded vector of Fourier coefficients, size is \ref n_total fftw_complex */\
  C *g1; /**< Input of fftw */\
  C *g2; /**< Output of fftw */\
\
  R *spline_coeffs; /**< Input for de Boor algorithm if B_SPLINE or SINC_POWER is defined */\
\
  NFFT_INT *index_x; /**< Index array for nodes x used when flag \ref NFFT_SORT_NODES is set. */\
} X(plan); \
\
NFFT_EXTERN void X(trafo_direct)(const X(plan) *ths);\
NFFT_EXTERN void X(adjoint_direct)(const X(plan) *ths);\
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
NFFT_EXTERN void X(init_guru)(X(plan) *ths, int d, int *N, int M, int *n, \
  int m, unsigned flags, unsigned fftw_flags);\
NFFT_EXTERN void X(init_lin)(X(plan) *ths, int d, int *N, int M, int *n, \
  int m, int K, unsigned flags, unsigned fftw_flags); \
NFFT_EXTERN void X(precompute_one_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_full_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_fg_psi)(X(plan) *ths); \
NFFT_EXTERN void X(precompute_lin_psi)(X(plan) *ths);\
NFFT_EXTERN const char* X(check)(X(plan) *ths);\
NFFT_EXTERN void X(finalize)(X(plan) *ths);

/* Nfft module API. */
NFFT_DEFINE_API(NFFT_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,float,fftwf_complex)
NFFT_DEFINE_API(NFFT_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,double,fftw_complex)
NFFT_DEFINE_API(NFFT_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

/* Flags for init routines. */
#define PRE_PHI_HUT                (1U<<0)
#define FG_PSI                     (1U<<1)
#define PRE_LIN_PSI                (1U<<2)
#define PRE_FG_PSI                 (1U<<3)
#define PRE_PSI                    (1U<<4)
#define PRE_FULL_PSI               (1U<<5)
#define MALLOC_X                   (1U<<6)
#define MALLOC_F_HAT               (1U<<7)
#define MALLOC_F                   (1U<<8)
#define FFT_OUT_OF_PLACE           (1U<<9)
#define FFTW_INIT                  (1U<<10)
#define NFFT_SORT_NODES            (1U<<11)
#define NFFT_OMP_BLOCKWISE_ADJOINT (1U<<12)
#define PRE_ONE_PSI (PRE_LIN_PSI| PRE_FG_PSI| PRE_PSI| PRE_FULL_PSI)

/* nfct */

/* name mangling macros */
#define NFCT_MANGLE_DOUBLE(name) NFFT_CONCAT(nfct_, name)
#define NFCT_MANGLE_FLOAT(name) NFFT_CONCAT(nfctf_, name)
#define NFCT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nfctl_, name)

/* huge second-order macro that defines prototypes for all nfct API functions.
 * We expand this macro for each supported precision.
 *   X: nfct name-mangling macro
 *   Y: fftw name-mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define NFCT_DEFINE_API(X,Y,R,C) \
/** data structure for an NFCT (nonequispaced fast cosine transform) plan with R precision */ \
typedef struct\
{\
  /* api */\
  MACRO_MV_PLAN(R)\
\
  NFFT_INT d; /**< dimension, rank */\
  NFFT_INT *N; /**< cut-off-frequencies (kernel) */\
  NFFT_INT *n; /**< length of DCT-I */\
  NFFT_INT n_total; /**< Combined total length of FFTW transforms. */\
  R *sigma; /**< oversampling-factor */\
  NFFT_INT m; /**< cut-off parameter in time-domain */\
\
  R *b; /**< shape parameters */\
  NFFT_INT K; /**< Number of equispaced samples of window function. Used for flag
               PRE_LIN_PSI. */\
\
  unsigned flags; /**< flags for precomputation, malloc */\
  unsigned fftw_flags; /**< flags for the fftw */\
\
  R *x; /**< nodes (in time/spatial domain)   */\
\
  double MEASURE_TIME_t[3]; /**< measured time for each step */\
\
  /* internal use only */\
  Y(plan)  my_fftw_r2r_plan; /**< fftw_plan */\
  Y(r2r_kind) *r2r_kind; /**< r2r transform type (DCT-I) */\
\
  R **c_phi_inv; /**< precomputed data, matrix D */\
  R *psi; /**< precomputed data, matrix B */\
  NFFT_INT size_psi; /**< only for thin B */\
  NFFT_INT *psi_index_g; /**< only for thin B */\
  NFFT_INT *psi_index_f; /**< only for thin B */\
\
  R *g;\
  R *g_hat;\
  R *g1; /**< input of fftw */\
  R *g2; /**< output of fftw */\
\
  R *spline_coeffs; /**< input for de Boor algorithm, if B_SPLINE or SINC_2m is defined   */\
} X(plan);\
\
NFFT_EXTERN void X(init_1d)(X(plan) *ths_plan, int N0, int M_total); \
NFFT_EXTERN void X(init_2d)(X(plan) *ths_plan, int N0, int N1, int M_total); \
NFFT_EXTERN void X(init_3d)(X(plan) *ths_plan, int N0, int N1, int N2, int M_total); \
NFFT_EXTERN void X(init)(X(plan) *ths_plan, int d, int *N, int M_total); \
NFFT_EXTERN void X(init_guru)(X(plan) *ths_plan, int d, int *N, int M_total, int *n, \
  int m, unsigned flags, unsigned fftw_flags); \
NFFT_EXTERN void X(precompute_one_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_full_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_fg_psi)(X(plan) *ths); \
NFFT_EXTERN void X(precompute_lin_psi)(X(plan) *ths);\
NFFT_EXTERN void X(trafo)(X(plan) *ths_plan); \
NFFT_EXTERN void X(trafo_direct)(const X(plan) *ths_plan); \
NFFT_EXTERN void X(adjoint)(X(plan) *ths_plan); \
NFFT_EXTERN void X(adjoint_direct)(const X(plan) *ths_plan); \
NFFT_EXTERN const char* X(check)(X(plan) *ths);\
NFFT_EXTERN void X(finalize)(X(plan) *ths_plan); \

/* nfct api */
NFCT_DEFINE_API(NFCT_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,float,fftwf_complex)
NFCT_DEFINE_API(NFCT_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,double,fftw_complex)
NFCT_DEFINE_API(NFCT_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

/* nfst */

/* name mangling macros */
#define NFST_MANGLE_DOUBLE(name) NFFT_CONCAT(nfst_, name)
#define NFST_MANGLE_FLOAT(name) NFFT_CONCAT(nfstf_, name)
#define NFST_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nfstl_, name)

/* huge second-order macro that defines prototypes for all nfct API functions.
 * We expand this macro for each supported precision.
 *   X: nfst name-mangling macro
 *   Y: fftw name-mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define NFST_DEFINE_API(X,Y,R,C) \
/** data structure for an NFST (nonequispaced fast sine transform) plan with R precision */ \
typedef struct\
{\
  /* api */\
  MACRO_MV_PLAN(R)\
\
  NFFT_INT d; /**< dimension, rank */\
  NFFT_INT *N; /**< bandwidth */\
  NFFT_INT *n; /**< length of DST-I */\
  NFFT_INT n_total; /**< Combined total length of FFTW transforms. */\
  R *sigma; /**< oversampling-factor */\
  NFFT_INT m; /**< cut-off parameter in time-domain */\
\
  R *b; /**< shape parameters */\
  NFFT_INT K; /**< Number of equispaced samples of window function. Used for flag
               PRE_LIN_PSI. */\
\
  unsigned flags; /**< flags for precomputation, malloc */\
  unsigned fftw_flags; /**< flags for the fftw */\
\
  R *x; /**< nodes (in time/spatial domain) */\
\
  double MEASURE_TIME_t[3]; /**< measured time for each step */\
\
  /* internal use only */\
  Y(plan)  my_fftw_r2r_plan; /**< fftw_plan forward */\
  Y(r2r_kind) *r2r_kind; /**< r2r transform type (dct-i) */\
\
  R **c_phi_inv; /**< precomputed data, matrix D */\
  R *psi; /**< precomputed data, matrix B */\
  NFFT_INT size_psi; /**< only for thin B */\
  NFFT_INT *psi_index_g; /**< only for thin B */\
  NFFT_INT *psi_index_f; /**< only for thin B */\
\
  R *g;\
  R *g_hat;\
  R *g1; /**< input of fftw */\
  R *g2; /**< output of fftw */\
\
  R *spline_coeffs; /**< input for de Boor algorithm, if B_SPLINE or SINC_2m is defined */\
\
  R X(full_psi_eps);\
} X(plan);\
\
NFFT_EXTERN void X(init_1d)(X(plan) *ths_plan, int N0, int M_total); \
NFFT_EXTERN void X(init_2d)(X(plan) *ths_plan, int N0, int N1, int M_total); \
NFFT_EXTERN void X(init_3d)(X(plan) *ths_plan, int N0, int N1, int N2, int M_total); \
NFFT_EXTERN void X(init)(X(plan) *ths_plan, int d, int *N, int M_total); \
NFFT_EXTERN void X(init_guru)(X(plan) *ths_plan, int d, int *N, int M_total, int *n, \
  int m, unsigned flags, unsigned fftw_flags); \
NFFT_EXTERN void X(precompute_one_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_full_psi)(X(plan) *ths);\
NFFT_EXTERN void X(precompute_fg_psi)(X(plan) *ths); \
NFFT_EXTERN void X(precompute_lin_psi)(X(plan) *ths);\
NFFT_EXTERN void X(trafo)(X(plan) *ths_plan); \
NFFT_EXTERN void X(trafo_direct)(const X(plan) *ths_plan); \
NFFT_EXTERN void X(adjoint)(X(plan) *ths_plan); \
NFFT_EXTERN void X(adjoint_direct)(const X(plan) *ths_plan); \
NFFT_EXTERN const char* X(check)(X(plan) *ths);\
NFFT_EXTERN void X(finalize)(X(plan) *ths_plan); \

/* nfst api */
NFST_DEFINE_API(NFST_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,float,fftwf_complex)
NFST_DEFINE_API(NFST_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,double,fftw_complex)
NFST_DEFINE_API(NFST_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

/* nnfft */

/* name mangling macros */
#define NNFFT_MANGLE_DOUBLE(name) NFFT_CONCAT(nnfft_, name)
#define NNFFT_MANGLE_FLOAT(name) NFFT_CONCAT(nnfftf_, name)
#define NNFFT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nnfftl_, name)

/* huge second-order macro that defines prototypes for all nfst API functions.
 * We expand this macro for each supported precision.
 *   X: nnfft name-mangling macro
 *   Y: fftw name-mangling macro
 *   Z: nfft name mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define NNFFT_DEFINE_API(X,Y,Z,R,C) \
/** data structure for an NNFFT (nonequispaced in time and frequency fast Fourier transform) plan with R precision */ \
typedef struct\
{\
  /* api */\
  MACRO_MV_PLAN(C)\
\
  int d; /**< dimension, rank */\
  R *sigma; /**< oversampling-factor */\
  R *a; /**< 1 + 2*m/N1 */\
  int *N; /**< cut-off-frequencies */\
  int *N1; /**< sigma*N */\
  int *aN1; /**< sigma*a*N */\
  int m; /**< cut-off parameter in time-domain*/\
  R *b; /**< shape parameters */\
  int K; /**< number of precomp. uniform psi */\
  int aN1_total; /**< aN1_total=aN1[0]* ... *aN1[d-1] */\
  Z(plan) *direct_plan; /**< plan for the nfft */\
  unsigned nnfft_flags; /**< flags for precomputation, malloc*/\
  int *n; /**< n=N1, for the window function */\
  R *x; /**< nodes (in time/spatial domain) */\
  R *v; /**< nodes (in fourier domain) */\
  R *c_phi_inv; /**< precomputed data, matrix D */\
  R *psi; /**< precomputed data, matrix B */\
  int size_psi; /**< only for thin B */\
  int *psi_index_g; /**< only for thin B */\
  int *psi_index_f; /**< only for thin B */\
  C *F;\
  R *spline_coeffs; /**< input for de Boor algorithm, if B_SPLINE or SINC_2m is defined */\
} X(plan);\
\
NFFT_EXTERN void X(init)(X(plan) *ths_plan, int d, int N_total, int M_total, int *N); \
NFFT_EXTERN void X(init_1d)(X(plan) *ths_plan, int N, int M_total); \
NFFT_EXTERN void X(init_guru)(X(plan) *ths_plan, int d, int N_total, int M_total, \
  int *N, int *N1, int m, unsigned nnfft_flags); \
NFFT_EXTERN void X(trafo_direct)(X(plan) *ths_plan); \
NFFT_EXTERN void X(adjoint_direct)(X(plan) *ths_plan); \
NFFT_EXTERN void X(trafo)(X(plan) *ths_plan); \
NFFT_EXTERN void X(adjoint)(X(plan) *ths_plan); \
NFFT_EXTERN void X(precompute_lin_psi)(X(plan) *ths_plan); \
NFFT_EXTERN void X(precompute_psi)(X(plan) *ths_plan); \
NFFT_EXTERN void X(precompute_full_psi)(X(plan) *ths_plan); \
NFFT_EXTERN void X(precompute_phi_hut)(X(plan) *ths_plan); \
NFFT_EXTERN void X(precompute_one_psi)(X(plan) *ths);\
NFFT_EXTERN void X(finalize)(X(plan) *ths_plan);

/* nnfft api */
NNFFT_DEFINE_API(NNFFT_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,float,fftwf_complex)
NNFFT_DEFINE_API(NNFFT_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,double,fftw_complex)
NNFFT_DEFINE_API(NNFFT_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

/* additional init flags */
#define MALLOC_V         (1U<< 11)

/* nsfft */

#define NSFFT_MANGLE_DOUBLE(name) NFFT_CONCAT(nsfft_, name)
#define NSFFT_MANGLE_FLOAT(name) NFFT_CONCAT(nsfftf_, name)
#define NSFFT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nsfftl_, name)

/* huge second-order macro that defines prototypes for all nnfft API functions.
 * We expand this macro for each supported precision.
 *   X: nnfft name-mangling macro
 *   Y: fftw name-mangling macro
 *   Z: nfft name mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define NSFFT_DEFINE_API(X,Y,Z,R,C) \
/** data structure for an NSFFT (nonequispaced sparse fast Fourier transform) plan with R precision */ \
typedef struct\
{\
  MACRO_MV_PLAN(C)\
\
  int d; /**< dimension, rank; d = 2, 3 */\
  int J; /**< problem size, i.e.,
                d=2: N_total=(J+4) 2^(J+1)
                d=3: N_total=2^J 6(2^((J+1)/2+1)-1)+2^(3(J/2+1)) */\
  int sigma; /**< oversampling-factor */\
  unsigned flags; /**< flags for precomputation, malloc*/\
  int *index_sparse_to_full; /**< index conversation, overflow for d=3, J=9! */\
  int r_act_nfft_plan; /**< index of current nfft block */\
  Z(plan) *act_nfft_plan; /**< current nfft block */\
  Z(plan) *center_nfft_plan; /**< central nfft block */\
  Y(plan) *set_fftw_plan1; /**< fftw plan for the nfft blocks */\
  Y(plan) *set_fftw_plan2; /**< fftw plan for the nfft blocks */\
  Z(plan) *set_nfft_plan_1d; /**< nfft plans for short nffts */\
  Z(plan) *set_nfft_plan_2d; /**< nfft plans for short nffts */\
  R *x_transposed; /**< coordinate exchanged nodes, d = 2 */\
  R *x_102,*x_201,*x_120,*x_021; /**< coordinate exchanged nodes, d=3 */\
} X(plan);\
\
NFFT_EXTERN void X(trafo_direct)(X(plan) *ths); \
NFFT_EXTERN void X(adjoint_direct)(X(plan) *ths); \
NFFT_EXTERN void X(trafo)(X(plan) *ths); \
NFFT_EXTERN void X(adjoint)(X(plan) *ths); \
NFFT_EXTERN void X(cp)(X(plan) *ths, Z(plan) *ths_nfft); \
NFFT_EXTERN void X(init_random_nodes_coeffs)(X(plan) *ths); \
NFFT_EXTERN void X(init)(X(plan) *ths, int d, int J, int M, int m, unsigned flags); \
NFFT_EXTERN void X(finalize)(X(plan) *ths);

/* nsfft api */
NSFFT_DEFINE_API(NSFFT_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,float,fftwf_complex)
NSFFT_DEFINE_API(NSFFT_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,double,fftw_complex)
NSFFT_DEFINE_API(NSFFT_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

/* additional init flags */
#define NSDFT            (1U<< 12)

/* mri */

/* name mangling macros */
#define MRI_MANGLE_DOUBLE(name) NFFT_CONCAT(mri_, name)
#define MRI_MANGLE_FLOAT(name) NFFT_CONCAT(mrif_, name)
#define MRI_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(mril_, name)

/* huge second-order macro that defines prototypes for all mri API functions.
 * We expand this macro for each supported precision.
 *   X: mri name-mangling macro
 *   Z: nfft name mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define MRI_DEFINE_API(X,Z,R,C) \
typedef struct\
{\
  MACRO_MV_PLAN(C)\
  Z(plan) plan;\
  int N3;\
  R sigma3;\
  R *t;\
  R *w;\
} X(inh_2d1d_plan);\
\
typedef struct\
{\
  MACRO_MV_PLAN(C)\
  Z(plan) plan;\
  int N3;\
  R sigma3;\
  R *t;\
  R *w;\
} X(inh_3d_plan);\
\
void X(inh_2d1d_trafo)(X(inh_2d1d_plan) *ths); \
void X(inh_2d1d_adjoint)(X(inh_2d1d_plan) *ths); \
void X(inh_2d1d_init_guru)(X(inh_2d1d_plan) *ths, int *N, int M, int *n, \
  int m, R sigma, unsigned nfft_flags, unsigned fftw_flags); \
void X(inh_2d1d_finalize)(X(inh_2d1d_plan) *ths); \
void X(inh_3d_trafo)(X(inh_3d_plan) *ths); \
void X(inh_3d_adjoint)(X(inh_3d_plan) *ths); \
void X(inh_3d_init_guru)(X(inh_3d_plan) *ths, int *N, int M, int *n, \
  int m, R sigma, unsigned nfft_flags, unsigned fftw_flags); \
void X(inh_3d_finalize)(X(inh_3d_plan) *ths);

  /* mri api */
MRI_DEFINE_API(MRI_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,float,fftwf_complex)
MRI_DEFINE_API(MRI_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,double,fftw_complex)
MRI_DEFINE_API(MRI_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

/* nfsft */

/* name mangling macros */
#define NFSFT_MANGLE_DOUBLE(name) NFFT_CONCAT(nfsft_, name)
#define NFSFT_MANGLE_FLOAT(name) NFFT_CONCAT(nfsftf_, name)
#define NFSFT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nfsftl_, name)

/* huge second-order macro that defines prototypes for all nfsft API functions.
 * We expand this macro for each supported precision.
 *   X: nfsft name-mangling macro
 *   Z: nfft name mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define NFSFT_DEFINE_API(X,Z,R,C) \
/** data structure for an NFSFT (nonequispaced fast spherical Fourier transform) plan with R precision */ \
typedef struct\
{\
  MACRO_MV_PLAN(C)\
  int N; /**< the bandwidth \f$N\f$ */\
  R *x; /**< the nodes \f$\mathbf{x}(m) = \left(x_1,x_2\right) \in
    [-\frac{1}{2},\frac{1}{2}) \times [0,\frac{1}{2}]\f$ for \f$m=0,\ldots,
    M-1\f$,\f$M \in \mathbb{N},\f$ */\
  /* internal use only */\
  int t; /**< the logarithm of NPT with respect to the basis 2 */\
  unsigned int flags; /**< the planner flags */\
  Z(plan) plan_nfft; /**< the internal NFFT plan */\
  C *f_hat_intern; /**< Internally used pointer to spherical Fourier
    coefficients */\
  double MEASURE_TIME_t[3]; /**< Measured time for each step if MEASURE_TIME is
    set */\
} X(plan);\
\
NFFT_EXTERN void X(init)(X(plan) *plan, int N, int M); \
NFFT_EXTERN void X(init_advanced)(X(plan)* plan, int N, int M, unsigned int \
  nfsft_flags); \
NFFT_EXTERN void X(init_guru)(X(plan) *plan, int N, int M, \
  unsigned int nfsft_flags, unsigned int nfft_flags, int nfft_cutoff); \
NFFT_EXTERN void X(precompute)(int N, R kappa, unsigned int nfsft_flags, \
  unsigned int fpt_flags); \
NFFT_EXTERN void X(forget)(void); \
NFFT_EXTERN void X(trafo_direct)(X(plan)* plan); \
NFFT_EXTERN void X(adjoint_direct)(X(plan)* plan); \
NFFT_EXTERN void X(trafo)(X(plan)* plan); \
NFFT_EXTERN void X(adjoint)(X(plan)* plan); \
NFFT_EXTERN void X(finalize)(X(plan) *plan); \
NFFT_EXTERN void X(precompute_x)(X(plan) *plan);

/* nfsft api */
NFSFT_DEFINE_API(NFSFT_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,float,fftwf_complex)
NFSFT_DEFINE_API(NFSFT_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,double,fftw_complex)
NFSFT_DEFINE_API(NFSFT_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

/* init flags */
#define NFSFT_NORMALIZED     (1U << 0)
#define NFSFT_USE_NDFT       (1U << 1)
#define NFSFT_USE_DPT        (1U << 2)
#define NFSFT_MALLOC_X       (1U << 3)
#define NFSFT_MALLOC_F_HAT   (1U << 5)
#define NFSFT_MALLOC_F       (1U << 6)
#define NFSFT_PRESERVE_F_HAT (1U << 7)
#define NFSFT_PRESERVE_X     (1U << 8)
#define NFSFT_PRESERVE_F     (1U << 9)
#define NFSFT_DESTROY_F_HAT  (1U << 10)
#define NFSFT_DESTROY_X      (1U << 11)
#define NFSFT_DESTROY_F      (1U << 12)
#define NFSFT_EQUISPACED     (1U << 17)

/* precompute flags */
#define NFSFT_NO_DIRECT_ALGORITHM (1U << 13)
#define NFSFT_NO_FAST_ALGORITHM   (1U << 14)
#define NFSFT_ZERO_F_HAT          (1U << 16)

/* helper macros */
#define NFSFT_INDEX(k,n,plan) ((2*(plan)->N+2)*((plan)->N-n+1)+(plan)->N+k+1)
#define NFSFT_F_HAT_SIZE(N) ((2*N+2)*(2*N+2))

/* fpt */

/* name mangling macros */
#define FPT_MANGLE_DOUBLE(name) NFFT_CONCAT(fpt_, name)
#define FPT_MANGLE_FLOAT(name) NFFT_CONCAT(fptf_, name)
#define FPT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(fptl_, name)

/* huge second-order macro that defines prototypes for all fpt API functions.
 * We expand this macro for each supported precision.
 *   X: fpt name-mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define FPT_DEFINE_API(X,Y,R,C) \
typedef struct X(set_s_) *X(set); /**< A set of precomputed data for a set of
  DPT transforms of equal maximum length. */\
\
NFFT_EXTERN X(set) X(init)(const int M, const int t, const unsigned int flags); \
NFFT_EXTERN void X(precompute)(X(set) set, const int m, R *alpha, R *beta, \
  R *gam, int k_start, const R threshold); \
NFFT_EXTERN void X(trafo_direct)(X(set) set, const int m, const C *x, C *y, \
  const int k_end, const unsigned int flags); \
NFFT_EXTERN void X(trafo)(X(set) set, const int m, const C *x, C *y, \
  const int k_end, const unsigned int flags); \
NFFT_EXTERN void X(transposed_direct)(X(set) set, const int m, C *x, \
  C *y, const int k_end, const unsigned int flags); \
NFFT_EXTERN void X(transposed)(X(set) set, const int m, C *x, \
  C *y, const int k_end, const unsigned int flags); \
NFFT_EXTERN void X(finalize)(X(set) set);

/* fpt api */
FPT_DEFINE_API(FPT_MANGLE_FLOAT,FFTW_MANGLE_FLOAT,float,fftwf_complex)
FPT_DEFINE_API(FPT_MANGLE_DOUBLE,FFTW_MANGLE_DOUBLE,double,fftw_complex)
FPT_DEFINE_API(FPT_MANGLE_LONG_DOUBLE,FFTW_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

/* init flags */
#define FPT_NO_STABILIZATION    (1U << 0)
#define FPT_NO_FAST_ALGORITHM   (1U << 2)
#define FPT_NO_DIRECT_ALGORITHM (1U << 3)
#define FPT_PERSISTENT_DATA     (1U << 4)
#define FPT_NO_INIT_FPT_DATA    (1U << 7)

/* transform flags */
#define FPT_FUNCTION_VALUES     (1U << 5)
#define FPT_AL_SYMMETRY         (1U << 6)

/* nfsoft*/

/* name mangling macros */
#define NFSOFT_MANGLE_DOUBLE(name) NFFT_CONCAT(nfsoft_, name)
#define NFSOFT_MANGLE_FLOAT(name) NFFT_CONCAT(nfsoftf_, name)
#define NFSOFT_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(nfsoftl_, name)

/* huge second-order macro that defines prototypes for all nfsoft API functions.
 * We expand this macro for each supported precision.
 *   X: nfsoft name-mangling macro
 *   Y: nfft name-mangling macro
 *   Z: fpt name-mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define NFSOFT_DEFINE_API(X,Y,Z,R,C) \
typedef struct X(plan_)\
{\
  MACRO_MV_PLAN(C) \
  R *x; /**< input nodes */\
  /* internal use only */\
  C *wig_coeffs; /**< deprecated variable */\
  C *cheby; /**< deprecated variable */\
  C *aux; /**< deprecated variable */\
  int t; /**< the logarithm of NPT with respect to the basis 2 */\
  unsigned int flags; /**< the planner flags  */\
  Y(plan) p_nfft; /**< the internal NFFT plan */\
  Z(set) *internal_fpt_set; /**< the internal FPT plan */\
  int nthreads; /**< the number of threads */\
} X(plan);\
\
NFFT_EXTERN void X(precompute)(X(plan) *plan); \
NFFT_EXTERN Z(set) X(SO3_single_fpt_init)(int l, int k, int m, unsigned int flags, int kappa); \
NFFT_EXTERN void X(SO3_fpt)(C *coeffs, Z(set) set, int l, int k, int m, unsigned int nfsoft_flags); \
NFFT_EXTERN void X(SO3_fpt_transposed)(C *coeffs, Z(set) set,int l, int k, int m,unsigned int nfsoft_flags); \
NFFT_EXTERN void X(init)(X(plan) *plan, int N, int M); \
NFFT_EXTERN void X(init_advanced)(X(plan) *plan, int N, int M,unsigned int nfsoft_flags); \
NFFT_EXTERN void X(init_guru)(X(plan) *plan, int N, int M,unsigned int nfsoft_flags,unsigned int nfft_flags,int nfft_cutoff,int fpt_kappa); \
NFFT_EXTERN void X(init_guru_advanced)(X(plan) *plan, int N, int M,unsigned int nfsoft_flags,unsigned int nfft_flags,int nfft_cutoff,int fpt_kappa, int nn_oversampled); \
NFFT_EXTERN void X(trafo)(X(plan) *plan_nfsoft); \
NFFT_EXTERN void X(adjoint)(X(plan) *plan_nfsoft); \
NFFT_EXTERN void X(finalize)(X(plan) *plan); \
NFFT_EXTERN int X(posN)(int n,int m, int B);

/* nfsoft api */
NFSOFT_DEFINE_API(NFSOFT_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,FPT_MANGLE_FLOAT,float,fftwf_complex)
NFSOFT_DEFINE_API(NFSOFT_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,FPT_MANGLE_DOUBLE,double,fftw_complex)
NFSOFT_DEFINE_API(NFSOFT_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,FPT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

/* init flags */
#define NFSOFT_NORMALIZED     (1U << 0)
#define NFSOFT_USE_NDFT       (1U << 1)
#define NFSOFT_USE_DPT        (1U << 2)
#define NFSOFT_MALLOC_X       (1U << 3)
#define NFSOFT_REPRESENT      (1U << 4)
#define NFSOFT_MALLOC_F_HAT   (1U << 5)
#define NFSOFT_MALLOC_F       (1U << 6)
#define NFSOFT_PRESERVE_F_HAT (1U << 7)
#define NFSOFT_PRESERVE_X     (1U << 8)
#define NFSOFT_PRESERVE_F     (1U << 9)
#define NFSOFT_DESTROY_F_HAT  (1U << 10)
#define NFSOFT_DESTROY_X      (1U << 11)
#define NFSOFT_DESTROY_F      (1U << 12)

/* precompute flags */
#define NFSOFT_NO_STABILIZATION (1U << 13)
#define NFSOFT_CHOOSE_DPT       (1U << 14)
#define NFSOFT_SOFT             (1U << 15)
#define NFSOFT_ZERO_F_HAT       (1U << 16)

/* helper macros */
#define NFSOFT_INDEX(m,n,l,B) (((l)+((B)+1))+(2*(B)+2)*(((n)+((B)+1))+(2*(B)+2)*((m)+((B)+1))))
#define NFSOFT_F_HAT_SIZE(B) (((B)+1)*(4*((B)+1)*((B)+1)-1)/3)

/* solver */

/* name mangling macros */
#define SOLVER_MANGLE_DOUBLE(name) NFFT_CONCAT(solver_, name)
#define SOLVER_MANGLE_FLOAT(name) NFFT_CONCAT(solverf_, name)
#define SOLVER_MANGLE_LONG_DOUBLE(name) NFFT_CONCAT(solverl_, name)

/* huge second-order macro that defines prototypes for all nfsoft API functions.
 * We expand this macro for each supported precision.
 *   X: nfsoft name-mangling macro
 *   Y: nfft name-mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define SOLVER_DEFINE_API(X,Y,R,C)\
/** data structure for an inverse NFFT plan with R precision */ \
typedef struct\
{\
  Y(mv_plan_complex) *mv; /**< matrix vector multiplication   */\
  unsigned flags; /**< iteration type */\
  R *w; /**< weighting factors */\
  R *w_hat; /**< damping factors */\
  C *y; /**< right hand side, samples */\
  C *f_hat_iter; /**< iterative solution */\
  C *r_iter; /**< iterated residual vector */\
  C *z_hat_iter; /**< residual of normal equation of first kind */\
  C *p_hat_iter; /**< search direction */\
  C *v_iter; /**< residual vector update */\
  R alpha_iter; /**< step size for search direction */\
  R beta_iter; /**< step size for search correction*/\
  R dot_r_iter; /**< weighted dotproduct of r_iter */\
  R dot_r_iter_old; /**< previous dot_r_iter */\
  R dot_z_hat_iter; /**< weighted dotproduct of z_hat_iter */\
  R dot_z_hat_iter_old; /**< previous dot_z_hat_iter */\
  R dot_p_hat_iter; /**< weighted dotproduct of p_hat_iter */\
  R dot_v_iter; /**< weighted dotproduct of v_iter */\
} X(plan_complex);\
\
NFFT_EXTERN void X(init_advanced_complex)(X(plan_complex)* ths, Y(mv_plan_complex) *mv, unsigned flags);\
NFFT_EXTERN void X(init_complex)(X(plan_complex)* ths, Y(mv_plan_complex) *mv);\
NFFT_EXTERN void X(before_loop_complex)(X(plan_complex)* ths);\
NFFT_EXTERN void X(loop_one_step_complex)(X(plan_complex) *ths);\
NFFT_EXTERN void X(finalize_complex)(X(plan_complex) *ths);\
\
/** data structure for an inverse NFFT plan with R precision */ \
typedef struct\
{\
  Y(mv_plan_double) *mv; /**< matrix vector multiplication   */\
  unsigned flags; /**< iteration type */\
  R *w; /**< weighting factors */\
  R *w_hat; /**< damping factors */\
  R *y; /**< right hand side, samples */\
  R *f_hat_iter; /**< iterative solution */\
  R *r_iter; /**< iterated residual vector */\
  R *z_hat_iter; /**< residual of normal equation of first kind */\
  R *p_hat_iter; /**< search direction */\
  R *v_iter; /**< residual vector update */\
  R alpha_iter; /**< step size for search direction */\
  R beta_iter; /**< step size for search correction */\
  R dot_r_iter; /**< weighted dotproduct of r_iter */\
  R dot_r_iter_old; /**< previous dot_r_iter */\
  R dot_z_hat_iter; /**< weighted dotproduct of z_hat_iter */\
  R dot_z_hat_iter_old; /**< previous dot_z_hat_iter */\
  R dot_p_hat_iter; /**< weighted dotproduct of p_hat_iter */\
  R dot_v_iter; /**< weighted dotproduct of v_iter */\
} X(plan_double);\
\
NFFT_EXTERN void X(init_advanced_double)(X(plan_double)* ths, Y(mv_plan_double) *mv, unsigned flags);\
NFFT_EXTERN void X(init_double)(X(plan_double)* ths, Y(mv_plan_double) *mv);\
NFFT_EXTERN void X(before_loop_double)(X(plan_double)* ths);\
NFFT_EXTERN void X(loop_one_step_double)(X(plan_double) *ths);\
NFFT_EXTERN void X(finalize_double)(X(plan_double) *ths);

/* solver api */
SOLVER_DEFINE_API(SOLVER_MANGLE_FLOAT,NFFT_MANGLE_FLOAT,float,fftwf_complex)
SOLVER_DEFINE_API(SOLVER_MANGLE_DOUBLE,NFFT_MANGLE_DOUBLE,double,fftw_complex)
SOLVER_DEFINE_API(SOLVER_MANGLE_LONG_DOUBLE,NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

/* init flags */
#define LANDWEBER             (1U<< 0)
#define STEEPEST_DESCENT      (1U<< 1)
#define CGNR                  (1U<< 2)
#define CGNE                  (1U<< 3)
#define NORMS_FOR_LANDWEBER   (1U<< 4)
#define PRECOMPUTE_WEIGHT     (1U<< 5)
#define PRECOMPUTE_DAMP       (1U<< 6)

/* util */

/* huge second-order macro that defines prototypes for all utility API functions.
 * We expand this macro for each supported precision.
 *   Y: nfft name-mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define NFFT_DEFINE_UTIL_API(Y,R,C) \
/* rand.c */ \
R Y(drand48)(void); \
void Y(srand48)(long int seed); \
\
/** Inits a vector of random complex numbers in \f$[0,1]\times[0,1]{\rm i}\f$. \
 */ \
void Y(vrand_unit_complex)(C *x, const NFFT_INT n); \
\
/** Inits a vector of random double numbers in \f$[-1/2,1/2]\f$. \
 */ \
void Y(vrand_shifted_unit_double)(R *x, const NFFT_INT n); \
\
void Y(vrand_real)(R *x, const NFFT_INT n, const R a, const R b); \
\
/* print.c */ \
/** Print real vector to standard output. */ \
void Y(vpr_double)(R *x, const NFFT_INT n, const char *text); \
\
/** Print complex vector to standard output. */ \
void Y(vpr_complex)(C *x, const NFFT_INT n, const char *text); \
/* thread.c */ \
NFFT_INT Y(get_num_threads)(void); \
void Y(set_num_threads)(NFFT_INT nthreads); \
NFFT_INT Y(has_threads_enabled)(void); \
/* time.c */ \
R Y(clock_gettime_seconds)(void); \
/* error.c: */ \
R Y(error_l_infty_complex)(const C *x, const C *y, const NFFT_INT n); \
R Y(error_l_infty_1_complex)(const C *x, const C *y, const NFFT_INT n, \
  const C *z, const NFFT_INT m); \
/* int.c: */ \
NFFT_INT Y(exp2i)(const NFFT_INT a); \
NFFT_INT Y(next_power_of_2)(const NFFT_INT N); \
/* vector1.c */ \
/** Computes the inner/dot product \f$x^H x\f$. */ \
R Y(dot_complex)(C *x, NFFT_INT n); \
/* vector3.c */ \
/** Updates \f$x \leftarrow a x + y\f$. */ \
void Y(upd_axpy_complex)(C *x, R a, C *y, NFFT_INT n); \
/** Swaps each half over N[d]/2. */ \
void Y(fftshift_complex)(C *x, NFFT_INT d, NFFT_INT* N); \
void Y(fftshift_complex_int)(C *x, int d, int* N); \
/** Return library version. */ \
void Y(get_version)(unsigned *major, unsigned *minor, unsigned *patch); \
/** \
 * Return name of window function. \
 * \
 * The window function to be used is configured at compile time. \
 */ \
const char *Y(get_window_name)(); \
NFFT_INT Y(get_default_window_cut_off)();

NFFT_DEFINE_UTIL_API(NFFT_MANGLE_FLOAT,float,fftwf_complex)
NFFT_DEFINE_UTIL_API(NFFT_MANGLE_DOUBLE,double,fftw_complex)
NFFT_DEFINE_UTIL_API(NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* defined(__NFFT3_H__) */
