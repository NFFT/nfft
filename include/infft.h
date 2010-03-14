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

/* NFFT internal header file */
#ifndef __INFFT_H__
#define __INFFT_H__

#include "config.h"

#include <math.h>
#include <float.h>
#include <complex.h>
#include <stdio.h>
#include <string.h>

#include <stdlib.h> /* size_t */
#include <stdarg.h> /* va_list */
#include <stddef.h> /* ptrdiff_t */

#if HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#if HAVE_STDINT_H
#include <stdint.h> /* uintptr_t, maybe */
#endif

#if HAVE_INTTYPES_H
#include <inttypes.h> /* uintptr_t, maybe */
#endif

#if defined(HAVE_LIBDISPATCH)
#include <dispatch/dispatch.h> /* libdispatch */
#endif

/** Macros for window functions. */
#if defined(DIRAC_DELTA)
  #define PHI_HUT(k,d) K(1.0)
  #define PHI(x,d) (FABS((x))<K(10e-8))? K(1.0) : K(0.0)
  #define WINDOW_HELP_INIT(d)
  #define WINDOW_HELP_FINALIZE
  #define WINDOW_HELP_ESTIMATE_m {ths->m = 0;}
#elif defined(GAUSSIAN)
  #define PHI_HUT(k,d) ((R)EXP(-(POW(PI*(k)/ths->n[d],K(2.0))*ths->b[d])))
  #define PHI(x,d) ((R)EXP(-POW((x)*((R)ths->n[d]),K(2.0))/ ths->b[d])/SQRT(PI*ths->b[d]))
  #define WINDOW_HELP_INIT \
    {                                                                          \
      int WINDOW_idx;                                                          \
      ths->b = (R*) nfft_malloc(ths->d*sizeof(R));                   \
      for(WINDOW_idx=0; WINDOW_idx<ths->d; WINDOW_idx++)                       \
      ths->b[WINDOW_idx]=(K(2.0)*ths->sigma[WINDOW_idx])/                   \
        (K(2.0)*ths->sigma[WINDOW_idx]-K(1.0))*(((R)ths->m) / PI);                  \
      }
  #define WINDOW_HELP_FINALIZE {nfft_free(ths->b);}
  #define WINDOW_HELP_ESTIMATE_m {ths->m = 12;}
#elif defined(B_SPLINE)
  #define PHI_HUT(k,d) ((R)(((k)==0)? K(1.0)/ths->n[(d)] :                   \
    POW(SIN((k)*PI/ths->n[(d)])/((k)*PI/ths->n[(d)]),K(2.0)*ths->m)/ths->n[(d)]))
  #define PHI(x,d) (nfft_bspline(2*ths->m,((x)*ths->n[(d)])+                   \
    (R)ths->m,ths->spline_coeffs)/ths->n[(d)])
  #define WINDOW_HELP_INIT \
    {                                                                          \
      ths->spline_coeffs= (R*)nfft_malloc(2*ths->m*sizeof(R));       \
    }
  #define WINDOW_HELP_FINALIZE {nfft_free(ths->spline_coeffs);}
  #define WINDOW_HELP_ESTIMATE_m {ths->m = 11;}
#elif defined(SINC_POWER)
  #define PHI_HUT(k,d) (nfft_bspline(2*ths->m,(K(2.0)*ths->m*(k))/          \
    ((K(2.0)*ths->sigma[(d)]-1)*ths->n[(d)]/ths->sigma[(d)])+ (R)ths->m,       \
    ths->spline_coeffs))
  #define PHI(x,d) ((R)(ths->n[(d)]/ths->sigma[(d)]*(K(2.0)*ths->sigma[(d)]-K(1.0))/\
    (K(2.0)*ths->m)*POW(nfft_sinc(PI*ths->n[(d)]/ths->sigma[(d)]*(x)*               \
    (K(2.0)*ths->sigma[(d)]-1)/(K(2.0)*ths->m)),2*ths->m)/ths->n[(d)]))
  #define WINDOW_HELP_INIT \
    {                                                                          \
      ths->spline_coeffs= (R*)nfft_malloc(2*ths->m*sizeof(R));       \
    }
  #define WINDOW_HELP_FINALIZE {nfft_free(ths->spline_coeffs);}
  #define WINDOW_HELP_ESTIMATE_m {ths->m = 9;}
#else /* Kaiser-Bessel is the default. */
  #define PHI_HUT(k,d) ((R)nfft_i0( ths->m*SQRT(\
    POW((R)(ths->b[d]),K(2.0)) - POW(K(2.0)*PI*(k)/ths->n[d],K(2.0)))))
  #define PHI(x,d) ((R)((POW((R)(ths->m),K(2.0))\
    -POW((x)*ths->n[d],2.0))>0)? \
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

/* precision and name-mangling scheme */
#define CONCAT(prefix, name) prefix ## name
#if defined(NFFT_SINGLE)
typedef float R;
typedef float _Complex C;
#define X(name) CONCAT(nfftf_, name)
#elif defined(NFFT_LDOUBLE)
typedef long double R;
typedef long double _Complex C;
#define X(name) CONCAT(nfftl_, name)
#else
typedef double R;
typedef double _Complex C;
#define X(name) CONCAT(nfft_, name)
#endif

#ifdef NFFT_LDOUBLE
#  define K(x) ((R) x##L)
#else
#  define K(x) ((R) x)
#endif
#define DK(name, value) const R name = K(value)

/**
 * Integral type large enough to contain a stride (what ``int'' should
 * have been in the first place)
 */
typedef ptrdiff_t INT;

#define IF(x,a,b) ((x)?(a):(b))

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define ABS(x) (((x)>K(0.0))?(x):(-(x)))
#define SIGN(a) (((a)>=0)?1:-1)
#define SIGN(a) (((a)>=0)?1:-1)
#define SIGNF(a) IF((a)<K(0.0),K(-1.0),K(1.0))

#if defined(NFFT_LDOUBLE)
#define LOG1P log1pl
#if HAVE_DECL_LOG1PL == 0
extern long double log1pl(const long double);
#endif
#define LOG10 log10l
#if HAVE_DECL_LOG10L == 0
extern long double log10l(const long double);
#endif
#define SQRT sqrtl
#if HAVE_DECL_SQRTL == 0
extern long double sqrtl(const long double);
#endif
#define TGAMMA tgammal
#if HAVE_DECL_TGAMMAL == 0
extern long double tgammal(const long double);
#endif
#define LGAMMA lgammal
#if HAVE_DECL_LGAMMAL == 0
extern long double lgammal(const long double);
#endif
#define LOG logl
#if HAVE_DECL_LOGL == 0
extern long double logl(const long double);
#endif
#define LOG2 log2l
#if HAVE_DECL_LOG2L == 0
extern long double log2l(const long double);
#endif
#define SIN sinl
#if HAVE_DECL_SINL == 0
extern long double sinl(const long double);
#endif
#define COS cosl
#if HAVE_DECL_COSL == 0
extern long double cosl(const long double);
#endif
#define ACOS acosl
#if HAVE_DECL_ACOSL == 0
extern long double acosl(const long double);
#endif
#define POW powl
#if HAVE_DECL_POWL == 0
extern long double powl(const long double, const long double);
#endif
#define EXP expl
#if HAVE_DECL_EXPL == 0
extern long double expl(const long double);
#endif
#define CEXP cexpl
#if HAVE_DECL_CEXPL == 0
extern long double _Complex cexpl(const long double _Complex);
#endif
#define COPYSIGN copysignl
#if HAVE_DECL_COPYSIGNL == 0
extern long double copysignl(const long double, const long double);
#endif
#define LRINT lrintl
#if HAVE_DECL_LRINTL == 0
extern long int lrintl(const long double);
#endif
#define FMAX fmaxl
#if HAVE_DECL_FMAXL == 0
extern double fmaxl(const long double x, const long double y);
#endif
#elif defined(NFFT_SINGLE)
#define LOG1P log1pf
#if HAVE_DECL_LOG1PF == 0
extern float log1pf(const float);
#endif
#define LOG10 log10f
#if HAVE_DECL_LOG10F == 0
extern float log10f(const float);
#endif
#define SQRT sqrtf
#if HAVE_DECL_SQRTF == 0
extern float sqrtf(const float);
#endif
#define >TGAMMA tgammaf
#if HAVE_DECL_TGAMMAF == 0
extern float tgammaf(const float);
#endif
#define LGAMMA lgammaf
#if HAVE_DECL_LGAMMAF == 0
extern float lgammaf(const float);
#endif
#define LOG logf
#if HAVE_DECL_LOGF == 0
extern float logf(const float);
#endif
#define LOG2 log2f
#if HAVE_DECL_LOG2F == 0
extern float log2f(const float);
#endif
#define SIN sinf
#if HAVE_DECL_SINF == 0
extern float sinf(const float);
#endif
#define COS cosf
#if HAVE_DECL_COSF == 0
extern float cosf(const float);
#endif
#define ACOS acosf
#if HAVE_DECL_ACOSF == 0
extern float acosf(const float);
#endif
#define POW powf
#if HAVE_DECL_POWF == 0
extern float powf(const float, const float);
#endif
#define EXP expf
#if HAVE_DECL_EXPF == 0
extern float expf(const float);
#endif
#define CEXP cexpf
#if HAVE_DECL_CEXPF == 0
extern float _Complex cexpf(const float _Complex);
#endif
#define COPYSIGN copysignf
#if HAVE_DECL_COPYSIGNF == 0
extern float copysignf(const float, const float);
#endif
#define LRINT lrintf
#if HAVE_DECL_LRINTF == 0
extern long int lrintf(const float);
#endif
#define FMAX fmaxf
#if HAVE_DECL_FMAXF == 0
extern double fmaxf(const float x, const float y);
#endif
#else
#define LOG1P log1p
#if HAVE_DECL_LOG1P == 0
extern double log1p(const double);
#endif
#define LOG10 log10
#if HAVE_DECL_LOG10 == 0
extern double log10(const double);
#endif
#define SQRT sqrt
#if HAVE_DECL_SQRT == 0
extern double sqrt(const double);
#endif
#define TGAMMA tgamma
#if HAVE_DECL_TGAMMA == 0
extern double tgamma(const double);
#endif
#define LGAMMA lgamma
#if HAVE_DECL_LGAMMA == 0
extern double lgamma(const double);
#endif
#define LOG log
#if HAVE_DECL_LOG == 0
extern double log(const double);
#endif
#define LOG2 log2
#if HAVE_DECL_LOG2 == 0
extern double log2(const double);
#endif
#define SIN sin
#if HAVE_DECL_SIN == 0
extern double sin(const double);
#endif
#define COS cos
#if HAVE_DECL_COS == 0
extern double cos(const double);
#endif
#define ACOS acos
#if HAVE_DECL_ACOS == 0
extern double acos(const double);
#endif
#define POW pow
#if HAVE_DECL_POW == 0
extern double pow(const double, const double);
#endif
#define EXP exp
#if HAVE_DECL_EXP == 0
extern double exp(const double);
#endif
#define CEXP cexp
#if HAVE_DECL_CEXP == 0
extern double _Complex cexp(const double _Complex);
#endif
#define COPYSIGN copysign
#if HAVE_DECL_COPYSIGN == 0
extern double copysign(const double, const double);
#endif
#define LRINT lrint
#if HAVE_DECL_LRINT == 0
extern long int lrint(const double);
#endif
#define FMAX fmax
#if HAVE_DECL_FMAX == 0
extern double fmax(const double x, const double y);
#endif
#endif

#if defined(NFFT_LDOUBLE)
#  define CEIL ceill
#  define FLOOR floorl
#  define ROUND roundl
#  define FABS fabsl
#  define FMAX fmaxl
#  define FMIN fminl
#  define FREXP frexpl
#  define LDEXP ldexpl
#  define EPS LDBL_EPSILON
#  define R_MIN LDBL_MIN
#  define R_MAX LDBL_MAX
#  define R_MIN_EXP LDBL_MIN_EXP
#  define R_MAX_EXP LDBL_MAX_EXP
#  define R_MIN_10_EXP LDBL_MIN_10_EXP
#  define R_MAX_10_EXP LDBL_MAX_10_EXP
#  define R_DIG LDBL_DIG
#  define CREAL creall
#  define CIMAG cimagl
#elif defined(NFFT_SINGLE)
#  define CEIL ceilf
#  define FLOOR floorf
#  define ROUND roundf
#  define FABS fabsf
#  define FMAX fmaxf
#  define FMIN fminf
#  define FREXP frexpf
#  define LDEXP ldexpf
#  define EPS FLT_EPSILON
#  define R_MIN FLT_MIN
#  define R_MAX FLT_MAX
#  define R_MIN_EXP FLT_MIN_EXP
#  define R_MAX_EXP FLT_MAX_EXP
#  define R_MIN_10_EXP FLT_MIN_10_EXP
#  define R_MAX_10_EXP FLT_MAX_10_EXP
#  define R_DIG FLT_DIG
#  define CREAL crealf
#  define CIMAG cimagf
#else
#  define CEIL ceil
#  define FLOOR floor
#  define ROUND round
#  define FABS fabs
#  define FMAX fmax
#  define FMIN fmin
#  define FREXP frexp
#  define LDEXP ldexp
#  define EPS DBL_EPSILON
#  define R_MIN DBL_MIN
#  define R_MAX DBL_MAX
#  define R_MIN_EXP DBL_MIN_EXP
#  define R_MAX_EXP DBL_MAX_EXP
#  define R_MAX_10_EXP DBL_MAX_10_EXP
#  define R_MIN_10_EXP DBL_MIN_10_EXP
#  define R_DIG DBL_DIG
#  define CREAL creal
#  define CIMAG cimag
#endif

#if HAVE_DECL_DRAND48 == 0
  extern double drand48(void);
#endif
#if HAVE_DECL_SRAND48 == 0
  extern void srand48(long int);
#endif
#define RAND ((R)drand48())
#define R_RADIX FLT_RADIX
#define II _Complex_I

/* format strings */
#if defined(NFFT_LDOUBLE)
#  define FE "LE"
#elif defined(NFFT_SINGLE)
#  define FE "E"
#else
#  define FE "lE"
#endif

#define TRUE 1
#define FALSE 0

#define KPI K(3.1415926535897932384626433832795028841971693993751)
#define K2PI K(6.2831853071795864769252867665590057683943387987502)

/** Dummy use of unused parameters to silence compiler warnings */
#define UNUSED(x) (void)x

#if defined(HAVE_LIBDISPATCH)
#define FOR(VAR,VAL) \
  dispatch_apply((size_t)VAL, dispatch_get_global_queue((long)DISPATCH_QUEUE_PRIORITY_HIGH, (unsigned long)0), ^(size_t VAR)
#else
#define FOR(VAR,VAL) \
  { \
    int VAR; \
    for (VAR = 0; VAR < VAL; VAR++)
#endif

#if defined(HAVE_LIBDISPATCH)
#define END_FOR );
#else
#define END_FOR }
#endif

extern void nfft_assertion_failed(const char *s, int line, const char *file);

/* always check */
#define CK(ex) \
  (void)((ex) || (nfft_assertion_failed(#ex, __LINE__, __FILE__), 0))

#ifdef NFFT_DEBUG
  /* check only if debug enabled */
  #define A(ex) \
    (void)((ex) || (nfft_assertion_failed(#ex, __LINE__, __FILE__), 0))
#else
  #define A(ex) /* nothing */
#endif

#ifdef HAVE_ALLOCA
  /* Use alloca if available. */
  #ifndef alloca
    #ifdef __GNUC__
      /* No alloca defined but can use GCC's builtin version. */
      #define alloca __builtin_alloca
    #else
      /* No alloca defined and not using GCC. */
      #ifdef _MSC_VER
        /* Using Microsoft's C compiler. Include header file and use _alloca
         * defined therein. */
        #include <malloc.h>
        #define alloca _alloca
      #else
        /* Also not using Microsoft's C compiler. */
        #if HAVE_ALLOCA_H
          /* Alloca header is available. */
          #include <alloca.h>
        #else
          /* No alloca header available. */
          #ifdef _AIX
            /* We're using the AIX C compiler. Use pragma. */
            #pragma alloca
          #else
            /* Not using AIX compiler. */
            #ifndef alloca /* HP's cc +Olibcalls predefines alloca. */
              void *alloca(size_t);
            #endif
          #endif
        #endif
      #endif
    #endif
  #endif
  /* So we have alloca. */
  #define STACK_MALLOC(T, p, x) p = (T)alloca(x)
  #define STACK_FREE(x) /* Nothing. Cleanup done automatically. */
#else /* ! HAVE_ALLOCA */
  /* Use malloc instead of alloca. So we allocate memory on the heap instead of
   * on the stack which is slower. */
  #define STACK_MALLOC(T, p, x) p = (T)nfft_malloc(x)
  #define STACK_FREE(x) nfft_free(x)
#endif /* ! HAVE_ALLOCA */

#endif
