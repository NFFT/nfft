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

/*! \file infft.h
 *  \brief Internal header file for auxiliary definitions and functions.
 *
 *  This header builds on nfft3util.h (which provides the precision-agnostic
 *  types, name-mangling, math macros, constants and shared util prototypes)
 *  and adds the machinery that is internal to the library object itself:
 *  window-function macros, plan-coupled timing macros, stack allocation,
 *  internal linear-algebra kernels and the assertion helpers.
 *
 *  Library sources (kernel/, solver, ...) include this header. Programs that
 *  are not part of the library object (examples/, applications/, benchmarks/)
 *  should include nfft3util.h instead.
 */
#ifndef __INFFT_H__
#define __INFFT_H__

/* Shared, precision-agnostic helpers (types, macros, math, util prototypes). */
#include "nfft3util.h"

/**
 * @addtogroup nfftutil
 * @{
 */

/* Window function macros (coupled to the transform plan structs). */
/* macros for window functions */

#if defined(DIRAC_DELTA)
  #define PHI_HUT(n,k,d) K(1.0)
  #define PHI(n,x,d) IF(FABS((x)) < K(10E-8),K(1.0),K(0.0))
  #define WINDOW_HELP_INIT(d)
  #define WINDOW_HELP_FINALIZE
  #define WINDOW_HELP_ESTIMATE_m 0
#elif defined(GAUSSIAN)
  #define PHI_HUT(n,k,d) ((R)EXP(-(POW(KPI*(k)/n,K(2.0))*ths->b[d])))
  #define PHI(n,x,d) ((R)EXP(-POW((x)*((R)n),K(2.0)) / \
    ths->b[d])/SQRT(KPI*ths->b[d]))
  #define WINDOW_HELP_INIT \
    { \
      int WINDOW_idx; \
      ths->b = (R*) Y(malloc)(ths->d*sizeof(R)); \
      for (WINDOW_idx = 0; WINDOW_idx < ths->d; WINDOW_idx++) \
        ths->b[WINDOW_idx]=(K(2.0)*ths->sigma[WINDOW_idx]) / \
          (K(2.0)*ths->sigma[WINDOW_idx] - K(1.0)) * (((R)ths->m) / KPI); \
    }
  #define WINDOW_HELP_FINALIZE {Y(free)(ths->b);}
  #if MANT_DIG == 113
    // IEEE 754 quadruple precision, 128 bits.
    // TODO: Set good value for quadruple precision.
    #define WINDOW_HELP_ESTIMATE_m 17
  #elif MANT_DIG == 64
    // Intel double extended, 80 bits.
    #define WINDOW_HELP_ESTIMATE_m 17
  #elif MANT_DIG == 53
    // IEEE 754 double precision, 64 bits.
    #define WINDOW_HELP_ESTIMATE_m 13
  #elif MANT_DIG == 24
    // IEEE 754 single precision, 32 bits.
    #define WINDOW_HELP_ESTIMATE_m 5
  #else
    // Unknown floating-point type.
    // Assume IEEE 754 double precision, 64 bits.
    #define WINDOW_HELP_ESTIMATE_m 13
  #endif
#elif defined(B_SPLINE)
  #define PHI_HUT(n,k,d) ((R)(((k) == 0) ? K(1.0) / n : \
    POW(SIN((k) * KPI / n) / ((k) * KPI / n), \
      K(2.0) * ths->m)/n))
  #define PHI(n,x,d) (Y(bsplines)(2*ths->m,((x)*n) + \
    (R)ths->m) / n)
  #define WINDOW_HELP_INIT
  #define WINDOW_HELP_FINALIZE
  #if MANT_DIG == 113
    // IEEE 754 quadruple precision, 128 bits.
    // TODO: Set good value for quadruple precision.
    #define WINDOW_HELP_ESTIMATE_m 11
  #elif MANT_DIG == 64
    // Intel double extended, 80 bits.
    #define WINDOW_HELP_ESTIMATE_m 11
  #elif MANT_DIG == 53
    // IEEE 754 double precision, 64 bits.
    #define WINDOW_HELP_ESTIMATE_m 11
  #elif MANT_DIG == 24
    // IEEE 754 single precision, 32 bits.
    #define WINDOW_HELP_ESTIMATE_m 11
  #else
    // Unknown floating-point type.
    // Assume IEEE 754 double precision, 64 bits.
    #define WINDOW_HELP_ESTIMATE_m 11
  #endif
#elif defined(SINC_POWER)
  #define PHI_HUT(n,k,d) (Y(bsplines)(2 * ths->m, (K(2.0) * ths->m*(k)) / \
    ((K(2.0) * ths->sigma[(d)] - 1) * n / \
      ths->sigma[(d)]) + (R)ths->m))
  #define PHI(n,x,d) ((R)(n / ths->sigma[(d)] * \
    (K(2.0) * ths->sigma[(d)] - K(1.0))/ (K(2.0)*ths->m) * \
    POW(Y(sinc)(KPI * n / ths->sigma[(d)] * (x) * \
    (K(2.0) * ths->sigma[(d)] - K(1.0)) / (K(2.0)*ths->m)) , 2*ths->m) / \
    n))
  #define WINDOW_HELP_INIT
  #define WINDOW_HELP_FINALIZE
  #if MANT_DIG == 113
    // IEEE 754 quadruple precision, 128 bits.
    // TODO: Set good value for quadruple precision.
    #define WINDOW_HELP_ESTIMATE_m 13
  #elif MANT_DIG == 64
    // Intel double extended, 80 bits.
    #define WINDOW_HELP_ESTIMATE_m 13
  #elif MANT_DIG == 53
    // IEEE 754 double precision, 64 bits.
    #define WINDOW_HELP_ESTIMATE_m 11
  #elif MANT_DIG == 24
    // IEEE 754 single precision, 32 bits.
    #define WINDOW_HELP_ESTIMATE_m 11
  #else
    // Unknown floating-point type.
    // Assume IEEE 754 double precision, 64 bits.
    #define WINDOW_HELP_ESTIMATE_m 11
  #endif
#else /* Kaiser-Bessel is the default. */
  #define PHI_HUT(n,k,d) (Y(bessel_i0)((R)(ths->m) * SQRT(ths->b[d] * ths->b[d] - (K(2.0) * KPI * (R)(k) / (R)(n)) * (K(2.0) * KPI * (R)(k) / (R)(n)))))
  #define PHI(n,x,d) (  (((R)(ths->m) * (R)(ths->m) - (x) * (R)(n) * (x) * (R)(n)) > K(0.0)) \
                      ?   SINH(ths->b[d] * SQRT((R)(ths->m) * (R)(ths->m) - (x) * (R)(n) * (x) * (R)(n))) \
                        / (KPI * SQRT((R)(ths->m) * (R)(ths->m) - (x) * (R)(n) * (x) * (R)(n))) \
                      :   ((((R)(ths->m) * (R)(ths->m) - (x) * (R)(n) * (x) * (R)(n)) < K(0.0)) \
                        ?   SIN(ths->b[d] * SQRT((x) * (R)(n) * (x) * (R)(n) - (R)(ths->m) * (R)(ths->m))) \
                          / (KPI * SQRT((x) * (R)(n) * (x) * (R)(n) - (R)(ths->m) * (R)(ths->m))) \
                        : ths->b[d] / KPI))
  #define WINDOW_HELP_INIT \
    { \
      int WINDOW_idx; \
      ths->b = (R*) Y(malloc)((size_t)(ths->d) * sizeof(R)); \
      for (WINDOW_idx = 0; WINDOW_idx < ths->d; WINDOW_idx++) \
        ths->b[WINDOW_idx] = (KPI * (K(2.0) - K(1.0) / ths->sigma[WINDOW_idx])); \
  }
  #define WINDOW_HELP_FINALIZE {Y(free)(ths->b);}
  #if MANT_DIG == 113
    // IEEE 754 quadruple precision, 128 bits.
    // TODO: Set good value for quadruple precision.
    #define WINDOW_HELP_ESTIMATE_m 10 
  #elif MANT_DIG == 64
    // Intel double extended, 80 bits.
    #define WINDOW_HELP_ESTIMATE_m 9
  #elif MANT_DIG == 53
    // IEEE 754 double precision, 64 bits.
    #define WINDOW_HELP_ESTIMATE_m 8
  #elif MANT_DIG == 24
      #define WINDOW_HELP_ESTIMATE_m 4
  #else
    // Unknown floating-point type.
    // Assume IEEE 754 double precision, 64 bits.
    #define WINDOW_HELP_ESTIMATE_m 8
  #endif
#endif

/* window.c */
INT Y(m2K)(const INT m);

/* Stack allocation helpers. */
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
  #define STACK_MALLOC(T, p, x) p = (T)Y(malloc)(x)
  #define STACK_FREE(x) Y(free)(x)
#endif /* ! HAVE_ALLOCA */

/* Plan-coupled timing macros (require MEASURE_TIME and a plan pointer). */
/** Timing, method works since the inaccurate timer is updated mostly in the
 *  measured function. For small times not every call of the measured function
 *  will also produce a 'unit' time step.
 *  Measuring the fftw might cause a wrong output vector due to the repeated
 *  ffts.
 */
#ifdef MEASURE_TIME
 int MEASURE_TIME_r;
 double MEASURE_TIME_tt;
 ticks MEASURE_TIME_t0, MEASURE_TIME_t1;

#define TIC(a)                                                                \
  ths->MEASURE_TIME_t[(a)]=0;                                                 \
  MEASURE_TIME_r=0;                                                           \
  /* DISABLED LOOP due to code blocks causing segfault when repeatedly run */ \
  /*while(ths->MEASURE_TIME_t[(a)]<0.01)*/                                    \
    {                                                                         \
      MEASURE_TIME_r++;                                                       \
      MEASURE_TIME_t0 = getticks();                                           \

/* THE MEASURED FUNCTION IS CALLED REPEATEDLY */

#define TOC(a)                                                                \
      MEASURE_TIME_t1 = getticks();                                           \
      MEASURE_TIME_tt = Y(elapsed_seconds)(MEASURE_TIME_t1,MEASURE_TIME_t0);\
      ths->MEASURE_TIME_t[(a)]+=MEASURE_TIME_tt;                              \
    }                                                                         \
  ths->MEASURE_TIME_t[(a)]/=MEASURE_TIME_r;                                   \

#else
#define TIC(a)
#define TOC(a)
#endif

#ifdef MEASURE_TIME_FFTW
#define TIC_FFTW(a) TIC(a)
#define TOC_FFTW(a) TOC(a)
#else
#define TIC_FFTW(a)
#define TOC_FFTW(a)
#endif

/* Internal-only util prototypes (used by the library implementation). */
/* sort.c: */
void Y(sort_node_indices_radix_msdf)(INT n, INT *keys0, INT *keys1, INT rhigh);
void Y(sort_node_indices_radix_lsdf)(INT n, INT *keys0, INT *keys1, INT rhigh);

/* assert.c */
void Y(assertion_failed)(const char *s, int line, const char *file);

/* vector1.c */
/** Computes the inner/dot product \f$x^H x\f$. */
R Y(dot_double)(R *x, INT n);
/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$. */
R Y(dot_w_complex)(C *x, R *w, INT n);
/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$. */
R Y(dot_w_double)(R *x, R *w, INT n);
/** Computes the weighted inner/dot product \f$x^H (w\odot w2\odot w2 \odot x)\f$. */
R Y(dot_w_w2_complex)(C *x, R *w, R *w2, INT n);
/** Computes the weighted inner/dot product \f$x^H (w2\odot w2 \odot x)\f$. */
R Y(dot_w2_complex)(C *x, R *w2, INT n);

/* vector2.c */
/** Copies \f$x \leftarrow y\f$. */
void Y(cp_complex)(C *x, C *y, INT n);
/** Copies \f$x \leftarrow y\f$. */
void Y(cp_double)(R *x, R *y, INT n);
/** Copies \f$x \leftarrow a y\f$. */
void Y(cp_a_complex)(C *x, R a, C *y, INT n);
/** Copies \f$x \leftarrow a y\f$. */
void Y(cp_a_double)(R *x, R a, R *y, INT n);
/** Copies \f$x \leftarrow w\odot y\f$. */
void Y(cp_w_complex)(C *x, R *w, C *y, INT n);
/** Copies \f$x \leftarrow w\odot y\f$. */
void Y(cp_w_double)(R *x, R *w, R *y, INT n);

/* vector3.c */
/** Updates \f$x \leftarrow a x + y\f$. */
void Y(upd_axpy_double)(R *x, R a, R *y, INT n);
/** Updates \f$x \leftarrow x + a y\f$. */
void Y(upd_xpay_complex)(C *x, R a, C *y, INT n);
/** Updates \f$x \leftarrow x + a y\f$. */
void Y(upd_xpay_double)(R *x, R a, R *y, INT n);
/** Updates \f$x \leftarrow a x + b y\f$. */
void Y(upd_axpby_complex)(C *x, R a, C *y, R b, INT n);
/** Updates \f$x \leftarrow a x + b y\f$. */
void Y(upd_axpby_double)(R *x, R a, R *y, R b, INT n);
/** Updates \f$x \leftarrow x + a w\odot y\f$. */
void Y(upd_xpawy_complex)(C *x, R a, R *w, C *y, INT n);
/** Updates \f$x \leftarrow x + a w\odot y\f$. */
void Y(upd_xpawy_double)(R *x, R a, R *w, R *y, INT n);
/** Updates \f$x \leftarrow a x +  w\odot y\f$. */
void Y(upd_axpwy_complex)(C *x, R a, R *w, C *y, INT n);
/** Updates \f$x \leftarrow a x +  w\odot y\f$. */
void Y(upd_axpwy_double)(R *x, R a, R *w, R *y, INT n);

/* always check */
#define CK(ex) \
  (void)((ex) || (Y(assertion_failed)(#ex, __LINE__, __FILE__), 0))

#ifdef NFFT_DEBUG
  /* check only if debug enabled */
  #define A(ex) \
    (void)((ex) || (Y(assertion_failed)(#ex, __LINE__, __FILE__), 0))
#else
  #define A(ex) /* nothing */
#endif

/** @}
 */
#endif /* __INFFT_H__ */
