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
 */
#ifndef __INFFT_H__
#define __INFFT_H__

#include "config.h"

#include <math.h>
#include <float.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
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

#include <fftw3.h>

#include "ticks.h"

/**
 * @defgroup nfftutil Util - Auxiliary functions
 * @{
 *
 * This module implements frequently used utility functions.
 * In particular, this includes simple measurement of resources, evaluation of
 * window functions, vector routines for basic linear algebra tasks, and
 * computation of weights for the inverse transforms.
 *
 */

/* Determine precision and name-mangling scheme. */
#define CONCAT(prefix, name) prefix ## name
#if defined(NFFT_SINGLE)
typedef float R;
typedef float _Complex C;
#define Y(name) CONCAT(nfftf_,name)
#define FFTW(name) CONCAT(fftwf_,name)
#define NFFT(name) CONCAT(nfftf_,name)
#define NFCT(name) CONCAT(nfctf_,name)
#define NFST(name) CONCAT(nfstf_,name)
#define NFSFT(name) CONCAT(nfsftf_,name)
#define SOLVER(name) CONCAT(solverf_,name)
#elif defined(NFFT_LDOUBLE)
typedef long double R;
typedef long double _Complex C;
#define Y(name) CONCAT(nfftl_,name)
#define FFTW(name) CONCAT(fftwl_,name)
#define NFFT(name) CONCAT(nfftl_,name)
#define NFCT(name) CONCAT(nfctl_,name)
#define NFST(name) CONCAT(nfstl_,name)
#define NFSFT(name) CONCAT(nfsftl_,name)
#define SOLVER(name) CONCAT(solverl_,name)
#else
typedef double R;
typedef double _Complex C;
#define Y(name) CONCAT(nfft_,name)
#define FFTW(name) CONCAT(fftw_,name)
#define NFFT(name) CONCAT(nfft_,name)
#define NFCT(name) CONCAT(nfct_,name)
#define NFST(name) CONCAT(nfst_,name)
#define NFSFT(name) CONCAT(nfsft_,name)
#define SOLVER(name) CONCAT(solver_,name)
#endif
#define X(name) Y(name)

#define STRINGIZEx(x) #x
#define STRINGIZE(x) STRINGIZEx(x)

#ifdef NFFT_LDOUBLE
#  define K(x) ((R) x##L)
#else
#  define K(x) ((R) x)
#endif
#define DK(name, value) const R name = K(value)

#if defined __CYGWIN32__ && !defined __CYGWIN__
   /* For backwards compatibility with Cygwin b19 and
      earlier, we define __CYGWIN__ here, so that
      we can rely on checking just for that macro. */
#  define __CYGWIN__  __CYGWIN32__
#endif

/* Integral type large enough to contain a stride (what ``int'' should have been
 * in the first place) */
typedef ptrdiff_t INT;

#if defined(NFFT_LDOUBLE)
  #define MANT_DIG LDBL_MANT_DIG
  #define MIN_EXP LDBL_MIN_EXP
  #define MAX_EXP LDBL_MAX_EXP
  #define EPSILON LDBL_EPSILON
#elif defined(NFFT_SINGLE)
  #define MANT_DIG FLT_MANT_DIG
  #define MIN_EXP FLT_MIN_EXP
  #define MAX_EXP FLT_MAX_EXP
  #define EPSILON FLT_EPSILON
#else
  #define MANT_DIG DBL_MANT_DIG
  #define MIN_EXP DBL_MIN_EXP
  #define MAX_EXP DBL_MAX_EXP
  #define EPSILON DBL_EPSILON
#endif

#define KPI K(3.1415926535897932384626433832795028841971693993751)
#define K2PI K(6.2831853071795864769252867665590057683943387987502)
#define K4PI K(12.5663706143591729538505735331180115367886775975004)
#define KE K(2.7182818284590452353602874713526624977572470937000)

#define IF(x,a,b) ((x)?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define ABS(x) (((x)>K(0.0))?(x):(-(x)))
#define SIGN(a) (((a)>=0)?1:-1)
#define SIGN(a) (((a)>=0)?1:-1)
#define SIGNF(a) IF((a)<K(0.0),K(-1.0),K(1.0))

/* Size of array. */
#define SIZE(x) sizeof(x)/sizeof(x[0])

/* How far Clenshaw on a Chebyshev piece of B_k may lose relative accuracy
 * before the evaluation falls back to de Boor. */
#define NFFT_CHEB_GUARD K(256.0)

/** Swap two vectors. */
#define CSWAP(x,y) {C* NFFT_SWAP_temp__; \
  NFFT_SWAP_temp__=(x); (x)=(y); (y)=NFFT_SWAP_temp__;}

/** Swap two vectors. */
#define RSWAP(x,y) {R* NFFT_SWAP_temp__; NFFT_SWAP_temp__=(x); \
  (x)=(y); (y)=NFFT_SWAP_temp__;}

/* macros for window functions */

/* Half-width, in grid spacings, that the 2m+2 point run reaches: the distance
 * to the nearest point uo() leaves out. The default is the floor(n x)
 * centring, whose offset lies in [0,1); a module whose uo() centres the run
 * differently defines this before including. Only the Gaussian reads it, the
 * other windows being zero past |x| <= m/n. */
#ifndef WINDOW_STENCIL_REACH
  #define WINDOW_STENCIL_REACH (((R)ths->m) + K(1.0))
#endif

#if defined(DIRAC_DELTA)
  #define PHI_HUT(n,k,d) K(1.0)
  #define PHI(n,x,d) IF(FABS((x)) < K(10E-8),K(1.0),K(0.0))
  #define WINDOW_HELP_INIT(d)
  #define WINDOW_HELP_FINALIZE
  #define WINDOW_HELP_ESTIMATE_m 0
#elif defined(GAUSSIAN)
  /* ths->b holds 3 * ths->d reals: b per axis, then 1/b and 1/sqrt(pi b). b
   * keeps index ax, where the FG_PSI paths read it. The other two are what an
   * evaluation would otherwise recompute per point: with them, phi is a
   * multiply, an exponential and a multiply. */
  #define GA_B(ax) (ths->b[(ax)])
  #define GA_B_INV(ax) (ths->b[(ths->d) + (ax)])
  #define GA_NORM(ax) (ths->b[2 * (ths->d) + (ax)])
  #define PHI_HUT(n,k,ax) (Y(gaussian_phi_hut)(GA_B(ax), (R)(k) * KPI / (R)(n)))
  #define PHI(n,x,ax) (Y(gaussian_phi)(GA_B_INV(ax), GA_NORM(ax), \
                          (R)(n) * (x)))
  #define PHI_RUN(dst,n,x,u,ax) \
    Y(gaussian_phi_run)((dst), GA_B_INV(ax), GA_NORM(ax), (ths->m), \
        (R)(n) * (x) - (R)(u))
  /* b balances the two error terms (frequency and time/space) of a Gaussian
   * truncated at half-width WINDOW_STENCIL_REACH. */
  #define WINDOW_HELP_INIT \
    { \
      int WINDOW_idx; \
      ths->b = (R*) Y(malloc)((size_t)(3 * ths->d) * sizeof(R)); \
      for (WINDOW_idx = 0; WINDOW_idx < ths->d; WINDOW_idx++) \
      { \
        const R WINDOW_b = (K(2.0)*ths->sigma[WINDOW_idx]) / \
          (K(2.0)*ths->sigma[WINDOW_idx] - K(1.0)) * \
          (WINDOW_STENCIL_REACH / KPI); \
        ths->b[WINDOW_idx] = WINDOW_b; \
        ths->b[ths->d + WINDOW_idx] = K(1.0) / WINDOW_b; \
        ths->b[2 * ths->d + WINDOW_idx] = K(1.0) / SQRT(KPI * WINDOW_b); \
      } \
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
  /* ths->b holds 2m reals of 1/(j+1) (PHI_RUN's de Boor triangular sweep),
   * then the (2m)x(2m) Chebyshev table Y(bspline_cheb_init) builds (PHI's
   * single-point evaluator) -- the same table-building function the
   * sinc-power window's phi_hut uses. PHI is a rarer single-point call (fg_psi
   * table builds, precompute_lin_psi) where the O(m) Chebyshev evaluator
   * beats de Boor's O(m^2) outright. PHI_RUN already amortizes to O(m) per
   * point via one shared triangular sweep over the whole run -- an
   * independent per-point Chebyshev evaluation was measured slower there, so
   * it keeps de Boor. */
  #define BS_CHEB (ths->b + 2 * ths->m)
  #define PHI_HUT(n,k,ax) (Y(bspline_phi_hut)((R)ths->m, (R)(n), (R)(k)))
  #define PHI(n,x,ax) (Y(bspline_cheb)(BS_CHEB, 2 * ths->m, \
                          (x) * (R)(n) + (R)ths->m) / (R)(n))
  #define PHI_RUN(dst,n,x,u,ax) \
    Y(bspline_phi_run)((dst), ths->b, (ths->m), \
        (R)(n) * (x) - (R)(u) + (R)ths->m, K(1.0) / (R)(n))
  #define WINDOW_HELP_INIT \
    { \
      int WINDOW_idx; \
      const INT WINDOW_k = 2 * ths->m; \
      ths->b = (R*) Y(malloc)((size_t)(WINDOW_k + WINDOW_k * WINDOW_k) \
          * sizeof(R)); \
      for (WINDOW_idx = 0; WINDOW_idx < WINDOW_k; WINDOW_idx++) \
        ths->b[WINDOW_idx] = K(1.0) / (R)(WINDOW_idx + 1); \
      Y(bspline_cheb_init)(BS_CHEB, WINDOW_k); \
    }
  #define WINDOW_HELP_FINALIZE {Y(free)(ths->b);}
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
  /* phi(x) = w sinc(pi n w x)^(2m) and phi_hut(k) = B_{2m}(k/(w n) + m), with
   * w = (2 sigma - 1)/(2 m sigma). n cancels out of both prefactors, which it
   * must: NFCT and NFST pass a doubled n to the same macros. ths->b holds 2d
   * reals, w per axis then 1/w. */
  #define SP_W(ax) (ths->b[(ax)])
  #define SP_W_INV(ax) (ths->b[(ths->d) + (ax)])
  #define SP_CHEB (ths->b + 2 * (ths->d))
  #define SP_CHEB_I0 ((INT)(ths->b[2 * (ths->d) + 4 * (ths->m) * (ths->m)]))
  /* Deconvolution divides by phi_hut, so this one needs relative accuracy and
   * takes the guarded evaluation. */
  #define PHI_HUT(n,k,ax) (Y(bspline_cheb_rel)(SP_CHEB, 2 * ths->m, \
                              SP_CHEB_I0, \
                              (R)(k) * SP_W_INV(ax) / (R)(n) + (R)ths->m))
  #define PHI(n,x,ax) (Y(sincpow_phi)(SP_W(ax), (R)ths->m, \
                           KPI * (R)(n) * SP_W(ax) * (x)))
  /* The sinc arguments across a run are one step of KPI * w apart. */
  #define PHI_RUN(dst,n,x,u,ax) \
    Y(sincpow_phi_run)((dst), SP_W(ax), (R)ths->m, (ths->m), \
        KPI * SP_W(ax), (R)(n) * (x) - (R)(u))
  #define WINDOW_HELP_INIT \
    { \
      int WINDOW_idx; \
      const INT WINDOW_k = 2 * ths->m; \
      ths->b = (R*) Y(malloc)((size_t)(2 * ths->d + WINDOW_k * WINDOW_k + 1) \
          * sizeof(R)); \
      for (WINDOW_idx = 0; WINDOW_idx < ths->d; WINDOW_idx++) \
      { \
        const R WINDOW_s = ths->sigma[WINDOW_idx]; \
        const R WINDOW_w = (K(2.0) * WINDOW_s - K(1.0)) \
            / (K(2.0) * ((R)ths->m) * WINDOW_s); \
        ths->b[WINDOW_idx] = WINDOW_w; \
        ths->b[ths->d + WINDOW_idx] = K(1.0) / WINDOW_w; \
      } \
      Y(bspline_cheb_init)(ths->b + 2 * ths->d, WINDOW_k); \
      ths->b[2 * ths->d + WINDOW_k * WINDOW_k] = \
          (R)Y(bspline_cheb_guard)(ths->b + 2 * ths->d, WINDOW_k, \
              NFFT_CHEB_GUARD); \
    }
  #define WINDOW_HELP_FINALIZE {Y(free)(ths->b);}
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
  /* PHI and PHI_HUT both carry the factor exp(-log I0(m b)) so that PHI_HUT(n,0,ax) 
   * is normalized to 1. Deconvolution divides by PHI_HUT and convolution multiplies 
   * by PHI, so the factor cancels and the transform is unchanged.
   *
   * ths->b holds 4 * ths->d reals: b per axis, log I0(m b) - m b,
   * 1 / (exp(-m b) * I0(m b)) and 1 / I0(m b). The last two are the peak in the
   * two forms the window needs it in; both are constant per axis, so neither
   * costs a Bessel evaluation inside a loop. The second follows from the first,
   * since log I0(m b) is m b + lg_tail. */
  #define KB_B(ax) (ths->b[(ax)])
  #define KB_LG_TAIL(ax) (ths->b[(ths->d) + (ax)])
  #define KB_I0E_PEAK_INV(ax) (ths->b[2 * (ths->d) + (ax)])
  #define KB_PEAK_INV(ax) (ths->b[3 * (ths->d) + (ax)])
  #define PHI_HUT(n,k,ax) (Y(kb_phi_hut)(KB_B(ax), KB_I0E_PEAK_INV(ax), \
                             (R)(ths->m), (R)(n), (R)(k)))
  #define PHI(n,x,ax) (Y(kb_phi)(KB_B(ax), KB_LG_TAIL(ax), KB_PEAK_INV(ax), \
                         (R)(ths->m), (R)(n) * (R)(x)))
  /* Fills dst[0 .. 2m+1] with phi(x - (u + l)/n), the run every psi table is
   * built from. The window sees the whole run, so it can hoist its constants,
   * step the argument by one grid cell rather than dividing per point, and put
   * the points that need the guarded evaluation in their own branch. */
  #define PHI_RUN(dst,n,x,u,ax) \
    Y(kb_phi_run)((dst), KB_B(ax), KB_LG_TAIL(ax), KB_PEAK_INV(ax), \
        (R)(ths->m), (ths->m), (R)(n) * (R)(x) - (R)(u))
  #define WINDOW_HELP_INIT \
    { \
      int WINDOW_idx; \
      ths->b = (R*) Y(malloc)((size_t)(4 * ths->d) * sizeof(R)); \
      for (WINDOW_idx = 0; WINDOW_idx < ths->d; WINDOW_idx++) \
      { \
        R WINDOW_b = (KPI * (K(2.0) - K(1.0) / ths->sigma[WINDOW_idx])); \
        R WINDOW_xpk = ((R)(ths->m)) * WINDOW_b; \
        R WINDOW_lg = Y(bessel_i0_logtail)(WINDOW_xpk); \
        ths->b[WINDOW_idx] = WINDOW_b; \
        ths->b[ths->d + WINDOW_idx] = WINDOW_lg; \
        ths->b[2 * ths->d + WINDOW_idx] = K(1.0) / \
            Y(bessel_i0_exp_scaled)(WINDOW_xpk); \
        ths->b[3 * ths->d + WINDOW_idx] = EXP(-WINDOW_xpk - WINDOW_lg); \
      } \
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

/* Generic run fill: dst[0 .. 2m+1] = phi(x - (u + l)/n). Windows that gain from
 * seeing the whole run at once override PHI_RUN above. */
#ifndef PHI_RUN
  #define PHI_RUN(dst,n,x,u,ax) \
    do { \
      INT PHI_RUN_l; \
      for (PHI_RUN_l = 0; PHI_RUN_l < 2 * ths->m + 2; PHI_RUN_l++) \
        (dst)[PHI_RUN_l] = PHI((n), \
            (x) - ((R)(PHI_RUN_l + (u))) / ((R)(n)), (ax)); \
    } while (0)
#endif

/* Gaussian run fill for the FG_PSI paths: dst[0 .. len-1], stepping outward
 * from the peak. seed is the window at the grid point nearest the node, kc
 * that point's index in the run, s is exp(2 d / b) for its offset d, e is
 * exp(-1/b) and q is e * e. Stepping outward keeps each step on a value no
 * larger than the one before, holding the run to about an ulp of its peak;
 * stepping up from the tail instead costs an ulp per grid point. */
#define FG_RUN(dst,len,seed,s,e,q,kc) \
  do { \
    const R FG_RUN_v0 = (seed), FG_RUN_s = (s), FG_RUN_e = (e), \
        FG_RUN_q = (q); \
    const INT FG_RUN_c = (kc), FG_RUN_n = (len); \
    R FG_RUN_v = FG_RUN_v0, FG_RUN_r = FG_RUN_s * FG_RUN_e; \
    INT FG_RUN_k; \
    (dst)[FG_RUN_c] = FG_RUN_v0; \
    for (FG_RUN_k = FG_RUN_c + 1; FG_RUN_k < FG_RUN_n; FG_RUN_k++) \
    { \
      FG_RUN_v *= FG_RUN_r; \
      (dst)[FG_RUN_k] = FG_RUN_v; \
      FG_RUN_r *= FG_RUN_q; \
    } \
    FG_RUN_v = FG_RUN_v0; \
    FG_RUN_r = FG_RUN_e / FG_RUN_s; \
    for (FG_RUN_k = FG_RUN_c - 1; FG_RUN_k >= 0; FG_RUN_k--) \
    { \
      FG_RUN_v *= FG_RUN_r; \
      (dst)[FG_RUN_k] = FG_RUN_v; \
      FG_RUN_r *= FG_RUN_q; \
    } \
  } while (0)

/* window.c */
INT Y(m2K)(const INT m);

#if defined(NFFT_LDOUBLE)
#if HAVE_DECL_COPYSIGNL == 0
extern long double copysignl(long double, long double);
#endif
#if HAVE_DECL_NEXTAFTERL == 0
extern long double nextafterl(long double, long double);
#endif
#if HAVE_DECL_NANL == 0
extern long double nanl(const char *tag);
#endif
#if HAVE_DECL_CEILL == 0
extern long double ceill(long double);
#endif
#if HAVE_DECL_FLOORL == 0
extern long double floorl(long double);
#endif
#if HAVE_DECL_NEARBYINTL == 0
extern long double nearbyintl(long double);
#endif
#if HAVE_DECL_RINTL == 0
extern long double rintl(long double);
#endif
#if HAVE_DECL_ROUNDL == 0
extern long double roundl(long double);
#endif
#if HAVE_DECL_LRINTL == 0
extern long int lrintl(long double);
#endif
#if HAVE_DECL_LROUNDL == 0
extern long int lroundl(long double);
#endif
#if HAVE_DECL_LLRINTL == 0
extern long long int llrintl(long double);
#endif
#if HAVE_DECL_LLROUNDL == 0
extern long long int llroundl(long double);
#endif
#if HAVE_DECL_TRUNCL == 0
extern long double truncl(long double);
#endif
#if HAVE_DECL_FMODL == 0
extern long double fmodl(long double, long double);
#endif
#if HAVE_DECL_REMAINDERL == 0
extern long double remainderl(long double, long double);
#endif
#if HAVE_DECL_REMQUOL == 0
extern long double remquol(long double x, long double y, int *);
#endif
#if HAVE_DECL_FDIML == 0
extern long double fdiml(long double, long double);
#endif
#if HAVE_DECL_FMAXL == 0
extern long double fmaxl(long double, long double);
#endif
#if HAVE_DECL_FMINL == 0
extern long double fminl(long double, long double);
#endif
#if HAVE_DECL_FMAL == 0
extern long double fmal(long double x, long double y, long double z);
#endif
#if HAVE_DECL_FABSL == 0
extern long double fabsl(long double);
#endif
#if HAVE_DECL_SQRTL == 0
extern long double sqrtl(long double);
#endif
#if HAVE_DECL_CBRTL == 0
extern long double cbrtl(long double);
#endif
#if HAVE_DECL_HYPOTL == 0
extern long double hypotl(long double, long double);
#endif
#if HAVE_DECL_EXPL == 0
extern long double expl(long double);
#endif
#if HAVE_DECL_EXP2L == 0
extern long double exp2l(long double);
#endif
#if HAVE_DECL_EXPM1L == 0
extern long double expm1l(long double);
#endif
#if HAVE_DECL_LOGL == 0
extern long double logl(long double);
#endif
#if HAVE_DECL_LOG2L == 0
extern long double log2l(long double);
#endif
#if HAVE_DECL_LOG10L == 0
extern long double log10l(long double);
#endif
#if HAVE_DECL_LOG1PL == 0
extern long double log1pl(long double);
#endif
#if HAVE_DECL_LOGBL == 0
extern long double logbl(long double);
#endif
#if HAVE_DECL_ILOGBL == 0
extern int ilogbl(long double);
#endif
#if HAVE_DECL_MODFL == 0
extern long double modfl(long double, long double *);
#endif
#if HAVE_DECL_FREXPL == 0
extern long double frexpl(long double, int *);
#endif
#if HAVE_DECL_LDEXPL == 0
extern long double ldexpl(long double, int);
#endif
#if HAVE_DECL_SCALBNL == 0
extern long double scalbnl(long double, int);
#endif
#if HAVE_DECL_SCALBLNL == 0
extern long double scalblnl(long double, long int);
#endif
#if HAVE_DECL_POWL == 0
extern long double powl(long double, long double);
#endif
#if HAVE_DECL_COSL == 0
extern long double cosl(long double);
#endif
#if HAVE_DECL_SINL == 0
extern long double sinl(long double);
#endif
#if HAVE_DECL_TANL == 0
extern long double tanl(long double);
#endif
#if HAVE_DECL_COSHL == 0
extern long double coshl(long double);
#endif
#if HAVE_DECL_SINHL == 0
extern long double sinhl(long double);
#endif
#if HAVE_DECL_TANHL == 0
extern long double tanhl(long double);
#endif
#if HAVE_DECL_ACOSL == 0
extern long double acosl(long double);
#endif
#if HAVE_DECL_ASINL == 0
extern long double asinl(long double);
#endif
#if HAVE_DECL_ATANL == 0
extern long double atanl(long double);
#endif
#if HAVE_DECL_ATAN2L == 0
extern long double atan2l(long double, long double);
#endif
#if HAVE_DECL_ACOSHL == 0
extern long double acoshl(long double);
#endif
#if HAVE_DECL_ASINHL == 0
extern long double asinhl(long double);
#endif
#if HAVE_DECL_ATANHL == 0
extern long double atanhl(long double);
#endif
#if HAVE_DECL_TGAMMAL == 0
extern long double tgammal(long double);
#endif
#if HAVE_DECL_LGAMMAL == 0
extern long double lgammal(long double);
#endif
#if HAVE_DECL_J0L == 0
extern long double j0l(long double);
#endif
#if HAVE_DECL_J1L == 0
extern long double j1l(long double);
#endif
#if HAVE_DECL_JNL == 0
extern long double jnl(int, long double);
#endif
#if HAVE_DECL_Y0L == 0
extern long double y0l(long double);
#endif
#if HAVE_DECL_Y1L == 0
extern long double y1l(long double);
#endif
#if HAVE_DECL_YNL == 0
extern long double ynl(int, long double);
#endif
#if HAVE_DECL_ERFL == 0
extern long double erfl(long double);
#endif
#if HAVE_DECL_ERFCL == 0
extern long double erfcl(long double);
#endif
#if HAVE_DECL_CREALL == 0
extern long double creall(long double _Complex z);
#endif
#if HAVE_DECL_CIMAGL == 0
extern long double cimagl(long double _Complex z);
#endif
#if HAVE_DECL_CABSL == 0
extern long double cabsl(long double _Complex z);
#endif
#if HAVE_DECL_CARGL == 0
extern long double cargl(long double _Complex z);
#endif
#if HAVE_DECL_CONJL == 0
extern long double _Complex conjl(long double _Complex z);
#endif
#if HAVE_DECL_CPROJL == 0
extern long double _Complex cprojl(long double _Complex z);
#endif
#if HAVE_DECL_CSQRTL == 0
extern long double _Complex csqrtl(long double _Complex z);
#endif
#if HAVE_DECL_CEXPL == 0
extern long double _Complex cexpl(long double _Complex z);
#endif
#if HAVE_DECL_CLOGL == 0
extern long double _Complex clogl(long double _Complex z);
#endif
#if HAVE_DECL_CPOWL == 0
extern long double _Complex cpowl(long double _Complex z, long double _Complex w);
#endif
#if HAVE_DECL_CSINL == 0
extern long double _Complex csinl(long double _Complex z);
#endif
#if HAVE_DECL_CCOSL == 0
extern long double _Complex ccosl(long double _Complex z);
#endif
#if HAVE_DECL_CTANL == 0
extern long double _Complex ctanl(long double _Complex z);
#endif
#if HAVE_DECL_CASINL == 0
extern long double _Complex casinl(long double _Complex z);
#endif
#if HAVE_DECL_CACOSL == 0
extern long double _Complex cacosl(long double _Complex z);
#endif
#if HAVE_DECL_CATANL == 0
extern long double _Complex catanl(long double _Complex z);
#endif
#if HAVE_DECL_CSINHL == 0
extern long double _Complex csinhl(long double _Complex z);
#endif
#if HAVE_DECL_CCOSHL == 0
extern long double _Complex ccoshl(long double _Complex z);
#endif
#if HAVE_DECL_CTANHL == 0
extern long double _Complex ctanhl(long double _Complex z);
#endif
#if HAVE_DECL_CASINHL == 0
extern long double _Complex casinhl(long double _Complex z);
#endif
#if HAVE_DECL_CACOSHL == 0
extern long double _Complex cacoshl(long double _Complex z);
#endif
#if HAVE_DECL_CATANHL == 0
extern long double _Complex catanhl(long double _Complex z);
#endif
#define COPYSIGN copysignl
#define NEXTAFTER  nextafterl
#define MKNAN nanl
#define CEIL ceill
#define FLOOR floorl
#define NEARBYINT nearbyintl
#define RINT rintl
#define ROUND roundl
#define LRINT lrintl
#define LROUND lroundl
#define LLRINT llrintl
#define LLROUND llroundl
#define TRUNC truncl
#define FMOD fmodl
#define REMAINDER remainderl
#define REMQUO remquol
#define FDIM fdiml
#define FMAX fmaxl
#define FMIN fminl
#define FFMA fmal
#define FABS fabsl
#define SQRT sqrtl
#define CBRT cbrtl
#define HYPOT hypotl
#define EXP expl
#define EXP2 exp2l
#define EXPM1 expm1l
#define LOG logl
#define LOG2 log2l
#define LOG10 log10l
#define LOG1P log1pl
#define LOGB logbl
#define ILOGB ilogbl
#define MODF modfl
#define FREXP frexpl
#define LDEXP ldexpl
#define SCALBN scalbnl
#define SCALBLN scalblnl
#define POW powl
#define COS cosl
#define SIN sinl
#define TAN tanl
#define COSH coshl
#define SINH sinhl
#define TANH tanhl
#define ACOS acosl
#define ASIN asinl
#define ATAN atanl
#define ATAN2 atan2l
#define ACOSH acoshl
#define ASINH asinhl
#define ATANH atanhl
#define TGAMMA tgammal
#define LGAMMA lgammal
#define J0 j0l
#define J1 j1l
#define JN jnl
#define Y0 y0l
#define Y1 y1l
#define YN ynl
#define ERF erfl
#define ERFC erfcl
#define CREAL creall
#define CIMAG cimagl
#define CABS cabsl
#define CARG cargl
#define CONJ conjl
#define CPROJ cprojl
#define CSQRT csqrtl
#define CEXP cexpl
#define CLOG clogl
#define CPOW cpowl
#define CSIN csinl
#define CCOS ccosl
#define CTAN ctanl
#define CASIN casinl
#define CACOS cacosl
#define CATAN catanl
#define CSINH csinhl
#define CCOSH ccoshl
#define CTANH ctanhl
#define CASINH casinhl
#define CACOSH cacoshl
#define CATANH catanhl
#elif defined(NFFT_SINGLE)
#if HAVE_DECL_COPYSIGNF == 0
extern float copysignf(float, float);
#endif
#if HAVE_DECL_NEXTAFTERF == 0
extern float nextafterf(float, float);
#endif
#if HAVE_DECL_NANF == 0
extern float nanf(const char *tag);
#endif
#if HAVE_DECL_CEILF == 0
extern float ceilf(float);
#endif
#if HAVE_DECL_FLOORF == 0
extern float floorf(float);
#endif
#if HAVE_DECL_NEARBYINTF == 0
extern float nearbyintf(float);
#endif
#if HAVE_DECL_RINTF == 0
extern float rintf(float);
#endif
#if HAVE_DECL_ROUNDF == 0
extern float roundf(float);
#endif
#if HAVE_DECL_LRINTF == 0
extern long int lrintf(float);
#endif
#if HAVE_DECL_LROUNDF == 0
extern long int lroundf(float);
#endif
#if HAVE_DECL_LLRINTF == 0
extern long long int llrintf(float);
#endif
#if HAVE_DECL_LLROUNDF == 0
extern long long int llroundf(float);
#endif
#if HAVE_DECL_TRUNCF == 0
extern float truncf(float);
#endif
#if HAVE_DECL_FMODF == 0
extern float fmodf(float, float);
#endif
#if HAVE_DECL_REMAINDERF == 0
extern float remainderf(float, float);
#endif
#if HAVE_DECL_REMQUOF == 0
extern float remquof(float x, float y, int *);
#endif
#if HAVE_DECL_FDIMF == 0
extern float fdimf(float, float);
#endif
#if HAVE_DECL_FMAXF == 0
extern float fmaxf(float, float);
#endif
#if HAVE_DECL_FMINF == 0
extern float fminf(float, float);
#endif
#if HAVE_DECL_FMAF == 0
extern float fmaf(float x, float y, float z);
#endif
#if HAVE_DECL_FABSF == 0
extern float fabsf(float);
#endif
#if HAVE_DECL_SQRTF == 0
extern float sqrtf(float);
#endif
#if HAVE_DECL_CBRTF == 0
extern float cbrtf(float);
#endif
#if HAVE_DECL_HYPOTF == 0
extern float hypotf(float, float);
#endif
#if HAVE_DECL_EXPF == 0
extern float expf(float);
#endif
#if HAVE_DECL_EXP2F == 0
extern float exp2f(float);
#endif
#if HAVE_DECL_EXPM1F == 0
extern float expm1f(float);
#endif
#if HAVE_DECL_LOGF == 0
extern float logf(float);
#endif
#if HAVE_DECL_LOG2F == 0
extern float log2f(float);
#endif
#if HAVE_DECL_LOG10F == 0
extern float log10f(float);
#endif
#if HAVE_DECL_LOG1PF == 0
extern float log1pf(float);
#endif
#if HAVE_DECL_LOGBF == 0
extern float logbf(float);
#endif
#if HAVE_DECL_ILOGBF == 0
extern int ilogbf(float);
#endif
#if HAVE_DECL_MODFF == 0
extern float modff(float, float *);
#endif
#if HAVE_DECL_FREXPF == 0
extern float frexpf(float, int *);
#endif
#if HAVE_DECL_LDEXPF == 0
extern float ldexpf(float, int);
#endif
#if HAVE_DECL_SCALBNF == 0
extern float scalbnf(float, int);
#endif
#if HAVE_DECL_SCALBLNF == 0
extern float scalblnf(float, long int);
#endif
#if HAVE_DECL_POWF == 0
extern float powf(float, float);
#endif
#if HAVE_DECL_COSF == 0
extern float cosf(float);
#endif
#if HAVE_DECL_SINF == 0
extern float sinf(float);
#endif
#if HAVE_DECL_TANF == 0
extern float tanf(float);
#endif
#if HAVE_DECL_COSHF == 0
extern float coshf(float);
#endif
#if HAVE_DECL_SINHF == 0
extern float sinhf(float);
#endif
#if HAVE_DECL_TANHF == 0
extern float tanhf(float);
#endif
#if HAVE_DECL_ACOSF == 0
extern float acosf(float);
#endif
#if HAVE_DECL_ASINF == 0
extern float asinf(float);
#endif
#if HAVE_DECL_ATANF == 0
extern float atanf(float);
#endif
#if HAVE_DECL_ATAN2F == 0
extern float atan2f(float, float);
#endif
#if HAVE_DECL_ACOSHF == 0
extern float acoshf(float);
#endif
#if HAVE_DECL_ASINHF == 0
extern float asinhf(float);
#endif
#if HAVE_DECL_ATANHF == 0
extern float atanhf(float);
#endif
#if HAVE_DECL_TGAMMAF == 0
extern float tgammaf(float);
#endif
#if HAVE_DECL_LGAMMAF == 0
extern float lgammaf(float);
#endif
#if HAVE_DECL_J0F == 0
extern float j0f(float);
#endif
#if HAVE_DECL_J1F == 0
extern float j1f(float);
#endif
#if HAVE_DECL_JNF == 0
extern float jnf(int, float);
#endif
#if HAVE_DECL_Y0F == 0
extern float y0f(float);
#endif
#if HAVE_DECL_Y1F == 0
extern float y1f(float);
#endif
#if HAVE_DECL_YNF == 0
extern float ynf(int, float);
#endif
#if HAVE_DECL_ERFF == 0
extern float erff(float);
#endif
#if HAVE_DECL_ERFCF == 0
extern float erfcf(float);
#endif
#if HAVE_DECL_CREALF == 0
extern float crealf(float _Complex z);
#endif
#if HAVE_DECL_CIMAGF == 0
extern float cimagf(float _Complex z);
#endif
#if HAVE_DECL_CABSF == 0
extern float cabsf(float _Complex z);
#endif
#if HAVE_DECL_CARGF == 0
extern float cargf(float _Complex z);
#endif
#if HAVE_DECL_CONJF == 0
extern float _Complex conjf(float _Complex z);
#endif
#if HAVE_DECL_CPROJF == 0
extern float _Complex cprojf(float _Complex z);
#endif
#if HAVE_DECL_CSQRTF == 0
extern float _Complex csqrtf(float _Complex z);
#endif
#if HAVE_DECL_CEXPF == 0
extern float _Complex cexpf(float _Complex z);
#endif
#if HAVE_DECL_CLOGF == 0
extern float _Complex clogf(float _Complex z);
#endif
#if HAVE_DECL_CPOWF == 0
extern float _Complex cpowf(float _Complex z, float _Complex w);
#endif
#if HAVE_DECL_CSINF == 0
extern float _Complex csinf(float _Complex z);
#endif
#if HAVE_DECL_CCOSF == 0
extern float _Complex ccosf(float _Complex z);
#endif
#if HAVE_DECL_CTANF == 0
extern float _Complex ctanf(float _Complex z);
#endif
#if HAVE_DECL_CASINF == 0
extern float _Complex casinf(float _Complex z);
#endif
#if HAVE_DECL_CACOSF == 0
extern float _Complex cacosf(float _Complex z);
#endif
#if HAVE_DECL_CATANF == 0
extern float _Complex catanf(float _Complex z);
#endif
#if HAVE_DECL_CSINHF == 0
extern float _Complex csinhf(float _Complex z);
#endif
#if HAVE_DECL_CCOSHF == 0
extern float _Complex ccoshf(float _Complex z);
#endif
#if HAVE_DECL_CTANHF == 0
extern float _Complex ctanhf(float _Complex z);
#endif
#if HAVE_DECL_CASINHF == 0
extern float _Complex casinhf(float _Complex z);
#endif
#if HAVE_DECL_CACOSHF == 0
extern float _Complex cacoshf(float _Complex z);
#endif
#if HAVE_DECL_CATANHF == 0
extern float _Complex catanhf(float _Complex z);
#endif
#define COPYSIGN copysignf
#define NEXTAFTER  nextafterf
#define MKNAN nanf
#define CEIL ceilf
#define FLOOR floorf
#define NEARBYINT nearbyintf
#define RINT rintf
#define ROUND roundf
#define LRINT lrintf
#define LROUND lroundf
#define LLRINT llrintf
#define LLROUND llroundf
#define TRUNC truncf
#define FMOD fmodf
#define REMAINDER remainderf
#define REMQUO remquof
#define FDIM fdimf
#define FMAX fmaxf
#define FMIN fminf
#define FFMA fmaf
#define FABS fabsf
#define SQRT sqrtf
#define CBRT cbrtf
#define HYPOT hypotf
#define EXP expf
#define EXP2 exp2f
#define EXPM1 expm1f
#define LOG logf
#define LOG2 log2f
#define LOG10 log10f
#define LOG1P log1pf
#define LOGB logbf
#define ILOGB ilogbf
#define MODF modff
#define FREXP frexpf
#define LDEXP ldexpf
#define SCALBN scalbnf
#define SCALBLN scalblnf
#define POW powf
#define COS cosf
#define SIN sinf
#define TAN tanf
#define COSH coshf
#define SINH sinhf
#define TANH tanhf
#define ACOS acosf
#define ASIN asinf
#define ATAN atanf
#define ATAN2 atan2f
#define ACOSH acoshf
#define ASINH asinhf
#define ATANH atanhf
#define TGAMMA tgammaf
#define LGAMMA lgammaf
#define J0 j0f
#define J1 j1f
#define JN jnf
#define Y0 y0f
#define Y1 y1f
#define YN ynf
#define ERF erff
#define ERFC erfcf
#define CREAL crealf
#define CIMAG cimagf
#define CABS cabsf
#define CARG cargf
#define CONJ conjf
#define CPROJ cprojf
#define CSQRT csqrtf
#define CEXP cexpf
#define CLOG clogf
#define CPOW cpowf
#define CSIN csinf
#define CCOS ccosf
#define CTAN ctanf
#define CASIN casinf
#define CACOS cacosf
#define CATAN catanf
#define CSINH csinhf
#define CCOSH ccoshf
#define CTANH ctanhf
#define CASINH casinhf
#define CACOSH cacoshf
#define CATANH catanhf
#else
#if HAVE_DECL_COPYSIGN == 0
extern double copysign(double, double);
#endif
#if HAVE_DECL_NEXTAFTER == 0
extern double nextafter(double, double);
#endif
#if HAVE_DECL_NAN == 0
extern double nan(const char *tag);
#endif
#if HAVE_DECL_CEIL == 0
extern double ceil(double);
#endif
#if HAVE_DECL_FLOOR == 0
extern double floor(double);
#endif
#if HAVE_DECL_NEARBYINT == 0
extern double nearbyint(double);
#endif
#if HAVE_DECL_RINT == 0
extern double rint(double);
#endif
#if HAVE_DECL_ROUND == 0
extern double round(double);
#endif
#if HAVE_DECL_LRINT == 0
extern long int lrint(double);
#endif
#if HAVE_DECL_LROUND == 0
extern long int lround(double);
#endif
#if HAVE_DECL_LLRINT == 0
extern long long int llrint(double);
#endif
#if HAVE_DECL_LLROUND == 0
extern long long int llround(double);
#endif
#if HAVE_DECL_TRUNC == 0
extern double trunc(double);
#endif
#if HAVE_DECL_FMOD == 0
extern double fmod(double, double);
#endif
#if HAVE_DECL_REMAINDER == 0
extern double remainder(double, double);
#endif
#if HAVE_DECL_REMQUO == 0
extern double remquo(double x, double y, int *);
#endif
#if HAVE_DECL_FDIM == 0
extern double fdim(double, double);
#endif
#if HAVE_DECL_FMAX == 0
extern double fmax(double, double);
#endif
#if HAVE_DECL_FMIN == 0
extern double fmin(double, double);
#endif
#if HAVE_DECL_FMA == 0
extern double fma(double x, double y, double z);
#endif
#if HAVE_DECL_FABS == 0
extern double fabs(double);
#endif
#if HAVE_DECL_SQRT == 0
extern double sqrt(double);
#endif
#if HAVE_DECL_CBRT == 0
extern double cbrt(double);
#endif
#if HAVE_DECL_HYPOT == 0
extern double hypot(double, double);
#endif
#if HAVE_DECL_EXP == 0
extern double exp(double);
#endif
#if HAVE_DECL_EXP2 == 0
extern double exp2(double);
#endif
#if HAVE_DECL_EXPM1 == 0
extern double expm1(double);
#endif
#if HAVE_DECL_LOG == 0
extern double log(double);
#endif
#if HAVE_DECL_LOG2 == 0
extern double log2(double);
#endif
#if HAVE_DECL_LOG10 == 0
extern double log10(double);
#endif
#if HAVE_DECL_LOG1P == 0
extern double log1p(double);
#endif
#if HAVE_DECL_LOGB == 0
extern double logb(double);
#endif
#if HAVE_DECL_ILOGB == 0
extern int ilogb(double);
#endif
#if HAVE_DECL_MODF == 0
extern double modf(double, double *);
#endif
#if HAVE_DECL_FREXP == 0
extern double frexp(double, int *);
#endif
#if HAVE_DECL_LDEXP == 0
extern double ldexp(double, int);
#endif
#if HAVE_DECL_SCALBN == 0
extern double scalbn(double, int);
#endif
#if HAVE_DECL_SCALBLN == 0
extern double scalbln(double, long int);
#endif
#if HAVE_DECL_POW == 0
extern double pow(double, double);
#endif
#if HAVE_DECL_COS == 0
extern double cos(double);
#endif
#if HAVE_DECL_SIN == 0
extern double sin(double);
#endif
#if HAVE_DECL_TAN == 0
extern double tan(double);
#endif
#if HAVE_DECL_COSH == 0
extern double cosh(double);
#endif
#if HAVE_DECL_SINH == 0
extern double sinh(double);
#endif
#if HAVE_DECL_TANH == 0
extern double tanh(double);
#endif
#if HAVE_DECL_ACOS == 0
extern double acos(double);
#endif
#if HAVE_DECL_ASIN == 0
extern double asin(double);
#endif
#if HAVE_DECL_ATAN == 0
extern double atan(double);
#endif
#if HAVE_DECL_ATAN2 == 0
extern double atan2(double, double);
#endif
#if HAVE_DECL_ACOSH == 0
extern double acosh(double);
#endif
#if HAVE_DECL_ASINH == 0
extern double asinh(double);
#endif
#if HAVE_DECL_ATANH == 0
extern double atanh(double);
#endif
#if HAVE_DECL_TGAMMA == 0
extern double tgamma(double);
#endif
#if HAVE_DECL_LGAMMA == 0
extern double lgamma(double);
#endif
#if HAVE_DECL_J0 == 0
extern double j0(double);
#endif
#if HAVE_DECL_J1 == 0
extern double j1(double);
#endif
#if HAVE_DECL_JN == 0
extern double jn(int, double);
#endif
#if HAVE_DECL_Y0 == 0
extern double y0(double);
#endif
#if HAVE_DECL_Y1 == 0
extern double y1(double);
#endif
#if HAVE_DECL_YN == 0
extern double yn(int, double);
#endif
#if HAVE_DECL_ERF == 0
extern double erf(double);
#endif
#if HAVE_DECL_ERFC == 0
extern double erfc(double);
#endif
#if HAVE_DECL_CREAL == 0
extern double creal(double _Complex z);
#endif
#if HAVE_DECL_CIMAG == 0
extern double cimag(double _Complex z);
#endif
#if HAVE_DECL_CABS == 0
extern double cabs(double _Complex z);
#endif
#if HAVE_DECL_CARG == 0
extern double carg(double _Complex z);
#endif
#if HAVE_DECL_CONJ == 0
extern double _Complex conj(double _Complex z);
#endif
#if HAVE_DECL_CPROJ == 0
extern double _Complex cproj(double _Complex z);
#endif
#if HAVE_DECL_CSQRT == 0
extern double _Complex csqrt(double _Complex z);
#endif
#if HAVE_DECL_CEXP == 0
extern double _Complex cexp(double _Complex z);
#endif
#if HAVE_DECL_CLOG == 0
extern double _Complex clog(double _Complex z);
#endif
#if HAVE_DECL_CPOW == 0
extern double _Complex cpow(double _Complex z, double _Complex w);
#endif
#if HAVE_DECL_CSIN == 0
extern double _Complex csin(double _Complex z);
#endif
#if HAVE_DECL_CCOS == 0
extern double _Complex ccos(double _Complex z);
#endif
#if HAVE_DECL_CTAN == 0
extern double _Complex ctan(double _Complex z);
#endif
#if HAVE_DECL_CASIN == 0
extern double _Complex casin(double _Complex z);
#endif
#if HAVE_DECL_CACOS == 0
extern double _Complex cacos(double _Complex z);
#endif
#if HAVE_DECL_CATAN == 0
extern double _Complex catan(double _Complex z);
#endif
#if HAVE_DECL_CSINH == 0
extern double _Complex csinh(double _Complex z);
#endif
#if HAVE_DECL_CCOSH == 0
extern double _Complex ccosh(double _Complex z);
#endif
#if HAVE_DECL_CTANH == 0
extern double _Complex ctanh(double _Complex z);
#endif
#if HAVE_DECL_CASINH == 0
extern double _Complex casinh(double _Complex z);
#endif
#if HAVE_DECL_CACOSH == 0
extern double _Complex cacosh(double _Complex z);
#endif
#if HAVE_DECL_CATANH == 0
extern double _Complex catanh(double _Complex z);
#endif
#define COPYSIGN copysign
#define NEXTAFTER  nextafter
#define MKNAN nan
#define CEIL ceil
#define FLOOR floor
#define NEARBYINT nearbyint
#define RINT rint
#define ROUND round
#define LRINT lrint
#define LROUND lround
#define LLRINT llrint
#define LLROUND llround
#define TRUNC trunc
#define FMOD fmod
#define REMAINDER remainder
#define REMQUO remquo
#define FDIM fdim
#define FMAX fmax
#define FMIN fmin
#define FFMA fma
#define FABS fabs
#define SQRT sqrt
#define CBRT cbrt
#define HYPOT hypot
#define EXP exp
#define EXP2 exp2
#define EXPM1 expm1
#define LOG log
#define LOG2 log2
#define LOG10 log10
#define LOG1P log1p
#define LOGB logb
#define ILOGB ilogb
#define MODF modf
#define FREXP frexp
#define LDEXP ldexp
#define SCALBN scalbn
#define SCALBLN scalbln
#define POW pow
#define COS cos
#define SIN sin
#define TAN tan
#define COSH cosh
#define SINH sinh
#define TANH tanh
#define ACOS acos
#define ASIN asin
#define ATAN atan
#define ATAN2 atan2
#define ACOSH acosh
#define ASINH asinh
#define ATANH atanh
#define TGAMMA tgamma
#define LGAMMA lgamma
#define J0 j0
#define J1 j1
#define JN jn
#define Y0 y0
#define Y1 y1
#define YN yn
#define ERF erf
#define ERFC erfc
#define CREAL creal
#define CIMAG cimag
#define CABS cabs
#define CARG carg
#define CONJ conj
#define CPROJ cproj
#define CSQRT csqrt
#define CEXP cexp
#define CLOG clog
#define CPOW cpow
#define CSIN csin
#define CCOS ccos
#define CTAN ctan
#define CASIN casin
#define CACOS cacos
#define CATAN catan
#define CSINH csinh
#define CCOSH ccosh
#define CTANH ctanh
#define CASINH casinh
#define CACOSH cacosh
#define CATANH catanh
#endif

#if defined(FLT_ROUND)
  #if FLT_ROUND != -1
    #define FLTROUND 1.0
  #else
    #define FLTROUND 0.0
  #endif
#else
  #define FLTROUND 0.0
#endif

#if HAVE_DECL_DRAND48 == 0
  extern double drand48(void);
#endif
#if HAVE_DECL_SRAND48 == 0
  extern void srand48(long int);
#endif
#define R_RADIX FLT_RADIX
#define II _Complex_I

/* format strings */
#if defined(NFFT_LDOUBLE)
#  define __FGS__ "Lg"
#  define __FES__ "LE"
#  define __FI__ "%Lf"
#  define __FIS__ "Lf"
#  define __FR__ "%Le"
#elif defined(NFFT_SINGLE)
#  define __FGS__ "g"
#  define __FES__ "E"
#  define __FI__ "%f"
#  define __FIS__ "f"
#  define __FR__ "%e"
#else // double precision
#  define __FGS__ "lg"
#  define __FES__ "lE"
#  define __FI__ "%lf"
#  define __FIS__ "lf"
#  define __FR__ "%le"
#endif

#if MANT_DIG == 113
  // IEEE 754 quadruple precision, 128 bits.
  #define __FE__ "% 36.32LE"
#elif MANT_DIG == 64
  // Intel double extended, 80 bits.
  #define __FE__ "% 24.20LE"
#elif MANT_DIG == 53
  // IEEE 754 double precision, 64 bits.
  #define __FE__ "% 20.16lE"
#elif MANT_DIG == 24
  #define __FE__ "% 12.8E"
#else
  // Unknown floating-point type.
  // Assume IEEE 754 double precision, 64 bits.
  #define __FE__ "% 20.16LE"
#endif

#define TRUE 1
#define FALSE 0

#if defined(_WIN32) || defined(_WIN64)
#  define __D__ "%Id"
#else
#  define __D__ "%td"
#endif

/** Dummy use of unused parameters to silence compiler warnings */
#define UNUSED(x) (void)x

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

/** Return number of elapsed seconds between two time points. */
R Y(elapsed_seconds)(ticks t1, ticks t0);

/** Dummy use of unused parameters to silence compiler warnings */
#define UNUSED(x) (void)x

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

/* sinc.c: */

/* Sinus cardinalis. */
R Y(sinc)(R x);

/* log|sinc(x)|; even, 0 at x = 0, -inf at every other multiple of pi */
R Y(log_sinc)(R x);

/* lambda.c: */

/* lambda(z, eps) = gamma(z + eps) / gamma(z + 1) */
R Y(lambda)(R z, R eps);

/* lambda2(mu, nu) = sqrt(gamma(mu + nu + 1) / (gamma(mu + 1) * gamma(nu + 1))) */
R Y(lambda2)(R mu, R nu);

/* bessel_i0.c: I0 in the two scaled forms the Kaiser-Bessel window needs.
 * Neither exponentiates a large argument, so neither inherits the relative
 * error of EXP, which grows with the argument on some targets. I0 itself is not
 * provided: nothing needs it, and every form of it would carry that error. */
R Y(bessel_i0_logtail)(R x); /* log I0(x) - x, <= 0 */
R Y(bessel_i0_exp_scaled)(R x); /* exp(-x) * I0(x), in (0, 1] */

/* Kaiser-Bessel window, scaled by 1/I0(m b) so that phi_hut(0) == 1.
 *
 *   phi(nx) = sinh(b ra) / (pi ra I0(m b)),  ra = sqrt(m^2 - nx^2),  |nx| < m
 *             sin (b ri) / (pi ri I0(m b)),  ri = sqrt(nx^2 - m^2),  |nx| > m
 *             b / (pi I0(m b)),                                      |nx| = m
 *
 * The argument is nx = n x rather than x: every caller has n x to hand, and the
 * points of a run differ by exactly one in nx, so the run form needs no divide.
 *
 * Constants per axis: lg_tail = log I0(m b) - m b, peak_inv = 1 / I0(m b).
 */

/* |nx| >= m. */
static inline R Y(kb_phi_out)(R b, R peak_inv, R m, R nx)
{
  const R a = (nx - m) * (nx + m);

  if (a <= K(0.0))
    return (b / KPI) * peak_inv;

  {
    const R ri = SQRT(a);
    return (SIN(b * ri) / (KPI * ri)) * peak_inv;
  }
}

/* |nx| <= m, accurate down to ra = 0 and free of overflow at any m b.
 *
 * sinh(b ra)/I0(m b) = exp(b ra - m b - lg_tail) * (1 - exp(-2 b ra)) / 2, with
 * b ra - m b = -b nx^2/(ra + m) so that the peak does not cancel. One
 * reciprocal serves both 1/(ra + m) and 1/ra. */
static inline R Y(kb_phi_in)(R b, R lg_tail, R peak_inv, R m, R nx)
{
  const R a = (m - nx) * (m + nx);

  if (a <= K(0.0))
    return (b / KPI) * peak_inv;

  {
    const R ra = SQRT(a);
    const R w = K(1.0) / (ra * (ra + m));

    return EXP(-b * (nx * nx) * (w * ra) - lg_tail)
      * (-EXPM1(K(-2.0) * b * ra)) * (w * (ra + m)) * (K(0.5) / KPI);
  }
}

/* |nx| <= m - 1, so b ra >= b sqrt(2m - 1) and exp(-2 b ra) is far below 1.
 *
 * With E = exp(b ra - m b - lg_tail) the second exponential is redundant:
 * exp(-2 b ra) = peak_inv^2 / E^2, so 1 - exp(-2 b ra) costs one reciprocal
 * instead of an EXPM1, and E - peak_inv^2/E cannot cancel while b ra is bounded
 * away from zero. peak_inv^2 underflows only past m b ~ 357, where exp(-2 b ra)
 * is already far below the format epsilon and zero is the right value. */
static inline R Y(kb_phi_in_interior)(R b, R lg_tail, R peak_inv_sq, R m, R nx)
{
  const R ra = SQRT((m - nx) * (m + nx));
  const R w = K(1.0) / (ra * (ra + m));
  const R e = EXP(-b * (nx * nx) * (w * ra) - lg_tail);

  return (e - peak_inv_sq / e) * (w * (ra + m)) * (K(0.5) / KPI);
}

static inline R Y(kb_phi)(R b, R lg_tail, R peak_inv, R m, R nx)
{
  if ((m - nx) * (m + nx) > K(0.0))
    return Y(kb_phi_in)(b, lg_tail, peak_inv, m, nx);

  return Y(kb_phi_out)(b, peak_inv, m, nx);
}

/* dst[l] = phi(nx0 - l) for l = 0 .. 2m+1.
 *
 * A run straddles the support: its ends lie outside, and next to them ra comes
 * arbitrarily close to zero. The stretch with |nx| <= m - 1 keeps b ra at or
 * above b sqrt(2m - 1), which is what the interior kernel needs; it is found by
 * inverting |nx0 - l| <= m - 1 rather than assumed, so a caller whose run sits
 * differently still gets the guarded evaluation everywhere it is due. The
 * interior loop has no branch and one transcendental, so it vectorises. */
static inline void Y(kb_phi_run)(R *dst, R b, R lg_tail, R peak_inv, R m,
  INT mi, R nx0)
{
  const INT last = 2 * mi + 1;
  const R lim = m - K(1.0);
  INT lo = (INT)CEIL(nx0 - lim);
  INT hi = (INT)FLOOR(nx0 + lim);
  INT l;

  if (lo < 0)
    lo = 0;
  if (hi > last)
    hi = last;

  for (l = 0; l < lo; l++)
    dst[l] = Y(kb_phi)(b, lg_tail, peak_inv, m, nx0 - (R)l);

  if (lo <= hi)
  {
    const R peak_inv_sq = peak_inv * peak_inv;

    for (l = lo; l <= hi; l++)
      dst[l] = Y(kb_phi_in_interior)(b, lg_tail, peak_inv_sq, m, nx0 - (R)l);
  }
  else
    hi = lo - 1;

  for (l = hi + 1; l <= last; l++)
    dst[l] = Y(kb_phi)(b, lg_tail, peak_inv, m, nx0 - (R)l);
}

/* I0(a)/I0(m b) with a = m sqrt(b^2 - t^2), t = 2 pi k / n. Both exponentially
 * scaled Bessel values lie in (0, 1] and a - m b = -m t^2/(ra + b) is formed
 * without the cancellation, so nothing overflows at any m b.
 * i0e_peak_inv = 1/(exp(-m b) I0(m b)) is constant per axis. */
static inline R Y(kb_phi_hut)(R b, R i0e_peak_inv, R m, R n, R k)
{
  const R t = K(2.0) * KPI * k / n;
  const R s = (b - t) * (b + t); /* not b*b - t*t; >= 0 in band, so SQRT is safe */
  const R ra = SQRT(s);

  return EXP(-m * t * t / (ra + b)) * Y(bessel_i0_exp_scaled)(m * ra)
    * i0e_peak_inv;
}

/* bspline.c: */
R Y(bsplines)(const INT, const R x);

/* Chebyshev series of every piece of B_k; tab holds k rows of k reals. */
void Y(bspline_cheb_init)(R *tab, const INT k);

/* First piece whose Clenshaw evaluation keeps relative accuracy within thresh;
 * the admissible pieces are i0 through k-1-i0. */
INT Y(bspline_cheb_guard)(const R *tab, const INT k, const R thresh);

/* B_k(x) from the Chebyshev table, in O(k) rather than de Boor's O(k^2).
 *
 * Clenshaw's error is bounded by the coefficient sum, so this holds accuracy
 * against the peak of B_k everywhere but loses relative accuracy in the tail,
 * where B_k runs below the peak by any number of orders. Callers that divide by
 * the result want Y(bspline_cheb_rel) instead. */
static inline R Y(bspline_cheb)(const R *tab, const INT k, const R x)
{
  const INT i = (INT)FLOOR(x);
  const R *c;
  R s, s2, b1 = K(0.0), b2 = K(0.0);
  INT l;

  if (x <= K(0.0) || x >= (R)k)
    return K(0.0);

  c = tab + i * k;
  s = K(2.0) * (x - (R)i) - K(1.0);
  s2 = K(2.0) * s;

  for (l = k - 1; l >= 1; l--)
  {
    const R b = s2 * b1 - b2 + c[l];

    b2 = b1;
    b1 = b;
  }

  return s * b1 - b2 + c[0];
}

/* B_k(x) with relative accuracy, for callers that divide by it. Outside the
 * pieces Y(bspline_cheb_guard) admits it falls back to de Boor, which costs
 * O(k^2) but is reached only where the deconvolution is already the transform's
 * limit. */
static inline R Y(bspline_cheb_rel)(const R *tab, const INT k, const INT i0,
  const R x)
{
  const INT i = (INT)FLOOR(x);

  if (x <= K(0.0) || x >= (R)k)
    return K(0.0);

  if (i < i0 || i >= k - i0)
    return Y(bsplines)(k, x);

  return Y(bspline_cheb)(tab, k, x);
}

/* sinc(t)^(2m)/n with t = k pi / n. Formed as exp(2m log|sinc t|), because 2m
 * multiplies the rounding error of sinc but only the absolute error of its
 * logarithm. Here |t| <= pi/(2 sigma), so log_sinc never leaves its polynomial
 * branch, and k = 0 needs no special case: log_sinc(0) is zero. */
static inline R Y(bspline_phi_hut)(R m, R n, R k)
{
  return EXP(K(2.0) * m * Y(log_sinc)(k * KPI / n)) / n;
}

/* dst[l] = B_{2m}(t - l) / n for l = 0 .. 2m+1, with inv[j] = 1/(j+1).
 *
 * The run arguments are one grid cell apart, so the nonzero entries are the 2m
 * cardinal B-splines of order 2m at a single point, and de Boor's scheme gives
 * all of them in one triangular sweep of O(m^2) rather than one sweep each.
 * The recurrence is linear in its seed, so the 1/n enters there instead of
 * costing a rounding per entry.
 *
 * The nonzero stretch is derived from t, so a run that sits differently against
 * the support still gets one evaluation per point. */
static inline void Y(bspline_phi_run)(R *dst, const R *inv, INT mi, R t,
  R inv_n)
{
  const INT k = 2 * mi;
  const INT last = 2 * mi + 1;
  const INT i = (INT)FLOOR(t);
  const INT lo = i - k + 1;
  INT j, l, r;

  if (lo < 0 || lo + k > last + 1)
  {
    for (l = 0; l <= last; l++)
      dst[l] = Y(bsplines)(k, t - (R)l) * inv_n;

    return;
  }

  {
    const R g = t - (R)i;
    R *b = dst + lo;

    for (l = 0; l < lo; l++)
      dst[l] = K(0.0);

    for (l = lo + k; l <= last; l++)
      dst[l] = K(0.0);

    b[0] = inv_n;

    for (j = 0; j < k - 1; j++)
    {
      const R iv = inv[j];
      R saved = K(0.0);

      for (r = 0; r <= j; r++)
      {
        const R term = b[r] * iv;

        b[r] = saved + ((R)(r + 1) - g) * term;
        saved = (g + (R)(j - r)) * term;
      }

      b[j + 1] = saved;
    }
  }
}

/* w sinc(a)^(2m). The exponent is even, so the sign of sinc past its first zero
 * does not matter and the logarithm may take the absolute value. */
static inline R Y(sincpow_phi)(R w, R m, R a)
{
  return w * EXP(K(2.0) * m * Y(log_sinc)(a));
}

/* dst[l] = phi(x - (u + l)/n) for l = 0 .. 2m+1, with the sinc argument
 * a_l = h (nx0 - l). Both h and the prefactor are constant per axis. */
static inline void Y(sincpow_phi_run)(R *dst, R w, R m, INT mi, R h, R nx0)
{
  const INT last = 2 * mi + 1;
  const R e = K(2.0) * m;
  INT l;

  for (l = 0; l <= last; l++)
    dst[l] = w * EXP(e * Y(log_sinc)(h * (nx0 - (R)l)));
}

/* Gaussian window, phi(nx) = exp(-nx^2/b)/sqrt(pi b) and phi_hut(t) =
 * exp(-t^2 b) with t = k pi / n.
 *
 * b_inv = 1/b and norm = 1/sqrt(pi b) are constant per axis, so neither
 * evaluation divides or takes a square root, and the argument is nx = n x
 * rather than x: every caller has n x to hand, and the points of a run differ
 * by exactly one in nx. */
static inline R Y(gaussian_phi_hut)(R b, R t)
{
  return EXP(-(t * t) * b);
}

static inline R Y(gaussian_phi)(R b_inv, R norm, R nx)
{
  return EXP(-(nx * nx) * b_inv) * norm;
}

/* dst[l] = phi(nx0 - l) for l = 0 .. 2m+1.
 *
 * One exponential per point, the same as the generic fill, but the argument is
 * stepped by one in nx rather than divided by n per point, and the two
 * constants are already formed. The loop carries nothing, so it vectorises,
 * which is what makes it the right shape here: a Gaussian run can also be had
 * from two exponentials and a geometric recurrence outward from its peak, the
 * way FG_RUN does it, and that measures faster only past m of about 12 -- the
 * serial dependency costs more than the exponentials it saves while the run is
 * short. Measured at sigma = 2 on one 1-D transform of 20000 nodes, the
 * recurrence is 12% slower at m = 4 and 13% faster at m = 24 against this
 * loop, and the crossover moves with the vector width and with whether the
 * libm has a vector exp at all; the cutoffs nfft_init picks straddle it, so
 * the form that does not depend on where it falls is the one to have.
 *
 * A run needs no clamping or special case: phi is entire and positive, so
 * every point of every run is an ordinary evaluation. */
static inline void Y(gaussian_phi_run)(R *dst, R b_inv, R norm, INT mi, R nx0)
{
  const INT last = 2 * mi + 1;
  INT l;

  for (l = 0; l <= last; l++)
  {
    const R t = nx0 - (R)l;

    dst[l] = EXP(-(t * t) * b_inv) * norm;
  }
}

/* float.c: */
typedef enum {NFFT_EPSILON = 0, NFFT_SAFE__MIN = 1, NFFT_BASE = 2,
  NFFT_PRECISION = 3, NFFT_MANT_DIG = 4, NFFT_FLTROUND = 5, NFFT_E_MIN = 6,
  NFFT_R_MIN = 7, NFFT_E_MAX = 8, NFFT_R_MAX = 9 } float_property;

R Y(float_property)(float_property);
R Y(prod_real)(R *vec, INT d);

/* int.c: */
INT Y(log2i)(const INT m);
void Y(next_power_of_2_exp)(const INT N, INT *N2, INT *t);
void Y(next_power_of_2_exp_int)(const int N, int *N2, int *t);

/* error.c: */
/* not used */ R Y(error_l_infty_double)(const R *x, const R *y, const INT n);
/* not used */ R Y(error_l_infty_1_double)(const R *x, const R *y, const INT n, const R *z,
  const INT m);
R Y(error_l_2_complex)(const C *x, const C *y, const INT n);
/* not used */ R Y(error_l_2_double)(const R *x, const R *y, const INT n);

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

/* voronoi.c */
void Y(voronoi_weights_1d)(R *w, R *x, const INT M);

/* damp.c */
/**
 * Compute damping factor for modified Fejer kernel:
 * /f$\frac{2}{N}\left(1-\frac{\left|2k+1\right|}{N}\right)/f$
 */
R Y(modified_fejer)(const INT N, const INT kk);
/** Compute damping factor for modified Jackson kernel. */
R Y(modified_jackson2)(const INT N, const INT kk);
/** Compute damping factor for modified generalised Jackson kernel. */
R Y(modified_jackson4)(const INT N, const INT kk);
/** Compute damping factor for modified Sobolev kernel. */
R Y(modified_sobolev)(const R mu, const INT kk);
/** Comput damping factor for modified multiquadric kernel. */
R Y(modified_multiquadric)(const R mu, const R c, const INT kk);

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

#endif
