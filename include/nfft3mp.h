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

#ifndef __NFFT3MP_H__
#define __NFFT3MP_H__

/* Canonical include order. <complex.h> must precede <fftw3.h> if the C99 complex 
 * type should be used. The header nfft3.h pulls it in. Otherwise, fftw_complex 
 * (and hence NFFT_C) is a 2-element real array instead. In a program:
 *
 *   #include <complex.h>
 *   #include <nfft3.h>
 *   #define NFFT_PRECISION_SINGLE  // or NFFT_PRECISION_DOUBLE/NFFT_PRECISION_LONG_DOUBLE, depending on desired type
 *   #include <nfft3mp.h>
 *
 * NFFT_R/NFFT_C are then the real/complex types for the selected floating-point type
 * so code can be written type-agnostic.
 */
#ifndef FFTW3_H
#error include nfft3.h (or <fftw3.h>) before nfft3mp.h
#endif

#include <float.h>
#include <math.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#define NFFT_MP_CONCAT(prefix, name) prefix ## name

#if defined(NFFT_PRECISION_SINGLE)
typedef float NFFT_R;
#define NFFT_C fftwf_complex
#define NFFT_K(x) ((NFFT_R) x)
#define NFFT_R_MANT_DIG FLT_MANT_DIG
#define NFFT_R_EPSILON FLT_EPSILON
#define NFFT_M(name) NFFT_MP_CONCAT(name,f)
#define FFTW(name) NFFT_MP_CONCAT(fftwf_,name)
#define NFFT(name) NFFT_MP_CONCAT(nfftf_,name)
#define NFCT(name) NFFT_MP_CONCAT(nfctf_,name)
#define NFST(name) NFFT_MP_CONCAT(nfstf_,name)
#define NFSFT(name) NFFT_MP_CONCAT(nfsftf_,name)
#define SOLVER(name) NFFT_MP_CONCAT(solverf_,name)
#elif defined(NFFT_PRECISION_LONG_DOUBLE)
typedef long double NFFT_R;
#define NFFT_C fftwl_complex
#define NFFT_K(x) ((NFFT_R) x##L)
#define NFFT_R_MANT_DIG LDBL_MANT_DIG
#define NFFT_R_EPSILON LDBL_EPSILON
#define NFFT_M(name) NFFT_MP_CONCAT(name,l)
#define FFTW(name) NFFT_MP_CONCAT(fftwl_,name)
#define NFFT(name) NFFT_MP_CONCAT(nfftl_,name)
#define NFCT(name) NFFT_MP_CONCAT(nfctl_,name)
#define NFST(name) NFFT_MP_CONCAT(nfstl_,name)
#define NFSFT(name) NFFT_MP_CONCAT(nfsftl_,name)
#define SOLVER(name) NFFT_MP_CONCAT(solverl_,name)
#elif defined(NFFT_PRECISION_DOUBLE)
typedef double NFFT_R;
#define NFFT_C fftw_complex
#define NFFT_K(x) ((NFFT_R) x)
#define NFFT_R_MANT_DIG DBL_MANT_DIG
#define NFFT_R_EPSILON DBL_EPSILON
#define NFFT_M(name) name
#define FFTW(name) NFFT_MP_CONCAT(fftw_,name)
#define NFFT(name) NFFT_MP_CONCAT(nfft_,name)
#define NFCT(name) NFFT_MP_CONCAT(nfct_,name)
#define NFST(name) NFFT_MP_CONCAT(nfst_,name)
#define NFSFT(name) NFFT_MP_CONCAT(nfsft_,name)
#define SOLVER(name) NFFT_MP_CONCAT(solver_,name)
#else
#error Either define macro NFFT_PRECISION_SINGLE, NFFT_PRECISION_DOUBLE or NFFT_PRECISION_LONG_DOUBLE for single, double or long double precision
#endif

/* math functions (C99/C11 <math.h> and <complex.h>), precision-dispatched */
#if defined(NFFT_PRECISION_SINGLE)
#define NFFT_MKNAN nanf
#define NFFT_CEIL ceilf
#define NFFT_FLOOR floorf
#define NFFT_ROUND roundf
#define NFFT_LRINT lrintf
#define NFFT_FMAX fmaxf
#define NFFT_FABS fabsf
#define NFFT_SQRT sqrtf
#define NFFT_EXP expf
#define NFFT_LOG logf
#define NFFT_POW powf
#define NFFT_COS cosf
#define NFFT_SIN sinf
#define NFFT_TAN tanf
#define NFFT_TGAMMA tgammaf
#define NFFT_CREAL crealf
#define NFFT_CIMAG cimagf
#define NFFT_CABS cabsf
#define NFFT_CEXP cexpf
#define NFFT_CPOW cpowf
#elif defined(NFFT_PRECISION_LONG_DOUBLE)
#define NFFT_MKNAN nanl
#define NFFT_CEIL ceill
#define NFFT_FLOOR floorl
#define NFFT_ROUND roundl
#define NFFT_LRINT lrintl
#define NFFT_FMAX fmaxl
#define NFFT_FABS fabsl
#define NFFT_SQRT sqrtl
#define NFFT_EXP expl
#define NFFT_LOG logl
#define NFFT_POW powl
#define NFFT_COS cosl
#define NFFT_SIN sinl
#define NFFT_TAN tanl
#define NFFT_TGAMMA tgammal
#define NFFT_CREAL creall
#define NFFT_CIMAG cimagl
#define NFFT_CABS cabsl
#define NFFT_CEXP cexpl
#define NFFT_CPOW cpowl
#elif defined(NFFT_PRECISION_DOUBLE)
#define NFFT_MKNAN nan
#define NFFT_CEIL ceil
#define NFFT_FLOOR floor
#define NFFT_ROUND round
#define NFFT_LRINT lrint
#define NFFT_FMAX fmax
#define NFFT_FABS fabs
#define NFFT_SQRT sqrt
#define NFFT_EXP exp
#define NFFT_LOG log
#define NFFT_POW pow
#define NFFT_COS cos
#define NFFT_SIN sin
#define NFFT_TAN tan
#define NFFT_TGAMMA tgamma
#define NFFT_CREAL creal
#define NFFT_CIMAG cimag
#define NFFT_CABS cabs
#define NFFT_CEXP cexp
#define NFFT_CPOW cpow
#else
#error Either define macro NFFT_PRECISION_SINGLE, NFFT_PRECISION_DOUBLE or NFFT_PRECISION_LONG_DOUBLE for single, double or long double precision
#endif

/* format strings */
#if defined(NFFT_PRECISION_LONG_DOUBLE)
#  define NFFT__FGS__ "Lg"
#  define NFFT__FES__ "LE"
#  define NFFT__FE__ "% 36.32LE"
#  define NFFT__FI__ "%Lf"
#  define NFFT__FIS__ "Lf"
#  define NFFT__FR__ "%Le"
#elif defined(NFFT_PRECISION_SINGLE)
#  define NFFT__FGS__ "g"
#  define NFFT__FES__ "E"
#  define NFFT__FE__ "% 12.8E"
#  define NFFT__FI__ "%f"
#  define NFFT__FIS__ "f"
#  define NFFT__FR__ "%e"
#elif defined(NFFT_PRECISION_DOUBLE)
#  define NFFT__FGS__ "lg"
#  define NFFT__FES__ "lE"
#  define NFFT__FE__ "% 20.16lE"
#  define NFFT__FI__ "%lf"
#  define NFFT__FIS__ "lf"
#  define NFFT__FR__ "%le"
#else
#error Either define macro NFFT_PRECISION_SINGLE, NFFT_PRECISION_DOUBLE or NFFT_PRECISION_LONG_DOUBLE for single, double or long double precision
#endif

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

/** Swap two vectors. */
#define NFFT_CSWAP(x,y) {NFFT_C* NFFT_SWAP_temp__; \
  NFFT_SWAP_temp__=(x); (x)=(y); (y)=NFFT_SWAP_temp__;}

#define NFFT_KPI NFFT_K(3.1415926535897932384626433832795028841971693993751)
#define NFFT_K2PI NFFT_K(6.2831853071795864769252867665590057683943387987502)

#define NFFT_MIN(a,b) (((a)<(b))?(a):(b))
#define NFFT_MAX(a,b) (((a)>(b))?(a):(b))
#define NFFT_ABS(x) (((x)>NFFT_K(0.0))?(x):(-(x)))

#define NFFT_II _Complex_I

/** Radix of the floating-point representation. */
#define NFFT_R_RADIX FLT_RADIX

#define NFFT_UNUSED(x) (void)x

#if defined(_WIN32) || defined(_WIN64)
#  define NFFT__D__ "%Id"
#else
#  define NFFT__D__ "%td"
#endif

#endif /* defined(__NFFT3MP_H__) */
