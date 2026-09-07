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

/* Helper functions (random data, printing, timing, error norms, ...) useful
 * to callers and to the example/application programs, but not themselves
 * part of any transform's API. See nfft3.h for the transform APIs. */

#ifndef __NFFT3UTIL_H__
#define __NFFT3UTIL_H__

#include "nfft3.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/* float.c: machine-precision parameters, LAPACK DLAMCH-style. */
typedef enum {NFFT_EPSILON = 0, NFFT_SAFE__MIN = 1, NFFT_BASE = 2,
  NFFT_PRECISION = 3, NFFT_MANT_DIG = 4, NFFT_FLTROUND = 5, NFFT_E_MIN = 6,
  NFFT_R_MIN = 7, NFFT_E_MAX = 8, NFFT_R_MAX = 9} float_property;

/* huge second-order macro that defines prototypes for all utility API functions.
 * We expand this macro for each supported precision.
 *   Y: nfft name-mangling macro
 *   R: real data type
 *   C: complex data type
 */
#define NFFT_DEFINE_UTIL_API(Y,R,C) \
/* rand.c */ \
NFFT_EXTERN R Y(drand48)(void); \
NFFT_EXTERN void Y(srand48)(long int seed); \
\
/** Inits a vector of random complex numbers in \f$[0,1]\times[0,1]{\rm i}\f$. \
 */ \
NFFT_EXTERN void Y(vrand_unit_complex)(C *x, const NFFT_INT n); \
\
/** Inits a vector of random double numbers in \f$[-1/2,1/2]\f$. \
 */ \
NFFT_EXTERN void Y(vrand_shifted_unit_double)(R *x, const NFFT_INT n); \
\
NFFT_EXTERN void Y(vrand_real)(R *x, const NFFT_INT n, const R a, const R b); \
\
/* print.c */ \
/** Print real vector to standard output. */ \
NFFT_EXTERN void Y(vpr_double)(R *x, const NFFT_INT n, const char *text); \
\
/** Print complex vector to standard output. */ \
NFFT_EXTERN void Y(vpr_complex)(C *x, const NFFT_INT n, const char *text); \
/* thread.c */ \
NFFT_EXTERN NFFT_INT Y(get_num_threads)(void); \
NFFT_EXTERN void Y(set_num_threads)(NFFT_INT nthreads); \
NFFT_EXTERN NFFT_INT Y(has_threads_enabled)(void); \
/* time.c */ \
/** Monotonic wall clock in seconds. The origin is arbitrary; take \
 * differences to measure an interval. Always double: an absolute timestamp \
 * is not a transform datum, and float quantises it to tens of ms. */ \
NFFT_EXTERN double Y(clock_gettime_seconds)(void); \
/* float.c */ \
NFFT_EXTERN R Y(float_property)(float_property p); \
/* error.c: */ \
NFFT_EXTERN R Y(error_l_infty_complex)(const C *x, const C *y, const NFFT_INT n); \
NFFT_EXTERN R Y(error_l_infty_1_complex)(const C *x, const C *y, const NFFT_INT n, \
  const C *z, const NFFT_INT m); \
NFFT_EXTERN R Y(error_l_2_complex)(const C *x, const C *y, const NFFT_INT n); \
/* int.c: */ \
NFFT_EXTERN NFFT_INT Y(exp2i)(const NFFT_INT a); \
NFFT_EXTERN NFFT_INT Y(next_power_of_2)(const NFFT_INT N); \
/* vector1.c */ \
/** Computes the inner/dot product \f$x^H x\f$. */ \
NFFT_EXTERN R Y(dot_complex)(C *x, NFFT_INT n); \
/* vector3.c */ \
/** Updates \f$x \leftarrow a x + y\f$. */ \
NFFT_EXTERN void Y(upd_axpy_complex)(C *x, R a, C *y, NFFT_INT n); \
/** Swaps each half over N[d]/2. */ \
NFFT_EXTERN void Y(fftshift_complex)(C *x, NFFT_INT d, NFFT_INT* N); \
NFFT_EXTERN void Y(fftshift_complex_int)(C *x, int d, int* N); \
/** Return library version. */ \
NFFT_EXTERN void Y(get_version)(unsigned *major, unsigned *minor, unsigned *patch); \
/** \
 * Return name of window function. \
 * \
 * The window function to be used is configured at compile time. \
 */ \
NFFT_EXTERN const char *Y(get_window_name)(void); \
NFFT_EXTERN NFFT_INT Y(get_default_window_cut_off)(void);

NFFT_DEFINE_UTIL_API(NFFT_MANGLE_FLOAT,float,fftwf_complex)
NFFT_DEFINE_UTIL_API(NFFT_MANGLE_DOUBLE,double,fftw_complex)
NFFT_DEFINE_UTIL_API(NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* defined(__NFFT3UTIL_H__) */
