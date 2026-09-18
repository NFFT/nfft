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

/* Helpers that belong to no transform: random data, printing, timing, error
 * norms. */

#ifndef __NFFT3UTIL_H__
#define __NFFT3UTIL_H__

#include "nfft3.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/* Utility API prototypes, expanded once per precision. */
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
/* time.c */ \
/** Wall clock in seconds, monotonic where available, arbitrary origin. */ \
NFFT_EXTERN double Y(clock_gettime_seconds)(void); \
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
NFFT_EXTERN void Y(fftshift_complex_int)(C *x, int d, int* N);

NFFT_DEFINE_UTIL_API(NFFT_MANGLE_FLOAT,float,fftwf_complex)
NFFT_DEFINE_UTIL_API(NFFT_MANGLE_DOUBLE,double,fftw_complex)
NFFT_DEFINE_UTIL_API(NFFT_MANGLE_LONG_DOUBLE,long double,fftwl_complex)

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* defined(__NFFT3UTIL_H__) */
