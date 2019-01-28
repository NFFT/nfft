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

#include "nfft3.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#if defined(NFFT_PRECISION_SINGLE)
typedef float NFFT_R;
typedef fftwf_complex NFFT_C;
#define NFFT_K(x) ((NFFT_R) x)
#define NFFT_M(name) NFFT_CONCAT(name,f)
#define FFTW(name) NFFT_CONCAT(fftwf_,name)
#define NFFT(name) NFFT_CONCAT(nfftf_,name)
#define NFCT(name) NFFT_CONCAT(nfctf_,name)
#define NFST(name) NFFT_CONCAT(nfstf_,name)
#define NFSFT(name) NFFT_CONCAT(nfsftf_,name)
#define SOLVER(name) NFFT_CONCAT(solverf_,name)
#elif defined(NFFT_PRECISION_LONG_DOUBLE)
typedef long double NFFT_R;
typedef fftwl_complex NFFT_C;
#define NFFT_K(x) ((NFFT_R) x##L)
#define NFFT_M(name) NFFT_CONCAT(name,l)
#define FFTW(name) NFFT_CONCAT(fftwl_,name)
#define NFFT(name) NFFT_CONCAT(nfftl_,name)
#define NFCT(name) NFFT_CONCAT(nfctl_,name)
#define NFST(name) NFFT_CONCAT(nfstl_,name)
#define NFSFT(name) NFFT_CONCAT(nfsftl_,name)
#define SOLVER(name) NFFT_CONCAT(solverl_,name)
#elif defined(NFFT_PRECISION_DOUBLE)
typedef double NFFT_R;
typedef fftw_complex NFFT_C;
#define NFFT_K(x) ((NFFT_R) x)
#define NFFT_M(name) name
#define FFTW(name) NFFT_CONCAT(fftw_,name)
#define NFFT(name) NFFT_CONCAT(nfft_,name)
#define NFCT(name) NFFT_CONCAT(nfct_,name)
#define NFST(name) NFFT_CONCAT(nfst_,name)
#define NFSFT(name) NFFT_CONCAT(nfsft_,name)
#define SOLVER(name) NFFT_CONCAT(solver_,name)
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

#if defined(_WIN32) || defined(_WIN64)
#  define NFFT__D__ "%Id"
#else
#  define NFFT__D__ "%td"
#endif

#endif /* defined(__NFFT3MP_H__) */
