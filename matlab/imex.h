/*
 * Copyright (c) 2002, 2016 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* NFFT mex internal header file */
#ifndef MEXUTIL_H
#define MEXUTIL_H

#include <math.h>

#ifdef HAVE_MATLAB_GCC_REQUIRE_UNDEF_STDC_UTF_16
  #undef __STDC_UTF_16__
#endif

#include <mex.h>
#include <matrix.h>

#ifdef HAVE_MATLAB_GCC_REQUIRE_UNDEF_STDC_UTF_16
  #define __STDC_UTF_16__
#endif

/*----------------------------------------------------------------------------*/
/* Replacements for nfft_malloc and nfft_free plus install routine */
extern void *nfft_mex_malloc(size_t n);
extern void nfft_mex_free(void *p);
extern void nfft_mex_install_mem_hooks(void);

int nfft_mex_get_int(const mxArray *p, const char *errmsg);
double nfft_mex_get_double(const mxArray *p, const char *errmsg);

#ifdef MATLAB_ARGCHECKS
#define DM(Y) Y
#else
#define DM(Y)
#endif

/*----------------------------------------------------------------------------*/
/* Checks if argument is a scalar. */
#define ARG_CHECK_SCALAR(x,y) \
if (mxGetM(prhs[x]) != 1 || mxGetN(prhs[x]) != 1 || mxIsDouble(prhs[x]) != 1) \
  mexErrMsgTxt(#y " argument must be a scalar.");

/* Gets and stores pointer to argument data. */
#define ARG_GET_PTR(x,y) \
y = mxGetPr(prhs[x]);

/* Checks and get argument as nonnegative integer. */
#define ARG_GET_NONNEG_INT(x,y,z) \
ARG_CHECK_SCALAR(x,y) \
ARG_GET_PTR(x,z) \
if (z[0] != round(z[0]) || z[0] < 0) \
  mexErrMsgTxt(#y " argument must be a nonnegative integer.");

/* Checks and get argument as positive integer. */
#define ARG_GET_POS_INT(x,y,z) \
ARG_CHECK_SCALAR(x,y) \
ARG_GET_PTR(x,z) \
if (z[0] != round(z[0]) || z[0] < 1) \
  mexErrMsgTxt(#y " argument must be a nonnegative integer.");

/* Checks and get argument as nonnegative double. */
#define ARG_GET_NONNEG_DOUBLE(x,y,z) \
ARG_CHECK_SCALAR(x,y) \
ARG_GET_PTR(x,z) \
if (z[0] < 0) \
  mexErrMsgTxt(#y " argument must be a nonnegative number.");

#endif
