/*
 * Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis
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

/* NFFT mex internal header file */
#ifndef MEXUTIL_H
#define MEXUTIL_H

#include <math.h>
#include <mex.h>
#include <matrix.h>

/*----------------------------------------------------------------------------*/
/* Replacements for nfft_malloc and nfft_free plus install routine */
extern void *nfft_mex_malloc(size_t n);
extern void nfft_mex_free(void *p);
extern void install_mem_hooks(void);

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
