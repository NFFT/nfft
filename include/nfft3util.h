/*
 * Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
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

/*! \file nfft3util.h
 *  \brief Header file for utility functions used by the nfft3 library.
 */
#ifndef __UTIL_H__
#define __UTIL_H__

/** Include header for FFTW3 library for its complex type. */
#include <fftw3.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/**
 * @defgroup nfftutil Util - Auxilliary functions
 * @{
 *
 * This module implements frequently used utility functions.
 * In particular, this includes simple measurement of resources, evaluation of
 * window functions, vector routines for basic linear algebra tasks, and
 * computation of weights for the inverse transforms.
 *
 */

/** Swapping of two vectors.
 */
#define NFFT_SWAP_complex(x,y) {fftw_complex* NFFT_SWAP_temp__; \
  NFFT_SWAP_temp__=(x); (x)=(y); (y)=NFFT_SWAP_temp__;}

/** Swapping of two vectors.
 */
#define NFFT_SWAP_double(x,y) {double* NFFT_SWAP_temp__; NFFT_SWAP_temp__=(x); \
  (x)=(y); (y)=NFFT_SWAP_temp__;}

/** Maximum of its two arguments.
 */
#define NFFT_MAX(a,b) ((a)>(b)? (a) : (b))

/** Mimimum of its two arguments.
 */
#define NFFT_MIN(a,b) ((a)<(b)? (a) : (b))

/* ######################################################################### */
/* ########## Helpers for inverse transforms ############################### */
/* ######################################################################### */

int nfft_smbi(const double x, const double alpha, const int nb, const int ize,
  double *b);

int nfft_get_num_threads(void);

#ifdef _OPENMP
int nfft_get_omp_num_threads(void);
#endif

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

/** @}
 */
#endif
