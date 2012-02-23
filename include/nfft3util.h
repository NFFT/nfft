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

/** Formerly known to be an irrational number.
 */
#define PI 3.141592653589793238462643383279502884197169399375105820974944592
#define PI2 6.283185307179586476925286766559005768394338798750211641949889185
#define PI4 12.56637061435917295385057353311801153678867759750042328389977837

/** Maximum of its two arguments.
 */
#define NFFT_MAX(a,b) ((a)>(b)? (a) : (b))

/** Mimimum of its two arguments.
 */
#define NFFT_MIN(a,b) ((a)<(b)? (a) : (b))

/* ######################################################################### */
/* ########## Window function related ###################################### */
/* ######################################################################### */

/** To test the new one
 */
double nfft_bspline_old(int k,double x,double *A);

/** Computes the B-spline \f$M_{k,0}\left(x\right)\f$,
    scratch is used for de Boor's scheme
 */
double nfft_bspline(int k, double x, double *scratch);

/* ######################################################################### */
/* ########## Vector routines ############################################## */
/* ######################################################################### */

/** Computes integer \f$\prod_{t=0}^{d-1} v_t\f$.
 */
int nfft_prod_int(int *vec, int d);

/** Computes integer \f$\prod_{t=0}^{d-1} v_t-a\f$.
 */
int nfst_prod_minus_a_int(int *vec, int a, int d);

/** Computes \f$\sum_{t=0}^{d-1} i_t \prod_{t'=t+1}^{d-1} N_{t'}\f$.
 */
int nfft_plain_loop(int *idx,int *N,int d);

/** Computes double \f$\prod_{t=0}^{d-1} v_t\f$.
 */
double nfft_prod_real(double *vec,int d);

/** Computes the inner/dot product \f$x^H x\f$.
 */
double nfft_dot_complex(fftw_complex* x, int n);

/** Computes the inner/dot product \f$x^H x\f$.
 */
double nfft_dot_double( double*  x, int n);

/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$.
 */
double nfft_dot_w_complex(fftw_complex* x, double* w, int n);

/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$.
 */
double nfft_dot_w_double( double*  x, double* w, int n);

/** Computes the weighted inner/dot product
    \f$x^H (w1\odot w2\odot w2 \odot x)\f$.
*/
double nfft_dot_w_w2_complex(fftw_complex* x, double* w, double* w2, int n);

/** Computes the weighted inner/dot product
    \f$x^H (w2\odot w2 \odot x)\f$.
 */
double nfft_dot_w2_complex(fftw_complex* x, double* w2, int n);

/** Copies \f$x \leftarrow y\f$.
 */
void nfft_cp_complex(fftw_complex* x, fftw_complex* y, int n);

/** Copies \f$x \leftarrow y\f$.
 */
void nfft_cp_double( double*  x, double*  y, int n);

/** Copies \f$x \leftarrow a y\f$.
 */
void nfft_cp_a_complex(fftw_complex* x, double a, fftw_complex* y, int n);

/** Copies \f$x \leftarrow a y\f$.
 */
void nfft_cp_a_double(double *x, double a, double *y, int n);

/** Copies \f$x \leftarrow w\odot y\f$.
 */
void nfft_cp_w_complex(fftw_complex* x, double* w, fftw_complex* y, int n);

/** Copies \f$x \leftarrow w\odot y\f$.
 */
void nfft_cp_w_double( double*  x, double* w, double*  y, int n);

/** Updates \f$x \leftarrow a x + y\f$.
 */
void nfft_upd_axpy_complex(fftw_complex* x, double a, fftw_complex* y, int n);

/** Updates \f$x \leftarrow a x + y\f$.
 */
void nfft_upd_axpy_double( double*  x, double a, double*  y, int n);

/** Updates \f$x \leftarrow x + a y\f$.
 */
void nfft_upd_xpay_complex(fftw_complex* x, double a, fftw_complex* y, int n);

/** Updates \f$x \leftarrow x + a y\f$.
 */
void nfft_upd_xpay_double( double*  x, double a, double*  y, int n);

/** Updates \f$x \leftarrow a x + b y\f$.
 */
void nfft_upd_axpby_complex(fftw_complex* x, double a, fftw_complex* y, double b, int n);

/** Updates \f$x \leftarrow a x + b y\f$.
 */
void nfft_upd_axpby_double(  double* x, double a, double*  y, double b, int n);

/** Updates \f$x \leftarrow x + a w\odot y\f$.
 */
void nfft_upd_xpawy_complex(fftw_complex* x, double a, double* w, fftw_complex* y, int n);

/** Updates \f$x \leftarrow x + a w\odot y\f$.
 */
void nfft_upd_xpawy_double( double*  x, double a, double* w, double*  y, int n);

/** Updates \f$x \leftarrow a x +  w\odot y\f$.
 */
void nfft_upd_axpwy_complex(fftw_complex* x, double a, double* w, fftw_complex* y, int n);

/** Updates \f$x \leftarrow a x +  w\odot y\f$.
 */
void nfft_upd_axpwy_double( double*  x, double a, double* w, double*  y, int n);

/** Swaps each half over N[d]/2.
 */
void nfft_fftshift_complex(fftw_complex *x, int d, int* N);

/** Prints a vector of integer numbers.
 */
void nfft_vpr_int(int *x, int n, char *text);

/** Prints a vector of doubles numbers.
 */
void nfft_vpr_double(double *x, int n, const char *text);

/** Prints a vector of complex numbers.
 */
void nfft_vpr_complex(fftw_complex *x, int n, const char *text);

/** Inits a vector of random complex numbers in \f$[0,1]\times[0,1]{\rm i}\f$.
 */
void nfft_vrand_unit_complex(fftw_complex *x, int n);

/** Inits a vector of random double numbers in \f$[-1/2,1/2]\f$.
 */
void nfft_vrand_shifted_unit_double(double *x, int n);

/* ######################################################################### */
/* ########## Helpers for inverse transforms ############################### */
/* ######################################################################### */

/** Computes non periodic voronoi weights, assumes ordered nodes \f$x_j\f$ */
void nfft_voronoi_weights_1d(double *w, double *x, int M);

/** Computes voronoi weights for nodes on the sphere S^2. */
void nfft_voronoi_weights_S2(double *w, double *xi, int M);

/** Computes the damping factor for the modified Fejer kernel, ie
    \f$\frac{2}{N}\left(1-\frac{\left|2k+1\right|}{N}\right)\f$
 */
double nfft_modified_fejer(int N,int kk);

/** Computes the damping factor for the modified Jackson kernel.
 */
double nfft_modified_jackson2(int N,int kk);

/** Computes the damping factor for the modified generalised Jackson kernel.
 */
double nfft_modified_jackson4(int N,int kk);

/** Computes the damping factor for the modified Sobolev kernel.
 */
double nfft_modified_sobolev(double mu,int kk);

/** Computes the damping factor for the modified multiquadric kernel.
 */
double nfft_modified_multiquadric(double mu,double c,int kk);

int nfft_smbi(const double x, const double alpha, const int nb, const int ize,
  double *b);

double nfft_drand48(void);

void nfft_srand48(long int seed);

/** Radix sort for node indices.
 */
void nfft_sort_node_indices_radix_lsdf(int n, int *keys0, int *keys1, int rhigh);

/** Radix sort for node indices.
 */
void nfft_sort_node_indices_radix_msdf(int n, int *keys0, int *keys1, int rhigh);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

/** @}
 */
#endif
