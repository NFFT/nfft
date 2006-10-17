/*! \file util.h
 *  \brief Header file for utility functions used by the nfft3 library.
 */
#ifndef utils_h_inc
#define utils_h_inc

/** Include header for C99 complex datatype. */
#include <complex.h>

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/**
 * @defgroup nfftutil UTIL - Utility functions for the NFFT library
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
#define SWAP_complex(x,y) {double complex* temp; temp=(x); (x)=(y); (y)=temp;}

/** Swapping of two vectors.
 */
#define SWAP_double(x,y) {double* temp; temp=(x); (x)=(y); (y)=temp;}

/** Formerly known to be an irrational number.
 */
#define PI 3.1415926535897932384

/** Maximum of its two arguments.
 */
#define MAX(a,b) ((a)>(b)? (a) : (b))

/** Mimimum of its two arguments.
 */
#define MIN(a,b) ((a)<(b)? (a) : (b))

/* ######################################################################### */
/* ########## Little helpers ############################################### */
/* ######################################################################### */

/** Actual used CPU time in seconds; calls getrusage, limited accuracy.
 */
double second();

/** Actual used memory in bytes; calls mallinfo if define HAVE_MALLOC_H.
 */
int total_used_memory();

/** Integer logarithm of 2.
 */
int ld(int m);

/** Integer power of 2.
 */
int int_2_pow(int a);

/** Computes \f$n\ge N\f$ such that \f$n=2^j,\, j\in\mathbb{N}_0\f$.
 */
int next_power_of_2(int N);

/** Computes ?
 */
void next_power_of_2_exp(int N, int *N2, int *t);

/* ######################################################################### */
/* ########## Window function related ###################################### */
/* ######################################################################### */

/** Computes the sinus cardinalis \f$\frac{sin\left(x\right)}{x}\f$.
 */
double sinc(double x);

/** To test the new one
 */
double bspline_old(int k,double x,double *A);

/** Computes the B-spline \f$M_{k,0}\left(x\right)\f$, 
    scratch is used for de Boor's scheme
 */
double bspline(int k, double x, double *scratch);

/** Modified Bessel function of order zero; 
    adapted from Stephen Moshier's Cephes Math Library Release 2.8
 */
double i0(double x);

/* ######################################################################### */
/* ########## Vector routines ############################################## */
/* ######################################################################### */

/** Computes integer \f$\prod_{t=0}^{d-1} v_t\f$.
 */
int nfft_prod_int(int *vec, int d);

/** Computes integer \f$\prod_{t=0}^{d-1} v_t\f$.
 */
int nfct_prod_int(int *vec, int d);

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
double dot_complex(double complex* x, int n);

/** Computes the inner/dot product \f$x^H x\f$.
 */
double dot_double( double*  x, int n);

/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$.
 */
double dot_w_complex(double complex* x, double* w, int n);

/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$.
 */
double dot_w_double( double*  x, double* w, int n);

/** Computes the weighted inner/dot product
    \f$x^H (w1\odot w2\odot w2 \odot x)\f$.
*/
double dot_w_w2_complex(double complex* x, double* w, double* w2, int n);

/** Computes the weighted inner/dot product
    \f$x^H (w2\odot w2 \odot x)\f$.
 */
double dot_w2_complex(double complex* x, double* w2, int n);

/** Copies \f$x \leftarrow y\f$.
 */
void cp_complex(double complex* x, double complex* y, int n);

/** Copies \f$x \leftarrow y\f$.
 */
void cp_double( double*  x, double*  y, int n);

/** Copies \f$x \leftarrow a y\f$.
 */
void cp_a_complex(double complex* x, double a, double complex* y, int n);

/** Copies \f$x \leftarrow w\odot y\f$.
 */
void cp_w_complex(double complex* x, double* w, double complex* y, int n);

/** Copies \f$x \leftarrow w\odot y\f$.
 */
void cp_w_double( double*  x, double* w, double*  y, int n);

/** Updates \f$x \leftarrow a x + y\f$.
 */
void upd_axpy_complex(double complex* x, double a, double complex* y, int n);

/** Updates \f$x \leftarrow a x + y\f$.
 */
void upd_axpy_double( double*  x, double a, double*  y, int n);

/** Updates \f$x \leftarrow x + a y\f$.
 */
void upd_xpay_complex(double complex* x, double a, double complex* y, int n);

/** Updates \f$x \leftarrow x + a y\f$.
 */
void upd_xpay_double( double*  x, double a, double*  y, int n);

/** Updates \f$x \leftarrow a x + b y\f$.
 */
void upd_axpby_complex(double complex* x, double a, double complex* y, double b, int n);

/** Updates \f$x \leftarrow a x + b y\f$.
 */
void upd_axpby_double(  double* x, double a, double*  y, double b, int n);

/** Updates \f$x \leftarrow x + a w\odot y\f$.
 */
void upd_xpawy_complex(double complex* x, double a, double* w, double complex* y, int n);

/** Updates \f$x \leftarrow x + a w\odot y\f$.
 */
void upd_xpawy_double( double*  x, double a, double* w, double*  y, int n);

/** Updates \f$x \leftarrow a x +  w\odot y\f$.
 */
void upd_axpwy_complex(double complex* x, double a, double* w, double complex* y, int n);

/** Updates \f$x \leftarrow a x +  w\odot y\f$.
 */
void upd_axpwy_double( double*  x, double a, double* w, double*  y, int n);

/** Swaps each half over N[d]/2.
 */
void fftshift_complex(double complex *x, int d, int* N);

/** Computes \f$\frac{\|x-y\|_{\infty}}{\|x\|_{\infty}} \f$
 */
double error_l_infty_complex(double complex *x, double complex *y, int n);

/** Computes \f$\frac{\|x-y\|_{\infty}}{\|x\|_{\infty}} \f$
 */
double error_l_infty_double(double *x, double *y, int n);

/** Computes \f$\frac{\|x-y\|_{\infty}}{\|z\|_1} \f$
 */
double error_l_infty_1_complex(double complex *x, double complex *y, int n, double complex *z,
                               int m);

/** Computes \f$\frac{\|x-y\|_{\infty}}{\|z\|_1} \f$
 */
double error_l_infty_1_double(double *x, double *y, int n, double *z,
			      int m);

/** Computes \f$\frac{\|x-y\|_2}{\|x\|_2} \f$
 */
double error_l_2_complex(double complex *x, double complex *y, int n);

/** Computes \f$\frac{\|x-y\|_2}{\|x\|_2} \f$
 */
double error_l_2_double(double *x, double *y, int n);

/** Prints a vector of integer numbers.
 */
void vpr_int(int *x, int n, char *text);

/** Prints a vector of doubles numbers.
 */
void vpr_double(double *x, int n, char *text);

/** Prints a vector of complex numbers.
 */
void vpr_complex(double complex *x, int n, char *text);

/** Inits a vector of random complex numbers in \f$[0,1]\times[0,1]{\rm i}\f$.
 */
void vrand_unit_complex(double complex *x, int n);

/** Inits a vector of random double numbers in \f$[-1/2,1/2]\f$.
 */
void vrand_shifted_unit_double(double *x, int n);

/* ######################################################################### */
/* ########## Helpers for inverse transforms ############################### */
/* ######################################################################### */

/** Computes non periodic voronoi weights, assumes ordered nodes \f$x_j\f$ */
void voronoi_weights_1d(double *w, double *x, int M);

/** Computes the damping factor for the modified Fejer kernel, ie
    \f$\frac{2}{N}\left(1-\frac{\left|2k+1\right|}{N}\right)\f$
 */
double modified_fejer(int N,int kk);

/** Computes the damping factor for the modified Jackson kernel.
 */
double modified_jackson2(int N,int kk);

/** Computes the damping factor for the modified generalised Jackson kernel.
 */
double modified_jackson4(int N,int kk);

/** Computes the damping factor for the modified Sobolev kernel.
 */
double modified_sobolev(double mu,int kk);

/** Computes the damping factor for the modified multiquadric kernel.
 */
double modified_multiquadric(double mu,double c,int kk);

/** @}
 */
#endif
