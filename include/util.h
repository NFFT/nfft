/** Header for utilities.
 *  functions for vectors, window functions, ...
 */

#ifndef utils_h_inc
#define utils_h_inc

#include <complex.h>

/** Macros
 */
#define SWAPC(x,y) {complex* temp; temp=(x); (x)=(y); (y)=temp;}
#define SWAP_double(x,y) {double* temp; temp=(x); (x)=(y); (y)=temp;}
#define SWAP_complex(x,y) {complex* temp; temp=(x); (x)=(y); (y)=temp;}
#define PI 3.1415926535897932384
#define MAX(a,b) ((a)>(b)? (a) : (b))
#define MIN(a,b) ((a)<(b)? (a) : (b))

/** @defgroup little_group Group for little helpers
 * This group contains little helpers, e.g. for
 * timing, ...
 * @{ 
 */

/** Actual used CPU time in seconds.
 *  Calls getrusage, limited accuracy
 */
double second();

/** Computes /f$n\ge N/f$ such that /f$n=2^j,\, j\in\mathhb{N}_0/f$.
 */
int next_power_of_2(int N);

/** Computes integer /f$\prod_{t=0}^{d-1} v_t/f$.
 */
int nfft_prod_int(int *vec, int d);
int nfct_prod_int(int *vec, int d);

/** Computes integer /f$\prod_{t=0}^{d-1} v_t-a/f$.
 */
int nfst_prod_minus_a_int(int *vec, int a, int d);

/** Computes /f$\sum_{t=0}^{d-1} i_t \prod_{t'=t+1}^{d-1} N_{t'}/f$.
 */
int nfft_plain_loop(int *idx,int *N,int d);

/** Computes double /f$\prod_{t=0}^{d-1} v_t/f$.
 */
double nfft_prod_real(double *vec,int d);

/** @} 
 */ 

/** @defgroup window_group Group for window function routines
 * This group contains routines regarding window functions
 * @{ 
 */

/** Sinus cardinalis
 *  \f$\frac{sin\left(x\right)}{x}$ 
 */
double sinc(double x);

/** To test the new one
 */
double bspline_old(int k,double x,double *A);

/** Computes \f$M_{k,0}\left(x\right)\f$
 *  scratch is used for de Boor's scheme
 */
double bspline(int k, double x, double *scratch);

/** Modified Bessel function of order zero.
 *  Cephes Math Library Release 2.8:  June, 2000
 *  Copyright 1984, 1987, 2000 by Stephen L. Moshier
 */
double i0(double x);

/** @} 
 */ 

/** @defgroup vector_group Group for vector routines
 * This group contains routines for handling vectors
 * @{ 
 */

/** Computes the inner/dot product \f$x^H x\f$.
 */
double dot_complex(complex* x, int n);
double dot_double( double*  x, int n);

/** Computes the weighted inner/dot product \f$x^H (w \odot x)\f$.
 */
double dot_w_complex(complex* x, double* w, int n);
double dot_w_double( double*  x, double* w, int n);

/** Computes the weighted inner/dot product 
    \f$x^H (w1\odot w2\odot w2 \odot x)\f$.
 */
double dot_w_w2_complex(complex* x, double* w, double* w2, int n);

/** Computes the weighted inner/dot product 
    \f$x^H (w2\odot w2 \odot x)\f$.
 */
double dot_w2_complex(complex* x, double* w2, int n);

/** Copies \f$x \leftarrow y\f$.
 */
void cp_complex(complex* x, complex* y, int n);
void cp_double( double*  x, double*  y, int n);

/** Copies \f$x \leftarrow a y\f$.
 */
void cp_a_complex(complex* x, double a, complex* y, int n);

/** Copies \f$x \leftarrow w\odot y\f$.
 */
void cp_w_complex(complex* x, double* w, complex* y, int n);
void cp_w_double( double*  x, double* w, double*  y, int n);

/** Updates \f$x \leftarrow a x + y\f$.
 */
void upd_axpy_complex(complex* x, double a, complex* y, int n);
void upd_axpy_double( double*  x, double a, double*  y, int n);

/** Updates \f$x \leftarrow x + a y\f$.
 */
void upd_xpay_complex(complex* x, double a, complex* y, int n);
void upd_xpay_double( double*  x, double a, double*  y, int n);

/** Updates \f$x \leftarrow a x + b y\f$.
 */
void upd_axpby_complex(complex* x, double a, complex* y, double b, int n);
void upd_axpby_double(  double* x, double a, double*  y, double b, int n);

/** Updates \f$x \leftarrow x + a w\odot y\f$.
 */
void upd_xpawy_complex(complex* x, double a, double* w, complex* y, int n);
void upd_xpawy_double( double*  x, double a, double* w, double*  y, int n);

/** Updates \f$x \leftarrow a x +  w\odot y\f$.
 */
void upd_axpwy_complex(complex* x, double a, double* w, complex* y, int n);
void upd_axpwy_double( double*  x, double a, double* w, double*  y, int n);

/** Swaps each half over N[d]/2.
 */
void fftshift_complex(complex *x, int d, int* N);

/** computes \f$\frac{\|x-y\|_{\infty}}{\|x\|_{\infty}} \f$
 */
double error_l_infty_complex(complex *x, complex *y, int n);

/** computes \f$\frac{\|x-y\|_{\infty}}{\|x\|_{\infty}} \f$
 */
double error_l_infty_double(double *x, double *y, int n);

/** computes \f$\frac{\|x-y\|_{\infty}}{\|z\|_1} \f$
 */
double error_l_infty_1_complex(complex *x, complex *y, int n, complex *z,
                               int m);

/** computes \f$\frac{\|x-y\|_{\infty}}{\|z\|_1} \f$
 */
double error_l_infty_1_double(double *x, double *y, int n, double *z,
                               int m);

/** computes \f$\frac{\|x-y\|_2}{\|x\|_2} \f$
 */
double error_l_2_complex(complex *x, complex *y, int n);

/** computes \f$\frac{\|x-y\|_2}{\|x\|_2} \f$
 */
double error_l_2_double(double *x, double *y, int n);


void vpr_int(int *x, int n, char *text);

void vpr_double(double *x, int n, char *text);

void vpr_complex(complex *x, int n, char *text);

void vrand_unit_complex(complex *x, int n);

void vrand_shifted_unit_double(double *x, int n);
/** @} 
 */ 

/** @defgroup dampinggroup Group for damping factors in iterative reconstruction
 * This group contains sub routines for different one dimensional damping
 * factors.
 * @{ 
 */

/** Computes non periodic voronoi weights 
 *  assumes ordered x_j */
void voronoi_weights_1d(double *w, double *x, int M);

/** Computes the damping factor for the modified Fejer kernel.
 *  /f$\frac{2}{N}\left(1-\frac{\left|2k+1\right|}{N}\right)/f$
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

#ifdef HAVE_MALLINFO
#define HAVE_TOTAL_USED_MEMORY
/**
 * 
 */
int total_used_memory();
#endif
/** @} 
 */ 
#endif
