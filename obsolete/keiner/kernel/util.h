#ifndef UTIL_H
#define UTIL_H

#include <complex.h>

/**
 * Compute the exponent of the next greater power of two.
 *
 * \param An integer n
 *
 * \return The exponent of the next greater power of two with respect to n
 */
inline int ngpt(int n);

double mysecond();
void myvpr (double* x, int n, char * text);
void myvprc (complex *x, int n, char * text);
int pow2(const int t);
double error_complex_inf(complex *x0, complex *x, int n);
double error_complex_inf_r(complex *x0, complex *x, int n);
double error_complex_1(complex *a, complex *b, int n);
double error_complex_2(complex *a, complex *b, int n);
double error_complex3(complex *a, complex *b, int n);
double err_f_hat_infty(complex **f_hat, complex **f_hat2, int M);
double err_f_hat_1(complex **f_hat, complex **f_hat2, int M);
double err_f_hat_2(complex **f_hat, complex **f_hat2, int M);

#define PI 3.1415926535897932385
#define QUADRATIC_KERNEL(k) k==0?1.0:((2*k+1)/((double)(k*(k+1))*(k*(k+1))))
#define ABEL_POISSON_KERNEL(k,h) (2*k+1)*pow(h,k)*((1-h)*(1-h)/(1+h))
#define GAUSS_WEIERSTRASS(k,p) exp(-k*(k+1)*p)*(2*k+1)/(4*PI) 
#define max(a,b) a<b?b:a


#endif
