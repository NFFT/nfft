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

#define SWAPCT(x,y,temp) {temp=(x); (x)=(y); (y)=temp;}

void vpr_c_hat (complex **x, int n, char * text);

double error_complex_inf(complex *x0, complex *x, int n);
double error_complex_inf_r(complex *x0, complex *x, int n);
double error_complex_1(complex *a, complex *b, int n);
double error_complex_2(complex *a, complex *b, int n);
double error_complex3(complex *a, complex *b, int n);
double err_f_hat_infty(complex **f_hat, complex **f_hat2, int M);
double err_f_hat_1(complex **f_hat, complex **f_hat2, int M);
double err_f_hat_2(complex **f_hat, complex **f_hat2, int M);
double norm_complex_inf(complex *x, int n);
double norm_complex_1(complex *x, int n);
double norm_complex_2(complex *x, int n);
double norm_f_hat_1(complex **f_hat, int M);
int min(const int a, const int b);
int max(const int a, const int b);
long int maxl(const long int a, const long int b);

void copyc_hat(complex **x,complex **y, int M);
void copyc_w_hat(complex **x, double **w, complex **y, int M);
double dotproductc_w_hat(complex **x, double **w, int M);
double dotproductc_hat(complex** x, int M);
void updatec_xpay_hat(complex **x,double a, complex **y, int M);
void updatec_xpawy_hat(complex **x, double a, double **w, complex **y, int M);
void updatec_axpy_hat(complex **x,double a, complex **y, int M);
void updatec_axpby_hat(complex **x, double a, complex **y, double b, int M);
void updatec_axpy_2(complex* x,double a, complex* y, int n);
void updatec_xpay_2(complex* x,double a, complex* y, int n);
void updatec_axpby_2(complex* x, double a, complex* y, double b, int n);
void updatec_xpawy_2(complex* x,double a, double* w, complex* y, int n);
void uvxpwy(complex* u, complex* x, double* v, complex* y, double* w, int n);
void auvxpwy(double a, complex* u, complex* x, double* v, complex* y, double* w, int n);
void abuvxpwy(double a, double b, complex* u, complex* x, double* v, complex* y, double* w, int n);
void normalize_f_hat(complex **f_hat, int M);

#endif
