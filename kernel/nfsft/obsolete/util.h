#ifndef UTIL_H
#define UTIL_H

#include <complex.h>
#include "../../../include/util.h"

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

double error_complex_inf_old(complex *x0, complex *x, int n);
double error_complex_inf_r_old(complex *x0, complex *x, int n);
double error_complex_1_old(complex *a, complex *b, int n);
double error_complex_2_old(complex *a, complex *b, int n);
double error_complex3_old(complex *a, complex *b, int n);
double err_f_hat_infty_old(complex **f_hat, complex **f_hat2, int M);
double err_f_hat_1_old(complex **f_hat, complex **f_hat2, int M);
double err_f_hat_2_old(complex **f_hat, complex **f_hat2, int M);
double norm_complex_inf_old(complex *x, int n);
double norm_complex_1_old(complex *x, int n);
double norm_complex_2_old(complex *x, int n);
double norm_f_hat_1_old(complex **f_hat, int M);
double norm_nfft_1_old(complex *f_hat, int M);
double norm_f_hat_1s_old(complex **f_hat, int M);
int min_old(const int a, const int b);
int max_old(const int a, const int b);
long int max_oldl(const long int a, const long int b);

void copyc_hat_old(complex **x,complex **y, int M);
void copyc_w_hat_old(complex **x, double **w, complex **y, int M);
double dotproductc_w_hat_old(complex **x, double **w, int M);
double dotproductc_hat_old(complex** x, int M);
void updatec_xpay_hat_old(complex **x,double a, complex **y, int M);
void updatec_xpawy_hat_old(complex **x, double a, double **w, complex **y, int M);
void updatec_axpy_hat_old(complex **x,double a, complex **y, int M);
void updatec_axpby_hat_old(complex **x, double a, complex **y, double b, int M);
void updatec_axpy_2_old(complex* x,double a, complex* y, int n);
void updatec_xpay_2_old(complex* x,double a, complex* y, int n);
void updatec_axpby_2(complex* x, double a, complex* y, double b, int n);
void updatec_xpawy_2_old(complex* x,double a, double* w, complex* y, int n);
void uvxpwy_old(complex* u, complex* x, double* v, complex* y, double* w, int n);
void auvxpwy_old(double a, complex* u, complex* x, double* v, complex* y, double* w, int n);
void abuvxpwy_old(double a, double b, complex* u, complex* x, double* v, complex* y, double* w, int n);
void normalize_f_hat_old(complex **f_hat, int M);
int abs(const int a);
#endif
