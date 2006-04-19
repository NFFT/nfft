/*! \file kernels.h
 *  \brief Header file with predefined kernels for the fast summation algorithm.
 */
#ifndef KERNELS_H
#define KERNELS_H
//#include "kernels.c"

complex gaussian(double x, int der, const double *param);              /* K(x)=exp(-x^2/c^2) */
complex multiquadric(double x, int der, const double *param);          /* K(x)=sqrt(x^2+c^2) */
complex inverse_multiquadric(double x, int der, const double *param);  /* K(x)=1/sqrt(x^2+c^2) */
complex logarithm(double x, int der, const double *param);             /* K(x)=log |x| */
complex thinplate_spline(double x, int der, const double *param);      /* K(x) = x^2 log |x| */
complex one_over_square(double x, int der, const double *param);       /* K(x) = 1/x^2 */
complex one_over_modulus(double x, int der, const double *param);      /* K(x) = 1/|x| */
complex one_over_x(double x, int der, const double *param);            /* K(x) = 1/x */
complex inverse_multiquadric3(double x, int der, const double *param); /* K(x) = 1/sqrt(x^2+c^2)^3 */
complex sinc_kernel(double x, int der, const double *param);           /* K(x) = sin(cx)/x */
complex cosc(double x, int der, const double *param);                  /* K(x) = cos(cx)/x */
complex cot(double x, int der, const double *param);                   /* K(x) = cot(cx) */

#endif
/* kernels.h */
