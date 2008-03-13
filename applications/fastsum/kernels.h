/*! \file kernels.h
 *  \brief Header file with predefined kernels for the fast summation algorithm.
 */
#ifndef KERNELS_H
#define KERNELS_H

#include <complex.h>

/** 
 * \addtogroup applications_fastsum
 * \{
 */

double _Complex gaussian(double x, int der, const double *param);              /* K(x)=exp(-x^2/c^2) */
double _Complex multiquadric(double x, int der, const double *param);          /* K(x)=sqrt(x^2+c^2) */
double _Complex inverse_multiquadric(double x, int der, const double *param);  /* K(x)=1/sqrt(x^2+c^2) */
double _Complex logarithm(double x, int der, const double *param);             /* K(x)=log |x| */
double _Complex thinplate_spline(double x, int der, const double *param);      /* K(x) = x^2 log |x| */
double _Complex one_over_square(double x, int der, const double *param);       /* K(x) = 1/x^2 */
double _Complex one_over_modulus(double x, int der, const double *param);      /* K(x) = 1/|x| */
double _Complex one_over_x(double x, int der, const double *param);            /* K(x) = 1/x */
double _Complex inverse_multiquadric3(double x, int der, const double *param); /* K(x) = 1/sqrt(x^2+c^2)^3 */
double _Complex sinc_kernel(double x, int der, const double *param);           /* K(x) = sin(cx)/x */
double _Complex cosc(double x, int der, const double *param);                  /* K(x) = cos(cx)/x */
double _Complex kcot(double x, int der, const double *param);                   /* K(x) = cot(cx) */
double _Complex one_over_cube(double x, int der, const double *param);                /* K(x) = 1/x^3 */
/* \} */

#endif
/* kernels.h */
