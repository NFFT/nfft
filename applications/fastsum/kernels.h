/*! \file kernels.h
 *  \brief Header file with predefined kernels for the fast summation algorithm.
 */
#ifndef KERNELS_H
#define KERNELS_H

/** 
 * \addtogroup applications_fastsum
 * \{
 */

complex double gaussian(double x, int der, const double *param);              /* K(x)=exp(-x^2/c^2) */
complex double multiquadric(double x, int der, const double *param);          /* K(x)=sqrt(x^2+c^2) */
complex double inverse_multiquadric(double x, int der, const double *param);  /* K(x)=1/sqrt(x^2+c^2) */
complex double logarithm(double x, int der, const double *param);             /* K(x)=log |x| */
complex double thinplate_spline(double x, int der, const double *param);      /* K(x) = x^2 log |x| */
complex double one_over_square(double x, int der, const double *param);       /* K(x) = 1/x^2 */
complex double one_over_modulus(double x, int der, const double *param);      /* K(x) = 1/|x| */
complex double one_over_x(double x, int der, const double *param);            /* K(x) = 1/x */
complex double inverse_multiquadric3(double x, int der, const double *param); /* K(x) = 1/sqrt(x^2+c^2)^3 */
complex double sinc_kernel(double x, int der, const double *param);           /* K(x) = sin(cx)/x */
complex double cosc(double x, int der, const double *param);                  /* K(x) = cos(cx)/x */
complex double cot(double x, int der, const double *param);                   /* K(x) = cot(cx) */
complex one_over_cube(double x, int der, const double *param);                /* K(x) = 1/x^3 */
/* \} */

#endif
/* kernels.h */
