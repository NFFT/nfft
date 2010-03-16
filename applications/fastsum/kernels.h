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

/*! \file kernels.h
 *  \brief Header file with predefined kernels for the fast summation algorithm.
 */
#ifndef KERNELS_H
#define KERNELS_H

#include <complex.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

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

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
/* kernels.h */
