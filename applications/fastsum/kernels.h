/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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

/*! \file kernels.h
 *  \brief Header file with predefined kernels for the fast summation algorithm.
 */
#ifndef KERNELS_H
#define KERNELS_H

#include "config.h"

#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "nfft3.h"
#include "infft.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/**
 * \addtogroup applications_fastsum
 * \{
 */

C gaussian(R x, int der, const R *param);              /**< K(x)=exp(-x^2/c^2) */
C multiquadric(R x, int der, const R *param);          /**< K(x)=sqrt(x^2+c^2) */
C inverse_multiquadric(R x, int der, const R *param);  /**< K(x)=1/sqrt(x^2+c^2) */
C logarithm(R x, int der, const R *param);             /**< K(x)=log |x| */
C thinplate_spline(R x, int der, const R *param);      /**< K(x) = x^2 log |x| */
C one_over_square(R x, int der, const R *param);       /**< K(x) = 1/x^2 */
C one_over_modulus(R x, int der, const R *param);      /**< K(x) = 1/|x| */
C one_over_x(R x, int der, const R *param);            /**< K(x) = 1/x */
C inverse_multiquadric3(R x, int der, const R *param); /**< K(x) = 1/sqrt(x^2+c^2)^3 */
C sinc_kernel(R x, int der, const R *param);           /**< K(x) = sin(cx)/x */
C cosc(R x, int der, const R *param);                  /**< K(x) = cos(cx)/x */
C kcot(R x, int der, const R *param);                  /**< K(x) = cot(cx) */
C one_over_cube(R x, int der, const R *param);         /**< K(x) = 1/x^3 */
C log_sin(R x, int der, const R *param);               /**< K(x) = log(|sin(cx)|) */
C laplacian_rbf(R x, int der, const R *param);         /**< K(x) = exp(-|x|/c) */
/* \} */

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
/* kernels.h */
