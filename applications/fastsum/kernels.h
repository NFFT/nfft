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

#include <complex.h>

#include "nfft3.h"
#include "nfft3mp.h"
#include "nfft3util.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/**
 * \addtogroup applications_fastsum
 * \{
 */

NFFT_C gaussian(NFFT_R x, int der, const NFFT_R *param);              /**< K(x)=exp(-x^2/c^2) */
NFFT_C multiquadric(NFFT_R x, int der, const NFFT_R *param);          /**< K(x)=sqrt(x^2+c^2) */
NFFT_C inverse_multiquadric(NFFT_R x, int der, const NFFT_R *param);  /**< K(x)=1/sqrt(x^2+c^2) */
NFFT_C logarithm(NFFT_R x, int der, const NFFT_R *param);             /**< K(x)=log |x| */
NFFT_C thinplate_spline(NFFT_R x, int der, const NFFT_R *param);      /**< K(x) = x^2 log |x| */
NFFT_C one_over_square(NFFT_R x, int der, const NFFT_R *param);       /**< K(x) = 1/x^2 */
NFFT_C one_over_modulus(NFFT_R x, int der, const NFFT_R *param);      /**< K(x) = 1/|x| */
NFFT_C one_over_x(NFFT_R x, int der, const NFFT_R *param);            /**< K(x) = 1/x */
NFFT_C inverse_multiquadric3(NFFT_R x, int der, const NFFT_R *param); /**< K(x) = 1/sqrt(x^2+c^2)^3 */
NFFT_C sinc_kernel(NFFT_R x, int der, const NFFT_R *param);           /**< K(x) = sin(cx)/x */
NFFT_C cosc(NFFT_R x, int der, const NFFT_R *param);                  /**< K(x) = cos(cx)/x */
NFFT_C kcot(NFFT_R x, int der, const NFFT_R *param);                  /**< K(x) = cot(cx) */
NFFT_C one_over_cube(NFFT_R x, int der, const NFFT_R *param);         /**< K(x) = 1/x^3 */
NFFT_C log_sin(NFFT_R x, int der, const NFFT_R *param);               /**< K(x) = log(|sin(cx)|) */
NFFT_C laplacian_rbf(NFFT_R x, int der, const NFFT_R *param);         /**< K(x) = exp(-|x|/c) */
NFFT_C der_laplacian_rbf(NFFT_R x, int der, const NFFT_R *param);     /**< K(x) = |x|/c exp(-|x|/c) */
NFFT_C xx_gaussian(NFFT_R x, int der, const NFFT_R *param);           /**< K(x) = x^2/c^2 exp(-x^2/c^2) */
NFFT_C absx(NFFT_R x, int der, const NFFT_R *param);                  /**< K(x) = |x| */
/* \} */

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
/* kernels.h */
