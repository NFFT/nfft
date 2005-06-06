/** 
 * \file c2f.h
 * \brief 
 * \author Jens Keiner
 */
#ifndef C2F_H
#define C2F_H

#include <complex.h>

/**
 * Converts Chebyshev coefficients to Fourier coefficients.
 * 
 * \arg f_hat Fourier coeffcients
 * \arg cheb Chebyshev coeffcients
 * \arg M Bandwidth 
 * \arg N
 *
 * \author Jens Keiner
 * \ingroup nfsft_internal
 */
void cheb2exp(complex *f_hat, complex **cheb, int M, int N);

/**
 * Adjoint algorithm for cheb2exp.
 * 
 * \arg f_hat Fourier coeffcients
 * \arg cheb Chebyshev coeffcients
 * \arg M Bandwidth 
 * \arg N
 *
 * \author Jens Keiner
 * \ingroup nfsft_internal
 */
void cheb2exp_adjoint(complex *f_hat, complex **cheb, int M, int N);

#endif
