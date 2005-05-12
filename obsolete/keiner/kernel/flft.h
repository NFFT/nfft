#ifndef FLFT_H
#define FLFT_H

#include "api.h"

/**
 * \brief Fast Legendre function transform
 *
 * Computes the fast Legendre function transform (FLFT) for a given order 
 * \f$n\f$.
 *
 * \arg M The bandwidth \f$M\f$
 * \arg t The exponent of the next greater power of 2 relative to M
 * \arg n The order \f$n\f$
 * \arg f_hat The vector containing the Fourier coeffcients with f_hat[k] = 
 *   \f$a_k^n (k=|n|,\dots,M)\f$
 * \arg U The three-dimensional array containing the matrices \f$U_{n,\tau,l}\f$
 * \arg tw The transform wisdom
 */
void flft(const int M, const int t, const int n, complex *f_hat, 
          const struct U_type ***U, const struct nfsft_transform_wisdom *tw);

/**
 * \brief AdjointfFast Legendre function transform
 *
 * Computes the adjoint fast Legendre function transform (FLFT) for a given 
 * order \f$n\f$.
 *
 * \arg M The bandwidth \f$M\f$
 * \arg t The exponent of the next greater power of 2 relative to M
 * \arg n The order \f$n\f$
 * \arg f_hat The vector containing the Fourier coeffcients with f_hat[k] = 
 *   \f$a_k^n (k=|n|,\dots,M)\f$
 * \arg U The three-dimensional array containing the matrices \f$U_{n,\tau,l}\f$
 * \arg tw The transform wisdom
 */
void flft_adjoint(const int M, const int t, const int n, complex *f_hat, const 
                  struct U_type ***U, const struct nfsft_transform_wisdom *tw);

#endif