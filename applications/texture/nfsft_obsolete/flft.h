#ifndef FLFT_H
#define FLFT_H

/* Internal API */
#include "api_old.h"

/**
 * \brief Fast Legendre function transform
 *
 * Computes the fast Legendre function transform (FLFT) for a given order 
 * \f$n\f$.
 *
 * \arg M The bandwidth \f$M \in \mathbb{N}, M \ge 4\f$
 * \arg t The exponent of the next greater power of 2 relative to M, hence 
 *        \f$t := \left\lceil \log_2 M \right\rceil\f$
 * \arg n The order \f$n, -M \le n \le M\f$
 * \arg f_hat The vector containing the Fourier coeffcients with f_hat[n+M][k] = 
 *        \f$a_k^n \quad (k=|n|,\dots,M)\f$
 * \arg wisdom The wisdom structure with precomputed data
 * 
 * 
 */
void flft_old(const int M, const int t, const int n, complex *const f_hat, 
          struct nfsft_wisdom_old *const wisdom, int *nstab, int *ntotal);

/**
 * \brief Adjoint fast Legendre function transform
 *
 * Computes the adjoint fast Legendre function transform (FLFT) for a given 
 * order \f$n\f$.
 *
 * \arg M The bandwidth \f$M \in \mathbb{N}, M \ge 4\f$
 * \arg t The exponent of the next greater power of 2 relative to M, hence 
 *        \f$t := \left\lceil \log_2 M \right\rceil\f$
 * \arg n The order \f$n, -M \le n \le M\f$
 * \arg f_hat The vector containing the Chebyshev coeffcients with 
 *        f_hat[n+M][k] = \f$b_k^n \quad (k=0\dots,M)\f$
 * \arg wisdom The wisdom structure with precomputed data
 * 
 * 
 */
void flft_adjoint_old(const int M, const int t, const int n, complex *const f_hat, 
                  struct nfsft_wisdom_old *const wisdom);

void flft_stab_old(const int M, const int t, const int n, complex *const f_hat, 
               struct nfsft_wisdom_old *const wisdom, int *nstab, int *ntotal);

#endif
