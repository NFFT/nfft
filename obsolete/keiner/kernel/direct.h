/**
 * Header file for direct NDSFT algorithms
 */
#ifndef DIRECT_H
#define DIRECT_H

#include "api.h"

/** \defgroup nfsft_direct NFSFT: Direct algorithms */

/*
 * Spherical Fourier transform and adjoint spherical Fourier transform 
 * implemented by means of the direct algorithms.
 */

/**
 * \brief Direct spherical Fourier transform
 *
 * Direct implementation of the forward spherical Fourier transform. Given nodes
 * \f$(\vartheta_d,\varphi_d)_{d=0}^{D-1}\f$ (\f$D \in \mathbb{N}\f$) and 
 * Fourier-coefficients \f$(a_k^n)\f$, (\f$k = 0,...,M\f$, \f$n = -k,...,k\f$) 
 * the function evaluates the Fourier-sum at the specified nodes.
 *
 * \param D The number of nodes
 * \param angles Pointer to an array containing the nodes 
 *   \f$(\vartheta_d,\varphi_d)_{d=0}^{D-1}\f$ (in spherical coordinates) 
 *   ordered as \f$\varphi_0,\vartheta_0,\varphi_1,\dots,\vartheta_{D-1}\f$
 * \param f Pointer to an array where the function values will be stored
 * \param M The bandwidth
 * \param N Next greater power of two with respect to the bandwidth \f$M\f$
 * \param f_hat Pointer to an array containing the spherical Fourier 
 *  coefficients \f$(a_k^n)\f$, \f$k = 0,...,M\f$, \f$n = -k,...,k\f$ in the 
 *  ordering \f$a_0^0,a_1^{-1},a_1^0,...,a_M^M\f$
 * \param wisdom Structure containing precomputed values of associated Legendre 
 *   functions.
 *
 * \ingroup nfsft_direct
 */
void ndsft(int D, double *angles, complex *f, int M, complex **f_hat, 
  struct nfsft_wisdom *wisdom);
	
/**
 * \brief Adjoint direct spherical Fourier transform
 *
 * Direct implementation of the adjoint forward spherical Fourier transform. The
 * spherical Fourier transform can be regarded as a multiplication of a vector
 * \f$a \in \mathbb{C}^{(M+1)^2}\f$ (\f$M \in \mathbb{N}\f$) containing Fourier-
 * coefficents with a matrix \f$Y \in \mathbb{C}^{D \times (M+1)^2}\f$ 
 * (\f$D \in \mathbb{N}\f$) resulting in a vector \f$f \in \mathbb{C}^{D}\f$ 
 * containing the Fourier-sum evaluated at nodes 
 * \f$(\vartheta_d,\varphi_d)_{d=0}^{D-1}\f$. This function
 * implements the matrix-vector-multiplication with the adjoint matrix 
 * \f$Y^{\text{H}}\f$, hence 
 * \f[\tilde{\mathbf{a}} = Y^{\text{H}} \tilde{\mathbf{f}}.\f]
 *
 * \param angles Array containing the "nodes" 
 *   \f$(\vartheta_d,\varphi_d)_{d=0}^{D-1}\f$ (in spherical coordinates) 
 *   aligned as \f$\varphi_0,\vartheta_0,\varphi_1,\dots,\vartheta_{D-1}\f$
 * \param D The number of "nodes"
 * \param f_hat Array containing the "Fourier coefficients" 
 *  \f$(a_k^n)\f$, \f$k = 0,...,M\f$, \f$n = -k,...,k\f$ in the ordering 
 *  \f$a_0^0,a_1^{-1},a_1^0,...,a_M^M\f$
 * \param M The bandwidth
 * \param N Next greater power of two with respect to the bandwidth \f$M\f$
 * \param wisdom Structure containing precomputed values of associated Legendre 
 *   functions.
 *
 * \ingroup nfsft_direct
 */	
void adjoint_ndsft(int D, double *angles, complex *f, int M, complex **f_hat, 
  struct nfsft_wisdom *wisdom);
#endif
