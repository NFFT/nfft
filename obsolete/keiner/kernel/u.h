/**
 * Header file for functions related to the matrices \f$U_{n,\tau,l}\f$
 */
#ifndef U_H
#define U_H

/**
 * Precomputation of the matrices \f$U_{n,\tau,l}\f$
 * 
 * \param t Exponent of the maximum bandwidth \f$N = 2^t\f$
 * \param threshold The threshold for stabilization steps
 * 
 * \return 4-dimensional array containing the matrics 
 *
 * \ingroup nfsft_internal
 */
struct U_type**** precomputeU(int t, double threshold, double *walpha, 
                             double *wbeta, double *wgamma, bool window);

/**
 * Forget U
 *
 * \ingroup nfsft_internal
 */
void forgetU(struct U_type**** U, int M, int t, bool window);

/**
 * Fast matrix Multiplication with matrices \f$U_{n,\tau,l}\f$.
 *
 * \param a The first array of Chebyshev-coefficients
 * \param b The second array of Chebyshev-coefficients
 * \param u The \f$2 \times 2\f$-matrix \f$U_{n,\tau,l}\f$
 * \param tau The parameter \f$\tau\f$
 * \param n The parameter n
 * \param l The parameter l
 *
 * \ingroup nfsft_internal
 */
void multiplyU(complex  *a, complex *b, struct U_type u, int tau, int n, int l, 
               struct nfsft_wisdom *tw, double gamma);

/**
 * Fast adjoint matrix Multiplication with matrices \f$U_{n,\tau,l}\f$.
 *
 * \param a The first array of Chebyshev-coefficients
 * \param b The second array of Chebyshev-coefficients
 * \param u The \f$2 \times 2\f$-matrix \f$U_{n,\tau,l}\f$
 * \param tau The parameter \f$\tau\f$
 * \param n The parameter n
 * \param l The parameter l
 *
 * \ingroup nfsft_internal
 */
void multiplyU_adjoint(complex  *a, complex *b, 
                       struct U_type u, int tau, int n, int l, 
                       struct nfsft_wisdom *tw, double gamma);

struct U_type**** precomputeU_guru(int t, double threshold, double *walpha, 
                                   double *wbeta, double *wgamma, int n);

void forgetU_guru(struct U_type**** U, int M, int t, int n);

struct U_type**** precomputeU_stab(int t, double threshold, double *walpha, 
                                   double *wbeta, double *wgamma);

void forgetU_stab(struct U_type**** U, int M, int t);


#endif
