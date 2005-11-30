#ifndef DPT_H_
#define DPT_H_

#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#include <fftw3.h>

/* Flags for dpt_init() */
#define DPT_NO_STABILIZATION  (1U << 0) /**< If set, no stabilization will be 
                                             used.                                 */
#define DPT_BANDWIDTH_WINDOW  (1U << 1) /**< If set, \TODO complete comment here   */
#define DPT_NO_FAST_TRANSFORM (1U << 2) /**< If set, \TODO complete comment here   */
#define DPT_NO_SLOW_TRANSFORM (1U << 3) /**< If set, \TODO complete comment here   */
#define DPT_PERSISTENT_DATA   (1U << 4) /**< If set, \TODO complete comment here   */  

/* Flags for fpt_trafo(), dpt_transposed(), fpt_trafo(), fpt_transposed() */
#define DPT_FUNCTION_VALUES   (1U << 5) /**< If set, the output are function values 
                                             at Chebyshev nodes rather than 
                                             Chebyshev coefficients.               */
                                             
/* Data structures */
typedef struct dpt_set_s_ *dpt_set;    /**< A set of precomputed data for a set 
                                            of DPT transforms of equal maximum 
                                            length.                                */

/**
 * Initializes a set of precomputed data for DPT transforms of equal length.
 * 
 * \arg M The maximum DPT transform index \f$M \in \mathbb{N}_0\f$. The individual
 *        transforms are addressed by and index number \f$m \in \mathbb{N}_0\f$ with 
 *        range \f$m = 0,\ldots,M\f$. The total number of transforms is therefore 
 *        \f$M+1\f$.
 * \arg t The exponent \f$t \in \mathbb{N}, t \ge 2\f$ of the transform length 
 *        \f$N = 2^t \in \mathbb{N}, N \ge 4\f$ 
 * \arg flags A bitwise combination of the flags DPT_NO_STABILIZATION, 
 *            DPT_BANDWIDTH_WINDOW
 *
 * \author Jens Keiner
 */
dpt_set dpt_init(const int M, const int t, const unsigned int flags);

/**
 * Computes the data required for a single DPT transform.
 * 
 * \arg set The set of DPT transform data where the computed data will be stored.
 * \arg m The transform index \f$m \in \mathbb{N}, 0 \le m \le M\f$.
 * \arg alpha The three-term recurrence coefficients \f$\alpha_k \in \mathbb{R}\f$ 
 *            for \f$k=0,\ldots,N\f$ such that \verbatim alpha[k] \endverbatim 
 *            \f$=\alpha_k\f$.
 * \arg beta The three-term recurrence coefficients \f$\beta_k \in \mathbb{R}\f$ 
 *            for \f$k=0,\ldots,N\f$ such that \verbatim beta[k] \endverbatim 
 *            \f$=\beta_k\f$.
 * \arg gamma The three-term recurrence coefficients \f$\gamma_k \in \mathbb{R}\f$ 
 *            for \f$k=0,\ldots,N\f$ such that \verbatim gamma[k] \endverbatim 
 *            \f$=\gamma_k\f$.
 * \arg k_start The index \f$k_{\text{start}} \in \mathbb{N}_0, 
 *              0 \le k_{\text{start}} \le N\f$
 * \arg threshold The treshold \f$\kappa \in \mathbb{R}, \kappa > 0\f$.
 *
 * \author Jens Keiner
 */
void dpt_precompute(dpt_set set, const int m, const double *alpha, 
                    const double *beta, const double *gamma, int k_start,
                    const double threshold);

/**
 * Computes a single DPT transform.
 * 
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void dpt_trafo(dpt_set set, const int m, const complex *x, complex *y, 
  const int k_end, const unsigned int flags);

/**
 * Computes a single DPT transform.
 * 
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void fpt_trafo(dpt_set set, const int m, const complex *x, complex *y, 
  const int k_end, const unsigned int flags);

/**
 * Computes a single DPT transform.
 * 
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void dpt_transposed(dpt_set set, const int m, complex *x, const complex *y, 
  const int k_end, const unsigned int flags);

/**
 * Computes a single DPT transform.
 * 
 * \arg set
 * \arg m
 * \arg x
 * \arg y
 * \arg k_end
 * \arg flags
 */
void fpt_transposed(dpt_set set, const int m, complex *x, const complex *y, 
  const int k_end, const unsigned int flags);

void dpt_finalize(dpt_set set); 

#endif /*DPT_H_*/

