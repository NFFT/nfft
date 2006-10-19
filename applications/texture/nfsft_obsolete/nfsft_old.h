/** 
 * \file nfsft.h
 * \brief Header file for the nfsft library
 * \author Jens Keiner
 */
#ifndef NFSFT_H
#define NFSFT_H

#include "nfft3.h"

/** \defgroup nfsft NFSFT */

/** 
 * If set, all computations are carried out with spherical harmonics normalized 
 * with respect to the \f$\text{L}^2\left(\mathbb{S}^2\right)\f$ standard inner 
 * product
 * \f[ 
 *   <f,g>_{\mathbb{S}^2} = \int_{0}^{\pi} \int_{-\pi}^{\pi} 
 *   f(\vartheta,\varphi) \overline{g(\vartheta,\varphi)} \; 
 *   d\varphi \; d\vartheta.
 * \f]
 * 
 * \see nfsft_init
 * \see nfsft_trafo_old
 * \see nfsft_adjoint_old
 * \see ndsft_trafo_old
 * \see ndsft_adjoint_old
 * \author Jens Keiner
 *  
 */
#define NFSFT_NORMALIZED_OLD 1<<0

/**
 * If set, only the fast transformations work.
 * 
 * \see nfsft_precompute_old
 * \see nfsft_trafo_old
 * \see nfsft_adjoint_old
 * \author Jens Keiner
 *  
 */
#define NFSFT_FAST_ONLY_OLD  1<<1

/**
 * If set, only fast transformations for bandwidth in a certain bandwidth
 * window will work. If \f$N\f$ is the power of two up to which precomputations 
 * is done, only transformations for bandwidth \f$M\f$ with \f$N/2 < M \le N\f$ 
 * will work.
 *
 * \see nfsft_precompute_old
 * \see nfsft_trafo_old
 * \see nfsft_adjoint_old
 * \author Jens Keiner
 *  
 */
#define NFSFT_BW_WINDOW_OLD  1<<2

/**
 * If set, the direct NDFT algorithm will be used.
 *
 * \see nfsft_trafo_old
 * \see nfsft_adjoint_old
 * \author Jens Keiner
 *  
 */
#define NFSFT_USE_NDFT_OLD 1<<3


/** 
 * Typedef for transform plans
 *
 * \author Jens Keiner
 *  
 */
typedef struct nfsft_plan_s_old *nfsft_plan_old;

/** 
 * Typedef for precomputation flags 
 * \author Jens Keiner
 *  
 */
typedef int nfsft_precompute_flags_old;

/** 
 * Typedef for transform flags 
 * \author Jens Keiner
 *  
 */
typedef int nfsft_flags_old;

/** 
 * Typedef for inverse transform plans 
 * \author Jens Keiner
 *  
 */
typedef struct infsft_plan_s_old *infsft_plan_old;

/** 
 * Typedef for inverse transform flags 
 * \author Jens Keiner
 *  
 */
typedef int infsft_flags_old;


/**
 * Creates a transform plan.
 *
 * \arg M The bandwidth \f$M\f$
 * \arg D The number of nodes \f$D\f$
 * \arg f_hat Fourier coefficients \f$\left(a_k^n\right)_{(k,n) \in 
 *   \mathcal{I}^M}\f$
 * \arg x The nodes \f$\left(\mathbf{\xi}_d\right)_{d = 0}^{D - 1}\f$
 * \arg f Function values \f$\left(f_d\right)_{d = 0}^{D - 1}\f$
 * \arg flags Flags
 *
 * \return The plan
 *
 * \author Jens Keiner
 *  
 */
nfsft_plan_old nfsft_init_old(int M, int D, complex **f_hat, double *x, complex *f, 
                      nfsft_flags_old flags);

nfsft_plan_old nfsft_init_guru_old(int M, int D, complex **f_hat, double *x, complex *f, 
                           nfsft_flags_old flags, int cutoff);

/**
 * Precomputes wisdom up to the next power of two with respect to a given 
 * bandwidth. Stabilization steps are precomputed for the given threshold.
 *
 * \arg M The bandwidth \F$M\f$
 * \arg threshold The threshold
 * \arg flags Flags
 *
 * \author Jens Keiner
 *  
 */
void nfsft_precompute_old(int M, double threshold, nfsft_precompute_flags_old flags);

/**
 * Forget all wisdom.
 *
 * \author Jens Keiner
 *  
 */
void nfsft_forget_old();

/**
 * Executes a direct NDSFT, i.e. computes for \f$d = 0,\ldots,D-1\f$
 * \f[
 *   f_d = f\left(\vartheta_d,\varphi_d\right) = 
 *         \sum_{(k,n) \in \mathcal{I}^M} a_k^n 
 *         Y_k^n\left(\vartheta_d,\varphi_d\right).  
 * \f]
 *
 * \arg plan The plan
 *
 * \author Jens Keiner
 *  
 */
void ndsft_trafo_old(nfsft_plan_old plan);

/**
 * Executes a direct adjoint NDSFT, i.e. computes for \f$(k,n) \in 
 * \mathcal{I}^M\f$
 * \f[
 *   a_k^n = \sum_{d = 0}^{D-1} f_d 
 *           Y_k^n\left(\vartheta_d,\varphi_d\right).  
 * \f]
 *
 * \arg plan The plan
 *
 * \author Jens Keiner
 *  
 */
void ndsft_adjoint_old(nfsft_plan_old plan);

/**
 * Executes a NFSFT, i.e. computes for \f$d = 0,\ldots,D-1\f$
 * \f[
 *   f_d = f\left(\vartheta_d,\varphi_d\right) = 
 *         \sum_{(k,n) \in \mathcal{I}^M} a_k^n 
 *         Y_k^n\left(\vartheta_d,\varphi_d\right).  
 * \f]
 *
 * \arg plan The plan
 *
 * \author Jens Keiner
 *  
 */
void nfsft_trafo_old(nfsft_plan_old plan);

/**
 * Executes an adjoint NFSFT, i.e. computes for \f$(k,n) \in 
 * \mathcal{I}^M\f$
 * \f[
 *   a_k^n = \sum_{d = 0}^{D-1} f_d 
 *           Y_k^n\left(\vartheta_d,\varphi_d\right).  
 * \f]
 *
 * \arg plan The plan
 *
 * \author Jens Keiner
 *  
 */
void nfsft_adjoint_old(nfsft_plan_old plan);

/**
 * Destroys a plan.
 *
 * \arg plan The plan
 *
 * \author Jens Keiner
 *  
 */
void nfsft_finalize_old(nfsft_plan_old plan);
#endif
