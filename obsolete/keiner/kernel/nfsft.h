/** 
 * \file nfsft.h
 * \brief Header file for the nfsft library
 * \author Jens Keiner
 */
#ifndef NFSFT_H
#define NFSFT_H

#include <complex.h>
#include <fftw3.h>

/** \defgroup nfsft NFSFT */

/** \defgroup nfsft_public_api Public API
 *  \ingroup nfsft
 */

/** Flag for normalization */
#define NFSFT_NORMALIZED 1<<0
#define NFSFT_FAST_ONLY  1<<1
#define NFSFT_BW_WINDOW  1<<2

/** 
 * Typedef for transform plans 
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
typedef struct nfsft_plan_s *nfsft_plan;

/** 
 * Typedef for transform flags 
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
typedef int nfsft_flags;

/** 
 * Typedef for inverse transform plans 
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
typedef struct infsft_plan_s *infsft_plan;

/** 
 * Typedef for inverse transform flags 
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
typedef int infsft_flags;


/**
 * Creates a transform plan.
 *
 * \arg D The number of nodes
 * \arg M The bandwidth \f$M\f$
 * \arg x The nodes
 * \arg f_hat Fourier coefficients
 * \arg f Function values
 * \arg flags NFSFT flags
 *
 * \return The plan
 *
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
nfsft_plan nfsft_init(int M, int N, double *x, complex **f_hat, complex *f, 
                      nfsft_flags flags);

/**
 * Precomputes wisdom up to the next power of two relative to a  given bandwidth.
 *
 * \arg m The bandwidth
 * \arg threshold The threshold
 *
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
void nfsft_precompute(int m, double threshold, int flags);

/**
 * Forgets all wisdom computed.
 *
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
void nfsft_forget();

/**
 * Executes a direct NDSFT
 *
 * \arg plan The plan
 *
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
void ndsft_trafo(nfsft_plan plan);

/**
 * Executes a direct adjoint NDSFT
 *
 * \arg plan The plan
 *
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
void ndsft_adjoint(nfsft_plan plan);

/**
 * Executes a NFSFT.
 *
 * \arg plan The plan
 *
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
void nfsft_trafo(nfsft_plan plan);

/**
 * Executes an adjoint NFSFT.
 *
 * \arg plan The plan
 *
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
void nfsft_adjoint(nfsft_plan plan);

/**
 * Destroys a plan.
 *
 * \arg plan The plan
 *
 * \author Jens Keiner
 * \ingroup nfsft_public_api 
 */
void nfsft_finalize(nfsft_plan plan);
#endif
