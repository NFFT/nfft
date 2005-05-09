/**
 * Header file with the public interface of the 
 * Nonequispaced Fast Fourier Transform (NFSFT) library
 */
#ifndef NFSFT_H
#define NFSFT_H

/* Include the header file for the complex datatype. */
#include <complex.h>

/** \defgroup nfsft NFSFT: public API */

/** 
 * Typedef for transform plans 
 * \ingroup nfsft 
 */
typedef struct nfsft_plan_s *nfsft_plan;

/** 
 * Typedef for transform flags 
 * \ingroup nfsft 
 */
typedef int nfsft_flags;

/**
 * Creates a transform plan.
 *
 * \arg type The type of transform, i.e. NFSFT_ITERATIVE, NFSFT_FORWARD or 
 *   NFSFT_ADJOINT
 * \arg D The number of nodes
 * \arg M The bandwidth
 * \arg phi Angles phi of nodes
 * \arg theta Angles theta of nodes
 * \arg f_hat Fourier coefficients
 * \arg f Function values
 * \arg flags Flags
 *
 * \return The plan
 *
 * \ingroup nfsft 
 */
nfsft_plan nfsft_init(int D, int M, double *angles, complex **f_hat, complex *f, 
           nfsft_flags flags);

/**
 * Executes a NFSFT.
 *
 * \arg plan The plan
 *
 * \ingroup nfsft 
 */
void nfsft_trafo(nfsft_plan plan);

/**
 * Executes a direct NDSFT
 *
 * \arg plan The plan
 *
 * \ingroup nfsft 
 */
void ndsft_trafo(nfsft_plan plan);

/**
 * Executes a adjoint NFSFT
 *
 * \arg plan The plan
 *
 * \ingroup nfsft 
 */
void nfsft_adjoint(nfsft_plan plan);

/**
 * Executes a direct adjoint NDSFT
 *
 * \arg plan The plan
 *
 * \ingroup nfsft 
 */
void ndsft_adjoint(nfsft_plan plan);

/**
 * Destroys a plan.
 *
 * \arg plan The plan
 *
 * \ingroup nfsft 
 *
 * \ingroup nfsft 
 */
void nfsft_finalize(nfsft_plan plan);

/**
 * Precomputes wisdom for a given bandwidth. Wisdom is computed in 
 * steps of powers of 2.
 *
 * \arg m The bandwidth
 *
 * \ingroup nfsft 
 */
void nfsft_compute_wisdom(int m, int threshold);

/**
 * Forgets all wisdom computed.
 *
 * \ingroup nfsft 
 */
void nfsft_forget_wisdom();

/**
 * Exports wisdom to a file
 *
 * \arg filename The filename
 *
 * \ingroup nfsft 
 */
/*void nfsft_export_wisdom(const char* filename);*/
#endif
