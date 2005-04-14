#ifndef NFSFT_H
#define NFSFT_H

#include <complex.h>

/** 
* Constant symbols for transform types.
*/
#define NFSFT_FORWARD    (1U<< 0)
#define NFSFT_BACKWARD   (1U<< 1)
#define NFSFT_ADJOINT    (1U<< 2)

/** Transform plan */
typedef struct nfsft_plan_s *nfsft_plan;

typedef int nfsft_flags;

/**
 * Creates a transform plan-
 *
 * \arg type The type of transform, i.e. NFSFT_FORWARD, NFSFT_BACKWARD or 
 *   NFSFT_ADJOINT
 * \arg D The number of nodes
 * \arg M The bandwidth
 * \arg phi Angles phi of nodes
 * \arg theta Angles theta of nodes
 * \arg f_hat Fourier coefficients
 * \arg f Function values
 * \arg flags Flags
 */
nfsft_plan nfsft_create_plan(int type, int d, int m, double *angles, 
  complex **f_hat, complex *f, nfsft_flags flags);

/**
 * Destroys a plan.
 *
 * \arg plan THe plan
 */
void nfsft_destroy_plan(nfsft_plan plan);

/**
 * Executes a plan using the fast algorithm.
 *
 * \arg plan The plan
 */
void nfsft_execute(nfsft_plan plan);

/**
 * Executes a plan using the slow algorithm.
 *
 * \arg plan The plan
 */
void ndsft_execute(nfsft_plan plan);

/**
 * Precomputes wisdom for a given bandwidth. In general, wisdom is computed in 
 * steps of powers of 2.
 *
 * \arg m The bandwidth
 */
void nfsft_compute_wisdom(int m);

/**
 * Forgets all wisdom computed.
 */
void nfsft_forget_wisdom();

/**
 * Exports wisdom to a file
 *
 * \arg filename The filename
 */
//void nfsft_export_wisdom(const char* filename);

void nfsft_test();
#endif
