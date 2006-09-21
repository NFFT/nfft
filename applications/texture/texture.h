#ifndef TRANSFORM_H
#define TRANSFORM_H

#include<assert.h>
#include"nfft3_texture.h"

/**
 * @addtogroup texture_private
 * @{
 */

/** @file texture.h
 * Containts declarations and datatype definitions used by the implementation
 * of the texture transform.
 * They are not part of the public interface.
 * texture.h may only be included by texture.c.
 */

/** Datatype for a variable storing the current state of precomputation.
 * @see texture_precompute
 * @see texture_precompute_advanced
 */
typedef struct precompute_params_type_ {

	/** Stores the precomputed bandwidth.
	 */
	int N;

	/** Stores the last value of texture_precompute_flags given to 
	 * ::texture_precomputed_advanced.
	 */
	unsigned int texture_precompute_flags;

	/** Stores the last value of nfsft_precompute_flags given to 
	 * ::texture_precompute_advanced.
	 */
	unsigned int nfsft_precompute_flags;

	/** Stores the last value of nfsft_threshold given to 
	 * ::texture_precompute_advanced.
	 */
	double nfsft_threshold;
} precompute_params_type;

/** Stores the state of precomputation.
 * @see texture_precompute
 * @see texture_precompute_advanced
 */
static precompute_params_type precompute_params;

/** Stores if the precomputed data is available, i.e. if the transformation 
 * can be carried out.
 * @see texture_precompute
 * @see texture_precompute_advanced
 * @see texture_forget
 */
static int is_precomputed = 0;

/** Helper method for texture_trafo. 
 * 
 * Calling nfsft_trafo on the nfsft_plan nfsft calculates the sample values
 * for some fixed pole figure and puts the result in ths. The values in ths->x
 * belonging to different pole figures and the remaining fields of ths remain
 * unchanged. 
 * 
 * @par ths gives the input for the nfsft_plan returned and provides the
 * storage for the result of the nfsft_trafo.
 * @par nfsft will be prepared.
 * @par i specifies the pole figure.
 */
static void prepare_nfsft_plan(texture_plan * ths, nfsft_plan * nfsft, int i);

/** Helper method for texture_adjoint.
 * Adds the result of the adjoined transformation restricted to a single pole
 * figure to the frequencies in a plan.
 *
 * @par ths is the plan giving the sample values and providing the storage for
 * the frequencies.
 * @par nfsft is the nfsft plan used to calculate ths->nfsft_f_hat.
 * @par i specifies the pole figure.
 *
 * @pre ths->nfsft_f_hat has to contain the result of nfsft_adjoint applied to
 * the 
 * samples of the i-th pole figure.
 */
static void process_results(texture_plan * ths, nfsft_plan * nfsft, int i);

/**@}
 */
#endif
