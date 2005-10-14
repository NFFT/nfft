#ifndef TRANSFORM_H
#define TRANSFORM_H

#include<nfft3.h>
#include<nfsft.h>

/**
 * \addtogroup texture_private
 * @{
 */

typedef struct precompute_params_type_ {
	int N;
	unsigned int texture_precompute_flags; 
	unsigned int nfsft_precompute_flags;
	double nfsft_threshold;
} precompute_params_type;

static precompute_params_type precompute_params;

static int is_precomputed = 0;

static nfsft_plan_old prepare_nfsft_plan(texture_plan *ths, int i);

static void process_results(texture_plan *ths, int i);

/**@}
 */
#endif
