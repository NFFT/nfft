#ifndef TRANSFORM_H
#define TRANSFORM_H

#include<nfft3.h>
#include<nfsft.h>

static inline void precompute(texture_plan *ths);

static nfsft_plan_old prepare_nfsft_plan(texture_plan *ths, int i);

static void process_results(texture_plan *ths, int i);

#endif
