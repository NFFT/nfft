#ifndef INFSFT_H
#define INFSFT_H

#include "nfsft.h"

#define CGNR_E            1<<0
#define CGNE_R            1<<1
#define LANDWEBER         1<<2
#define STEEPEST_DESCENT  1<<3
#define ITERATE_2nd       1<<4
#define PRECOMPUTE_WEIGHT 1<<5
#define PRECOMPUTE_DAMP   1<<6

void infsft_init(infsft_plan this_iplan, nfsft_plan direct_plan);

void infsft_init_guru(infsft_plan this_iplan, nfsft_plan direct_plan, 
                      int infsft_flags);

infsft_plan infsft_make_plan();

void infsft_before_loop(infsft_plan this_iplan);

void infsft_loop_one_step(infsft_plan this_iplan);

void infsft_finalize(infsft_plan this_iplan);

#endif