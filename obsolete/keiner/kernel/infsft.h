#ifndef INFSFT_H
#define INFSFT_H

#include "nfsft.h"

#define NFSFT_CGNR_E            (1U<<0)
#define NFSFT_CGNE_R            (1U<<1)
#define NFSFT_LANDWEBER         (1U<<2)
#define NFSFT_STEEPEST_DESCENT  (1U<<3)
#define NFSFT_ITERATE_2nd       (1U<<4)
#define NFSFT_PRECOMPUTE_WEIGHT (1U<<5)
#define NFSFT_PRECOMPUTE_DAMP   (1U<<6)
#define NFSFT_NORMS_FOR_LANDWEBER    (1U<<7)

void infsft_init(infsft_plan this_iplan, nfsft_plan direct_plan);

void infsft_init_guru(infsft_plan this_iplan, nfsft_plan direct_plan, 
                      int infsft_flags);

infsft_plan infsft_make_plan();

void infsft_before_loop(infsft_plan this_iplan);

void infsft_loop_one_step(infsft_plan this_iplan);

void infsft_finalize(infsft_plan this_iplan);

#endif