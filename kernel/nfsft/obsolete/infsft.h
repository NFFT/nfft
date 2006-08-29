#ifndef INFSFT_H
#define INFSFT_H

#include "nfsft_old.h"

#define NFSFT_CGNR_E            (1U<<0)
#define NFSFT_CGNE_R            (1U<<1)
#define NFSFT_LANDWEBER         (1U<<2)
#define NFSFT_STEEPEST_DESCENT  (1U<<3)
#define NFSFT_ITERATE_2nd       (1U<<4)
#define NFSFT_PRECOMPUTE_WEIGHT (1U<<5)
#define NFSFT_PRECOMPUTE_DAMP   (1U<<6)
#define NFSFT_NORMS_FOR_LANDWEBER    (1U<<7)

void infsft_init_old(infsft_plan_old this_iplan, nfsft_plan_old direct_plan);

void infsft_init_guru_old(infsft_plan_old this_iplan, nfsft_plan_old direct_plan, 
                      int infsft_flags_old);

infsft_plan_old infsft_make_plan_old();

void infsft_before_loop_old(infsft_plan_old this_iplan);

void infsft_loop_one_step_old(infsft_plan_old this_iplan);

void infsft_finalize_old(infsft_plan_old this_iplan);

#endif
