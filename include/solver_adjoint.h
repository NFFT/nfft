/*
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* $Id$ */

/*! \file solver_adjoint.h
 *  \brief Header file for adjoint solver.
 */
#ifndef SOLVER_ADJOINT_H
#define SOLVER_ADJOINT_H

/** Include NFFT3 header. */
#include "nfft3.h"

/*
 * Macro for mangling an adjoint transform.
 * temporary added 01.2007 by tim becker
 */
                                                                               \
#define MACRO_SOLVER_ADJOINT_PLAN(MV, FLT, FLT_TYPE)                           \
                                                                               \
/** Structure for an adjoint transform plan  */                                \
typedef struct                                                                 \
{                                                                              \
  MV ## _plan *mv;                      /**< matrix vector multiplication    */\
  unsigned flags;                       /**< iteration type                  */\
                                                                               \
  double *w;                            /**< weighting factors               */\
  double *w_hat;                        /**< damping factors                 */\
                                                                               \
  FLT_TYPE *y_hat;                      /**< right hand side, samples        */\
                                                                               \
  FLT_TYPE *f_iter;                     /**< iterative solution              */\
  FLT_TYPE *r_hat_iter;                 /**< iterated residual vector        */\
  FLT_TYPE *z_iter;                     /**< residual of normal equation of    \
                                             first kind                      */\
  FLT_TYPE *p_iter;                     /**< search direction                */\
  FLT_TYPE *v_hat_iter;                 /**< residual vector update          */\
                                                                               \
  double alpha_iter;                    /**< step size for search direction  */\
  double beta_iter;                     /**< step size for search correction */\
                                                                               \
  double dot_r_hat_iter;                /**< weighted dotproduct of r_iter   */\
  double dot_r_hat_iter_old;            /**< previous dot_r_iter             */\
  double dot_z_iter;                    /**< weighted dotproduct of            \
                                             z_hat_iter                      */\
  double dot_z_iter_old;                /**< previous dot_z_hat_iter         */\
  double dot_p_iter;                    /**< weighted dotproduct of            \
                                             p_hat_iter                      */\
  double dot_v_hat_iter;                /**< weighted dotproduct of v_iter   */\
} i ## MV ## _adjoint_plan;                                                    \
                                                                               \
/** Simple initialisation. */                                                  \
void i ## MV ## _adjoint_init(adjoint ## MV ## _plan *ths, MV ## _plan *mv);   \
/** Advanced initialisation. */                                               \
void i ## MV ## _adjoint_init_advanced(adjoint ## MV ## _plan *ths, MV ## _plan,\
*mv, unsigned adjoint ## MV ## _flags);                                       \
/** Setting up residuals before the actual iteration.                       */\
void i ## MV ## _adjoint_before_loop(adjoint ## MV ## _plan *ths);              \
/** Doing one step in the iteration.                                        */\
void i ## MV ## _adjoint_loop_one_step(adjoint ## MV ## _plan *ths);            \
/** Destroys the plan for the adjoint transform.                            */\
void i ## MV ## _adjoint_finalize(adjoint ## MV ## _plan *ths);                 \

/** TODO: different solvers */
MACRO_SOLVER_ADJOINT_PLAN(nfsft, complex, double _Complex)
MACRO_SOLVER_ADJOINT_PLAN(nfft, complex, double _Complex)
MACRO_SOLVER_ADJOINT_PLAN(nfct, double, double)
MACRO_SOLVER_ADJOINT_PLAN(nfst, double, double)
MACRO_SOLVER_ADJOINT_PLAN(nnfft, complex, double _Complex)
MACRO_SOLVER_ADJOINT_PLAN(mri_inh_2d1d, complex, double _Complex)
MACRO_SOLVER_ADJOINT_PLAN(mri_inh_3d, complex, double _Complex)

#endif
/* solver_adjoint.h */
