/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
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

/*! \file solver.c
 *  \brief Implementation file for the solver module
 *  \author Stefan Kunis
 */
#include "config.h"

#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include "nfft3.h"
#include "infft.h"

#undef X
#if defined(NFFT_SINGLE)
#define X(name) CONCAT(solverf_,name)
#elif defined(NFFT_LDOUBLE)
#define X(name) CONCAT(solverl_,name)
#else
#define X(name) CONCAT(solver_,name)
#endif

void X(init_advanced_complex)(X(plan_complex)* ths, Y(mv_plan_complex) *mv,
    unsigned flags)
{
  ths->mv = mv;
  ths->flags = flags;

  ths->y          = (C*)Y(malloc)((size_t)(ths->mv->M_total) * sizeof(C));
  ths->r_iter     = (C*)Y(malloc)((size_t)(ths->mv->M_total) * sizeof(C));
  ths->f_hat_iter = (C*)Y(malloc)((size_t)(ths->mv->N_total) * sizeof(C));
  ths->p_hat_iter = (C*)Y(malloc)((size_t)(ths->mv->N_total) * sizeof(C));

  if(ths->flags & LANDWEBER)
    ths->z_hat_iter = ths->p_hat_iter;

  if(ths->flags & STEEPEST_DESCENT)
    {
      ths->z_hat_iter = ths->p_hat_iter;
      ths->v_iter     = (C*)Y(malloc)((size_t)(ths->mv->M_total) * sizeof(C));
    }

  if(ths->flags & CGNR)
    {
      ths->z_hat_iter = (C*)Y(malloc)((size_t)(ths->mv->N_total) * sizeof(C));
      ths->v_iter     = (C*)Y(malloc)((size_t)(ths->mv->M_total) * sizeof(C));
    }

  if(ths->flags & CGNE)
    ths->z_hat_iter = ths->p_hat_iter;

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->w = (R*) Y(malloc)((size_t)(ths->mv->M_total) * sizeof(R));

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->w_hat = (R*) Y(malloc)((size_t)(ths->mv->N_total) * sizeof(R));
}

void X(init_complex)(X(plan_complex)* ths, Y(mv_plan_complex) *mv)
{
  X(init_advanced_complex)(ths, mv, CGNR);
}

void X(before_loop_complex)(X(plan_complex)* ths)
{
  Y(cp_complex)(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);

  CSWAP(ths->r_iter, ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  CSWAP(ths->r_iter, ths->mv->f);

  Y(upd_axpy_complex)(ths->r_iter, K(-1.0), ths->y, ths->mv->M_total);

  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))
    {
      if(ths->flags & PRECOMPUTE_WEIGHT)
	ths->dot_r_iter = Y(dot_w_complex)(ths->r_iter, ths->w,
					     ths->mv->M_total);
      else
	ths->dot_r_iter = Y(dot_complex)(ths->r_iter, ths->mv->M_total);
    }

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(cp_w_complex)(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    Y(cp_complex)(ths->mv->f, ths->r_iter, ths->mv->M_total);

  CSWAP(ths->z_hat_iter, ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  CSWAP(ths->z_hat_iter, ths->mv->f_hat);

  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))
    {
      if(ths->flags & PRECOMPUTE_DAMP)
	ths->dot_z_hat_iter = Y(dot_w_complex)(ths->z_hat_iter, ths->w_hat,
						 ths->mv->N_total);
      else
	ths->dot_z_hat_iter = Y(dot_complex)(ths->z_hat_iter,
					       ths->mv->N_total);
    }

  if(ths->flags & CGNE)
    ths->dot_p_hat_iter = ths->dot_z_hat_iter;

  if(ths->flags & CGNR)
    Y(cp_complex)(ths->p_hat_iter, ths->z_hat_iter, ths->mv->N_total);
} /* void solver_before_loop */

/** void solver_loop_one_step_landweber */
static void solver_loop_one_step_landweber_complex(X(plan_complex)* ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(upd_xpawy_complex)(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->z_hat_iter, ths->mv->N_total);
  else
    Y(upd_xpay_complex)(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  Y(cp_complex)(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);

  CSWAP(ths->r_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  CSWAP(ths->r_iter,ths->mv->f);

  Y(upd_axpy_complex)(ths->r_iter, K(-1.0), ths->y, ths->mv->M_total);

  if(ths->flags & NORMS_FOR_LANDWEBER)
    {
      if(ths->flags & PRECOMPUTE_WEIGHT)
	ths->dot_r_iter = Y(dot_w_complex)(ths->r_iter,ths->w,
					     ths->mv->M_total);
      else
	ths->dot_r_iter = Y(dot_complex)(ths->r_iter, ths->mv->M_total);
    }

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(cp_w_complex)(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    Y(cp_complex)(ths->mv->f, ths->r_iter, ths->mv->M_total);

  CSWAP(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  CSWAP(ths->z_hat_iter,ths->mv->f_hat);

  if(ths->flags & NORMS_FOR_LANDWEBER)
    {
      if(ths->flags & PRECOMPUTE_DAMP)
	ths->dot_z_hat_iter = Y(dot_w_complex)(ths->z_hat_iter, ths->w_hat,
						 ths->mv->N_total);
      else
	ths->dot_z_hat_iter = Y(dot_complex)(ths->z_hat_iter,
					       ths->mv->N_total);
    }
} /* void solver_loop_one_step_landweber */

/** void solver_loop_one_step_steepest_descent */
static void solver_loop_one_step_steepest_descent_complex(X(plan_complex) *ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(cp_w_complex)(ths->mv->f_hat, ths->w_hat, ths->z_hat_iter,
		      ths->mv->N_total);
  else
    Y(cp_complex)(ths->mv->f_hat, ths->z_hat_iter, ths->mv->N_total);

  CSWAP(ths->v_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  CSWAP(ths->v_iter,ths->mv->f);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_v_iter = Y(dot_w_complex)(ths->v_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_v_iter = Y(dot_complex)(ths->v_iter, ths->mv->M_total);

  /*-----------------*/
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(upd_xpawy_complex)(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->z_hat_iter, ths->mv->N_total);
  else
    Y(upd_xpay_complex)(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  Y(upd_xpay_complex)(ths->r_iter, -ths->alpha_iter, ths->v_iter,
			ths->mv->M_total);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = Y(dot_w_complex)(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = Y(dot_complex)(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(cp_w_complex)(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    Y(cp_complex)(ths->mv->f, ths->r_iter, ths->mv->M_total);

  CSWAP(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  CSWAP(ths->z_hat_iter,ths->mv->f_hat);

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_z_hat_iter = Y(dot_w_complex)(ths->z_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_z_hat_iter = Y(dot_complex)(ths->z_hat_iter, ths->mv->N_total);
} /* void solver_loop_one_step_steepest_descent */

/** void solver_loop_one_step_cgnr */
static void solver_loop_one_step_cgnr_complex(X(plan_complex) *ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(cp_w_complex)(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,
		      ths->mv->N_total);
  else
    Y(cp_complex)(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);

  CSWAP(ths->v_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  CSWAP(ths->v_iter,ths->mv->f);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_v_iter = Y(dot_w_complex)(ths->v_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_v_iter = Y(dot_complex)(ths->v_iter, ths->mv->M_total);

  /*-----------------*/
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(upd_xpawy_complex)(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->p_hat_iter, ths->mv->N_total);
  else
    Y(upd_xpay_complex)(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  Y(upd_xpay_complex)(ths->r_iter, -ths->alpha_iter, ths->v_iter,
			ths->mv->M_total);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = Y(dot_w_complex)(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = Y(dot_complex)(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(cp_w_complex)(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    Y(cp_complex)(ths->mv->f, ths->r_iter, ths->mv->M_total);

  CSWAP(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  CSWAP(ths->z_hat_iter,ths->mv->f_hat);

  ths->dot_z_hat_iter_old = ths->dot_z_hat_iter;
  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_z_hat_iter = Y(dot_w_complex)(ths->z_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_z_hat_iter = Y(dot_complex)(ths->z_hat_iter, ths->mv->N_total);

  /*-----------------*/
  ths->beta_iter = ths->dot_z_hat_iter / ths->dot_z_hat_iter_old;

  /*-----------------*/
  Y(upd_axpy_complex)(ths->p_hat_iter, ths->beta_iter, ths->z_hat_iter,
			ths->mv->N_total);
} /* void solver_loop_one_step_cgnr */

/** void solver_loop_one_step_cgne */
static void solver_loop_one_step_cgne_complex(X(plan_complex) *ths)
{
  ths->alpha_iter = ths->dot_r_iter / ths->dot_p_hat_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(upd_xpawy_complex)(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->p_hat_iter, ths->mv->N_total);
  else
    Y(upd_xpay_complex)(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,
                          ths->mv->N_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(cp_w_complex)(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,
		      ths->mv->N_total);
  else
    Y(cp_complex)(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);

  ths->mv->mv_trafo(ths->mv);

  Y(upd_xpay_complex)(ths->r_iter, -ths->alpha_iter, ths->mv->f,
			ths->mv->M_total);

  ths->dot_r_iter_old = ths->dot_r_iter;
  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = Y(dot_w_complex)(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = Y(dot_complex)(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  ths->beta_iter = ths->dot_r_iter / ths->dot_r_iter_old;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(cp_w_complex)(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    Y(cp_complex)(ths->mv->f, ths->r_iter, ths->mv->M_total);

  ths->mv->mv_adjoint(ths->mv);

  Y(upd_axpy_complex)(ths->p_hat_iter, ths->beta_iter, ths->mv->f_hat,
			ths->mv->N_total);

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_p_hat_iter = Y(dot_w_complex)(ths->p_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_p_hat_iter = Y(dot_complex)(ths->p_hat_iter, ths->mv->N_total);
} /* void solver_loop_one_step_cgne */

/** void solver_loop_one_step */
void X(loop_one_step_complex)(X(plan_complex) *ths)
{
  if(ths->flags & LANDWEBER)
    solver_loop_one_step_landweber_complex(ths);

  if(ths->flags & STEEPEST_DESCENT)
    solver_loop_one_step_steepest_descent_complex(ths);

  if(ths->flags & CGNR)
    solver_loop_one_step_cgnr_complex(ths);

  if(ths->flags & CGNE)
    solver_loop_one_step_cgne_complex(ths);
} /* void solver_loop_one_step */

/** void solver_finalize */
void X(finalize_complex)(X(plan_complex) *ths)
{
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(free)(ths->w);

  if(ths->flags & PRECOMPUTE_DAMP)
    Y(free)(ths->w_hat);

  if(ths->flags & CGNR)
    {
      Y(free)(ths->v_iter);
      Y(free)(ths->z_hat_iter);
    }

  if(ths->flags & STEEPEST_DESCENT)
    Y(free)(ths->v_iter);

  Y(free)(ths->p_hat_iter);
  Y(free)(ths->f_hat_iter);

  Y(free)(ths->r_iter);
  Y(free)(ths->y);
} /** void solver_finalize */


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void X(init_advanced_double)(X(plan_double)* ths, Y(mv_plan_double) *mv, unsigned flags)
{
  ths->mv = mv;
  ths->flags = flags;

  ths->y          = (R*)Y(malloc)((size_t)(ths->mv->M_total) * sizeof(R));
  ths->r_iter     = (R*)Y(malloc)((size_t)(ths->mv->M_total) * sizeof(R));
  ths->f_hat_iter = (R*)Y(malloc)((size_t)(ths->mv->N_total) * sizeof(R));
  ths->p_hat_iter = (R*)Y(malloc)((size_t)(ths->mv->N_total) * sizeof(R));

  if(ths->flags & LANDWEBER)
    ths->z_hat_iter = ths->p_hat_iter;

  if(ths->flags & STEEPEST_DESCENT)
    {
      ths->z_hat_iter = ths->p_hat_iter;
      ths->v_iter     = (R*)Y(malloc)((size_t)(ths->mv->M_total) * sizeof(R));
    }

  if(ths->flags & CGNR)
    {
      ths->z_hat_iter = (R*)Y(malloc)((size_t)(ths->mv->N_total) * sizeof(R));
      ths->v_iter     = (R*)Y(malloc)((size_t)(ths->mv->M_total) * sizeof(R));
    }

  if(ths->flags & CGNE)
    ths->z_hat_iter = ths->p_hat_iter;

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->w = (R*) Y(malloc)((size_t)(ths->mv->M_total) * sizeof(R));

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->w_hat = (R*) Y(malloc)((size_t)(ths->mv->N_total) * sizeof(R));
}

void X(init_double)(X(plan_double)* ths, Y(mv_plan_double) *mv)
{
  X(init_advanced_double)(ths, mv, CGNR);
}

void X(before_loop_double)(X(plan_double)* ths)
{
  Y(cp_double)(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);

  RSWAP(ths->r_iter, ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  RSWAP(ths->r_iter, ths->mv->f);

  Y(upd_axpy_double)(ths->r_iter, K(-1.0), ths->y, ths->mv->M_total);

  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))
    {
      if(ths->flags & PRECOMPUTE_WEIGHT)
	ths->dot_r_iter = Y(dot_w_double)(ths->r_iter, ths->w,
					     ths->mv->M_total);
      else
	ths->dot_r_iter = Y(dot_double)(ths->r_iter, ths->mv->M_total);
    }

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(cp_w_double)(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    Y(cp_double)(ths->mv->f, ths->r_iter, ths->mv->M_total);

  RSWAP(ths->z_hat_iter, ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  RSWAP(ths->z_hat_iter, ths->mv->f_hat);

  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))
    {
      if(ths->flags & PRECOMPUTE_DAMP)
	ths->dot_z_hat_iter = Y(dot_w_double)(ths->z_hat_iter, ths->w_hat,
						 ths->mv->N_total);
      else
	ths->dot_z_hat_iter = Y(dot_double)(ths->z_hat_iter,
					       ths->mv->N_total);
    }

  if(ths->flags & CGNE)
    ths->dot_p_hat_iter = ths->dot_z_hat_iter;

  if(ths->flags & CGNR)
    Y(cp_double)(ths->p_hat_iter, ths->z_hat_iter, ths->mv->N_total);
} /* void solver_before_loop */

/** void solver_loop_one_step_landweber */
static void solver_loop_one_step_landweber_double(X(plan_double)* ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(upd_xpawy_double)(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->z_hat_iter, ths->mv->N_total);
  else
    Y(upd_xpay_double)(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  Y(cp_double)(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);

  RSWAP(ths->r_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  RSWAP(ths->r_iter,ths->mv->f);

  Y(upd_axpy_double)(ths->r_iter, K(-1.0), ths->y, ths->mv->M_total);

  if(ths->flags & NORMS_FOR_LANDWEBER)
    {
      if(ths->flags & PRECOMPUTE_WEIGHT)
	ths->dot_r_iter = Y(dot_w_double)(ths->r_iter,ths->w,
					     ths->mv->M_total);
      else
	ths->dot_r_iter = Y(dot_double)(ths->r_iter, ths->mv->M_total);
    }

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(cp_w_double)(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    Y(cp_double)(ths->mv->f, ths->r_iter, ths->mv->M_total);

  RSWAP(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  RSWAP(ths->z_hat_iter,ths->mv->f_hat);

  if(ths->flags & NORMS_FOR_LANDWEBER)
    {
      if(ths->flags & PRECOMPUTE_DAMP)
	ths->dot_z_hat_iter = Y(dot_w_double)(ths->z_hat_iter, ths->w_hat,
						 ths->mv->N_total);
      else
	ths->dot_z_hat_iter = Y(dot_double)(ths->z_hat_iter,
					       ths->mv->N_total);
    }
} /* void solver_loop_one_step_landweber */

/** void solver_loop_one_step_steepest_descent */
static void solver_loop_one_step_steepest_descent_double(X(plan_double) *ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(cp_w_double)(ths->mv->f_hat, ths->w_hat, ths->z_hat_iter,
		      ths->mv->N_total);
  else
    Y(cp_double)(ths->mv->f_hat, ths->z_hat_iter, ths->mv->N_total);

  RSWAP(ths->v_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  RSWAP(ths->v_iter,ths->mv->f);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_v_iter = Y(dot_w_double)(ths->v_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_v_iter = Y(dot_double)(ths->v_iter, ths->mv->M_total);

  /*-----------------*/
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(upd_xpawy_double)(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->z_hat_iter, ths->mv->N_total);
  else
    Y(upd_xpay_double)(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  Y(upd_xpay_double)(ths->r_iter, -ths->alpha_iter, ths->v_iter,
			ths->mv->M_total);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = Y(dot_w_double)(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = Y(dot_double)(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(cp_w_double)(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    Y(cp_double)(ths->mv->f, ths->r_iter, ths->mv->M_total);

  RSWAP(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  RSWAP(ths->z_hat_iter,ths->mv->f_hat);

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_z_hat_iter = Y(dot_w_double)(ths->z_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_z_hat_iter = Y(dot_double)(ths->z_hat_iter, ths->mv->N_total);
} /* void solver_loop_one_step_steepest_descent */

/** void solver_loop_one_step_cgnr */
static void solver_loop_one_step_cgnr_double(X(plan_double) *ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(cp_w_double)(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,
		      ths->mv->N_total);
  else
    Y(cp_double)(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);

  RSWAP(ths->v_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  RSWAP(ths->v_iter,ths->mv->f);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_v_iter = Y(dot_w_double)(ths->v_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_v_iter = Y(dot_double)(ths->v_iter, ths->mv->M_total);

  /*-----------------*/
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(upd_xpawy_double)(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->p_hat_iter, ths->mv->N_total);
  else
    Y(upd_xpay_double)(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  Y(upd_xpay_double)(ths->r_iter, -ths->alpha_iter, ths->v_iter,
			ths->mv->M_total);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = Y(dot_w_double)(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = Y(dot_double)(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(cp_w_double)(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    Y(cp_double)(ths->mv->f, ths->r_iter, ths->mv->M_total);

  RSWAP(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  RSWAP(ths->z_hat_iter,ths->mv->f_hat);

  ths->dot_z_hat_iter_old = ths->dot_z_hat_iter;
  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_z_hat_iter = Y(dot_w_double)(ths->z_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_z_hat_iter = Y(dot_double)(ths->z_hat_iter, ths->mv->N_total);

  /*-----------------*/
  ths->beta_iter = ths->dot_z_hat_iter / ths->dot_z_hat_iter_old;

  /*-----------------*/
  Y(upd_axpy_double)(ths->p_hat_iter, ths->beta_iter, ths->z_hat_iter,
			ths->mv->N_total);
} /* void solver_loop_one_step_cgnr */

/** void solver_loop_one_step_cgne */
static void solver_loop_one_step_cgne_double(X(plan_double) *ths)
{
  ths->alpha_iter = ths->dot_r_iter / ths->dot_p_hat_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(upd_xpawy_double)(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->p_hat_iter, ths->mv->N_total);
  else
    Y(upd_xpay_double)(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,
                          ths->mv->N_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    Y(cp_w_double)(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,
		      ths->mv->N_total);
  else
    Y(cp_double)(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);

  ths->mv->mv_trafo(ths->mv);

  Y(upd_xpay_double)(ths->r_iter, -ths->alpha_iter, ths->mv->f,
			ths->mv->M_total);

  ths->dot_r_iter_old = ths->dot_r_iter;
  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = Y(dot_w_double)(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = Y(dot_double)(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  ths->beta_iter = ths->dot_r_iter / ths->dot_r_iter_old;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(cp_w_double)(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    Y(cp_double)(ths->mv->f, ths->r_iter, ths->mv->M_total);

  ths->mv->mv_adjoint(ths->mv);

  Y(upd_axpy_double)(ths->p_hat_iter, ths->beta_iter, ths->mv->f_hat,
			ths->mv->N_total);

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_p_hat_iter = Y(dot_w_double)(ths->p_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_p_hat_iter = Y(dot_double)(ths->p_hat_iter, ths->mv->N_total);
} /* void solver_loop_one_step_cgne */

/** void solver_loop_one_step */
void X(loop_one_step_double)(X(plan_double) *ths)
{
  if(ths->flags & LANDWEBER)
    solver_loop_one_step_landweber_double(ths);

  if(ths->flags & STEEPEST_DESCENT)
    solver_loop_one_step_steepest_descent_double(ths);

  if(ths->flags & CGNR)
    solver_loop_one_step_cgnr_double(ths);

  if(ths->flags & CGNE)
    solver_loop_one_step_cgne_double(ths);
} /* void solver_loop_one_step */

/** void solver_finalize */
void X(finalize_double)(X(plan_double) *ths)
{
  if(ths->flags & PRECOMPUTE_WEIGHT)
    Y(free)(ths->w);

  if(ths->flags & PRECOMPUTE_DAMP)
    Y(free)(ths->w_hat);

  if(ths->flags & CGNR)
    {
      Y(free)(ths->v_iter);
      Y(free)(ths->z_hat_iter);
    }

  if(ths->flags & STEEPEST_DESCENT)
    Y(free)(ths->v_iter);

  Y(free)(ths->p_hat_iter);
  Y(free)(ths->f_hat_iter);

  Y(free)(ths->r_iter);
  Y(free)(ths->y);
} /** void solver_finalize */
