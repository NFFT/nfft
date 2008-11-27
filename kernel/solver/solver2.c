/* $Id: nfft3.h 2168 2008-03-13 14:55:11Z kunis $
 *
 * Copyright (c) 2008  Jens Keiner, Stefan Kunis, Daniel Potts
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

#include <complex.h>

#include "util.h"
#include "nfft3.h"

void solver_init_advanced_complex(solver_plan_complex* ths, mv_plan_complex *mv, unsigned flags)
{
  ths->mv = mv;
  ths->flags = flags;

  ths->y          = (fftw_complex*)nfft_malloc(ths->mv->M_total*sizeof(fftw_complex));
  ths->r_iter     = (fftw_complex*)nfft_malloc(ths->mv->M_total*sizeof(fftw_complex));
  ths->f_hat_iter = (fftw_complex*)nfft_malloc(ths->mv->N_total*sizeof(fftw_complex));
  ths->p_hat_iter = (fftw_complex*)nfft_malloc(ths->mv->N_total*sizeof(fftw_complex));

  if(ths->flags & LANDWEBER)
    ths->z_hat_iter = ths->p_hat_iter;

  if(ths->flags & STEEPEST_DESCENT)
    {
      ths->z_hat_iter = ths->p_hat_iter;
      ths->v_iter     = (fftw_complex*)nfft_malloc(ths->mv->M_total*sizeof(fftw_complex));
    }

  if(ths->flags & CGNR)
    {
      ths->z_hat_iter = (fftw_complex*)nfft_malloc(ths->mv->N_total*sizeof(fftw_complex));
      ths->v_iter     = (fftw_complex*)nfft_malloc(ths->mv->M_total*sizeof(fftw_complex));
    }

  if(ths->flags & CGNE)
    ths->z_hat_iter = ths->p_hat_iter;

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->w = (double*) nfft_malloc(ths->mv->M_total*sizeof(double));

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->w_hat = (double*) nfft_malloc(ths->mv->N_total*sizeof(double));
}

void solver_init_complex(solver_plan_complex* ths, mv_plan_complex *mv)
{
  solver_init_advanced_complex(ths, mv, CGNR);
}

void solver_before_loop_complex(solver_plan_complex* ths)
{
  nfft_cp_complex(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);

  NFFT_SWAP_complex(ths->r_iter, ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  NFFT_SWAP_complex(ths->r_iter, ths->mv->f);

  nfft_upd_axpy_complex(ths->r_iter, -1.0, ths->y, ths->mv->M_total);

  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))
    {
      if(ths->flags & PRECOMPUTE_WEIGHT)
	ths->dot_r_iter = nfft_dot_w_complex(ths->r_iter, ths->w,
					     ths->mv->M_total);
      else
	ths->dot_r_iter = nfft_dot_complex(ths->r_iter, ths->mv->M_total);
    }

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_cp_w_complex(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    nfft_cp_complex(ths->mv->f, ths->r_iter, ths->mv->M_total);

  NFFT_SWAP_complex(ths->z_hat_iter, ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  NFFT_SWAP_complex(ths->z_hat_iter, ths->mv->f_hat);

  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))
    {
      if(ths->flags & PRECOMPUTE_DAMP)
	ths->dot_z_hat_iter = nfft_dot_w_complex(ths->z_hat_iter, ths->w_hat,
						 ths->mv->N_total);
      else
	ths->dot_z_hat_iter = nfft_dot_complex(ths->z_hat_iter,
					       ths->mv->N_total);
    }

  if(ths->flags & CGNE)
    ths->dot_p_hat_iter = ths->dot_z_hat_iter;

  if(ths->flags & CGNR)
    nfft_cp_complex(ths->p_hat_iter, ths->z_hat_iter, ths->mv->N_total);
} /* void solver_before_loop */

/** void solver_loop_one_step_landweber */
void solver_loop_one_step_landweber_complex(solver_plan_complex* ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_upd_xpawy_complex(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->z_hat_iter, ths->mv->N_total);
  else
    nfft_upd_xpay_complex(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  nfft_cp_complex(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);

  NFFT_SWAP_complex(ths->r_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  NFFT_SWAP_complex(ths->r_iter,ths->mv->f);

  nfft_upd_axpy_complex(ths->r_iter, -1.0, ths->y, ths->mv->M_total);

  if(ths->flags & NORMS_FOR_LANDWEBER)
    {
      if(ths->flags & PRECOMPUTE_WEIGHT)
	ths->dot_r_iter = nfft_dot_w_complex(ths->r_iter,ths->w,
					     ths->mv->M_total);
      else
	ths->dot_r_iter = nfft_dot_complex(ths->r_iter, ths->mv->M_total);
    }

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_cp_w_complex(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    nfft_cp_complex(ths->mv->f, ths->r_iter, ths->mv->M_total);

  NFFT_SWAP_complex(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  NFFT_SWAP_complex(ths->z_hat_iter,ths->mv->f_hat);

  if(ths->flags & NORMS_FOR_LANDWEBER)
    {
      if(ths->flags & PRECOMPUTE_DAMP)
	ths->dot_z_hat_iter = nfft_dot_w_complex(ths->z_hat_iter, ths->w_hat,
						 ths->mv->N_total);
      else
	ths->dot_z_hat_iter = nfft_dot_complex(ths->z_hat_iter,
					       ths->mv->N_total);
    }
} /* void solver_loop_one_step_landweber */

/** void solver_loop_one_step_steepest_descent */
void solver_loop_one_step_steepest_descent_complex(solver_plan_complex *ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_cp_w_complex(ths->mv->f_hat, ths->w_hat, ths->z_hat_iter,
		      ths->mv->N_total);
  else
    nfft_cp_complex(ths->mv->f_hat, ths->z_hat_iter, ths->mv->N_total);

  NFFT_SWAP_complex(ths->v_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  NFFT_SWAP_complex(ths->v_iter,ths->mv->f);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_v_iter = nfft_dot_w_complex(ths->v_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_v_iter = nfft_dot_complex(ths->v_iter, ths->mv->M_total);

  /*-----------------*/
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_upd_xpawy_complex(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->z_hat_iter, ths->mv->N_total);
  else
    nfft_upd_xpay_complex(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  nfft_upd_xpay_complex(ths->r_iter, -ths->alpha_iter, ths->v_iter,
			ths->mv->M_total);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = nfft_dot_w_complex(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = nfft_dot_complex(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_cp_w_complex(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    nfft_cp_complex(ths->mv->f, ths->r_iter, ths->mv->M_total);

  NFFT_SWAP_complex(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  NFFT_SWAP_complex(ths->z_hat_iter,ths->mv->f_hat);

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_z_hat_iter = nfft_dot_w_complex(ths->z_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_z_hat_iter = nfft_dot_complex(ths->z_hat_iter, ths->mv->N_total);
} /* void solver_loop_one_step_steepest_descent */

/** void solver_loop_one_step_cgnr */
void solver_loop_one_step_cgnr_complex(solver_plan_complex *ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_cp_w_complex(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,
		      ths->mv->N_total);
  else
    nfft_cp_complex(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);

  NFFT_SWAP_complex(ths->v_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  NFFT_SWAP_complex(ths->v_iter,ths->mv->f);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_v_iter = nfft_dot_w_complex(ths->v_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_v_iter = nfft_dot_complex(ths->v_iter, ths->mv->M_total);

  /*-----------------*/
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_upd_xpawy_complex(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->p_hat_iter, ths->mv->N_total);
  else
    nfft_upd_xpay_complex(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  nfft_upd_xpay_complex(ths->r_iter, -ths->alpha_iter, ths->v_iter,
			ths->mv->M_total);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = nfft_dot_w_complex(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = nfft_dot_complex(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_cp_w_complex(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    nfft_cp_complex(ths->mv->f, ths->r_iter, ths->mv->M_total);

  NFFT_SWAP_complex(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  NFFT_SWAP_complex(ths->z_hat_iter,ths->mv->f_hat);

  ths->dot_z_hat_iter_old = ths->dot_z_hat_iter;
  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_z_hat_iter = nfft_dot_w_complex(ths->z_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_z_hat_iter = nfft_dot_complex(ths->z_hat_iter, ths->mv->N_total);

  /*-----------------*/
  ths->beta_iter = ths->dot_z_hat_iter / ths->dot_z_hat_iter_old;

  /*-----------------*/
  nfft_upd_axpy_complex(ths->p_hat_iter, ths->beta_iter, ths->z_hat_iter,
			ths->mv->N_total);
} /* void solver_loop_one_step_cgnr */

/** void solver_loop_one_step_cgne */
void solver_loop_one_step_cgne_complex(solver_plan_complex *ths)
{
  ths->alpha_iter = ths->dot_r_iter / ths->dot_p_hat_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_upd_xpawy_complex(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->p_hat_iter, ths->mv->N_total);
  else
    nfft_upd_xpay_complex(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,
                          ths->mv->N_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_cp_w_complex(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,
		      ths->mv->N_total);
  else
    nfft_cp_complex(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);

  ths->mv->mv_trafo(ths->mv);

  nfft_upd_xpay_complex(ths->r_iter, -ths->alpha_iter, ths->mv->f,
			ths->mv->M_total);

  ths->dot_r_iter_old = ths->dot_r_iter;
  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = nfft_dot_w_complex(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = nfft_dot_complex(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  ths->beta_iter = ths->dot_r_iter / ths->dot_r_iter_old;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_cp_w_complex(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    nfft_cp_complex(ths->mv->f, ths->r_iter, ths->mv->M_total);

  ths->mv->mv_adjoint(ths->mv);

  nfft_upd_axpy_complex(ths->p_hat_iter, ths->beta_iter, ths->mv->f_hat,
			ths->mv->N_total);

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_p_hat_iter = nfft_dot_w_complex(ths->p_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_p_hat_iter = nfft_dot_complex(ths->p_hat_iter, ths->mv->N_total);
} /* void solver_loop_one_step_cgne */

/** void solver_loop_one_step */
void solver_loop_one_step_complex(solver_plan_complex *ths)
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
void solver_finalize_complex(solver_plan_complex *ths)
{
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_free(ths->w);

  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_free(ths->w_hat);

  if(ths->flags & CGNR)
    {
      nfft_free(ths->v_iter);
      nfft_free(ths->z_hat_iter);
    }

  if(ths->flags & STEEPEST_DESCENT)
    nfft_free(ths->v_iter);

  nfft_free(ths->p_hat_iter);
  nfft_free(ths->f_hat_iter);

  nfft_free(ths->r_iter);
  nfft_free(ths->y);
} /** void solver_finalize */


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void solver_init_advanced_double(solver_plan_double* ths, mv_plan_double *mv, unsigned flags)
{
  ths->mv = mv;
  ths->flags = flags;

  ths->y          = (double*)nfft_malloc(ths->mv->M_total*sizeof(double));
  ths->r_iter     = (double*)nfft_malloc(ths->mv->M_total*sizeof(double));
  ths->f_hat_iter = (double*)nfft_malloc(ths->mv->N_total*sizeof(double));
  ths->p_hat_iter = (double*)nfft_malloc(ths->mv->N_total*sizeof(double));

  if(ths->flags & LANDWEBER)
    ths->z_hat_iter = ths->p_hat_iter;

  if(ths->flags & STEEPEST_DESCENT)
    {
      ths->z_hat_iter = ths->p_hat_iter;
      ths->v_iter     = (double*)nfft_malloc(ths->mv->M_total*sizeof(double));
    }

  if(ths->flags & CGNR)
    {
      ths->z_hat_iter = (double*)nfft_malloc(ths->mv->N_total*sizeof(double));
      ths->v_iter     = (double*)nfft_malloc(ths->mv->M_total*sizeof(double));
    }

  if(ths->flags & CGNE)
    ths->z_hat_iter = ths->p_hat_iter;

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->w = (double*) nfft_malloc(ths->mv->M_total*sizeof(double));

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->w_hat = (double*) nfft_malloc(ths->mv->N_total*sizeof(double));
}

void solver_init_double(solver_plan_double* ths, mv_plan_double *mv)
{
  solver_init_advanced_double(ths, mv, CGNR);
}

void solver_before_loop_double(solver_plan_double* ths)
{
  nfft_cp_double(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);

  NFFT_SWAP_double(ths->r_iter, ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  NFFT_SWAP_double(ths->r_iter, ths->mv->f);

  nfft_upd_axpy_double(ths->r_iter, -1.0, ths->y, ths->mv->M_total);

  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))
    {
      if(ths->flags & PRECOMPUTE_WEIGHT)
	ths->dot_r_iter = nfft_dot_w_double(ths->r_iter, ths->w,
					     ths->mv->M_total);
      else
	ths->dot_r_iter = nfft_dot_double(ths->r_iter, ths->mv->M_total);
    }

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_cp_w_double(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    nfft_cp_double(ths->mv->f, ths->r_iter, ths->mv->M_total);

  NFFT_SWAP_double(ths->z_hat_iter, ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  NFFT_SWAP_double(ths->z_hat_iter, ths->mv->f_hat);

  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))
    {
      if(ths->flags & PRECOMPUTE_DAMP)
	ths->dot_z_hat_iter = nfft_dot_w_double(ths->z_hat_iter, ths->w_hat,
						 ths->mv->N_total);
      else
	ths->dot_z_hat_iter = nfft_dot_double(ths->z_hat_iter,
					       ths->mv->N_total);
    }

  if(ths->flags & CGNE)
    ths->dot_p_hat_iter = ths->dot_z_hat_iter;

  if(ths->flags & CGNR)
    nfft_cp_double(ths->p_hat_iter, ths->z_hat_iter, ths->mv->N_total);
} /* void solver_before_loop */

/** void solver_loop_one_step_landweber */
void solver_loop_one_step_landweber_double(solver_plan_double* ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_upd_xpawy_double(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->z_hat_iter, ths->mv->N_total);
  else
    nfft_upd_xpay_double(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  nfft_cp_double(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);

  NFFT_SWAP_double(ths->r_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  NFFT_SWAP_double(ths->r_iter,ths->mv->f);

  nfft_upd_axpy_double(ths->r_iter, -1.0, ths->y, ths->mv->M_total);

  if(ths->flags & NORMS_FOR_LANDWEBER)
    {
      if(ths->flags & PRECOMPUTE_WEIGHT)
	ths->dot_r_iter = nfft_dot_w_double(ths->r_iter,ths->w,
					     ths->mv->M_total);
      else
	ths->dot_r_iter = nfft_dot_double(ths->r_iter, ths->mv->M_total);
    }

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_cp_w_double(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    nfft_cp_double(ths->mv->f, ths->r_iter, ths->mv->M_total);

  NFFT_SWAP_double(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  NFFT_SWAP_double(ths->z_hat_iter,ths->mv->f_hat);

  if(ths->flags & NORMS_FOR_LANDWEBER)
    {
      if(ths->flags & PRECOMPUTE_DAMP)
	ths->dot_z_hat_iter = nfft_dot_w_double(ths->z_hat_iter, ths->w_hat,
						 ths->mv->N_total);
      else
	ths->dot_z_hat_iter = nfft_dot_double(ths->z_hat_iter,
					       ths->mv->N_total);
    }
} /* void solver_loop_one_step_landweber */

/** void solver_loop_one_step_steepest_descent */
void solver_loop_one_step_steepest_descent_double(solver_plan_double *ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_cp_w_double(ths->mv->f_hat, ths->w_hat, ths->z_hat_iter,
		      ths->mv->N_total);
  else
    nfft_cp_double(ths->mv->f_hat, ths->z_hat_iter, ths->mv->N_total);

  NFFT_SWAP_double(ths->v_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  NFFT_SWAP_double(ths->v_iter,ths->mv->f);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_v_iter = nfft_dot_w_double(ths->v_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_v_iter = nfft_dot_double(ths->v_iter, ths->mv->M_total);

  /*-----------------*/
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_upd_xpawy_double(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->z_hat_iter, ths->mv->N_total);
  else
    nfft_upd_xpay_double(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  nfft_upd_xpay_double(ths->r_iter, -ths->alpha_iter, ths->v_iter,
			ths->mv->M_total);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = nfft_dot_w_double(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = nfft_dot_double(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_cp_w_double(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    nfft_cp_double(ths->mv->f, ths->r_iter, ths->mv->M_total);

  NFFT_SWAP_double(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  NFFT_SWAP_double(ths->z_hat_iter,ths->mv->f_hat);

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_z_hat_iter = nfft_dot_w_double(ths->z_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_z_hat_iter = nfft_dot_double(ths->z_hat_iter, ths->mv->N_total);
} /* void solver_loop_one_step_steepest_descent */

/** void solver_loop_one_step_cgnr */
void solver_loop_one_step_cgnr_double(solver_plan_double *ths)
{
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_cp_w_double(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,
		      ths->mv->N_total);
  else
    nfft_cp_double(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);

  NFFT_SWAP_double(ths->v_iter,ths->mv->f);
  ths->mv->mv_trafo(ths->mv);
  NFFT_SWAP_double(ths->v_iter,ths->mv->f);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_v_iter = nfft_dot_w_double(ths->v_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_v_iter = nfft_dot_double(ths->v_iter, ths->mv->M_total);

  /*-----------------*/
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_upd_xpawy_double(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->p_hat_iter, ths->mv->N_total);
  else
    nfft_upd_xpay_double(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,
			  ths->mv->N_total);

  /*-----------------*/
  nfft_upd_xpay_double(ths->r_iter, -ths->alpha_iter, ths->v_iter,
			ths->mv->M_total);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = nfft_dot_w_double(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = nfft_dot_double(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_cp_w_double(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    nfft_cp_double(ths->mv->f, ths->r_iter, ths->mv->M_total);

  NFFT_SWAP_double(ths->z_hat_iter,ths->mv->f_hat);
  ths->mv->mv_adjoint(ths->mv);
  NFFT_SWAP_double(ths->z_hat_iter,ths->mv->f_hat);

  ths->dot_z_hat_iter_old = ths->dot_z_hat_iter;
  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_z_hat_iter = nfft_dot_w_double(ths->z_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_z_hat_iter = nfft_dot_double(ths->z_hat_iter, ths->mv->N_total);

  /*-----------------*/
  ths->beta_iter = ths->dot_z_hat_iter / ths->dot_z_hat_iter_old;

  /*-----------------*/
  nfft_upd_axpy_double(ths->p_hat_iter, ths->beta_iter, ths->z_hat_iter,
			ths->mv->N_total);
} /* void solver_loop_one_step_cgnr */

/** void solver_loop_one_step_cgne */
void solver_loop_one_step_cgne_double(solver_plan_double *ths)
{
  ths->alpha_iter = ths->dot_r_iter / ths->dot_p_hat_iter;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_upd_xpawy_double(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
			   ths->p_hat_iter, ths->mv->N_total);
  else
    nfft_upd_xpay_double(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,
                          ths->mv->N_total);

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_cp_w_double(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,
		      ths->mv->N_total);
  else
    nfft_cp_double(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);

  ths->mv->mv_trafo(ths->mv);

  nfft_upd_xpay_double(ths->r_iter, -ths->alpha_iter, ths->mv->f,
			ths->mv->M_total);

  ths->dot_r_iter_old = ths->dot_r_iter;
  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = nfft_dot_w_double(ths->r_iter,ths->w,ths->mv->M_total);
  else
    ths->dot_r_iter = nfft_dot_double(ths->r_iter, ths->mv->M_total);

  /*-----------------*/
  ths->beta_iter = ths->dot_r_iter / ths->dot_r_iter_old;

  /*-----------------*/
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_cp_w_double(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    nfft_cp_double(ths->mv->f, ths->r_iter, ths->mv->M_total);

  ths->mv->mv_adjoint(ths->mv);

  nfft_upd_axpy_double(ths->p_hat_iter, ths->beta_iter, ths->mv->f_hat,
			ths->mv->N_total);

  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_p_hat_iter = nfft_dot_w_double(ths->p_hat_iter, ths->w_hat,
					     ths->mv->N_total);
  else
    ths->dot_p_hat_iter = nfft_dot_double(ths->p_hat_iter, ths->mv->N_total);
} /* void solver_loop_one_step_cgne */

/** void solver_loop_one_step */
void solver_loop_one_step_double(solver_plan_double *ths)
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
void solver_finalize_double(solver_plan_double *ths)
{
  if(ths->flags & PRECOMPUTE_WEIGHT)
    nfft_free(ths->w);

  if(ths->flags & PRECOMPUTE_DAMP)
    nfft_free(ths->w_hat);

  if(ths->flags & CGNR)
    {
      nfft_free(ths->v_iter);
      nfft_free(ths->z_hat_iter);
    }

  if(ths->flags & STEEPEST_DESCENT)
    nfft_free(ths->v_iter);

  nfft_free(ths->p_hat_iter);
  nfft_free(ths->f_hat_iter);

  nfft_free(ths->r_iter);
  nfft_free(ths->y);
} /** void solver_finalize */
