/**
 * Library.
 * Includes iterative solution to the inverse problem.
 * authors: D. Potts, S. Kunis (c) 2002-2005
 */

#include "solver.h"
#include "util.h"


void infft_init_specific(infft_plan *ths, nfft_plan *mv, int infft_flags)
{
  ths->mv = mv;
  ths->flags = infft_flags;
  
  ths->y = (complex*)fftw_malloc(ths->mv->M_total*sizeof(complex));

  ths->r_iter = (complex*)fftw_malloc(ths->mv->M_total*sizeof(complex));

  ths->f_hat_iter = (complex*)fftw_malloc(ths->mv->N_total*sizeof(complex));

  ths->p_hat_iter = (complex*)fftw_malloc(ths->mv->N_total*sizeof(complex));

  if(ths->flags & LANDWEBER)
    ths->z_hat_iter = ths->p_hat_iter;

  if(ths->flags & STEEPEST_DESCENT)
    {
      ths->z_hat_iter = ths->p_hat_iter;

      ths->v_iter = (complex*)fftw_malloc(ths->mv->M_total*sizeof(complex));
    }

  if(ths->flags & CGNR)
    {
      ths->z_hat_iter =(complex*)fftw_malloc(ths->mv->N_total*sizeof(complex));

      ths->v_iter = (complex*)fftw_malloc(ths->mv->M_total*sizeof(complex));
    }

  if(ths->flags & CGNE)
    ths->z_hat_iter = ths->p_hat_iter;

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->w = (double*) fftw_malloc(ths->mv->M_total*sizeof(double));
  
  if(ths->flags & PRECOMPUTE_DAMP)
    ths->w_hat = (double*) fftw_malloc(ths->mv->N_total*sizeof(double));
}

void infft_init(infft_plan *ths, nfft_plan *mv)
{
  infft_init_specific(ths, mv, CGNR);
}

void infft_before_loop(infft_plan *ths)
{
  /** step 2
   *  overwrites ths->mv->f_hat 
   *  overwrites ths->r_iter
   */
  copyc(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);
  
  SWAPC(ths->r_iter, ths->mv->f);
  nfft_trafo(ths->mv);
  SWAPC(ths->r_iter, ths->mv->f);

  updatec_axpy(ths->r_iter, -1.0, ths->y, ths->mv->M_total);

  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))
    {
      if(ths->flags & PRECOMPUTE_WEIGHT)
	ths->dot_r_iter =
	  dotproductc_w(ths->r_iter, ths->w, ths->mv->M_total);
      else
	ths->dot_r_iter = dotproductc(ths->r_iter, ths->mv->M_total);
    }
  
  /** step 3
   *  overwrites ths->mv->f
   *  overwrites ths->z_hat_iter resp. ths->z_hat_iter
   */ 

  if(ths->flags & PRECOMPUTE_WEIGHT)
    copyc_w(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    copyc(ths->mv->f, ths->r_iter, ths->mv->M_total);
  
  SWAPC(ths->z_hat_iter, ths->mv->f_hat);
  nfft_adjoint(ths->mv);
  SWAPC(ths->z_hat_iter, ths->mv->f_hat);
  
  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))
    {
      if(ths->flags & PRECOMPUTE_DAMP)
	ths->dot_z_hat_iter = 
	  dotproductc_w(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);
      else
	ths->dot_z_hat_iter = dotproductc(ths->z_hat_iter, ths->mv->N_total);
    }

  if(ths->flags & CGNE)
    ths->dot_p_hat_iter = ths->dot_z_hat_iter;

  if(ths->flags & CGNR)
    copyc(ths->p_hat_iter, ths->z_hat_iter, ths->mv->N_total);
} /* void infft_before_loop */

void infft_loop_one_step_landweber(infft_plan *ths)
{
  /** step 5
   *  updates ths->f_hat_iter
   */
  if(ths->flags & PRECOMPUTE_DAMP)
    updatec_xpawy(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
		  ths->z_hat_iter, ths->mv->N_total);
  else
    updatec_xpay(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,
		 ths->mv->N_total);
  
  /** step 6
   *  original residual, not the updated residual,
   *  overwrites ths->r_iter
   *  overwrites ths->mv->f_hat 
   */
  copyc(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);

  SWAPC(ths->r_iter,ths->mv->f);
  nfft_trafo(ths->mv);
  SWAPC(ths->r_iter,ths->mv->f);
  
  updatec_axpy(ths->r_iter, -1.0, ths->y, ths->mv->M_total);

  if(ths->flags & NORMS_FOR_LANDWEBER)
    {
      if(ths->flags & PRECOMPUTE_WEIGHT)
	ths->dot_r_iter = dotproductc_w(ths->r_iter,ths->w, ths->mv->M_total);
      else
	ths->dot_r_iter = dotproductc(ths->r_iter, ths->mv->M_total);
    }

  /** step 7
   *  overwrites ths->mv->f 
   *  overwrites ths->z_hat_iter
   */
  if(ths->flags & PRECOMPUTE_WEIGHT)
    copyc_w(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    copyc(ths->mv->f, ths->r_iter, ths->mv->M_total);
    
  SWAPC(ths->z_hat_iter,ths->mv->f_hat);
  nfft_adjoint(ths->mv);
  SWAPC(ths->z_hat_iter,ths->mv->f_hat);

  if(ths->flags & NORMS_FOR_LANDWEBER)
    {
      if(ths->flags & PRECOMPUTE_DAMP)
	ths->dot_z_hat_iter = 
	  dotproductc_w(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);
      else
	ths->dot_z_hat_iter = dotproductc(ths->z_hat_iter, ths->mv->N_total);
    }
} /* void infft_loop_one_step_landweber */

void infft_loop_one_step_steepest_descent(infft_plan *ths)
{
  /** step 5
   *  overwrites ths->mv->f_hat 
   *  overwrites ths->v_iter
   */
  if(ths->flags & PRECOMPUTE_DAMP)
    copyc_w(ths->mv->f_hat, ths->w_hat, ths->z_hat_iter, ths->mv->N_total);
  else
    copyc(ths->mv->f_hat, ths->z_hat_iter, ths->mv->N_total);
  
  SWAPC(ths->v_iter,ths->mv->f);
  nfft_trafo(ths->mv);
  SWAPC(ths->v_iter,ths->mv->f);
  
  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_v_iter = dotproductc_w(ths->v_iter, ths->w, ths->mv->M_total);
  else
    ths->dot_v_iter = dotproductc(ths->v_iter, ths->mv->M_total);
  
  /** step 6
   */
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;

  /** step 7
   *  updates ths->f_hat_iter
   */
  if(ths->flags & PRECOMPUTE_DAMP)
    updatec_xpawy(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
		  ths->z_hat_iter, ths->mv->N_total);
  else
    updatec_xpay(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,
		 ths->mv->N_total);

  /** step 8
   *  updates ths->r_iter
   */
  updatec_xpay(ths->r_iter, -ths->alpha_iter, ths->v_iter, ths->mv->M_total);

  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = dotproductc_w(ths->r_iter, ths->w, ths->mv->M_total);
  else
    ths->dot_r_iter = dotproductc(ths->r_iter, ths->mv->M_total);

  /** step 9
   *  overwrites ths->mv->f
   *  overwrites ths->z_hat_iter
   */
  if(ths->flags & PRECOMPUTE_WEIGHT)
    copyc_w(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    copyc(ths->mv->f, ths->r_iter, ths->mv->M_total);
  
  SWAPC(ths->z_hat_iter,ths->mv->f_hat);
  nfft_adjoint(ths->mv);
  SWAPC(ths->z_hat_iter,ths->mv->f_hat);
  
  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_z_hat_iter =
      dotproductc_w(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);
  else
    ths->dot_z_hat_iter = dotproductc(ths->z_hat_iter, ths->mv->N_total);
} /* void infft_loop_one_step_steepest_descent */

void infft_loop_one_step_cgnr(infft_plan *ths)
{
  /** step 9
   *  overwrites ths->mv->f_hat 
   *  overwrites ths->v_iter
   */
  if(ths->flags & PRECOMPUTE_DAMP)
    copyc_w(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter, ths->mv->N_total);
  else
    copyc(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);
  
  SWAPC(ths->v_iter,ths->mv->f);
  nfft_trafo(ths->mv);
  SWAPC(ths->v_iter,ths->mv->f);
  
  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_v_iter = dotproductc_w(ths->v_iter, ths->w, ths->mv->M_total);
  else
    ths->dot_v_iter = dotproductc(ths->v_iter, ths->mv->M_total);
  
  /** step 10
   */
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;

  /** step 11
   *  updates ths->f_hat_iter
   */
  if(ths->flags & PRECOMPUTE_DAMP)
    updatec_xpawy(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
		  ths->p_hat_iter, ths->mv->N_total);
  else
    updatec_xpay(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,
		 ths->mv->N_total);
  
  /** step 16
   *  updates ths->r_iter
   */
  updatec_xpay(ths->r_iter, -ths->alpha_iter, ths->v_iter, ths->mv->M_total);
  
  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = dotproductc_w(ths->r_iter, ths->w, ths->mv->M_total);
  else
    ths->dot_r_iter = dotproductc(ths->r_iter, ths->mv->M_total);

  /** step 17
   *  overwrites ths->mv->f
   *  overwrites ths->mv->r_iter
   */
  if(ths->flags & PRECOMPUTE_WEIGHT)
    copyc_w(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    copyc(ths->mv->f, ths->r_iter, ths->mv->M_total);
  
  SWAPC(ths->z_hat_iter,ths->mv->f_hat);
  nfft_adjoint(ths->mv);
  SWAPC(ths->z_hat_iter,ths->mv->f_hat);

  ths->dot_z_hat_iter_old = ths->dot_z_hat_iter;
  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_z_hat_iter =
      dotproductc_w(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);
  else
    ths->dot_z_hat_iter = dotproductc(ths->z_hat_iter, ths->mv->N_total);
  
  /** step 18
   */
  ths->beta_iter = ths->dot_z_hat_iter / ths->dot_z_hat_iter_old;
  
  /** step 19
   *  updates ths->p_hat_iter
   */
  updatec_axpy(ths->p_hat_iter, ths->beta_iter, ths->z_hat_iter,
	       ths->mv->N_total);
} /* void infft_loop_one_step_cgnr */

void infft_loop_one_step_cgne(infft_plan *ths)
{
  /** step 9
   */
  ths->alpha_iter = ths->dot_r_iter / ths->dot_p_hat_iter;
  
  /** step 10
   *  updates ths->f_hat_iter
   */
  if(ths->flags & PRECOMPUTE_DAMP)
    updatec_xpawy(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,
		  ths->p_hat_iter, ths->mv->N_total);
  else
    updatec_xpay(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,
		 ths->mv->N_total);
  
  /** step 11
   *  overwrites ths->mv->f_hat 
   *  overwrites ths->mv->f
   *  updates ths->r_iter
   */
  if(ths->flags & PRECOMPUTE_DAMP)
    copyc_w(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter, ths->mv->N_total);
  else
    copyc(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);
  
  nfft_trafo(ths->mv);

  updatec_xpay(ths->r_iter, -ths->alpha_iter, ths->mv->f, ths->mv->M_total);
  
  ths->dot_r_iter_old = ths->dot_r_iter;
  if(ths->flags & PRECOMPUTE_WEIGHT)
    ths->dot_r_iter = dotproductc_w(ths->r_iter, ths->w, ths->mv->M_total);
  else
    ths->dot_r_iter = dotproductc(ths->r_iter, ths->mv->M_total);
  
  /** step 12
   */
  ths->beta_iter = ths->dot_r_iter / ths->dot_r_iter_old;
    
  /** step 16
   *  overwrites ths->mv->f
   *  overwrites ths->mv->f_hat
   *  updates ths->p_hat_iter
   */
  if(ths->flags & PRECOMPUTE_WEIGHT)
    copyc_w(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);
  else
    copyc(ths->mv->f, ths->r_iter, ths->mv->M_total); 
  
  nfft_adjoint(ths->mv);
  
  updatec_axpy(ths->p_hat_iter, ths->beta_iter, ths->mv->f_hat,
	       ths->mv->N_total);
  
  if(ths->flags & PRECOMPUTE_DAMP)
    ths->dot_p_hat_iter = 
      dotproductc_w(ths->p_hat_iter, ths->w_hat, ths->mv->N_total);
  else
    ths->dot_p_hat_iter = dotproductc(ths->p_hat_iter, ths->mv->N_total);
} /* void infft_loop_one_step_cgne */

void infft_loop_one_step(infft_plan *ths)
{
  if(ths->flags & LANDWEBER)
    infft_loop_one_step_landweber(ths);

  if(ths->flags & STEEPEST_DESCENT)
    infft_loop_one_step_steepest_descent(ths);
  
  if(ths->flags & CGNR)
    infft_loop_one_step_cgnr(ths);
  
  if(ths->flags & CGNE)
    infft_loop_one_step_cgne(ths);
}

void infft_finalize(infft_plan *ths)
{
  if(ths->flags & PRECOMPUTE_WEIGHT)
    fftw_free(ths->w);
  
  if(ths->flags & PRECOMPUTE_DAMP)
    fftw_free(ths->w_hat);
  
  if(ths->flags & CGNR)
    {
      fftw_free(ths->v_iter);
      fftw_free(ths->z_hat_iter);
    }
  
  if(ths->flags & STEEPEST_DESCENT)
    fftw_free(ths->v_iter);
  
  fftw_free(ths->p_hat_iter);
  fftw_free(ths->f_hat_iter);

  fftw_free(ths->r_iter);
  fftw_free(ths->y);
}
