/**
 * Library.
 * Includes iterative solution to the inverse problem.
 * authors: D. Potts, S. Kunis (c) 2002-2005
 */

#include "util.h"
#include "nfft3.h"


#define MACRO_SOLVER_IMPL(MV, FLT)			                      \
                                                                              \
F(MV, FLT, init_advanced, i ## MV ## _plan *ths, MV ## _plan *mv,             \
		     	  unsigned i ## MV ## _flags)		              \
{                                                                             \
  ths->mv = mv;                                                               \
  ths->flags = i ## MV ## _flags;                                             \
                                                                              \
  ths->y          = (FLT*)fftw_malloc(ths->mv->M_total*sizeof(FLT));          \
  ths->r_iter     = (FLT*)fftw_malloc(ths->mv->M_total*sizeof(FLT));          \
  ths->f_hat_iter = (FLT*)fftw_malloc(ths->mv->N_total*sizeof(FLT));          \
  ths->p_hat_iter = (FLT*)fftw_malloc(ths->mv->N_total*sizeof(FLT));          \
                                                                              \
  if(ths->flags & LANDWEBER)                                                  \
    ths->z_hat_iter = ths->p_hat_iter;                                        \
                                                                              \
  if(ths->flags & STEEPEST_DESCENT)                                           \
    {                                                                         \
      ths->z_hat_iter = ths->p_hat_iter;                                      \
      ths->v_iter     = (FLT*)fftw_malloc(ths->mv->M_total*sizeof(FLT));      \
    }                                                                         \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    {                                                                         \
      ths->z_hat_iter = (FLT*)fftw_malloc(ths->mv->N_total*sizeof(FLT));      \
      ths->v_iter     = (FLT*)fftw_malloc(ths->mv->M_total*sizeof(FLT));      \
    }                                                                         \
                                                                              \
  if(ths->flags & CGNE)                                                       \
    ths->z_hat_iter = ths->p_hat_iter;                                        \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->w = (double*) fftw_malloc(ths->mv->M_total*sizeof(double));          \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->w_hat = (double*) fftw_malloc(ths->mv->N_total*sizeof(double));      \
}                                                                             \
                                                                              \
/** void i<mv>_init */                                                        \
F(MV, FLT, init, i ## MV ## _plan *ths, MV ## _plan *mv)                      \
{                                                                             \
  i ## MV ## _init_advanced(ths, mv, CGNR);                                   \
} /* void i<mv>_init */                                                       \
                                                                              \
/** void i<mv>_before_loop */                                                 \
F(MV, FLT, before_loop,   i ## MV ## _plan *ths)                              \
{                                                                             \
  cp_ ## FLT(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);              \
                                                                              \
  SWAP_ ## FLT(ths->r_iter, ths->mv->f);                                      \
  MV ## _trafo(ths->mv);                                                      \
  SWAP_ ## FLT(ths->r_iter, ths->mv->f);                                      \
                                                                              \
  upd_axpy_ ## FLT(ths->r_iter, -1.0, ths->y, ths->mv->M_total);              \
                                                                              \
  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))       \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_WEIGHT)                                      \
	ths->dot_r_iter =                                                     \
	    dot_w_ ## FLT(ths->r_iter, ths->w, ths->mv->M_total);             \
      else                                                                    \
	ths->dot_r_iter = dot_ ## FLT(ths->r_iter, ths->mv->M_total);         \
    }                                                                         \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);          \
  else                                                                        \
    cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);                    \
                                                                              \
  SWAP_ ## FLT(ths->z_hat_iter, ths->mv->f_hat);                              \
  MV ## _adjoint(ths->mv);                                                    \
  SWAP_ ## FLT(ths->z_hat_iter, ths->mv->f_hat);                              \
                                                                              \
  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))       \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_DAMP)                                        \
	ths->dot_z_hat_iter =                                                 \
	  dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);       \
      else                                                                    \
	ths->dot_z_hat_iter = dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total); \
    }                                                                         \
                                                                              \
  if(ths->flags & CGNE)                                                       \
    ths->dot_p_hat_iter = ths->dot_z_hat_iter;                                \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    cp_ ## FLT(ths->p_hat_iter, ths->z_hat_iter, ths->mv->N_total);           \
} /* void i<mv>_before_loop */                                                \
                                                                              \
/** void i<mv>_loop_one_step_landweber */                                     \
F(MV, FLT, loop_one_step_landweber, i ## MV ## _plan *ths)                    \
{                                                                             \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,           \
		      ths->z_hat_iter, ths->mv->N_total);                     \
  else                                                                        \
    upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,       \
		     ths->mv->N_total);                                       \
                                                                              \
  /*-----------------*/                                                       \
  cp_ ## FLT(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);              \
                                                                              \
  SWAP_ ## FLT(ths->r_iter,ths->mv->f);                                       \
  MV ## _trafo(ths->mv);                                                      \
  SWAP_ ## FLT(ths->r_iter,ths->mv->f);                                       \
                                                                              \
  upd_axpy_ ## FLT(ths->r_iter, -1.0, ths->y, ths->mv->M_total);              \
                                                                              \
  if(ths->flags & NORMS_FOR_LANDWEBER)                                        \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_WEIGHT)                                      \
	ths->dot_r_iter = dot_w_ ## FLT(ths->r_iter,ths->w, ths->mv->M_total);\
      else                                                                    \
	ths->dot_r_iter = dot_ ## FLT(ths->r_iter, ths->mv->M_total);         \
    }                                                                         \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);          \
  else                                                                        \
    cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);                    \
                                                                              \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
  MV ## _adjoint(ths->mv);                                                    \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
                                                                              \
  if(ths->flags & NORMS_FOR_LANDWEBER)                                        \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_DAMP)                                        \
	ths->dot_z_hat_iter =                                                 \
	  dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);       \
      else                                                                    \
	ths->dot_z_hat_iter = dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total); \
    }                                                                         \
} /* void i<mv>_loop_one_step_landweber */                                    \
                                                                              \
/** void i<mv>_loop_one_step_steepest_descent */                              \
F(MV, FLT, loop_one_step_steepest_descent, i ## MV ## _plan *ths)             \
{                                                                             \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    cp_w_ ## FLT(ths->mv->f_hat, ths->w_hat, ths->z_hat_iter,                 \
		 ths->mv->N_total);                                           \
  else                                                                        \
    cp_ ## FLT(ths->mv->f_hat, ths->z_hat_iter, ths->mv->N_total);            \
                                                                              \
  SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                       \
  MV ## _trafo(ths->mv);                                                      \
  SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                       \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_v_iter = dot_w_ ## FLT(ths->v_iter, ths->w, ths->mv->M_total);   \
  else                                                                        \
    ths->dot_v_iter = dot_ ## FLT(ths->v_iter, ths->mv->M_total);             \
                                                                              \
  /*-----------------*/                                                       \
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;                    \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,           \
		      ths->z_hat_iter, ths->mv->N_total);                     \
  else                                                                        \
    upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,       \
		     ths->mv->N_total);                                       \
                                                                              \
  /*-----------------*/                                                       \
  upd_xpay_ ## FLT(ths->r_iter, -ths->alpha_iter, ths->v_iter,                \
		   ths->mv->M_total);                                         \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_r_iter = dot_w_ ## FLT(ths->r_iter, ths->w, ths->mv->M_total);   \
  else                                                                        \
    ths->dot_r_iter = dot_ ## FLT(ths->r_iter, ths->mv->M_total);             \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);          \
  else                                                                        \
    cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);                    \
                                                                              \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
  MV ## _adjoint(ths->mv);                                                    \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->dot_z_hat_iter =                                                     \
      dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);           \
  else                                                                        \
    ths->dot_z_hat_iter = dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total);     \
} /* void i<mv>_loop_one_step_steepest_descent */                             \
                                                                              \
/** void i<mv>_loop_one_step_cgnr */                                          \
F(MV, FLT, loop_one_step_cgnr, i ## MV ## _plan *ths)                         \
{                                                                             \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    cp_w_ ## FLT(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,                 \
		 ths->mv->N_total);                                           \
  else                                                                        \
    cp_ ## FLT(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);            \
                                                                              \
  SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                       \
  MV ## _trafo(ths->mv);                                                      \
  SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                       \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_v_iter = dot_w_ ## FLT(ths->v_iter, ths->w, ths->mv->M_total);   \
  else                                                                        \
    ths->dot_v_iter = dot_ ## FLT(ths->v_iter, ths->mv->M_total);             \
                                                                              \
  /*-----------------*/                                                       \
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;                    \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,           \
		      ths->p_hat_iter, ths->mv->N_total);                     \
  else                                                                        \
    upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,       \
		     ths->mv->N_total);                                       \
                                                                              \
  /*-----------------*/                                                       \
  upd_xpay_ ## FLT(ths->r_iter, -ths->alpha_iter, ths->v_iter,                \
		   ths->mv->M_total);                                         \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_r_iter = dot_w_ ## FLT(ths->r_iter, ths->w, ths->mv->M_total);   \
  else                                                                        \
    ths->dot_r_iter = dot_ ## FLT(ths->r_iter, ths->mv->M_total);             \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);          \
  else                                                                        \
    cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);                    \
                                                                              \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
  MV ## _adjoint(ths->mv);                                                    \
  SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                               \
                                                                              \
  ths->dot_z_hat_iter_old = ths->dot_z_hat_iter;                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->dot_z_hat_iter =                                                     \
      dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);           \
  else                                                                        \
    ths->dot_z_hat_iter = dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total);     \
                                                                              \
  /*-----------------*/                                                       \
  ths->beta_iter = ths->dot_z_hat_iter / ths->dot_z_hat_iter_old;             \
                                                                              \
  /*-----------------*/                                                       \
  upd_axpy_ ## FLT(ths->p_hat_iter, ths->beta_iter, ths->z_hat_iter,          \
		   ths->mv->N_total);                                         \
} /* void i<mv>_loop_one_step_cgnr */                                         \
                                                                              \
/** void i<mv>_loop_one_step_cgne */                                          \
F(MV, FLT, loop_one_step_cgne, i ## MV ## _plan *ths)                         \
{                                                                             \
  ths->alpha_iter = ths->dot_r_iter / ths->dot_p_hat_iter;                    \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,           \
		      ths->p_hat_iter, ths->mv->N_total);                     \
  else                                                                        \
    upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,       \
		     ths->mv->N_total);                                       \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    cp_w_ ## FLT(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,                 \
		 ths->mv->N_total);                                           \
  else                                                                        \
    cp_ ## FLT(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);            \
                                                                              \
  MV ## _trafo(ths->mv);                                                      \
                                                                              \
  upd_xpay_ ## FLT(ths->r_iter, -ths->alpha_iter, ths->mv->f,                 \
		   ths->mv->M_total);                                         \
                                                                              \
  ths->dot_r_iter_old = ths->dot_r_iter;                                      \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_r_iter = dot_w_ ## FLT(ths->r_iter, ths->w, ths->mv->M_total);   \
  else                                                                        \
    ths->dot_r_iter = dot_ ## FLT(ths->r_iter, ths->mv->M_total);             \
                                                                              \
  /*-----------------*/                                                       \
  ths->beta_iter = ths->dot_r_iter / ths->dot_r_iter_old;                     \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);          \
  else                                                                        \
    cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);                    \
                                                                              \
  MV ## _adjoint(ths->mv);                                                    \
                                                                              \
  upd_axpy_ ## FLT(ths->p_hat_iter, ths->beta_iter, ths->mv->f_hat,           \
		   ths->mv->N_total);                                         \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->dot_p_hat_iter =                                                     \
      dot_w_ ## FLT(ths->p_hat_iter, ths->w_hat, ths->mv->N_total);           \
  else                                                                        \
    ths->dot_p_hat_iter = dot_ ## FLT(ths->p_hat_iter, ths->mv->N_total);     \
} /* void i<mv>_loop_one_step_cgne */                                         \
                                                                              \
/** void i<mv>_loop_one_step */                                               \
F(MV, FLT, loop_one_step, i ## MV ## _plan *ths)                              \
{                                                                             \
  if(ths->flags & LANDWEBER)                                                  \
    i ## MV ## _loop_one_step_landweber(ths);                                 \
                                                                              \
  if(ths->flags & STEEPEST_DESCENT)                                           \
    i ## MV ## _loop_one_step_steepest_descent(ths);                          \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    i ## MV ## _loop_one_step_cgnr(ths);                                      \
                                                                              \
  if(ths->flags & CGNE)                                                       \
    i ## MV ## _loop_one_step_cgne(ths);                                      \
} /* void i<mv>_loop_one_step */                                              \
                                                                              \
/** void i<mv>_finalize */                                                    \
F(MV, FLT, finalize, i ## MV ## _plan *ths)                                   \
{                                                                             \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    fftw_free(ths->w);                                                        \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    fftw_free(ths->w_hat);                                                    \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    {                                                                         \
      fftw_free(ths->v_iter);                                                 \
      fftw_free(ths->z_hat_iter);                                             \
    }                                                                         \
                                                                              \
  if(ths->flags & STEEPEST_DESCENT)                                           \
    fftw_free(ths->v_iter);                                                   \
                                                                              \
  fftw_free(ths->p_hat_iter);                                                 \
  fftw_free(ths->f_hat_iter);                                                 \
                                                                              \
  fftw_free(ths->r_iter);                                                     \
  fftw_free(ths->y);                                                          \
} /** void i<mv>_finalize */

MACRO_SOLVER_IMPL(nfft, complex)
MACRO_SOLVER_IMPL(nnfft, complex)
MACRO_SOLVER_IMPL(mri_inh_2d1d, complex)
MACRO_SOLVER_IMPL(mri_inh_3d, complex)
//MACRO_SOLVER_IMPL(texture, complex)
