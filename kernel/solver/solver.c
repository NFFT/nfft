/* $Id$
 *
 * Copyright (c) 2005, 2008 Jens Keiner, Stefan Kunis, Daniel Potts
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

/**
 * Library.
 * Includes iterative solution to the inverse problem.
 * authors: D. Potts, S. Kunis (c) 2002-2006
 */

#include <complex.h>

#include "util.h"
#include "nfft3.h"

/**
 * Function mangling macro.
 *
 * \arg MV Matrix vector multiplication type (nfft, nfct, ...)
 * \arg FLT Float used as prefix for function names (double or double _Complex)
 * \arg FLT_TYPE Float type (double or double _Complex)
 * \arg NAME Name of the functions suffix
 * \arg ...  argument list of the function
 *
 * \author Stefan Kunis
 */
#define F(MV, FLT, FLT_TYPE, NAME, ...) void i ## MV ## _ ## NAME(__VA_ARGS__)


#define MACRO_SOLVER_IMPL(MV, FLT, FLT_TYPE)                                  \
                                                                              \
F(MV, FLT, FLT_TYPE, init_advanced, i ## MV ## _plan *ths, MV ## _plan *mv,   \
                                    unsigned i ## MV ## _flags)               \
{                                                                             \
  ths->mv = mv;                                                               \
  ths->flags = i ## MV ## _flags;                                             \
                                                                              \
  ths->y          = (FLT_TYPE*)fftw_malloc(ths->mv->M_total*sizeof(FLT_TYPE));\
  ths->r_iter     = (FLT_TYPE*)fftw_malloc(ths->mv->M_total*sizeof(FLT_TYPE));\
  ths->f_hat_iter = (FLT_TYPE*)fftw_malloc(ths->mv->N_total*sizeof(FLT_TYPE));\
  ths->p_hat_iter = (FLT_TYPE*)fftw_malloc(ths->mv->N_total*sizeof(FLT_TYPE));\
                                                                              \
  if(ths->flags & LANDWEBER)                                                  \
    ths->z_hat_iter = ths->p_hat_iter;                                        \
                                                                              \
  if(ths->flags & STEEPEST_DESCENT)                                           \
    {                                                                         \
      ths->z_hat_iter = ths->p_hat_iter;                                      \
      ths->v_iter     = (FLT_TYPE*)fftw_malloc(ths->mv->M_total*              \
                                               sizeof(FLT_TYPE));             \
    }                                                                         \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    {                                                                         \
      ths->z_hat_iter = (FLT_TYPE*)fftw_malloc(ths->mv->N_total*              \
                                               sizeof(FLT_TYPE));             \
      ths->v_iter     = (FLT_TYPE*)fftw_malloc(ths->mv->M_total*              \
                                               sizeof(FLT_TYPE));             \
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
F(MV, FLT, FLT_TYPE, init, i ## MV ## _plan *ths, MV ## _plan *mv)            \
{                                                                             \
  i ## MV ## _init_advanced(ths, mv, CGNR);                                   \
} /* void i<mv>_init */                                                       \
                                                                              \
/** void i<mv>_before_loop */                                                 \
F(MV, FLT, FLT_TYPE, before_loop,   i ## MV ## _plan *ths)                    \
{                                                                             \
  nfft_cp_ ## FLT(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);         \
                                                                              \
  NFFT_SWAP_ ## FLT(ths->r_iter, ths->mv->f);                                 \
  MV ## _trafo(ths->mv);                                                      \
  NFFT_SWAP_ ## FLT(ths->r_iter, ths->mv->f);                                 \
                                                                              \
  nfft_upd_axpy_ ## FLT(ths->r_iter, -1.0, ths->y, ths->mv->M_total);         \
                                                                              \
  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))       \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_WEIGHT)                                      \
  ths->dot_r_iter =                                                           \
      nfft_dot_w_ ## FLT(ths->r_iter, ths->w, ths->mv->M_total);              \
      else                                                                    \
  ths->dot_r_iter = nfft_dot_ ## FLT(ths->r_iter, ths->mv->M_total);          \
    }                                                                         \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    nfft_cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);     \
  else                                                                        \
    nfft_cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);               \
                                                                              \
  NFFT_SWAP_ ## FLT(ths->z_hat_iter, ths->mv->f_hat);                         \
  MV ## _adjoint(ths->mv);                                                    \
  NFFT_SWAP_ ## FLT(ths->z_hat_iter, ths->mv->f_hat);                         \
                                                                              \
  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))       \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_DAMP)                                        \
  ths->dot_z_hat_iter =                                                       \
    nfft_dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);        \
      else                                                                    \
  ths->dot_z_hat_iter = nfft_dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total);  \
    }                                                                         \
                                                                              \
  if(ths->flags & CGNE)                                                       \
    ths->dot_p_hat_iter = ths->dot_z_hat_iter;                                \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    nfft_cp_ ## FLT(ths->p_hat_iter, ths->z_hat_iter, ths->mv->N_total);      \
} /* void i<mv>_before_loop */                                                \
                                                                              \
/** void i<mv>_loop_one_step_landweber */                                     \
F(MV, FLT, FLT_TYPE, loop_one_step_landweber, i ## MV ## _plan *ths)          \
{                                                                             \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    nfft_upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,      \
          ths->z_hat_iter, ths->mv->N_total);                                 \
  else                                                                        \
    nfft_upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,  \
         ths->mv->N_total);                                                   \
                                                                              \
  /*-----------------*/                                                       \
  nfft_cp_ ## FLT(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total);         \
                                                                              \
  NFFT_SWAP_ ## FLT(ths->r_iter,ths->mv->f);                                  \
  MV ## _trafo(ths->mv);                                                      \
  NFFT_SWAP_ ## FLT(ths->r_iter,ths->mv->f);                                  \
                                                                              \
  nfft_upd_axpy_ ## FLT(ths->r_iter, -1.0, ths->y, ths->mv->M_total);         \
                                                                              \
  if(ths->flags & NORMS_FOR_LANDWEBER)                                        \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_WEIGHT)                                      \
  ths->dot_r_iter = nfft_dot_w_ ## FLT(ths->r_iter,ths->w, ths->mv->M_total); \
      else                                                                    \
  ths->dot_r_iter = nfft_dot_ ## FLT(ths->r_iter, ths->mv->M_total);          \
    }                                                                         \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    nfft_cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);     \
  else                                                                        \
    nfft_cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);               \
                                                                              \
  NFFT_SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                          \
  MV ## _adjoint(ths->mv);                                                    \
  NFFT_SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                          \
                                                                              \
  if(ths->flags & NORMS_FOR_LANDWEBER)                                        \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_DAMP)                                        \
  ths->dot_z_hat_iter =                                                       \
    nfft_dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);        \
      else                                                                    \
  ths->dot_z_hat_iter = nfft_dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total);  \
    }                                                                         \
} /* void i<mv>_loop_one_step_landweber */                                    \
                                                                              \
/** void i<mv>_loop_one_step_steepest_descent */                              \
F(MV, FLT, FLT_TYPE, loop_one_step_steepest_descent, i ## MV ## _plan *ths)   \
{                                                                             \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    nfft_cp_w_ ## FLT(ths->mv->f_hat, ths->w_hat, ths->z_hat_iter,            \
     ths->mv->N_total);                                                       \
  else                                                                        \
    nfft_cp_ ## FLT(ths->mv->f_hat, ths->z_hat_iter, ths->mv->N_total);       \
                                                                              \
  NFFT_SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                  \
  MV ## _trafo(ths->mv);                                                      \
  NFFT_SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                  \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_v_iter = nfft_dot_w_ ## FLT(ths->v_iter, ths->w,                 \
                                         ths->mv->M_total);                   \
  else                                                                        \
    ths->dot_v_iter = nfft_dot_ ## FLT(ths->v_iter, ths->mv->M_total);        \
                                                                              \
  /*-----------------*/                                                       \
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;                    \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    nfft_upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,      \
          ths->z_hat_iter, ths->mv->N_total);                                 \
  else                                                                        \
    nfft_upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->z_hat_iter,  \
         ths->mv->N_total);                                                   \
                                                                              \
  /*-----------------*/                                                       \
  nfft_upd_xpay_ ## FLT(ths->r_iter, -ths->alpha_iter, ths->v_iter,           \
       ths->mv->M_total);                                                     \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_r_iter = nfft_dot_w_ ## FLT(ths->r_iter, ths->w,                 \
                                         ths->mv->M_total);                   \
  else                                                                        \
    ths->dot_r_iter = nfft_dot_ ## FLT(ths->r_iter, ths->mv->M_total);        \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    nfft_cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);     \
  else                                                                        \
    nfft_cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);               \
                                                                              \
  NFFT_SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                          \
  MV ## _adjoint(ths->mv);                                                    \
  NFFT_SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                          \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->dot_z_hat_iter =                                                     \
      nfft_dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);      \
  else                                                                        \
    ths->dot_z_hat_iter = nfft_dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total);\
} /* void i<mv>_loop_one_step_steepest_descent */                             \
                                                                              \
/** void i<mv>_loop_one_step_cgnr */                                          \
F(MV, FLT, FLT_TYPE, loop_one_step_cgnr, i ## MV ## _plan *ths)               \
{                                                                             \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    nfft_cp_w_ ## FLT(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,            \
     ths->mv->N_total);                                                       \
  else                                                                        \
    nfft_cp_ ## FLT(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);       \
                                                                              \
  NFFT_SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                  \
  MV ## _trafo(ths->mv);                                                      \
  NFFT_SWAP_ ## FLT(ths->v_iter,ths->mv->f);                                  \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_v_iter = nfft_dot_w_ ## FLT(ths->v_iter, ths->w,                 \
                                         ths->mv->M_total);                   \
  else                                                                        \
    ths->dot_v_iter = nfft_dot_ ## FLT(ths->v_iter, ths->mv->M_total);        \
                                                                              \
  /*-----------------*/                                                       \
  ths->alpha_iter = ths->dot_z_hat_iter / ths->dot_v_iter;                    \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    nfft_upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,      \
          ths->p_hat_iter, ths->mv->N_total);                                 \
  else                                                                        \
    nfft_upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,  \
         ths->mv->N_total);                                                   \
                                                                              \
  /*-----------------*/                                                       \
  nfft_upd_xpay_ ## FLT(ths->r_iter, -ths->alpha_iter, ths->v_iter,           \
       ths->mv->M_total);                                                     \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_r_iter = nfft_dot_w_ ## FLT(ths->r_iter, ths->w,                 \
                                         ths->mv->M_total);                   \
  else                                                                        \
    ths->dot_r_iter = nfft_dot_ ## FLT(ths->r_iter, ths->mv->M_total);        \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    nfft_cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);     \
  else                                                                        \
    nfft_cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);               \
                                                                              \
  NFFT_SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                          \
  MV ## _adjoint(ths->mv);                                                    \
  NFFT_SWAP_ ## FLT(ths->z_hat_iter,ths->mv->f_hat);                          \
                                                                              \
  ths->dot_z_hat_iter_old = ths->dot_z_hat_iter;                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->dot_z_hat_iter =                                                     \
      nfft_dot_w_ ## FLT(ths->z_hat_iter, ths->w_hat, ths->mv->N_total);      \
  else                                                                        \
    ths->dot_z_hat_iter = nfft_dot_ ## FLT(ths->z_hat_iter, ths->mv->N_total);\
                                                                              \
  /*-----------------*/                                                       \
  ths->beta_iter = ths->dot_z_hat_iter / ths->dot_z_hat_iter_old;             \
                                                                              \
  /*-----------------*/                                                       \
  nfft_upd_axpy_ ## FLT(ths->p_hat_iter, ths->beta_iter, ths->z_hat_iter,     \
       ths->mv->N_total);                                                     \
} /* void i<mv>_loop_one_step_cgnr */                                         \
                                                                              \
/** void i<mv>_loop_one_step_cgne */                                          \
F(MV, FLT, FLT_TYPE, loop_one_step_cgne, i ## MV ## _plan *ths)               \
{                                                                             \
  ths->alpha_iter = ths->dot_r_iter / ths->dot_p_hat_iter;                    \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    nfft_upd_xpawy_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->w_hat,      \
          ths->p_hat_iter, ths->mv->N_total);                                 \
  else                                                                        \
    nfft_upd_xpay_ ## FLT(ths->f_hat_iter, ths->alpha_iter, ths->p_hat_iter,  \
                          ths->mv->N_total);                                  \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    nfft_cp_w_ ## FLT(ths->mv->f_hat, ths->w_hat, ths->p_hat_iter,            \
                      ths->mv->N_total);                                      \
  else                                                                        \
    nfft_cp_ ## FLT(ths->mv->f_hat, ths->p_hat_iter, ths->mv->N_total);       \
                                                                              \
  MV ## _trafo(ths->mv);                                                      \
                                                                              \
  nfft_upd_xpay_ ## FLT(ths->r_iter, -ths->alpha_iter, ths->mv->f,            \
                        ths->mv->M_total);                                    \
                                                                              \
  ths->dot_r_iter_old = ths->dot_r_iter;                                      \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_r_iter = nfft_dot_w_ ## FLT(ths->r_iter, ths->w,                 \
                                         ths->mv->M_total);                   \
  else                                                                        \
    ths->dot_r_iter = nfft_dot_ ## FLT(ths->r_iter, ths->mv->M_total);        \
                                                                              \
  /*-----------------*/                                                       \
  ths->beta_iter = ths->dot_r_iter / ths->dot_r_iter_old;                     \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    nfft_cp_w_ ## FLT(ths->mv->f, ths->w, ths->r_iter, ths->mv->M_total);     \
  else                                                                        \
    nfft_cp_ ## FLT(ths->mv->f, ths->r_iter, ths->mv->M_total);               \
                                                                              \
  MV ## _adjoint(ths->mv);                                                    \
                                                                              \
  nfft_upd_axpy_ ## FLT(ths->p_hat_iter, ths->beta_iter, ths->mv->f_hat,      \
       ths->mv->N_total);                                                     \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->dot_p_hat_iter =                                                     \
      nfft_dot_w_ ## FLT(ths->p_hat_iter, ths->w_hat, ths->mv->N_total);      \
  else                                                                        \
    ths->dot_p_hat_iter = nfft_dot_ ## FLT(ths->p_hat_iter, ths->mv->N_total);\
} /* void i<mv>_loop_one_step_cgne */                                         \
                                                                              \
/** void i<mv>_loop_one_step */                                               \
F(MV, FLT, FLT_TYPE, loop_one_step, i ## MV ## _plan *ths)                    \
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
F(MV, FLT, FLT_TYPE, finalize, i ## MV ## _plan *ths)                         \
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

MACRO_SOLVER_IMPL(nfft, complex, double _Complex)
#ifdef HAVE_NFCT
MACRO_SOLVER_IMPL(nfct, double, double)
#endif
#ifdef HAVE_NFST
MACRO_SOLVER_IMPL(nfst, double, double)
#endif
#ifdef HAVE_NNFFT
MACRO_SOLVER_IMPL(nnfft, complex, double _Complex)
#endif
#ifdef HAVE_MRI
MACRO_SOLVER_IMPL(mri_inh_2d1d, complex, double _Complex)
MACRO_SOLVER_IMPL(mri_inh_3d, complex, double _Complex)
#endif
#ifdef HAVE_NFSFT
MACRO_SOLVER_IMPL(nfsft, complex, double _Complex)
#endif
