/*
 * $Id$
 *
 * Copyright (c) 2002, 2009 Jens Keiner, Daniel Potts, Stefan Kunis
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
* Library to solve the adjoint problem tim becker
*/

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
 * \author Stefan Kunis
 */

 /**
  *  TAKE CARE
  *  M_total  =  number of points
  *  N_total  =  bandwidth
  */



#define F(MV, FLT, FLT_TYPE, NAME, ...) void i ## MV ## _ ## NAME(__VA_ARGS__)


#define MACRO_SOLVER_ADJOINT_IMPL(MV, FLT, FLT_TYPE)                          \
                                                                              \
F(MV, FLT, FLT_TYPE, init_advanced, i ## MV ## _adjoint_plan *ths,             \
MV ## _adjoint_plan *mv, unsigned adjoint ## MV ## _flags)                     \
{                                                                             \
  ths->mv = mv;                                                               \
  ths->flags = i ## MV ## _flags;                                             \
                                                                              \
  ths->y_hat      = (FLT_TYPE*)nfft_malloc(ths->mv->N_total*sizeof(FLT_TYPE));\
  ths->r_hat_iter = (FLT_TYPE*)nfft_malloc(ths->mv->N_total*sizeof(FLT_TYPE));\
  ths->f_iter = (FLT_TYPE*)nfft_malloc(ths->mv->M_total*sizeof(FLT_TYPE));    \
  ths->p_iter = (FLT_TYPE*)nfft_malloc(ths->mv->M_total*sizeof(FLT_TYPE));    \
                                                                              \
  if(ths->flags & LANDWEBER)                                                  \
    ths->z_iter = ths->p_iter;                                                \
                                                                              \
  if(ths->flags & STEEPEST_DESCENT)                                           \
    {                                                                         \
      ths->z_iter = ths->p_iter;                                              \
      ths->v_hat_iter     = (FLT_TYPE*)nfft_malloc(ths->mv->N_total*          \
                                               sizeof(FLT_TYPE));             \
    }                                                                         \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    {                                                                         \
      ths->z_iter = (FLT_TYPE*)nfft_malloc(ths->mv->M_total*                  \
                                               sizeof(FLT_TYPE));             \
      ths->v_hat_iter     = (FLT_TYPE*)nfft_malloc(ths->mv->N_total*          \
                                               sizeof(FLT_TYPE));             \
    }                                                                         \
                                                                              \
  if(ths->flags & CGNE)                                                       \
    ths->z_iter = ths->p_iter;                                                \
                                                                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
/** nicht sicher... w und w_hat sind hier nicht vertauscht...               */\
    ths->w = (double*) nfft_malloc(ths->mv->M_total*sizeof(double));          \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->w_hat = (double*) nfft_malloc(ths->mv->N_total*sizeof(double));      \
}                                                                             \
                                                                              \
                                                                              \
/** void adjoint <mv>_init */                                                 \
F(MV, FLT, FLT_TYPE, init, i ## MV ## _adjoint_plan *ths, MV ## _plan *mv)     \
{                                                                             \
  adjoint ## MV ## _init_advanced(ths, mv, CGNR);                             \
} /* void adjoint<mv>_init */                                                 \
                                                                              \
/** void adjoint <mv>_before_loop */                                          \
F(MV, FLT, FLT_TYPE, before_loop,   i ## MV ## _adjoint_plan *ths)             \
{                                                                             \
/**  nfft_cp_ ## FLT(ths->mv->f_hat, ths->f_hat_iter, ths->mv->N_total); */   \
                                                                              \
   nfft_cp_ ## FLT(ths->mv->f, ths->f_iter, ths->mv->M_total);                \
   /** Vertauschung der Zeiger */                                             \
                                                                              \
   NFFT_SWAP_ ## FLT(ths->r_hat_iter, ths->mv->f_hat);                        \
   MV ## _adjoint(ths->mv);                                                    \
   NFFT_SWAP_ ## FLT(ths->r_hat_iter, ths->mv->f_hat);                        \
                                                                              \
                                                                              \
   nfft_upd_axpy_ ## FLT(ths->r_hat_iter, -1.0, ths->y_hat, ths->mv->N_total);\
                                                                              \
 if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))        \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_WEIGHT)                                      \
  ths->dot_r_hat_iter =                                                       \
      nfft_dot_w_ ## FLT(ths->r_hat_iter, ths->w, ths->mv->N_total);          \
      else                                                                    \
  ths->dot_r_hat_iter = nfft_dot_ ## FLT(ths->r_hat_iter, ths->mv->N_total);  \
    }                                                                         \
                                                                              \
                                                                              \
/*-----------------*/                                                         \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    nfft_cp_w_ ## FLT(ths->mv->f_hat, ths->w_hat, ths->r_hat_iter, ths->mv->N_total);     \
  else                                                                          \
    nfft_cp_ ## FLT(ths->mv->f_hat, ths->r_hat_iter, ths->mv->N_total);         \
                                                                                \
/** NFFT_SWAP_ ## FLT(ths->z_hat_iter, ths->mv->f_hat);                         \
    MV ## _adjoint(ths->mv);                                                    \
    NFFT_SWAP_ ## FLT(ths->z_hat_iter, ths->mv->f_hat);                       */\
                                                                                \
  if((!(ths->flags & LANDWEBER)) || (ths->flags & NORMS_FOR_LANDWEBER))       \
    {                                                                         \
      if(ths->flags & PRECOMPUTE_DAMP)                                        \
        ths->dot_z_iter =                                                     \
          nfft_dot_w_ ## FLT(ths->z_iter, ths->w, ths->mv->M_total);          \
      else                                                                    \
         ths->dot_z_iter = nfft_dot_ ## FLT(ths->z_iter, ths->mv->M_total);   \
    }                                                                         \
                                                                              \
  if(ths->flags & CGNE)                                                       \
    ths->dot_p_iter = ths->dot_z_iter;                                        \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    nfft_cp_ ## FLT(ths->p_iter, ths->z_iter, ths->mv->M_total);              \
} /* void i<mv>_before_loop */                                                \
                                                                              \
/** neu */                                                                    \
                                                                              \
/** void i<mv>_loop_one_step_cgne */                                          \
F(MV, FLT, FLT_TYPE, loop_one_step_cgne, i ## MV ## _adjoint_plan *ths)       \
{                                                                             \
  /**    alpha berechnen                            */            \
  ths->alpha_iter = ths->dot_r_hat_iter / ths->dot_p_iter;                    \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    nfft_upd_xpawy_ ## FLT(ths->f_iter, ths->alpha_iter, ths->w,              \
          ths->p_iter, ths->mv->M_total);                                     \
  else                                                                        \
    nfft_upd_xpay_ ## FLT(ths->f_ter, ths->alpha_iter, ths->p_iter,           \
                          ths->mv->M_total);                                  \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    nfft_cp_w_ ## FLT(ths->mv->f, ths->w, ths->p_iter,                        \
                      ths->mv->M_total);                                      \
  else                                                                        \
    nfft_cp_ ## FLT(ths->mv->f, ths->p_iter, ths->mv->M_total);               \
                                                                              \
  /** anstatt trafo adjoint.... MV ## _trafo(ths->mv);                      */\
   MV ## _adjoint(ths->mv);                                                   \
                                                                              \
  nfft_upd_xpay_ ## FLT(ths->r_hat_iter, -ths->alpha_iter, ths->mv->f_hat,    \
                        ths->mv->N_total);                                    \
                                                                              \
  ths->dot_r_hat_iter_old = ths->dot_r_hat_iter;                              \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    ths->dot_r_hat_iter = nfft_dot_w_ ## FLT(ths->r_hat_iter, ths->w_hat, \
                                         ths->mv->N_total);                   \
  else                                                                        \
    ths->dot_r_hat_iter = nfft_dot_ ## FLT(ths->r_hat_iter, ths->mv->N_total);\
                                                                              \
  /*-----------------*/                                                       \
  ths->beta_iter = ths->dot_r_hat_iter / ths->dot_r_hat_iter_old;             \
                                                                              \
  /*-----------------*/                                                       \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    nfft_cp_w_ ## FLT(ths->mv->f_hat, ths->w_hat, ths->r_hat_iter, ths->mv->N_total);     \
  else                                                                        \
    nfft_cp_ ## FLT(ths->mv->f_hat, ths->r_hat_iter, ths->mv->N_total);       \
                                                                              \
  /** diesmal nicht adjoint, sondern trafo  MV ## _adjoint(ths->mv);        */\
    MV ## _trafo(ths->mv);                                                    \
                                                                              \
                                                                              \
                                                                              \
                                                                              \
                                                                              \
  nfft_upd_axpy_ ## FLT(ths->p_iter, ths->beta_iter, ths->mv->f,      \
       ths->mv->M_total);                                                     \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    ths->dot_p_hat_iter =                                                     \
      nfft_dot_w_ ## FLT(ths->p_hat_iter, ths->w_hat, ths->mv->N_total);      \
  else                                                                        \
    ths->dot_p_hat_iter = nfft_dot_ ## FLT(ths->p_hat_iter, ths->mv->N_total);\
} /* void i<mv>_adjoint_loop_one_step_cgne */                                 \
                                                                              \
                                                                              \
                                                                              \
                                                                              \
                                                                              \
                                                                              \
/** void i<mv>_adjoint_loop_one_step */                                       \
F(MV, FLT, FLT_TYPE, adjoint_loop_one_step, i ## MV ## _adjoint_plan *ths)    \
{                                                                             \
/**  angepasst, aber noch auskommentiert, bis die methoden auch angepasst     \
 *   wurden                                                                   \
 * if(ths->flags & LANDWEBER)                                                 \
 *   i ## MV ## _adjoint_loop_one_step_landweber(ths);                        \
 *                                                                            \
 * if(ths->flags & STEEPEST_DESCENT)                                          \
 *   i ## MV ## _adjoint_loop_one_step_steepest_descent(ths);                 \
 *                                                                            \
 * if(ths->flags & CGNR)                                                      \
 *   i ## MV ## _adjoint_loop_one_step_cgnr(ths);                             \
 **/                                                                          \
  if(ths->flags & CGNE)                                                       \
    i ## MV ## _adjoint_loop_one_step_cgne(ths);                              \
} /* void i<mv>_adjoint_loop_one_step */                                      \
                                                                              \
/** void i<mv>_adjoint_finalize */                                            \
F(MV, FLT, FLT_TYPE, adjoint_finalize, i ## MV ## _adjoint_plan *ths)         \
{                                                                             \
  if(ths->flags & PRECOMPUTE_WEIGHT)                                          \
    nfft_free(ths->w);                                                        \
                                                                              \
  if(ths->flags & PRECOMPUTE_DAMP)                                            \
    nfft_free(ths->w_hat);                                                    \
                                                                              \
  if(ths->flags & CGNR)                                                       \
    {                                                                         \
      nfft_free(ths->v_hatiter);                                              \
      nfft_free(ths->z_iter);                                                 \
    }                                                                         \
                                                                              \
  if(ths->flags & STEEPEST_DESCENT)                                           \
    nfft_free(ths->v_hat_iter);                                               \
                                                                              \
  nfft_free(ths->p_iter);                                                     \
  nfft_free(ths->f_iter);                                                     \
                                                                              \
  nfft_free(ths->r_hat_iter);                                                 \
  nfft_free(ths->y);                                                          \
} /** void i<mv>adjoint_finalize */                                           \
MACRO_SOLVER_ADJOINT_IMPL(nfsft, complex, double _Complex)

