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

/**
 * Library.
 * Includes simple and fast computation of the NFST (direct problem)
 * author: S. Klatt 2004-2006
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "nfft3util.h"
#include "nfft3.h"
#include "infft.h"

/**
 *  handy shortcuts
 **/
#define NFST_DEFAULT_FLAGS   PRE_PHI_HUT|\
                             PRE_PSI|\
                             MALLOC_X|\
                             MALLOC_F_HAT|\
                             MALLOC_F|\
                             FFTW_INIT|\
                             FFT_OUT_OF_PLACE

#define FFTW_DEFAULT_FLAGS   FFTW_ESTIMATE|\
                             FFTW_DESTROY_INPUT

#define NFST_SUMMANDS ( 2 * ths->m + 2)
#define NODE(p,r) ( ths->x[(p) * ths->d + (r)])

#define MACRO_ndst_init_result_trafo      \
  memset( f, 0, ths->M_total * sizeof( double));
#define MACRO_ndst_init_result_adjoint \
  memset( f_hat, 0, ths->N_total * sizeof( double));


#define MACRO_nfst_D_init_result_A        \
  memset(g_hat, 0, nfst_prod_minus_a_int( ths->n, 0, ths->d) * sizeof( double));
#define MACRO_nfst_D_init_result_T        \
  memset(f_hat, 0, ths->N_total * sizeof( double));

#define MACRO_nfst_B_init_result_A        \
  memset(f, 0,     ths->M_total * sizeof( double));
#define MACRO_nfst_B_init_result_T        \
  memset(g, 0, nfst_prod_minus_a_int( ths->n, 0, ths->d) * sizeof( double));


#define NFST_PRE_WINFUN( d)  ths->N[d] = 2 * ths->N[d];         \
                             ths->n[d] = nfst_fftw_2N( ths->n[d]);

#define NFST_POST_WINFUN( d) ths->N[d] = (LRINT(0.5 * ths->N[d]));   \
                             ths->n[d] = nfst_fftw_2N_rev( ths->n[d]);


#define NFST_WINDOW_HELP_INIT  WINDOW_HELP_INIT


double nfst_phi_hut( nfst_plan *ths, int k, int d)
{
  NFST_PRE_WINFUN( d);
  double phi_hut_tmp = PHI_HUT( k, d);
  NFST_POST_WINFUN( d);

  return phi_hut_tmp;
}

double nfst_phi( nfst_plan *ths, double x, int d)
{
  NFST_PRE_WINFUN( d);
  double phi_tmp = PHI( x, d);
  NFST_POST_WINFUN( d);

  return phi_tmp;
}

int nfst_fftw_2N( int n)
{
  return 2 * ( n + 1);
}

int nfst_fftw_2N_rev( int n)
{
  div_t n_div;

  n_div = div(n, 2);
  return n_div.quot - 1;
}

#define MACRO_with_sin_vec     sin_vec[t][ka[t]]
#define MACRO_without_sin_vec  sin( 2.0 * PI * (ka[t]+1) * NODE(j,t))


#define MACRO_with_PRE_PHI_HUT     ths->c_phi_inv[t][kg[t]];
#define MACRO_compute_PHI_HUT_INV  (1.0 / (nfst_phi_hut( ths, kg[t]+1, t)))

#define MACRO_with_PRE_PSI     ths->psi[(j * ths->d + t) * NFST_SUMMANDS + lc[t]];
#define MACRO_compute_PSI      \
  nfst_phi( ths, NODE(j,t) - (( double)(lc[t] + lb[t])) / nfst_fftw_2N( ths->n[t]), t)



/** direct computation of non equispaced sine transforms
 *  ndst_trafo,  ndst_adjoint_
 *  require O(M N^d) arithemtical operations
 *
 * direct computation of the ndst_trafo, formula (1.1)
 * ndst_trafo:
 * for j=0,...,M-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * sin(2 (pi) k x[j])
 *
 * direct computation of the ndft_adjoint, formula (1.2)
 * ndst_adjoint:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * sin(2 (pi) k x[j])
 */

#define MACRO_ndst_malloc__sin_vec                                              \
                                                                                \
  double **sin_vec;                                                             \
  sin_vec = (double**)nfft_malloc( ths->d * sizeof( double*));                       \
  for( t = 0; t < ths->d; t++)                                                  \
    sin_vec[t] = (double*)nfft_malloc( ( ths->N[t] - 1) * sizeof( double));          \




#define MACRO_ndst_free__sin_vec                                                \
{                                                                               \
  /* free allocated memory */                                                   \
  for( t = 0; t < ths->d; t++)                                                  \
  nfft_free( sin_vec[t]);                                                          \
  nfft_free( sin_vec);                                                               \
}



#define MACRO_ndst_init__sin_vec                                                \
{                                                                               \
  for( t = 0; t < ths->d; t++)                                                  \
  {                                                                             \
    cos_x[t] = cos( 2.0 * PI * NODE(j,t));                                      \
    sin_vec[t][0] = sin( 2.0 * PI * NODE(j,t));                                 \
    sin_vec[t][1] = sin( 4.0 * PI * NODE(j,t));                                 \
    for( k = 2; k < ths->N[t] - 1; k++)                                         \
      sin_vec[t][k] = 2.0 * cos_x[t] * sin_vec[t][k-1]                          \
                      - sin_vec[t][k-2];                                        \
  }                                                                             \
}


#define MACRO_ndst_init__k__sin_k( which_one)                                   \
{                                                                               \
  sin_k[0] = 1.0;                                                               \
  for( t = 0; t < ths->d; t++)                                                  \
    ka[t] = 0;                                                                  \
                                                                                \
  for( t = 0; t < ths->d; t++)                                                  \
  {                                                                             \
    sin_k[t+1] = sin_k[t] * MACRO_ ##which_one;                                 \
  }                                                                             \
}


#define MACRO_ndst_count__k__sin_k( which_one)                                  \
{                                                                               \
  ka[ths->d-1]++;                                                               \
  i = ths->d - 1;                                                               \
  while( ( ka[i] == ths->N[i] - 1) && ( i > 0))                                 \
  {                                                                             \
    ka[i - 1]++;                                                                \
    ka[i] = 0;                                                                  \
                                                                                \
    i--;                                                                        \
  }                                                                             \
  for( t = i; t < ths->d; t++)                                                  \
    sin_k[t+1] = sin_k[t] * MACRO_ ##which_one;                                 \
}


#define MACRO_ndst_compute__trafo                                               \
{                                                                               \
  f[j] += f_hat[k] * sin_k[ths->d];                                             \
}

#define MACRO_ndst_compute__adjoint                                             \
{                                                                               \
  f_hat[k] += f[j] * sin_k[ths->d];                                             \
}


/* slow (trafo) transform */
#define MACRO_ndst( which_one)                                                  \
  void ndst_ ## which_one ( nfst_plan *ths)                                     \
  {                                                                             \
    int j, k, t, i;                                                             \
    int ka[ths->d];                                                             \
    double sin_k[ths->d+1];                                                     \
    double cos_x[ths->d];                                                       \
                                                                                \
    double *f     = ths->f;                                                     \
    double *f_hat = ths->f_hat;                                                 \
                                                                                \
    MACRO_ndst_init_result_ ## which_one;                                       \
                                                                                \
    if( ths->d == 1)                                                            \
      for( j = 0; j < ths->M_total; j++)                                        \
        for( k = 0; k < ths->N_total; k++)                                      \
        {                                                                       \
          sin_k[ths->d] = sin( 2.0 * PI * (k+1) * NODE(j,0));                   \
          MACRO_ndst_compute__ ## which_one;                                    \
        }                                                                       \
    else                                                                        \
      if( 1 == 0) /*FIXME: remove slow slow ... */                              \
        /* slow ndst */                                                         \
        for( j = 0; j < ths->M_total; j++)                                      \
        {                                                                       \
          MACRO_ndst_init__k__sin_k(without_sin_vec);                           \
                                                                                \
          for( k = 0; k < ths->N_total; k++)                                    \
          {                                                                     \
            MACRO_ndst_compute__ ## which_one;                                  \
                                                                                \
            MACRO_ndst_count__k__sin_k(without_sin_vec);                        \
          }                                                                     \
        }                                                                       \
      else                                                                      \
      {                                                                         \
        /* fast ndst_trafo */                                                   \
        MACRO_ndst_malloc__sin_vec;                                             \
                                                                                \
        for( j = 0; j < ths->M_total; j++)                                      \
        {                                                                       \
          MACRO_ndst_init__sin_vec;                                             \
                                                                                \
          MACRO_ndst_init__k__sin_k(with_sin_vec);                              \
                                                                                \
          for( k = 0; k < ths->N_total; k++)                                    \
          {                                                                     \
            MACRO_ndst_compute__ ## which_one;                                  \
                                                                                \
            MACRO_ndst_count__k__sin_k(with_sin_vec);                           \
          }                                                                     \
        }                                                                       \
        MACRO_ndst_free__sin_vec;                                               \
      }                                                                         \
  } /* ndst_{trafo, adjoint} */


MACRO_ndst(trafo)
MACRO_ndst(adjoint)




/** fast computation of non equispaced sine transforms
 *  require O(N^d log(N) + M) arithemtical operations
 *
 * fast computation of the nfst_trafo, formula (1.1)
 * nfst_trafo:
 * for j=0,...,M-1
 *  f[j] = sum_{k in I_N^d} f_hat[k] * sin(2 (pi) k x[j])
 *
 * direct computation of the nfst_adjoint, formula (1.2)
 * nfst_adjoint:
 * for k in I_N^d
 *  f_hat[k] = sum_{j=0}^{M-1} f[j] * sin(2 (pi) k x[j])
 */

#define MACRO_nfst__lower_boundary( j,act_dim)                                  \
{                                                                               \
  lb[(act_dim)] =                                                               \
    (LRINT(NODE((j),(act_dim)) * nfst_fftw_2N( ths->n[(act_dim)]))) - ths->m;   \
}

#define MACRO_nfst_D_compute_A                                                  \
{                                                                               \
  g_hat[kg_plain[ths->d]] = f_hat[k_L] * c_phi_inv_k[ths->d];                   \
}

#define MACRO_nfst_D_compute_T                                                  \
{                                                                               \
  f_hat[k_L] = g_hat[kg_plain[ths->d]] * c_phi_inv_k[ths->d];                   \
}


#define MACRO_init__kg                                                          \
{                                                                               \
  for( t = 0; t < ths->d; t++)                                                  \
    kg[t] = 0;                                                                  \
                                                                                \
  i = 0;                                                                        \
}


#define MACRO_count__kg                                                         \
{                                                                               \
  kg[ths->d - 1]++;                                                             \
  i = ths->d - 1;                                                               \
  while( ( kg[i] == ths->N[i] - 1) && ( i > 0))                                 \
  {                                                                             \
    kg[i - 1]++;                                                                \
    kg[i] = 0;                                                                  \
                                                                                \
    i--;                                                                        \
  }                                                                             \
}


#define MACRO_update__c_phi_inv_k__lg_plain( which_one, which_phi)              \
{                                                                               \
  for( t = i; t < ths->d; t++) {                                                \
    MACRO__c_phi_inv_k( which_phi);                                             \
    kg_plain[t+1] = kg_plain[t] * ths->n[t] + kg[t];                            \
  }                                                                             \
}


#define MACRO__c_phi_inv_k( which_phi)                                          \
{                                                                               \
  c_phi_inv_k[t+1] = 0.5 * c_phi_inv_k[t] * MACRO_ ## which_phi;                \
}


#define MACRO_nfst_D(which_one)                                                 \
static inline void nfst_D_ ## which_one (nfst_plan *ths)                               \
{                                                                               \
  int k_L;                              /**< plain index                */      \
                                                                                \
  int i, t;                                                                     \
  int kg[ths->d];                       /**< multi index in g_hat,c_phi */      \
  double c_phi_inv_k[ths->d+1];         /**< postfix product of PHI_HUT */      \
  int kg_plain[ths->d+1];               /**< postfix plain index        */      \
                                                                                \
  double *g_hat, *f_hat;                /**< local copy                 */      \
                                                                                \
  g_hat = ths->g_hat;                                                           \
  f_hat = ths->f_hat;                                                           \
                                                                                \
  MACRO_nfst_D_init_result_ ## which_one                                        \
                                                                                \
  c_phi_inv_k[0] = 1;                                                           \
  kg_plain[0]    = 0;                                                           \
                                                                                \
  MACRO_init__kg;                                                               \
                                                                                \
  if( ths->nfst_flags & PRE_PHI_HUT)                                            \
                                                                                \
    for( k_L = 0; k_L < ths->N_total; k_L++)                                    \
    {                                                                           \
      MACRO_update__c_phi_inv_k__lg_plain( which_one, with_PRE_PHI_HUT);        \
                                                                                \
      MACRO_nfst_D_compute_ ## which_one;                                       \
                                                                                \
      MACRO_count__kg;                                                          \
                                                                                \
    } /* for(k_L) */                                                            \
                                                                                \
  else                                                                          \
                                                                                \
    for( k_L = 0; k_L < ths->N_total; k_L++)                                    \
    {                                                                           \
      MACRO_update__c_phi_inv_k__lg_plain( which_one, compute_PHI_HUT_INV);     \
                                                                                \
      MACRO_nfst_D_compute_ ## which_one;                                       \
                                                                                \
      MACRO_count__kg                                                           \
                                                                                \
    } /* for(k_L) */                                                            \
} /* nfst_D */

MACRO_nfst_D(A)
MACRO_nfst_D(T)







/** sub routines for the fast transforms
 *  matrix vector multiplication with \f$B, B^{\rm T}\f$
 */
#define MACRO_nfst_B_PRE_FULL_PSI_compute_A                                     \
{                                                                               \
  (*fj) += ths->psi[ix] * g[ths->psi_index_g[ix]];                              \
}

#define MACRO_nfst_B_PRE_FULL_PSI_compute_T                                     \
{                                                                               \
  g[ths->psi_index_g[ix]] += ths->psi[ix] * (*fj);                              \
}



#define MACRO_nfst_B_compute_A                                                  \
{                                                                               \
  (*fj) += phi_tilde[ths->d] * g[lg_plain[ths->d]];                             \
}

#define MACRO_nfst_B_compute_T                                                  \
{                                                                               \
  g[lg_plain[ths->d]] += phi_tilde[ths->d] * (*fj);                             \
}



#define MACRO_compute_lg_offset__count_lg( i0)                                  \
{                                                                               \
  /* determine index in g-array corresponding to lb[(i0)] */                    \
  if( lb[(i0)] < 0)                                                             \
  {                                                                             \
    lg_offset[(i0)] =                                                           \
      (lb[(i0)] % nfst_fftw_2N( ths->n[(i0)])) + nfst_fftw_2N( ths->n[(i0)]);   \
  }                                                                             \
  else                                                                          \
  {                                                                             \
    lg_offset[(i0)] = lb[(i0)] % nfst_fftw_2N( ths->n[(i0)]);                   \
  }                                                                             \
                                                                                \
  if( lg_offset[(i0)] > ths->n[(i0)]+1)                                         \
      lg_offset[(i0)] = -( nfst_fftw_2N( ths->n[(i0)]) - lg_offset[(i0)]);      \
}



#define MACRO_set__lg__to__lg_offset                                            \
{                                                                               \
  if( lg_offset[i] <= 0)                                                        \
  {                                                                             \
    lg[i] = -lg_offset[i];                                                      \
    count_lg[i] = -1;                                                           \
  }                                                                             \
  else                                                                          \
  {                                                                             \
    lg[i] = +lg_offset[i];                                                      \
    count_lg[i] = +1;                                                           \
  }                                                                             \
}



#define MACRO_count__lg(dim)                                                    \
{                                                                               \
  /* turn around when we hit one of the boundaries */                           \
  if( ((lg[(dim)] == 0) || (lg[(dim)] == (ths->n[(dim)] + 1))) )                \
    count_lg[(dim)] *= -1;                                                      \
                                                                                \
  lg[(dim)] += count_lg[(dim)];                                                 \
}



#define MACRO_init_lb_lg_lc_phi_tilde_lg_plain( which_psi)                      \
{                                                                               \
  for( i = 0; i < ths->d; i++)                                                  \
  {                                                                             \
    MACRO_nfst__lower_boundary( j, i);                                          \
                                                                                \
    MACRO_compute_lg_offset__count_lg( i);                                      \
    MACRO_set__lg__to__lg_offset;                                               \
                                                                                \
    /* counter for lg */                                                        \
    lc[i] = 0;                                                                  \
  }                                                                             \
                                                                                \
  for( t = 0; t < ths->d; t++)                                                  \
  {                                                                             \
    if( lg[t] == 0)                                                             \
    {                                                                           \
      lg_plain[t+1]  = lg_plain[t] * ths->n[t];                                 \
      phi_tilde[t+1] = 0.0;                                                     \
    }                                                                           \
    else                                                                        \
      if( lg[t] == ths->n[t]+1)                                                 \
      {                                                                         \
        lg_plain[t+1]  = lg_plain[t] * ths->n[t] + ths->n[t]-1;                 \
        phi_tilde[t+1] = 0.0;                                                   \
      }                                                                         \
      else                                                                      \
      {                                                                         \
        MACRO__phi_tilde( which_psi);                                           \
        lg_plain[t+1]  = lg_plain[t] * ths->n[t] + lg[t]-1;                     \
      }                                                                         \
  }                                                                             \
                                                                                \
  i = 0;                                                                        \
}



#define MACRO_count__lg_lc                                                      \
{                                                                               \
  MACRO_count__lg( ths->d-1);                                                   \
                                                                                \
  lc[ths->d - 1]++;                                                             \
  i = ths->d - 1;                                                               \
                                                                                \
  while( (lc[i] == NFST_SUMMANDS) && (i > 0))                                   \
  {                                                                             \
    lc[i - 1]++;                                                                \
    lc[i] = 0;                                                                  \
                                                                                \
    /* ansonsten lg[i-1] verschieben */                                         \
    MACRO_count__lg( i - 1);                                                    \
    /* lg[i] = anfangswert */                                                   \
    MACRO_set__lg__to__lg_offset;                                               \
                                                                                \
    i--;                                                                        \
  }                                                                             \
}


#define MACRO_update__phi_tilde__lg_plain( which_psi)                           \
{                                                                               \
  for( t = i; t < ths->d; t++)                                                  \
  {                                                                             \
    if( (lg[t] != 0) && (lg[t] != ths->n[t]+1))                                 \
    {                                                                           \
      MACRO__phi_tilde( which_psi);                                             \
      lg_plain[t+1] = lg_plain[t] * ths->n[t] + lg[t]-1;                        \
    }                                                                           \
    else                                                                        \
      phi_tilde[t+1] = 0.0;                                                     \
  }                                                                             \
}



#define MACRO__phi_tilde( which_psi)                                            \
{                                                                               \
  phi_tilde[t+1] = (double)count_lg[t] * phi_tilde[t] * MACRO_ ## which_psi;    \
}




#define MACRO_nfst_B( which_one)                                                \
  static inline void nfst_B_ ## which_one ( nfst_plan *ths)                            \
  { /* MACRO_nfst_B */                                                          \
    int lb[ths->d];                     /**< multi band with respect to x_j */  \
    int j, t, i;                        /**< index nodes, help vars         */  \
    int lprod, l_L, ix;                 /**< index one row of B             */  \
    int lc[ths->d];                     /**< multi index 0<=lj<u+o+1        */  \
    int lg[ths->d];                     /**< real index of g in array       */  \
    int lg_offset[ths->d];              /**< offset in g according to u     */  \
    int count_lg[ths->d];               /**< count summands (2m+2)          */  \
    int lg_plain[ths->d+1];             /**< index of g in multi_array      */  \
    double *f, *g;                      /**< local copy                     */  \
    double phi_tilde[ths->d+1];         /**< holds values for psi           */  \
    double *fj;                         /**< pointer to final result        */  \
                                                                                \
    f = ths->f; g = ths->g;                                                     \
                                                                                \
    MACRO_nfst_B_init_result_ ## which_one                                      \
                                                                                \
    /* both flags are set */                                                    \
    if( (ths->nfst_flags & PRE_PSI) && (ths->nfst_flags & PRE_FULL_PSI))        \
    {                                                                           \
      for( ix = 0, j = 0, fj = &f[0]; j < ths->M_total; j++, fj += 1)           \
        for( l_L = 0; l_L < ths->psi_index_f[j]; l_L++, ix++)                   \
        {                                                                       \
          MACRO_nfst_B_PRE_FULL_PSI_compute_ ## which_one;                      \
        }                                                                       \
    }                                                                           \
    else                                                                        \
    {                                                                           \
      phi_tilde[0] = 1;                                                         \
      lg_plain[0]  = 0;                                                         \
                                                                                \
      for( t = 0, lprod = 1; t < ths->d; t++)                                   \
        lprod *= NFST_SUMMANDS;                                                 \
                                                                                \
      /* PRE_PSI flag is set */                                                 \
      if( ths->nfst_flags & PRE_PSI)                                            \
      {                                                                         \
        for( j = 0, fj = &f[0]; j < ths->M_total; j++, fj += 1)                 \
        {                                                                       \
          MACRO_init_lb_lg_lc_phi_tilde_lg_plain( with_PRE_PSI);                \
                                                                                \
          for( l_L = 0; l_L < lprod; l_L++)                                     \
          {                                                                     \
            MACRO_update__phi_tilde__lg_plain( with_PRE_PSI);                   \
                                                                                \
            MACRO_nfst_B_compute_ ## which_one;                                 \
                                                                                \
            MACRO_count__lg_lc;                                                 \
                                                                                \
          } /* for( l_L) */                                                     \
        } /* for( j) */                                                         \
      } /* if( PRE_PSI) */                                                      \
                                                                                \
      /* no PSI flag is set */                                                  \
      else                                                                      \
      {                                                                         \
        for( j = 0, fj = &f[0]; j < ths->M_total; j++, fj += 1)                 \
        {                                                                       \
          MACRO_init_lb_lg_lc_phi_tilde_lg_plain( compute_PSI);                 \
                                                                                \
          for( l_L = 0; l_L < lprod; l_L++)                                     \
          {                                                                     \
            MACRO_update__phi_tilde__lg_plain( compute_PSI);                    \
                                                                                \
            MACRO_nfst_B_compute_ ## which_one;                                 \
                                                                                \
            MACRO_count__lg_lc;                                                 \
                                                                                \
          } /* for(l_L) */                                                      \
        } /* for(j) */                                                          \
      } /* else(PRE_PSI) */                                                     \
    }/* else( PRE_PRE && FULL_PRE_PSI) */                                       \
} /* nfst_B */

MACRO_nfst_B(A)
MACRO_nfst_B(T)






/**
 * user routines
 *
 */
void nfst_trafo( nfst_plan *ths)
{
  /**
   * use ths->my_fftw_r2r_plan
   *
   */
  ths->g_hat = ths->g1;
  ths->g     = ths->g2;


  /**
   * form \f$ \hat g_k = \frac{\hat f_k}{c_k\left(\phi\right)} \text{ for }
   * k \in I_N \f$
   *
   */
  TIC(0)
  nfst_D_A( ths);
  TOC(0)


  /**
   * compute by d-variate discrete Fourier transform
   * \f$ g_l = \sum_{k \in I_N} \hat g_k {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
   * \text{ for } l \in I_n \f$
   *
   */
  TIC(1)
  fftw_execute( ths->my_fftw_r2r_plan);
  TOC(1)


  /**
   * set \f$ f_j = \sum_{l \in I_n,m(x_j)} g_l \psi\left(x_j-\frac{l}{n}\right)
   * \text{ for } j=0,\hdots,M-1 \f$
   *
   */
  TIC(2)
  nfst_B_A( ths);
  TOC(2)

} /* nfst_trafo */




void nfst_adjoint( nfst_plan *ths)
{
  /**
   * use ths->my_fftw_plan
   *
   **/
  ths->g_hat = ths->g2;
  ths->g     = ths->g1;


  /**
   * set \f$ g_l = \sum_{j=0}^{M-1} f_j \psi\left(x_j-\frac{l}{n}\right)
   * \text{ for } l \in I_n,m(x_j) \f$
   *
   */
  TIC(2)
  nfst_B_T( ths);
  TOC(2)


  /**
   * compute by d-variate discrete cosine transform
   * \f$ \hat g_k = \sum_{l \in I_n} g_l {\rm e}^{-2\pi {\rm i} \frac{kl}{n}}
   * \text{ for }  k \in I_N\f$
   *
   */
  TIC(1)
  fftw_execute( ths->my_fftw_r2r_plan);
  TOC(1)


  /**
   * form \f$ \hat f_k = \frac{\hat g_k}{c_k\left(\phi\right)} \text{ for }
   * k \in I_N \f$
   *
   */
  TIC(0)
  nfst_D_T( ths);
  TOC(0)

} /* nfst_adjoint */



/**
 * initialization of direct transform
 *
 */
void nfst_precompute_phi_hut( nfst_plan *ths)
{
  int kg[ths->d];                      /**< index over all frequencies       */
  int t;                               /**< index over all dimensions        */

  ths->c_phi_inv = (double**)nfft_malloc( ths->d * sizeof( double*));

  for( t = 0; t < ths->d; t++)
  {
    ths->c_phi_inv[t] = (double*)nfft_malloc( ( ths->N[t] - 1) * sizeof( double));

    for( kg[t] = 0; kg[t] < ths->N[t] - 1; kg[t]++)
    {
      ths->c_phi_inv[t][kg[t]] = MACRO_compute_PHI_HUT_INV;
    }
  }
} /* nfst_phi_hut */



void nfst_precompute_psi( nfst_plan *ths)
{
  int t;                                /**< index over all dimensions        */
  int j;                                /**< index over all nodes             */
  int lc[ths->d];                       /**< index 0<=lj<u+o+1                */
  int lb[ths->d];                       /**< depends on x_j                   */

  for (t = 0; t < ths->d; t++)
  {
    for(j = 0; j < ths->M_total; j++)
    {
      MACRO_nfst__lower_boundary( j, t);

      for( lc[t] = 0; lc[t] < NFST_SUMMANDS; lc[t]++)
        ths->psi[(j * ths->d + t) * NFST_SUMMANDS + lc[t]] = MACRO_compute_PSI;

    } /* for(j) */
  }  /* for(t) */

  /* full precomputation of psi */
  if ( ths->nfst_flags & PRE_FULL_PSI)
    nfst_full_psi( ths, ths->nfst_full_psi_eps);

} /* nfst_precompute_psi */



/** more memory usage, a bit faster */
void nfst_full_psi(nfst_plan *ths, double eps)
{
  int t, i;                             /**< index over all dimensions        */
  int j;                                /**< index over all nodes             */
  int l_L;                              /**< plain index 0<=l_L<lprod         */
  int lc[ths->d];                       /**< multi index 0<=lj<u+o+1          */
  int lg_plain[ths->d+1];               /**< postfix plain index              */
  int count_lg[ths->d];
  int lg_offset[ths->d];
  int lg[ths->d];
  int lprod;                            /**< 'bandwidth' of matrix B          */
  int lb[ths->d];                       /**< depends on x_j                   */

  double phi_tilde[ths->d+1];

  int *index_g, *index_f;
  double *new_psi;
  int ix, ix_old, size_psi;

  phi_tilde[0] = 1.0;
  lg_plain[0]  =   0;

  if(ths->nfst_flags & PRE_PSI)
  {
    size_psi = ths->M_total;
    index_f  =    (int*)nfft_malloc( ths->M_total  * sizeof( int));
    index_g  =    (int*)nfft_malloc( size_psi * sizeof( int));
    new_psi  = (double*)nfft_malloc( size_psi * sizeof( double));

    for( t = 0,lprod = 1; t < ths->d; t++)
    {
      lprod *= NFST_SUMMANDS;
      eps *= PHI( 0, t);
    }

    for( ix = 0, ix_old = 0, j = 0; j < ths->M_total; j++)
    {
      MACRO_init_lb_lg_lc_phi_tilde_lg_plain( with_PRE_PSI);

      for( l_L = 0; l_L < lprod; l_L++)
      {
        MACRO_update__phi_tilde__lg_plain( with_PRE_PSI);

        if( fabs(phi_tilde[ths->d]) > eps)
        {
          index_g[ix] =  lg_plain[ths->d];
          new_psi[ix] = phi_tilde[ths->d];

          ix++;
          if( ix == size_psi)
          {
            size_psi += ths->M_total;
            index_g   =    (int*)realloc( index_g, size_psi * sizeof( int));
            new_psi   = (double*)realloc( new_psi, size_psi * sizeof( double));
          }
        }
        MACRO_count__lg_lc;

      } /* for(l_L) */

      index_f[j] = ix - ix_old;
      ix_old     = ix;

    } /* for(j) */

    nfft_free( ths->psi);

    size_psi      = ix;
    ths->size_psi = size_psi;
    index_g       = (int*)realloc( index_g, size_psi * sizeof( int));
    new_psi       = (double*)realloc( new_psi, size_psi * sizeof( double));

    ths->psi         = new_psi;
    ths->psi_index_g = index_g;
    ths->psi_index_f = index_f;

  } /* if(PRE_PSI) */
} /* nfst_full_psi */




void nfst_init_help( nfst_plan *ths)
{
  int t;                                /**< index over all dimensions        */

  ths->N_total = nfst_prod_minus_a_int( ths->N, 1, ths->d);

  ths->sigma   = (double*)nfft_malloc( ths->d * sizeof( double));

  for( t = 0; t < ths->d; t++)
    /* FIXME: n/N or (n+1)/N */
    ths->sigma[t] = ((double)ths->n[t] + 1) / ths->N[t];

  /* assign r2r transform kinds for each dimension */
  ths->r2r_kind = (fftw_r2r_kind*) nfft_malloc ( ths->d * sizeof( fftw_r2r_kind));
  for (t = 0; t < ths->d; t++)
    ths->r2r_kind[t] = FFTW_RODFT00;


  WINDOW_HELP_INIT;

  if(ths->nfst_flags & MALLOC_X)
    ths->x = (double*)nfft_malloc( ths->d * ths->M_total * sizeof( double));

  if(ths->nfst_flags & MALLOC_F_HAT)
    ths->f_hat = (double*)nfft_malloc( ths->N_total * sizeof( double));

  if(ths->nfst_flags & MALLOC_F)
    ths->f = (double*)nfft_malloc( ths->M_total * sizeof( double));

  if(ths->nfst_flags & PRE_PHI_HUT)
    nfst_precompute_phi_hut( ths);

  /* NO FFTW_MALLOC HERE */
  if(ths->nfst_flags & PRE_PSI)
  {
    ths->psi =
      (double*)nfft_malloc( ths->M_total * ths->d * NFST_SUMMANDS * sizeof( double));

    /**
     * set default for full_psi_eps
     **/
    ths->nfst_full_psi_eps = pow(10, -10);
  }

  if(ths->nfst_flags & FFTW_INIT)
  {
      ths->g1 =
        (double*)nfft_malloc( nfst_prod_minus_a_int( ths->n, 0, ths->d) * sizeof( double));

      if(ths->nfst_flags & FFT_OUT_OF_PLACE)
        ths->g2 =
          (double*)nfft_malloc( nfst_prod_minus_a_int( ths->n, 0, ths->d) * sizeof( double));
      else
        ths->g2 = ths->g1;

      ths->my_fftw_r2r_plan =
        fftw_plan_r2r( ths->d, ths->n, ths->g1, ths->g2, ths->r2r_kind, ths->fftw_flags);
  }

  ths->mv_trafo = (void (*) (void* ))nfst_trafo;
  ths->mv_adjoint = (void (*) (void* ))nfst_adjoint;
}

void nfst_init( nfst_plan *ths, int d, int *N, int M_total)
{
  int t;                                /**< index over all dimensions        */

  ths->d = d;

  ths->N      = (int*)nfft_malloc( ths->d * sizeof( int));

  for(t = 0;t < d; t++)
    ths->N[t] = N[t];

  ths->n      = (int*)nfft_malloc( ths->d * sizeof( int));

  for( t = 0; t < d; t++)
    ths->n[t] = 2 * nfft_next_power_of_2( ths->N[t]) - 1;

  ths->M_total = M_total;

  WINDOW_HELP_ESTIMATE_m;

  ths->nfst_flags = NFST_DEFAULT_FLAGS;
  ths->fftw_flags = FFTW_DEFAULT_FLAGS;

  nfst_init_help( ths);
}


void nfst_init_m( nfst_plan *ths, int d, int *N, int M_total, int m)
{
  int t, n[d];

  for( t = 0; t < d; t++)
    n[t] = nfst_fftw_2N( nfft_next_power_of_2( N[t]));

  nfst_init_guru( ths, d, N, M_total, n, m, NFST_DEFAULT_FLAGS, FFTW_DEFAULT_FLAGS);
}


void nfst_init_guru( nfst_plan *ths, int d, int *N,
     int M_total, int *n, int m,
     unsigned nfst_flags, unsigned fftw_flags)
{
  int t;             /**< index over all dimensions */

  ths->d = d;
  ths->M_total = M_total;

  ths->N      = (int*)nfft_malloc( ths->d * sizeof( int));

  for( t = 0; t < d; t++)
    ths->N[t]      = N[t];

  ths->n      = (int*)nfft_malloc( ths->d * sizeof( int));

  for( t = 0; t < d; t++)
    ths->n[t]      = n[t];

  ths->m = m;

  ths->nfst_flags = nfst_flags;
  ths->fftw_flags = fftw_flags;

  nfst_init_help( ths);
}


void nfst_init_1d( nfst_plan *ths, int N0, int M_total)
{
  int N[1];

  N[0] = N0;
  nfst_init( ths, 1, N, M_total);
}

void nfst_init_2d( nfst_plan *ths, int N0, int N1, int M_total)
{
  int N[2];

  N[0] = N0;
  N[1] = N1;
  nfst_init( ths, 2, N, M_total);
}

void nfst_init_3d( nfst_plan *ths, int N0, int N1, int N2, int M_total)
{
  int N[3];

  N[0] = N0;
  N[1] = N1;
  N[2] = N2;
  nfst_init( ths, 3, N, M_total);
}

void nfst_finalize( nfst_plan *ths)
{
  int t; /* index over dimensions */

  if( ths->nfst_flags & FFTW_INIT)
  {
    fftw_destroy_plan( ths->my_fftw_r2r_plan);

    if( ths->nfst_flags & FFT_OUT_OF_PLACE)
      nfft_free( ths->g2);

    nfft_free( ths->g1);
  }

  /* NO FFTW_FREE HERE */
  if( ths->nfst_flags & PRE_PSI)
  {
    if( ths->nfst_flags & PRE_FULL_PSI)
    {
      nfft_free( ths->psi_index_g);
      nfft_free( ths->psi_index_f);
    }

    nfft_free( ths->psi);
  }

  if( ths->nfst_flags & PRE_PHI_HUT) {
    for( t = 0; t < ths->d; t++)
      nfft_free( ths->c_phi_inv[t]);
    nfft_free( ths->c_phi_inv);
  }

  if( ths->nfst_flags & MALLOC_F)
    nfft_free( ths->f);

  if( ths->nfst_flags & MALLOC_F_HAT)
    nfft_free( ths->f_hat);

  if( ths->nfst_flags & MALLOC_X)
    nfft_free( ths->x);

  WINDOW_HELP_FINALIZE;

  nfft_free( ths->N);
  nfft_free( ths->n);
  nfft_free( ths->sigma);

  nfft_free(ths->r2r_kind);
} /* nfst_finalize */

